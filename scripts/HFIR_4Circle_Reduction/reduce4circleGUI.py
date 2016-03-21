#pylint: disable=invalid-name,relative-import,W0611,R0921,R0902,R0904,R0921,C0302
################################################################################
#
# MainWindow application for reducing HFIR 4-circle
#
################################################################################
import os
import sys
import math
import csv
import time
import random
import numpy

from PyQt4 import QtCore, QtGui
try:
    from mantidqtpython import MantidQt
except ImportError as e:
    NO_SCROLL = True
else:
    NO_SCROLL = False

import reduce4circleControl as r4c
import guiutility as gutil
from peakinfo import PeakInfo
import fourcircle_utility as fcutil
import plot3dwindow

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

# import line for the UI python class
from ui_MainWindow import Ui_MainWindow


class MainWindow(QtGui.QMainWindow):
    """ Class of Main Window (top)
    """
    TabPage = {'View Raw Data': 1,
               'Calculate UB': 2,
               'UB Matrix': 3,
               'Peak Integration': 5}

    def __init__(self, parent=None):
        """ Initialization and set up
        """
        # Base class
        QtGui.QMainWindow.__init__(self,parent)

        # UI Window (from Qt Designer)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # Make UI scrollable
        if NO_SCROLL is False:
            self._scrollbars = MantidQt.API.WidgetScrollbarDecorator(self)
            self._scrollbars.setEnabled(True) # Must follow after setupUi(self)!

        # Mantid configuration
        self._instrument = str(self.ui.comboBox_instrument.currentText())
        # config = ConfigService.Instance()
        # self._instrument = config["default.instrument"]

        # Event handling definitions
        # Top
        self.connect(self.ui.pushButton_setExp, QtCore.SIGNAL('clicked()'),
                     self.do_set_experiment)

        # Tab 'Data Access'
        self.connect(self.ui.pushButton_applySetup, QtCore.SIGNAL('clicked()'),
                     self.do_apply_setup)
        self.connect(self.ui.pushButton_browseLocalDataDir, QtCore.SIGNAL('clicked()'),
                     self.do_browse_local_spice_data)
        self.connect(self.ui.pushButton_testURLs, QtCore.SIGNAL('clicked()'),
                     self.do_test_url)
        self.connect(self.ui.pushButton_ListScans, QtCore.SIGNAL('clicked()'),
                     self.do_list_scans)
        self.connect(self.ui.pushButton_downloadExpData, QtCore.SIGNAL('clicked()'),
                     self.do_download_spice_data)
        self.connect(self.ui.comboBox_mode, QtCore.SIGNAL('currentIndexChanged(int)'),
                     self.change_data_access_mode)

        # Tab 'View Raw Data'
        self.connect(self.ui.pushButton_setScanInfo, QtCore.SIGNAL('clicked()'),
                     self.do_load_scan_info)
        self.connect(self.ui.pushButton_plotRawPt, QtCore.SIGNAL('clicked()'),
                     self.do_plot_pt_raw)
        self.connect(self.ui.pushButton_prevPtNumber, QtCore.SIGNAL('clicked()'),
                     self.do_plot_prev_pt_raw)
        self.connect(self.ui.pushButton_nextPtNumber, QtCore.SIGNAL('clicked()'),
                     self.do_plot_next_pt_raw)
        self.connect(self.ui.pushButton_showPtList, QtCore.SIGNAL('clicked()'),
                     self.show_scan_pt_list)
        self.connect(self.ui.pushButton_usePt4UB, QtCore.SIGNAL('clicked()'),
                     self.do_add_peak_to_find)
        self.connect(self.ui.pushButton_addPeakNoIndex, QtCore.SIGNAL('clicked()'),
                     self.do_add_peak_no_index)

        # Tab 'calculate ub matrix'
        self.connect(self.ui.pushButton_findPeak, QtCore.SIGNAL('clicked()'),
                     self.do_find_peak)
        self.connect(self.ui.pushButton_addPeakToCalUB, QtCore.SIGNAL('clicked()'),
                     self.do_add_ub_peak)
        self.connect(self.ui.pushButton_calUB, QtCore.SIGNAL('clicked()'),
                     self.do_cal_ub_matrix)
        self.connect(self.ui.pushButton_acceptUB, QtCore.SIGNAL('clicked()'),
                     self.doAcceptCalUB)
        self.connect(self.ui.pushButton_indexUBPeaks, QtCore.SIGNAL('clicked()'),
                     self.do_index_ub_peaks)
        self.connect(self.ui.pushButton_deleteUBPeak, QtCore.SIGNAL('clicked()'),
                     self.do_del_ub_peaks)
        self.connect(self.ui.pushButton_clearUBPeakTable, QtCore.SIGNAL('clicked()'),
                     self.do_clear_ub_peaks)
        self.connect(self.ui.pushButton_resetPeakHKLs, QtCore.SIGNAL('clicked()'),
                     self.do_reset_ub_peaks_hkl)
        self.connect(self.ui.pushButton_selectAllPeaks, QtCore.SIGNAL('clicked()'),
                     self.do_select_all_peaks)
        self.connect(self.ui.pushButton_viewScan3D, QtCore.SIGNAL('clicked()'),
                     self.do_view_data_3d)
        self.connect(self.ui.pushButton_plotSelectedData, QtCore.SIGNAL('clicked()'),
                     self.do_view_data_set_3d)

        self.connect(self.ui.pushButton_refineUB, QtCore.SIGNAL('clicked()'),
                     self.do_refine_ub)
        self.connect(self.ui.pushButton_refineUBFFT, QtCore.SIGNAL('clicked()'),
                     self.do_refine_ub_fft)

        # Tab 'Setup'
        self.connect(self.ui.pushButton_useDefaultDir, QtCore.SIGNAL('clicked()'),
                     self.do_setup_dir_default)
        self.connect(self.ui.pushButton_browseLocalCache, QtCore.SIGNAL('clicked()'),
                     self.do_browse_local_cache_dir)
        self.connect(self.ui.pushButton_browseWorkDir, QtCore.SIGNAL('clicked()'),
                     self.do_browse_working_dir)
        self.connect(self.ui.comboBox_instrument, QtCore.SIGNAL('currentIndexChanged(int)'),
                     self.change_instrument_name)

        # Tab 'UB Matrix'
        self.connect(self.ui.pushButton_showUB2Edit, QtCore.SIGNAL('clicked()'),
                     self.do_show_ub_in_box)
        self.connect(self.ui.pushButton_syncUB, QtCore.SIGNAL('clicked()'),
                     self.do_sync_ub)
        self.connect(self.ui.pushButton_setUBSliceView, QtCore.SIGNAL('clicked()'),
                     self.do_set_ub_from_text)

        # Tab 'Merge'
        self.connect(self.ui.pushButton_addScanSliceView, QtCore.SIGNAL('clicked()'),
                     self.do_add_scans_merge)
        self.connect(self.ui.pushButton_mergeScans, QtCore.SIGNAL('clicked()'),
                     self.do_merge_scans)
        self.connect(self.ui.pushButton_integratePeaks, QtCore.SIGNAL('clicked()'),
                     self.do_integrate_peaks)
        self.connect(self.ui.pushButton_setupPeakIntegration, QtCore.SIGNAL('clicked()'),
                     self.do_switch_tab_peak_int)
        self.connect(self.ui.pushButton_refreshMerged, QtCore.SIGNAL('clicked()'),
                     self.do_refresh_merged_scans_table)
        self.connect(self.ui.pushButton_plotMergedScans, QtCore.SIGNAL('clicked()'),
                     self.do_view_merged_scans_3d)
        self.connect(self.ui.pushButton_showUB, QtCore.SIGNAL('clicked()'),
                     self.do_view_ub)

        # Tab 'Integrate Peaks'
        self.connect(self.ui.pushButton_integratePt, QtCore.SIGNAL('clicked()'),
                     self.do_integrate_per_pt)
        self.connect(self.ui.pushButton_integratePeak, QtCore.SIGNAL('clicked()'),
                     self.do_integrate_peaks)

        # Tab survey
        self.connect(self.ui.pushButton_survey, QtCore.SIGNAL('clicked()'),
                     self.do_survey)
        self.connect(self.ui.pushButton_saveSurvey, QtCore.SIGNAL('clicked()'),
                     self.do_save_survey)
        self.connect(self.ui.pushButton_loadSurvey, QtCore.SIGNAL('clicked()'),
                     self.do_load_survey)
        self.connect(self.ui.pushButton_viewSurveyPeak, QtCore.SIGNAL('clicked()'),
                     self.do_view_survey_peak)
        self.connect(self.ui.pushButton_addPeaksToRefine, QtCore.SIGNAL('clicked()'),
                     self.do_add_peaks_for_ub)
        self.connect(self.ui.pushButton_selectAllSurveyPeaks, QtCore.SIGNAL('clicked()'),
                     self.do_select_all_survey)
        self.connect(self.ui.pushButton_sortInfoTable, QtCore.SIGNAL('clicked()'),
                     self.do_filter_sort_survey_table)
        self.connect(self.ui.pushButton_clearSurvey, QtCore.SIGNAL('clicked()'),
                     self.do_clear_survey)

        # Menu
        self.connect(self.ui.actionExit, QtCore.SIGNAL('triggered()'),
                     self.menu_quit)

        self.connect(self.ui.actionSave_Session, QtCore.SIGNAL('triggered()'),
                     self.save_current_session)
        self.connect(self.ui.actionLoad_Session, QtCore.SIGNAL('triggered()'),
                     self.load_session)

        # Validator ... (NEXT)

        # Declaration of class variable
        # some configuration
        self._homeSrcDir = os.getcwd()
        self._homeDir = os.getcwd()

        # Control
        self._myControl = r4c.CWSCDReductionControl(self._instrument)
        self._allowDownload = True
        self._dataAccessMode = 'Download'
        self._surveyTableFlag = True
        self._ubPeakTableFlag = True

        # Sub window
        self._my3DWindow = None

        # Initial setup
        self.ui.tabWidget.setCurrentIndex(0)
        self._init_table_widgets()
        self.ui.radioButton_ubMantidStyle.setChecked(True)
        self.ui.lineEdit_numSurveyOutput.setText('50')
        self.ui.checkBox_loadHKLfromFile.setChecked(True)
        self.ui.checkBox_sortDescending.setChecked(False)
        self.ui.radioButton_sortByCounts.setChecked(True)

        # Tab 'Access'
        self.ui.lineEdit_url.setText('http://neutron.ornl.gov/user_data/hb3a/')
        self.ui.comboBox_mode.setCurrentIndex(0)
        self.ui.lineEdit_localSpiceDir.setEnabled(True)
        self.ui.pushButton_browseLocalDataDir.setEnabled(True)

        # QSettings
        self.load_settings()

        return

    def closeEvent(self, QCloseEvent):
        """
        Close event
        :param QCloseEvent:
        :return:
        """
        print '[QCloseEvent=]', str(QCloseEvent)
        self.menu_quit()

    def _init_table_widgets(self):
        """ Initialize the table widgets
        :return:
        """
        # UB-peak table
        # NOTE: have to call this because pyqt set column and row to 0 after __init__
        #       thus a 2-step initialization has to been adopted
        self.ui.tableWidget_peaksCalUB.setup()
        self.ui.tableWidget_ubMatrix.setup()
        self.ui.tableWidget_surveyTable.setup()
        self.ui.tableWidget_peakIntegration.setup()
        self.ui.tableWidget_mergeScans.setup()

        return

    def do_advance_to_integrate_peaks(self):
        """
        Advance from 'merge'-tab to peak integration tab
        :return:
        """
        # Check whether there is any scan merged and selected
        ret_list = self.ui.tableWidget_mergeScans.get_rows_by_state('Done')
        if len(ret_list) == 0:
            self.pop_one_button_dialog('No scan is selected for integration!')
            return
        else:
            print '[DB] Total %d rows are selected for peak integration.' % len(ret_list)

        # Switch tab
        self.ui.tabWidget.setCurrentIndex(MainWindow.TabPage['Peak Integration'])

        # Add table
        for row_index in ret_list:
            # merged_ws_name = self.ui.tableWidget_mergeScans.get_merged_ws_name(row_index)
            # status, merged_info = self._myControl.get_merged_scan_info()
            status = False
            merged_info = None
            if status is False:
                err_msg = merged_info
                self.pop_one_button_dialog(err_msg)
                return
            else:
                print '[DB] Add selected row %d. ' % row_index
                status, msg = self.ui.tableWidget_peakIntegration.append_scan(merged_info)
                if status is False:
                    self.pop_one_button_dialog(msg)
        # END-FOR

        return

    def do_integrate_per_pt(self):
        """
        Use the peak center of the merged scan. Then on each Pt. with specified radius,
        calculate the integrated peak intensity and number of counts
        :return:
        """
        # get experiment and scan number
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_exp,
                                                        self.ui.lineEdit_scanIntegratePeak])
        if not status:
            self.pop_one_button_dialog('Unable to integrate peak due to %s.' % ret_obj)
            return
        else:
            exp_number, scan_number = ret_obj

        # call myController to get peak center from all pt in scans
        if self._myControl.has_peak_info(exp_number, scan_number) is False:
            # merge scan and find peak
            if not self._myControl.has_merged_data(exp_number, scan_number):
                # merge all Pts. in scan if necessary
                status, pt_number_list = self._myControl.get_pt_numbers(exp_number, scan_number)
                assert status
                self._myControl.merge_pts_in_scan(exp_number, scan_number, pt_number_list,
                                                  'q-sample')
            # find peak
            self._myControl.find_peak(exp_number, scan_number)
        scan_peak_info = self._myControl.get_peak_info(exp_number, scan_number)
        weighted_peak_center = scan_peak_info.get_peak_centre()
        print '[INFO] Exp %d Scan %d: weighted peak center = %s.' % (exp_number, scan_number,
                                                                     str(weighted_peak_center))

        # get peak radius
        status, peak_radius = gutil.parse_float_editors(self.ui.lineEdit_peakRadius)
        if not status:
            self.pop_one_button_dialog('Unable to get valid peak radius from GUI.')
            return
        else:
            assert peak_radius > 0, 'Peak radius cannot be zero or negative!'

        # call myController to integrate on each individual Pt.
        self._myControl.integrate_scan_peaks(exp_number, scan_number, peak_radius,
                                             peak_centre=weighted_peak_center,
                                             merge=False)

        # plot all the result to the table
        vec_x, vec_y = self._myControl.get_peaks_integrated_intensities(exp_number, scan_number, None)

        self.ui.graphicsView_integratedPeakView.clear_all_lines()
        self.ui.graphicsView_integratedPeakView.add_plot_1d(vec_x, vec_y)

        return

    def do_integrate_peaks(self):
        """ Integrate selected peaks with.  If any scan is not merged, then it will merge the scan first
        Integrate peaks
        :return:
        """
        # get rows to merge
        row_number_list = self.ui.tableWidget_mergeScans.get_selected_rows(True)
        if len(row_number_list) == 0:
            self.pop_one_button_dialog('No scan is selected for scan')
            return

        # get the parameters for integration
        status, peak_radius = gutil.parse_float_editors([self.ui.lineEdit_peakRadius], allow_blank=False)
        if status is False:
            self.pop_one_button_dialog('Peak radius is not given right')
            return

        # merged workspace base name
        use_default_merge_name = self.ui.checkBox_useDefaultMergedName.isChecked()

        # integrate peak
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_exp], allow_blank=False)
        if status:
            exp_number = ret_obj[0]
        else:
            self.pop_one_button_dialog('Unable to get valid experiment number due to %s.' % ret_obj)
            return

        for row_number in row_number_list:
            scan_number = self.ui.tableWidget_mergeScans.get_scan_number(row_number)
            pt_number_list = self._myControl.get_pt_numbers(exp_number, scan_number)

            # merge if required
            merged = self.ui.tableWidget_mergeScans.get_merged_status(row_number)
            if merged is False:
                self._myControl.merge_pts_in_scan(exp_no=exp_number, scan_no=scan_number, pt_num_list=pt_number_list,
                                                  target_frame='qsample')
                self.ui.tableWidget_mergeScans.set_status_by_row(row_number, 'Done')
            # END-IF

            # integrate peak
            peak_intensity = self._myControl.integrate_peak(exp_number, scan_number, pt_number_list, peak_radius)
            self.ui.tableWidget_mergeScans.set_peak_intensity(row_number, peak_intensity)

        # END-FOR

        return

    def _how_to_deal_with_this(self):
        # Get peak integration parameters
        line_editors = [self.ui.lineEdit_peakRadius,
                        self.ui.lineEdit_bkgdInnerR,
                        self.ui.lineEdit_bkgdOuterR]
        status, value_list = gutil.parse_float_editors(line_editors)
        if status is False:
            err_msg = value_list
            print '[DB] Error message: %s' % err_msg
            return
        else:
            peak_radius = value_list[0]
            bkgd_inner_radius = value_list[1]
            bkgd_outer_radius = value_list[2]

        # Get peak integration options
        adapt_q_bkgd = self.ui.checkBox_adaptQBkgd.isChecked()
        integrate_on_edge = self.ui.checkBox_integrateOnEdge.isChecked()
        print '[DB-NEXT] Unused widgets: AdaptQBackground (%s) and IntegrateOnEdge (%s).' % (
            str(adapt_q_bkgd), str(integrate_on_edge)
        )
        is_cylinder = self.ui.checkBox_cylinder.isChecked()

        # Choose the peaks to be integrated
        row_index_list = self.ui.tableWidget_peakIntegration.get_selected_rows()

        for i_row in row_index_list:
            md_ws_name = self.ui.tableWidget_peakIntegration.get_md_ws_name(i_row)
            exp_num = None
            scan_num = self.ui.tableWidget_peakIntegration.get_scan_number(i_row)
            pt_list = None
            self._myControl.integrate_peaks(exp_num, scan_num, pt_list, md_ws_name,
                                            peak_radius, bkgd_inner_radius, bkgd_outer_radius,
                                            is_cylinder)
        # END-FOR

    def do_refine_ub(self):
        """
        Refine UB matrix
        :return:
        """
        peak_info_list = self._build_peak_info_list()
        set_hkl_int = self.ui.checkBox_roundHKLInt.isChecked()

        # Refine UB matrix
        try:
            self._myControl.refine_ub_matrix_indexed_peaks(peak_info_list, set_hkl_int)
        except AssertionError as error:
            self.pop_one_button_dialog(str(error))
            return

        # show result
        self._show_refined_ub_result()

        return

    def do_refine_ub_fft(self):
        """
        Refine UB matrix by calling FFT method
        :return:
        """
        # get PeakInfo list and check
        peak_info_list = self._build_peak_info_list()
        assert isinstance(peak_info_list, list), \
            'PeakInfo list must be a list but not %s.' % str(type(peak_info_list))
        assert len(peak_info_list) >= 3, \
            'PeakInfo must be larger or equal to 3 (.now given %d) to refine UB matrix' % len(peak_info_list)

        # get lattice range information
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_minD,
                                                        self.ui.lineEdit_maxD])
        if status is False:
            self.pop_one_button_dialog('Must specify Min D and max D to refine UB using FFT.')
            return
        min_d, max_d = ret_obj
        if (0 < min_d < max_d) is False:
            self.pop_one_button_dialog('Range of d is not correct!')
            return

        # friendly suggestion
        if len(peak_info_list) <= 9:
            self.pop_one_button_dialog('It is recommended to use at least 9 reflections'
                                       'to refine UB matrix without prior knowledge.')

        # refine
        self._myControl.refine_ub_matrix_least_info(peak_info_list, min_d, max_d)

        # set value
        self._show_refined_ub_result()

        return

    def _show_refined_ub_result(self):
        """
        Show the result from refined UB matrix
        :return:
        """
        # Deal with result
        ub_matrix, lattice, lattice_error = self._myControl.get_refined_ub_matrix()
        # ub matrix
        self.ui.tableWidget_ubMatrix.set_from_matrix(ub_matrix)

        # lattice parameter
        assert isinstance(lattice, list)
        assert len(lattice) == 6
        self.ui.lineEdit_aUnitCell.setText('%.5f' % lattice[0])
        self.ui.lineEdit_bUnitCell.setText('%.5f' % lattice[1])
        self.ui.lineEdit_cUnitCell.setText('%.5f' % lattice[2])
        self.ui.lineEdit_alphaUnitCell.setText('%.5f' % lattice[3])
        self.ui.lineEdit_betaUnitCell.setText('%.5f' % lattice[4])
        self.ui.lineEdit_gammaUnitCell.setText('%.5f' % lattice[5])

        assert isinstance(lattice_error, list)
        assert len(lattice_error) == 6
        self.ui.lineEdit_aError.setText('%.5f' % lattice_error[0])
        self.ui.lineEdit_bError.setText('%.5f' % lattice_error[1])
        self.ui.lineEdit_cError.setText('%.5f' % lattice_error[2])
        self.ui.lineEdit_alphaError.setText('%.5f' % lattice_error[3])
        self.ui.lineEdit_betaError.setText('%.5f' % lattice_error[4])
        self.ui.lineEdit_gammaError.setText('%.5f' % lattice_error[5])

        return

    def _build_peak_info_list(self):
        """ Build a list of PeakInfo to build peak workspace
        :return:
        """
        # Collecting all peaks that will be used to refine UB matrix
        row_index_list = self.ui.tableWidget_peaksCalUB.get_selected_rows(True)
        if len(row_index_list) < 3:
            err_msg = 'At least 3 peaks must be selected to refine UB matrix.' \
                      'Now it is only %d selected.' % len(row_index_list)
            self.pop_one_button_dialog(err_msg)
            return

        # loop over all peaks for peak information
        peak_info_list = list()
        status, exp_number = gutil.parse_integers_editors(self.ui.lineEdit_exp)
        assert status
        for i_row in row_index_list:
            scan_num, pt_num = self.ui.tableWidget_peaksCalUB.get_exp_info(i_row)
            try:
                if pt_num < 0:
                    pt_num = None
                peak_info = self._myControl.get_peak_info(exp_number, scan_num, pt_num)
            except AssertionError as ass_err:
                raise RuntimeError('Unable to retrieve PeakInfo due to %s.' % str(ass_err))
            assert isinstance(peak_info, r4c.PeakInfo)
            peak_info_list.append(peak_info)
        # END-FOR

        return peak_info_list

    def change_data_access_mode(self):
        """ Change data access mode between downloading from server and local
        Event handling methods
        :return:
        """
        new_mode = str(self.ui.comboBox_mode.currentText())
        self._dataAccessMode = new_mode

        if new_mode.startswith('Local') is True:
            self.ui.lineEdit_localSpiceDir.setEnabled(True)
            self.ui.pushButton_browseLocalDataDir.setEnabled(True)
            self.ui.lineEdit_url.setEnabled(False)
            self.ui.lineEdit_localSrcDir.setEnabled(False)
            self.ui.pushButton_browseLocalCache.setEnabled(False)
            self._allowDownload = False
        else:
            self.ui.lineEdit_localSpiceDir.setEnabled(False)
            self.ui.pushButton_browseLocalDataDir.setEnabled(False)
            self.ui.lineEdit_url.setEnabled(True)
            self.ui.lineEdit_localSrcDir.setEnabled(True)
            self.ui.pushButton_browseLocalCache.setEnabled(True)
            self._allowDownload = True

        return

    def change_instrument_name(self):
        """ Handing the event as the instrument name is changed
        :return:
        """
        new_instrument = str(self.ui.comboBox_instrument.currentText())
        self.pop_one_button_dialog('Change of instrument during data processing is dangerous.')
        status, error_message = self._myControl.set_instrument_name(new_instrument)
        if status is False:
            self.pop_one_button_dialog(error_message)

        return

    def do_add_peak_no_index(self):
        """
        Purpose: add a peak from 'View Raw Data' tab to the UB peak table
        without indexing it
        :return:
        """
        # Get exp, scan and Pt information
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_exp,
                                                        self.ui.lineEdit_run,
                                                        self.ui.lineEdit_rawDataPtNo])
        if not status:
            err_msg = ret_obj
            self.pop_one_button_dialog(err_msg)
            return
        else:
            int_list = ret_obj
            exp_no, scan_no, pt_no = int_list

        # Switch tab
        self.ui.tabWidget.setCurrentIndex(MainWindow.TabPage['Calculate UB'])

        # Find and index peak
        self._myControl.find_peak(exp_no, scan_no, [pt_no])
        try:
            peak_info = self._myControl.get_peak_info(exp_no, scan_no, pt_no)
            self.set_ub_peak_table(peak_info)
        except AssertionError as ass_err:
            self.pop_one_button_dialog(str(ass_err))
            return

        return

    def do_add_ub_peak(self):
        """ Add current to ub peaks
        :return:
        """
        # Add peak
        status, int_list = gutil.parse_integers_editors([self.ui.lineEdit_exp,
                                                         self.ui.lineEdit_scanNumber])
        if status is False:
            self.pop_one_button_dialog(int_list)
            return
        exp_no, scan_no = int_list

        # Get HKL from GUI
        status, float_list = gutil.parse_float_editors([self.ui.lineEdit_H,
                                                        self.ui.lineEdit_K,
                                                        self.ui.lineEdit_L])
        if status is False:
            err_msg = float_list
            self.pop_one_button_dialog(err_msg)
            return
        h, k, l = float_list

        try:
            peak_info_obj = self._myControl.get_peak_info(exp_no, scan_no)
        except AssertionError as ass_err:
            self.pop_one_button_dialog(str(ass_err))
            return

        assert isinstance(peak_info_obj, r4c.PeakInfo)
        if self.ui.checkBox_roundHKLInt.isChecked():
            h = math.copysign(1, h)*int(abs(h)+0.5)
            k = math.copysign(1, k)*int(abs(k)+0.5)
            l = math.copysign(1, l)*int(abs(l)+0.5)
        peak_info_obj.set_user_hkl(h, k, l)
        self.set_ub_peak_table(peak_info_obj)

        # Clear
        self.ui.lineEdit_scanNumber.setText('')

        self.ui.lineEdit_sampleQx.setText('')
        self.ui.lineEdit_sampleQy.setText('')
        self.ui.lineEdit_sampleQz.setText('')

        self.ui.lineEdit_H.setText('')
        self.ui.lineEdit_K.setText('')
        self.ui.lineEdit_L.setText('')

        return

    def doAcceptCalUB(self):
        """ Accept the calculated UB matrix
        """
        raise RuntimeError('ASAP')

    def doAddScanPtToRefineUB(self):
        """ Add scan/pt numbers to the list of data points for refining ub matrix

        And the added scan number and pt numbers will be reflected in the (left sidebar)

        """
        raise RuntimeError("ASAP")

    def do_add_peak_to_find(self):
        """
        Add the scan/pt to the next
        :return:
        """
        # peak finding will use all points in the selected scan.
        scan_no = self.ui.lineEdit_run.text()
        self.ui.tabWidget.setCurrentIndex(MainWindow.TabPage['Calculate UB'])
        self.ui.lineEdit_scanNumber.setText(scan_no)

        return

    def do_add_peaks_for_ub(self):
        """ In tab-survey, merge selected scans, find peaks in merged data and
         switch to UB matrix calculation tab and add to table
        :return:
        """
        # get selected scans
        selected_row_index_list = self.ui.tableWidget_surveyTable.get_selected_rows(True)
        scan_number_list = self.ui.tableWidget_surveyTable.get_scan_numbers(selected_row_index_list)
        if len(scan_number_list) == 0:
            self.pop_one_button_dialog('No scan is selected.')
            return

        # get experiment number
        status, exp_number = gutil.parse_integers_editors(self.ui.lineEdit_exp)
        assert status

        # switch to tab-3
        self.ui.tabWidget.setCurrentIndex(MainWindow.TabPage['Calculate UB'])

        # find peak and add peak
        failed_list = list()
        for scan_number in scan_number_list:
            # merge peak
            status, err_msg = self._myControl.merge_pts_in_scan(exp_number, scan_number, [], 'q-sample')

            # continue to the next scan if there is something wrong
            if status is False:
                failed_list.append((scan_number, err_msg))
                continue

            # find peak
            self._myControl.find_peak(exp_number, scan_number)

            # get PeakInfo
            peak_info = self._myControl.get_peak_info(exp_number, scan_number)
            assert isinstance(peak_info, r4c.PeakInfo)

            # retrieve and set HKL from spice table
            peak_info.retrieve_hkl_from_spice_table()

            # add to table
            self.set_ub_peak_table(peak_info)
        # END-FOR

        # pop error if there is any scan that is not reduced right
        if len(failed_list) > 0:
            failed_scans_str = 'Unable to merge scans: '
            sum_error_str = ''
            for fail_tup in failed_list:
                failed_scans_str += '%d, ' % fail_tup[0]
                sum_error_str += '%s\n' % fail_tup[1]
            # END-FOR

            self.pop_one_button_dialog(failed_scans_str)
            self.pop_one_button_dialog(sum_error_str)
        # END-FOR

        return

    def do_browse_local_cache_dir(self):
        """ Browse local cache directory
        :return:
        """
        local_cache_dir = str(QtGui.QFileDialog.getExistingDirectory(self,
                                                                     'Get Local Cache Directory',
                                                                     self._homeSrcDir))

        # Set local directory to control
        status, error_message = self._myControl.set_local_data_dir(local_cache_dir)
        if status is False:
            self.pop_one_button_dialog(error_message)
            return

        # Synchronize to local data/spice directory and local cache directory
        if str(self.ui.lineEdit_localSpiceDir.text()) != '':
            prev_dir = str(self.ui.lineEdit_localSrcDir.text())
            self.pop_one_button_dialog('Local data directory was set up as %s' %
                                       prev_dir)
        self.ui.lineEdit_localSrcDir.setText(local_cache_dir)
        self.ui.lineEdit_localSpiceDir.setText(local_cache_dir)

        return

    def do_browse_local_spice_data(self):
        """ Browse local source SPICE data directory
        """
        src_spice_dir = str(QtGui.QFileDialog.getExistingDirectory(self, 'Get Directory',
                                                                   self._homeSrcDir))
        # Set local data directory to controller
        status, error_message = self._myControl.set_local_data_dir(src_spice_dir)
        if status is False:
            self.pop_one_button_dialog(error_message)
            return

        self._homeSrcDir = src_spice_dir
        self.ui.lineEdit_localSpiceDir.setText(src_spice_dir)

        return

    def do_browse_working_dir(self):
        """
        Browse and set up working directory
        :return:
        """
        work_dir = str(QtGui.QFileDialog.getExistingDirectory(self, 'Get Working Directory', self._homeDir))
        status, error_message = self._myControl.set_working_directory(work_dir)
        if status is False:
            self.pop_one_button_dialog(error_message)
        else:
            self.ui.lineEdit_workDir.setText(work_dir)

        return

    def do_cal_ub_matrix(self):
        """ Calculate UB matrix by 2 or 3 reflections
        """
        # Get reflections selected to calculate UB matrix
        num_rows = self.ui.tableWidget_peaksCalUB.rowCount()
        peak_info_list = list()
        status, exp_number = gutil.parse_integers_editors(self.ui.lineEdit_exp)
        assert status
        for i_row in xrange(num_rows):
            if self.ui.tableWidget_peaksCalUB.is_selected(i_row) is True:
                scan_num, pt_num = self.ui.tableWidget_peaksCalUB.get_exp_info(i_row)
                if pt_num < 0:
                    pt_num = None
                peak_info = self._myControl.get_peak_info(exp_number, scan_num, pt_num)
                assert isinstance(peak_info, r4c.PeakInfo)
                peak_info_list.append(peak_info)
        # END-FOR

        # Get lattice
        status, ret_obj = self._get_lattice_parameters()
        if status is True:
            a, b, c, alpha, beta, gamma = ret_obj
        else:
            err_msg = ret_obj
            self.pop_one_button_dialog(err_msg)
            return

        # Calculate UB matrix
        status, ub_matrix = self._myControl.calculate_ub_matrix(peak_info_list, a, b, c,
                                                                alpha, beta, gamma)

        # Deal with result
        if status is True:
            self.ui.tableWidget_ubMatrix.set_from_matrix(ub_matrix)

        else:
            err_msg = ub_matrix
            self.pop_one_button_dialog(err_msg)

        return

    def do_clear_survey(self):
        """
        Clear survey and survey table.
        As myController does not store any survey information,
        there is no need to clear any data structure in myController
        :return:
        """
        # Clear table
        self.ui.tableWidget_surveyTable.remove_all_rows()
        self.ui.tableWidget_surveyTable.reset()

        return

    def do_clear_ub_peaks(self):
        """
        Clear all peaks in UB-Peak table
        :return:
        """
        num_rows = self.ui.tableWidget_peaksCalUB.rowCount()
        row_number_list = range(num_rows)
        self.ui.tableWidget_peaksCalUB.delete_rows(row_number_list)

        return

    def do_del_ub_peaks(self):
        """
        Delete a peak in UB-Peak table
        :return:
        """
        # Find out the lines to get deleted
        row_num_list = self.ui.tableWidget_peaksCalUB.get_selected_rows()
        print '[DB] Row %s are selected' % str(row_num_list)

        # Delete
        self.ui.tableWidget_peaksCalUB.delete_rows(row_num_list)

        return

    def do_download_spice_data(self):
        """ Download SPICE data
        :return:
        """
        # Check scans to download
        scan_list_str = str(self.ui.lineEdit_downloadScans.text())
        if len(scan_list_str) > 0:
            # user specifies scans to download
            valid, scan_list = fcutil.parse_int_array(scan_list_str)
            if valid is False:
                error_message = scan_list
                self.pop_one_button_dialog(error_message)
        else:
            # Get all scans
            status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_exp])
            if status is False:
                self.pop_one_button_dialog(ret_obj)
                return
            exp_no = ret_obj
            assert isinstance(exp_no, int)
            server_url = str(self.ui.lineEdit_url.text())
            scan_list = fcutil.get_scans_list(server_url, exp_no, return_list=True)
        self.pop_one_button_dialog('Going to download scans %s.' % str(scan_list))

        # Check location
        destination_dir = str(self.ui.lineEdit_localSrcDir.text())
        status, error_message = self._myControl.set_local_data_dir(destination_dir)
        if status is False:
            self.pop_one_button_dialog(error_message)
        else:
            self.pop_one_button_dialog('Spice files will be downloaded to %s.' % destination_dir)

        # Set up myControl for downloading data
        exp_no = int(self.ui.lineEdit_exp.text())
        self._myControl.set_exp_number(exp_no)

        server_url = str(self.ui.lineEdit_url.text())
        status, error_message = self._myControl.set_server_url(server_url)
        if status is False:
            self.pop_one_button_dialog(error_message)
            return

        # Download
        self._myControl.download_data_set(scan_list)

        return

    def do_find_peak(self):
        """ Find peak in a given scan and record it
        """
        # Get experiment, scan and pt
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_exp,
                                                        self.ui.lineEdit_scanNumber])
        if status is True:
            exp_no, scan_no = ret_obj
        else:
            self.pop_one_button_dialog(ret_obj)
            return

        # merge peak if necessary
        if self._myControl.has_merged_data(exp_no, scan_no) is False:
            status, err_msg = self._myControl.merge_pts_in_scan(exp_no, scan_no, [], 'q-sample')
            if status is False:
                self.pop_one_button_dialog(err_msg)

        # Find peak
        status, err_msg = self._myControl.find_peak(exp_no, scan_no)
        if status is False:
            self.pop_one_button_dialog(ret_obj)
            return

        # Get information from the latest (integrated) peak
        if self.ui.checkBox_loadHKLfromFile.isChecked() is True:
            # This is the first time that in the workflow to get HKL from MD workspace
            peak_info = self._myControl.get_peak_info(exp_no, scan_no)
            try:
                peak_info.retrieve_hkl_from_spice_table()
            except RuntimeError as run_err:
                self.pop_one_button_dialog('Unable to locate peak info due to %s.' % str(run_err))
        # END-IF

        # Set up correct values to table tableWidget_peaksCalUB
        peak_info = self._myControl.get_peak_info(exp_no, scan_no)
        h, k, l = peak_info.get_user_hkl()
        self.ui.lineEdit_H.setText('%.2f' % h)
        self.ui.lineEdit_K.setText('%.2f' % k)
        self.ui.lineEdit_L.setText('%.2f' % l)

        q_x, q_y, q_z = peak_info.get_peak_centre()
        self.ui.lineEdit_sampleQx.setText('%.5E' % q_x)
        self.ui.lineEdit_sampleQy.setText('%.5E' % q_y)
        self.ui.lineEdit_sampleQz.setText('%.5E' % q_z)

        return

    def do_index_ub_peaks(self):
        """ Index the peaks in the UB matrix peak table
        :return:
        """
        # Get UB matrix
        ub_matrix = self.ui.tableWidget_ubMatrix.get_matrix()
        print '[Info] Get UB matrix from table ', ub_matrix

        # Index all peaks
        num_peaks = self.ui.tableWidget_peaksCalUB.rowCount()
        err_msg = ''
        for i_peak in xrange(num_peaks):
            scan_no = self.ui.tableWidget_peaksCalUB.get_exp_info(i_peak)[0]
            status, ret_obj = self._myControl.index_peak(ub_matrix, scan_number=scan_no)
            if status is True:
                hkl_value = ret_obj[0]
                hkl_error = ret_obj[1]
                self.ui.tableWidget_peaksCalUB.set_hkl(i_peak, hkl_value, hkl_error)
            else:
                err_msg += ret_obj + '\n'
        # END-FOR

        if len(err_msg) > 0:
            self.pop_one_button_dialog(err_msg)

        return

    def do_list_scans(self):
        """ List all scans available
        :return:
        """
        # Experiment number
        exp_no = int(self.ui.lineEdit_exp.text())

        access_mode = str(self.ui.comboBox_mode.currentText())
        if access_mode == 'Local':
            spice_dir = str(self.ui.lineEdit_localSpiceDir.text())
            message = fcutil.get_scans_list_local_disk(spice_dir, exp_no)
        else:
            url = str(self.ui.lineEdit_url.text())
            message = fcutil.get_scans_list(url, exp_no)

        self.pop_one_button_dialog(message)

        return

    def do_load_scan_info(self):
        """ Load SIICE's scan file
        :return:
        """
        # Get scan number
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_run])
        if status is True:
            scan_no = ret_obj[0]
        else:
            err_msg = ret_obj
            self.pop_one_button_dialog('Unable to get scan number in raw data tab due to %s.' % err_msg)
            return

        status, err_msg = self._myControl.load_spice_scan_file(exp_no=None, scan_no=scan_no)
        if status is False:
            self.pop_one_button_dialog(err_msg)

        return

    def do_load_survey(self):
        """ Load csv file containing experiment-scan survey's result.
        :return:
        """
        # check validity
        num_rows = int(self.ui.lineEdit_numSurveyOutput.text())

        # get the csv file
        file_filter = 'CSV Files (*.csv);;All Files (*.*)'
        csv_file_name = str(QtGui.QFileDialog.getOpenFileName(self, 'Open Exp-Scan Survey File', self._homeDir,
                                                              file_filter))
        if csv_file_name is None or len(csv_file_name) == 0:
            # return if file selection is cancelled
            return

        # call controller to load
        survey_tuple = self._myControl.load_scan_survey_file(csv_file_name)
        scan_sum_list = survey_tuple[1]
        assert isinstance(scan_sum_list, list), 'Returned value from load scan survey file must be a dictionary.'

        # set the table
        self.ui.tableWidget_surveyTable.set_survey_result(scan_sum_list)
        self.ui.tableWidget_surveyTable.remove_all_rows()
        self.ui.tableWidget_surveyTable.show_reflections(num_rows)

        return

    def do_plot_pt_raw(self):
        """ Plot the Pt.
        """
        # Get measurement pt and the file number
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_exp,
                                                        self.ui.lineEdit_run,
                                                        self.ui.lineEdit_rawDataPtNo])
        if status is True:
            exp_no = ret_obj[0]
            scan_no = ret_obj[1]
            pt_no = ret_obj[2]
        else:
            self.pop_one_button_dialog(ret_obj)
            return

        # Call to plot 2D
        self._plot_raw_xml_2d(exp_no, scan_no, pt_no)

        return

    def do_plot_prev_pt_raw(self):
        """ Plot the Pt.
        """
        # Get measurement pt and the file number
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_exp,
                                                        self.ui.lineEdit_run,
                                                        self.ui.lineEdit_rawDataPtNo])
        if status is True:
            exp_no = ret_obj[0]
            scan_no = ret_obj[1]
            pt_no = ret_obj[2]
        else:
            self.pop_one_button_dialog(ret_obj)
            return

        # Previous one
        pt_no -= 1
        if pt_no <= 0:
            self.pop_one_button_dialog('Pt. = 1 is the first one.')
            return
        else:
            self.ui.lineEdit_rawDataPtNo.setText('%d' % pt_no)

        # Plot
        self._plot_raw_xml_2d(exp_no, scan_no, pt_no)

        return

    def do_plot_next_pt_raw(self):
        """ Plot the Pt.
        """
        # Get measurement pt and the file number
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_exp,
                                                        self.ui.lineEdit_run,
                                                        self.ui.lineEdit_rawDataPtNo])
        if status is True:
            exp_no = ret_obj[0]
            scan_no = ret_obj[1]
            pt_no = ret_obj[2]
        else:
            self.pop_one_button_dialog(ret_obj)
            return

        # Previous one
        pt_no += 1
        # get last Pt. number
        status, last_pt_no = self._myControl.get_pt_numbers(exp_no, scan_no)
        if status is False:
            error_message = last_pt_no
            self.pop_one_button_dialog('Unable to access Spice table for scan %d. Reason" %s.' % (
                scan_no, error_message))
        if pt_no > last_pt_no:
            self.pop_one_button_dialog('Pt. = %d is the last one of scan %d.' % (pt_no, scan_no))
            return
        else:
            self.ui.lineEdit_rawDataPtNo.setText('%d' % pt_no)

        # Plot
        self._plot_raw_xml_2d(exp_no, scan_no, pt_no)

        return

    def do_add_scans_merge(self):
        """ Add scans to merge
        :return:
        """
        # Get list of scans
        scan_list = gutil.parse_integer_list(str(self.ui.lineEdit_listScansSliceView.text()))
        if len(scan_list) == 0:
            self.pop_one_button_dialog('Scan list is empty.')

        # Set table
        self.ui.tableWidget_mergeScans.append_scans(scans=scan_list)

        return

    def do_merge_scans(self):
        """ Process data for slicing view
        :return:
        """
        # Get UB matrix
        ub_matrix = self.ui.tableWidget_ubInUse.get_matrix()
        self._myControl.set_ub_matrix(exp_number=None, ub_matrix=ub_matrix)

        # Warning
        self.pop_one_button_dialog('Data processing is long. Be patient!')

        # Process
        scan_row_list = self.ui.tableWidget_mergeScans.get_scan_list()
        print '[DB] %d scans have been selected to merge.' % len(scan_row_list)
        frame = str(self.ui.comboBox_mergeScanFrame.currentText())
        for tup2 in scan_row_list:
            #
            scan_no, i_row = tup2

            # Download/check SPICE file
            self._myControl.download_spice_file(None, scan_no, over_write=False)

            # Get some information
            status, pt_list = self._myControl.get_pt_numbers(None, scan_no)
            if status is False:
                err_msg = pt_list
                self.pop_one_button_dialog('Failed to get Pt. number: %s' % err_msg)
                return
            else:
                # Set information to table
                err_msg = self.ui.tableWidget_mergeScans.set_scan_pt(scan_no, pt_list)
                if len(err_msg) > 0:
                    self.pop_one_button_dialog(err_msg)

            self.ui.tableWidget_mergeScans.set_status_by_row(i_row, 'In Processing')
            merge_status = 'UNKNOWN'
            merged_name = '???'
            group_name = '???'

            status, ret_tup = self._myControl.merge_pts_in_scan(exp_no=None, scan_no=scan_no,
                                                                pt_num_list=[], target_frame=frame)
            merge_status = 'Done'
            merged_name = ret_tup[0]
            group_name = ret_tup[1]

            if status is False:
                merge_status = 'Failed. Reason: %s' % str(e)
                merged_name = ''
                group_name = ''
                print merge_status
            else:
                self.ui.tableWidget_mergeScans.set_status_by_row(i_row, merge_status)
                self.ui.tableWidget_mergeScans.set_ws_names_by_row(i_row, merged_name, group_name)

            # Sleep for a while
            time.sleep(0.1)
        # END-FOR

        return

    def do_refresh_merged_scans_table(self):
        """ Find the merged
        :return:
        """
        # find out the merged runs
        scan_info_tup_list = self._myControl.get_merged_scans()
        print '[DB-BAT] Scan Info List: ', scan_info_tup_list

        # append the row to the merged scan table
        for scan_info_tup in scan_info_tup_list:
            exp_number = scan_info_tup[0]
            scan_number = scan_info_tup[1]
            # pt_number_list = scan_info_tup[2]
            ws_name = 'whatever'
            ws_group = -1
            self.ui.tableWidget_mergeScans.add_new_merged_data(exp_number, scan_number, ws_name, ws_group)

        return

    def do_reset_ub_peaks_hkl(self):
        """
        Reset user specified HKL value to peak table
        :return:
        """
        # get experiment number
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_exp])
        assert status, ret_obj
        exp_number = ret_obj[0]

        # reset all rows back to SPICE HKL
        num_rows = self.ui.tableWidget_peaksCalUB.rowCount()
        for i_row in xrange(num_rows):
            scan, pt = self.ui.tableWidget_peaksCalUB.get_scan_pt(i_row)
            if pt < 0:
                pt = None
            peak_info = self._myControl.get_peak_info(exp_number, scan, pt)
            h, k, l = peak_info.get_user_hkl()
            self.ui.tableWidget_peaksCalUB.update_hkl(i_row, h, k, l)
        # END-FOR

        return

    def do_save_survey(self):
        """
        Save the survey to a file
        :return:
        """
        # Get file name
        file_filter = 'CSV Files (*.csv);;All Files (*.*)'
        out_file_name = str(QtGui.QFileDialog.getSaveFileName(self, 'Save scan survey result',
                                                              self._homeDir, file_filter))

        # Save file
        self._myControl.save_scan_survey(out_file_name)

        return

    def do_select_all_peaks(self):
        """
        Purpose: select all peaks in table tableWidget_peaksCalUB
        :return:
        """
        self.ui.tableWidget_peaksCalUB.select_all_rows(self._ubPeakTableFlag)
        self._ubPeakTableFlag = not self._ubPeakTableFlag

        return

    def do_select_all_survey(self):
        """
        Select or de-select all rows in survey items
        :return:
        """
        self.ui.tableWidget_surveyTable.select_all_rows(self._surveyTableFlag)
        self._surveyTableFlag = not self._surveyTableFlag

        return

    def do_set_experiment(self):
        """ Set experiment
        :return:
        """
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_exp])
        if status is True:
            exp_number = ret_obj[0]
            curr_exp_number = self._myControl.get_experiment()
            if curr_exp_number is not None and exp_number != curr_exp_number:
                self.pop_one_button_dialog('Changing experiment to %d.  Clean previous experiment %d\'s result'
                                           ' in Mantid manually.' % (exp_number, curr_exp_number))
            self._myControl.set_exp_number(exp_number)
            self.ui.lineEdit_exp.setStyleSheet('color: black')
        else:
            err_msg = ret_obj
            self.pop_one_button_dialog('Unable to set experiment as %s' % err_msg)
            self.ui.lineEdit_exp.setStyleSheet('color: red')

        self.ui.tabWidget.setCurrentIndex(0)

        return

    def do_setup_dir_default(self):
        """
        Set up default directory for storing data and working
        :return:
        """
        home_dir = os.path.expanduser('~')

        # Data cache directory
        data_cache_dir = os.path.join(home_dir, 'Temp/HB3ATest')
        self.ui.lineEdit_localSpiceDir.setText(data_cache_dir)
        self.ui.lineEdit_localSrcDir.setText(data_cache_dir)

        work_dir = os.path.join(data_cache_dir, 'Workspace')
        self.ui.lineEdit_workDir.setText(work_dir)

        return

    def do_apply_setup(self):
        """
        Purpose:
         - Apply the setup to controller.
        Requirements:
         - data directory, working directory must be given; but not necessarily correct
         - URL must be given; but not necessary to be correct
        :return:
        """
        # get data directory, working directory and data server URL from GUI
        local_data_dir = str(self.ui.lineEdit_localSpiceDir.text()).strip()
        working_dir = str(self.ui.lineEdit_workDir.text()).strip()
        data_server = str(self.ui.lineEdit_url.text()).strip()

        # set to my controller
        self._myControl.set_local_data_dir(local_data_dir)
        self._myControl.set_working_directory(working_dir)
        self._myControl.set_server_url(data_server, check_link=False)

        # check
        error_message = ''

        # local data dir
        if local_data_dir == '':
            error_message += 'Local data directory is not specified!\n'
        elif os.path.exists(local_data_dir) is False:
            try:
                os.mkdir(local_data_dir)
            except OSError as os_error:
                error_message += 'Unable to create local data directory %s due to %s.\n' % (
                    local_data_dir, str(os_error))
                self.ui.lineEdit_localSpiceDir.setStyleSheet("color: red;")
            else:
                self.ui.lineEdit_localSpiceDir.setStyleSheet("color: green;")
        else:
            self.ui.lineEdit_localSpiceDir.setStyleSheet("color: green;")
        # END-IF-ELSE

        # working directory
        if working_dir == '':
            error_message += 'Working directory is not specified!\n'
        elif os.path.exists(working_dir) is False:
            try:
                os.mkdir(working_dir)
                self.ui.lineEdit_workDir.setStyleSheet("color: green;")
            except OSError as os_error:
                error_message += 'Unable to create working directory %s due to %s.\n' % (
                    working_dir, str(os_error))
                self.ui.lineEdit_workDir.setStyleSheet("color: red;")
        else:
            self.ui.lineEdit_workDir.setStyleSheet("color: green;")
        # END-IF-ELSE

        # Set the URL red as it is better not check at this stage. Leave it to user
        self.ui.lineEdit_url.setStyleSheet("color: black;")

        if len(error_message) > 0:
            self.pop_one_button_dialog(error_message)

        return

    def do_filter_sort_survey_table(self):
        """
        Sort and filter survey table by specified field
        Requirements:
        Guarantees: the original table is cleared and a new list is appended
        :return:
        """
        # Get column name
        if self.ui.radioButton_sortByScan.isChecked():
            column_name = 'Scan'
        elif self.ui.radioButton_sortByCounts.isChecked():
            column_name = 'Max Counts'
        else:
            self.pop_one_button_dialog('No column is selected to sort.')
            return

        # Get filters
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_filterScanLower,
                                                        self.ui.lineEdit_filterScanUpper],
                                                       allow_blank=True)

        # return with error
        if status is False:
            self.pop_one_button_dialog(ret_obj)
            return

        # set up default with return as None
        start_scan_number = ret_obj[0]
        if start_scan_number is None:
            start_scan_number = 0
        end_scan_number = ret_obj[1]
        if end_scan_number is None:
            end_scan_number = sys.maxint

        status, ret_obj = gutil.parse_float_editors([self.ui.lineEdit_filterCountsLower,
                                                     self.ui.lineEdit_filterCountsUpper],
                                                    allow_blank=True)
        if status is False:
            # return with error message
            self.pop_one_button_dialog(ret_obj)
            return

        # set up default with return as None
        min_counts = ret_obj[0]
        if min_counts is None:
            min_counts = -0.0
        max_counts = ret_obj[1]
        if max_counts is None:
            max_counts = sys.float_info.max

        # filter and sort
        ascending_order = not self.ui.checkBox_sortDescending.isChecked()
        if ascending_order:
            sort_order = 0
        else:
            sort_order = 1
        self.ui.tableWidget_surveyTable.filter_and_sort(start_scan_number, end_scan_number,
                                                        min_counts, max_counts,
                                                        column_name, sort_order)

        return

    def do_set_ub_from_text(self):
        """ Purpose: Set UB matrix in use from plain text edit plainTextEdit_ubInput.
        Requirements:
          1. the string in the plain text edit must be able to be split to 9 floats by ',', ' ', '\t' and '\n'
        Guarantees: the matrix will be set up the UB matrix in use
        :return:
        """
        ub_str = str(self.ui.plainTextEdit_ubInput.toPlainText())
        status, ret_obj = gutil.parse_float_array(ub_str)
        if status is False:
            # unable to parse to float arrays
            self.pop_one_button_dialog(ret_obj)
            return

        elif len(ret_obj) != 9:
            # number of floats is not 9
            self.pop_one_button_dialog('Requiring 9 floats for UB matrix.  Only %d are given.' % len(ret_obj))
            return

        # in good UB matrix format
        ub_str = ret_obj
        if self.ui.radioButton_ubMantidStyle.isChecked():
            # UB matrix in mantid style
            self.ui.tableWidget_ubInUse.set_from_list(ub_str)

        elif self.ui.radioButton_ubSpiceStyle.isChecked():
            # UB matrix in SPICE style
            spice_ub = gutil.convert_str_to_matrix(ub_str, (3, 3))
            mantid_ub = r4c.convert_spice_ub_to_mantid(spice_ub)
            self.ui.tableWidget_ubInUse.set_from_matrix(mantid_ub)

        else:
            # not defined
            self.pop_one_button_dialog('Neither Mantid or SPICE-styled UB is checked!')

        return

    def do_show_ub_in_box(self):
        """ Get UB matrix in table tableWidget_ubMergeScan and write to plain text edit plainTextEdit_ubInput
        :return:
        """
        ub_matrix = self.ui.tableWidget_ubInUse.get_matrix()

        text = ''
        for i in xrange(3):
            for j in xrange(3):
                text += '%.7f, ' % ub_matrix[i][j]
            text += '\n'

        self.ui.plainTextEdit_ubInput.setPlainText(text)

        return

    def do_survey(self):
        """
        Purpose: survey for the strongest reflections
        :return:
        """
        # Get experiment number
        exp_number = int(self.ui.lineEdit_exp.text())
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_surveyStartPt,
                                                        self.ui.lineEdit_surveyEndPt])
        if status is False:
            err_msg = ret_obj
            self.pop_one_button_dialog(err_msg)
        start_scan = ret_obj[0]
        end_scan = ret_obj[1]

        max_number = int(self.ui.lineEdit_numSurveyOutput.text())

        # Get value
        status, ret_obj = self._myControl.survey(exp_number, start_scan, end_scan)
        if status is False:
            self.pop_one_button_dialog(ret_obj)
            return
        scan_sum_list = ret_obj
        self.ui.tableWidget_surveyTable.set_survey_result(scan_sum_list)
        self.ui.tableWidget_surveyTable.show_reflections(max_number)

        return

    def do_switch_tab_peak_int(self):
        """ Switch to tab 'Peak Integration' to set up and learn how to do peak integration
        :return:
        """
        # switch tab
        self.ui.tabWidget.setCurrentIndex(MainWindow.TabPage['Peak Integration'])

        # set up value
        selected_scans = self.ui.tableWidget_mergeScans.get_scan_list()
        if len(selected_scans) > 0:
            self.ui.lineEdit_scanIntegratePeak.setText(str(selected_scans[0]))

        return

    def do_sync_ub(self):
        """ Purpose: synchronize UB matrix in use with UB matrix calculated.
        :return:
        """
        self.ui.tableWidget_ubInUse.set_from_matrix(self.ui.tableWidget_ubMatrix.get_matrix())

        return

    def do_test_url(self):
        """ Test whether the root URL provided specified is good
        """
        url = str(self.ui.lineEdit_url.text())

        url_is_good, err_msg = fcutil.check_url(url)
        if url_is_good is True:
            self.pop_one_button_dialog("URL %s is valid." % url)
            self.ui.lineEdit_url.setStyleSheet("color: green;")
        else:
            self.pop_one_button_dialog(err_msg)
            self.ui.lineEdit_url.setStyleSheet("color: read;")

        return url_is_good

    def do_view_data_set_3d(self):
        """
        Launch the sub window to view merged data in 3D.
        :return:
        """
        # get a list of rows that are selected
        row_index_list = self.ui.tableWidget_peaksCalUB.get_selected_rows(True)
        assert len(row_index_list) > 0, 'There is no row that is selected to view.'

        # get the information from the selected row
        status, exp_number = gutil.parse_integers_editors(self.ui.lineEdit_exp)
        assert status, 'Experiment number is not a valid integer.'

        # create window
        if self._my3DWindow is None:
            self._my3DWindow = plot3dwindow.Plot3DWindow(self)

        for i_row in row_index_list:
            exp_info = self.ui.tableWidget_peaksCalUB.get_exp_info(i_row)
            scan_number = exp_info[0]

            # prepare data
            ret_obj = self._prepare_view_merged(exp_number, scan_number)
            assert len(ret_obj) == 5
            md_file_name, weight_peak_centers, weight_peak_intensities, avg_peak_centre, avg_peak_intensity = ret_obj

            # add the 3D window
            self._my3DWindow.initialize_group('%s' % scan_number)
            self._my3DWindow.add_plot_by_file(md_file_name)
            self._my3DWindow.add_plot_by_array(weight_peak_centers, weight_peak_intensities)
            self._my3DWindow.add_plot_by_array(avg_peak_centre, avg_peak_intensity)
            self._my3DWindow.close_session()

        # END-FOR

        # Show
        self._my3DWindow.show()

        return

    def do_view_data_3d(self):
        """
        View merged scan data in 3D after FindPeaks
        :return:
        """
        # get experiment and scan number
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_exp,
                                                        self.ui.lineEdit_scanNumber])
        if status:
            exp_number = ret_obj[0]
            scan_number = ret_obj[1]
        else:
            self.pop_one_button_dialog(ret_obj)
            return

        # Check
        if self._myControl.has_merged_data(exp_number, scan_number) is False:
            self.pop_one_button_dialog('No merged MD workspace found for %d, %d' % (exp_number, scan_number))
            return

        # Generate data by writing out temp file
        ret_obj = self._prepare_view_merged(exp_number, scan_number)
        assert len(ret_obj) == 5
        md_file_name, weight_peak_centers, weight_peak_intensities, avg_peak_centre, avg_peak_intensity = ret_obj

        print 'Write file to %s' % md_file_name
        for i_peak in xrange(len(weight_peak_centers)):
            peak_i = weight_peak_centers[i_peak]
            print '%f, %f, %f' % (peak_i[0], peak_i[1], peak_i[2])
        print
        print avg_peak_centre

        # Plot
        if self._my3DWindow is None:
            self._my3DWindow = plot3dwindow.Plot3DWindow(self)

        self._my3DWindow.add_plot_by_file(md_file_name)
        self._my3DWindow.add_plot_by_array(weight_peak_centers, weight_peak_intensities)
        self._my3DWindow.add_plot_by_array(avg_peak_centre, avg_peak_intensity)

        # Show
        self._my3DWindow.show()

        return

    def _prepare_view_merged(self, exp_number, scan_number):
        """
        Prepare the data to view 3D for merged scan
        :return:
        """
        # check
        assert isinstance(exp_number, int) and isinstance(scan_number, int)

        # get file name for 3D data
        base_file_name = 'md_%d.dat' % random.randint(1, 1001)
        md_file_name = self._myControl.export_md_data(exp_number, scan_number, base_file_name)
        peak_info = self._myControl.get_peak_info(exp_number, scan_number)

        # peak center
        weight_peak_centers, weight_peak_intensities = peak_info.get_weighted_peak_centres()
        qx, qy, qz = peak_info.get_peak_centre()
        # [NEXT/FUTURE] Use a real peak intensity other than 100000
        intensity = 100000

        # convert from list to ndarray
        num_pt_peaks = len(weight_peak_centers)
        assert num_pt_peaks == len(weight_peak_intensities)

        pt_peak_centre_array = numpy.ndarray(shape=(num_pt_peaks, 3), dtype='float')
        pt_peak_intensity_array = numpy.ndarray(shape=(num_pt_peaks,), dtype='float')
        for i_pt in xrange(num_pt_peaks):
            for j in xrange(3):
                pt_peak_centre_array[i_pt][j] = weight_peak_centers[i_pt][j]
            pt_peak_intensity_array[i_pt] = weight_peak_intensities[i_pt]

        avg_peak_centre = numpy.ndarray(shape=(1, 3), dtype='float')
        avg_peak_intensity = numpy.ndarray(shape=(1,), dtype='float')
        avg_peak_centre[0][0] = qx
        avg_peak_centre[0][1] = qy
        avg_peak_centre[0][2] = qz
        avg_peak_intensity[0] = intensity

        return md_file_name, pt_peak_centre_array, pt_peak_intensity_array, avg_peak_centre, avg_peak_intensity

    def do_view_merged_scans_3d(self):
        """ Get selected merged scan and plot them individually in 3D
        :return:
        """
        # collect the selected rows and thus workspace names
        row_index_list = self.ui.tableWidget_mergeScans.get_selected_rows(True)
        exp_scan_list = list()
        for row_index in row_index_list:
            exp_number, scan_number = self.ui.tableWidget_mergeScans.get_exp_scan(row_index)
            pt_number_list = self._myControl.get_pt_numbers(exp_number, scan_number)
            md_data = self._myControl.get_merged_data(exp_number, scan_number, pt_number_list)
            exp_scan_list.append((scan_number, md_data))

        # initialize 3D plot
        if self._my3DWindow is None:
            self._my3DWindow = plot3dwindow.Plot3DWindow(self)
        self._my3DWindow.set_merged_data_set(exp_scan_list)
        self._my3DWindow.show()

        return

    def do_view_ub(self):
        """ View UB matrix in tab 'UB matrix'
        :return:
        """
        self.ui.tabWidget.setCurrentIndex(MainWindow.TabPage['UB Matrix'])

        return

    def do_view_survey_peak(self):
        """ View selected peaks from survey table
        Requirements: one and only 1 run is selected
        Guarantees: the scan number and pt number that are selected will be set to
            tab 'View Raw' and the tab is switched.
        :return:
        """
        # get values
        try:
            scan_num, pt_num = self.ui.tableWidget_surveyTable.get_selected_run_surveyed()
        except RuntimeError, err:
            self.pop_one_button_dialog(str(err))
            return

        # clear selection
        self.ui.tableWidget_surveyTable.select_all_rows(False)

        # switch tab
        self.ui.tabWidget.setCurrentIndex(MainWindow.TabPage['View Raw Data'])
        self.ui.lineEdit_run.setText(str(scan_num))
        self.ui.lineEdit_rawDataPtNo.setText(str(pt_num))

        return

    def pop_one_button_dialog(self, message):
        """ Pop up a one-button dialog
        :param message:
        :return:
        """
        assert isinstance(message, str)
        QtGui.QMessageBox.information(self, '4-circle Data Reduction', message)

        return

    def save_current_session(self, filename=None):
        """ Save current session/value setup to
        :return:
        """
        # Set up dictionary
        save_dict = dict()

        # Setup
        save_dict['lineEdit_localSpiceDir'] = str(self.ui.lineEdit_localSpiceDir.text())
        save_dict['lineEdit_url'] = str(self.ui.lineEdit_url.text())
        save_dict['lineEdit_workDir']= str(self.ui.lineEdit_workDir.text())

        # Experiment
        save_dict['lineEdit_exp'] = str(self.ui.lineEdit_exp.text())
        save_dict['lineEdit_scanNumber'] = self.ui.lineEdit_scanNumber.text()

        # Lattice
        save_dict['lineEdit_a'] = str(self.ui.lineEdit_a.text())
        save_dict['lineEdit_b'] = str(self.ui.lineEdit_b.text())
        save_dict['lineEdit_c'] = str(self.ui.lineEdit_c.text())
        save_dict['lineEdit_alpha'] = str(self.ui.lineEdit_alpha.text())
        save_dict['lineEdit_beta'] = str(self.ui.lineEdit_beta.text())
        save_dict['lineEdit_gamma'] = str(self.ui.lineEdit_gamma.text())

        # Merge scan
        # FIXME : use self._myUBMatrix instead
        #         save_dict['plainTextEdit_ubInput'] = str(self.ui.plainTextEdit_ubInput.toPlainText())
        save_dict['lineEdit_listScansSliceView'] = str(self.ui.lineEdit_listScansSliceView.text())
        save_dict['lineEdit_baseMergeMDName'] = str(self.ui.lineEdit_baseMergeMDName.text())

        # Save to csv file
        if filename is None:
            # default
            save_dir = os.path.expanduser('~/.mantid/')
            if os.path.exists(save_dir) is False:
                os.mkdir(save_dir)
            filename = os.path.join(save_dir, 'session_backup.csv')
        # END-IF

        ofile = open(filename, 'w')
        writer = csv.writer(ofile)
        for key, value in save_dict.items():
            writer.writerow([key, value])
        ofile.close()

        return

    def load_session(self, filename=None):
        """
        To load a session, i.e., read it back:
        :param filename:
        :return:
        """
        if filename is None:
            filename = 'session_backup.csv'
            filename = os.path.join(os.path.expanduser('~/.mantid/'), filename)

        in_file = open(filename, 'r')
        reader = csv.reader(in_file)
        my_dict = dict(x for x in reader)

        # set the data from saved file
        for key, value in my_dict.items():
            if key.startswith('lineEdit') is True:
                self.ui.__getattribute__(key).setText(value)
            elif key.startswith('plainText') is True:
                self.ui.__getattribute__(key).setPlainText(value)
            elif key.startswith('comboBox') is True:
                self.ui.__getattribute__(key).setCurrentIndex(int(value))
            else:
                self.pop_one_button_dialog('Error! Widget name %s is not supported' % key)
        # END-FOR

        # set the experiment
        self._myControl.set_local_data_dir(str(self.ui.lineEdit_localSpiceDir.text()))
        self._myControl.set_working_directory(str(self.ui.lineEdit_workDir.text()))
        self._myControl.set_server_url(str(self.ui.lineEdit_url.text()))

        return

    def menu_quit(self):
        """

        :return:
        """
        self.save_settings()
        self.close()

    def show_scan_pt_list(self):
        """ Show the range of Pt. in a scan
        :return:
        """
        # Get parameters
        status, inp_list = gutil.parse_integers_editors([self.ui.lineEdit_exp, self.ui.lineEdit_run])
        if status is False:
            self.pop_one_button_dialog(inp_list)
            return
        else:
            exp_no = inp_list[0]
            scan_no = inp_list[1]

        status, ret_obj = self._myControl.get_pt_numbers(exp_no, scan_no)

        # Form message
        if status is False:
            # Failed to get Pt. list
            error_message = ret_obj
            self.pop_one_button_dialog(error_message)
        else:
            # Form message
            pt_list = sorted(ret_obj)
            num_pts = len(pt_list)
            info = 'Exp %d Scan %d has %d Pt. ranging from %d to %d.\n' % (exp_no, scan_no, num_pts,
                                                                           pt_list[0], pt_list[-1])
            num_miss_pt = pt_list[-1] - pt_list[0] + 1 - num_pts
            if num_miss_pt > 0:
                info += 'There are %d Pt. skipped.\n' % num_miss_pt

            self.pop_one_button_dialog(info)

        return

    def set_ub_peak_table(self, peak_info):
        """
        Set up the table of peaks to calculate UB matrix
        Requirements: peak_info is a valid PeakInfo instance
        :param peak_info:
        :return:
        """
        # Check requirements
        assert isinstance(peak_info, r4c.PeakInfo)

        # Get data
        exp_number, scan_number = peak_info.get_experiment_info()
        h, k, l = peak_info.get_user_hkl()
        q_x, q_y, q_z = peak_info.get_peak_centre()
        m1 = self._myControl.get_sample_log_value(exp_number, scan_number, 1, '_m1')

        # Set to table
        status, err_msg = self.ui.tableWidget_peaksCalUB.append_row(
            [scan_number, -1, h, k, l, q_x, q_y, q_z, False, m1, ''])
        if status is False:
            self.pop_one_button_dialog(err_msg)

        return

    def save_settings(self):
        """
        Save settings (parameter set) upon quiting
        :return:
        """
        settings = QtCore.QSettings()

        # directories
        local_spice_dir = str(self.ui.lineEdit_localSpiceDir.text())
        settings.setValue("local_spice_dir", local_spice_dir)
        work_dir = str(self.ui.lineEdit_workDir.text())
        settings.setValue('work_dir', work_dir)

        # experiment number
        exp_num = str(self.ui.lineEdit_exp.text())
        settings.setValue('exp_number', exp_num)

        # lattice parameters
        lattice_a = str(self.ui.lineEdit_a.text())
        settings.setValue('a', lattice_a)
        lattice_b = str(self.ui.lineEdit_b.text())
        settings.setValue('b', lattice_b)
        lattice_c = str(self.ui.lineEdit_c.text())
        settings.setValue('c', lattice_c)
        lattice_alpha = str(self.ui.lineEdit_alpha.text())
        settings.setValue('alpha', lattice_alpha)
        lattice_beta = str(self.ui.lineEdit_beta.text())
        settings.setValue('beta', lattice_beta)
        lattice_gamma = str(self.ui.lineEdit_gamma.text())
        settings.setValue('gamma', lattice_gamma)

        return

    def load_settings(self):
        """
        Load QSettings from previous saved file
        :return:
        """
        settings = QtCore.QSettings()

        # directories
        try:
            spice_dir = settings.value('local_spice_dir', '')
            self.ui.lineEdit_localSpiceDir.setText(str(spice_dir))
            work_dir = settings.value('work_dir')
            self.ui.lineEdit_workDir.setText(str(work_dir))

            # experiment number
            exp_num = settings.value('exp_number')
            self.ui.lineEdit_exp.setText(str(exp_num))

            # lattice parameters
            lattice_a = settings.value('a')
            self.ui.lineEdit_a.setText(str(lattice_a))
            lattice_b = settings.value('b')
            self.ui.lineEdit_b.setText(str(lattice_b))
            lattice_c = settings.value('c')
            self.ui.lineEdit_c.setText(str(lattice_c))
            lattice_alpha = settings.value('alpha')
            self.ui.lineEdit_alpha.setText(str(lattice_alpha))
            lattice_beta = settings.value('beta')
            self.ui.lineEdit_beta.setText(str(lattice_beta))
            lattice_gamma = settings.value('gamma')
            self.ui.lineEdit_gamma.setText(str(lattice_gamma))
        except TypeError as e:
            self.pop_one_button_dialog(str(e))
            return

        return

    def _get_lattice_parameters(self):
        """
        Get lattice parameters from GUI
        :return: (Boolean, Object).  True, 6-tuple as a, b, c, alpha, beta, gamm
                                     False: error message
        """
        status, ret_list = gutil.parse_float_editors([self.ui.lineEdit_a,
                                                      self.ui.lineEdit_b,
                                                      self.ui.lineEdit_c,
                                                      self.ui.lineEdit_alpha,
                                                      self.ui.lineEdit_beta,
                                                      self.ui.lineEdit_gamma])
        if status is False:
            err_msg = ret_list
            err_msg = 'Unable to parse unit cell due to %s' % err_msg
            return False, err_msg

        a, b, c, alpha, beta, gamma = ret_list

        return True, (a, b, c, alpha, beta, gamma)

    def _plot_raw_xml_2d(self, exp_no, scan_no, pt_no):
        """ Plot raw workspace from XML file for a measurement/pt.
        """
        # Check and load SPICE table file
        does_exist = self._myControl.does_spice_loaded(exp_no, scan_no)
        if does_exist is False:
            # Download data
            status, error_message = self._myControl.download_spice_file(exp_no, scan_no, over_write=False)
            if status is True:
                status, error_message = self._myControl.load_spice_scan_file(exp_no, scan_no)
                if status is False and self._allowDownload is False:
                    self.pop_one_button_dialog(error_message)
                    return
            else:
                self.pop_one_button_dialog(error_message)
                return
        # END-IF(does_exist)

        # Load Data for Pt's xml file
        does_exist = self._myControl.does_raw_loaded(exp_no, scan_no, pt_no)

        if does_exist is False:
            # Check whether needs to download
            status, error_message = self._myControl.download_spice_xml_file(scan_no, pt_no, exp_no=exp_no)
            if status is False:
                self.pop_one_button_dialog(error_message)
                return
            # Load SPICE xml file
            status, error_message = self._myControl.load_spice_xml_file(exp_no, scan_no, pt_no)
            if status is False:
                self.pop_one_button_dialog(error_message)
                return

        # Convert a list of vector to 2D numpy array for imshow()
        # Get data and plot
        raw_det_data = self._myControl.get_raw_detector_counts(exp_no, scan_no, pt_no)
        self.ui.graphicsView.clear_canvas()
        self.ui.graphicsView.add_plot_2d(raw_det_data, x_min=0, x_max=256, y_min=0, y_max=256,
                                         hold_prev_image=False)

        # Information
        info = '%-10s: %d\n%-10s: %d\n%-10s: %d\n' % ('Exp', exp_no,
                                                      'Scan', scan_no,
                                                      'Pt', pt_no)
        self.ui.plainTextEdit_rawDataInformation.setPlainText(info)

        return
