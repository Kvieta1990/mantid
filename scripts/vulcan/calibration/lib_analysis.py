# Zoo of methods that are develooped for analyze the calibration
from mantid.simpleapi import (AlignDetectors, FitPeaks, FindPeakBackground, DiffractionFocussing, Rebin,
                              ConvertToMatrixWorkspace, EditInstrumentGeometry, SaveNexusProcessed,
                              MaskDetectors, ConvertUnits)
from mantid.simpleapi import mtd
import numpy as np
from typing import Union, Tuple


class FindDiamondPeaks(object):
    """
    Class to handle diamond peak positions for calibration diagnosis.

    It is supposed to replace diagnose.get_peak_centers()
    """

    def __init(self):

        self._diamond_ws_name = None

    def find_peak(self, peak_name, exp_pos, d_min, d_max, start_ws_index, end_ws_index):

        # fit peaks
        peak_pos_ws_name = self._fit_peaks(exp_pos, d_min, d_max, start_ws_index, end_ws_index, peak_name)

        # separate good fit/bad fit/no fit
        self._analyze_fitted_peaks(peak_pos_ws_name, start_ws_index, end_ws_index)

        # do statistics on fit result (average, standard deviation)

        # ...

    def _fit_peaks(self, expected_peak_pos, min_d, max_d, start_ws_index, end_ws_index, suffix):

        # base workspace name (output)
        tag = f'{self._diamond_ws_name}_{suffix}'

        # output workspace names
        peak_pos_ws_name = f'{tag}_peak_positions'

        # Default to Gaussian and Linear background
        FitPeaks(InputWorkspace=self._diamond_ws_name,
                 OutputWorkspace=peak_pos_ws_name,
                 StartWorkspaceIndex=int(start_ws_index),
                 StopWorkspaceIndex=int(end_ws_index),
                 PeakCenters=f'{expected_peak_pos}',
                 FitWindowBoundaryList=f'{min_d}, {max_d}',
                 FittedPeaksWorkspace=f'{tag}_model',
                 OutputPeakParametersWorkspace=f'{tag}_params',
                 OutputParameterFitErrorsWorkspace=f'{tag}_errors')

        return peak_pos_ws_name

    def _analyze_fitted_peaks(self, peak_pos_ws_name: str,
                              start_ws_index: int,
                              end_ws_index: int):
        # . -1 for data is zero;
        # -2 for maximum value is smaller than specified minimum value.and
        # -3 for non-converged fitting.

        peak_pos_ws = mtd[peak_pos_ws_name]

        ws_index_vec = np.arange(start_ws_index, end_ws_index)
        # ws_index_vec = np.arange(6468, 24900)

        # Sanity check
        print(f'Number of histograms in {peak_pos_ws_name} = {peak_pos_ws.getNumberHistograms()}, '
              f'Spectrum vector shape = {ws_index_vec.shape}')
        # peak positions
        peak_center_vec = peak_pos_ws.extractY().flatten()
        assert ws_index_vec.shape == peak_center_vec.shape

        # Type 1: zero counts
        zero_counts_pixels = np.where(np.abs(peak_center_vec - (-1)) < 1E-5)[0]
        print(f'Zero counts = {len(zero_counts_pixels)}')

        # Type 2: low counts
        low_counts_pixels = np.where(np.abs(peak_center_vec - (-2)) < 1E-5)[0]
        print(f'Low counts = {len(low_counts_pixels)}')

        bad_fit_pixels = np.where(np.abs(peak_center_vec - (-3)) < 1E-5)[0]
        print(f'Type 1 bad counts = {len(bad_fit_pixels)}')
        # print(f'  they are ... {bad_fit_pixels + 6468}')

        bad_fit_pixels = np.where(np.abs(peak_center_vec - (-4)) < 1E-5)[0]
        print(f'Type 2 bad counts = {len(bad_fit_pixels)}')
        # print(f'  they are ... {bad_fit_pixels + 6468}')

        # Check again with counts
        diamond_ws = mtd[self._diamond_ws_name]
        counts_vec = diamond_ws.extractY().sum(axis=1)[start_ws_index:end_ws_index]
        # sanity check
        assert counts_vec.shape == peak_center_vec.shape

        # Check zero cunts
        assert np.max(counts_vec[peak_center_vec - (-1)]) < 1E-6

    def get_peak_height_max_y(self):
        # Get peak height as maximum Y value by removing background

        bkgd_param_ws_name = 'allmultiple'
        FindPeakBackground(InputWorkspace='VULCAN_164960_matrix', WorkspaceIndex=10000, FitWindow='1,1.2',
                           OutputWorkspace=bkgd_param_ws_name)

        bkgd_param_ws = mtd[bkgd_param_ws_name]
        ws_index = bkgd_param_ws.cell(0, 0)
        peak_min_index = bkgd_param_ws.cell(0, 1)
        peak_max_index = bkgd_param_ws.cell(0, 2)
        bgkd_0 = bkgd_param_ws.cell(0, 3)
        bkgd_1 = bkgd_param_ws.cell(0, 4)
        bkgd_2 = bkgd_param_ws.cell(0, 5)

        print(f'Spectrum (index) {ws_index} from d-index = {peak_min_index}, {peak_max_index}: '
              f'B = {bgkd_0} + {bkgd_1} * d + {bkgd_2} * d**2')


def report_masked_pixels(data_workspace, mask_ws, wi_start, wi_stop):
    """Generate a report for all masked pixels

    Parameters
    ----------
    data_workspace: str, MatrixWorkspace
        A diamond workspace.  It can be either the raw EventWorkspace or a Workspace2D with 1 bin (count)
    mask_ws: str, MaskWorkspace
        Mask workspace
    wi_start: int, None
        starting workspace index
    wi_stop: int, None
        stopping workspace index (excluded)

    Returns
    -------

    """
    # Get input
    if isinstance(data_workspace, str):
        data_workspace = mtd[data_workspace]
    if isinstance(mask_ws, str):
        mask_ws = mtd[mask_ws]
    if wi_start is None:
        wi_start = 0
    if wi_stop is None:
        wi_stop = mask_ws.getNumberHistograms()

    # Check input
    # get the number of events per spectrum
    events_number_vec = data_workspace.extractY().sum(axis=1).flatten()
    assert mask_ws.getNumberHistograms() == events_number_vec.shape[0],\
        f'Number of histograms does not match.  counts vector shape = {events_number_vec.shape}'

    # Set initial value for statistics
    num_masked = 0
    zero_masked = 0
    event_spectrum_list = list()
    for ws_index in range(wi_start, wi_stop):
        # skip non-masked pixels
        if mask_ws.readY(ws_index)[0] < 0.1:
            continue
        else:
            num_masked += 1

        # analyze masking information
        if events_number_vec[ws_index] == 0:
            zero_masked += 1
        else:
            event_spectrum_list.append((events_number_vec[ws_index], ws_index))

    # Make report
    report = '[REPORT]'
    report += '\nFrom {} to {}: Number of masked pixels = {}, including '.format(wi_start, wi_stop-1, num_masked)
    report += f'\n  (1) {zero_masked} pixels with zero counts'
    report += f'\n  (2) {num_masked - zero_masked} pixels with non-zero counts and they are ... '
    event_spectrum_list.sort(reverse=True)
    for i in range(min(100, len(event_spectrum_list))):
        num_events_i, ws_index = event_spectrum_list[i]
        report += '\n      ws-index = {}, num of events = {}'.format(ws_index, num_events_i)

    return report


def align_focus_event_ws(event_ws_name,
                         calib_ws_name: Union[str, None],
                         group_ws_name: str,
                         mask_ws_name: Union[str, None]) -> Tuple[str, str]:
    """
    overwrite the input
    """
    # determine tag
    file_tag = ''

    # Align detector or not
    print(f'Event workspace: {event_ws_name}.  X unit = {mtd[event_ws_name].getAxis(0).getUnit().unitID()}')

    if calib_ws_name:
        # align detectors and convert unit to dSpacing
        AlignDetectors(InputWorkspace=event_ws_name, OutputWorkspace=event_ws_name,
                       CalibrationWorkspace=calib_ws_name)
        file_tag += '_Cal'

    else:
        # optionally not align detectors: convert to dSpacing
        ConvertUnits(InputWorkspace=event_ws_name, OutputWorkspace=event_ws_name, Target='dSpacing')
        file_tag += '_Raw'

    # Rebin
    Rebin(InputWorkspace=event_ws_name, OutputWorkspace=event_ws_name, Params='0.3,-0.0003,1.5')
    # Convert to matrix workspace
    matrix_ws_name = f'{event_ws_name}_matrix'
    ConvertToMatrixWorkspace(InputWorkspace=event_ws_name, OutputWorkspace=matrix_ws_name)

    # Save nexus for 2D alignment view
    SaveNexusProcessed(InputWorkspace=matrix_ws_name, Filename=f'{event_ws_name}{file_tag}.nxs')
    print(f'[CHECK] saved aligned workspace size: {mtd[matrix_ws_name].extractY().shape}')

    # Mask group workspace
    if mask_ws_name:
        MaskDetectors(Workspace=group_ws_name, MaskedWorkspace=mask_ws_name)
        file_tag += 'Masked'
    else:
        file_tag += '_Nomask'

    # Diffraction focus
    DiffractionFocussing(InputWorkspace=event_ws_name, OutputWorkspace=event_ws_name,
                         GroupingWorkspace=group_ws_name)

    # Convert from event workspace to workspace 2D
    ConvertToMatrixWorkspace(InputWorkspace=event_ws_name, OutputWorkspace=event_ws_name)

    # Edit instrument geometry
    # NOTE: Disable EditInstrumentGeometry as
    #   1.  The geometry information won't be saved to processed NeXus
    #   2.  It destroys the geometry information that can be used for FitPeaks with instrument parameters
    # EditInstrumentGeometry(Workspace=event_ws_name, PrimaryFlightPath=42, SpectrumIDs='1-3', L2='2,2,2',
    #                        Polar='89.9284,90.0716,150.059', Azimuthal='0,0,0', DetectorIDs='1-3',
    #                        InstrumentName='vulcan_3bank')

    # TODO - need to relax to allow user to  determine
    focused_run_nxs = f'{event_ws_name}{file_tag}_3banks.nxs'

    SaveNexusProcessed(InputWorkspace=event_ws_name, Filename=focused_run_nxs)

    return event_ws_name, focused_run_nxs


def get_masked_ws_indexes(mask_ws):
    """
    get the workspace indexes that are masked
    :param mask_ws:
    :return:
    """
    if isinstance(mask_ws, str):
        mask_ws = mtd[mask_ws]

    masked_list = list()
    for iws in range(mask_ws.getNumberHistograms()):
        if mask_ws.readY(iws)[0] > 0.5:
            masked_list.append(iws)

    return masked_list
