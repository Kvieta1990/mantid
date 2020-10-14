# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +

from copy import deepcopy
import numpy as np
import re
from typing import Any, Callable, Dict, Optional, Tuple
from corelli.calibration.utils import InputTable, WorkspaceTypes  # custom type aliases
# imports from Mantid
from mantid import AnalysisDataService, mtd
from mantid.api import TextAxis, WorkspaceGroup
from mantid.dataobjects import TableWorkspace, Workspace2D
from mantid.simpleapi import (CloneWorkspace, CreateEmptyTableWorkspace, CreateWorkspace,
                              DeleteTableRows, DeleteWorkspaces, GroupWorkspaces, MaskBTP, RenameWorkspace)
from Calibration import tube
from Calibration.tube_spec import TubeSpec
from Calibration.tube_calib_fit_params import TubeCalibFitParams
from corelli.calibration.utils import bank_numbers, PIXELS_PER_TUBE, TUBES_IN_BANK, wire_positions


def sufficient_intensity(input_workspace: WorkspaceTypes, bank_name: str, minimum_intensity:float = 10000) -> bool:
    r"""
    Assert if the average intensity per pixel in the bank surpasses a minimum threshold.

    :param input_workspace: input Workspace2D containing total neutron counts per pixel
    :param bank_name: a string of the form 'bankI' where 'I' is a bank number
    :minimum_intensity: neutron counts averaged over all pixels
    """
    workspace = mtd[str(input_workspace)]
    tube_set = TubeSpec(workspace)
    tube_set.setTubeSpecByString(bank_name)
    workspace_indexes = list()
    for tube_index in range(tube_set.getNumTubes()):
        workspace_indexes_in_tube, skipped = tube_set.getTube(tube_index)
        workspace_indexes.extend(list(workspace_indexes_in_tube))
    return bool(np.mean(workspace.extractY()[workspace_indexes].flatten()) > minimum_intensity)


def fit_bank(workspace: WorkspaceTypes, bank_name: str, shadow_height: float = 1000, shadow_width: float = 4,
             fit_domain: float = 7, minimum_intensity: float = 1000,
             calibration_table: str = 'CalibTable', peak_pixel_positions_table: str = 'PeakTable') -> None:
    r"""
    Find the position of the wire shadow on each tube in the bank, in units of pixel positions

    :param workspace: input Workspace2D containing total neutron counts per pixel
    :param bank_name: a string of the form 'bankI' where 'I' is a bank number
    :param shadow_height: initial guess for the decrease in neutron counts caused by the calibrating wire
    :param shadow_width: initial guess for the number of pixels shadowed by the calibrating wire
    :param fit_domain: number of pixels over which a calibrating wire may have any significant influence
    :param minimum_intensity: mininum number of neutron counts per pixel to warrant a significant fit
    session. This number is compared against the neutron counts per pixel, averaged over all pixels in the bank
    :param calibration_table: output TableWorkspace containing one column for detector ID and one column
    for its calibrated XYZ coordinates, in meters
    :param peak_pixel_positions_table: output table containing the positions of the wire shadows for each tube, in
    units of pixel coordinates

    :raises AssertionError: the input workspace is of incorrect type or cannot be found
    :raises AssertionError: the input `shadow_height` is a non-positive value
    :raises AssertionError: the `bank_name` if not of the form 'bankI'
    :raises AssertionError: insufficient neutron counts per pixel

    :return: workspace handles to the calibration and peak table
    """
    message = f'Cannot process workspace {workspace}. Pass the name of an existing workspace or a workspace handle'
    assert isinstance(workspace, (str, Workspace2D)), message
    workspace_name = str(workspace)
    assert AnalysisDataService.doesExist(workspace_name), f'Input workspace {workspace_name} does not exists'
    assert shadow_height > 0, 'shadow height must be positive'
    peak_height, peak_width = -shadow_height, shadow_width
    assert re.match(r'^bank\d+$', bank_name), 'The bank name must be of the form "bankI" where "I" in an integer'
    message = f'Insufficient counts per pixel in workspace {workspace_name} for a confident calibration'
    assert sufficient_intensity(workspace, bank_name, minimum_intensity=minimum_intensity), message
    # Fit only the inner 14 dips because the extrema wires are too close to the tube tips.
    # The dead zone in the tube tips interferes with the shadow cast by the extrema  wires
    # preventing a good fitting
    wire_positions_pixels = wire_positions(units='pixels')[1: -1]
    wire_count = len(wire_positions_pixels)
    peaks_form = [1] * wire_count  # signals we'll be fitting dips (peaks with negative heights)

    fit_par = TubeCalibFitParams(wire_positions_pixels, height=peak_height, width=peak_width, margin=fit_domain)
    fit_par.setAutomatic(True)

    tube.calibrate(workspace_name, bank_name, wire_positions(units='meters')[1: -1],
                   peaks_form, fitPar=fit_par, outputPeak=True)
    if calibration_table != 'CalibTable':
        RenameWorkspace(InputWorkspace='CalibTable', OutputWorkspace=calibration_table)
    if peak_pixel_positions_table != 'PeakTable':
        RenameWorkspace(InputWorkspace='PeakTable', OutputWorkspace=peak_pixel_positions_table)


def criterium_peak_pixel_position(peak_table: InputTable, summary: Optional[str] = None,
                                  zscore_threshold: float = 2.5, deviation_threshold: float = 3.0) -> np.ndarray:
    r"""
    Flag tubes whose peak pixel positions deviate considerably from the peak pixel positions when
    averaged for all tubes in the bank.

    .. math::

      <p_i> = \frac{1}{n_t} \Sum_{j=1}^{n_t} p_{ij}
      \delta_j^2 = \frac{1}{n_w} \Sum (p_{ij} - <p_i>)^2
      assert d_j < threshold

    :param peak_table: pixel positions of the peaks for each tube
    :param summary: name of output Workspace2D containing deviations and Z-score for each tube.
    :param zscore_threshold: maximum Z-score for the pixels positions of a tube.
    :param deviation_threshold: maximum deviation (in pixels) for the pixels positions of a tube.
    :return: array of booleans, one per tube. `True` is the tube passes the acceptance criterium, `False` otherwise.
    """
    table = mtd[str(peak_table)]  # handle to the peak table
    peak_count = table.columnCount() - 1  # the first column contains the names of the tubes
    # `positions_average` stores the pixel position for each peak, averaged for all tubes
    positions_average = [np.mean(table.column(column_number)) for column_number in range(1, 1 + peak_count)]

    deviations = list()  # a measure of how much the peak positions in a tube deviate from the mean positions
    tube_count = table.rowCount()  # number of tubes in the bank
    for tube_index in range(tube_count):
        positions = np.array(list(table.row(tube_index).values())[1:])  # peak positions for the current tube
        deviations.append(np.sqrt(np.mean(np.square(positions - positions_average))))

    # find tubes with a large Z-score
    outlier_values = list()
    values = deepcopy(deviations)
    z_score = 1000
    outlier_value = 1000
    deviation_threshold = 3.0  # three pixels
    while z_score > zscore_threshold and outlier_value > deviation_threshold and len(values) > 0:
        # find the tube with the highest Z-score, possibly signaling a large deviation from the mean
        mean, std = np.mean(values), np.std(values)
        outlier_index = np.argmax(np.abs((values - mean) / std))
        outlier_value = values[outlier_index]
        # recalculate the Z-score of the tube, but removing it from the pool of values. This removes
        # any skewing effects from including the aberrant tube in the calculation of its Z-score
        del values[outlier_index]
        mean, std = np.mean(values), np.std(values)
        z_score = np.abs((outlier_value - mean) / std)
        if z_score > zscore_threshold and outlier_value > deviation_threshold:
            outlier_values.append(outlier_value)

    # flag the outlier tubes as failing the criterium
    criterium_pass = np.tile(True, tube_count)  # initialize as all tubes passing the criterium
    if len(outlier_values) > 0:
        failure_indexes = [deviations.index(value) for value in outlier_values]
        criterium_pass[failure_indexes] = False

    # create an analysis summary if so requested
    if isinstance(summary, str) and len(summary) > 0:
        success = [1 if criterium else 0 for criterium in criterium_pass]
        x_values = list(range(1, 1 + TUBES_IN_BANK))
        mean, std = np.mean(deviations), np.std(deviations)
        z_scores = np.abs((deviations - mean) / std)
        y_values = np.array([success, deviations, z_scores]).flatten()
        workspace = CreateWorkspace(x_values, y_values, NSpec=3, OutputWorkspace=summary,
                                    WorkspaceTitle='Tube deviations from averages taken over the bank',
                                    YUnitLabel='Pixel Units', EnableLogging=False)
        labels = ('success', 'deviation', 'Z-score')
        axis = TextAxis.create(len(labels))
        [axis.setLabel(index, label) for index, label in enumerate(labels)]
        workspace.replaceAxis(1, axis)

    return criterium_pass


def purge_table(workspace: WorkspaceTypes, calibration_table: TableWorkspace,
                tubes_fit_success: np.ndarray,  output_table: str = None) -> None:
    r"""
    Remove the detectorID's corresponding to the failing tubes from the calibration table

    Assumptions:
    - Each tube has PIXELS_PER_TUBE number of pixels
    - detector ID's in the calibration table are sorted according to tube number

    :param workspace: input Workspace2D containing total neutron counts per pixel
    :param calibration_table: input TableWorkspace containing one column for detector ID and one column
    for its calibrated XYZ coordinates, in meters
    :param tubes_fit_success: array of booleans of length TUBES_IN_BANK. `False` if a tube was unsuccessfully fitted.
    :param output_table: name of the purged table. If `None`, the input `calibration_table` is purged.
    """
    # validate input arguments
    if False not in tubes_fit_success:
        return  # nothing to do
    # validate the input workspace
    message = f'Cannot process workspace {workspace}. Pass the name of an existing workspace or a workspace handle'
    assert isinstance(workspace, (str, Workspace2D)), message
    workspace_name = str(workspace)
    assert AnalysisDataService.doesExist(workspace_name), f'Input workspace {workspace_name} does not exists'
    # validate the input calibraton table
    message = f'Cannot process table {calibration_table}. Pass the name of an existing TableWorkspace' \
              ' or a TableWorkspace handle'
    assert isinstance(calibration_table, (str, TableWorkspace)), message
    assert AnalysisDataService.doesExist(str(calibration_table)), f'Input table {calibration_table} does not exists'
    if output_table is not None:
        CloneWorkspace(InputWorkspace=calibration_table, OutputWorkspace=output_table)
    else:
        output_table = str(calibration_table)
    tube_fail_indexes = np.where(tubes_fit_success == False)[0]  # noqa E712 indexes of tubes unsuccessfully fitted
    row_indexes_of_first_tube = np.arange(PIXELS_PER_TUBE)  # 0, 1, ... 255
    fail_rows = [row_indexes_of_first_tube + (i * PIXELS_PER_TUBE) for i in tube_fail_indexes]
    fail_rows = np.array(fail_rows, dtype=int).flatten().tolist()
    DeleteTableRows(output_table, fail_rows)


def mask_bank(bank_name: str, tubes_fit_success: np.ndarray, output_table: str) -> Optional[TableWorkspace]:
    r"""
    Creates a single-column `TableWorkspace` object containing the detector ID's of the
    unsuccessfully fitted tubes

    If all tubes were fit successfully, no `TableWorkspace` is created, and `None` is returned.

    :param bank_name: a string of the form 'bankI' where 'I' is a bank number
    :param tubes_fit_success: array of boolean containing a True/False entry for each tube, indicating wether
    the tube was successfully calibrated.
    :param output_table: name of the output TableWorkspace containing one column for detector ID from tubes
    not successfully calibrated.

    :raise AssertionError: the string bank_name does not follow the pattern 'bankI' where 'I' in an integer
    :return: name of the mask TableWorkspace. Returns `None` if no TableWorkspace is created.
    """
    assert re.match(r'^bank\d+$', bank_name), 'The bank name must be of the form "bankI" where "I" in an integer'
    if False not in tubes_fit_success:
        return None  # al tubes were fit successfully
    bank_number = bank_name[4:]  # drop 'bank' from bank_name
    tube_numbers = 1 + np.where(tubes_fit_success == False)[0]  # noqa E712 unsuccessfully fitted tube numbers
    tube_numbers = ','.join([str(n) for n in tube_numbers])  # failing tubes as a string
    detector_ids = MaskBTP(Instrument='CORELLI', Bank=bank_number, Tube=tube_numbers)
    table = CreateEmptyTableWorkspace(OutputWorkspace=output_table)
    table.addColumn('long64', 'DETECTOR ID')
    [table.addRow([detector_id]) for detector_id in detector_ids.tolist()]
    if AnalysisDataService.doesExist('CORELLIMaskBTP'):
        DeleteWorkspaces(['CORELLIMaskBTP'])
    return mtd[output_table]


def calibrate_bank(workspace: WorkspaceTypes, bank_name: str,
                   calibration_table: str, mask_table: str = 'MaskTable',
                   acceptance_criterium: Callable[[Any], np.ndarray] = criterium_peak_pixel_position,
                   criterium_kwargs: Dict[str, Any] = {},
                   shadow_height: float = 1000, shadow_width: float = 4, fit_domain: float = 7,
                   minimum_intensity: float = 1000) -> Tuple[TableWorkspace, Optional[TableWorkspace]]:
    r"""
    Calibrate the tubes in a bank and assess their goodness-of-fit with an acceptance function. Creates a
    table of calibrated detector IDs and a table of non-calibrated detector IDs

    :param workspace: input Workspace2D containing total neutron counts per pixel
    :param bank_name: a string of the form 'bankI' where 'I' is a bank number
    :param shadow_height: initial guess for the decrease in neutron counts caused by the calibrating wire
    :param shadow_width: initial guess for the number of pixels shadowed by the calibrating wire
    :param fit_domain: number of pixels over which a calibrating wire may have any significant influence
    :param minimum_intensity: mininum number of neutron counts per pixel to warrant a significant fit
    session. This number is compared against the neutron counts per pixel, averaged over all pixels in the bank
    :param acceptance_criterium: a function the determines wether each pixel in the tube was calibrated correctly
    :param criterium_kwargs: optional arguments to the `acceptance_criterium` function.
    :param calibration_table: output TableWorkspace containing one column for detector ID and one column
    for its calibrated XYZ coordinates, in meters
    :param mask_table: output TableWorkspace containing containing the detector ID's of the
    unsuccessfully fitted tubes

    :return: Workspace2D handles for the calibration and mask tables
    """
    # Validate inputs are taken care in function fit_bank
    # Fit the tubes in the bank
    fit_bank(workspace, bank_name, shadow_height, shadow_width, fit_domain, minimum_intensity,
             calibration_table=calibration_table, peak_pixel_positions_table='PeakTable')
    # Run the acceptance criterium to determine the failing tubes
    tubes_fit_success = acceptance_criterium('PeakTable', **criterium_kwargs)
    # purge the calibration table of detector ID's with failing tubes
    purge_table(workspace, calibration_table, tubes_fit_success)
    # Create table of masked detector ID's, or None if all tubes were successfully fitted
    mask_table_workspace = mask_bank(bank_name, tubes_fit_success, mask_table)
    DeleteWorkspaces(['PeakTable'])
    return mtd[calibration_table], mask_table_workspace


def calibrate_banks(workspace: WorkspaceTypes, bank_selection: str,
                    calibration_group: str = 'calibrations', mask_group: str = 'masks',
                    acceptance_group: str = 'acceptances') -> Tuple[WorkspaceGroup, Optional[WorkspaceGroup]]:
    r"""
    Calibrate the tubes in a selection of banks, and assess their goodness-of-fit with an acceptance function.

    For each bank, creates a table of calibrated detector IDs, a table of non-calibrated detector IDs, and one
    'acceptance' Workspace2D object containing, for each tube, the typical deviation for the position of the wire shadows
    with respect to the average positions for the while bank. A Z-score measure of these deviations is also stored.

    :param workspace: input Workspace2D containing total neutron counts per pixel
    :param bank_selection: selection string, such as '32-36,38,40-43'
    :param calibration_group: name of the output WorkspaceGroup containing the calibration tables for each bank
    :param mask_group: name of the output WorkspaceGroup containing the mask tables for each bank
    :param acceptance_group: name of the output WorkspaceGroup containing the acceptance workspaces

    :return: handles to the calibrations and masks WorkspaceGroup objects
    """
    # create a list of banks names

    # Calibrate each bank
    calibrations, masks, acceptances = list(), list(), list()
    for n in bank_numbers(bank_selection):
        calibration, mask = calibrate_bank(workspace, 'bank' + n, 'calib' + n, 'mask' + n,
                                           criterium_kwargs={'summary': 'acceptance' + n})
        acceptances.append(mtd['acceptance' + n])
        calibrations.append(calibration)
        if mask is not None:
            masks.append(mask)
    # Create the group workspaces
    GroupWorkspaces(InputWorkspaces=calibrations, OutputWorkspace=calibration_group)
    GroupWorkspaces(InputWorkspaces=acceptances, OutputWorkspace=acceptance_group)
    if len(masks) > 0:
        GroupWorkspaces(InputWorkspaces=masks, OutputWorkspace=mask_group)

    return mtd[calibration_group], None if len(masks) == 0 else mtd[mask_group]
