# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +

from mantid.api import PythonAlgorithm, MatrixWorkspaceProperty, MultipleFileProperty, Progress
from mantid.kernel import Direction, IntBoundedValidator
from mantid.simpleapi import *


class IndirectILLReductionDIFF(PythonAlgorithm):
    """
    Performs reduction on IN16B's diffraction data. It can be on mode Doppler or BATS.
    """

    runs = None
    mode = None
    scan_parameter = None
    mask_start_pixels = None
    mask_end_pixels = None
    output = None
    progress = None

    def category(self):
        return "ILL\\Indirect"

    def summary(self):
        return "Performs reduction on IN16B's diffraction data. Mode is either Doppler or BATS."

    def name(self):
        return "IndirectILLReductionDIFF"

    def setUp(self):
        self.runs = self.getPropertyValue('SampleRuns').split(',')
        self.scan_parameter = self.getPropertyValue('Observable')
        self.mask_start_pixels = self.getProperty('MaskPixelsFromStart').value
        self.mask_end_pixels = self.getProperty('MaskPixelsFromEnd').value
        self.output = self.getPropertyValue('OutputWorkspace')
        self.progress = Progress(self, start=0.0, end=1.0, nreports=10)

    def PyInit(self):
        self.declareProperty(MultipleFileProperty('SampleRuns', extensions=['nxs']), doc="File path for run(s).")

        self.declareProperty(MatrixWorkspaceProperty('OutputWorkspace', '', direction=Direction.Output),
                             doc='The output workspace group containing reduced data.')

        self.declareProperty("MaskPixelsFromStart", 10, validator=IntBoundedValidator(lower=0),
                             doc="Number of pixels to mask at the start of each tube")
        self.declareProperty("MaskPixelsFromEnd", 10, validator=IntBoundedValidator(lower=0),
                             doc="Number of pixels to mask at the end of each tube")

        self.declareProperty("Observable", "sample.temperature",
                             doc="If multiple files, the parameter from SampleLog to use as an index when conjoined.")

    def _normalize_by_monitor(self, ws):
        """
            Normalizes the workspace by monitor values (ID is 0 for IN16B)
            @param ws : the input workspace
        """
        monitorID = 0
        monitor_ws = ws + '_mon'
        ExtractSpectra(InputWorkspace=ws, DetectorList=monitorID, OutputWorkspace=monitor_ws)

        # in case of 0 counts monitors, replace 0s by 1s so the division becomes neutral
        # (since it's generally division of 0 detector's counts by 0 monitor's counts,
        # they weren't very useful to begin with)

        ReplaceSpecialValues(InputWorkspace=monitor_ws, OutputWorkspace=monitor_ws, SmallNumberThreshold=0.00000001,
                             SmallNumberValue=1)

        Divide(LHSWorkspace=ws, RHSWorkspace=monitor_ws, OutputWorkspace=ws, WarnOnZeroDivide=True)

        cache = list(range(1, self.mask_start_pixels)) + list(range(257 - self.mask_end_pixels, 257))
        to_mask = [i + 256 * j for i in cache for j in range(8)]
        MaskDetectors(Workspace=ws, DetectorList=to_mask)
        DeleteWorkspace(monitor_ws)

    def _treat_doppler(self, ws):
        """
            Reduce Doppler diffraction data presents in workspace.
            @param ws: the input workspace
        """
        run = None
        if len(self.runs) > 1:
            number_of_channels = mtd[mtd[ws].getNames()[0]].blocksize()
            run = mtd[mtd[ws].getNames()[0]].getRun()
        else:
            number_of_channels = mtd[ws].blocksize()
            run = mtd[ws].getRun()

        if run.hasProperty('Doppler.incident_energy'):
            energy = run.getLogData('Doppler.incident_energy').value / 1000
        else:
            raise RuntimeError("Unable to find incident energy for Doppler mode")

        Rebin(InputWorkspace=ws, OutputWorkspace=self.output, Params=[0, number_of_channels, number_of_channels])
        self._normalize_by_monitor(self.output)

        ConvertSpectrumAxis(InputWorkspace=self.output,
                            OutputWorkspace=self.output,
                            Target='ElasticQ',
                            EMode="Direct",
                            EFixed=energy)

        ConvertToPointData(InputWorkspace=self.output, OutputWorkspace=self.output)

        ConjoinXRuns(InputWorkspaces=self.output,
                     SampleLogAsXAxis=self.scan_parameter,
                     FailBehaviour="Stop",
                     OutputWorkspace="conjoined_" + self.output)
        mtd[self.output].delete()

        ExtractUnmaskedSpectra(InputWorkspace="conjoined_" + self.output, OutputWorkspace=self.output)
        mtd["conjoined_" + self.output].delete()

        Transpose(InputWorkspace=self.output, OutputWorkspace=self.output)

    def _treat_BATS(self, ws):
        self.log().warning("BATS treatment not implemented yet.")
        pass

    def PyExec(self):
        self.setUp()
        LoadAndMerge(Filename=self.getPropertyValue('SampleRuns'), OutputWorkspace=self.output,
                     LoaderOptions={"LoadDetectors": "Diffractometer"}, startProgress=0, endProgress=0.4)

        if len(self.runs) > 1:
            run = mtd[mtd[self.output].getNames()[0]].getRun()
        else:
            run = mtd[self.output].getRun()

        if run.hasProperty('acquisition_mode') and run.getLogData('acquisition_mode').value == 1:
            self.mode = "BATS"
            self.log().information("Mode recognized as BATS.")
        else:
            self.mode = "Doppler"
            self.log().information("Mode recognized as Doppler.")

        self.progress.report(4, "Treating data")
        if self.mode == "Doppler":
            self._treat_doppler(self.output)
        elif self.mode == "BATS":
            self._treat_BATS(self.output)

        self.setProperty("OutputWorkspace", mtd[self.output])


AlgorithmFactory.subscribe(IndirectILLReductionDIFF)
