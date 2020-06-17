# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantid workbench.
#
from enum import Enum
from typing import Dict, List

from mantid.api import MatrixWorkspace, MultipleExperimentInfos, SpecialCoordinateSystem
from mantid.plots.datafunctions import get_indices
from mantid.simpleapi import BinMD
import numpy as np

from .sliceinfo import SliceInfo
from .transform import NonOrthogonalTransform

# Constants
PROJ_MATRIX_LOG_NAME = "W_MATRIX"


class WS_TYPE(Enum):
    MDE = 0
    MDH = 1
    MATRIX = 2


class SliceViewerModel:
    """Store the workspace to be plotted. Can be MatrixWorkspace, MDEventWorkspace or MDHistoWorkspace"""
    def __init__(self, ws):
        if isinstance(ws, MatrixWorkspace):
            if ws.getNumberHistograms() < 2:
                raise ValueError("workspace must contain at least 2 spectrum")
            if ws.getDimension(0).getNBins() < 2:  # Check number of x bins
                raise ValueError("workspace must contain at least 2 bin")
        elif isinstance(ws, MultipleExperimentInfos):
            if ws.isMDHistoWorkspace():
                if len(ws.getNonIntegratedDimensions()) < 2:
                    raise ValueError("workspace must have at least 2 non-integrated dimensions")
            else:
                if ws.getNumDims() < 2:
                    raise ValueError("workspace must have at least 2 dimensions")
        else:
            raise ValueError("only works for MatrixWorkspace and MDWorkspace")

        self._ws = ws

        if self.get_ws_type() == WS_TYPE.MDE:
            self.get_ws = self.get_ws_MDE
            self.get_data = self.get_data_MDE
        else:
            self.get_ws = self._get_ws
            self.get_data = self.get_data_MDH

    def can_normalize_workspace(self) -> bool:
        if self.get_ws_type() == WS_TYPE.MATRIX and not self._get_ws().isDistribution():
            return True
        return False

    def get_ws_name(self) -> str:
        return self._ws.name()

    def can_support_nonorthogonal_axes(self) -> bool:
        """
        Query if the workspace can support non-orthogonal axes.
        Workspace must be an MD workspace, in HKL coordinate system
        and have a defined oriented lattice on the first experiment info.
        """
        ws_type = self.get_ws_type()
        if ws_type == WS_TYPE.MDE or ws_type == WS_TYPE.MDH:
            if self.get_frame() != SpecialCoordinateSystem.HKL:
                return False

            sample = self._get_ws().getExperimentInfo(0).sample()
            if not sample.hasOrientedLattice():
                return False

            return True
        else:
            return False

    def can_support_peaks_overlays(self) -> bool:
        """
        Check if the given workspace can support peaks overlay
        """
        try:
            if self._get_ws().getNumDims() > 2:
                return True
        except AttributeError:
            pass

        return False

    def get_frame(self) -> SpecialCoordinateSystem:
        """Return the coordinate system of the workspace"""
        return self._ws.getSpecialCoordinateSystem()

    def get_ws_MDE(self, slicepoint, bin_params):
        params = {}
        for n in range(self._get_ws().getNumDims()):
            if slicepoint[n] is None:
                params['AlignedDim{}'.format(n)] = '{},{},{},{}'.format(
                    self._get_ws().getDimension(n).name,
                    self._get_ws().getDimension(n).getMinimum(),
                    self._get_ws().getDimension(n).getMaximum(), bin_params[n])
            else:
                params['AlignedDim{}'.format(n)] = '{},{},{},{}'.format(
                    self._get_ws().getDimension(n).name, slicepoint[n] - bin_params[n] / 2,
                    slicepoint[n] + bin_params[n] / 2, 1)
        return BinMD(InputWorkspace=self._get_ws(),
                     OutputWorkspace=self._get_ws().name() + '_rebinned',
                     **params)

    def get_data_MDH(self, slicepoint, transpose=False):
        indices, _ = get_indices(self.get_ws(), slicepoint=slicepoint)
        if transpose:
            return np.ma.masked_invalid(self.get_ws().getSignalArray()[indices]).T
        else:
            return np.ma.masked_invalid(self.get_ws().getSignalArray()[indices])

    def get_data_MDE(self, slicepoint, bin_params, transpose=False):
        if transpose:
            return np.ma.masked_invalid(
                self.get_ws_MDE(slicepoint, bin_params).getSignalArray().squeeze()).T
        else:
            return np.ma.masked_invalid(
                self.get_ws_MDE(slicepoint, bin_params).getSignalArray().squeeze())

    def get_dim_limits(self, slicepoint, transpose):
        """
        Return a xlim, ylim) for the display dimensions where xlim, ylim are tuples
        :param slicepoint: Sequence containing either a float or None where None indicates a display dimension
        :param transpose: A boolean flag indicating if the display dimensions are transposed
        """
        workspace = self._get_ws()
        assert len(slicepoint) == workspace.getNumDims(
        ), "Expected len(slicepoint) to match number of workspace dimensions"
        limits = []
        for index, pt in enumerate(slicepoint):
            if pt is None:
                dimension = workspace.getDimension(index)
                limits.append((dimension.getMinimum(), dimension.getMaximum()))
        assert len(
            limits) == 2, f"There should be exactly 2 display dimensions, found {len(limits)}"
        xlim, ylim = limits
        if transpose:
            ylim, xlim = xlim, ylim

        return xlim, ylim

    def get_dim_info(self, n: int) -> dict:
        """
        returns dict of (minimum, maximun, number_of_bins, width, name, units) for dimension n
        """
        workspace = self._get_ws()
        dim = workspace.getDimension(n)
        return {
            'minimum': dim.getMinimum(),
            'maximum': dim.getMaximum(),
            'number_of_bins': dim.getNBins(),
            'width': dim.getBinWidth(),
            'name': dim.name,
            'units': dim.getUnits(),
            'type': self.get_ws_type().name,
            'qdim': dim.getMDFrame().isQ()
        }

    def get_dimensions_info(self) -> List[Dict]:
        """
        returns a list of dict for each dimension containing dim_info
        """
        return [self.get_dim_info(n) for n in range(self._get_ws().getNumDims())]

    def get_ws_type(self) -> WS_TYPE:
        if isinstance(self._get_ws(), MatrixWorkspace):
            return WS_TYPE.MATRIX
        elif isinstance(self._get_ws(), MultipleExperimentInfos):
            if self._get_ws().isMDHistoWorkspace():
                return WS_TYPE.MDH
            else:
                return WS_TYPE.MDE
        else:
            raise ValueError("Unsupported workspace type")

    def create_nonorthogonal_transform(self, slice_info: SliceInfo):
        """
        Calculate the transform object for nonorthogonal axes.
        :param slice_info: SliceInfoType object giving information on the slice
        :raises: RuntimeError if model does not support nonorthogonal axes
        """
        if not self.can_support_nonorthogonal_axes():
            raise RuntimeError("Workspace cannot support non-orthogonal axes view")

        expt_info = self._get_ws().getExperimentInfo(0)
        lattice = expt_info.sample().getOrientedLattice()
        try:
            proj_matrix = np.array(expt_info.run().get(PROJ_MATRIX_LOG_NAME).value)
            proj_matrix = proj_matrix.reshape(3, 3)
        except KeyError:
            # assume orthogonal projection
            proj_matrix = np.diag([1., 1., 1.])

        display_indices = slice_info.transform([0, 1, 2])
        return NonOrthogonalTransform.from_lattice(lattice,
                                                   x_proj=proj_matrix[:, display_indices[0]],
                                                   y_proj=proj_matrix[:, display_indices[1]])

    # private api
    def _get_ws(self):
        return self._ws
