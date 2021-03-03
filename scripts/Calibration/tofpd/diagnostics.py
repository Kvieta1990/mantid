from mantid.plots.resampling_image.samplingimage import imshow_sampling
from mantid.plots.datafunctions import get_axes_labels
from mantid.simpleapi import mtd
from mantid.api import WorkspaceFactory
import matplotlib.pyplot as plt
import numpy as np


# Diamond peak positions in d-space which may differ from actual sample
DIAMOND = np.asarray((2.06,1.2615,1.0758,0.892,0.8186,0.7283,0.6867,0.6307,0.5642,0.5441,0.515,0.4996,0.4768,
                      0.4645,0.4205,0.3916,0.3499,0.3257,0.3117))


def _get_xrange(wksp, xmarkers, tolerance, xmin, xmax):
    # start with the range determined by the markers and tolerance
    x_mins = [xmarkers.min() * (1. - 3 * tolerance)]
    x_maxs = [xmarkers.max() * (1. + 3 * tolerance)]
    # add range based on user supplied parameters
    if not np.isnan(xmin):
        x_mins.append(xmin)
    if not np.isnan(xmax):
        x_maxs.append(xmax)
    # add data range if possible
    try:
        temp = wksp.getTofMin()  # only a method on EventWorkspace
        if not np.isnan(temp):
            x_mins.append(temp)
        temp = wksp.getTofMax()
        if not np.isnan(temp):
            x_maxs.append(temp)
    except:
        pass  # don't use the data range
    xmin = np.max(x_mins)
    xmax = np.min(x_maxs)

    return xmin, xmax


def plot2d(workspace, tolerance: float=0.001, peakpositions: np.ndarray=DIAMOND,
           xmin: float=np.nan, xmax: float=np.nan, horiz_markers=[]):
    TOLERANCE_COLOR = 'w'  # color to mark the area within tolerance
    MARKER_COLOR = 'b'  # color for peak markers and horizontal (between bank) markers

    # convert workspace to pointer to workspace
    wksp = mtd[str(workspace)]

    # determine x-range
    xmin, xmax = _get_xrange(wksp, peakpositions, tolerance, xmin, xmax)

    # create a new figure
    fig, ax = plt.subplots()
    # show the image - this uses the same function used by "colorfill" plot
    # and slice viewer. The data is subsampled according to how may screen
    # pixels are available
    im = imshow_sampling(ax, wksp, aspect='auto', origin='lower')
    # add the colorbar
    fig.colorbar(im, ax=ax)
    # set the x-range to show
    ax.set_xlim(xmin, xmax)

    # annotate expected peak positions and tolerances
    for peak in peakpositions:
        if peak > 0.4 and peak < 1.5:
            shade = ((1-tolerance) * peak, (1+tolerance) * peak)
            ax.axvspan(shade[0], shade[1], color=TOLERANCE_COLOR, alpha=.2)
            ax.axvline(x=peak, color=MARKER_COLOR, alpha=0.4)

    # annotate areas of the detector
    for position in horiz_markers:
        ax.axhline(position, color=MARKER_COLOR, alpha=0.4)

    # add axes labels
    labels = get_axes_labels(wksp)
    ax.set_xlabel(labels[1])
    ax.set_ylabel(labels[2])

    # show the figure
    fig.show()

    # return the figure so others can customize it
    return fig, fig.axes


def collect_peaks(wksp, outputname: str, donor=None):
    wksp = mtd[str(wksp)]

    numSpec = int(wksp.rowCount())
    peak_names = [item for item in wksp.getColumnNames()
                  if item not in ['detid', 'chisq', 'normchisq']]
    peaks = np.asarray([float(item[1:]) for item in peak_names])
    numPeaks = len(peaks)

    # convert the d-space table to a Workspace2d
    if donor:
        donor = mtd[str(donor)]
    else:
        donor = 'Workspace2D'
    output = WorkspaceFactory.create(donor, NVectors=numSpec,
                                     XLength=numPeaks, YLength=numPeaks)
    output.getAxis(0).setUnit('dSpacing')
    for i in range(numSpec):  # TODO get the detID correct
        output.setX(i, peaks)
        output.setY(i, list(wksp.row(i).values())[1:-2] / peaks)

    # add the workspace to the AnalysisDataService
    mtd.addOrReplace(outputname, output)
    return mtd[outputname]


def extract_peak_info(wksp, outputname: str, peak_position: float):
    '''
    Extract information about a single peak from a Workspace2D. The input workspace is expected to have
    common x-axis of observed d-spacing. The y-values and errors are extracted.

    The output workspace will be a single spectra with the x-axis being the detector-id. The y-values
    and errors are extracted from the input workspace.
    '''
    # confirm that the input is a workspace pointer
    wksp = mtd[str(wksp)]
    numSpec = wksp.getNumberHistograms()

    # get the index into the x/y arrays of the peak position
    peak_index = wksp.readX(0).searchsorted(peak_position)

    # create a workspace to put the result into
    single = WorkspaceFactory.create('Workspace2D', NVectors=1,
                                     XLength=wksp.getNumberHistograms(),
                                     YLength=wksp.getNumberHistograms())
    single.setTitle('d-spacing={}\\A'.format(wksp.readX(0)[peak_index]))

    # get a handle to map the detector positions
    detids = wksp.detectorInfo().detectorIDs()
    have_detids = bool(len(detids) > 0)

    # fill in the data values
    x = single.dataX(0)
    y = single.dataY(0)
    e = single.dataE(0)
    for wksp_index in range(numSpec):
        if have_detids:
            x[wksp_index] = detids[wksp_index]
        else:
            x[wksp_index] = wksp_index
        y[wksp_index] = wksp.readY(wksp_index)[peak_index]
        e[wksp_index] = wksp.readE(wksp_index)[peak_index]

    # add the workspace to the AnalysisDataService
    mtd.addOrReplace(outputname, single)
    return mtd[outputname]
