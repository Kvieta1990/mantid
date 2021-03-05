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


def __calculate_strain(obs, exp: np.ndarray):
    return np.asarray(list(obs.values())[1:-2]) / exp


def __calculate_rel_diff(obs, exp: np.ndarray):
    obs_ndarray = np.asarray(list(obs.values())[1:-2])
    return (obs_ndarray - exp) / exp


def collect_peaks(wksp, outputname: str, donor=None, infotype: str = 'strain'):
    if infotype not in ['strain', 'reldiff']:
        raise ValueError('Do not know how to calculate "{}"'.format(infotype))

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
        if infotype == 'strain':
            output.setY(i, __calculate_strain(wksp.row(i), peaks))
        elif infotype == 'reldiff':
            output.setY(i, __calculate_rel_diff(wksp.row(i), peaks))

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


def plot_peakd(wksp, peak_positions):
    """
    Plots peak d spacing value for each peak position in peaks
    :param wksp: Workspace returned from collect_peaks
    :param peak_positions: List of peak positions
    :return: plot, plot axes
    """

    peaks = peak_positions
    if isinstance(peak_positions, float):
        peaks = [peak_positions]

    if len(peaks) == 0:
        raise ValueError("Expected one or more peak positions")

    if not mtd.doesExist(str(wksp)):
        raise ValueError("Could not find provided workspace in ADS")

    fig, ax = plt.subplots()
    ax.set_xlabel("det IDs")

    # Hold the mean and stddev for each peak to compute total at end
    means = []
    stddevs = []

    cm = plt.get_cmap("jet")
    ax.set_prop_cycle(color=[cm(1.*i/len(peaks)) for i in range(len(peaks))])

    # Plot data for each peak position
    for peak in peaks:
        print("Processing peak position {}".format(peak))
        single = extract_peak_info(wksp, 'single', peak)

        # get x and y arrays from single peak ws
        x = single.dataX(0)
        y = single.dataY(0)

        # filter out any nans
        y_val = y[~np.isnan(y)]

        # skip if y was entirely nans
        if len(y_val) == 0:
            continue

        means.append(np.mean(y_val))
        stddevs.append(np.std(y_val))

        ax.plot(x, y, marker="x", linestyle="None", label="{:0.6f}".format(peak))
        ax.legend(bbox_to_anchor=(1, 1), loc="upper left")

    # If every peak had nans, raise error
    if len(means) == 0 or len(stddevs) == 0:
        raise RuntimeError("No valid peak data was found for provided peak positions")

    # Calculate total mean and stddev of all peaks
    total_mean = np.mean(means)
    total_stddev = np.std(stddevs)

    # Draw solid line at mean
    ax.axhline(total_mean, color="black", lw=2.5)  # default lw=1.5

    # Upper and lower lines calibration lines (1 percent of mean)
    band = total_mean * 0.01
    ax.axhline(total_mean + band, color="black", ls="--")
    ax.axhline(total_mean - band, color="black", ls="--")

    # Add mean and stddev text annotations
    stat_str = "Mean = {:0.6f} Stdev = {:0.6f}".format(total_mean, total_stddev)
    plt_text = ax.text(0, 0, stat_str, transform=ax.transAxes)

    # Center text in the figure
    box = plt_text.get_window_extent(renderer=fig.canvas.get_renderer()).inverse_transformed(ax.transAxes)
    text_x = 0.5 - box.width * 0.5
    plt_text.set_position((text_x, 0.875))

    plt.show()

    return fig, ax
