/*WIKI*
 Integrates arbitrary shaped single crystal peaks defined on an [[MDHistoWorkspace]] using connected component analysis to determine
 regions of interest around each peak of the [[PeaksWorkspace]]. The output is an integrated [[PeaksWorkspace]] as well as an image
 containing the labels assigned to each cluster for diagnostic and visualisation purposes.

 '''The algorithm makes no assumptions about Peak shape or size''' and can therefore be used where integration over defined shapes
 [[IntegratePeaksMD]] and [[IntegrateEllipsoids]], for example, will not work.

 [[File:ClusterImage.png|400px]]

 ''Cluster Label region displayed in the [[SliceViewer]]. Peak centre is marked with an X. The green circle illustrates the integration region used by [[IntegratePeaksMD]]''
 
 A threshold for the Peak should be defined below which, parts of the image are treated as background. The normalization method in combination with the
 threshold may both be used to define a background. We suggest keeping the default of VolumeNormalization so that changes in the effective bin size
 do not affect the background filtering.

 This algorithm uses an imaging technique, and it is therefore important that the MDHistoWorkspace you are using is binned to a sufficient
 resolution via [[BinMD]]. You can overlay the integrated peaks workspace in the [[MantidPlot:_SliceViewer#Viewing_Peaks_Workspaces|Slice Viewer]] over
 the generated Cluster Labeled OutputWorkspaceMD to see what the integration region used for each peak amounts to.

 == Notes for running ==

 It is suggested that you '''initially run the algorithm on a coarse image'''. This will help you tune the Threshold parameters. The algorithm generates
 a large memory footprint, so it is suggested that you keep the initial image small, and run on hardware with sufficient memory to store multiple workspace
 of equal size to the input MDWorkspace (generated as part of the connected component analysis).

 == Warnings and Logging ==
 The algorithm will generate warning. There are three main warning to know about.
 === Off the Image Edge ===
 The algorithm will warn about unreachable peaks (off the image). This may be because the peaks detected were off
 the edge of the detector, or because the image was cropped in BinMD in such a way that that part of the detector/TOF space is no
 longer accessible.
 === No Cluster Corresponding to Peak ===
 This is because the input [[PeaksWorkspace]] has peaks that do not align with peaks in the image. The error could either
 be on the side of the input PeaksWorkspace (spurious peaks), or of the [[MDHistoWorkspace]] generated as part of processing. One thing to verify
 is that the combination of Threshold and Normalization input parameters are not so low that they are treating genuine peaks in the image
 as background.
 === Multiple Peaks Assigned to the same Cluster ===
 This means overlapping peaks in the image. This is a problem because both peaks will be given an integrated value that is the sum
 of the entire cluster. You may need to increase the Threshold parameter to resolve this problem.

 For more in-depth analysis, the algorithm will produce debug log messages.

 *WIKI*/

#include "MantidCrystal/IntegratePeaksUsingClusters.h"
#include "MantidCrystal/ICluster.h"
#include "MantidCrystal/ConnectedComponentLabeling.h"
#include "MantidCrystal/HardThresholdBackground.h"
#include "MantidCrystal/PeakClusterProjection.h"
#include "MantidAPI/IMDHistoWorkspace.h"
#include "MantidAPI/WorkspaceProperty.h"
#include "MantidAPI/IMDIterator.h"
#include "MantidAPI/AlgorithmManager.h"
#include "MantidAPI/Progress.h"
#include "MantidKernel/MultiThreaded.h"
#include "MantidKernel/CompositeValidator.h"
#include "MantidKernel/MandatoryValidator.h"
#include "MantidKernel/BoundedValidator.h"
#include "MantidKernel/ListValidator.h"
#include "MantidKernel/Utils.h"
#include "MantidDataObjects/PeaksWorkspace.h"

#include <boost/make_shared.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <map>
#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <cmath>

using namespace Mantid::API;
using namespace Mantid::Kernel;
using namespace Mantid::DataObjects;
using namespace Mantid::Crystal::ConnectedComponentMappingTypes;

namespace Mantid
{
  namespace Crystal
  {

    // Register the algorithm into the AlgorithmFactory
    DECLARE_ALGORITHM(IntegratePeaksUsingClusters)

    //----------------------------------------------------------------------------------------------
    /** Constructor
     */
    IntegratePeaksUsingClusters::IntegratePeaksUsingClusters()
    {
    }

    //----------------------------------------------------------------------------------------------
    /** Destructor
     */
    IntegratePeaksUsingClusters::~IntegratePeaksUsingClusters()
    {
    }

    //----------------------------------------------------------------------------------------------
    /// Algorithm's name for identification. @see Algorithm::name
    const std::string IntegratePeaksUsingClusters::name() const
    {
      return "IntegratePeaksUsingClusters";
    }
    ;

    /// Algorithm's version for identification. @see Algorithm::version
    int IntegratePeaksUsingClusters::version() const
    {
      return 1;
    }
    ;

    /// Algorithm's category for identification. @see Algorithm::category
    const std::string IntegratePeaksUsingClusters::category() const
    {
      return "MDAlgorithms";
    }

    //----------------------------------------------------------------------------------------------
    /// Sets documentation strings for this algorithm
    void IntegratePeaksUsingClusters::initDocs()
    {
      this->setWikiSummary("Integrate single crystal peaks using connected component analysis");
      this->setOptionalMessage(this->getWikiSummary());
    }

    //----------------------------------------------------------------------------------------------
    /** Initialize the algorithm's properties.
     */
    void IntegratePeaksUsingClusters::init()
    {
      declareProperty(new WorkspaceProperty<IMDHistoWorkspace>("InputWorkspace", "", Direction::Input),
          "Input md workspace.");
      declareProperty(new WorkspaceProperty<IPeaksWorkspace>("PeaksWorkspace", "", Direction::Input),
          "A PeaksWorkspace containing the peaks to integrate.");

      auto positiveValidator = boost::make_shared<BoundedValidator<double> >();
      positiveValidator->setExclusive(true);
      positiveValidator->setLower(0);

      auto compositeValidator = boost::make_shared<CompositeValidator>();
      compositeValidator->add(positiveValidator);
      compositeValidator->add(boost::make_shared<MandatoryValidator<double> >());

      declareProperty(
          new PropertyWithValue<double>("Threshold", 0, compositeValidator, Direction::Input),
          "Threshold signal above which to consider peaks");

      std::vector<std::string> normalizations(3);
      normalizations[0] = "NoNormalization";
      normalizations[1] = "VolumeNormalization";
      normalizations[2] = "NumEventsNormalization";

      declareProperty("Normalization", normalizations[1],
          Kernel::IValidator_sptr(new Kernel::ListValidator<std::string>(normalizations)),
          "Normalization to use with Threshold. Defaults to VolumeNormalization to account for different binning.");

      declareProperty(new WorkspaceProperty<IPeaksWorkspace>("OutputWorkspace", "", Direction::Output),
          "An output integrated peaks workspace.");
      declareProperty(
          new WorkspaceProperty<IMDHistoWorkspace>("OutputWorkspaceMD", "", Direction::Output),
          "MDHistoWorkspace containing the labeled clusters used by the algorithm.");
    }

    /**
     * Get the normalization. For use with iterators + background strategies.
     * @return Chosen normalization
     */
    MDNormalization IntegratePeaksUsingClusters::getNormalization()
    {
      std::string normProp = getPropertyValue("Normalization");
      Mantid::API::MDNormalization normalization;
      if (normProp == "NoNormalization")
      {
        normalization = NoNormalization;
      }
      else if (normProp == "VolumeNormalization")
      {
        normalization = VolumeNormalization;
      }
      else
      {
        normalization = NumEventsNormalization;
      }
      return normalization;
    }

    //----------------------------------------------------------------------------------------------
    /** Execute the algorithm.
     */
    void IntegratePeaksUsingClusters::exec()
    {
      IMDHistoWorkspace_sptr mdWS = getProperty("InputWorkspace");
      IPeaksWorkspace_sptr inPeakWS = getProperty("PeaksWorkspace");
      IPeaksWorkspace_sptr peakWS = getProperty("OutputWorkspace");
      if (peakWS != inPeakWS)
      {
        peakWS = IPeaksWorkspace_sptr(dynamic_cast<IPeaksWorkspace*>(inPeakWS->clone()));
      }

      {
        const SpecialCoordinateSystem mdCoordinates = mdWS->getSpecialCoordinateSystem();
        if (mdCoordinates == None)
        {
          throw std::invalid_argument(
              "The coordinate system of the input MDWorkspace cannot be established. Run SetSpecialCoordinates on InputWorkspace.");
        }
      }

      const double threshold = getProperty("Threshold");
      // Make a background strategy for the CCL analysis to use.
      HardThresholdBackground backgroundStrategy(threshold, this->getNormalization());
      // CCL. Multi-processor version.
      ConnectedComponentLabeling analysis;

      Progress progress(this, 0, 1, 1);
      // Perform CCL.
      ClusterTuple clusters = analysis.executeAndFetchClusters(mdWS, &backgroundStrategy, progress);
      // Extract the clusters
      ConnectedComponentMappingTypes::ClusterMap& clusterMap = clusters.get<1>();
      // Extract the labeled image
      IMDHistoWorkspace_sptr outHistoWS = clusters.get<0>();
      // Labels taken by peaks.
      std::map<size_t, size_t> labelsTakenByPeaks;
      // Make a peak transform so that we can understand a peak in the context of the mdworkspace coordinate setup.
      PeakClusterProjection projection(outHistoWS); // Projection of PeaksWorkspace over labelled cluster workspace.

      progress.doReport("Performing Peak Integration");
      g_log.information("Starting Integration");
      progress.resetNumSteps(peakWS->getNumberPeaks(), 0.9, 1);
      PARALLEL_FOR1(peakWS)
      for (int i = 0; i < peakWS->getNumberPeaks(); ++i)
      {
        PARALLEL_START_INTERUPT_REGION
        IPeak& peak = peakWS->getPeak(i);
        const Mantid::signal_t signalValue = projection.signalAtPeakCenter(peak); // No normalization when extracting label ids!
        if (boost::math::isnan(signalValue))
        {
          g_log.warning() << "Warning: image for integration is off edge of detector for peak " << i
              << std::endl;
        }
        else if (signalValue < static_cast<Mantid::signal_t>(analysis.getStartLabelId()))
        {
          g_log.information() << "Peak: " << i
              << " Has no corresponding cluster/blob detected on the image. This could be down to your Threshold settings."
              << std::endl;
        }
        else
        {
          const size_t labelIdAtPeak = static_cast<size_t>(signalValue);
          ICluster * const cluster = clusterMap[labelIdAtPeak].get();
          ICluster::ClusterIntegratedValues integratedValues = cluster->integrate(mdWS);
          peak.setIntensity(integratedValues.get<0>());
          peak.setSigmaIntensity(std::sqrt(integratedValues.get<1>()));

          PARALLEL_CRITICAL(IntegratePeaksUsingClusters)
          {
            auto it = labelsTakenByPeaks.find(labelIdAtPeak);
            if (it != labelsTakenByPeaks.end())
            {
              g_log.warning() << "Overlapping Peaks. Peak: " << i << " overlaps with another Peak: "
                  << it->second << " and shares label id: " << it->first << std::endl;
            }
            labelsTakenByPeaks.insert(std::make_pair(labelIdAtPeak, i));
          }
          progress.report();
        }
      PARALLEL_END_INTERUPT_REGION
    }
    PARALLEL_CHECK_INTERUPT_REGION

    setProperty("OutputWorkspace", peakWS);
    setProperty("OutputWorkspaceMD", outHistoWS);
  }

} // namespace Crystal
} // namespace Mantid
