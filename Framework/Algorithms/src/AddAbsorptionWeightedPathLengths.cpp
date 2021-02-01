// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2020 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidAlgorithms/AddAbsorptionWeightedPathLengths.h"
#include "MantidAPI/ExperimentInfo.h"
#include "MantidAPI/Sample.h"
#include "MantidAPI/WorkspaceProperty.h"
#include "MantidAlgorithms/SampleCorrections/CircularBeamProfile.h"
#include "MantidAlgorithms/SampleCorrections/MCAbsorptionStrategy.h"
#include "MantidAlgorithms/SampleCorrections/MCInteractionStatistics.h"
#include "MantidAlgorithms/SampleCorrections/MCInteractionVolume.h"
#include "MantidAlgorithms/SampleCorrections/RectangularBeamProfile.h"
#include "MantidDataObjects/PeaksWorkspace.h"
#include "MantidGeometry/Instrument.h"
#include "MantidGeometry/Instrument/ReferenceFrame.h"
#include "MantidGeometry/Instrument/SampleEnvironment.h"
#include "MantidKernel/BoundedValidator.h"
#include "MantidKernel/EnabledWhenProperty.h"
#include "MantidKernel/Material.h"
#include "MantidKernel/MersenneTwister.h"

using namespace Mantid::API;
using namespace Mantid::DataObjects;
using namespace Mantid::Geometry;
using namespace Mantid::Kernel;

namespace Mantid {
namespace Algorithms {

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(AddAbsorptionWeightedPathLengths)

//----------------------------------------------------------------------------------------------

namespace {

constexpr int DEFAULT_NEVENTS = 1000;
constexpr int DEFAULT_SEED = 123456789;
constexpr int NLAMBDA = 1;

} // namespace

//----------------------------------------------------------------------------------------------
/** Initialize the algorithm's properties.
 */
void AddAbsorptionWeightedPathLengths::init() {
  declareProperty(
      std::make_unique<WorkspaceProperty<PeaksWorkspace_sptr::element_type>>("InputWorkspace", "", Direction::InOut),
      "An input/output peaks workspace that the path distances will be added "
      "to.");
  declareProperty("UseSinglePath", false, "Use a single path with a scatter point at the sample position");
  auto positiveInt = std::make_shared<Kernel::BoundedValidator<int>>();
  positiveInt->setLower(1);
  declareProperty("EventsPerPoint", DEFAULT_NEVENTS, positiveInt,
                  "The number of \"neutron\" events to generate per peak");
  declareProperty("SeedValue", DEFAULT_SEED, positiveInt, "Seed the random number generator with this value");
  declareProperty("MaxScatterPtAttempts", 5000, positiveInt,
                  "Maximum number of tries made to generate a scattering point "
                  "within the sample. Objects with "
                  "holes in them, e.g. a thin annulus can cause problems "
                  "if this number is too low.\n"
                  "If a scattering point cannot be generated by increasing "
                  "this value then there is most likely a problem with "
                  "the sample geometry.");
  setPropertySettings("SeedValue",
                      std::make_unique<EnabledWhenProperty>("UseSinglePath", ePropertyCriterion::IS_DEFAULT));
  setPropertySettings("EventsPerPoint",
                      std::make_unique<EnabledWhenProperty>("UseSinglePath", ePropertyCriterion::IS_DEFAULT));
  setPropertySettings("MaxScatterPtAttempts",
                      std::make_unique<EnabledWhenProperty>("UseSinglePath", ePropertyCriterion::IS_DEFAULT));
}

std::map<std::string, std::string> AddAbsorptionWeightedPathLengths::validateInputs() {
  PeaksWorkspace_sptr inputWS = getProperty("InputWorkspace");
  std::map<std::string, std::string> issues;
  Geometry::IComponent_const_sptr sample = inputWS->getInstrument()->getSample();
  if (!sample) {
    issues["InputWorkspace"] = "Input workspace does not have a Sample";
  } else {
    if (inputWS->sample().hasEnvironment()) {
      issues["InputWorkspace"] = "Sample must not have a sample environment";
    }

    if (inputWS->sample().getMaterial().numberDensity() == 0) {
      issues["InputWorkspace"] = "Sample must have a material set up with a non-zero number density";
    }
  }
  return issues;
}

//----------------------------------------------------------------------------------------------
/** Execute the algorithm.
 */
void AddAbsorptionWeightedPathLengths::exec() {

  const PeaksWorkspace_sptr inputWS = getProperty("InputWorkspace");
  const int nevents = getProperty("EventsPerPoint");
  const int maxScatterPtAttempts = getProperty("MaxScatterPtAttempts");

  auto instrument = inputWS->getInstrument();
  auto beamProfile = createBeamProfile(*instrument, inputWS->sample());

  const auto npeaks = inputWS->getNumberPeaks();

  // Configure progress
  Progress prog(this, 0.0, 1.0, npeaks);
  prog.setNotifyStep(0.01);
  const std::string reportMsg = "Computing path lengths";

  // Configure strategy
  MCInteractionVolume interactionVol(inputWS->sample(), maxScatterPtAttempts);
  MCAbsorptionStrategy strategy(interactionVol, *beamProfile, DeltaEMode::Elastic, nevents, maxScatterPtAttempts, true);

  const int seed = getProperty("SeedValue");

  PARALLEL_FOR_IF(Kernel::threadSafe(*inputWS))
  for (int i = 0; i < npeaks; ++i) {
    PARALLEL_START_INTERUPT_REGION
    IPeak &peak = inputWS->getPeak(i);
    auto peakWavelength = peak.getWavelength();

    std::vector<double> lambdas{peakWavelength}, absFactors(NLAMBDA), absFactorErrors(NLAMBDA);

    bool useSinglePath = getProperty("UseSinglePath");
    if (useSinglePath) {
      auto inst = inputWS->getInstrument();
      const auto sourcePos = inst->getSource()->getPos();
      const auto samplePos = inst->getSample()->getPos();
      const auto reverseBeamDir = normalize(samplePos - sourcePos);

      const IObject *sampleShape = &(inputWS->sample().getShape());

      Track beforeScatter(samplePos, reverseBeamDir);
      sampleShape->interceptSurface(beforeScatter);
      const auto detDir = normalize(peak.getDetPos() - samplePos);
      Track afterScatter(samplePos, detDir);
      sampleShape->interceptSurface(afterScatter);

      absFactors[0] = beforeScatter.calculateAttenuation(lambdas[0]) * afterScatter.calculateAttenuation(lambdas[0]);

    } else {
      MersenneTwister rng(seed + int(i));
      MCInteractionStatistics detStatistics(peak.getDetectorID(), inputWS->sample());

      strategy.calculate(rng, peak.getDetectorPosition(), lambdas, peakWavelength, absFactors, absFactorErrors,
                         detStatistics);

      if (g_log.is(Kernel::Logger::Priority::PRIO_DEBUG)) {
        g_log.debug(detStatistics.generateScatterPointStats());
      }
    }

    double mu = inputWS->sample().getMaterial().attenuationCoefficient(peakWavelength); // m-1
    double absWeightedPathLength = -log(absFactors[0]) / mu;                            // metres
    peak.setAbsorptionWeightedPathLength(absWeightedPathLength * 100);                  // cm

    prog.report(reportMsg);
    PARALLEL_END_INTERUPT_REGION
  }
  PARALLEL_CHECK_INTERUPT_REGION
}

/**
 * Create the beam profile. Currently only supports Rectangular and Circular.
 * The dimensions are either specified by those provided by `SetBeam` algorithm
 * or default to the width and height of the samples bounding box
 * @param instrument A reference to the instrument object
 * @param sample A reference to the sample object
 * @return A new IBeamProfile object
 */
std::unique_ptr<IBeamProfile> AddAbsorptionWeightedPathLengths::createBeamProfile(const Instrument &instrument,
                                                                                  const Sample &sample) const {
  const auto frame = instrument.getReferenceFrame();
  const auto source = instrument.getSource();
  double beamWidth(-1.0), beamHeight(-1.0), beamRadius(-1.0);

  std::string beamShapeParam = source->getParameterAsString("beam-shape");
  if (beamShapeParam.compare("Slit") == 0) {
    auto beamWidthParam = source->getNumberParameter("beam-width");
    auto beamHeightParam = source->getNumberParameter("beam-height");
    if (beamWidthParam.size() == 1 && beamHeightParam.size() == 1) {
      beamWidth = beamWidthParam[0];
      beamHeight = beamHeightParam[0];
      return std::make_unique<RectangularBeamProfile>(*frame, source->getPos(), beamWidth, beamHeight);
    }
  } else if (beamShapeParam.compare("Circle") == 0) {
    auto beamRadiusParam = source->getNumberParameter("beam-radius");
    if (beamRadiusParam.size() == 1) {
      beamRadius = beamRadiusParam[0];
      return std::make_unique<CircularBeamProfile>(*frame, source->getPos(), beamRadius);
    }
  } // revert to sample dimensions if no return by this point
  const auto bbox = sample.getShape().getBoundingBox().width();
  beamWidth = bbox[frame->pointingHorizontal()];
  beamHeight = bbox[frame->pointingUp()];
  return std::make_unique<RectangularBeamProfile>(*frame, source->getPos(), beamWidth, beamHeight);
}

} // namespace Algorithms
} // namespace Mantid