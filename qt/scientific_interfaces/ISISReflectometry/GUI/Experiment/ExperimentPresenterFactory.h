// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
#ifndef MANTID_ISISREFLECTOMETRY_EXPERIMENTPRESENTERFACTORY_H
#define MANTID_ISISREFLECTOMETRY_EXPERIMENTPRESENTERFACTORY_H
#include "../../Reduction/Experiment.h"
#include "Common/DllConfig.h"
#include "ExperimentPresenter.h"
#include "IExperimentPresenter.h"
#include "IExperimentView.h"
#include <memory>

namespace MantidQt {
namespace CustomInterfaces {

class ExperimentPresenterFactory {
public:
  explicit ExperimentPresenterFactory(double thetaTolerance)
      : m_thetaTolerance(thetaTolerance) {}

  std::unique_ptr<IExperimentPresenter> make(IExperimentView *view) {
    return std::make_unique<ExperimentPresenter>(view, makeModel(),
                                                 m_thetaTolerance);
  }

private:
  double m_thetaTolerance;

  Experiment makeModel() {
    // TODO inject model instead of creating it here
    auto polarizationCorrections =
        PolarizationCorrections(PolarizationCorrectionType::None);
    auto floodCorrections(FloodCorrectionType::Workspace);
    auto stitchParameters = std::map<std::string, std::string>();
    auto perThetaDefaults = std::vector<PerThetaDefaults>(
        {PerThetaDefaults(boost::none, std::pair<std::string, std::string>(),
                          RangeInQ(), boost::none, ProcessingInstructions())});
    return Experiment(AnalysisMode::PointDetector, ReductionType::Normal,
                      SummationType::SumInLambda, false, false,
                      std::move(polarizationCorrections),
                      std::move(floodCorrections), boost::none,
                      std::move(stitchParameters), std::move(perThetaDefaults));
  }
};
} // namespace CustomInterfaces
} // namespace MantidQt
#endif // MANTID_ISISREFLECTOMETRY_EXPERIMENTPRESENTERFACTORY_H
