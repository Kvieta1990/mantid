// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2009 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidDataHandling/ReadMaterial.h"
#include "MantidAPI/Algorithm.h"
#include "MantidKernel/Material.h"

namespace Mantid {
namespace DataHandling {

/**
 * Validate the parameters to build the material from, this returns
 * any errors in the inputs.
 *
 * @param params A struct containing all the parameters to be set.
 * @returns A map containing the relevent failure messages, if any.
 */
ValidationErrors
ReadMaterial::validateInputs(const MaterialParameters &params) {
  ValidationErrors result;
  const bool chemicalSymbol{!params.chemicalSymbol.empty()};
  const bool atomicNumber{params.atomicNumber != 0};
  if (!chemicalSymbol && !atomicNumber) {
    if (isEmpty(params.coherentXSection)) {
      result["CoherentXSection"] = "The cross section must be specified when "
                                   "no ChemicalFormula or AtomicNumber is "
                                   "given.";
    }
    if (isEmpty(params.incoherentXSection)) {
      result["IncoherentXSection"] = "The cross section must be specified when "
                                     "no ChemicalFormula or AtomicNumber is "
                                     "given.";
    }
    if (isEmpty(params.attenuationXSection) &&
        params.attenuationProfileFileName.empty()) {
      result["AttenuationXSection"] = "The cross section must be specified "
                                      "when no ChemicalFormula or AtomicNumber "
                                      "is given.";
    }
    if (isEmpty(params.scatteringXSection)) {
      result["ScatteringXSection"] = "The cross section must be specified when "
                                     "no ChemicalFormula or AtomicNumber is "
                                     "given.";
    }
    if (isEmpty(params.numberDensity)) {
      result["NumberDensity"] =
          "The number density must be specified with a use-defined material.";
    }
  } else if (chemicalSymbol && atomicNumber) {
    result["AtomicNumber"] =
        "Cannot specify both ChemicalFormula and AtomicNumber";
  }

  if (params.massNumber > 0 && params.atomicNumber <= 0)
    result["AtomicNumber"] = "Specified MassNumber without AtomicNumber";

  if (!isEmpty(params.zParameter)) {
    if (isEmpty(params.unitCellVolume)) {
      result["UnitCellVolume"] =
          "UnitCellVolume must be provided with ZParameter";
    }
    if (!isEmpty(params.numberDensity)) {
      result["ZParameter"] = "Cannot give ZParameter with NumberDensity set";
    }
    if (!isEmpty(params.massDensity)) {
      result["MassDensity"] = "Cannot give MassDensity with ZParameter set";
    }
  } else if (!isEmpty(params.numberDensity)) {
    if (!isEmpty(params.massDensity)) {
      result["MassDensity"] = "Cannot give MassDensity with NumberDensity set";
    }
    bool canCalculateMassDensity =
        ((!isEmpty(params.mass)) && (!isEmpty(params.volume)));
    if (canCalculateMassDensity) {
      result["MassDensity"] = "Cannot give MassDensity with NumberDensity set";
    }
  }
  return result;
}

/**
 * Set the parameters to build the material to the builder,
 * taking into account which values were and weren't set.
 *
 * @param params A struct containing all the parameters to be set.
 */
void ReadMaterial::setMaterialParameters(const MaterialParameters &params) {
  setMaterial(params.chemicalSymbol, params.atomicNumber, params.massNumber);

  // calculate the mass density if it wasn't provided
  double massDensity = params.massDensity;
  if (isEmpty(massDensity)) {
    if (!(isEmpty(params.mass) || isEmpty(params.volume)))
      massDensity = params.mass / params.volume;
  }

  setNumberDensity(massDensity, params.numberDensity, params.numberDensityUnit,
                   params.zParameter, params.unitCellVolume);
  setScatteringInfo(params.coherentXSection, params.incoherentXSection,
                    params.attenuationXSection, params.scatteringXSection,
                    params.attenuationProfileFileName);
}

/**
 * Construct the material,
 *
 *  @returns A unique pointer to the newly made material
 */
std::unique_ptr<Kernel::Material> ReadMaterial::buildMaterial() {
  return std::make_unique<Kernel::Material>(builder.build());
}

void ReadMaterial::setMaterial(const std::string &chemicalSymbol,
                               const int atomicNumber, const int massNumber) {
  if (!chemicalSymbol.empty()) {
    builder.setFormula(chemicalSymbol);
  } else if (atomicNumber != 0) {
    builder.setAtomicNumber(atomicNumber);
    builder.setMassNumber(massNumber);
  }
}

void ReadMaterial::setNumberDensity(
    const double rho_m, const double rho,
    Kernel::MaterialBuilder::NumberDensityUnit rhoUnit, const double zParameter,
    const double unitCellVolume) {
  if (!isEmpty(rho_m))
    builder.setMassDensity(rho_m);
  if (isEmpty(rho)) {
    if (!isEmpty(zParameter)) {
      builder.setZParameter(zParameter);
      builder.setUnitCellVolume(unitCellVolume);
    }
  } else {
    builder.setNumberDensity(rho);
    builder.setNumberDensityUnit(rhoUnit);
  }
}

void ReadMaterial::setScatteringInfo(double coherentXSection,
                                     double incoherentXSection,
                                     double attenuationXSection,
                                     double scatteringXSection,
                                     std::string attenuationProfileFileName) {
  builder.setCoherentXSection(coherentXSection);       // in barns
  builder.setIncoherentXSection(incoherentXSection);   // in barns
  builder.setAbsorptionXSection(attenuationXSection);  // in barns
  builder.setTotalScatterXSection(scatteringXSection); // in barns
  builder.setAttenuationProfileFilename(attenuationProfileFileName);
}

bool ReadMaterial::isEmpty(const double toCheck) {
  return std::abs((toCheck - EMPTY_DBL()) / (EMPTY_DBL())) < 1e-8;
}
} // namespace DataHandling
} // namespace Mantid
