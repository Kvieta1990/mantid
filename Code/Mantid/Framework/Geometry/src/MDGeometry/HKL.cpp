#include "MantidGeometry/MDGeometry/HKL.h"

namespace {

void checkUnitCompatibility(Mantid::Kernel::MDUnit const * const unit){
    if(!unit->isQUnit()){
        throw std::invalid_argument("HKL unit must be a QUnit");
    }
}

}

namespace Mantid {
namespace Geometry {

const std::string HKL::HKLName = "HKL";

/**
 * Constructor
 * @param unit : Unit to use
 */
HKL::HKL(std::unique_ptr<Kernel::MDUnit>& unit){
    checkUnitCompatibility(unit.get());
    // Only change ownership once we are happy. Gives exception safety for input unit.
    m_unit.swap(unit);
}

/**
 * Constructor
 * @param unit : Unit to use
 */
HKL::HKL(Kernel::MDUnit* unit) : m_unit(unit){
    checkUnitCompatibility(unit);
}

/**
 * Assignment
 * @param other : Other unit
 * @return assigned unit.
 */
HKL::HKL(const HKL &other) : m_unit(other.getMDUnit().clone())
{

}

HKL& HKL::operator=(const HKL& other){
    if(this == &other){
        this->m_unit = std::unique_ptr<Mantid::Kernel::MDUnit>(other.getMDUnit().clone());
    }
    return *this;
}

//----------------------------------------------------------------------------------------------
/** Destructor
 */
HKL::~HKL() {}

Kernel::UnitLabel HKL::getUnitLabel() const
{
    return m_unit->getUnitLabel();
}

const Kernel::MDUnit& HKL::getMDUnit() const
{
    return *m_unit;
}

bool HKL::canConvertTo(const Kernel::MDUnit &otherUnit) const
{
    return this->m_unit->canConvertTo(otherUnit);
}

std::string HKL::name() const
{
    return HKLName;
}

} // namespace Geometry
} // namespace Mantid
