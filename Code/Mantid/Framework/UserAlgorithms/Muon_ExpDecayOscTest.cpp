/*  This is another User test fitting function.

    For a general description of how to create a test
    fitting function see LorentzianTest.h and LorentzianTest.h

    This function was originally provided by the ISIS Muon group.
*/

//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "Muon_ExpDecayOscTest.h"
#include <cmath>

namespace Mantid
{
namespace CurveFitting
{

using namespace Kernel;
using namespace API;

DECLARE_FUNCTION(Muon_ExpDecayOscTest)

void Muon_ExpDecayOscTest::init()
{
  declareParameter("A", 0.2);
  declareParameter("lambda", 0.2);
  declareParameter("frequency", 0.5);
  declareParameter("phi", 0.0);
}


void Muon_ExpDecayOscTest::functionLocal(double* out, const double* xValues, const size_t nData)const
{
  const double& gA0 = getParameter("A");
  const double& gs = getParameter("lambda");
	const double& gf = getParameter("frequency");
	const double& gphi = getParameter("phi");
  
 
  for (int i = 0; i < nData; i++) 
  {
    if ( gs*xValues[i] < 20 )
      out[i] = gA0*exp(-gs*xValues[i])*cos(2*3.1415926536*gf*xValues[i]+gphi);
    else
      out[i] = 0.0;
  }
}
void Muon_ExpDecayOscTest::functionDerivLocal(API::Jacobian* out, const double* xValues, const size_t nData)
{
  calNumericalDeriv(out, xValues, nData);
}




} // namespace CurveFitting
} // namespace Mantid
