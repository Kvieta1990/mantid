#ifndef MANTID_API_REFAXIS_H_
#define MANTID_API_REFAXIS_H_

//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidAPI/Axis.h"

namespace Mantid
{
namespace API
{
/** A class to represent the axis of a 2D (or more) workspace where the value at
    a given point on the axis varies along the other dimension.
    This class does not hold the axis values itself; they are held by the X vectors
    of the workspace itself.
    An axis of this kind is always of numeric type.

    @author Russell Taylor, Tessella Support Services plc
    @date 18/05/2008
    
    Copyright &copy; 2008 STFC Rutherford Appleton Laboratory

    This file is part of Mantid.
  
    Mantid is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.
  
    Mantid is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
  
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
    File change history is stored at: <https://svn.mantidproject.org/mantid/trunk/Code/Mantid>.
    Code Documentation is available at: <http://doxygen.mantidproject.org>
*/
class DLLExport RefAxis : public Axis
{
public:
	explicit RefAxis(const Workspace* const parentWorkspace);
	virtual ~RefAxis();

	Axis* clone(const Workspace* const parentWorkspace);
	
  virtual const double operator()(const int index, const int verticalIndex) const;

private:
  RefAxis(const RefAxis& right, const Workspace* const parentWorkspace);
  /// Private, undefined 'regular' copy constructor
  RefAxis(const RefAxis&);  
  /// Private, undefined copy assignment operator
  const RefAxis& operator=(const RefAxis&);

  /// A pointer to the workspace holding the axis
  const Workspace* const m_parentWS;
};

} // namespace API
} // namespace Mantid

#endif /*MANTID_API_REFAXIS_H_*/
