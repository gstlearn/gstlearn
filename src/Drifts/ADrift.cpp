/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "Drifts/ADrift.hpp"

#include "Db/Db.hpp"
#include "Space/SpacePoint.hpp"
#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"

ADrift::ADrift(const ASpace* space)
: ASpaceObject(space)
{
}

ADrift::ADrift(const ADrift &r)
: ASpaceObject(r)
{
}

ADrift& ADrift::operator=(const ADrift &r)
{
  if (this != &r)
  {
    ASpaceObject::operator=(r);
  }
  return *this;
}

ADrift::~ADrift()
{
}
