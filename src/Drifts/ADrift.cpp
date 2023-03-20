/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Drifts/ADrift.hpp"

#include "Db/Db.hpp"
#include "Space/SpacePoint.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorNumT.hpp"

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
