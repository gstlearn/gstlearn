/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Covariances/ACor.hpp"
#include "Covariances/CovCalcMode.hpp"

#include "Db/Db.hpp"

#include "Space/ASpace.hpp"
#include "geoslib_define.h"

#include <math.h>

ACor::ACor(const ASpace *space,int nvar)
    : ASpaceObject(space),
      _nvar(nvar)
{
}

ACor::ACor(const ACor &r)
    : ASpaceObject(r),
      _nvar(r._nvar)
{
}

ACor& ACor::operator=(const ACor &r)
{
  if (this != &r)
  {
    ASpaceObject::operator=(r);
    _nvar = r._nvar;
  }
  return *this;
}

ACor::~ACor()
{
}


double ACor::eval0(int ivar,
                   int jvar,
                   const CovCalcMode* mode) const
{
  DECLARE_UNUSED(ivar,jvar,mode)
  return 1.;
}
