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

#include "Covariances/CovGneiting.hpp"
#include "Covariances/ACov.hpp"
#include "Space/SpacePoint.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "geoslib_define.h"

CovGneiting::CovGneiting()
{

}

CovGneiting::CovGneiting(const CovGneiting& r):
ACov(r)
{

}

CovGneiting& CovGneiting::operator=(const CovGneiting &r)
{
  if (this != &r)
  {
    ACov::operator =(r);
    _ctxt = r._ctxt;
  }
  return *this;
}


CovGneiting::~CovGneiting()
{

}


double CovGneiting::eval(const SpacePoint& p1,
                    const SpacePoint& p2,
                    int ivar,
                    int jvar,
                    const CovCalcMode* mode) const
{

  auto p1_0 = p1.projection(0);
  auto p2_0 = p2.projection(0);
  auto p1_1 = p1.projection(1);
  auto p2_1 = p2.projection(1);
  return _covS->eval(p1_0, p2_0, ivar, jvar, mode) * 
         _covTemp->eval(p1_1, p2_1, ivar, jvar, mode);

}
