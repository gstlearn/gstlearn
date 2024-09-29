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
#pragma once

#include "Covariances/CovGneiting.hpp"
#include "Covariances/ACov.hpp"
#include "Space/SpacePoint.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "geoslib_define.h"

CovGneiting::CovGneiting()
{

}

CovGneiting::CovGneiting(const CovGneiting& r)
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

double CovGneiting::eval0(int ivar,
                          int jvar,
                          const CovCalcMode* mode) const
{
    DECLARE_UNUSED(ivar,jvar,mode)
    return 1.;
}

void CovGneiting::eval0MatInPlace(MatrixSquareGeneral &mat,
                                 const CovCalcMode *mode) const 
{
    DECLARE_UNUSED(mat,mode)
}

CovGneiting::~CovGneiting()
{

}


double CovGneiting::eval(const SpacePoint& p1,
                    const SpacePoint& p2,
                    int ivar = 0,
                    int jvar = 0,
                    const CovCalcMode* mode = nullptr) const
{
    constvect p1;
    constvect p2;
    p1 = constvect(p1.getCoords(0));
    p2 = constvect(p2.getCoords(0));
    _covS->eval(p1, const SpacePoint &p2)
}
