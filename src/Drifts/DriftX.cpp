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
#include "Drifts/DriftX.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"
#include "geoslib_enum.h"

DriftX::DriftX(const CovContext& ctxt)
    : ADriftElem(DRIFT_X, ctxt)
{
}

DriftX::DriftX(const DriftX &r)
    : ADriftElem(r)
{
}

DriftX& DriftX::operator=(const DriftX &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftX::~DriftX()
{
}

double DriftX::eval(const Db* db, int iech) const
{
  return db->getCoordinate(iech,0);
}

IClonable* DriftX::clone() const
{
  return new DriftX(*this);
}
