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
#include "Drifts/DriftX2Y.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"
#include "geoslib_enum.h"

DriftX2Y::DriftX2Y(const CovContext& ctxt)
    : ADriftElem(DRIFT_X2Y, ctxt)
{
}

DriftX2Y::DriftX2Y(const DriftX2Y &r)
    : ADriftElem(r)
{
}

DriftX2Y& DriftX2Y::operator=(const DriftX2Y &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftX2Y::~DriftX2Y()
{
}

double DriftX2Y::eval(const Db* db, int iech) const
{
  double valx = db->getCoordinate(iech,0);
  double valy = db->getCoordinate(iech,1);
  return valx * valx * valy;
}

IClonable* DriftX2Y::clone() const
{
  return new DriftX2Y(*this);
}
