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
#include "Drifts/DriftXY2.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"
#include "Drifts/EDrift.hpp"

DriftXY2::DriftXY2(const CovContext& ctxt)
    : ADriftElem(EDrift::XY2, ctxt)
{
}

DriftXY2::DriftXY2(const DriftXY2 &r)
    : ADriftElem(r)
{
}

DriftXY2& DriftXY2::operator=(const DriftXY2 &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftXY2::~DriftXY2()
{
}

double DriftXY2::eval(const Db* db, int iech) const
{
  double valx = db->getCoordinate(iech,0);
  double valy = db->getCoordinate(iech,1);
  return valx * valy * valy;
}

ICloneable* DriftXY2::clone() const
{
  return new DriftXY2(*this);
}
