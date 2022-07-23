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
#include "Drifts/DriftXY.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"
#include "Drifts/EDrift.hpp"

DriftXY::DriftXY(const CovContext& ctxt)
    : ADriftElem(EDrift::XY, ctxt)
{
}

DriftXY::DriftXY(const DriftXY &r)
    : ADriftElem(r)
{
}

DriftXY& DriftXY::operator=(const DriftXY &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftXY::~DriftXY()
{
}

double DriftXY::eval(const Db* db, int iech) const
{
  double valx = db->getCoordinate(iech,0);
  double valy = db->getCoordinate(iech,1);
  return valx * valy;
}

ICloneable* DriftXY::clone() const
{
  return new DriftXY(*this);
}
