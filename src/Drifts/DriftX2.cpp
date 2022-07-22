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
#include "Drifts/DriftX2.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"
#include "Drifts/EDrift.hpp"

DriftX2::DriftX2(const CovContext& ctxt)
    : ADriftElem(EDrift::X2, ctxt)
{
}

DriftX2::DriftX2(const DriftX2 &r)
    : ADriftElem(r)
{
}

DriftX2& DriftX2::operator=(const DriftX2 &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftX2::~DriftX2()
{
}

double DriftX2::eval(const Db* db, int iech) const
{
  double valx = db->getCoordinate(iech,0);
  return valx * valx;
}

ICloneable* DriftX2::clone() const
{
  return new DriftX2(*this);
}
