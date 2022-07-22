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
#include "Drifts/DriftY2.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"
#include "Drifts/EDrift.hpp"

DriftY2::DriftY2(const CovContext& ctxt)
    : ADriftElem(EDrift::Y2, ctxt)
{
}

DriftY2::DriftY2(const DriftY2 &r)
    : ADriftElem(r)
{
}

DriftY2& DriftY2::operator=(const DriftY2 &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftY2::~DriftY2()
{
}

double DriftY2::eval(const Db* db, int iech) const
{
  double valy = db->getCoordinate(iech,1);
  return valy * valy;
}

ICloneable* DriftY2::clone() const
{
  return new DriftY2(*this);
}
