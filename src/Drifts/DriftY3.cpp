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
#include "Drifts/DriftY3.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"
#include "Drifts/EDrift.hpp"

DriftY3::DriftY3(const CovContext& ctxt)
    : ADriftElem(EDrift::Y3, ctxt)
{
}

DriftY3::DriftY3(const DriftY3 &r)
    : ADriftElem(r)
{
}

DriftY3& DriftY3::operator=(const DriftY3 &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftY3::~DriftY3()
{
}

double DriftY3::eval(const Db* db, int iech) const
{
  double valy = db->getCoordinate(iech,1);
  return valy * valy * valy;
}

IClonable* DriftY3::clone() const
{
  return new DriftY3(*this);
}
