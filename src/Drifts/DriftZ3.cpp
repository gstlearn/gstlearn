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
#include "Drifts/DriftZ3.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"
#include "geoslib_enum.h"

DriftZ3::DriftZ3(const CovContext& ctxt)
    : ADriftElem(DRIFT_Z3, ctxt)
{
}

DriftZ3::DriftZ3(const DriftZ3 &r)
    : ADriftElem(r)
{
}

DriftZ3& DriftZ3::operator=(const DriftZ3 &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftZ3::~DriftZ3()
{
}

double DriftZ3::eval(const Db* db, int iech) const
{
  double valz = db->getCoordinate(iech,2);
  return valz * valz * valz;
}

IClonable* DriftZ3::clone() const
{
  return new DriftZ3(*this);
}
