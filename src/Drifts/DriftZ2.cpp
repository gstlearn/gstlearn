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
#include "Enum/EDrift.hpp"

#include "Drifts/DriftZ2.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"

DriftZ2::DriftZ2(const CovContext& ctxt)
    : ADriftElem(EDrift::Z2, ctxt)
{
}

DriftZ2::DriftZ2(const DriftZ2 &r)
    : ADriftElem(r)
{
}

DriftZ2& DriftZ2::operator=(const DriftZ2 &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftZ2::~DriftZ2()
{
}

double DriftZ2::eval(const Db* db, int iech) const
{
  double valz = db->getCoordinate(iech,2);
  return valz * valz;
}

