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

#include "Drifts/DriftZ.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"

DriftZ::DriftZ(const CovContext& ctxt)
    : ADriftElem(EDrift::Z, ctxt)
{
}

DriftZ::DriftZ(const DriftZ &r)
    : ADriftElem(r)
{
}

DriftZ& DriftZ::operator=(const DriftZ &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftZ::~DriftZ()
{
}

double DriftZ::eval(const Db* db, int iech) const
{
  return db->getCoordinate(iech,2);
}

