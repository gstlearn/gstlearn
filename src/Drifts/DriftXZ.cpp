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
#include "Drifts/DriftXZ.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"
#include "Drifts/EDrift.hpp"

DriftXZ::DriftXZ(const CovContext& ctxt)
    : ADriftElem(EDrift::XZ, ctxt)
{
}

DriftXZ::DriftXZ(const DriftXZ &r)
    : ADriftElem(r)
{
}

DriftXZ& DriftXZ::operator=(const DriftXZ &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftXZ::~DriftXZ()
{
}

double DriftXZ::eval(const Db* db, int iech) const
{
  double valx = db->getCoordinate(iech,0);
  double valz = db->getCoordinate(iech,2);
  return valx * valz;
}

IClonable* DriftXZ::clone() const
{
  return new DriftXZ(*this);
}
