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
#include "Drifts/DriftYZ.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"
#include "Drifts/EDrift.hpp"

DriftYZ::DriftYZ(const CovContext& ctxt)
    : ADriftElem(EDrift::YZ, ctxt)
{
}

DriftYZ::DriftYZ(const DriftYZ &r)
    : ADriftElem(r)
{
}

DriftYZ& DriftYZ::operator=(const DriftYZ &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftYZ::~DriftYZ()
{
}

double DriftYZ::eval(const Db* db, int iech) const
{
  double valy = db->getCoordinate(iech,1);
  double valz = db->getCoordinate(iech,2);
  return valy * valz;
}

ICloneable* DriftYZ::clone() const
{
  return new DriftYZ(*this);
}
