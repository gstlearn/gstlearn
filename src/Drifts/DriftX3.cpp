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
#include "Drifts/DriftX3.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"
#include "Drifts/EDrift.hpp"

DriftX3::DriftX3(const CovContext& ctxt)
    : ADriftElem(EDrift::X3 ,ctxt)
{
}

DriftX3::DriftX3(const DriftX3 &r)
    : ADriftElem(r)
{
}

DriftX3& DriftX3::operator=(const DriftX3 &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftX3::~DriftX3()
{
}

double DriftX3::eval(const Db* db, int iech) const
{
  double valx = db->getCoordinate(iech,0);
  return valx * valx * valx;
}

ICloneable* DriftX3::clone() const
{
  return new DriftX3(*this);
}
