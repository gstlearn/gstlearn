/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Enum/EDrift.hpp"

#include "Drifts/DriftX3.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"

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

