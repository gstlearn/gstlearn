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

#include "Drifts/DriftXY.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"

DriftXY::DriftXY(const CovContext& ctxt)
    : ADriftElem(EDrift::XY, ctxt)
{
}

DriftXY::DriftXY(const DriftXY &r)
    : ADriftElem(r)
{
}

DriftXY& DriftXY::operator=(const DriftXY &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftXY::~DriftXY()
{
}

double DriftXY::eval(const Db* db, int iech) const
{
  double valx = db->getCoordinate(iech,0);
  double valy = db->getCoordinate(iech,1);
  return valx * valy;
}

