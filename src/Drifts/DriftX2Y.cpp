/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Enum/EDrift.hpp"

#include "Drifts/DriftX2Y.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"

DriftX2Y::DriftX2Y(const CovContext& ctxt)
    : ADriftElem(EDrift::X2Y, ctxt)
{
}

DriftX2Y::DriftX2Y(const DriftX2Y &r)
    : ADriftElem(r)
{
}

DriftX2Y& DriftX2Y::operator=(const DriftX2Y &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftX2Y::~DriftX2Y()
{
}

double DriftX2Y::eval(const Db* db, int iech) const
{
  double valx = db->getCoordinate(iech,0);
  double valy = db->getCoordinate(iech,1);
  return valx * valx * valy;
}

