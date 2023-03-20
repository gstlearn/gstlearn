/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Enum/EDrift.hpp"

#include "Drifts/DriftX.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"

DriftX::DriftX(const CovContext& ctxt)
    : ADriftElem(EDrift::X, ctxt)
{
}

DriftX::DriftX(const DriftX &r)
    : ADriftElem(r)
{
}

DriftX& DriftX::operator=(const DriftX &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftX::~DriftX()
{
}

double DriftX::eval(const Db* db, int iech) const
{
  return db->getCoordinate(iech,0);
}

