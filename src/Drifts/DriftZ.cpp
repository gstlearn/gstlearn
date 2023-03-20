/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
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

