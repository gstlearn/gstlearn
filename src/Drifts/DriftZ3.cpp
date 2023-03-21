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

#include "Drifts/DriftZ3.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"

DriftZ3::DriftZ3(const CovContext& ctxt)
    : ADriftElem(EDrift::Z3, ctxt)
{
}

DriftZ3::DriftZ3(const DriftZ3 &r)
    : ADriftElem(r)
{
}

DriftZ3& DriftZ3::operator=(const DriftZ3 &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftZ3::~DriftZ3()
{
}

double DriftZ3::eval(const Db* db, int iech) const
{
  double valz = db->getCoordinate(iech,2);
  return valz * valz * valz;
}

