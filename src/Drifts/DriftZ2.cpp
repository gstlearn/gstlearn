/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "Enum/EDrift.hpp"

#include "Drifts/DriftZ2.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"

DriftZ2::DriftZ2(const CovContext& ctxt)
    : ADriftElem(EDrift::Z2, ctxt)
{
}

DriftZ2::DriftZ2(const DriftZ2 &r)
    : ADriftElem(r)
{
}

DriftZ2& DriftZ2::operator=(const DriftZ2 &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftZ2::~DriftZ2()
{
}

double DriftZ2::eval(const Db* db, int iech) const
{
  double valz = db->getCoordinate(iech,2);
  return valz * valz;
}

