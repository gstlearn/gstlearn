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

#include "Drifts/DriftXY2.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"

DriftXY2::DriftXY2(const CovContext& ctxt)
    : ADriftElem(EDrift::XY2, ctxt)
{
}

DriftXY2::DriftXY2(const DriftXY2 &r)
    : ADriftElem(r)
{
}

DriftXY2& DriftXY2::operator=(const DriftXY2 &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftXY2::~DriftXY2()
{
}

double DriftXY2::eval(const Db* db, int iech) const
{
  double valx = db->getCoordinate(iech,0);
  double valy = db->getCoordinate(iech,1);
  return valx * valy * valy;
}

