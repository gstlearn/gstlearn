/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Enum/EDrift.hpp"

#include "Drifts/DriftY2.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"

DriftY2::DriftY2(const CovContext& ctxt)
    : ADriftElem(EDrift::Y2, ctxt)
{
}

DriftY2::DriftY2(const DriftY2 &r)
    : ADriftElem(r)
{
}

DriftY2& DriftY2::operator=(const DriftY2 &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftY2::~DriftY2()
{
}

double DriftY2::eval(const Db* db, int iech) const
{
  double valy = db->getCoordinate(iech,1);
  return valy * valy;
}

