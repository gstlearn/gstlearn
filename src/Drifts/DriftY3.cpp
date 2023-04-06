/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "Enum/EDrift.hpp"

#include "Drifts/DriftY3.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"

DriftY3::DriftY3(const CovContext& ctxt)
    : ADriftElem(EDrift::Y3, ctxt)
{
}

DriftY3::DriftY3(const DriftY3 &r)
    : ADriftElem(r)
{
}

DriftY3& DriftY3::operator=(const DriftY3 &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftY3::~DriftY3()
{
}

double DriftY3::eval(const Db* db, int iech) const
{
  double valy = db->getCoordinate(iech,1);
  return valy * valy * valy;
}

