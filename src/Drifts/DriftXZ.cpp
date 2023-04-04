/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Enum/EDrift.hpp"

#include "Drifts/DriftXZ.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"

DriftXZ::DriftXZ(const CovContext& ctxt)
    : ADriftElem(EDrift::XZ, ctxt)
{
}

DriftXZ::DriftXZ(const DriftXZ &r)
    : ADriftElem(r)
{
}

DriftXZ& DriftXZ::operator=(const DriftXZ &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftXZ::~DriftXZ()
{
}

double DriftXZ::eval(const Db* db, int iech) const
{
  double valx = db->getCoordinate(iech,0);
  double valz = db->getCoordinate(iech,2);
  return valx * valz;
}

