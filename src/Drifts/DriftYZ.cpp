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

#include "Drifts/DriftYZ.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"

DriftYZ::DriftYZ(const CovContext& ctxt)
    : ADriftElem(EDrift::YZ, ctxt)
{
}

DriftYZ::DriftYZ(const DriftYZ &r)
    : ADriftElem(r)
{
}

DriftYZ& DriftYZ::operator=(const DriftYZ &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftYZ::~DriftYZ()
{
}

double DriftYZ::eval(const Db* db, int iech) const
{
  double valy = db->getCoordinate(iech,1);
  double valz = db->getCoordinate(iech,2);
  return valy * valz;
}

