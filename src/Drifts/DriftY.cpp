/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Enum/EDrift.hpp"

#include "Drifts/DriftY.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"

DriftY::DriftY(const CovContext& ctxt)
    : ADriftElem(EDrift::Y, ctxt)
{
}

DriftY::DriftY(const DriftY &r)
    : ADriftElem(r)
{
}

DriftY& DriftY::operator=(const DriftY &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftY::~DriftY()
{
}

double DriftY::eval(const Db* db, int iech) const
{
  return db->getCoordinate(iech,1);
}

