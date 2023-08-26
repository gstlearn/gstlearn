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

#include "Drifts/DriftX2.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"

DriftX2::DriftX2(const CovContext& ctxt)
    : ADriftElem(EDrift::X2, ctxt)
{
}

DriftX2::DriftX2(const DriftX2 &r)
    : ADriftElem(r)
{
}

DriftX2& DriftX2::operator=(const DriftX2 &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftX2::~DriftX2()
{
}

double DriftX2::eval(const Db* db, int iech) const
{
  double valx = db->getCoordinate(iech,0);
  return valx * valx;
}

