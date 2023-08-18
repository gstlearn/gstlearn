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

#include "Drifts/DriftF.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"

DriftF::DriftF(const CovContext& ctxt)
    : ADriftElem(EDrift::F, ctxt)
{
}

DriftF::DriftF(const DriftF &r)
    : ADriftElem(r)
{
}

DriftF& DriftF::operator=(const DriftF &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftF::~DriftF()
{
}

double DriftF::eval(const Db* db, int iech) const
{
  return db->getLocVariable(ELoc::F,iech,getRankFex());
}

