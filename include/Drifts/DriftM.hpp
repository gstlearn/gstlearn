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
#pragma once

#include "gstlearn_export.hpp"
#include "Drifts/ADriftElem.hpp"

/**
 * Monomial drift term
 */
class GSTLEARN_EXPORT DriftM : public ADriftElem
{
public:
  DriftM(const VectorInt &powers = VectorInt(),
         const CovContext &ctxt = CovContext());
  DriftM(const DriftM &r);
  DriftM& operator= (const DriftM &r);
  virtual ~DriftM();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftM)

  String getDriftName() const override;
  int    getOrderIRF()  const override;
  int    getOrderIRFIdim(int idim) const override;
  int    getNDim()      const override;
  double eval(const Db* db, int iech) const override;
  VectorInt getPowers() const override { return _monomialPower; }

private:
  VectorInt _monomialPower;
};
