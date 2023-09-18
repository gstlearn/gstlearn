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
         double coeff0 = 1.,
         const VectorDouble &coeffs = VectorDouble(),
         const CovContext &ctxt = CovContext());
  DriftM(const DriftM &r);
  DriftM& operator= (const DriftM &r);
  virtual ~DriftM();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftM)

  String getDriftName() const override;
  int    getOrderIRF()  const override;
  int    getNDim()      const override;
  double eval(const Db* db, int iech) const override;

private:
  VectorInt    _monomialPower;
  double       _coeff0;
  VectorDouble _monomialCoeffs;
};
