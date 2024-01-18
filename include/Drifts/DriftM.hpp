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
#include "Drifts/ADrift.hpp"

/**
 * Monomial drift term.
 * Examples:
 * - driftM() is the Universality condition
 * - driftM([2]) (where [] stands for a vector of integers) stands for x_1^2
 * - driftM([2,3]) stands for x_1^2 * x_2^3
 * Note: the size of the vector (when defined) must be smaller or equal to the space dimension
 */
class GSTLEARN_EXPORT DriftM : public ADrift
{
public:
  DriftM(const VectorInt &powers = VectorInt());
  DriftM(const DriftM &r);
  DriftM& operator= (const DriftM &r);
  virtual ~DriftM();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftM)

  String getDriftName() const override;
  int    getOrderIRF()  const override;
  int    getOrderIRFIdim(int idim) const override;
  int    getDriftNDimMax()      const override;
  double eval(const Db* db, int iech) const override;
  VectorInt getPowers() const override { return _monomialPower; }

  static DriftM* createByIdentifier(const String &driftname);

private:
  VectorInt _monomialPower;
};
