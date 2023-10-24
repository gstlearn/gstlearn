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
 * Monomial drift term
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
