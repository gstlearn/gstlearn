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

namespace gstlrn
{
class GSTLEARN_EXPORT DriftF : public ADrift
{
public:
  DriftF(Id rank_fex = 0);
  DriftF(const DriftF &r);
  DriftF& operator= (const DriftF &r);
  virtual ~DriftF();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftF)

  String getDriftName() const override;
  Id    getOrderIRF() const override { return -1; }
  Id    getOrderIRFIdim(Id idim) const override {
    DECLARE_UNUSED(idim);
    return -1;
  }
  bool   isDriftExternal() const override { return true; }
  double eval(const Db* db, Id iech) const override;
  Id    getRankFex() const override { return _rankFex; }

  static DriftF* createByIdentifier(const String &driftname);

private:
  Id _rankFex;       /* Rank of the external drift */
};
}