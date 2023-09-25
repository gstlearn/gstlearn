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

class GSTLEARN_EXPORT DriftF : public ADriftElem
{
public:
  DriftF(int rank_fex = 0, const CovContext& ctxt = CovContext());
  DriftF(const DriftF &r);
  DriftF& operator= (const DriftF &r);
  virtual ~DriftF();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftF)

  String getDriftName() const override;
  int    getOrderIRF() const override { return -1; }
  int    getOrderIRFIdim(int idim) const override { return -1; }
  bool   isDriftExternal() const override { return true; }
  double eval(const Db* db, int iech) const override;
  int    getRankFex() const override { return _rankFex; }

  static DriftF* createByIdentifier(const String &driftname,
                                    const CovContext &ctxt = CovContext());

private:
  int _rankFex;       /* Rank of the external drift */
};
