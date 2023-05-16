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
#pragma once

#include "gstlearn_export.hpp"
#include "Neigh/ANeigh.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Db/Db.hpp"
#include "geoslib_define.h"

class GSTLEARN_EXPORT NeighWork:  public ANeigh
{
public:
  NeighWork(const Db *dbin = nullptr,
            const ANeighParam *neighparam = nullptr,
            const Db *dbout = nullptr);
  NeighWork(const NeighWork& r);
  NeighWork& operator=(const NeighWork& r);
  virtual ~NeighWork();

  int initialize(const Db *dbin,
                 const ANeighParam *neighparam,
                 const Db *dbout = nullptr) override;
  bool hasChanged(int iech_out) const override;
  VectorInt getNeigh(int iech_out) override;

  void clear();
  VectorDouble summary(int iech_out);
  void setFlagSimu(bool flagSimu) { _flagSimu = flagSimu; }
  void setRankColCok(const VectorInt &rankColCok) { _rankColCok = rankColCok; }

private:
  void _unique(int iech_out, VectorInt& ranks);
  void _bench(int iech_out, VectorInt& ranks);
  int  _moving(int iech_out, VectorInt& ranks, double eps = EPSILON9);
  bool _discardUndefined(int iech);
  int  _xvalid(int iech_in, int iech_out, double eps = EPSILON9);
  int  _movingSectorDefine(double dx, double dy);
  void _movingSectorNsmax(int nsel, VectorInt& ranks);
  void _movingSelect(int nsel, VectorInt& ranks);
  void _display(const VectorInt& ranks);
  double _movingDist(int iech_in, int iech_out);
  bool _belongsToCell(int iech, int iech_out);
  void _clearMemoryNeigh();
  void _clearMemoryMoving();
  bool _isSameTargetBench(int iech_out) const;
  bool _hiddenByFault(int iech, int iech_out) const;

private:
  bool _flagSimu;
  mutable VectorInt    _movingInd;
  mutable VectorInt    _movingIsect;
  mutable VectorInt    _movingNsect;
  mutable VectorDouble _movingX1;
  mutable VectorDouble _movingX2;
  mutable VectorDouble _movingDst;
};
