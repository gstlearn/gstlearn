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
#include "Neigh/ANeighParam.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"

class Db;

class GSTLEARN_EXPORT ANeigh
{
public:
  ANeigh(const Db *dbin = nullptr,
         const ANeighParam *neighparam = nullptr,
         const Db *dbout = nullptr);
  ANeigh(const ANeigh& r);
  ANeigh& operator=(const ANeigh& r);
  virtual ~ANeigh();

  virtual int initialize(const Db *dbin,
                         const ANeighParam *neighparam,
                         const Db *dbout = nullptr);
  virtual VectorInt getNeigh(int iech_out) = 0;
  virtual bool hasChanged(int iech_out) const { return true; }

  VectorInt select(int iech_out);
  VectorInt selectBySP(const SpacePoint& pt_out);
  bool isUnchanged() const { return _flagIsUnchanged; }
  void setIsChanged();

protected:
  bool _isNbghMemoEmpty() const { return _nbghMemo.empty(); }

private:
  bool _isSameTarget(int iech_out);
  void _checkUnchanged(int iech_out, const VectorInt &ranks);
  void _updateColCok(VectorInt& ranks, int iech_out);

protected:
  const Db* _dbin; // compulsory
  const Db* _dbout; // optional
  const DbGrid* _dbgrid; // optional (is similar to dbout and defined only for grid
  const ANeighParam* _neighParam;
  VectorInt _rankColCok;
  int _iechMemo;

private:
  bool _flagIsUnchanged;
  mutable SpacePoint _ptOut;
  mutable VectorInt  _nbghMemo;
};
