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

#include "Skin/ISkinFunctions.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/VectorNumT.hpp"

namespace gstlrn
{
class GSTLEARN_EXPORT Skin
{
public:
  Skin(const ISkinFunctions* skf, DbGrid* dbgrid = nullptr);
  Skin(const Skin& r);
  Skin& operator=(const Skin& r);
  virtual ~Skin();

  Id gridShift(Id lec, Id dir);
  Id init(bool verbose = false);
  Id remains(bool verbose = false);
  void getNext(Id *rank, Id *ipos);
  Id unstack(Id rank0, Id ipos0);
  void skinPrint() const;

private:
  double _getWeight(Id ipos, Id idir);
  Id    _gridShift(const VectorInt& indg0, Id dir);
  void   _cellDelete(Id rank);
  Id    _cellAlreadyFilled(Id ipos);
  void   _cellModify(Id rank, double energy);
  Id    _cellAdd(Id ipos, double energy);
  Id    _getNDim() const;

private:
  const ISkinFunctions* _skf;
  DbGrid* _dbgrid;
  Id _nxyz;
  Id _nval;
  Id _date;
  Id _nvalMax;
  double  _total;
  double  _totalMax;
  VectorInt _address;
  VectorDouble _energy;
};
}