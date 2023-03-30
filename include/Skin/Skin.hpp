/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Skin/ISkinFunctions.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/VectorNumT.hpp"

class GSTLEARN_EXPORT Skin
{
public:
  Skin(const ISkinFunctions* skf, DbGrid* dbgrid = nullptr);
  Skin(const Skin& r);
  Skin& operator=(const Skin& r);
  virtual ~Skin();

  int gridShift(int lec, int dir);
  int init(bool verbose = false);
  int remains(bool verbose = false);
  void getNext(int *rank, int *ipos);
  int unstack(int rank0, int ipos0);
  void skinPrint();

private:
  double _getWeight(int ipos, int idir);
  int    _gridShift(const VectorInt& indg0, int dir);
  void   _cellDelete(int rank);
  int    _cellAlreadyFilled(int ipos);
  void   _cellModify(int rank, double energy);
  int    _cellAdd(int ipos, double energy);
  int    _getNDim() const;

private:
  const ISkinFunctions* _skf;
  DbGrid* _dbgrid;
  int _nxyz;
  int _nval;
  int _date;
  int _nvalMax;
  double  _total;
  double  _totalMax;
  VectorInt _address;
  VectorDouble _energy;
};
