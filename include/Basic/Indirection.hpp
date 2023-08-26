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
#include "geoslib_define.h"

#include "Basic/VectorNumT.hpp"
#include "Basic/AStringable.hpp"

#include <map>


class GSTLEARN_EXPORT Indirection: public AStringable
{
public:
  Indirection(int mode = 0);
  ~Indirection();
  Indirection(const Indirection &m);
  Indirection& operator=(const Indirection &m);

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void buildFromSel(const VectorDouble& sel);
  void buildFromRankRInA(const VectorInt& rels, int nabs);
  void buildFromMap(const std::map<int, int> &map, int nabs);
  int  getAToR(int iabs) const;
  int  getRToA(int irel) const;
  int  getAbsSize() const { return _nabs; }
  int  getRelSize() const { return _nrel; }
  void setMode(int mode);

  bool isDefined() const { return _defined; }

  VectorInt getRelRanks() const { return _vecRToA; }
  int getMode() const { return _mode; }

private:
  void _resetMap();
  int  _getMapAToR(int iabs) const;
  int  _getArrayAToR(int iabs) const;
  bool _isValidAbs(int iabs) const;
  bool _isValidRel(int irel) const;

private:
  bool _defined;
  int  _mode;     // 0 by array; 1 by MAP
  int  _nabs;     // Number of absolute elements
  int  _nrel;     // Number of relative elements
  VectorInt _vecRToA;
  VectorInt _vecAToR;
  std::map<int, int> _mapAToR;
};
