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

namespace gstlrn
{
class GSTLEARN_EXPORT Indirection: public AStringable
{
public:
  Indirection(Id mode = 0);
  ~Indirection();
  Indirection(const Indirection &m);
  Indirection& operator=(const Indirection &m);

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  void buildFromSel(const VectorDouble& sel);
  void buildFromRankRInA(const VectorInt& rels, Id nabs);
  void buildFromMap(const std::map<Id, Id> &map, Id nabs);
  Id  getAToR(Id iabs) const;
  Id  getRToA(Id irel) const;
  Id  getAbsSize() const { return _nabs; }
  Id  getRelSize() const { return _nrel; }
  void setMode(Id mode);

  bool isDefined() const { return _defined; }

  VectorInt getRelRanks() const { return _vecRToA; }
  Id getMode() const { return _mode; }

private:
  void _resetMap();
  Id  _getMapAToR(Id iabs) const;
  Id  _getArrayAToR(Id iabs) const;
  bool _isValidAbs(Id iabs) const;
  bool _isValidRel(Id irel) const;

private:
  bool _defined;
  Id  _mode;     // 0 by array; 1 by MAP
  Id  _nabs;     // Number of absolute elements
  Id  _nrel;     // Number of relative elements
  VectorInt _vecRToA;
  VectorInt _vecAToR;
  std::map<Id, Id> _mapAToR;
};
}