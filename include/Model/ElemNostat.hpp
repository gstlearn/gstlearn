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

#include "Enum/EConsElem.hpp"

#include "Basic/AStringable.hpp"

namespace gstlrn
{
class GSTLEARN_EXPORT ElemNostat: public AStringable
{
public:
  ElemNostat();
  ElemNostat(const ElemNostat &m);
  ElemNostat& operator= (const ElemNostat &m);
  virtual ~ElemNostat();

  String toString(const AStringFormat* strfmt = nullptr) const override;

  void init(const EConsElem& loctype, Id rank_grf, Id rank_str, Id rank_v1, Id rank_v2);

  const EConsElem& getLocType() const { return _locType; }
  Id getRankGrf() const { return _rankGRF; }
  Id getRankStr() const { return _rankStr; }
  Id getRankV1() const { return _rankV1; }
  Id getRankV2() const { return _rankV2; }
  double getVal1() const { return _val1; }
  double getVal2() const { return _val2; }
  void setVal1(double val1) { _val1 = val1; }
  void setVal2(double val2) { _val2 = val2; }

private:
  EConsElem _locType; /* Type of parameter (by its locator type) */
  Id _rankGRF; /* Rank of the GRF */
  Id _rankStr; /* Rank of the basic structure (from 0) */
  Id _rankV1; /* Rank of the first variable (from 0) */
  Id _rankV2; /* Rank of the second variable (from 0) */
  double _val1; /* Value at the first point */
  double _val2; /* Value at the second point */
};
}