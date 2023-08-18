/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/EConsElem.hpp"

#include "Basic/VectorNumT.hpp"
#include "Basic/AStringable.hpp"

class GSTLEARN_EXPORT ElemNostat: public AStringable
{
public:
  ElemNostat();
  ElemNostat(const ElemNostat &m);
  ElemNostat& operator= (const ElemNostat &m);
  virtual ~ElemNostat();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void init(const EConsElem& loctype, int rank_grf, int rank_str, int rank_v1, int rank_v2);

  const EConsElem& getLocType() const { return _locType; }
  int getRankGrf() const { return _rankGRF; }
  int getRankStr() const { return _rankStr; }
  int getRankV1() const { return _rankV1; }
  int getRankV2() const { return _rankV2; }
  double getVal1() const { return _val1; }
  double getVal2() const { return _val2; }
  void setVal1(double val1) { _val1 = val1; }
  void setVal2(double val2) { _val2 = val2; }

private:
  EConsElem _locType; /* Type of parameter (by its locator type) */
  int _rankGRF; /* Rank of the GRF */
  int _rankStr; /* Rank of the basic structure (from 0) */
  int _rankV1; /* Rank of the first variable (from 0) */
  int _rankV2; /* Rank of the second variable (from 0) */
  double _val1; /* Value at the first point */
  double _val2; /* Value at the second point */
};
