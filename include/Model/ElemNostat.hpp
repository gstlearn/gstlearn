/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"

class ElemNostat: public AStringable
{
public:
  ElemNostat();
  ElemNostat(const ElemNostat &m);
  ElemNostat& operator= (const ElemNostat &m);
  virtual ~ElemNostat();

  virtual String toString(int level) const override;

  void init(int loctype, int rank_grf, int rank_str, int rank_v1, int rank_v2);
  int getLocType() const
  {
    return _locType;
  }

  int getRankGrf() const
  {
    return _rankGRF;
  }

  int getRankStr() const
  {
    return _rankStr;
  }

  int getRankV1() const
  {
    return _rankV1;
  }

  int getRankV2() const
  {
    return _rankV2;
  }

  double getVal1() const
  {
    return _val1;
  }

  double getVal2() const
  {
    return _val2;
  }
  void setVal1(double val1) { _val1 = val1; }
  void setVal2(double val2) { _val2 = val2; }

private:
  int _locType; /* Type of parameter (by its locator type) */
  int _rankGRF; /* Rank of the GRF */
  int _rankStr; /* Rank of the basic structure (from 0) */
  int _rankV1; /* Rank of the first variable (from 0) */
  int _rankV2; /* Rank of the second variable (from 0) */
  double _val1; /* Value at the first point */
  double _val2; /* Value at the second point */
};
