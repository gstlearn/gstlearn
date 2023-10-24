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

#include "LinearOp/CGParam.hpp"
#include "Basic/VectorNumT.hpp"

class ALinearOp;

class GSTLEARN_EXPORT SPDEParam {

public:
  SPDEParam(int refineK = 11,
            int refineS = 18,
            int border = 8,
            double epsNugget = EPSILON2,
            const CGParam cgparams = CGParam());
  SPDEParam(const SPDEParam &m);
  SPDEParam& operator=(const SPDEParam &m);
  virtual ~SPDEParam();

  int getBorder() const { return _border; }
  const CGParam getCGparams() const { return _CGparams; }
  double getEpsNugget() const { return _epsNugget; }
  int getRefineK() const { return _refineK; }
  int getRefineS() const { return _refineS; }

  void setBorder(int border) { _border = border; }
  void setCGparams(const CGParam &CGparams) { _CGparams = CGparams; }
  void setEpsNugget(double epsNugget) { _epsNugget = epsNugget; }
  void setRefineK(int refineK) { _refineK = refineK; }
  void setRefineS(int refineS) { _refineS = refineS; }

private:
  int _refineK;
  int _refineS;
  int _border;
  double _epsNugget;

  CGParam _CGparams;
};
