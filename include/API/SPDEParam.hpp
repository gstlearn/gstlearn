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

#include "LinearOp/CGParam.hpp"

class ALinearOp;

class GSTLEARN_EXPORT SPDEParam {

public:
  SPDEParam(int            refineK   = 11,
            int            refineS   = 18,
            int            border    = 8,
            int            nxmax     = 300,
            double         epsNugget = EPSILON2,
            const CGParam& cgparams  = CGParam());
  SPDEParam(const SPDEParam& m);
  SPDEParam& operator=(const SPDEParam& m);
  virtual ~SPDEParam();

  int     getBorder() const { return _border; }
  CGParam getCGparams() const { return _CGparams; }
  double  getEpsNugget() const { return _epsNugget; }
  int     getRefineK() const { return _refineK; }
  int     getRefineS() const { return _refineS; }
  int     getNxMax() const { return _nxmax; }

  void setBorder(int border) { _border = border; }
  void setCGparams(const CGParam& CGparams) { _CGparams = CGparams; }
  void setEpsNugget(double epsNugget) { _epsNugget = epsNugget; }
  void setRefineK(int refineK) { _refineK = refineK; }
  void setRefineS(int refineS) { _refineS = refineS; }
  void setNxMax(int nxmax) { _nxmax = nxmax; }

private:
  int     _refineK;
  int     _refineS;
  int     _border;
  int     _nxmax;
  double  _epsNugget;
  CGParam _CGparams;
};
