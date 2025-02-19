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

/**
 * @brief Definition of the parameters used within SPDE
 *
 * refineK Discretization factor used for Kriging
 * refineS Discretization factor used for Simulation
 * border  Border size
 * nxmax   Maximum number of vertices in the internal mesh (0 : no limit)
 * epsNugget Nugget effect
 * cgparams Parameters for the Conjugate Gradient method
 */
class GSTLEARN_EXPORT SPDEParam {

public:
  SPDEParam(int            refineK   = 11,
            int            refineS   = 18,
            int            border    = 8,
            bool flag_polarized = true,
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
  bool isPolarized() const { return _flagPolarized; }
  void setPolarized(bool flagPolarized) { _flagPolarized = flagPolarized; }  int     getNxMax() const { return _nxmax; }

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
  bool _flagPolarized;
  int     _nxmax;
  double  _epsNugget;
  CGParam _CGparams;
};
