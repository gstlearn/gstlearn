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
#include "API/SPDEParam.hpp"

SPDEParam::SPDEParam(int            refineK,
                     int            refineS,
                     int            border,
                     bool flag_polarized,
                     int            nxmax,
                     double         epsNugget,
                     const CGParam& cgparams)
  : _refineK(refineK)
  , _refineS(refineS)
  , _border(border)
  ,_flagPolarized(flag_polarized)
  , _nxmax(nxmax)
  , _epsNugget(epsNugget)
  , _CGparams(cgparams) {
}

SPDEParam::SPDEParam(const SPDEParam& m)
  : _refineK(m._refineK)
  , _refineS(m._refineS)
  , _border(m._border)
  ,_flagPolarized(m._flagPolarized)
  , _nxmax(m._nxmax)
  , _epsNugget(m._epsNugget)
  , _CGparams(m._CGparams) {
}

SPDEParam& SPDEParam::operator=(const SPDEParam &m)
{
  if (this != &m) {
    _refineK   = m._refineK;
    _refineS   = m._refineS;
    _border    = m._border;
    _flagPolarized = m._flagPolarized;
    _nxmax     = m._nxmax;
    _epsNugget = m._epsNugget;
    _CGparams  = m._CGparams;
  }
  return *this;
}

SPDEParam::~SPDEParam()
{
}
