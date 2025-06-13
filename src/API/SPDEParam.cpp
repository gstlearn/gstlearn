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

SPDEParam::SPDEParam(int refineK,
                     int refineS,
                     int border,
                     bool flag_polarized,
                     int nxmax,
                     double epsNugget,
                     bool useStencil,
                     int nMC,
                     int seedMC,
                     const CGParam& cgparams)
  : _refineK(refineK)
  , _refineS(refineS)
  , _border(border)
  , _flagPolarized(flag_polarized)
  , _nxmax(nxmax)
  , _epsNugget(epsNugget)
  , _useStencil(useStencil)
  , _nMC(nMC)
  , _seedMC(seedMC)
  , _CGparams(cgparams)
{
}

SPDEParam::SPDEParam(const SPDEParam& m)
  : _refineK(m._refineK)
  , _refineS(m._refineS)
  , _border(m._border)
  ,_flagPolarized(m._flagPolarized)
  , _nxmax(m._nxmax)
  , _epsNugget(m._epsNugget)
  , _useStencil(m._useStencil)
  , _nMC(m._nMC)
  , _seedMC(m._seedMC)
  , _CGparams(m._CGparams) {
}

SPDEParam& SPDEParam::operator=(const SPDEParam &m)
{
  if (this != &m)
  {
    _refineK       = m._refineK;
    _refineS       = m._refineS;
    _border        = m._border;
    _flagPolarized = m._flagPolarized;
    _nxmax         = m._nxmax;
    _epsNugget     = m._epsNugget;
    _useStencil    = m._useStencil;
    _nMC           = m._nMC;
    _seedMC        = m._seedMC;
    _CGparams      = m._CGparams;
  }
  return *this;
}

SPDEParam::~SPDEParam()
{
}

SPDEParam* SPDEParam::create(int refineK,
                             int refineS,
                             int border,
                             bool flag_polarized,
                             int nxmax,
                             double epsNugget,
                             bool useStencil,
                             int nMC,
                             int seedMC,
                             const CGParam& cgparams)
{
  return new SPDEParam(refineK, refineS, border, flag_polarized, nxmax,
                       epsNugget, useStencil, nMC, seedMC, cgparams);
}