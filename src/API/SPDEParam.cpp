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

namespace gstlrn
{
SPDEParam::SPDEParam(Id refineK,
                     Id refineS,
                     Id border,
                     bool flag_polarized,
                     Id nxmax,
                     double epsNugget,
                     bool useStencil,
                     Id nMC,
                     Id seedMC,
                     const CGParam& cgparams)
  : AStringable()
  , _refineK(refineK)
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
  : AStringable(m)
  , _refineK(m._refineK)
  , _refineS(m._refineS)
  , _border(m._border)
  , _flagPolarized(m._flagPolarized)
  , _nxmax(m._nxmax)
  , _epsNugget(m._epsNugget)
  , _useStencil(m._useStencil)
  , _nMC(m._nMC)
  , _seedMC(m._seedMC)
  , _CGparams(m._CGparams)
{
}

SPDEParam& SPDEParam::operator=(const SPDEParam& m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
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

SPDEParam* SPDEParam::create(Id refineK,
                             Id refineS,
                             Id border,
                             bool flag_polarized,
                             Id nxmax,
                             double epsNugget,
                             bool useStencil,
                             Id nMC,
                             Id seedMC,
                             const CGParam& cgparams)
{
  return new SPDEParam(refineK, refineS, border, flag_polarized, nxmax,
                       epsNugget, useStencil, nMC, seedMC, cgparams);
}

String SPDEParam::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);
  std::stringstream sstr;

  sstr << "Discretization factor for Kriging = " << _refineK << std::endl;
  sstr << "Discretization factor for Simulation = " << _refineS << std::endl;
  sstr << "Border Size  = " << _border << std::endl;
  sstr << "Nugget effect = " << _epsNugget << std::endl;
  sstr << "Default option for no Cholesky (only used in stationary case for Turbo Meshing) = " << _useStencil << std::endl;
  sstr << "Number of Monte-Carlo simulations (used for Variance and logdet) = " << _nMC << std::endl;
  sstr << "Seed for the random number generator (used for Variance and logdet) = " << _seedMC << std::endl;

  sstr << _CGparams.toString() << std::endl;
  return sstr.str();
}
}