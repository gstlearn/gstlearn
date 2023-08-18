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
#include "Simulation/SimuSphericalParam.hpp"

SimuSphericalParam::SimuSphericalParam(int special,
                                       int nbf,
                                       int nfmax,
                                       int degmax,
                                       int ndisc,
                                       double tol)
    : AStringable(),
      _special(special),
      _nbf(nbf),
      _nfmax(nfmax),
      _degmax(degmax),
      _ndisc(ndisc),
      _tol(tol)
{
}

SimuSphericalParam::SimuSphericalParam(const SimuSphericalParam &r)
    : AStringable(r),
      _special(r._special),
      _nbf(r._nbf),
      _nfmax(r._nfmax),
      _degmax(r._degmax),
      _ndisc(r._ndisc),
      _tol(r._tol)
{
}

SimuSphericalParam& SimuSphericalParam::operator=(const SimuSphericalParam &r)
{
  if (this != &r)
  {
    AStringable::operator =(r);
    _special = r._special;
    _nbf = r._nbf;
    _nfmax = r._nfmax;
    _degmax = r._degmax;
    _ndisc = r._ndisc;
    _tol = r._tol;
  }
  return *this;
}

SimuSphericalParam::~SimuSphericalParam()
{
}

String SimuSphericalParam::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << toTitle(1, "Option for constructing the covariance spectrum");
  if (_special == 0)
    sstr << "For all standard covariances" << std::endl;
  else if (_special == 1)
    sstr << "For Chentsov construction" << std::endl;
  else if (_special == 2)
    sstr << "For particular Exponential Model" << std::endl;

  sstr << "Number of basic functions = " << _nbf << std::endl;
  if (_nfmax > 0)
    sstr << "Maximum number of frequencies = " << _nfmax << std::endl;

  sstr << "Number of discretization  = " << _ndisc << std::endl;
  sstr << "Spectrum Tolerance        = " << _tol << std::endl;

  return sstr.str();
}
