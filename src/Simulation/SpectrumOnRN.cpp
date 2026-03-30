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
#include "Simulation/SpectrumOnRN.hpp"
#include "Basic/Message.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorT.hpp"
#include "Db/Db.hpp"
#include "Matrix/MatrixDense.hpp"
#include "geoslib_define.h"

#include <cmath>

namespace gstlrn
{
///////////////////////////////////////////////////////////////////////////////////////////////////////
// Abstract class for the spectral simulation: SpectrumOnRN
///////////////////////////////////////////////////////////////////////////////////////////////////////

SpectrumOnRN::SpectrumOnRN(Id nv, Id nd, Id ns)
  : _nv(nv)
  , _nd(nd)
  , _ns(ns)
  , _gamma(ns, nv)
  , _u(ns)
  , _v(ns)
  , _SQRT2(sqrt(2.0))
{
}

SpectrumOnRN::SpectrumOnRN(const SpectrumOnRN& r)
  : _nv(r._nv)
  , _nd(r._nd)
  , _ns(r._ns)
  , _gamma(r._gamma)
  , _u(r._ns)
  , _v(r._ns)
  , _SQRT2(r._SQRT2)
{
}

SpectrumOnRN& SpectrumOnRN::operator=(const SpectrumOnRN& r)
{
  if (this != &r)
  {
    _nv    = r._nv;
    _nd    = r._nd;
    _ns    = r._ns;
    _gamma = r._gamma;
    _u     = r._u;
    _v     = r._v;
    _SQRT2 = r._SQRT2;
  }
  return *this;
}

SpectrumOnRN::~SpectrumOnRN()
{
}

bool SpectrumOnRN::setGamma(const MatrixDense& gamma)
{
  if ((_ns != gamma.getNRows()))
  {
    return false;
  }
  _nv = gamma.getNCols();
  _gamma.reset(_ns, _nv);
  for (Id iv = 0; iv < _nv; iv++)
  {
    VectorDouble tab = gamma.getColumn(iv);
    _gamma.setColumn(iv, tab);
  }
  return true;
}

void SpectrumOnRN::compute(Db* dbout, const VectorBool& activeArray, VectorVectorDouble& tab)
{
  // Loop on the samples
  auto nech = dbout->getNSample();
  VectorDouble x(_nd);
  VectorDouble values(_nv);
  _initWorkingVariables();
  for (Id iech = 0; iech < nech; iech++)
  {
    if (!activeArray[iech]) continue;
    dbout->getCoordinatesInPlace(x, iech);
    _computeValues(x, values);
    for (Id iv = 0; iv < _nv; iv++)
    {
      tab[iv][iech] = values[iv];
    }
  }
}
MatrixDense
SpectrumOnRN::computeToMatrix(Db* dbout)
{
  auto nvar = getNVar();
  auto nech = dbout->getNSample();
  VectorVectorDouble tab(nvar);
  for (Id ivar = 0; ivar < nvar; ivar++) tab[ivar].resize(nech);
  VectorBool activeArray = dbout->getActiveArray();
  compute(dbout, activeArray, tab);
  MatrixDense res(nech, nvar);
  for (Id ivar = 0; ivar < nvar; ivar++)
  {
    res.setColumn(ivar, tab[ivar]);
  }
  return res;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// concrete class defining a simple spectrum on RN
///////////////////////////////////////////////////////////////////////////////////////////////////////
SpectrumOnRNSimple::SpectrumOnRNSimple(Id nv, Id nd, Id ns)
  : SpectrumOnRN(nv, nd, ns)
  , _omega(ns, nd)
  , _phi(ns)
{
}

SpectrumOnRNSimple::SpectrumOnRNSimple(const SpectrumOnRNSimple& r)
  : SpectrumOnRN(r)
  , _omega(r._ns, _nd)
  , _phi(r._ns)
{
}

SpectrumOnRNSimple& SpectrumOnRNSimple::operator=(const SpectrumOnRNSimple& r)
{
  if (this != &r)
  {
    _nv    = r._nv;
    _nd    = r._nd;
    _ns    = r._ns;
    _gamma = r._gamma;
    _omega = r._omega;
    _phi   = r._phi;
  }
  return *this;
}

SpectrumOnRNSimple::~SpectrumOnRNSimple()
{
}

// virtual functions
bool SpectrumOnRNSimple::addFactor(
  const MatrixDense& omega,
  const VectorDouble& phi,
  const MatrixDense& proj,
  const MatrixDense& omega0,
  const VectorDouble& xi0)
{
  DECLARE_UNUSED(proj)
  if ((omega.getNRows() != _ns) || (omega.getNCols() != _nd))
  {
    messerr("Inconsistent dimensions of the frequency matrix: nrow = %d, ncol = %d (ns = %d, nv =  %d)",
            omega.getNRows(), omega.getNCols(), _ns, _nd);
    return false;
  }
  if ((phi.length() != _ns))
  {
    messerr("Inconsistent length of the phase vector: length = %d (ns = %d)",
            phi.length(), _ns);
    return false;
  }
  if ((omega0.getNRows() != 0) && ((omega0.getNRows() != _ns) || (omega0.getNCols() != _nd)))
  {
    messerr("Inconsistent dimensions of the raw frequency matrix: nrow = %d, ncol = %d (ns = %d, nv =  %d)",
            omega0.getNRows(), omega0.getNCols(), _ns, _nd);
    return false;
  }
  if ((xi0.length() != 0) && (xi0.length() != _ns))
  {
    messerr("Inconsistent length of the mixture coefficient vector: length = %d (ns = %d)",
            xi0.length(), _ns);
    return false;
  }
  _omega  = omega;
  _phi    = phi;
  _omega0 = omega0;
  _xi0    = xi0;
  return true;
}

MatrixDense SpectrumOnRNSimple::getOmega(Id ifac, Id is) const
{
  DECLARE_UNUSED(ifac)
  DECLARE_UNUSED(is)
  return _omega;
}

VectorDouble SpectrumOnRNSimple::getPhi(Id ifac, Id is) const
{
  DECLARE_UNUSED(ifac)
  DECLARE_UNUSED(is)
  return _phi;
}

MatrixDense SpectrumOnRNSimple::getProjection(Id ifac, Id is) const
{
  DECLARE_UNUSED(ifac)
  DECLARE_UNUSED(is)
  MatrixDense id(_nd, _nd);
  id.setDiagonalToConstant();
  return id;
}

MatrixDense SpectrumOnRNSimple::getOmega0(Id ifac, Id is) const
{
  DECLARE_UNUSED(ifac)
  DECLARE_UNUSED(is)
  return _omega0;
}

VectorDouble SpectrumOnRNSimple::getXi0(Id ifac, Id is) const
{
  DECLARE_UNUSED(ifac)
  DECLARE_UNUSED(is)
  return _xi0;
}

void SpectrumOnRNSimple::_initWorkingVariables()
{
  // Do nothing
}

void SpectrumOnRNSimple::_computeValues(const VectorDouble& coor, VectorDouble& values)
{
  MatrixDense::productInPlace(_u, _omega, coor, false);
  for (Id is = 0; is < _ns; is++)
  {
    _u[is] = (_SQRT2 * cos(_u[is] + _phi[is]));
  }
  MatrixDense::productInPlace(values, _gamma, _u, true);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// concrete class defining a factorized spectrum on RN for the simulation of CorFactorized and CorGneiting
///////////////////////////////////////////////////////////////////////////////////////////////////////
SpectrumOnRNFactorized::SpectrumOnRNFactorized(Id nv, Id nd, Id ns)
  : SpectrumOnRN(nv, nd, ns)
  , _nf(0)
  , _omegas()
  , _phis()
  , _projections()
  , _xproj()
{
}

SpectrumOnRNFactorized::SpectrumOnRNFactorized(const SpectrumOnRNFactorized& r)
  : SpectrumOnRN(r)
  , _nf(0)
  , _omegas()
  , _phis()
  , _projections()
  , _xproj()
{
}

SpectrumOnRNFactorized& SpectrumOnRNFactorized::operator=(const SpectrumOnRNFactorized& r)
{
  if (this != &r)
  {
    _nv          = r._nv;
    _nd          = r._nd;
    _ns          = r._ns;
    _gamma       = r._gamma;
    _nf          = r._nf;
    _omegas      = r._omegas;
    _phis        = r._phis;
    _projections = r._projections;
    _xproj       = r._xproj;
  }
  return *this;
}

SpectrumOnRNFactorized::~SpectrumOnRNFactorized()
{
}

bool SpectrumOnRNFactorized::addFactor(
  const MatrixDense& omega,
  const VectorDouble& phi,
  const MatrixDense& proj,
  const MatrixDense& omega0,
  const VectorDouble& xi0)
{
  DECLARE_UNUSED(omega0)
  DECLARE_UNUSED(xi0)
  if ((omega.getNRows() != _ns) || (omega.getNCols() > _nd) || (omega.getNCols() != proj.getNCols()))
  {
    messerr("Inconsistent dimensions of the frequency matrix: nrow = %d, ncol = %d (ns = %d, nv =  %d, nproj = %d)",
            omega.getNRows(), omega.getNCols(), _ns, _nv, proj.getNCols());
    return false;
  }
  if ((phi.length() != _ns))
  {
    messerr("Inconsistent length of the phase vector: length = %d (ns = %d)",
            phi.length(), _ns);
    return false;
  }
  if ((proj.getNRows() != _nd) || (proj.getNCols() > _nd))
  {
    messerr("Inconsistent dimensions of the projection matrix: nrow = %d, ncol = %d (nd = %d)",
            proj.getNRows(), proj.getNCols(), _nd);
    return false;
  }
  _omegas.push_back(omega);
  _phis.push_back(phi);
  _projections.push_back(proj);
  _nf++;
  return true;
}

Id SpectrumOnRNFactorized::getNFac() const
{
  return _nf;
}

MatrixDense SpectrumOnRNFactorized::getOmega(Id ifac, Id is) const
{
  DECLARE_UNUSED(is)
  if (!_isValidNf(ifac))
  {
    messerr("Invalid factor index: ifac = %d (nf = %d)", ifac, _nf);
    return MatrixDense();
  }
  return _omegas[ifac];
}

VectorDouble SpectrumOnRNFactorized::getPhi(Id ifac, Id is) const
{
  DECLARE_UNUSED(is)
  if (!_isValidNf(ifac))
  {
    messerr("Invalid factor index: ifac = %d (nf = %d)", ifac, _nf);
    return VectorDouble();
  }
  return _phis[ifac];
}

MatrixDense SpectrumOnRNFactorized::getOmega0(Id ifac, Id is) const
{
  DECLARE_UNUSED(is)
  DECLARE_UNUSED(ifac)
  messerr("Raw frequency matrix not defined for factorized spectrum");
  return MatrixDense();
}

VectorDouble SpectrumOnRNFactorized::getXi0(Id ifac, Id is) const
{
  DECLARE_UNUSED(is)
  DECLARE_UNUSED(ifac)
  messerr("Mixture coefficient vector not defined for factorized spectrum");
  return VectorDouble();
}

MatrixDense SpectrumOnRNFactorized::getProjection(Id ifac, Id is) const
{
  DECLARE_UNUSED(is)
  if (!_isValidNf(ifac))
  {
    messerr("Invalid factor index: ifac = %d (nf = %d)", ifac, _nf);
    return MatrixDense();
  }
  return _projections[ifac];
}

void SpectrumOnRNFactorized::_initWorkingVariables()
{
  _xproj.clear();
  for (Id i = 0; i < _nf; i++)
  {
    Id ndim_sub = _projections[i].getNCols();
    _xproj.push_back(VectorDouble(ndim_sub));
  }
}

void SpectrumOnRNFactorized::_computeValues(const VectorDouble& coor, VectorDouble& values)
{
  for (Id is = 0; is < _ns; is++)
  {
    _u[is] = 1.0;
  }
  // computing the monochromatic waves for each factor
  for (Id ifac = 0; ifac < _nf; ifac++)
  {
    MatrixDense::productInPlace(_xproj[ifac],_projections[ifac], coor , true);
    MatrixDense::productInPlace(_v, _omegas[ifac], _xproj[ifac], false);
    for (Id ib = 0; ib < _ns; ib++)
    {
      _u[ib] *= (_SQRT2 * cos(_v[ib] + _phis[ifac][ib]));
    }
  }
  // linear combination of the monochromatic waves
  MatrixDense::productInPlace(values, _gamma, _u, true);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// concrete class defining a list of spectrums on RN for the simulation of CovList
///////////////////////////////////////////////////////////////////////////////////////////////////////
SpectrumOnRNList::SpectrumOnRNList(Id nv, Id nd, Id ns)
  : SpectrumOnRN(nv, nd, 0)
  , _sps() {
      DECLARE_UNUSED(ns)}

  SpectrumOnRNList::SpectrumOnRNList(const SpectrumOnRNList& r)
  : SpectrumOnRN(r)
  , _sps()
{
  for (const auto& ptr: r._sps)
  {
    _sps.push_back(std::unique_ptr<SpectrumOnRN>(ptr->clone()));
  }
}

SpectrumOnRNList& SpectrumOnRNList::operator=(const SpectrumOnRNList& r)
{
  if (this != &r)
  {
    _nv    = r._nv;
    _nd    = r._nd;
    _ns    = r._ns;
    _gamma = r._gamma;
    for (const auto& ptr: r._sps)
    {
      _sps.push_back(std::unique_ptr<SpectrumOnRN>(ptr->clone()));
    }
  }
  return *this;
}

SpectrumOnRNList::~SpectrumOnRNList()
{
}

bool SpectrumOnRNList::addFactor(
  const MatrixDense& omega, 
  const VectorDouble& phi, 
  const MatrixDense& proj,
  const MatrixDense& omega0,
  const VectorDouble& xi0)
{
  DECLARE_UNUSED(omega)
  DECLARE_UNUSED(phi)
  DECLARE_UNUSED(proj)
  DECLARE_UNUSED(omega0)
  DECLARE_UNUSED(xi0)
  return false;
}

MatrixDense SpectrumOnRNList::getOmega(Id ifac, Id is) const
{
  if (_isValidNSpectrum(is) && (_sps[is]->_isValidNf(ifac)))
  {
    return _sps[is]->getOmega(ifac);
  }
  messerr("Omega not defined for ifac = %d and is = %d", ifac, is);
  return MatrixDense();
}

VectorDouble SpectrumOnRNList::getPhi(Id ifac, Id is) const
{
  if (_isValidNSpectrum(is) && (_sps[is]->_isValidNf(ifac)))
  {
    return _sps[is]->getPhi(ifac);
  }
  messerr("Phi not defined for ifac = %d and is = %d", ifac, is);
  return VectorDouble();
}

MatrixDense SpectrumOnRNList::getOmega0(Id ifac, Id is) const
{
  if (_isValidNSpectrum(is) && (_sps[is]->_isValidNf(ifac)))
  {
    return _sps[is]->getOmega0(ifac);
  }
  messerr("Omega0 not defined for ifac = %d and is = %d", ifac, is);
  return MatrixDense();
}

VectorDouble SpectrumOnRNList::getXi0(Id ifac, Id is) const
{
  if (_isValidNSpectrum(is) && (_sps[is]->_isValidNf(ifac)))
  {
    return _sps[is]->getXi0(ifac);
  }
  messerr("Xi0 not defined for ifac = %d and is = %d", ifac, is);
  return VectorDouble();
}

MatrixDense SpectrumOnRNList::getProjection(Id ifac, Id is) const
{
  if (_isValidNSpectrum(is) && (_sps[is]->_isValidNf(ifac)))
  {
    return _sps[is]->getProjection(ifac);
  }
  messerr("Projection not defined for ifac = %d and is = %d", ifac, is);
  return MatrixDense();
}

bool SpectrumOnRNList::setGamma(const MatrixDense& gamma) {
  DECLARE_UNUSED(gamma)
  messerr("Setting the Gamma matrix not possible for a SpectrumOnRNList...");
  return false;
}

MatrixDense SpectrumOnRNList::getGamma(Id is) const
{
  if (_isValidNSpectrum(is))
  {
    return _sps[is]->getGamma();
  }
  messerr("Gamma not defined for is = %d", is);
  return MatrixDense();
}

void SpectrumOnRNList::_initWorkingVariables()
{
  for (Id is = 0; is < getNSpectrum(); is++)
  {
    _sps[is]->_initWorkingVariables();
  }
}

bool SpectrumOnRNList::addSpectrum(std::unique_ptr<SpectrumOnRN> sp)
{
  if (sp->getNVar() != getNVar())
  {
    messerr("Inconsistent number of variables = %d (vs. nv = %d)", sp->getNVar(), getNVar());
    return false;
  }
  if (sp->getNDim() != getNDim())
  {
    messerr("Inconsistent space dimension = %d (vs. nd = %d)", sp->getNDim(), getNDim());
    return false;
  }
  _ns += sp->getNs();
  _sps.push_back(std::move(sp));
  return true;
};

void SpectrumOnRNList::_computeValues(const VectorDouble& coor, VectorDouble& values)
{
  VectorDouble val(_nv);
  for (Id iv = 0; iv < _nv; iv++)
  {
    values[iv] = 0.0;
  }
  for (Id is = 0; is < getNSpectrum(); is++)
  {
    _sps[is]->_computeValues(coor, val);
    for (Id iv = 0; iv < _nv; iv++)
    {
      values[iv] += val[iv];
    }
  }
}

} // namespace gstlrn