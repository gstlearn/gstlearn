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
#include "Covariances/CorFactorized.hpp"
#include "Basic/Message.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/CovFactory.hpp"
#include "Db/Db.hpp"
#include "Enum/ESimuType.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Model/ModelFitSillsVMap.hpp"
#include "Model/ModelFitSillsVario.hpp"
#include "Simulation/SpectrumOnRN.hpp"
#include "Simulation/SpectrumRN.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"

#include <cmath>
#include <memory>
#include <vector>

namespace gstlrn
{

CorFactorized::CorFactorized(const CovContext& ctxt)
  : ACov(ctxt)
  , _covs()
  , _projs()
//  , _subspaces()
{
}

CorFactorized::CorFactorized(const CorFactorized& r)
  : ACov(r)
{
  Id nfac = r.getNFac();
  for (Id ifac = 0; ifac < nfac; ifac++)
  {
    _covs.push_back(std::dynamic_pointer_cast<ACov>(r._covs[ifac]->cloneShared()));
    _projs.push_back(std::dynamic_pointer_cast<MatrixDense>(r._projs[ifac]->cloneShared()));
    //  _subspaces.addSpaceComponent(r._subspaces.getComponent(ifac));
  }
}

CorFactorized& CorFactorized::operator=(const CorFactorized& r)
{
  if (this != &r)
  {
    ACov::operator=(r);
    Id nfac = r.getNFac();
    for (Id ifac = 0; ifac < nfac; ifac++)
    {
      _covs.push_back(std::dynamic_pointer_cast<ACov>(r._covs[ifac]->cloneShared()));
      _projs.push_back(std::dynamic_pointer_cast<MatrixDense>(r._projs[ifac]->cloneShared()));
      //  _subspaces.addSpaceComponent(r._subspaces.getComponent(ifac));
    }
  }
  return *this;
}

CorFactorized::~CorFactorized()
{
  deleteAllFactors();
}

CorFactorized* CorFactorized::create(const CovContext& ctxt)
{
  auto* cov = new CorFactorized(ctxt);
  return cov;
}

double CorFactorized::eval0(Id ivar, Id jvar, const CovCalcMode* mode) const
{
  double cov = 1.0;
  Id nfac    = getNFac();
  for (Id j = 0; j < nfac; j++)
  {
    cov *= _covs[j]->eval0(ivar, jvar, mode);
  }
  return cov;
}

String CorFactorized::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  if (getNFac() <= 0) return sstr.str();

  for (Id ifac = 0, nfac = getNFac(); ifac < nfac; ifac++)
  {
    sstr << _covs[ifac]->toString();
  }
  sstr << std::endl;

  return sstr.str();
}

Id CorFactorized::getNFac() const
{
  Id nfac = static_cast<Id>(_covs.size());
  return nfac;
}

void CorFactorized::addFactor(const ACov& cov, const MatrixDense& proj)
{
  Id nvar     = getNVar();
  Id ndim     = getNDim();
  Id ndim_sub = proj.getNCols();
  if ((cov.getNVar() != 1) && (cov.getNVar() != nvar))
  {
    messerr("Error: the number of variables is inconsistent (nv = %d != %d)", cov.getNVar(), nvar);
    messerr("Operation is cancelled");
  }
  if (cov.getNDim() > ndim)
  {
    messerr("Error: dimension of covariance inconsistent (d = %d > %d)", cov.getNDim(), ndim);
    messerr("Operation is cancelled");
  }
  if (proj.getNRows() != ndim)
  {
    messerr("Error: subspace dimension inconsistent (d = %d != %d)", proj.getNRows(), ndim);
    messerr("Operation is cancelled");
  }
  if (ndim_sub > ndim)
  {
    messerr("Error: subspace dimension inconsistent (d = %d > %d)", proj.getNCols(), ndim);
    messerr("Operation is cancelled");
  }
  _covs.push_back(std::dynamic_pointer_cast<ACov>(cov.cloneShared()));
  _projs.push_back(std::dynamic_pointer_cast<MatrixDense>(proj.cloneShared()));
  //_subspaces.addSpaceComponent(cov.getSpace());
}

void CorFactorized::deleteAllFactors()
{
  _covs.clear();
  _projs.clear();
}

const ACov* CorFactorized::getCovariance(Id ifac) const
{
  if (!_isFactorIndexValid(ifac)) return nullptr;
  return _covs[ifac].get();
}

MatrixDense CorFactorized::getProjection(Id ifac) const
{
  if (!_isFactorIndexValid(ifac)) return MatrixDense();
  return MatrixDense(*_projs[ifac]);
}

bool CorFactorized::isValidForSimulation(const ESimuType& simuType) const
{
  /* Loop on the factors */
  for (int ifac = 0; ifac < getNFac(); ifac++)
  {
    if (!_covs[ifac]->isValidForSimulation(simuType))
    {
      messerr("The structure %d# is not valid for %s Simulation on RN", ifac, simuType.getKey());
      return false;
    }
  }
  return true;
}

SpectrumRN CorFactorized::simulateSpectrum(Id ns) const
{
  // TODO: SpectrumRN should be replaced by SpectrumOnRN
  DECLARE_UNUSED(ns)
  messerr("deprecated methods: simulateOnRN should be used (nfac = %d)", getNFac());
  return SpectrumRN();
}

SpectrumOnRN* CorFactorized::simulateOnRN(Id ns) const
{
  if (!isValidForSimulation(ESimuType::SPECTRAL))
  {
    messerr("Covariance not valid for spectral simulation on RN");
    return nullptr;
  }
  Id ndim = getNDim();
  Id nvar = getNVar();
  Id nfac = getNFac();
  MatrixDense gamma(ns, nvar);
  for (Id ivar = 0; ivar < nvar; ivar++)
    gamma.setColumnToConstant(ivar, 1.0);

  auto* res = new SpectrumOnRNFactorized(nvar, ndim, ns);
  for (Id ifac = 0; ifac < nfac; ifac++)
  {
    const ACov* cova        = getCovariance(ifac);
    MatrixDense proj = getProjection(ifac);
    if (cova->getNFac() > 1)
    {
      messerr("Non elementary covariance for factor #%d", ifac);
      continue;
    }
    if ((cova->getNVar() != 1) && (cova->getNVar() != getNVar()))
    {
      messerr("Inconsistent number of variables for factor #%d (nvar = %d != %d)", ifac, cova->getNVar(), getNVar());
      continue;
    }
    SpectrumOnRN* spf = cova->simulateOnRN(ns);
    res->addFactor(spf->getOmega(0), spf->getPhi(0), proj);
    const MatrixDense& gf = spf->getGamma();
    if (gf.getNCols() == 1)
    {
      gamma.multiplyRow(gf.getColumn(0));
    }
    else
    {
      // TODO: can we multiply term by term?
      for (Id i = 0; i < ns; i++)
      {
        for (Id j = 0; j < nvar; j++)
        {
          gamma.setValue(i, j, gamma.getValue(i, j) * gf.getValue(i, j));
        }
      }
    }
  }
  double coeff = (nfac - 1.) / 2.0;
  gamma.prodCst(pow(ns, coeff));
  res->setGamma(gamma);
  return res;
}

bool CorFactorized::_isFactorIndexValid(Id ifac) const
{
  return checkArg("Factor Index", ifac, getNFac());
}

double CorFactorized::_eval(const SpacePoint& p1,
                            const SpacePoint& p2,
                            Id ivar,
                            Id jvar,
                            const CovCalcMode* mode) const
{
  double cov = 1.0;
  Id nd      = getNDim();
  VectorDouble x1(p1.getCoords());
  VectorDouble x2(p2.getCoords());
  VectorDouble x1_proj(nd);
  VectorDouble x2_proj(nd);

  for (Id ifac = 0; ifac < getNFac(); ifac++)
  {
    SpacePoint p1_proj(_covs[ifac]->getSpace());
    SpacePoint p2_proj(_covs[ifac]->getSpace());
    for (Id j = 0; j < _covs[ifac]->getNDim(); j++)
    {
      x1_proj[j] = 0.0;
      x2_proj[j] = 0.0;
      for (Id k = 0; k < nd; k++)
      {
        x1_proj[j] += (x1[k] * _projs[ifac]->getValue(k, j));
        x2_proj[j] += (x2[k] * _projs[ifac]->getValue(k, j));
      }
      p1_proj.setCoord(j, x1_proj[j]);
      p2_proj.setCoord(j, x2_proj[j]);
    }
    cov *= _covs[ifac]->evalCov(p1_proj, p2_proj, ivar, jvar, mode);
  }
  return cov;
}

} // namespace gstlrn
