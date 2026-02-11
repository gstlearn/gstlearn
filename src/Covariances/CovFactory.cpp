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
#include "Covariances/CovFactory.hpp"
#include "Basic/AException.hpp"
#include "Basic/String.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/KernelBesselJ.hpp"
#include "Covariances/KernelCauchy.hpp"
#include "Covariances/KernelCauchyGen.hpp"
#include "Covariances/KernelCosExp.hpp"
#include "Covariances/KernelCosinus.hpp"
#include "Covariances/KernelCubic.hpp"
#include "Covariances/KernelExponential.hpp"
#include "Covariances/KernelGC1.hpp"
#include "Covariances/KernelGC3.hpp"
#include "Covariances/KernelGC5.hpp"
#include "Covariances/KernelGCspline.hpp"
#include "Covariances/KernelGCspline2.hpp"
#include "Covariances/KernelGamma.hpp"
#include "Covariances/KernelGaussian.hpp"
#include "Covariances/KernelGeometric.hpp"
#include "Covariances/KernelLinear.hpp"
#include "Covariances/KernelLinearSph.hpp"
#include "Covariances/KernelMarkov.hpp"
#include "Covariances/KernelMatern.hpp"
#include "Covariances/KernelNugget.hpp"
#include "Covariances/KernelPenta.hpp"
#include "Covariances/KernelPoisson.hpp"
#include "Covariances/KernelPower.hpp"
#include "Covariances/KernelReg1D.hpp"
#include "Covariances/KernelSincard.hpp"
#include "Covariances/KernelSpherical.hpp"
#include "Covariances/KernelStable.hpp"
#include "Covariances/KernelStorkey.hpp"
#include "Covariances/KernelTriangle.hpp"
#include "Covariances/KernelWendland0.hpp"
#include "Covariances/KernelWendland1.hpp"
#include "Covariances/KernelWendland2.hpp"
#include <cctype>

namespace gstlrn
{
bool _isValid(AKernel* cova, const CovContext& ctxt)
{
  return static_cast<Id>(cova->getMaxNDim()) <= 0 ||
         static_cast<Id>(ctxt.getNDim()) <= static_cast<Id>(cova->getMaxNDim());
}

AKernel* CovFactory::createCovFunc(const ECov& type, const CovContext& ctxt)
{
  switch (type.toEnum())
  {
    case ECov::E_NUGGET: return new KernelNugget(ctxt);
    case ECov::E_EXPONENTIAL: return new KernelExponential(ctxt);
    case ECov::E_SPHERICAL: return new KernelSpherical(ctxt);
    case ECov::E_GAUSSIAN: return new KernelGaussian(ctxt);
    case ECov::E_CUBIC: return new KernelCubic(ctxt);
    case ECov::E_SINCARD: return new KernelSincard(ctxt);
    case ECov::E_BESSELJ: return new KernelBesselJ(ctxt);
    case ECov::E_MATERN: return new KernelMatern(ctxt);
    case ECov::E_GAMMA: return new KernelGamma(ctxt);
    case ECov::E_CAUCHY: return new KernelCauchy(ctxt);
    case ECov::E_CAUCHY_GEN: return new KernelCauchyGen(ctxt);
    case ECov::E_STABLE: return new KernelStable(ctxt);
    case ECov::E_LINEAR: return new KernelLinear(ctxt);
    case ECov::E_POWER: return new KernelPower(ctxt);
    case ECov::E_ORDER1_GC: return new KernelGC1(ctxt);
    case ECov::E_SPLINE_GC: return new KernelGCspline(ctxt);
    case ECov::E_SPLINE2_GC: return new KernelGCspline2(ctxt);
    case ECov::E_ORDER3_GC: return new KernelGC3(ctxt);
    case ECov::E_ORDER5_GC: return new KernelGC5(ctxt);
    case ECov::E_COSINUS: return new KernelCosinus(ctxt);
    case ECov::E_TRIANGLE: return new KernelTriangle(ctxt);
    case ECov::E_COSEXP: return new KernelCosExp(ctxt);
    case ECov::E_REG1D: return new KernelReg1D(ctxt);
    case ECov::E_PENTA: return new KernelPenta(ctxt);
    case ECov::E_STORKEY: return new KernelStorkey(ctxt);
    case ECov::E_WENDLAND0: return new KernelWendland0(ctxt);
    case ECov::E_WENDLAND1: return new KernelWendland1(ctxt);
    case ECov::E_WENDLAND2: return new KernelWendland2(ctxt);
    case ECov::E_MARKOV: return new KernelMarkov(ctxt);
    case ECov::E_GEOMETRIC: return new KernelGeometric(ctxt);
    case ECov::E_POISSON: return new KernelPoisson(ctxt);
    case ECov::E_LINEARSPH: return new KernelLinearSph(ctxt);
    default: break;
  }
  return nullptr;
}

AKernel* CovFactory::duplicateCovFunc(const AKernel& cov)
{
  switch (cov.getType().toEnum())
  {
    // Warning : if a crash with "bad cast" occurs, please check the type of your CovFunc
    case ECov::E_NUGGET: return new KernelNugget(dynamic_cast<const KernelNugget&>(cov));
    case ECov::E_EXPONENTIAL: return new KernelExponential(dynamic_cast<const KernelExponential&>(cov));
    case ECov::E_SPHERICAL: return new KernelSpherical(dynamic_cast<const KernelSpherical&>(cov));
    case ECov::E_GAUSSIAN: return new KernelGaussian(dynamic_cast<const KernelGaussian&>(cov));
    case ECov::E_CUBIC: return new KernelCubic(dynamic_cast<const KernelCubic&>(cov));
    case ECov::E_SINCARD: return new KernelSincard(dynamic_cast<const KernelSincard&>(cov));
    case ECov::E_BESSELJ: return new KernelBesselJ(dynamic_cast<const KernelBesselJ&>(cov));
    case ECov::E_MATERN: return new KernelMatern(dynamic_cast<const KernelMatern&>(cov));
    case ECov::E_GAMMA: return new KernelGamma(dynamic_cast<const KernelGamma&>(cov));
    case ECov::E_CAUCHY: return new KernelCauchy(dynamic_cast<const KernelCauchy&>(cov));
    case ECov::E_CAUCHY_GEN: return new KernelCauchyGen(dynamic_cast<const KernelCauchyGen&>(cov));
    case ECov::E_STABLE: return new KernelStable(dynamic_cast<const KernelStable&>(cov));
    case ECov::E_LINEAR: return new KernelLinear(dynamic_cast<const KernelLinear&>(cov));
    case ECov::E_POWER: return new KernelPower(dynamic_cast<const KernelPower&>(cov));
    case ECov::E_ORDER1_GC: return new KernelGC1(dynamic_cast<const KernelGC1&>(cov));
    case ECov::E_SPLINE_GC: return new KernelGCspline(dynamic_cast<const KernelGCspline&>(cov));
    case ECov::E_ORDER3_GC: return new KernelGC3(dynamic_cast<const KernelGC3&>(cov));
    case ECov::E_ORDER5_GC: return new KernelGC5(dynamic_cast<const KernelGC5&>(cov));
    case ECov::E_COSINUS: return new KernelCosinus(dynamic_cast<const KernelCosinus&>(cov));
    case ECov::E_TRIANGLE: return new KernelTriangle(dynamic_cast<const KernelTriangle&>(cov));
    case ECov::E_COSEXP: return new KernelCosExp(dynamic_cast<const KernelCosExp&>(cov));
    case ECov::E_REG1D: return new KernelReg1D(dynamic_cast<const KernelReg1D&>(cov));
    case ECov::E_PENTA: return new KernelPenta(dynamic_cast<const KernelPenta&>(cov));
    case ECov::E_SPLINE2_GC: return new KernelGCspline2(dynamic_cast<const KernelGCspline2&>(cov));
    case ECov::E_STORKEY: return new KernelStorkey(dynamic_cast<const KernelStorkey&>(cov));
    case ECov::E_WENDLAND0: return new KernelWendland0(dynamic_cast<const KernelWendland0&>(cov));
    case ECov::E_WENDLAND1: return new KernelWendland1(dynamic_cast<const KernelWendland1&>(cov));
    case ECov::E_WENDLAND2: return new KernelWendland2(dynamic_cast<const KernelWendland2&>(cov));
    case ECov::E_MARKOV: return new KernelMarkov(dynamic_cast<const KernelMarkov&>(cov));
    case ECov::E_POISSON: return new KernelPoisson(dynamic_cast<const KernelPoisson&>(cov));
    case ECov::E_GEOMETRIC: return new KernelGeometric(dynamic_cast<const KernelGeometric&>(cov));
    case ECov::E_LINEARSPH: return new KernelLinearSph(dynamic_cast<const KernelLinearSph&>(cov));
    default: break;
  }
  my_throw("Covariance function not yet implemented!");
}

/**
 * Prints the list of covariances available for a given context
 *
 * @param ctxt  Context from which we want authorized covariances
 */
void CovFactory::displayCovList(const CovContext& ctxt)
{
  message("List of authorized covariance / variogram names:\n");
  auto it = ECov::getIterator();
  while (it.hasNext())
  {
    if (*it != ECov::FUNCTION)
    {
      AKernel* cova = createCovFunc(*it, ctxt);
      if (_isValid(cova, ctxt))
        message("%2d - %s\n", it.getValue(), cova->getCovName().c_str());
      delete cova;
    }
    it.toNext();
  }
}

/**
 * Return the list of covariances names available for a given context
 *
 * @param ctxt  Context from which we want authorized covariances
 * @param order Maximum order for the IRF
 */
VectorString CovFactory::getCovList(const CovContext& ctxt, Id order)
{
  VectorString names;
  auto it = ECov::getIterator();
  while (it.hasNext())
  {
    if (*it != ECov::FUNCTION)
    {
      AKernel* cova = createCovFunc(*it, ctxt);
      if (_isValid(cova, ctxt))
      {
        if (cova->getMinOrder() <= order)
          names.push_back(cova->getCovName());
      }
      delete cova;
    }
    it.toNext();
  }
  return names;
}

/**
 * Return the ECov object from the given covariance name.
 * The name must correspond to one of the getCovName().
 * If the name doesn't exists, this method returns ECov::NUGGET
 * and display available covariances for the given context.
 *
 * @param cov_name  Name of the required covariance
 * @param ctxt      Context from which we want authorized covariances
 */
ECov CovFactory::identifyCovariance(const String& cov_name,
                                    const CovContext& ctxt)
{
  auto it = ECov::getIterator();
  while (it.hasNext())
  {
    // Test covariance name using AKernel::getCovName (not the ECov keys!)
    // (This permits to ensure RGeostats scripts retro compatibility)
    if (*it != ECov::FUNCTION)
    {
      AKernel* cova = createCovFunc(*it, ctxt);
      String cn     = toUpper(cov_name);
      String ccn    = toUpper(cova->getCovName());
      delete cova;
      if (cn == ccn)
        return *it;
    }
    it.toNext();
  }
  messerr("Unknown covariance name:%s!", cov_name.c_str());
  displayCovList(ctxt);
  return ECov::NUGGET;
}

double CovFactory::getScaleFactor(const ECov& type, double param)
{
  CovContext ctxt(1, 1);
  AKernel* cova = CovFactory::createCovFunc(type, ctxt);
  cova->setParam(param);
  double scadef = cova->getScadef();
  delete cova;
  return scadef;
}
} // namespace gstlrn