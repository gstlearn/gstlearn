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
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/AException.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovBesselJ.hpp"
#include "Covariances/CovMatern.hpp"
#include "Covariances/CovCauchy.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/CovCosExp.hpp"
#include "Covariances/CovCosinus.hpp"
#include "Covariances/CovCubic.hpp"
#include "Covariances/CovExponential.hpp"
#include "Covariances/CovGamma.hpp"
#include "Covariances/CovGaussian.hpp"
#include "Covariances/CovGC1.hpp"
#include "Covariances/CovGC3.hpp"
#include "Covariances/CovGC5.hpp"
#include "Covariances/CovGCspline.hpp"
#include "Covariances/CovGCspline2.hpp"
#include "Covariances/CovLinear.hpp"
#include "Covariances/CovNugget.hpp"
#include "Covariances/CovPenta.hpp"
#include "Covariances/CovPower.hpp"
#include "Covariances/CovReg1D.hpp"
#include "Covariances/CovSincard.hpp"
#include "Covariances/CovSpherical.hpp"
#include "Covariances/CovStable.hpp"
#include "Covariances/CovStorkey.hpp"
#include "Covariances/CovTriangle.hpp"
#include "Covariances/CovWendland0.hpp"
#include "Covariances/CovWendland1.hpp"
#include "Covariances/CovWendland2.hpp"
#include "Covariances/CovMarkov.hpp"
#include "Covariances/CovGeometric.hpp"
#include "Covariances/CovPoisson.hpp"
#include "Covariances/CovLinearSph.hpp"

#include <iostream>
#include <cctype>

bool _isValid(ACovFunc* cova, const CovContext& ctxt)
{
  if ((int) cova->getMaxNDim() > 0   &&
      (int) ctxt.getNDim() > (int) cova->getMaxNDim())  return false;
  return true;
}

ACovFunc* CovFactory::createCovFunc(const ECov& type, const CovContext& ctxt)
{
  switch(type.toEnum())
  {
    case ECov::E_NUGGET:      return new CovNugget(ctxt);
    case ECov::E_EXPONENTIAL: return new CovExponential(ctxt);
    case ECov::E_SPHERICAL:   return new CovSpherical(ctxt);
    case ECov::E_GAUSSIAN:    return new CovGaussian(ctxt);
    case ECov::E_CUBIC:       return new CovCubic(ctxt);
    case ECov::E_SINCARD:     return new CovSincard(ctxt);
    case ECov::E_BESSELJ:     return new CovBesselJ(ctxt);
    case ECov::E_MATERN:      return new CovMatern(ctxt);
    case ECov::E_GAMMA:       return new CovGamma(ctxt);
    case ECov::E_CAUCHY:      return new CovCauchy(ctxt);
    case ECov::E_STABLE:      return new CovStable(ctxt);
    case ECov::E_LINEAR:      return new CovLinear(ctxt);
    case ECov::E_POWER:       return new CovPower(ctxt);
    case ECov::E_ORDER1_GC:   return new CovGC1(ctxt);
    case ECov::E_SPLINE_GC:   return new CovGCspline(ctxt);
    case ECov::E_SPLINE2_GC:  return new CovGCspline2(ctxt);
    case ECov::E_ORDER3_GC:   return new CovGC3(ctxt);
    case ECov::E_ORDER5_GC:   return new CovGC5(ctxt);
    case ECov::E_COSINUS:     return new CovCosinus(ctxt);
    case ECov::E_TRIANGLE:    return new CovTriangle(ctxt);
    case ECov::E_COSEXP:      return new CovCosExp(ctxt);
    case ECov::E_REG1D:       return new CovReg1D(ctxt);
    case ECov::E_PENTA:       return new CovPenta(ctxt);
    case ECov::E_STORKEY:     return new CovStorkey(ctxt);
    case ECov::E_WENDLAND0:   return new CovWendland0(ctxt);
    case ECov::E_WENDLAND1:   return new CovWendland1(ctxt);
    case ECov::E_WENDLAND2:   return new CovWendland2(ctxt);
    case ECov::E_MARKOV:      return new CovMarkov(ctxt);
    case ECov::E_GEOMETRIC:   return new CovGeometric(ctxt);
    case ECov::E_POISSON:     return new CovPoisson(ctxt);
    case ECov::E_LINEARSPH:   return new CovLinearSph(ctxt);
    default: break;
  }
  return nullptr;
}

ACovFunc* CovFactory::duplicateCovFunc(const ACovFunc& cov)
{
  switch(cov.getType().toEnum())
  {
    // Warning : if a crash with "bad cast" occurs, please check the type of your CovFunc
    case ECov::E_NUGGET:      return new CovNugget(     dynamic_cast<const CovNugget&>     (cov));
    case ECov::E_EXPONENTIAL: return new CovExponential(dynamic_cast<const CovExponential&>(cov));
    case ECov::E_SPHERICAL:   return new CovSpherical(  dynamic_cast<const CovSpherical&>  (cov));
    case ECov::E_GAUSSIAN:    return new CovGaussian(   dynamic_cast<const CovGaussian&>   (cov));
    case ECov::E_CUBIC:       return new CovCubic(      dynamic_cast<const CovCubic&>      (cov));
    case ECov::E_SINCARD:     return new CovSincard(    dynamic_cast<const CovSincard&>    (cov));
    case ECov::E_BESSELJ:     return new CovBesselJ(    dynamic_cast<const CovBesselJ&>    (cov));
    case ECov::E_MATERN:      return new CovMatern(     dynamic_cast<const CovMatern&>     (cov));
    case ECov::E_GAMMA:       return new CovGamma(      dynamic_cast<const CovGamma&>      (cov));
    case ECov::E_CAUCHY:      return new CovCauchy(     dynamic_cast<const CovCauchy&>     (cov));
    case ECov::E_STABLE:      return new CovStable(     dynamic_cast<const CovStable&>     (cov));
    case ECov::E_LINEAR:      return new CovLinear(     dynamic_cast<const CovLinear&>     (cov));
    case ECov::E_POWER:       return new CovPower(      dynamic_cast<const CovPower&>      (cov));
    case ECov::E_ORDER1_GC:   return new CovGC1(        dynamic_cast<const CovGC1&>        (cov));
    case ECov::E_SPLINE_GC:   return new CovGCspline(   dynamic_cast<const CovGCspline&>   (cov));
    case ECov::E_ORDER3_GC:   return new CovGC3(        dynamic_cast<const CovGC3&>        (cov));
    case ECov::E_ORDER5_GC:   return new CovGC5(        dynamic_cast<const CovGC5&>        (cov));
    case ECov::E_COSINUS:     return new CovCosinus(    dynamic_cast<const CovCosinus&>    (cov));
    case ECov::E_TRIANGLE:    return new CovTriangle(   dynamic_cast<const CovTriangle&>   (cov));
    case ECov::E_COSEXP:      return new CovCosExp(     dynamic_cast<const CovCosExp&>     (cov));
    case ECov::E_REG1D:       return new CovReg1D(      dynamic_cast<const CovReg1D&>      (cov));
    case ECov::E_PENTA:       return new CovPenta(      dynamic_cast<const CovPenta&>      (cov));
    case ECov::E_SPLINE2_GC:  return new CovGCspline2(  dynamic_cast<const CovGCspline2&>  (cov));
    case ECov::E_STORKEY:     return new CovStorkey(    dynamic_cast<const CovStorkey&>    (cov));
    case ECov::E_WENDLAND0:   return new CovWendland0(  dynamic_cast<const CovWendland0&>  (cov));
    case ECov::E_WENDLAND1:   return new CovWendland1(  dynamic_cast<const CovWendland1&>  (cov));
    case ECov::E_WENDLAND2:   return new CovWendland2(  dynamic_cast<const CovWendland2&>  (cov));
    case ECov::E_MARKOV:      return new CovMarkov(     dynamic_cast<const CovMarkov&>     (cov));
    case ECov::E_POISSON:     return new CovPoisson(    dynamic_cast<const CovPoisson&>    (cov));
    case ECov::E_GEOMETRIC:   return new CovGeometric(  dynamic_cast<const CovGeometric&>  (cov));
    case ECov::E_LINEARSPH:   return new CovLinearSph(  dynamic_cast<const CovLinearSph&>  (cov));
    default: break;
  }
  my_throw ("Covariance function not yet implemented!");
  return nullptr;
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
    if (*it != ECov::UNKNOWN && *it != ECov::FUNCTION)
    {
      ACovFunc* cova = createCovFunc(*it, ctxt);
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
VectorString CovFactory::getCovList(const CovContext& ctxt, int order)
{
  VectorString names;
  auto it = ECov::getIterator();
  while (it.hasNext())
  {
    if (*it != ECov::UNKNOWN && *it != ECov::FUNCTION)
    {
      ACovFunc* cova = createCovFunc(*it, ctxt);
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
 * If the name doesn't exists, this method returns ECov::UNKNOWN
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
    // Test covariance name using ACovFunc::getCovName (not the ECov keys!)
    // (This permits to ensure RGeostats scripts retro compatibility)
    if (*it != ECov::UNKNOWN && *it != ECov::FUNCTION)
    {
      ACovFunc* cova = createCovFunc(*it, ctxt);
      String cn = toUpper(cov_name);
      String ccn = toUpper(cova->getCovName());
      delete cova;
      if (cn == ccn)
        return *it;
    }
    it.toNext();
  }
  messerr("Unknown covariance name:%s!", cov_name.c_str());
  displayCovList(ctxt);
  return ECov::UNKNOWN;
}

double CovFactory::getScaleFactor(const ECov &type, double param)
{
  CovContext ctxt = CovContext(1, 1);
  ACovFunc *cova = CovFactory::createCovFunc(type, ctxt);
  cova->setParam(param);
  double scadef = cova->getScadef();
  delete cova;
  return scadef;
}
