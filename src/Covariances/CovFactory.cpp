/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "Covariances/CovFactory.hpp"
#include "geoslib_f.h"

#include <iostream>

#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Basic/Vector.hpp"
#include "Basic/AException.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovBesselJ.hpp"
#include "Covariances/CovBesselK.hpp"
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
#include "Covariances/CovP8.hpp"
#include "Covariances/CovPenta.hpp"
#include "Covariances/CovPower.hpp"
#include "Covariances/CovReg1D.hpp"
#include "Covariances/CovSincard.hpp"
#include "Covariances/CovSpherical.hpp"
#include "Covariances/CovStable.hpp"
#include "Covariances/CovStorkey.hpp"
#include "Covariances/CovTriangle.hpp"
#include "Covariances/CovWendland1.hpp"
#include "Covariances/CovWendland2.hpp"

bool _isValid(ACovFunc* cova, const CovContext& ctxt)
{
  if ((int) cova->getMaxNDim() > 0   &&
      (int) ctxt.getNDim() > (int) cova->getMaxNDim())  return false;
  if ((int) cova->getMinOrder() >= 0 &&
      ctxt.getIrfMaxDegree() < (int) cova->getMinOrder()) return false;
  return true;
}

ACovFunc* CovFactory::createCovFunc(const ENUM_COVS& type, const CovContext& ctxt)
{
  switch(type)
  {
    case COV_NUGGET:      return new CovNugget(ctxt);
    case COV_EXPONENTIAL: return new CovExponential(ctxt);
    case COV_SPHERICAL:   return new CovSpherical(ctxt);
    case COV_GAUSSIAN:    return new CovGaussian(ctxt);
    case COV_CUBIC:       return new CovCubic(ctxt);
    case COV_SINCARD:     return new CovSincard(ctxt);
    case COV_BESSEL_J:    return new CovBesselJ(ctxt);
    case COV_BESSEL_K:    return new CovBesselK(ctxt);
    case COV_GAMMA:       return new CovGamma(ctxt);
    case COV_CAUCHY:      return new CovCauchy(ctxt);
    case COV_STABLE:      return new CovStable(ctxt);
    case COV_LINEAR:      return new CovLinear(ctxt);
    case COV_POWER:       return new CovPower(ctxt);
    case COV_ORDER1_GC:   return new CovGC1(ctxt);
    case COV_SPLINE_GC:   return new CovGCspline(ctxt);
    case COV_SPLINE2_GC:  return new CovGCspline2(ctxt);
    case COV_ORDER3_GC:   return new CovGC3(ctxt);
    case COV_ORDER5_GC:   return new CovGC5(ctxt);
    case COV_COSINUS:     return new CovCosinus(ctxt);
    case COV_TRIANGLE:    return new CovTriangle(ctxt);
    case COV_COSEXP:      return new CovCosExp(ctxt);
    case COV_REG1D:       return new CovReg1D(ctxt);
    case COV_PENTA:       return new CovPenta(ctxt);
    case COV_STORKEY:     return new CovStorkey(ctxt);
    case COV_WENDLAND1:   return new CovWendland1(ctxt);
    case COV_WENDLAND2:   return new CovWendland2(ctxt);
    case COV_P8:          return new CovP8(ctxt);
    default:            std::cout << "Error unknown covariance !" << std::endl; break;
  }
  my_throw ("Covariance function not yet implemented !");
  return nullptr;
}

ACovFunc* CovFactory::duplicateCovFunc(const ACovFunc& cov)
{
  switch(cov.getType())
  {
    // Warning : if a crash with "bad cast" occurs, please check the type of your CovFunc
    case COV_NUGGET:      return new CovNugget(     dynamic_cast<const CovNugget&>     (cov));
    case COV_EXPONENTIAL: return new CovExponential(dynamic_cast<const CovExponential&>(cov));
    case COV_SPHERICAL:   return new CovSpherical(  dynamic_cast<const CovSpherical&>  (cov));
    case COV_GAUSSIAN:    return new CovGaussian(   dynamic_cast<const CovGaussian&>   (cov));
    case COV_CUBIC:       return new CovCubic(      dynamic_cast<const CovCubic&>      (cov));
    case COV_SINCARD:     return new CovSincard(    dynamic_cast<const CovSincard&>    (cov));
    case COV_BESSEL_J:    return new CovBesselJ(    dynamic_cast<const CovBesselJ&>    (cov));
    case COV_BESSEL_K:    return new CovBesselK(    dynamic_cast<const CovBesselK&>    (cov));
    case COV_GAMMA:       return new CovGamma(      dynamic_cast<const CovGamma&>      (cov));
    case COV_CAUCHY:      return new CovCauchy(     dynamic_cast<const CovCauchy&>     (cov));
    case COV_STABLE:      return new CovStable(     dynamic_cast<const CovStable&>     (cov));
    case COV_LINEAR:      return new CovLinear(     dynamic_cast<const CovLinear&>     (cov));
    case COV_POWER:       return new CovPower(      dynamic_cast<const CovPower&>      (cov));
    case COV_ORDER1_GC:   return new CovGC1(        dynamic_cast<const CovGC1&>        (cov));
    case COV_SPLINE_GC:   return new CovGCspline(   dynamic_cast<const CovGCspline&>   (cov));
    case COV_ORDER3_GC:   return new CovGC3(        dynamic_cast<const CovGC3&>        (cov));
    case COV_ORDER5_GC:   return new CovGC5(        dynamic_cast<const CovGC5&>        (cov));
    case COV_COSINUS:     return new CovCosinus(    dynamic_cast<const CovCosinus&>    (cov));
    case COV_TRIANGLE:    return new CovTriangle(   dynamic_cast<const CovTriangle&>   (cov));
    case COV_COSEXP:      return new CovCosExp(     dynamic_cast<const CovCosExp&>     (cov));
    case COV_REG1D:       return new CovReg1D(      dynamic_cast<const CovReg1D&>      (cov));
    case COV_PENTA:       return new CovPenta(      dynamic_cast<const CovPenta&>      (cov));
    case COV_SPLINE2_GC:  return new CovGCspline2(  dynamic_cast<const CovGCspline2&>  (cov));
    case COV_STORKEY:     return new CovStorkey(    dynamic_cast<const CovStorkey&>    (cov));
    case COV_WENDLAND1:   return new CovWendland1(  dynamic_cast<const CovWendland1&>  (cov));
    case COV_WENDLAND2:   return new CovWendland2(  dynamic_cast<const CovWendland2&>  (cov));
    case COV_P8:          return new CovP8(         dynamic_cast<const CovP8&>         (cov));
    default:              break;
  }
  message("Covariance type = %d\n",cov.getType());
  my_throw ("Covariance function not yet implemented !");
  return nullptr;
}

/**
 * Prints the list of Covariance available
 * @param ctxt  Pointer to CovContext structure
 */
void CovFactory::displayList(const CovContext& ctxt)
{
  message("List of authorized covariance / variogram names :\n");
  int navail = 0;
  for (int i = 0; i < COV_NUMBER; i++)
  {
    ACovFunc* cova = createCovFunc((ENUM_COVS) i, ctxt);
    if (_isValid(cova, ctxt))
      message("%2d - %s\n", navail + 1, cova->getCovName().c_str());
    delete cova;
    navail++;
  }
}

VectorString CovFactory::getCovList(const CovContext& ctxt)
{
  VectorString names;

  for (int i = 0; i < COV_NUMBER; i++)
  {
    ACovFunc* cova = createCovFunc((ENUM_COVS) i, ctxt);
    if (_isValid(cova, ctxt))
      names.push_back(cova->getCovName());
    delete cova;
  }
  return names;
}

int CovFactory::identifyCovariance(const String& cov_name,
                                   ENUM_COVS *rank,
                                   const CovContext& ctxt)
{
  for (int i = 0; i < COV_NUMBER; i++)
  {
    ACovFunc* cova = CovFactory::createCovFunc((ENUM_COVS) i, ctxt);
    if (matchRegexp(cova->getCovName(), cov_name, false))
    {
        *rank = (ENUM_COVS) i;
        return 0;
    }
    delete cova;
  }
  messerr("Unknown Covariance name : %s",cov_name.c_str());
  CovFactory::displayList(ctxt);
  return 1;
}

int CovFactory::getCovarianceNumber()
{
  return COV_NUMBER;
}
