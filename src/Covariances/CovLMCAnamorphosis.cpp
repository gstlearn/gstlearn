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
#include "Covariances/CovLMCAnamorphosis.hpp"

#include "Space/ASpace.hpp"
#include "Basic/AException.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/EConvType.hpp"
#include "Covariances/EConvDir.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Anamorphosis/AnamDiscreteDD.hpp"
#include "Anamorphosis/EAnam.hpp"
#include "geoslib_f.h"

#include <math.h>

CovLMCAnamorphosis::CovLMCAnamorphosis(const AAnam* anam,
                                       int anam_iclass,
                                       int anam_var,
                                       const VectorInt& strcnt,
                                       const ASpace* space)
    : CovLMC(space),
      _anamIClass(0),
      _anamPointBlock(0),
      _anamStrCount(),
      _anam(anam)
{
  init(anam_iclass,anam_var,strcnt);
}

CovLMCAnamorphosis::CovLMCAnamorphosis(const CovLMCAnamorphosis &r)
    : CovLMC(r),
      _anamIClass(r._anamIClass),
      _anamPointBlock(r._anamPointBlock),
      _anamStrCount(r._anamStrCount),
      _anam(r._anam)
{
}

CovLMCAnamorphosis& CovLMCAnamorphosis::operator=(const CovLMCAnamorphosis &r)
{
  if (this != &r)
  {
    CovLMC::operator=(r);
    _anamIClass = r._anamIClass;
    _anamPointBlock = r._anamPointBlock;
    _anamStrCount = r._anamStrCount;
    _anam = r._anam;
  }
  return *this;
}

CovLMCAnamorphosis::~CovLMCAnamorphosis()
{
}

int CovLMCAnamorphosis::init(int anam_iclass,
                             int anam_var,
                             const VectorInt& anam_strcnt)
{
  if (_anam == nullptr)
  {
    messerr("You must define 'anam'");
    return 1;
  }
  if (! (anam_iclass == 0 || anam_iclass < _anam->getNFactor()))
  {
    messerr("The rank of the active factor (%d) is incorrect",anam_iclass);
    messerr("It should lie between 1 and the number of factors (%d)",
            _anam->getNFactor()-1);
    messerr("or be set to 0 to estimate the whole discretized grade");
    return 1;
  }

  EAnam type = _anam->getType();
  if (type != EAnam::HERMITIAN && type != EAnam::DISCRETE_IR && type != EAnam::DISCRETE_DD)
  {
    messerr("Unknown Anamorphosis Definition of Model Transformation");
    return 1;
  }

  _anamIClass = anam_iclass;
  _anamPointBlock = FFFF(anam_var) ? 0 : anam_var - 1;
  _anamStrCount = anam_strcnt;

  return 0;
}

String CovLMCAnamorphosis::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  sstr << ACovAnisoList::toString(strfmt);

  sstr << _anam->toString(strfmt);

  return sstr.str();
}

double CovLMCAnamorphosis::eval0(int /*ivar*/,
                                 int /*jvar*/,
                                 const CovCalcMode& /*mode*/) const
{
  double cov0 = 0.;
  return cov0;
}

double CovLMCAnamorphosis::eval(int ivar,
                                int jvar,
                                const SpacePoint& p1,
                                const SpacePoint& p2,
                                const CovCalcMode& mode) const
{
  if (_anam == nullptr) return TEST;

  // The calculation flag 'as.Vario' must be treated here rather than relying on calculation
  // performed internally in 'eval' function
  CovCalcMode modeloc(mode);
  modeloc.setAsVario(false);

  int anam_var = getAnamPointBlock();
  /* 'anam_var' is negative if model evaluation is called from dk() */
  /* This modification is performed in model_anamorphosis_set_factor() */
  // TODO : this must be checked! anam_vario always equal to LHS, RHS or VAR ???
  if (anam_var >= 0) modeloc.setMember(ECalcMember::fromValue(anam_var));

  EAnam type = _anam->getType();
  if (type == EAnam::HERMITIAN)
  {
    return _evalHermite(ivar, jvar, p1, p2, mode);
  }

  if (type == EAnam::DISCRETE_DD)
  {
    return _evalDiscreteDD(ivar, jvar, p1, p2, mode);
  }

  if (type == EAnam::DISCRETE_IR)
  {
    return _evalDiscreteIR(ivar, jvar, p1, p2, mode);
  }

  return TEST;
}

double CovLMCAnamorphosis::_evalHermite(int ivar,
                                        int jvar,
                                        const SpacePoint& p1,
                                        const SpacePoint& p2,
                                        const CovCalcMode& mode) const
{
  double rho, coeff, psin2, rn, rhon;
  const AnamHermite *anamH = dynamic_cast<const AnamHermite*>(_anam);


  /* Check if the distance is zero */

  rho = coeff = 0.;
  double dist2 = getDistance(p1, p2);
  double cov = 0.;
  if (dist2 > 0) cov = CovLMC::eval(ivar, jvar, p1, p2, mode);

  /* Update the covariance */

  if (getAnamIClass() == 0)
  {

    /*********************************************/
    /* Structure for the whole discretized grade */
    /*********************************************/

    psin2 = 1.;
    if (dist2 <= 0.)
    {
      rn = 1.;
      for (int iclass = 1; iclass < anamH->getNFactor(); iclass++)
      {
//        psin2 = getAnamMeans(iclass) * getAnamMeans(iclass);
        rn *= anamH->getRCoef();
        switch (mode.getMember().toEnum())
        {
          case ECalcMember::E_LHS:
            coeff = 1. / (rn * rn);
            break;
          case ECalcMember::E_RHS:
            coeff = 1. / rn;
            break;
          case ECalcMember::E_VAR:
            coeff = 1.;
            break;
        }
        cov += coeff * psin2;
      }
    }
    else
    {
      rn = 1.;
      rhon = 1.;
      for (int iclass = 1; iclass < anamH->getNFactor(); iclass++)
      {
//        psin2 = getAnamMeans(iclass) * getAnamMeans(iclass);
        rn *= anamH->getRCoef();
        rhon *= rho;
        switch (mode.getMember().toEnum())
        {
          case ECalcMember::E_LHS:
            coeff = 1.;
            break;

          case ECalcMember::E_RHS:
            coeff = 1. / rn;
            break;

          case ECalcMember::E_VAR:
            coeff = 1.;
            break;
        }
        cov += coeff * psin2 * rhon;
      }
    }
  }
  else
  {

    /**************************************************/
    /* Structure for the factor 'modtrs.anam_iclass' */
    /**************************************************/

    int iclass = getAnamIClass();
    if (dist2 <= 0.)
    {
      rn = pow(anamH->getRCoef(), (double) iclass);
      switch (mode.getMember().toEnum())
      {
        case ECalcMember::E_LHS:
          coeff = 1.;
          break;

        case ECalcMember::E_RHS:
          coeff = rn;
          break;

        case ECalcMember::E_VAR:
          coeff = 1.;
          break;
      }
      cov = coeff;
    }
    else
    {
      rn = pow(anamH->getRCoef(), (double) iclass);
      rhon = pow(rho, (double) iclass);
      switch (mode.getMember().toEnum())
      {
        case ECalcMember::E_LHS:
          coeff = rn * rn;
          break;
        case ECalcMember::E_RHS:
          coeff = rn;
          break;
        case ECalcMember::E_VAR:
          coeff = 1.;
          break;
      }
      cov = coeff * rhon;
    }
  }
  return cov;
}

double CovLMCAnamorphosis::_evalDiscreteDD(int ivar,
                                           int jvar,
                                           const SpacePoint& p1,
                                           const SpacePoint& p2,
                                           const CovCalcMode& mode) const
{
  const AnamDiscreteDD *anamDD = dynamic_cast<const AnamDiscreteDD*>(_anam);
  double gamref, csi, li, mui, coeff;

  /* Check if the distance is zero */

  double cov = 0.;
  gamref = coeff = 0.;
  double dist2 = getDistance(p1, p2);
  if (dist2 > 0)
  {
    double cov1 = CovLMC::eval(ivar, jvar, p1, p2, mode);
    double cov0 = CovLMC::eval0(ivar, jvar, mode);
    gamref = cov0 - cov1;
  }

  /* Update the covariance */

  if (getAnamIClass() == 0)
  {

    /*********************************************/
    /* Structure for the whole discretized grade */
    /*********************************************/

    if (dist2 <= 0.)
      for (int iclass = 1; iclass < anamDD->getNClass(); iclass++)
      {
        csi = anamDD->getDDStatCnorm(iclass);
        mui = anamDD->getDDStatMul(iclass);
        switch (mode.getMember().toEnum())
        {
          case ECalcMember::E_LHS:
            coeff = csi * csi;
            break;
          case ECalcMember::E_RHS:
            coeff = csi * csi / mui;
            break;
          case ECalcMember::E_VAR:
            coeff = csi * csi;
            break;
        }
        cov += coeff;
      }
    else
      for (int iclass = 1; iclass < anamDD->getNClass(); iclass++)
      {
        li  = anamDD->getDDStatLambda(iclass);
        csi = anamDD->getDDStatCnorm(iclass);
        mui = anamDD->getDDStatMul(iclass);
        switch (mode.getMember().toEnum())
        {
          case ECalcMember::E_LHS:
            coeff = csi * csi;
            break;

          case ECalcMember::E_RHS:
            coeff = csi * csi / mui;
            break;

          case ECalcMember::E_VAR:
            coeff = csi * csi;
            break;
        }
        cov += coeff * exp(-li * gamref);
      }
  }
  else
  {

    /**************************************************/
    /* Structure for the factor 'modtrs.anam_iclass' */
    /**************************************************/

    int iclass = getAnamIClass();
    if (dist2 <= 0.)
    {
      mui = anamDD->getDDStatMul(iclass);
      switch (mode.getMember().toEnum())
      {
        case ECalcMember::E_LHS:
          coeff = 1.;
          break;

        case ECalcMember::E_RHS:
          coeff = mui;
          break;

        case ECalcMember::E_VAR:
          coeff = 1.;
          break;
      }
      cov = coeff;
    }
    else
    {
      mui = anamDD->getDDStatMul(iclass);
      li  = anamDD->getDDStatLambda(iclass);
      switch (mode.getMember().toEnum())
      {
        case ECalcMember::E_LHS:
          coeff = mui * mui;
          break;
        case ECalcMember::E_RHS:
          coeff = mui;
          break;
        case ECalcMember::E_VAR:
          coeff = 1.;
          break;
      }
      cov = coeff * exp(-li * gamref);
    }
  }
  return cov;
}

double CovLMCAnamorphosis::_evalDiscreteIR(int ivar,
                                           int jvar,
                                           const SpacePoint& p1,
                                           const SpacePoint& p2,
                                           const CovCalcMode& mode) const
{
  const AnamDiscreteIR *anamIR = dynamic_cast<const AnamDiscreteIR*>(_anam);
  CovCalcMode modeloc(mode);

  /* Initializations */

  int nclass = anamIR->getNClass();
  int ncut = nclass - 1;
  double r = anamIR->getRCoef();

  /* Calculate the generic variogram value */

  double cov = 0.;
  double dist2 = getDistance(p1, p2);
  if (dist2 > 0)
    cov = CovLMC::eval(ivar, jvar, p1, p2, mode);

  /* Modification of the covariance */

  int icut0 = getAnamIClass() - 1;
  if (icut0 < 0)
  {

    /* Structure for the whole discretized grade */

    double cov2 = 1.;
    for (int icut = 0; icut < ncut; icut++)
    {
      double cov1 = cov2;
      double bi = anamIR->getIRStatB(icut);
      cov2 = pow(__covSumResidualIR(modeloc, icut, ivar, jvar, p1, p2), r);
      cov += bi * bi * (cov2 - cov1);
    }
  }
  else
  {
    /* Check if the distance is zero */


    /* Structure for the factor 'modtrs.anam_iclass' */

    if (dist2 <= 0)
    {
      cov = anamIR->getIRStatR(icut0 + 1);
    }
    else
    {
      double cov1 = pow(__covSumResidualIR(modeloc, icut0 - 1, ivar, jvar, p1, p2), r);
      double cov2 = pow(__covSumResidualIR(modeloc, icut0, ivar, jvar, p1, p2), r);
      cov = cov2 - cov1;
    }
  }
  return cov;
}

double CovLMCAnamorphosis::_covResidualIR(CovCalcMode& mode,
                                            int icut0,
                                            int ivar,
                                            int jvar,
                                            const SpacePoint& p1,
                                            const SpacePoint& p2) const
{
  const AnamDiscreteIR *anamIR = dynamic_cast<const AnamDiscreteIR*>(_anam);

  /* Get the pointer of the basic structure for the current model */

  int icov = 0;
  for (int icut = 0; icut < icut0; icut++)
    icov += _anamStrCount[icut];

  /* Loop of the covariance basic structures */

  double cov = 0.;
  int number = _anamStrCount[icut0];
  for (int i = 0; i < number; i++, icov++)
  {
    mode.setKeepOnlyCovIdx(i);
    double covloc = CovLMC::eval(ivar, jvar, p1, p2, mode);
    cov += covloc;
  }
  cov *= anamIR->getIRStatR(icut0 + 1);

  return cov;
}

double CovLMCAnamorphosis::__covSumResidualIR(CovCalcMode& mode,
                                               int icut0,
                                               int ivar,
                                               int jvar,
                                               const SpacePoint& p1,
                                               const SpacePoint& p2) const
{
  double covsum = 0.;
  for (int icut = 0; icut <= icut0; icut++)
    covsum += _covResidualIR(mode, icut, ivar, jvar, p1, p2);
  return (1. + covsum);
}

int CovLMCAnamorphosis::setAnamIClass(int anam_iclass)
{
  if (! (anam_iclass == 0 || anam_iclass < _anam->getNFactor()))
  {
    messerr("The rank of the active factor (%d) is incorrect", anam_iclass);
    messerr("It should lie between 1 and the number of factors (%d)", _anam->getNFactor() - 1);
    messerr("or be set to 0 to estimate the whole discretized grade");
    return 1;
  }
  return 0;
}

const EAnam CovLMCAnamorphosis::getAnamType() const
{
  if (_anam == nullptr) return EAnam::UNKNOWN;
  return _anam->getType();
}
