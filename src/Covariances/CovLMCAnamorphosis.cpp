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

CovLMCAnamorphosis::CovLMCAnamorphosis(const EAnam& anam_type,
                                       int anam_nclass,
                                       int anam_iclass,
                                       int anam_var,
                                       double anam_coefr,
                                       double anam_coefs,
                                       VectorDouble& anam_strcnt,
                                       VectorDouble& anam_stats,
                                       const ASpace* space)
    : CovLMC(space),
      _anamType(),
      _anamNClass(0),
      _anamIClass(0),
      _anamPointBlock(0),
      _anamStrCount(),
      _anamMeans(),
      _anam(nullptr)
{
 init(anam_type,anam_nclass,anam_iclass,anam_var,anam_coefr,anam_coefs,anam_strcnt,anam_stats);
}

CovLMCAnamorphosis::CovLMCAnamorphosis(const CovLMCAnamorphosis &r)
    : CovLMC(r),
      _anamType(r._anamType),
      _anamNClass(r._anamNClass),
      _anamIClass(r._anamIClass),
      _anamPointBlock(r._anamPointBlock),
      _anamStrCount(r._anamStrCount),
      _anamMeans(r._anamMeans),
      _anam(r._anam)
{
}

CovLMCAnamorphosis& CovLMCAnamorphosis::operator=(const CovLMCAnamorphosis &r)
{
  if (this != &r)
  {
    CovLMC::operator=(r);
    _anamType = r._anamType;
    _anamNClass = r._anamNClass;
    _anamIClass = r._anamIClass;
    _anamPointBlock = r._anamPointBlock;
    _anamStrCount = r._anamStrCount;
    _anamMeans = r._anamMeans;
    _anam = r._anam;
  }
  return *this;
}

CovLMCAnamorphosis::~CovLMCAnamorphosis()
{
}

int CovLMCAnamorphosis::init(const EAnam& anam_type,
                             int anam_nclass,
                             int anam_iclass,
                             int anam_var,
                             double anam_coefr,
                             double anam_coefs,
                             VectorDouble& anam_strcnt,
                             VectorDouble& anam_stats)
{
  /* Preliminary checks */

  if (! (anam_iclass == 0 || anam_iclass < anam_nclass))
  {
    messerr("The rank of the active factor (%d) is incorrect",anam_iclass);
    messerr("It should lie between 1 and the number of factors (%d)",
            anam_nclass-1);
    messerr("or be set to 0 to estimate the whole discretized grade");
    return 1;
  }

  /* Load the parameters */

  if (anam_type == EAnam::HERMITIAN)
  {
    _anam = new AnamHermite();
    AnamHermite* anam_hermite = dynamic_cast<AnamHermite*>(_anam);
    anam_hermite->setRCoef(anam_coefr);
  }
  else if (anam_type == EAnam::DISCRETE_IR)
  {
    _anam = new AnamDiscreteIR();
    AnamDiscreteIR* anam_discrete_IR = dynamic_cast<AnamDiscreteIR*>(_anam);
    anam_discrete_IR->setNCut(anam_nclass);
    anam_discrete_IR->setRCoef(anam_coefr);
    anam_discrete_IR->setStats(anam_stats);
  }
  else if (anam_type == EAnam::DISCRETE_DD)
  {
    _anam = new AnamDiscreteDD();
    AnamDiscreteDD* anam_discrete_DD = dynamic_cast<AnamDiscreteDD*>(_anam);
    anam_discrete_DD->setNCut(anam_nclass);
    anam_discrete_DD->setSCoef(anam_coefs);
    anam_discrete_DD->setStats(anam_stats);
  }
  else
  {
    messerr("Unknown Anamorphosis type int Definition of Model Transformation");
    return 1;
  }

  _anamType = anam_type;
  _anamIClass = anam_iclass;
  _anamNClass = anam_nclass;
  _anamPointBlock = FFFF(anam_var) ? 0 : anam_var - 1;

  if (!anam_strcnt.empty())
  {
    _anamStrCount.resize(anam_nclass - 1);
    for (int i = 0; i < anam_nclass - 1; i++)
      _anamStrCount[i] = anam_strcnt[i];
  }

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
  bool asVario = mode.getAsVario(); // TODO: Pourquoi faire (voir ligne suivante)
  modeloc.setAsVario(false);

  int anam_var = getAnamPointBlock();
  /* 'anam_var' is negative if model evaluation is called from dk() */
  /* This modification is performed in model_anamorphosis_set_factor() */
  // TODO : this must be checked! anam_vario always equal to LHS, RHS or VAR ???
  if (anam_var >= 0) modeloc.setMember(ECalcMember::fromValue(anam_var));

  if (_anamType == EAnam::HERMITIAN)
  {
    AnamHermite *anam_hermite = dynamic_cast<AnamHermite*>(_anam);
    return _evalHermite(anam_hermite, ivar, jvar, p1, p2, mode);
  }

  if (_anamType == EAnam::DISCRETE_DD)
  {
    AnamDiscreteDD *anam_discreteDD = dynamic_cast<AnamDiscreteDD*>(_anam);
    return _evalDiscreteDD(anam_discreteDD, ivar, jvar, p1, p2, mode);
  }

  if (_anamType == EAnam::DISCRETE_IR)
  {
    AnamDiscreteIR *anam_discreteIR = dynamic_cast<AnamDiscreteIR*>(_anam);
    return _evalDiscreteIR(anam_discreteIR, ivar, jvar, p1, p2, mode);
  }

  return TEST;
}

double CovLMCAnamorphosis::_evalHermite(AnamHermite *anam,
                                        int ivar,
                                        int jvar,
                                        const SpacePoint& p1,
                                        const SpacePoint& p2,
                                        const CovCalcMode& mode) const
{
  double rho, coeff, psin2, rn, rhon;

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

    if (dist2 <= 0.)
    {
      rn = 1.;
      for (int iclass = 1; iclass < getAnamNClass(); iclass++)
      {
        psin2 = getAnamMeans(iclass) * getAnamMeans(iclass);
        rn *= anam->getRCoef();
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
      for (int iclass = 1; iclass < getAnamNClass(); iclass++)
      {
        psin2 = getAnamMeans(iclass) * getAnamMeans(iclass);
        rn *= anam->getRCoef();
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

    if (dist2 <= 0.)
    {
      rn = pow(anam->getRCoef(), (double) getAnamIClass());
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
      rn = pow(anam->getRCoef(), (double) getAnamIClass());
      rhon = pow(rho, (double) getAnamIClass());
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

double CovLMCAnamorphosis::_evalDiscreteDD(AnamDiscreteDD* anam,
                                           int ivar,
                                           int jvar,
                                           const SpacePoint& p1,
                                           const SpacePoint& p2,
                                           const CovCalcMode& mode) const
{
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
      for (int iclass = 1; iclass < getAnamNClass(); iclass++)
      {
        csi = anam->getDDStatCnorm(iclass);
        mui = anam->getDDStatMul(iclass);
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
      for (int iclass = 1; iclass < getAnamNClass(); iclass++)
      {
        li = anam->getDDStatLambda(iclass);
        csi = anam->getDDStatCnorm(iclass);
        mui = anam->getDDStatMul(iclass);
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

    if (dist2 <= 0.)
    {
      mui = anam->getDDStatMul(getAnamIClass());
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
      mui = anam->getDDStatMul(getAnamIClass());
      li = anam->getDDStatLambda(getAnamIClass());
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

double CovLMCAnamorphosis::_evalDiscreteIR(AnamDiscreteIR* anam,
                                           int ivar,
                                           int jvar,
                                           const SpacePoint& p1,
                                           const SpacePoint& p2,
                                           const CovCalcMode& mode) const
{
  CovCalcMode modeloc(mode);

  /* Initializations */

  int nclass = getAnamNClass();
  int ncut = nclass - 1;
  double r = anam->getRCoef();

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
      double bi = anam->getIRStatB(icut);
      cov2 = pow(_st_covsum_residual(anam, modeloc, icut, ivar, jvar, p1, p2), r);
      cov += bi * bi * (cov2 - cov1);
    }
  }
  else
  {
    /* Check if the distance is zero */


    /* Structure for the factor 'modtrs.anam_iclass' */

    if (dist2 <= 0)
    {
      cov = anam->getIRStatR(icut0 + 1);
    }
    else
    {
      double cov1 = pow(
          _st_covsum_residual(anam, modeloc, icut0 - 1, ivar, jvar, p1, p2), r);
      double cov2 = pow(
          _st_covsum_residual(anam, modeloc, icut0, ivar, jvar, p1, p2), r);
      cov = cov2 - cov1;
    }
  }
  return cov;
}

double CovLMCAnamorphosis::_st_cov_residual(AnamDiscreteIR *anam,
                                            CovCalcMode& mode,
                                            int icut0,
                                            int ivar,
                                            int jvar,
                                            const SpacePoint& p1,
                                            const SpacePoint& p2) const
{

  /* Get the pointer of the basic structure for the current model */

  int icov = 0;
  for (int icut = 0; icut < icut0; icut++)
    icov += (int) getAnamStrCount()[icut];

  /* Loop of the covariance basic structures */

  double cov = 0.;
  int number = (int) getAnamStrCount()[icut0];
  for (int i = 0; i < number; i++, icov++)
  {
    mode.setKeepOnlyCovIdx(i);
    double covloc = CovLMC::eval(ivar, jvar, p1, p2, mode);
    cov += covloc;
  }
  cov *= anam->getIRStatR(icut0 + 1);

  return cov;
}

double CovLMCAnamorphosis::_st_covsum_residual(AnamDiscreteIR* anam,
                                               CovCalcMode& mode,
                                               int icut0,
                                               int ivar,
                                               int jvar,
                                               const SpacePoint& p1,
                                               const SpacePoint& p2) const
{
  double covsum = 0.;
  for (int icut = 0; icut <= icut0; icut++)
    covsum += _st_cov_residual(anam, mode, icut, ivar, jvar, p1, p2);
  return (1. + covsum);
}

int CovLMCAnamorphosis::setAnamIClass(int anam_iclass)
{
  if (! (anam_iclass == 0 || anam_iclass < getAnamNClass()))
  {
    messerr("The rank of the active factor (%d) is incorrect", anam_iclass);
    messerr("It should lie between 1 and the number of factors (%d)", getAnamNClass() - 1);
    messerr("or be set to 0 to estimate the whole discretized grade");
    return 1;
  }
  return 0;
}
