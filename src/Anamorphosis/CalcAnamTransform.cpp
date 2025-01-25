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
#include "Calculators/ACalcDbVarCreator.hpp"
#include "geoslib_define.h"
#include "geoslib_old_f.h"

#include "Db/Db.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Law.hpp"
#include "Anamorphosis/AAnam.hpp"
#include "Anamorphosis/CalcAnamTransform.hpp"
#include "Anamorphosis/AnamContinuous.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Anamorphosis/AnamDiscreteDD.hpp"
#include "Polynomials/Hermite.hpp"
#include "Polynomials/MonteCarlo.hpp"

#include <math.h>

CalcAnamTransform::CalcAnamTransform(AAnam* anam)
    : ACalcDbVarCreator(),
      _iattVar(-1),
      _iattFac(-1),
      _iattSel(-1),
      _flagVars(false),
      _flagToFactors(false),
      _flagDisjKrig(false),
      _flagCondExp(false),
      _flagUniCond(false),
      _flagZToY(true),
      _flagNormalScore(false),
      _ifacs(),
      _iptrEst(),
      _iptrStd(),
      _nbsimu(0),
      _flagOK(false),
      _proba(TEST),
      _anam(anam),
      _selectivity(nullptr)
{
}

CalcAnamTransform::~CalcAnamTransform()
{
}

bool CalcAnamTransform::_hasAnam(const EAnam& anamType) const
{
  if (_anam == nullptr)
  {
    messerr("The argument 'anam' must be defined");
    return false;
  }
  if (anamType == EAnam::UNKNOWN) return true;
  if (anamType != _anam->getType())
  {
    messerr("The argument 'anam'  should be of type");
    return false;
  }
  return true;
}

bool CalcAnamTransform::_hasInputVarDefined(int mode) const
{
  // Check that the vector is not empty
  if (_iptrEst.empty())
  {
    messerr("'db' should contain an Estimate variable");
    return false;
  }

  // Check that the UID are valid
  for (int i = 0; i < (int) _iptrEst.size(); i++)
    if (_iptrEst[i] < 0)
    {
      messerr("An estimation variable is not correctly defined");
      return false;
    }

  // Check that the vector is not empty
  if (_iptrStd.empty())
  {
    if (mode == 0)
      messerr("'db' should contain an St.Dev of Estimation Error variable");
    else
      messerr("'db' should contain an Variance of Estimation Error variable");
    return false;
  }

  // Check that the UID are valid
  for (int i = 0; i < (int) _iptrStd.size(); i++)
    if (_iptrStd[i] < 0)
    {
      if (mode == 0)
        messerr("A St. Dev. variable is not correctly defined");
      else
        messerr("A Variance variable is not correctly defined");
      return false;
    }

  return true;
}

bool CalcAnamTransform::_hasSelectivity() const
{
  int ncuts = _selectivity->getNCuts();
  if (ncuts <= 0 && ! _selectivity->isOnlyZDefined())
  {
    messerr("You must define some cutoff values");
    return false;
  }
  int nvarout = _selectivity->getNVar();
  if (nvarout <= 0)
  {
    messerr("No recovery function is defined");
    return false;
  }
  return true;
}

bool CalcAnamTransform::_hasVariableNumber(bool equal1) const
{
  int number = getDb()->getNLoc(ELoc::Z);
  if (! equal1)
  {
    if (number <= 0)
    {
      messerr("The argument 'db'  must have some variable(s) defined");
      return false;
    }
  }
  else
  {
    if (number != 1)
    {
      messerr("The argument 'db'  must have a single variable defined");
      return false;
    }
  }
  return true;
}

bool CalcAnamTransform::_check()
{
  if (!ACalcDbVarCreator::_check()) return false;

  if (!hasDb()) return false;
  if (! _hasAnam()) return false;
  if (! _hasVariableNumber()) return false;

  if (_flagVars)
  {
    AnamContinuous* anamC = dynamic_cast<AnamContinuous*>(_anam);
    if (anamC == nullptr)
    {
      messerr("The argument 'anam'  must be of type AnamContinuous");
      return false;
    }
    return true;
  }

  if (_flagToFactors)
  {
    if (! _hasVariableNumber(true)) return false;
    int nmax = _anam->getNFactor();
    int nfact = _getNfact();
    for (int ifac = 0; ifac < nfact; ifac++)
      if (_ifacs[ifac] < 1 || _ifacs[ifac] > nmax)
      {
        messerr("Error in the rank of the factor(%d): it should lie in [1,%d]",
                _ifacs[ifac], nmax);
        return false;
      }
    return true;
  }

  if (_flagDisjKrig)
  {
    if (! _hasAnam(EAnam::HERMITIAN)) return false;
    if (! _hasInputVarDefined()) return false;
    if (! _hasSelectivity()) return false;
    return true;
  }

  if (_flagCondExp)
  {
    if (! _hasAnam(EAnam::HERMITIAN)) return false;
    if (! _hasInputVarDefined()) return false;
    if (! _hasSelectivity()) return false;
    return true;
  }

  if (_flagUniCond)
  {
    if (! _hasAnam(EAnam::HERMITIAN)) return false;
    if (! _hasInputVarDefined(1)) return false;
    if (! _hasSelectivity()) return false;
    if (_selectivity->isUsed(ESelectivity::Z))
    {
      messerr("The recovery option 'Z' is not available in this function");
      return false;
    }
    return true;
  }

  messerr("No Transformation option has been defined");
  return false;
}

bool CalcAnamTransform::_preprocess()
{
  if (!ACalcDbVarCreator::_preprocess()) return false;
  
  if (_flagVars)
  {
    int nvar = _getNVar();
    _iattVar = getDb()->addColumnsByConstant(nvar);
    return (_iattVar >= 0);
  }

  if (_flagToFactors)
  {
    int nfact = _getNfact();
    _iattFac = getDb()->addColumnsByConstant(nfact);
    return (_iattFac >= 0);
  }

  if (_flagDisjKrig)
  {
    int nvarout = _getNSel();
    _iattSel = getDb()->addColumnsByConstant(nvarout, TEST);
    if (_iattSel < 0) return 1;
    return true;
  }

  if (_flagCondExp)
  {
    int nvarout = _getNSel();
    _iattSel = getDb()->addColumnsByConstant(nvarout, TEST);
    if (_iattSel < 0) return 1;
    return true;
  }

  if (_flagUniCond)
  {
    int nvarout = _getNSel();
    _iattSel = getDb()->addColumnsByConstant(nvarout, TEST);
    if (_iattSel < 0) return 1;
    return true;
  }

  return false;
}

bool CalcAnamTransform::_postprocess()
{
  /* Free the temporary variables */
  _cleanVariableDb(2);

  if (_flagVars)
  {
    int nvar = _getNVar();
    _renameVariable(nvar, _iattVar, ELoc::Z, String(), 1);
    return true;
  }

  if (_flagToFactors)
  {
    int nfact = _getNfact();
    _renameVariable(1, _iattFac, ELoc::Z, String(), nfact);
    return true;
  }

  if (_flagDisjKrig)
  {
    int nsel = _getNSel();
    for (int i = 0; i < nsel; i++)
      _renameVariable(1, _iattSel+i, ELoc::UNKNOWN, _selectivity->getVariableName(i), 1);
    return true;
  }

  if (_flagCondExp)
  {
    int nsel = _getNSel();
    for (int i = 0; i < nsel; i++)
      _renameVariable(1, _iattSel+i, ELoc::Z, _selectivity->getVariableName(i), 1);
    return true;
  }

  if (_flagUniCond)
  {
    int nsel = _getNSel();
    for (int i = 0; i < nsel; i++)
      _renameVariable(1, _iattSel+i, ELoc::Z, _selectivity->getVariableName(i), 1);
    return true;
  }

  return true;
}

void CalcAnamTransform::_rollback()
{
  _cleanVariableDb(1);
}

bool CalcAnamTransform::_run()
{
  if (_flagVars)
  {
    if (_flagZToY)
    {

      // Transform from Raw to Gaussian
      if (_flagNormalScore)
      {
        // Transform using Normal Score
        if (_ZToYByNormalScore()) return true;
      }
      else
      {

        // Transform using Hermite expansion
        if (_ZToYByHermite()) return true;
      }
    }
    else
    {

      // Transform from Gaussian to Raw
      if (_YToZByHermite()) return true;
    }
    return false;
  }

  if (_flagToFactors)
  {
    return _ZToFactors();
  }

  if (_flagDisjKrig)
  {
    return _FactorsToSelectivity();
  }

  if (_flagCondExp)
  {
    return (!_conditionalExpectation(getDb(), _anam, _selectivity, _iattSel,
                                 _iptrEst[0], _iptrStd[0], _flagOK, _proba,
                                 _nbsimu));
  }

  if (_flagUniCond)
  {
    AnamHermite* anam_hermite = dynamic_cast<AnamHermite*>(_anam);
    return (!_uniformConditioning(getDb(), anam_hermite, _selectivity, _iattSel,
                                  _iptrEst[0], _iptrStd[0]));
  }

  return false;
}

bool CalcAnamTransform::_ZToYByNormalScore()
{
  int nech = getDb()->getNSample();
  int nvar = _getNVar();

  // Reading the optional weight

  VectorDouble wt;
  if (getDb()->hasLocVariable(ELoc::W))
    wt = getDb()->getColumnByLocator(ELoc::W);
  else
    wt.resize(nech, 1.);

  // Loop on the variables

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    VectorDouble data = getDb()->getColumnByLocator(ELoc::Z, ivar);
    VectorDouble vec = VH::normalScore(data, wt);
    if (! vec.empty())
      getDb()->setColumnByUID(vec, _iattVar + ivar);
  }
  return true;
}

bool CalcAnamTransform::_ZToYByHermite()
{
  int nvar = _getNVar();
  const AnamContinuous* anamC = dynamic_cast<const AnamContinuous*>(getAnam());

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    VectorDouble z = getDb()->getColumnByLocator(ELoc::Z, ivar, true);
    if (z.size() <= 0) continue;
    VectorDouble y = anamC->rawToGaussianVector(z);
    getDb()->setColumnByUID(y, _iattVar + ivar, true);
  }
  return true;
}

bool CalcAnamTransform::_YToZByHermite()
{
  int nvar = _getNVar();
  const AnamContinuous* anamC = dynamic_cast<const AnamContinuous*>(getAnam());

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    VectorDouble y = getDb()->getColumnByLocator(ELoc::Z, ivar, true);
    if (y.size() <= 0) continue;
    VectorDouble z = anamC->gaussianToRawVector(y);
    getDb()->setColumnByUID(z, _iattVar + ivar, true);
  }
  return true;
}

bool CalcAnamTransform::_ZToFactors()
{
  int nfact = _getNfact();

  for (int iech = 0; iech < getDb()->getNSample(); iech++)
  {
    if (! getDb()->isActive(iech)) continue;
    double zval = getDb()->getZVariable(iech, 0);
    if (FFFF(zval)) continue;
    VectorDouble factors = _anam->z2factor(zval, _ifacs);
    if (factors.empty()) continue;
    for (int ifac = 0; ifac < nfact; ifac++)
      getDb()->setArray(iech, _iattFac + ifac, factors[ifac]);
  }
  return true;
}

bool CalcAnamTransform::_FactorsToSelectivity()
{
  AnamHermite    *anam_hermite     = dynamic_cast<AnamHermite*>(_anam);
  AnamDiscreteDD *anam_discrete_DD = dynamic_cast<AnamDiscreteDD*>(_anam);
  AnamDiscreteIR *anam_discrete_IR = dynamic_cast<AnamDiscreteIR*>(_anam);

  /* Dispatch according to the type of Anamorphosis */

  switch (_anam->getType().toEnum())
  {
    case EAnam::E_HERMITIAN:
      anam_hermite->factor2Selectivity(getDb(), _selectivity, _iptrEst,
                                       _iptrStd, _iattSel);
      return true;

    case EAnam::E_DISCRETE_DD:
      anam_discrete_DD->factor2Selectivity(getDb(), _selectivity, _iptrEst,
                                           _iptrStd, _iattSel);
      return true;

    case EAnam::E_DISCRETE_IR:
      anam_discrete_IR->factor2Selectivity(getDb(), _selectivity, _iptrEst,
                                           _iptrStd, _iattSel);
      return true;

    default:
      messerr("This method is not programmed yet for this anamorphosis");
      return false;
  }
}


/*****************************************************************************/
/*!
 **  Calculate the recoveries (z,T,Q,m,B) starting from the factors
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         Point anamorphosis
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  name_est     Array of variable names for factor estimation
 ** \param[in]  name_std     Array of variable names for factor St. Dev.
 ** \param[in]  namconv      Naming convention
 **
 *****************************************************************************/
int DisjunctiveKriging(Db *db,
                        AAnam *anam,
                        Selectivity *selectivity,
                        const VectorString &name_est,
                        const VectorString &name_std,
                        const NamingConvention &namconv)
{
  CalcAnamTransform transfo(anam);
  transfo.setDb(db);
  transfo.setSelectivity(selectivity);
  transfo.setIptrEst(db->getUIDs(name_est));
  transfo.setIptrStd(db->getUIDs(name_std));
  transfo.setFlagDisjKrig(true);
  transfo.setNamingConvention(namconv);

  // Run the calculator
  int error = (transfo.run()) ? 0 : 1;
  return error;
}

/*****************************************************************************/
/*!
 **  Calculate the Conditional Expectation
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         Point anamorphosis
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  name_est     Name of the Kriging estimate
 ** \param[in]  name_std     Name of the Kriging St. deviation
 ** \param[in]  flag_OK      1 if kriging has ben performed in Ordinary Kriging
 ** \param[in]  proba        Probability
 ** \param[in]  nbsimu       Number of Simulation outcomes
 ** \param[in]  namconv      Naming convention
 **
 *****************************************************************************/
int ConditionalExpectation(Db *db,
                           AAnam *anam,
                           Selectivity *selectivity,
                           const String &name_est,
                           const String &name_std,
                           bool flag_OK,
                           double proba,
                           int nbsimu,
                           const NamingConvention &namconv)
{
  CalcAnamTransform transfo(anam);
  transfo.setDb(db);
  transfo.setSelectivity(selectivity);
  transfo.setIptrEst({db->getUID(name_est)});
  transfo.setIptrStd({db->getUID(name_std)});
  transfo.setFlagCondExp(true);
  transfo.setFlagOk(flag_OK);
  transfo.setNbsimu(nbsimu);
  transfo.setProba(proba);
  transfo.setNamingConvention(namconv);

  // Run the calculator
  int error = (transfo.run()) ? 0 : 1;
  return error;
}

/*****************************************************************************/
/*!
 **  Calculate the Uniform Conditioning
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         Point anamorphosis
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  name_est     Name of the Kriging estimate
 ** \param[in]  name_varz    Name of the Variance of Kriging estimate
 ** \param[in]  namconv      Naming Convention
 **
 ** \remark We need the variance of Estimation Error... even if it will be
 ** \remark temporarily stored in a member names iptrStd.
 **
 *****************************************************************************/
int UniformConditioning(Db *db,
                        AAnam *anam,
                        Selectivity *selectivity,
                        const String &name_est,
                        const String &name_varz,
                        const NamingConvention &namconv)
{
  CalcAnamTransform transfo(anam);
  transfo.setDb(db);
  transfo.setSelectivity(selectivity);
  transfo.setIptrEst({db->getUID(name_est)});
  transfo.setIptrStd({db->getUID(name_varz)});
  transfo.setFlagUniCond(true);
  transfo.setNamingConvention(namconv);

  // Run the calculator
  int error = (transfo.run()) ? 0 : 1;
  return error;
}

/*****************************************************************************/
/*!
 **  Calculate the Conditional Expectation
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         Anamorphosis structure
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  iptr0        Rank of the pointer for storage
 ** \param[in]  col_est      Rank of variable containing Kriging estimate
 ** \param[in]  col_std      Rank of Variable containing Kriging St. deviation
 ** \param[in]  flag_OK      1 if kriging has ben performed in Ordinary Kriging
 ** \param[in]  proba        Probability
 ** \param[in]  nbsimu       Number of Simulation outcomes
 **
 *****************************************************************************/
int CalcAnamTransform::_conditionalExpectation(Db* db,
                                               AAnam* anam,
                                               const Selectivity* selectivity,
                                               int iptr0,
                                               int col_est,
                                               int col_std,
                                               bool flag_OK,
                                               double proba,
                                               int nbsimu)
{
  AnamHermite* anam_hermite = dynamic_cast<AnamHermite*>(anam);

  /* Analyzing the codes */

  VectorDouble ycuts = anam_hermite->rawToTransformVec(selectivity->getZcut());
  int need_T         = selectivity->isNeededT();

  /* Computing the estimation */

  if (selectivity->isUsed(ESelectivity::Z))
  {
    if (_ceZ(db, anam_hermite, selectivity, iptr0, col_est, col_std, nbsimu,
                flag_OK))
      return 1;
  }

  /* Compute Conditional Expectation for Tonnage */

  if (selectivity->isUsed(ESelectivity::T))
  {
    if (_ceT(1, db, selectivity, iptr0, col_est, col_std, ycuts, nbsimu,
                flag_OK))
      return 1;
  }

  /* Compute Conditional Expectation for Metal Quantity */

  if (selectivity->isUsed(ESelectivity::Q))
  {
    if (_ceQ(db, anam_hermite, selectivity, iptr0, col_est, col_std, ycuts,
                nbsimu, flag_OK))
      return 1;
  }

  /* Compute Conditional Expectation for Conventional Benefit */

  if (selectivity->isUsed(ESelectivity::B))
  {
    if (_ceB(db, selectivity, iptr0, ycuts)) return 1;
  }

  /* Compute Conditional Expectation for Average recoveable grade */

  if (selectivity->isUsed(ESelectivity::M))
  {
    if (_ceM(db, selectivity, iptr0)) return 1;
  }

  /* Compute Conditional Expectation for Tonnage */

  if (selectivity->isUsed(ESelectivity::PROP) && need_T)
  {
    if (_ceT(2, db, selectivity, iptr0, col_est, col_std, ycuts, nbsimu,
                flag_OK))
      return 1;
  }

  /* Compute Conditional Expectation for Quantile */

  if (selectivity->isUsed(ESelectivity::QUANT))
  {
    if (_ceQuant(db, anam_hermite, selectivity, iptr0, col_est, col_std,
                    proba, flag_OK))
      return 1;
  }

  if (!selectivity->isUsed(ESelectivity::T))
    db->deleteColumnsByUIDRange(
      selectivity->getAddressQTEst(ESelectivity::T, iptr0),
      selectivity->getNumberQTEst(ESelectivity::T));
  if (!selectivity->isUsed(ESelectivity::Q))
    db->deleteColumnsByUIDRange(
      selectivity->getAddressQTEst(ESelectivity::Q, iptr0),
      selectivity->getNumberQTEst(ESelectivity::Q));
  return 0;
}

/*****************************************************************************/
/*!
 **  Calculate the Uniform Conditioning
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         Block Hermite anamorphosis
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  iptr0        Pointer for storage
 ** \param[in]  col_est      Rank of variable containing Kriging estimate
 ** \param[in]  col_var      Rank of Variable containing Variance of Kriging
 *estimate
 **
 *****************************************************************************/
int CalcAnamTransform::_uniformConditioning(Db* db,
                                            AnamHermite* anam,
                                            Selectivity* selectivity,
                                            int iptr0,
                                            int col_est,
                                            int col_var)
{
  //  anam->setFlagBound(1);
  int nbpoly = anam->getNbPoly();
  int ncuts  = selectivity->getNCuts();

  /* Core allocation */

  VectorVectorDouble phi_b_zc(ncuts);

  /* Memorize the punctual Hermite polynomials */

  VectorDouble psi_hn = anam->getPsiHns();

  /* Add variables for storage */

  int iptr_sV = db->addColumnsByConstant(1, TEST);
  if (iptr_sV < 0) return 1;
  int iptr_yV = db->addColumnsByConstant(1, TEST);
  if (iptr_yV < 0) return 1;

  // Calculate the change of support coefficient

  double r_coef = anam->getRCoef();

  /* Transform zcuts into gaussian equivalent */

  VectorDouble ycuts = anam->rawToTransformVec(selectivity->getZcut());

  /* Fill the array phi_b_zc */

  for (int icut = 0; icut < ncuts; icut++)
  {
    VectorDouble hn = hermitePolynomials(ycuts[icut], 1., nbpoly);
    phi_b_zc[icut]  = hermiteCoefMetal(ycuts[icut], hn);
  }

  /* Computing S and Y on panels */

  double vv_min = 1.e30;
  double zv_min = 1.e30;
  double vv_max = -1.e30;
  double zv_max = -1.e30;
  double sv_min = 1.e30;
  double yv_min = 1.e30;
  double sv_max = -1.e30;
  double yv_max = -1.e30;

  for (int iech = 0; iech < db->getNSample(); iech++)
  {
    if (!db->isActive(iech)) continue;
    anam->setPsiHns(psi_hn);
    anam->calculateMeanAndVariance();
    double zvstar = db->getArray(iech, col_est);
    double varv   = db->getArray(iech, col_var);
    if (anamPointToBlock(anam, 0, varv, TEST, TEST)) continue;
    db->setArray(iech, iptr_sV, anam->getRCoef());
    db->setArray(iech, iptr_yV, anam->rawToTransformValue(zvstar));

    if (varv < vv_min) vv_min = varv;
    if (varv > vv_max) vv_max = varv;
    if (zvstar < zv_min) zv_min = zvstar;
    if (zvstar > zv_max) zv_max = zvstar;
  }

  /* Loop on the panels to compute the grade-tonnage functions */

  for (int iech = 0; iech < db->getNSample(); iech++)
  {
    if (!db->isActive(iech)) continue;
    double sv = db->getArray(iech, iptr_sV);
    double yv = db->getArray(iech, iptr_yV);

    /* Loop on the cutoffs */

    for (int icut = 0; icut < ncuts; icut++)
    {
      double yc = ycuts[icut];
      double ore =
        1. - law_cdf_gaussian((yc - yv * sv / r_coef) /
                              sqrt(1. - (sv / r_coef) * (sv / r_coef)));
      VectorDouble hn = hermitePolynomials(yv, 1., nbpoly);
      for (int ih = 0; ih < nbpoly; ih++)
        hn[ih] *= pow(sv / r_coef, (double)ih);
      double metal;
      matrix_product_safe(1, nbpoly, 1, hn.data(), phi_b_zc[icut].data(),
                          &metal);

      if (sv < sv_min) sv_min = sv;
      if (sv > sv_max) sv_max = sv;
      if (yv < yv_min) yv_min = yv;
      if (yv > yv_max) yv_max = yv;

      /* Storing the grade-tonnage functions */

      selectivity->setTest(icut, ore);
      selectivity->setQest(icut, metal);
    }
    selectivity->calculateBenefitAndGrade();
    selectivity->storeInDb(db, iech, iptr0, TEST, TEST);
  }

  if (iptr_sV >= 0) db->deleteColumnByUID(iptr_sV);
  if (iptr_yV >= 0) db->deleteColumnByUID(iptr_yV);
  return 0;
}

/*****************************************************************************/
/*!
 **  Correct the estimation and st. deviation of estimation if KO
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  iech         Rank of the sample
 ** \param[in]  col_est      Rank of the Kriging estimate
 ** \param[in]  col_std      Rank of the St, Deviation of Kriging estimate
 ** \param[in]  flag_OK      true when used in Ordinary Kriging
 **
 ** \param[out] krigest      Kriging estimation
 ** \param[out] krigstd      Standard deviation of the estimation error
 **
 ** \remarks In SK: krigstd returns the standard deviation of the estimation
 *error
 ** \remarks In OK: krigstd reads the square root of the estimation variance
 ** \remarks and returns the standard deviation
 **
 *****************************************************************************/
void CalcAnamTransform::_correctForOK(Db* db,
                                      int iech,
                                      int col_est,
                                      int col_std,
                                      bool flag_OK,
                                      double* krigest,
                                      double* krigstd)
{
  DECLARE_UNUSED(flag_OK);
  *krigest = db->getArray(iech, col_est);
  *krigstd = db->getArray(iech, col_std);

  //  if (flag_OK)
  //  {
  //    double var2 = 1. - (*krigstd) * (*krigstd);
  //    double stdv = sqrt(var2);
  //    *krigstd  = stdv;
  //  }
}

/*****************************************************************************/
/*!
 **  Prepare the vectors of estimation and st. dev.
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  col_est      Rank of the Kriging estimate
 ** \param[in]  col_std      Rank of the St, Deviation of Kriging estimate
 ** \param[in]  flag_OK      1 if kriging is performed with OK
 **
 ** \param[out] krigest      Vector of estamations
 ** \param[out] krigstd      Vector of standard deviation
 **
 *****************************************************************************/
void CalcAnamTransform::_getVectorsForCE(Db* db,
                                          int col_est,
                                          int col_std,
                                          bool flag_OK,
                                          VectorDouble& krigest,
                                          VectorDouble& krigstd)
{
  int nech = db->getNSample();
  krigest.resize(nech);
  krigstd.resize(nech);
  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    _correctForOK(db, iech, col_est, col_std, flag_OK, &krigest[iech],
                  &krigstd[iech]);
  }
}

/*****************************************************************************/
/*!
 **  Calculate the Conditional value and variance in the Gaussian Model
 **
 ** \return Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         Anamorphosis Hermite
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  iptr0        Starting storage address
 ** \param[in]  col_est      Rank of the Kriging estimate
 ** \param[in]  col_std      Rank of the St, Deviation of Kriging estimate
 ** \param[in]  nbsimu       Number of Monte Carlo simulations (0: Hermite)
 ** \param[in]  flag_OK      1 if kriging is performed with OK
 **
 *****************************************************************************/
int CalcAnamTransform::_ceZ(Db* db,
                            const AnamHermite* anam,
                            const Selectivity* selectivity,
                            int iptr0,
                            int col_est,
                            int col_std,
                            int nbsimu,
                            bool flag_OK)
{
  VectorDouble krigest, krigstd, valest, valstd;

  _getVectorsForCE(db, col_est, col_std, flag_OK, krigest, krigstd);

  if (nbsimu <= 0)
  {
    valest = hermiteCondExp(krigest, krigstd, anam->getPsiHns());
    valstd = hermiteCondStd(krigest, krigstd, anam->getPsiHns());
  }
  else
  {
    valest = MCCondExp(krigest, krigstd, anam->getPsiHns(), nbsimu);
    valstd = MCCondStd(krigest, krigstd, anam->getPsiHns(), nbsimu);
  }

  int iptrEst = selectivity->getAddressQTEst(ESelectivity::Z, iptr0);
  int iptrStd = selectivity->getAddressQTStd(ESelectivity::Z, iptr0);
  for (int iech = 0; iech < db->getNSample(); iech++)
  {
    if (!db->isActive(iech)) continue;
    if (selectivity->isUsedEst(ESelectivity::Z))
      db->setArray(iech, iptrEst, valest[iech]);
    if (selectivity->isUsedStD(ESelectivity::Z))
      db->setArray(iech, iptrStd, valstd[iech]);
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Calculate the Tonnage by Conditional Expectation
 **
 ** \return Error return code
 **
 ** \param[in]  mode         1 for T (Proba. above); 2 for Proba. below
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  iptr0        Starting storage address
 ** \param[in]  col_est      Rank of the Kriging estimate
 ** \param[in]  col_std      Rank of the St, Deviation of Kriging estimate
 ** \param[in]  ycuts        Vector of Gaussian cutoffs
 ** \param[in]  nbsimu       Number of Monte Carlo simulations (0: Hermite)
 ** \param[in]  flag_OK      1 if kriging is performed with OK
 **
 *****************************************************************************/
int CalcAnamTransform::_ceT(int mode,
                            Db* db,
                            const Selectivity* selectivity,
                            int iptr0,
                            int col_est,
                            int col_std,
                            const VectorDouble& ycuts,
                            int nbsimu,
                            bool flag_OK)
{
  VectorDouble krigest, krigstd, valest, valstd;

  /* Loop on the samples */

  for (int icut = 0; icut < selectivity->getNCuts(); icut++)
  {
    _getVectorsForCE(db, col_est, col_std, flag_OK, krigest, krigstd);

    if (nbsimu <= 0)
    {
      valest = hermiteIndicator(ycuts[icut], krigest, krigstd);
      valstd = hermiteIndicatorStd(ycuts[icut], krigest, krigstd);
    }
    else
    {
      valest = MCIndicator(ycuts[icut], krigest, krigstd, nbsimu);
      valstd = MCIndicatorStd(ycuts[icut], krigest, krigstd, nbsimu);
    }

    int iptrEst = selectivity->getAddressQTEst(ESelectivity::T, iptr0, icut);
    int iptrStd = selectivity->getAddressQTStd(ESelectivity::T, iptr0, icut);
    int iptrProba =
      selectivity->getAddressQTEst(ESelectivity::PROP, iptr0, icut);
    for (int iech = 0; iech < db->getNSample(); iech++)
    {
      if (!db->isActive(iech)) continue;
      if (mode == 1)
      {
        if (selectivity->isUsedEst(ESelectivity::T))
          db->setArray(iech, iptrEst, valest[iech]);
        if (selectivity->isUsedStD(ESelectivity::T))
          db->setArray(iech, iptrStd, valstd[iech]);
      }
      else
      {
        if (selectivity->isUsedEst(ESelectivity::PROP))
          db->setArray(iech, iptrProba, 1. - valest[iech]);
      }
    }
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Calculate the Quantile by Conditional Expectation
 **
 ** \return Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         Hermite anamorphosis
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  iptr0        Staring storage address
 ** \param[in]  col_est      Rank of the Kriging estimate
 ** \param[in]  col_std      Rank of the Kriging St, Deviation
 ** \param[in]  proba        Probability threshold
 ** \param[in]  flag_OK      1 if kriging is performed with OK
 **
 *****************************************************************************/
int CalcAnamTransform::_ceQuant(Db* db,
                                const AnamHermite* anam,
                                const Selectivity* selectivity,
                                int iptr0,
                                int col_est,
                                int col_std,
                                double proba,
                                bool flag_OK)
{
  double krigest, krigstd;

  int iptrQuant = selectivity->getAddressQTEst(ESelectivity::QUANT, iptr0);

  for (int iech = 0; iech < db->getNSample(); iech++)
  {
    if (!db->isActive(iech)) continue;
    _correctForOK(db, iech, col_est, col_std, flag_OK, &krigest, &krigstd);
    double y = krigest + krigstd * law_invcdf_gaussian(proba);
    db->setArray(iech, iptrQuant, anam->transformToRawValue(y));
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Calculate the Metal Quantity by Conditional Expectation
 **
 ** \return Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         Hermite anamorphosis
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  iptr0        Starting storage address
 ** \param[in]  col_est      Rank of the Kriging estimate
 ** \param[in]  col_std      Rank of the St, Deviation of Kriging estimate
 ** \param[in]  ycuts        Array of the requested cutoffs (gaussian scale)
 ** \param[in]  nbsimu       Number of Monte Carlo simulations (0: Hermite)
 ** \param[in]  flag_OK      1 if kriging is performed with OK
 **
 *****************************************************************************/
int CalcAnamTransform::_ceQ(Db* db,
                            const AnamHermite* anam,
                            const Selectivity* selectivity,
                            int iptr0,
                            int col_est,
                            int col_std,
                            const VectorDouble& ycuts,
                            int nbsimu,
                            bool flag_OK)
{
  VectorDouble krigest, krigstd, valest, valstd;

  /* Loop on the cutoffs */

  for (int icut = 0; icut < selectivity->getNCuts(); icut++)
  {
    _getVectorsForCE(db, col_est, col_std, flag_OK, krigest, krigstd);

    if (nbsimu <= 0)
    {
      valest = hermiteMetal(ycuts[icut], krigest, krigstd, anam->getPsiHns());
      valstd =
        hermiteMetalStd(ycuts[icut], krigest, krigstd, anam->getPsiHns());
    }
    else
    {
      valest =
        MCMetal(ycuts[icut], krigest, krigstd, anam->getPsiHns(), nbsimu);
      valstd =
        MCMetalStd(ycuts[icut], krigest, krigstd, anam->getPsiHns(), nbsimu);
    }

    int iptrEst = selectivity->getAddressQTEst(ESelectivity::Q, iptr0, icut);
    int iptrStd = selectivity->getAddressQTStd(ESelectivity::Q, iptr0, icut);
    for (int iech = 0; iech < db->getNSample(); iech++)
    {
      if (!db->isActive(iech)) continue;
      if (selectivity->isUsedEst(ESelectivity::Q))
        db->setArray(iech, iptrEst, valest[iech]);
      if (selectivity->isUsedStD(ESelectivity::Q))
        db->setArray(iech, iptrStd, valstd[iech]);
    }
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Calculate the Conventional Benefit by Conditional Expectation
 **
 ** \return Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  iptr0        Starting storage address
 ** \param[in]  ycuts        Array of the requested cutoffs
 **
 *****************************************************************************/
int CalcAnamTransform::_ceB(Db* db,
                            const Selectivity* selectivity,
                            int iptr0,
                            const VectorDouble& ycuts)
{
  double t, q;

  /* Loop on the cutoffs */

  for (int icut = 0; icut < selectivity->getNCuts(); icut++)
  {
    int iptrT = selectivity->getAddressQTEst(ESelectivity::T, iptr0, icut);
    int iptrQ = selectivity->getAddressQTEst(ESelectivity::Q, iptr0, icut);
    int iptrB = selectivity->getAddressQTEst(ESelectivity::B, iptr0, icut);

    /* Loop on the samples */

    for (int iech = 0; iech < db->getNSample(); iech++)
    {
      if (!db->isActive(iech)) continue;

      t = db->getArray(iech, iptrT);
      q = db->getArray(iech, iptrQ);
      db->setArray(iech, iptrB, q - t * ycuts[icut]);
    }
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Calculate the Average Grade by Conditional Expectation
 **
 ** \return Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  iptr0        Starting storage address
 **
 *****************************************************************************/
int CalcAnamTransform::_ceM(Db* db, const Selectivity* selectivity, int iptr0)
{
  double t, q;

  /* Loop on the cutoffs */

  for (int icut = 0; icut < selectivity->getNCuts(); icut++)
  {
    int iptrT = selectivity->getAddressQTEst(ESelectivity::T, iptr0, icut);
    int iptrQ = selectivity->getAddressQTEst(ESelectivity::Q, iptr0, icut);
    int iptrM = selectivity->getAddressQTEst(ESelectivity::M, iptr0, icut);

    /* Loop on the samples */

    for (int iech = 0; iech < db->getNSample(); iech++)
    {
      if (!db->isActive(iech)) continue;

      t = db->getArray(iech, iptrT);
      q = db->getArray(iech, iptrQ);
      db->setArray(iech, iptrM, (t > EPSILON3) ? q / t : TEST);
    }
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Transform a point anamorphosis into a block anamorphosis
 **
 ** \return  Error return code
 **
 ** \param[in]  anam        Point anamorphosis -> Block anamorphosis [out]
 ** \param[in]  verbose     Verbose option
 ** \param[in]  cvv         Block variance
 ** \param[in]  coeff       Coefficient of change of support
 ** \param[in]  mu          Additional coefficient for Discrete case
 **
 ** \remark If 'coeff' is provided, it is used directly ('cvv' is ignored)
 ** \remark Otherwise, it is derived from 'cvv'
 **
 *****************************************************************************/
int anamPointToBlock(AAnam* anam, int verbose, double cvv, double coeff, double mu)
{
  if (anam == nullptr) return (1);
  double r_coef                    = 0.;
  AnamHermite* anam_hermite        = dynamic_cast<AnamHermite*>(anam);
  AnamDiscreteDD* anam_discrete_DD = dynamic_cast<AnamDiscreteDD*>(anam);
  AnamDiscreteIR* anam_discrete_IR = dynamic_cast<AnamDiscreteIR*>(anam);

  /* Preliminary check */

  if (!FFFF(coeff) && (coeff < 0 || coeff > 1.))
  {
    messerr("Change of support coefficient (%lf) must lie between 0 and 1.",
            coeff);
    return 1;
  }

  /* Dispatch according to the anamorphosis type */

  switch (anam->getType().toEnum())
  {
    case EAnam::E_HERMITIAN:

      if (FFFF(coeff))
      {
        r_coef = sqrt(anam->invertVariance(cvv));
        if (verbose)
        {
          mestitle(1, "Calculation of the Change of Support Coefficient");
          message("Average Block covariance      = %lf\n", cvv);
          message("Change of support coefficient = %lf\n",
                  anam_hermite->getRCoef());
        }
      }
      else
        r_coef = coeff;
      break;

    case EAnam::E_DISCRETE_DD:
      anam_discrete_DD->setMu(mu);
      if (FFFF(coeff))
      {
        r_coef = sqrt(anam->invertVariance(cvv));
        if (verbose)
        {
          mestitle(1, "Calculation of the Change of Support Coefficient");
          message("Point Variance                = %lf\n",
                  anam_discrete_DD->getVariance());
          message("Average Block covariance      = %lf\n", cvv);
          message("Coefficient mu                = %lf\n",
                  anam_discrete_DD->getMu());
          message("Change of support coefficient = %lf\n",
                  anam_discrete_DD->getSCoef());
        }
      }
      else
        r_coef = coeff;
      break;

    case EAnam::E_DISCRETE_IR:
      if (FFFF(coeff))
      {
        r_coef = sqrt(anam->invertVariance(cvv));
        if (verbose)
        {
          mestitle(1, "Calculation of the Change of Support Coefficient");
          message("Average Block covariance      = %lf\n", cvv);
          message("Change of support coefficient = %lf\n",
                  anam_discrete_IR->getRCoef());
        }
      }
      else
        r_coef = coeff;
      break;

    default:
      messerr("The change of support is not defined for this Anamorphosis");
      return 1;
  }

  /* Update the Point Anamorphosis into Block Anamorphosis */

  anam->updatePointToBlock(r_coef);

  return 0;
}
