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
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Basic/Utilities.hpp"
#include "Db/Db.hpp"

#include <cmath>

#define RESIDUALS(icut, iech) (residuals[iech * ncut + icut])

namespace gstlrn
{
AnamDiscreteIR::AnamDiscreteIR(double rcoef)
  : AnamDiscrete()
  , _sCoef(rcoef)
{
}

AnamDiscreteIR::AnamDiscreteIR(const AnamDiscreteIR& m)
  : AnamDiscrete(m)
  , _sCoef(m._sCoef)
{
}

AnamDiscreteIR& AnamDiscreteIR::operator=(const AnamDiscreteIR& m)
{
  if (this != &m)
  {
    AnamDiscrete::operator=(m);
    _sCoef = m._sCoef;
  }
  return *this;
}

AnamDiscreteIR::~AnamDiscreteIR()
{
}

AnamDiscreteIR* AnamDiscreteIR::createFromNF(const String& NFFilename, bool verbose)
{
  auto* anam = new AnamDiscreteIR();
  if (anam->_fileOpenAndDeserialize(NFFilename, verbose)) return anam;
  delete anam;
  return nullptr;
}

AnamDiscreteIR* AnamDiscreteIR::create(double rcoef)
{
  return new AnamDiscreteIR(rcoef);
}

void AnamDiscreteIR::reset(Id ncut,
                           double r_coef,
                           const VectorDouble& zcut,
                           const VectorDouble& stats)
{
  setNCut(ncut);
  setZCut(zcut);
  setRCoef(r_coef);
  setStats(stats);
  calculateMeanAndVariance();
}

String AnamDiscreteIR::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  if (getNCut() <= 0 && getNElem() <= 0) return sstr.str();

  sstr << "Indicator Residuals Anamorphosis" << std::endl;

  sstr << AnamDiscrete::toString(strfmt);

  if (!_isFitted()) return sstr.str();

  if (_sCoef != 1.)
  {
    sstr << "Change of Support = " << _sCoef << std::endl;
  }

  sstr << "In the following printout:" << std::endl;
  sstr << "[,1] : Tonnage          'T'" << std::endl;
  sstr << "[,2] : Metal Quantity   'Q'" << std::endl;
  sstr << "[,3] : Average in class 'z'" << std::endl;
  sstr << "[,4] : B coefficient    'b'" << std::endl;
  sstr << "[,5] : Residual Point   'R'" << std::endl;
  sstr << "[,5] : Residual Block   'Rv'" << std::endl;
  sstr << std::endl;

  return sstr.str();
}

void AnamDiscreteIR::calculateMeanAndVariance()
{
  double var, mean;

  /* Loop on the classes */

  auto nclass = getNClass();
  var = mean = 0.;
  for (Id iclass = 0; iclass < nclass; iclass++)
  {
    double zn = (iclass < nclass - 1) ? getIRStatR(iclass + 1) : 0.;
    double tn = (iclass < nclass - 1) ? getIRStatT(iclass + 1) : 0.;
    var += getIRStatB(iclass) * getIRStatB(iclass) * zn;
    mean += getIRStatZ(iclass) * (getIRStatT(iclass) - tn);
  }
  setMean(mean);
  setVariance(var);
}

Id AnamDiscreteIR::fitFromArray(const VectorDouble& tab,
                                const VectorDouble& /*wt*/)
{
  double mean, dt, dq, tnext, qnext, tcur, tprev;
  Id nsorted;

  /* Initializations */

  Id nech     = static_cast<Id>(tab.size());
  auto nclass = getNClass();
  auto ncut   = getNCut();

  /* Core allocation */

  VectorDouble residuals(nech * ncut);
  VectorDouble T(ncut);
  VectorDouble Q(ncut);

  /* Calculate the residuals */

  if (_stats_residuals(false, nech, tab, &nsorted, &mean,
                       residuals.data(), T.data(), Q.data()))
    return 1;

  /* Store the the statistics */

  setMean(mean);
  setIRStatT(0, 1.);
  setIRStatQ(0, mean);
  for (Id icut = 0; icut < ncut; icut++)
  {
    setIRStatT(icut + 1, T[icut]);
    setIRStatQ(icut + 1, Q[icut]);
  }

  /* Calculate the B coefficients */

  for (Id iclass = 0; iclass < nclass; iclass++)
  {
    tnext = (iclass < nclass - 1) ? getIRStatT(iclass + 1) : 0.;
    qnext = (iclass < nclass - 1) ? getIRStatQ(iclass + 1) : 0.;
    tcur  = getIRStatT(iclass);
    dt    = tcur - tnext;
    dq    = getIRStatQ(iclass) - qnext;
    setIRStatZ(iclass, (dt > 0) ? dq / dt : 0.);
    setIRStatB(iclass, getIRStatQ(iclass) - tcur * getIRStatZ(iclass));
    if (iclass <= 0)
      setIRStatR(iclass, 0.);
    else
    {
      tprev = getIRStatT(iclass - 1);
      setIRStatR(iclass, (tprev > 0 && tcur > 0) ? 1. / tcur - 1. / tprev : 0.);
    }
    setIRStatRV(iclass, getIRStatR(iclass));
  }

  /* Evaluate mean and variance */

  calculateMeanAndVariance();

  return 0;
}

Id AnamDiscreteIR::_stats_residuals(Id verbose,
                                    Id nech,
                                    const VectorDouble& tab,
                                    Id* nsorted,
                                    double* mean,
                                    double* residuals,
                                    double* T,
                                    double* Q)
{
  double value, moyenne;
  Id iech, icut, jcut, nactive;

  /* Initializations */

  auto ncut = getNCut();
  nactive = (*nsorted) = 0;
  moyenne              = 0.;
  for (icut = 0; icut < ncut; icut++)
  {
    T[icut] = Q[icut] = 0.;
    for (iech = 0; iech < nech; iech++)
      RESIDUALS(icut, iech) = 0.;
  }

  /* Loop on the samples to calculate the indicators */

  for (iech = 0; iech < nech; iech++)
  {
    value = tab[iech];
    if (FFFF(value)) continue;
    moyenne += value;
    nactive++;

    /* Loop on the cutoffs */

    for (icut = 0; icut < ncut; icut++)
    {
      if (value < getZCut(icut)) continue;
      RESIDUALS(icut, iech) = 1.;
      Q[icut] += value;
      T[icut] += 1.;
    }
  }
  if (nactive <= 0)
  {
    messerr("The calculation failed as there is no active sample");
    return (1);
  }

  /* Calculate the tonnage and meal quantity per class */

  moyenne /= static_cast<double>(nactive);
  for (icut = 0; icut < ncut; icut++)
  {
    T[icut] /= static_cast<double>(nactive);
    Q[icut] /= static_cast<double>(nactive);
  }

  /* Calculate the residuals */

  for (iech = 0; iech < nech; iech++)
  {
    value = tab[iech];
    if (FFFF(value)) continue;

    /* Loop on the cutoffs */

    for (icut = ncut - 1; icut >= 0; icut--)
    {
      value = RESIDUALS(icut, iech) / T[icut];
      if (icut > 0)
      {
        jcut = icut - 1;
        value -= RESIDUALS(jcut, iech) / T[jcut];
      }
      else
      {
        value -= 1.;
      }
      RESIDUALS(icut, iech) = value;
    }
  }

  /* Verbose optional option */

  if (verbose)
  {
    mestitle(0, "Building residuals");
    message("Number of sorted samples = %d\n", nactive);
    for (icut = 0; icut < ncut; icut++)
      message("Cutoff %2d (above %lf) - Tonnage = %lf - Metal = %lf\n",
              icut + 1, getZCut(icut), T[icut], Q[icut]);
  }

  (*nsorted) = nactive;
  (*mean)    = moyenne;
  return (0);
}

VectorDouble AnamDiscreteIR::z2factor(double z, const VectorInt& ifacs) const
{
  VectorDouble factors;
  Id nfact = static_cast<Id>(ifacs.size());
  factors.resize(nfact, 0);

  for (Id ifac = 0; ifac < nfact; ifac++)
  {
    Id iclass     = ifacs[ifac] - 1;
    double v1     = _getResidual(iclass, z);
    double v2     = _getResidual(iclass - 1, z);
    factors[ifac] = v1 - v2;
  }
  return factors;
}

/**
 *
 * @param iclass Rank of the class
 * @param z Input value
 * @return Calculate the normalized residual for a given value
 */
double AnamDiscreteIR::_getResidual(Id iclass, double z) const
{
  double seuil  = (iclass < 0) ? 0. : getZCut(iclass);
  double retval = (z >= seuil) ? 1. : 0.;
  retval /= getIRStatT(iclass + 1);

  return (retval);
}

bool AnamDiscreteIR::_serializeAscii(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret      = ret && AnamDiscrete::_serializeAscii(os, verbose);
  ret      = ret && _recordWrite<double>(os, "Change of support coefficient", getRCoef());
  return ret;
}

bool AnamDiscreteIR::_deserializeAscii(std::istream& is, bool verbose)
{
  double r = TEST;

  bool ret = true;
  ret      = ret && AnamDiscrete::_deserializeAscii(is, verbose);
  ret      = ret && _recordRead<double>(is, "Anamorphosis 'r' coefficient", r);
  if (ret) setRCoef(r);
  return ret;
}

double AnamDiscreteIR::computeVariance(double sval) const
{
  if (!allowChangeSupport()) return TEST;
  auto nclass = getNClass();

  double var = 0.;
  for (Id iclass = 0; iclass < nclass - 1; iclass++)
  {
    double b     = getIRStatB(iclass);
    double tcur  = getIRStatT(iclass + 1);
    double tprev = getIRStatT(iclass);
    double resid = (tprev > 0 && tcur > 0) ? 1. / pow(tcur, sval) - 1. / pow(tprev, sval) : 0.;
    var += b * b * resid;
  }
  return (var);
}

Id AnamDiscreteIR::updatePointToBlock(double r_coef)
{
  if (!allowChangeSupport()) return 1;
  setRCoef(r_coef);

  auto nclass = getNClass();
  double zie  = 0.;
  double zik  = 0.;
  double zje  = 0.;
  double zjk  = 0.;

  /* Loop on the classes */

  for (Id iclass = 0; iclass < nclass; iclass++)
  {
    double tcur = getIRStatT(iclass);
    if (iclass > 0)
    {
      zjk = getIRStatZ(iclass - 1);
      zie = getIRStatZ(iclass);
      zik = (tcur > 0) ? zjk + (zie - zje) * pow(tcur, 1. - r_coef) : 0.;
    }
    else
    {
      zik = getIRStatZ(iclass);
    }
    zje = getIRStatZ(iclass);

    double Tval = getIRStatT(iclass);
    double Bval = getIRStatB(iclass);
    setIRStatZ(iclass, zik);
    setIRStatT(iclass, pow(Tval, r_coef));
    setIRStatQ(iclass, Bval + zik * Tval);
    if (iclass <= 0)
      setIRStatRV(iclass, 0.);
    else
    {
      tcur         = getIRStatT(iclass);
      double tprev = getIRStatT(iclass - 1);
      setIRStatRV(
        iclass, (tprev > 0 && tcur > 0) ? 1. / tcur - 1 / tprev : 0.);
    }
  }

  /* Update mean and variance */

  calculateMeanAndVariance();

  return 0;
}

/****************************************************************************/
/*!
 **  Calculate the theoretical grade tonnage value (Discrete Indicator Residuals)
 **
 *****************************************************************************/
void AnamDiscreteIR::_globalSelectivity(Selectivity* selectivity)
{
  bool cutDefined = (selectivity->getNCuts() > 0);

  /* Define the working Selectivity structure */

  Selectivity* selloc;
  if (cutDefined)
    selloc = selectivity;
  else
  {
    selloc = selectivity->clone();
    selloc->resetCuts(getZCut());
  }
  Id ncuts = selloc->getNCuts();

  /* Calculate the Grade-Tonnage curves */

  for (Id icut = 0; icut < ncuts; icut++)
  {
    selloc->setZcut(icut, (icut == 0) ? 0. : getZCut(icut - 1));
    selloc->setTest(icut, getIRStatT(icut));
    selloc->setQest(icut, getIRStatQ(icut));
  }

  /* Correct order relationship */

  selloc->correctTonnageOrder();

  /* Store the results */

  if (cutDefined)
  {
    selectivity->interpolateSelectivity(selloc);
    selectivity->calculateBenefitAndGrade();
  }
  else
  {
    selloc->calculateBenefitAndGrade();
  }
}

/*****************************************************************************/
/*!
 **  Calculate Experimental Grade-Tonnage curves from factors
 **  Case of Discrete Indicator Residuals
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  cols_est     Array of UIDs for factor estimation
 ** \param[in]  cols_std     Array of UIDs for factor St. Dev.
 ** \param[in]  iptr0        Rank for storing the results
 **
 *****************************************************************************/
Id AnamDiscreteIR::factor2Selectivity(Db* db,
                                      Selectivity* selectivity,
                                      const VectorInt& cols_est,
                                      const VectorInt& cols_std,
                                      Id iptr0)
{
  Id nech         = db->getNSample();
  Id nb_est       = static_cast<Id>(cols_est.size());
  Id nb_std       = static_cast<Id>(cols_std.size());
  Id ncleff       = MAX(nb_est, nb_std);
  bool cutDefined = (selectivity->getNCuts() > 0);

  if (db == nullptr)
  {
    messerr("You must define a Db");
    return 1;
  }
  Id nfactor = MAX(nb_est, nb_std);
  if (nfactor > getNClass())
  {
    messerr("Number of factors (%d) must be smaller or equal to Number of Classes",
            nfactor, getNClass());
    return 1;
  }

  /* Define the working Selectivity structure */

  Selectivity* selloc;
  if (cutDefined)
    selloc = selectivity;
  else
  {
    selloc = selectivity->clone();
    selloc->resetCuts(getZCut());
  }

  /* Calculate the Recovery Functions from the factors */

  for (Id iech = 0; iech < nech; iech++)
  {
    if (_isSampleSkipped(db, iech, cols_est, cols_std)) continue;

    /* Calculate the tonnage and the recovered grade */

    double total = 0.;
    for (Id ivar = 0; ivar < ncleff; ivar++)
    {
      double value = db->getArray(iech, cols_est[ivar]);
      total += value;
      selloc->setTest(ivar, total * getIRStatT(ivar + 1));
    }

    /* Correct order relationship */

    selloc->correctTonnageOrder();

    /* Tonnage: Standard Deviation */

    if (selloc->isUsedStD(ESelectivity::T))
    {
      total = 0.;
      for (Id ivar = 0; ivar < ncleff; ivar++)
      {
        double value = db->getArray(iech, cols_std[ivar]);
        total += value * value;
        selloc->setTstd(ivar, sqrt(total) * getIRStatT(ivar + 1));
      }
    }

    /* Metal Quantity: Estimation */

    if (selloc->isUsedEst(ESelectivity::Q))
    {
      selloc->setQest(ncleff - 1, getIRStatZ(ncleff) * selloc->getTest(ncleff - 1));
      for (Id ivar = ncleff - 2; ivar >= 0; ivar--)
      {
        selloc->setQest(ivar,
                        selloc->getQest(ivar + 1) + getIRStatZ(ivar + 1) * (selloc->getTest(ivar) - selloc->getTest(ivar + 1)));
      }
    }

    /* Metal Quantity: Standard Deviation */

    if (selloc->isUsedStD(ESelectivity::Q))
    {
      for (Id ivar = 0; ivar < ncleff; ivar++)
      {
        total = 0.;
        for (Id jvar = 0; jvar < ivar; jvar++)
        {
          double value = db->getArray(iech, cols_std[jvar]);
          total += value * value;
        }
        double prod = getIRStatB(ivar + 1) +
                      getIRStatZ(ivar + 1) * getIRStatT(ivar + 1);
        total *= prod * prod;
        for (Id jvar = ivar + 1; jvar < ncleff; jvar++)
        {
          prod = db->getArray(iech, cols_std[jvar]) * getIRStatB(ivar + 1);
          total += prod * prod;
        }
        selloc->setQstd(ivar, sqrt(total));
      }
    }

    /* Z: Estimation */

    double zestim = 0.;
    if (selloc->isUsedEst(ESelectivity::Z))
    {
      zestim = getIRStatZ(ncleff) * selloc->getTest(ncleff - 1);
      for (Id ivar = 0; ivar < ncleff - 1; ivar++)
        zestim += getIRStatZ(ivar + 1) * (selloc->getTest(ivar) - selloc->getTest(ivar + 1));
    }

    /* Z: Standard Deviation */

    double zstdev = 0.;
    if (selloc->isUsedStD(ESelectivity::Z))
    {
      total = 0.;
      for (Id ivar = 0; ivar < ncleff; ivar++)
      {
        double prod = db->getArray(iech, cols_std[ivar]) * getIRStatB(ivar + 1);
        total += prod * prod;
      }
      zstdev = sqrt(total);
    }

    /* Store the results */

    if (cutDefined)
    {
      selectivity->interpolateSelectivity(selloc);
      selectivity->calculateBenefitAndGrade();
      selectivity->storeInDb(db, iech, iptr0, zestim, zstdev);
    }
    else
    {
      selloc->calculateBenefitAndGrade();
      selloc->storeInDb(db, iech, iptr0, zestim, zstdev);
    }
  }
  return (0);
}

#ifdef HDF5
bool AnamDiscreteIR::_deserializeH5(H5::Group& grp, [[maybe_unused]] bool verbose)
{
  auto anamG = SerializeHDF5::getGroup(grp, "AnamDiscreteIR");
  if (!anamG)
  {
    return false;
  }

  /* Read the grid characteristics */
  bool ret = true;
  double r = 0.;

  ret = ret && SerializeHDF5::readValue(*anamG, "R", r);

  ret = ret && AnamDiscrete::_deserializeH5(*anamG, verbose);

  if (ret)
    setRCoef(r);

  return ret;
}

bool AnamDiscreteIR::_serializeH5(H5::Group& grp, [[maybe_unused]] bool verbose) const
{
  auto anamG = grp.createGroup("AnamDiscreteIR");

  bool ret = true;

  ret = ret && SerializeHDF5::writeValue(anamG, "R", getRCoef());

  ret = ret && AnamDiscrete::_serializeH5(anamG, verbose);

  return ret;
}
#endif
} // namespace gstlrn