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
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Db/Db.hpp"
#include "Basic/Utilities.hpp"

#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "geoslib_enum.h"

#include <math.h>

#define RESIDUALS(icut,iech) (residuals[iech * ncut + icut])
#define QT_EST    0
#define QT_STD    1
#define QT_VARS(i,j)              (qt_vars[(i) + 2 * (j)])
#define QT_FLAG(j)                (QT_VARS(QT_EST,j) > 0 || \
                                   QT_VARS(QT_STD,j) > 0)

AnamDiscreteIR::AnamDiscreteIR(double rcoef)
    : AnamDiscrete(),
      _sCoef(rcoef)
{
}

AnamDiscreteIR::AnamDiscreteIR(const AnamDiscreteIR &m)
    : AnamDiscrete(m),
      _sCoef(m._sCoef)
{

}

AnamDiscreteIR& AnamDiscreteIR::operator=(const AnamDiscreteIR &m)
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

AnamDiscreteIR* AnamDiscreteIR::createFromNF(const String& neutralFilename, bool verbose)
{
  AnamDiscreteIR* anam = nullptr;
  std::ifstream is;
  anam = new AnamDiscreteIR();
  bool success = false;
  if (anam->_fileOpenRead(neutralFilename, is, verbose))
  {
    success =  anam->deserialize(is, verbose);
  }
  if (! success)
  {
    delete anam;
    anam = nullptr;
  }
  return anam;
}

AnamDiscreteIR* AnamDiscreteIR::create(double rcoef)
{
  return new AnamDiscreteIR(rcoef);
}

void AnamDiscreteIR::reset(int ncut,
                           double r_coef,
                           const VectorDouble &zcut,
                           const VectorDouble &stats)
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

  int nclass = getNClass();
  var = mean = 0.;
  for (int iclass = 0; iclass < nclass; iclass++)
  {
    double zn = (iclass < nclass - 1) ? getIRStatR(iclass + 1) :
                                        0.;
    double tn = (iclass < nclass - 1) ? getIRStatT(iclass + 1) :
                                        0.;
    var += getIRStatB(iclass) * getIRStatB(iclass) * zn;
    mean += getIRStatZ(iclass) * (getIRStatT(iclass) - tn);
  }
  setMean(mean);
  setVariance(var);
}

int AnamDiscreteIR::fit(const VectorDouble& tab, int verbose)
{
  double *residuals, *T, *Q, mean, dt, dq, tnext, qnext, tcur, tprev;
  int error, nsorted;

  /* Initializations */

  error = 1;
  int nech = static_cast<int> (tab.size());
  int nclass = getNClass();
  int ncut = getNCut();
  residuals = T = Q = nullptr;

  /* Core allocation */

  residuals = (double *) mem_alloc(sizeof(double) * nech * ncut, 1);
  T = (double *) mem_alloc(sizeof(double) * ncut, 1);
  Q = (double *) mem_alloc(sizeof(double) * ncut, 1);

  /* Calculate the residuals */

  if (_stats_residuals(verbose, nech, tab, &nsorted, &mean, residuals, T, Q))
    goto label_end;

  /* Store the the statistics */

  setMean(mean);
  setIRStatT(0, 1.);
  setIRStatQ(0, mean);
  for (int icut = 0; icut < ncut; icut++)
  {
    setIRStatT(icut + 1, T[icut]);
    setIRStatQ(icut + 1, Q[icut]);
  }

  /* Calculate the B coefficients */

  for (int iclass = 0; iclass < nclass; iclass++)
  {
    tnext = (iclass < nclass - 1) ? getIRStatT(iclass + 1) : 0.;
    qnext = (iclass < nclass - 1) ? getIRStatQ(iclass + 1) : 0.;
    tcur = getIRStatT(iclass);
    dt = tcur - tnext;
    dq = getIRStatQ(iclass) - qnext;
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

  /* Set the error return code */

  error = 0;

  label_end: residuals = (double *) mem_free((char * ) residuals);
  T = (double *) mem_free((char * ) T);
  Q = (double *) mem_free((char * ) Q);
  return (error);
}

int AnamDiscreteIR::_stats_residuals(int verbose,
                                       int nech,
                                       const VectorDouble& tab,
                                       int *nsorted,
                                       double *mean,
                                       double *residuals,
                                       double *T,
                                       double *Q)
{
  double value, moyenne;
  int iech, icut, jcut, nactive;

  /* Initializations */

  int ncut = getNCut();
  nactive = (*nsorted) = 0;
  moyenne = 0.;
  for (icut = 0; icut < ncut; icut++)
  {
    T[icut] = Q[icut] = 0.;
    for (iech = 0; iech < nech; iech++)
      RESIDUALS(icut,iech) = 0.;
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
      RESIDUALS(icut,iech) = 1.;
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

  moyenne /= (double) nactive;
  for (icut = 0; icut < ncut; icut++)
  {
    T[icut] /= (double) nactive;
    Q[icut] /= (double) nactive;
  }

  /* Calculate the residuals */

  for (iech = 0; iech < nech; iech++)
  {
    value = tab[iech];
    if (FFFF(value)) continue;

    /* Loop on the cutoffs */

    for (icut = ncut - 1; icut >= 0; icut--)
    {
      value = RESIDUALS(icut,iech) / T[icut];
      if (icut > 0)
      {
        jcut = icut - 1;
        value -= RESIDUALS(jcut,iech) / T[jcut];
      }
      else
      {
        value -= 1.;
      }
      RESIDUALS(icut,iech) = value;
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
  (*mean) = moyenne;
  return (0);
}

VectorDouble AnamDiscreteIR::z2factor(double z, const VectorInt& ifacs) const
{
  VectorDouble factors;
  int nfact = (int) ifacs.size();
  factors.resize(nfact, 0);

  for (int ifac = 0; ifac < nfact; ifac++)
  {
    int iclass = ifacs[ifac] - 1;
    double v1 = _getResidual(iclass, z);
    double v2 = _getResidual(iclass - 1, z);
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
double AnamDiscreteIR::_getResidual(int iclass, double z) const
{
  double seuil = (iclass < 0) ? 0. : getZCut(iclass);
  double retval = (z >= seuil) ? 1. : 0.;
  retval /= getIRStatT(iclass + 1);

  return (retval);
}

bool AnamDiscreteIR::_serialize(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret = ret && AnamDiscrete::_serialize(os, verbose);
  ret = ret && _recordWrite<double>(os, "Change of support coefficient", getRCoef());
  return ret;
}

bool AnamDiscreteIR::_deserialize(std::istream& is, bool verbose)
{
  double r = TEST;

  bool ret = true;
  ret = ret && AnamDiscrete::_deserialize(is, verbose);
  ret = ret && _recordRead<double>(is, "Anamorphosis 'r' coefficient", r);
  if (ret) setRCoef(r);
  return ret;
}

double AnamDiscreteIR::getBlockVariance(double sval, double /*power*/) const
{
  if (! allowChangeSupport()) return TEST;
  int nclass = getNClass();

  double var = 0.;
  for (int iclass = 0; iclass < nclass - 1; iclass++)
  {
    double b = getIRStatB(iclass);
    double tcur = getIRStatT(iclass + 1);
    double tprev = getIRStatT(iclass);
    double resid = (tprev > 0 && tcur > 0) ?
        1. / pow(tcur, sval) - 1. / pow(tprev, sval) : 0.;
    var += b * b * resid;
  }
  return (var);
}

int AnamDiscreteIR::updatePointToBlock(double r_coef)
{
  if (! allowChangeSupport()) return 1;
  setRCoef(r_coef);

  int nclass = getNClass();
  double zie = 0.;
  double zik = 0.;
  double zje = 0.;
  double zjk = 0.;

  /* Loop on the classes */

  for (int iclass = 0; iclass < nclass; iclass++)
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
      tcur = getIRStatT(iclass);
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
 ** \param[in] flag_correct 1 if Tonnage order relationship must be corrected
 **
 *****************************************************************************/
Selectivity AnamDiscreteIR::calculateSelectivity(bool flag_correct)
{
  int nclass = getNClass();
  Selectivity calest(nclass);

  /* Calculate the Grade-Tonnage curves */

  for (int iclass = 0; iclass < nclass; iclass++)
  {
    calest.setZcut(iclass, (iclass == 0) ? 0. : getZCut(iclass - 1));
    calest.setTest(iclass, getIRStatT(iclass));
    calest.setQest(iclass, getIRStatQ(iclass));
  }

  /* Correct order relationship */

  if (flag_correct) calest.correctTonnageOrder();

  /* Store the results */

  calest.calculateBenefitGrade();

  return calest;
}

/*****************************************************************************/
/*!
 **  Calculate Experimental Grade-Tonnage curves from factors
 **  Case of Discrete Indicator Residuals
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  cutmine      Array of the requested cutoffs
 ** \param[in]  z_max        Maximum grade array (only for QT interpolation)
 ** \param[in]  flag_correct 1 if Tonnage order relationship must be corrected
 ** \param[in]  cols_est     Array of columns for factor estimation
 ** \param[in]  cols_std     Array of columns for factor st. dev.
 ** \param[in]  iptr         Rank for storing the results
 ** \param[in]  codes        Array of codes for stored results
 ** \param[in]  qt_vars      Array of variables to be calculated
 **
 *****************************************************************************/
int AnamDiscreteIR::factor2QT(Db *db,
                              const VectorDouble& cutmine,
                              double z_max,
                              int flag_correct,
                              const VectorInt& cols_est,
                              const VectorInt& cols_std,
                              int iptr,
                              const VectorInt& codes,
                              VectorInt& qt_vars)
{
  int nech = db->getSampleNumber();
  int nb_est = (int) cols_est.size();
  int nb_std = (int) cols_std.size();
  int ncleff = MAX(nb_est, nb_std);
  int ncutmine = (int) cutmine.size();

  if (db == nullptr)
  {
    messerr("You must define a Db");
    return 1;
  }
  if (nb_est <= 0 && nb_std <= 0)
  {
    messerr("The number of factors is zero");
    return 1;
  }
  int nvar = MAX(nb_est, nb_std);
  int nmax = getNCut();
  if (nvar > nmax)
  {
    messerr("Number of factors (%d) must be smaller or equal to Number of cutoffs",
            nvar, nmax);
    return 1;
  }
  if (ncutmine <= 0) ncutmine = nmax;

  /* Core allocation */

  Selectivity calest(nmax);
  Selectivity calcut;
  if (ncutmine > 0)
    calcut = Selectivity(ncutmine);

  /* Calculate the Recovery Functions from the factors */

  for (int iech = 0; iech < nech; iech++)
  {
    if (_isSampleSkipped(db, iech, cols_est, cols_std)) continue;

    /* Calculate the tonnage and the recovered grade */

    double total = 0.;
    for (int ivar = 0; ivar < ncleff; ivar++)
    {
      double value = db->getArray(iech, cols_est[ivar]);
      total += value;
      calest.setTest(ivar, total * getIRStatT(ivar + 1));
    }

    /* Correct order relationship */

    if (flag_correct) calest.correctTonnageOrder();

    /* Tonnage: Standard Deviation */

    if (QT_VARS(QT_STD,ANAM_QT_T) > 0)
    {
      total = 0.;
      for (int ivar = 0; ivar < ncleff; ivar++)
      {
        double value = db->getArray(iech, cols_std[ivar]);
        total += value * value;
        calest.setTstd(ivar, sqrt(total) * getIRStatT(ivar + 1));
      }
    }

    /* Metal Quantity: Estimation */

    if (QT_VARS(QT_EST,ANAM_QT_Q) > 0)
    {
      calest.setQest(ncleff-1, getIRStatZ(ncleff) * calest.getTest(ncleff-1));
      for (int ivar = ncleff - 2; ivar >= 0; ivar--)
      {
        calest.setQest(ivar,
                       calest.getQest(ivar+1) + getIRStatZ(ivar + 1)
                       * (calest.getTest(ivar) - calest.getTest(ivar + 1)));
      }
    }

    /* Metal Quantity: Standard Deviation */

    if (QT_VARS(QT_STD,ANAM_QT_Q) > 0)
    {
      for (int ivar = 0; ivar < ncleff; ivar++)
      {
        total = 0.;
        for (int jvar = 0; jvar < ivar; jvar++)
        {
          double value = db->getArray(iech, cols_std[jvar]);
          total += value * value;
        }
        double prod = getIRStatB(ivar + 1) +
            getIRStatZ(ivar + 1) * getIRStatT(ivar + 1);
        total *= prod * prod;
        for (int jvar = ivar + 1; jvar < ncleff; jvar++)
        {
          prod = db->getArray(iech, cols_std[jvar]) * getIRStatB(ivar + 1);
          total += prod * prod;
        }
        calest.setQstd(ivar, sqrt(total));
      }
    }

    /* Z: Estimation */

    double zestim = 0.;
    if (QT_VARS(QT_EST,ANAM_QT_Z) > 0)
    {
      zestim = getIRStatZ(ncleff) * calest.getTest(ncleff-1);
      for (int ivar = 0; ivar < ncleff - 1; ivar++)
        zestim += getIRStatZ(ivar + 1)
            * (calest.getTest(ivar) - calest.getTest(ivar + 1));
    }

    /* Z: Standard Deviation */

    double zstdev = 0.;
    if (QT_VARS(QT_STD,ANAM_QT_Z) > 0)
    {
      total = 0.;
      for (int ivar = 0; ivar < ncleff; ivar++)
      {
        double prod = db->getArray(iech, cols_std[ivar])
            * getIRStatB(ivar + 1);
        total += prod * prod;
      }
      zstdev = sqrt(total);
    }

    /* Store the results */

    if (ncutmine > 0)
    {
      _interpolateQTLocal(z_max, cutmine, calest, calcut);
      calcut.calculateBenefitGrade();
      recoveryLocal(db, iech, iptr, codes, qt_vars, zestim, zstdev, calcut);
    }
    else
    {
      calest.calculateBenefitGrade();
      recoveryLocal(db, iech, iptr, codes, qt_vars, zestim, zstdev, calest);
    }
  }
  return (0);
}
