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

#include "Basic/Utilities.hpp"
#include "geoslib_f.h"

#define RESIDUALS(icut,iech) (residuals[iech * ncut + icut])

AnamDiscreteIR::AnamDiscreteIR()
    : AnamDiscrete(ANAM_DISCRETE_IR),
      _rCoef(0.)
{
}

AnamDiscreteIR::AnamDiscreteIR(const AnamDiscreteIR &m)
    : AnamDiscrete(m),
      _rCoef(m._rCoef)
{

}

AnamDiscreteIR& AnamDiscreteIR::operator=(const AnamDiscreteIR &m)
{
  if (this != &m)
  {
    _rCoef = m._rCoef;
  }
  return *this;
}

AnamDiscreteIR::~AnamDiscreteIR()
{

}

String AnamDiscreteIR::toString(int level) const
{
  std::stringstream sstr;
  sstr << Anam::toString(level);
  sstr << "Discrete Indicator Residuals Anamorphosis" << std::endl;
  if (_rCoef != 1.)
  {
    sstr << "Change of Support = " << _rCoef << std::endl;
  }

  sstr << std::endl;
  sstr << "In the following printout:" << std::endl;
  sstr << "[,1] : Tonnage          'T'" << std::endl;
  sstr << "[,2] : Metal Quantity   'Q'" << std::endl;
  sstr << "[,3] : Average in class 'z'" << std::endl;
  sstr << "[,4] : B coefficient    'b'" << std::endl;
  sstr << "[,5] : Residual Point   'R'" << std::endl;
  sstr << "[,5] : Residual Block   'Rv'" << std::endl;
  sstr << std::endl;
  sstr << toMatrix(String(), VectorString(), VectorString(), true, getNElem(),
                  getNClass(), getStats().getValues());

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
  residuals = T = Q = (double *) NULL;

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
    tnext = (iclass < nclass - 1) ? getIRStatT(iclass + 1) :
                                    0.;
    qnext = (iclass < nclass - 1) ? getIRStatQ(iclass + 1) :
                                    0.;
    tcur = getIRStatT(iclass);
    dt = tcur - tnext;
    dq = getIRStatQ(iclass) - qnext;
    setIRStatZ(iclass, (dt > 0) ? dq / dt :
                                  0.);
    setIRStatB(iclass, getIRStatQ(iclass) - tcur * getIRStatZ(iclass));
    if (iclass <= 0)
      setIRStatR(iclass, 0.);
    else
    {
      tprev = getIRStatT(iclass - 1);
      setIRStatR(iclass, (tprev > 0 && tcur > 0) ? 1. / tcur - 1. / tprev :
                                                   0.);
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

VectorDouble AnamDiscreteIR::z2f(int nfact,
                                   const VectorInt& ifacs,
                                   double z) const
{
  VectorDouble factors;
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
  double seuil = (iclass < 0) ? 0. :
                                getZCut(iclass);
  double retval = (z >= seuil) ? 1. :
                                 0.;
  retval /= getIRStatT(iclass + 1);

  return (retval);
}

