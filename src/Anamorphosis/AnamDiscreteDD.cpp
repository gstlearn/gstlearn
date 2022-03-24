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
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "Anamorphosis/AnamDiscreteDD.hpp"
#include "Stats/Selectivity.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"

#include <math.h>

#define EIGVEC(i,j)      eigvec[(i)*nclass+(j)]
#define CHI(i,j)         chi[(i)*nclass+(j)]
#define I2CHI(i,j)       _i2Chi[(i)*nclass+(j)]
#define PTAB(i,j)        ptab[(i)*nclass+(j)]
#define C_S(i,j)         c_s[(i)*nclass+(j)]
#define Q_S(i,j)         q_s[(i)*nclass+(j)]
#define MATI(i,j)        mati[(i)*nclass+(j)]

AnamDiscreteDD::AnamDiscreteDD(double mu, double scoef)
    : AnamDiscrete(),
      _mu(mu),
      _rCoef(scoef),
      _maf(),
      _i2Chi()
{
}

AnamDiscreteDD::AnamDiscreteDD(const AnamDiscreteDD &m)
    : AnamDiscrete(m),
      _mu(m._mu),
      _rCoef(m._rCoef),
      _maf(m._maf),
      _i2Chi(m._i2Chi)
{
}

AnamDiscreteDD& AnamDiscreteDD::operator=(const AnamDiscreteDD &m)
{
  if (this != &m)
  {
    _mu = m._mu;
    _rCoef = m._rCoef;
    _maf = m._maf;
    _i2Chi = m._i2Chi;
  }
  return *this;
}

AnamDiscreteDD::~AnamDiscreteDD()
{

}

int AnamDiscreteDD::dumpToNF(const String& neutralFilename, bool verbose) const
{
  std::ofstream os;
  int ret = 1;
  if (_fileOpenWrite(neutralFilename, "AnamDiscreteDD", os, verbose))
  {
    ret = _serialize(os, verbose);
    if (ret && verbose) messerr("Problem writing in the Neutral File.");
    os.close();
  }
  return ret;
}

AnamDiscreteDD* AnamDiscreteDD::createFromNF(const String& neutralFilename, bool verbose)
{
  AnamDiscreteDD* anam = nullptr;
  std::ifstream is;
  if (_fileOpenRead(neutralFilename, "AnamDiscreteDD", is, verbose))
  {
    anam = new AnamDiscreteDD();
    if (anam->_deserialize(is, verbose))
    {
      if (verbose) messerr("Problem reading the Neutral File");
      delete anam;
      anam = nullptr;
    }
    is.close();
  }
  return anam;
}

AnamDiscreteDD* AnamDiscreteDD::create(double mu, double scoef)
{
  return new AnamDiscreteDD(mu, scoef);
}

void AnamDiscreteDD::reset(int ncut,
                           double scoef,
                           double mu,
                           const VectorDouble &zcut,
                           const VectorDouble &pcaz2f,
                           const VectorDouble &pcaf2z,
                           const VectorDouble &stats)
{
  setNCut(ncut);
  setZCut(zcut);
  setRCoef(scoef);
  setMu(mu);
  setPcaF2Z(pcaf2z);
  setPcaZ2F(pcaz2f);
  setStats(stats);
  calculateMeanAndVariance();
}

String AnamDiscreteDD::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  sstr << "Discrete Diffusion Anamorphosis" << std::endl;

  sstr << AnamDiscrete::toString(strfmt);

  if (_rCoef != 0.)
  {
    sstr << "Mu Coefficient    = " << _mu << std::endl;
    sstr << "Change of Support = " << _rCoef << std::endl;
  }

  sstr << std::endl;
  sstr << "In the previous printout:" << std::endl;
  sstr << "[,1] : Class Probability 'w'" << std::endl;
  sstr << "[,2] : Class Mean Value 'zc'" << std::endl;
  sstr << "[,3] : Anamorphosis coefficient 'c_s'" << std::endl;
  sstr << "[,4] : Spectral Value 'lambda'" << std::endl;
  sstr << "[,5] : Spectral Weight 'U'" << std::endl;
  sstr << "[,6] : Terms pow(mu/(mu+li),s/2)" << std::endl;
  sstr << std::endl;

  return sstr.str();
}

void AnamDiscreteDD::calculateMeanAndVariance()
{
  double var,mean;

  mean = var = 0.;
  for (int iclass=0; iclass<getNClass(); iclass++)
  {
    double prop  = getDDStatProp(iclass);
    double zval  = getDDStatZmoy(iclass);
    mean += zval * prop;
    var  += zval * zval * prop;
  }
  var -= mean * mean;
  setMean(mean);
  setVariance(var);
}

int AnamDiscreteDD::fit(const VectorDouble& tab, bool verbose)
{
  VectorDouble chi;

  int nech = static_cast<int> (tab.size());

  // Calculate statistics on data

  _stats(nech,tab);

  // Modeling the diffusion process

  chi = factors_exp(verbose);
  if (chi.empty()) return 0;

  /* Invert the anamorphosis */

  VectorDouble chi2i = chi2I(chi, 1);
  matrix_invert_copy(chi2i.data(), getNClass(), _i2Chi.data());

  // Update statistics

  calculateMeanAndVariance();

  return 1;
}

int AnamDiscreteDD::_stats(int nech, const VectorDouble& tab)
{
  double zmin,zmax;
  int nclass = getNClass();

  /* Reset the statistics */

  for (int iclass=0; iclass<nclass; iclass++)
  {
    setDDStatProp(iclass,0.);
    setDDStatZmoy(iclass,0.);
  }

  /* Loop on the samples */

  int nactive = 0;
  for (int iech=0; iech<nech; iech++)
  {
    if (FFFF(tab[iech])) continue;
    nactive++;
    for (int iclass=0; iclass<nclass; iclass++)
    {
      zmin = (iclass ==        0) ? 0     : getZCut(iclass-1);
      zmax = (iclass == nclass-1) ? 1.e30 : getZCut(iclass);
      if (tab[iech] <  zmin || tab[iech] >= zmax) continue;
      setDDStatProp(iclass, getDDStatProp(iclass) + 1.);
      setDDStatZmoy(iclass, getDDStatZmoy(iclass) + tab[iech]);
    }
  }
  if (nactive <= 0)
  {
    messerr("No active sample");
    return(1);
  }
  for (int iclass = 0; iclass < nclass; iclass++)
  {
    setDDStatZmoy(iclass, getDDStatZmoy(iclass) / getDDStatProp(iclass));
    setDDStatProp(iclass, getDDStatProp(iclass) / nactive);
  }
  return(0);
}

VectorDouble AnamDiscreteDD::factors_exp(bool verbose)
{
  VectorDouble chi, maf, lambda, veca, vecb, vecc, f1, eigval, eigvec;

  /* Initializations */

  int nclass = getNClass();

  /* Core allocation */

  f1.resize(nclass);
  veca.resize(nclass);
  vecb.resize(nclass);
  vecc.resize(nclass);
  eigvec.resize(nclass * nclass);
  eigval.resize(nclass);

  /* Calculate the experimental MAF array */

  maf = factors_maf(verbose);

  /* Calculate the array 'F1' (based on the first MAF) */

  for (int iclass=0; iclass<nclass; iclass++)
    f1[iclass] = maf[iclass] / maf[0];

  /* Establish the tri-diagonal matrix */

  for (int iclass=0; iclass<nclass; iclass++)
  {
    veca[iclass] = vecb[iclass] = vecc[iclass] = 0.;
    if (iclass < nclass-1)
    {
      for (int jclass=0; jclass<=iclass; jclass++)
        veca[iclass] -= getDDStatProp(jclass) * f1[jclass];
      veca[iclass] /= getDDStatProp(iclass) * (f1[iclass+1] - f1[iclass]);
    }
    if (iclass > 0)
      vecb[iclass] =
        veca[iclass-1] * getDDStatProp(iclass-1) / getDDStatProp(iclass);
    vecc[iclass] = -(veca[iclass] + vecb[iclass]);
  }

  /* Calculate the infinitesimal generator */

  chi = _generator(vecc,veca,vecb,eigvec,eigval);
  if (chi.empty()) return chi;

  /* Calculate the lambda vector from eigen values */

  for (int iclass=0; iclass<nclass; iclass++)
    setDDStatLambda(iclass,-eigval[iclass]);
  lambda.resize(nclass);
  for (int iclass=0; iclass<nclass; iclass++)
    lambda[iclass] = getDDStatLambda(iclass);
  ut_sort_double(0,nclass,NULL,lambda.data());
  for (int iclass=0; iclass<nclass; iclass++)
    setDDStatLambda(iclass, lambda[iclass]);
  setDDStatLambda(0, 0.);
  setDDStatLambda(1, 1.);

  /* Derive the MUL Terms */

  _lambdaToMul();

  /* Calculate the spectrum weighting function */

  for (int iclass=0; iclass<nclass; iclass++)
  {
    double sum = 0.;
    for (int jclass=0; jclass<nclass; jclass++)
      sum += getDDStatProp(jclass) *
        EIGVEC(iclass,jclass) * EIGVEC(iclass,jclass);
    setDDStatU(iclass, getDDStatProp(0) / sum);
  }

  /* Calculate the array of point C_i (normalized polynomials) */

  for (int iclass=0; iclass<nclass; iclass++)
  {
    double value = 0.;
    for (int jclass=0; jclass<nclass; jclass++)
      value += getDDStatZmoy(jclass) * getDDStatProp(jclass) *
        CHI(iclass,jclass);
    setDDStatCnorm(iclass, value);
  }

  /* Verbose option */

  if (verbose)
    print_matrix("Factors",0,1,nclass,nclass,NULL,chi.data());

  return chi;
}

VectorDouble AnamDiscreteDD::factors_maf(bool verbose)
{
  VectorDouble maf, tab;
  int ncut   = getNCut();
  int nclass = getNClass();

  /* Core allocation */

  maf.resize(nclass * nclass,0);
  tab.resize(nclass * nclass,0);

  /* Calculate the experimental MAF array */

  int ecr = 0;
  for (int icut=0; icut<ncut; icut++)
    for (int iclass=0; iclass<nclass; iclass++,ecr++)
    {
      double bval = (iclass >=     icut) ? 1 : 0;
      double cval = (iclass >= (icut+1)) ? 1 : 0;
      double prop = getDDStatProp(icut);
      tab[ecr] = ((bval - cval) - prop) / sqrt(prop * (1. - prop));
    }
  matrix_product(nclass,ncut,ncut,tab.data(),getPcaZ2F().data(),maf.data());

  /* Verbose option */

  if (verbose)
    print_matrix("MAF",0,1,ncut,nclass,NULL,maf.data());

  return maf;
}

/**
 *
 * @param vecc Vector of the Tridiagonal matrix
 * @param veca Vector of the Tridiagonal matrix
 * @param vecb Vector of the Tridiagonal matrix
 * @param eigvec Returned Eigen vectors
 * @param eigval Returned Eigen values
 * @return Calculate the infinitesimal generator
 */
VectorDouble AnamDiscreteDD::_generator(const VectorDouble& vecc,
                                          const VectorDouble& veca,
                                          const VectorDouble& vecb,
                                          VectorDouble& eigvec,
                                          VectorDouble& eigval)
{
  VectorDouble hvar, chi;

  /* Initializations */

  int nclass = getNClass();

  /* Preliminary checks */

  for (int iclass=0; iclass<nclass; iclass++)
  {
    if (veca[iclass] < 0)
    {
      messerr("Diffusion hypothesis invalid: superdiagonal term (class=%d) is not positive",
              iclass+1);
      return chi;
    }
    if (vecb[iclass] < 0)
    {
      messerr("Diffusion hypothesis invalid: subdiagonal term (class=%d) is not positive",
              iclass+1);
      return chi;
    }
  }

  /* Core allocation */

  hvar.resize(nclass);
  chi.resize(nclass * nclass,0);

  /* Diagonalize the infinitesimal generator */

  matrix_eigen_tridiagonal(vecc.data(),vecb.data(),veca.data(),
                           nclass,eigvec.data(),eigval.data());

  /* Choose to set the Hn(0) = 1 */

  for (int iclass=0; iclass<nclass; iclass++)
    for (int jclass=nclass-1; jclass>=0; jclass--)
      EIGVEC(iclass,jclass) /= EIGVEC(iclass,0);

  /* Calculate the statistics on the factors */

  for (int iclass=0; iclass<nclass; iclass++)
  {
    double sum = 0.;
    for (int jclass=0; jclass<nclass; jclass++)
      sum += getDDStatProp(jclass) *
        EIGVEC(iclass,jclass) * EIGVEC(iclass,jclass);
    hvar[iclass] = sum;
  }

  /* Normalize the factors */

  for (int iclass=0; iclass<nclass; iclass++)
    for (int jclass=0; jclass<nclass; jclass++)
    {
      if (iclass == 0)
        CHI(iclass,jclass) = 1.;
      else
        CHI(iclass,jclass) = EIGVEC(iclass,jclass) / sqrt(hvar[iclass]);
    }
  return chi;
}


void AnamDiscreteDD::_lambdaToMul()
{
  int nclass   = getNClass();
  double scoef = getSCoef();
  double mu    = getMu();

  /* Loop on the classes */

  for (int iclass=0; iclass<nclass; iclass++)
    setDDStatMul(iclass, pow(mu / (mu + getDDStatLambda(iclass)),scoef/2.));
}

VectorDouble AnamDiscreteDD::z2factor(double z, const VectorInt& ifacs) const
{
  VectorDouble factors;
  int nfact = (int) ifacs.size();
  factors.resize(nfact,0);

  int nclass = getNClass();
  for (int ifac=0; ifac<nfact; ifac++)
  {
    double value = 0.;
    for (int iclass=0; iclass<nclass; iclass++)
    {
      double zmax   = (iclass == nclass-1) ? 1.e30 : getZCut(iclass);
      value += I2CHI(ifacs[ifac],iclass);
      if (zmax > z) break;
    }
    factors[ifac] = value;
  }
  return factors;
}

VectorDouble AnamDiscreteDD::factors_mod()
{
  VectorDouble chi;
  VectorDouble q2_s, q_s, tri1, tri2, c_s, ptab;
  VectorDouble stats, veca, vecb, vecc, eigval, eigvec;

  /* Initializations */

  int nclass = getNClass();
  int ntri = nclass * (nclass + 1) / 2;

  /* Core allocation */

  q2_s.resize(nclass);
  q_s.resize(nclass * nclass);
  tri1.resize(ntri);
  tri2.resize(ntri);
  ptab.resize(nclass * nclass);
  c_s.resize(nclass * nclass);

  veca.resize(nclass);
  vecb.resize(nclass);
  vecc.resize(nclass);
  eigvec.resize(nclass * nclass);
  eigval.resize(nclass);

  /* Calculate the monomials */

  for (int iclass=0; iclass<nclass; iclass++)
    for (int jclass=0; jclass<nclass; jclass++)
    {
      double value = 0.;
      for (int ic=0; ic<nclass; ic++)
        value += getDDStatU(ic) * pow(getDDStatLambda(ic),iclass);
      PTAB(iclass,jclass) = pow(getDDStatLambda(jclass),iclass) / sqrt(value);
    }

  /* Covariance of monomials in L2(R,u) */

  for (int iclass=0; iclass<nclass; iclass++)
    for (int jclass=0; jclass<nclass; jclass++)
    {
      double value = 0.;
      for (int ic=0; ic<nclass; ic++)
        value += PTAB(iclass,ic) * getDDStatU(ic) * PTAB(jclass,ic);
      C_S(iclass,jclass) = value;
    }

  if (matrix_cholesky_decompose(c_s.data(),tri1.data(),nclass)) return chi;
  matrix_cholesky_invert(nclass,tri1.data(),tri2.data());
  matrix_cholesky_product(2,nclass,nclass,tri2.data(),ptab.data(),q_s.data());

  for (int jclass=nclass-1; jclass>=0; jclass--)
    for (int iclass=0; iclass<nclass; iclass++)
      Q_S(iclass,jclass) /= Q_S(iclass,0);

  double sum = 0.;
  for (int iclass=0; iclass<nclass; iclass++)
  {
    double value = 0.;
    for (int ic=0; ic<nclass; ic++)
      value += Q_S(iclass,ic) * getDDStatU(ic) * Q_S(iclass,ic);
    q2_s[iclass] = value;
    sum += 1. / q2_s[iclass];
  }

  /* Derive the Stationary probabilities */

  for (int iclass=0; iclass<nclass; iclass++)
    setDDStatProp(iclass, (1. / q2_s[iclass]) / sum);

  /* Calculation of the diffusion coefficients */

  for (int iclass=0; iclass<nclass-1; iclass++)
  {
    double local = 0.;
    for (int ic = 0; ic < nclass; ic++)
      local -= (getDDStatLambda(ic)* Q_S(iclass,ic) * getDDStatU(ic) * Q_S(iclass+1,ic));
    veca[iclass] = local / q2_s[iclass + 1];
  }
  veca[nclass-1] = 0.;

  vecb[0] = 0.;
  for (int iclass=1; iclass<nclass; iclass++)
  {
    double local = 0.;
    for (int ic = 0; ic < nclass; ic++)
      local -= (getDDStatLambda(ic) * Q_S(iclass,ic) * getDDStatU(ic) * Q_S(iclass-1,ic));
    vecb[iclass] = local / q2_s[iclass - 1];
  }

  for (int iclass=0; iclass<nclass; iclass++)
    vecc[iclass] = -(veca[iclass] + vecb[iclass]);

  /* Calculate the infinitesimal generator */

  chi = _generator(vecc,veca,vecb,eigvec,eigval);
  return chi;
}

/**
 *
 * @param chi Chi Matrix
 * @param mode Type of recovery function
** \li                  1 : Indicator
** \li                  2 : Metal quantity
** \li                  3 : Benefit
 * @return Calculate the transition matrix from factor to a recovery item
**  i.e. Indicator, Metal Quantity or Benefit (according to mode)
**  (Diffusion Discrete)
 */
VectorDouble AnamDiscreteDD::chi2I(const VectorDouble& chi, int mode)
{
  VectorDouble stats, chi2i, mati;

  /* Initializations */

  int nclass = getNClass();
  chi2i.resize(nclass * nclass,0);
  mati.resize(nclass * nclass,0);

  /* MATI contains the matrix of indicators */

  for (int iclass=0; iclass<nclass; iclass++)
    for (int jclass=0; jclass<nclass; jclass++)
      switch (mode)
      {
        case 1:
          MATI(iclass,jclass) = (jclass >= iclass);
          break;

        case 2:
          MATI(iclass,jclass) = (jclass >= iclass) * getDDStatZmoy(jclass);
          break;

        case 3:
          MATI(iclass,jclass) = (jclass >= iclass) *
            (getDDStatZmoy(jclass) - getDDStatZmoy(iclass));
          break;
      }

  /* Calculate the matrix for CHI_2_I */

  int ecr = 0;
  for (int iclass=0; iclass<nclass; iclass++)
    for (int jclass=0; jclass<nclass; jclass++)
    {
      double value = 0;
      for (int ic=0; ic<nclass; ic++)
        value += MATI(iclass,ic) * getDDStatProp(ic) * CHI(jclass,ic);
      chi2i[ecr++] = value;
    }
  return chi2i;
}

int AnamDiscreteDD::_serialize(std::ostream& os, bool verbose) const
{

  if (AnamDiscrete::_serialize(os, verbose)) return 1;

  bool ret = _recordWrite<double>(os, "Change of support coefficient", getSCoef());
  ret = ret && _recordWrite<double>(os, "Additional Mu coefficient", getMu());
  ret = ret && _tableWrite(os, "PCA Z2Y", getNCut() * getNCut(), getPcaZ2F());
  ret = ret && _tableWrite(os, "PCA Y2Z", getNCut() * getNCut(), getPcaF2Z());

  return ret ? 0 : 1;
}

int AnamDiscreteDD::_deserialize(std::istream& is, bool verbose)
{
  VectorDouble pcaf2z, pcaz2f;
  double s = TEST;
  double mu = TEST;

  if (! AnamDiscrete::_deserialize(is, verbose)) return 1;

  bool ret = _recordRead<double>(is, "Anamorphosis 's' coefficient", s);
  ret = ret && _recordRead<double>(is, "Anamorphosis 'mu' coefficient", mu);
  if (! ret) return 1;
  pcaz2f.resize(getNCut() * getNCut());
  pcaf2z.resize(getNCut() * getNCut());

  if (_tableRead(is, getNCut() * getNCut(), pcaz2f.data())) return 1;
  if (_tableRead(is, getNCut() * getNCut(), pcaf2z.data())) return 1;

  setRCoef(s);
  setMu(mu);
  setPcaF2Z(pcaf2z);
  setPcaZ2F(pcaz2f);

  return 0;
}

double AnamDiscreteDD::modifyCov(const ECalcMember& member,
                                 int iclass,
                                 double dist,
                                 double cov0,
                                 double cov1,
                                 double /*cov2*/) const
{
  double cov;
  double coeff = 0.;
  double mui = getDDStatMul(iclass);
  double gamref = cov0 - cov1;
  if (dist <= 0.)
  {
    switch (member.toEnum())
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
    double li = getDDStatLambda(iclass);
    switch (member.toEnum())
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
  return cov;
}

double AnamDiscreteDD::getBlockVariance(double sval, double /*power*/) const
{
  if (! hasChangeSupport()) return TEST;
  int nclass = getNClass();

  // At this stage (point -> block)) cnorm designate the point C_i

  double var = 0.;
  for (int iclass = 1; iclass < nclass; iclass++)
  {
    double ci = getDDStatCnorm(iclass);
    var += ci * ci * pow(_mu / (_mu + getDDStatLambda(iclass)), sval);
  }
  return (var);
}

int AnamDiscreteDD::updatePointToBlock(double r_coef)
{
  if (! hasChangeSupport()) return 1;
  setRCoef(r_coef);
  int nclass = getNClass();

  /* Update the coefficients mul */

  _lambdaToMul();

  /* Spectral measure */

  double sum = 0.;
  for (int iclass = 0; iclass < nclass; iclass++)
  {
    double mul = getDDStatMul(iclass);
    double newU = getDDStatU(iclass) / (mul * mul);
    setDDStatU(iclass, newU);
    sum += newU;
  }

  for (int iclass = 0; iclass < nclass; iclass++)
  {
    double value = getDDStatU(iclass) / sum;
    setDDStatU(iclass, value);
  }

  /* Update the C_i from point to block */

  for (int iclass = 0; iclass < nclass; iclass++)
  {
    double cnorm = getDDStatCnorm(iclass);
    double mul = getDDStatMul(iclass);
    setDDStatCnorm(iclass, cnorm * mul);
  }

  /* Modeling the diffusion process */

  VectorDouble chi = factors_mod();
  if (chi.empty()) return 1;

  /* Establish the block anamorphosis */

  _blockAnamorphosis(chi);

  return 0;
}

/*****************************************************************************/
/*!
 **  Calculate the block anamorphosis (from point anamorphosis)
 **
 ** \param[in]  chi      Array containing the Chi factors
 **
 *****************************************************************************/
void AnamDiscreteDD::_blockAnamorphosis(const VectorDouble& chi)
{
  int nclass = getNClass();

  /* Block anamorphosis on the indicators */

  for (int iclass = 0; iclass < nclass; iclass++)
  {
    double sum = 0.;
    for (int jclass = 0; jclass < nclass; jclass++)
      sum += getDDStatCnorm(jclass) * CHI(jclass, iclass);
    setDDStatZmoy(iclass, sum);
  }

  /* Update mean and variance */

  calculateMeanAndVariance();

  return;
}

/****************************************************************************/
/*!
 **  Calculate the theoretical grade tonnage value (Discrete Diffusion case)
 **
 ** \param[in] flag_correct 1 if Tonnage order relationship must be corrected
 **
 ** \remark Can only calculate the grade-tonnage curve for the discretization
 ** \remark cutoffs.
 **
 *****************************************************************************/
Selectivity AnamDiscreteDD::calculateSelectivity(bool flag_correct)
{
  int nclass = getNClass();
  Selectivity calest(nclass);

  /* Calculate the Grade-Tonnage curves */

  for (int iclass = 0; iclass < nclass; iclass++)
  {
    double zval = (iclass == nclass - 1) ? 0. : getZCut(nclass - iclass - 2);
    double tval = 0.;
    double qval = 0.;
    for (int jclass = 0; jclass <= iclass; jclass++)
    {
      int ic = nclass - jclass - 1;
      tval += getDDStatProp(ic);
      qval += getDDStatProp(ic) * getDDStatZmoy(ic);
    }
    calest.setZcut(nclass-iclass-1, zval);
    calest.setTest(nclass-iclass-1, tval);
    calest.setQest(nclass-iclass-1, qval);
  }

  /* Correct order relationship */

  if (flag_correct) calest.correctTonnageOrder();

  /* Store the results */

  calest.calculateBenefitGrade();

  return calest;
}
