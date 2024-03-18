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
#include "geoslib_old_f.h"
#include "geoslib_enum.h"
#include "Anamorphosis/AnamDiscreteDD.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Db/Db.hpp"
#include "Stats/Selectivity.hpp"

#include <math.h>

#define EIGVEC(i,j)      eigvec[(i)*nclass+(j)]
#define CHI(i,j)         chi[(i)*nclass+(j)]
#define PTAB(i,j)        ptab[(i)*nclass+(j)]
#define C_S(i,j)         c_s[(i)*nclass+(j)]
#define Q_S(i,j)         q_s[(i)*nclass+(j)]

AnamDiscreteDD::AnamDiscreteDD(double mu, double scoef)
    : AnamDiscrete(),
      _mu(mu),
      _sCoef(scoef),
      _maf(),
      _i2Chi()
{
}

AnamDiscreteDD::AnamDiscreteDD(const AnamDiscreteDD &m)
    : AnamDiscrete(m),
      _mu(m._mu),
      _sCoef(m._sCoef),
      _maf(m._maf),
      _i2Chi(m._i2Chi)
{
}

AnamDiscreteDD& AnamDiscreteDD::operator=(const AnamDiscreteDD &m)
{
  if (this != &m)
  {
    AnamDiscrete::operator=(m);
    _mu = m._mu;
    _sCoef = m._sCoef;
    _maf = m._maf;
    _i2Chi = m._i2Chi;
  }
  return *this;
}

AnamDiscreteDD::~AnamDiscreteDD()
{

}

AnamDiscreteDD* AnamDiscreteDD::createFromNF(const String& neutralFilename, bool verbose)
{
  AnamDiscreteDD* anam = nullptr;
  std::ifstream is;
  anam = new AnamDiscreteDD();
  bool success = false;
  if (anam->_fileOpenRead(neutralFilename, is, verbose))
  {
    success = anam->deserialize(is, verbose);
  }
  if (! success)
  {
    delete anam;
    anam = nullptr;
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
                           const MatrixSquareGeneral &pcaz2f,
                           const MatrixSquareGeneral &pcaf2z,
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
  if (getNCut() <= 0 && getNElem() <= 0) return sstr.str();

  sstr << "Discrete Diffusion Anamorphosis" << std::endl;

  sstr << AnamDiscrete::toString(strfmt);

  if (! _isFitted()) return sstr.str();

  if (_sCoef != 0.)
  {
    sstr << "Mu Coefficient    = " << _mu << std::endl;
    sstr << "Change of Support = " << _sCoef << std::endl;
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

int AnamDiscreteDD::fitFromArray(const VectorDouble& tab,
                                 const VectorDouble& /*wt*/)
{
  VectorDouble chi;

  int nech = static_cast<int> (tab.size());

  // Calculate statistics on data

  _stats(nech,tab);

  // Modeling the diffusion process

  chi = factors_exp();
  if (chi.empty()) return 0;

  /* Invert the anamorphosis */

  _i2Chi = chi2I(chi, 1);
  _i2Chi.invert();

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
  VH::sortInPlace(lambda);
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
  matrix_product_safe(nclass,ncut,ncut,tab.data(),getPcaZ2Fs().getValues().data(),maf.data());

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
VectorDouble AnamDiscreteDD::_generator(const VectorDouble &vecc,
                                        const VectorDouble &veca,
                                        const VectorDouble &vecb,
                                        VectorDouble &eigvec,
                                        VectorDouble &eigval)
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
      value += _i2Chi.getValue(iclass,ifacs[ifac]);
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
MatrixSquareGeneral AnamDiscreteDD::chi2I(const VectorDouble& chi, int mode)
{
  VectorDouble stats;

  /* Initializations */

  int nclass = getNClass();
  MatrixSquareGeneral chi2i(nclass);
  MatrixSquareGeneral mati(nclass);
  chi2i.fill(0.);
  mati.fill(0.);

  /* 'mati' contains the matrix of indicators */

  for (int iclass=0; iclass<nclass; iclass++)
    for (int jclass=0; jclass<nclass; jclass++)
      switch (mode)
      {
        case 1:
          mati.setValue(iclass,jclass,jclass >= iclass);
          break;

        case 2:
          mati.setValue(iclass,jclass,(jclass >= iclass) * getDDStatZmoy(jclass));
          break;

        case 3:
          mati.setValue(iclass,jclass, (jclass >= iclass) *
            (getDDStatZmoy(jclass) - getDDStatZmoy(iclass)));
          break;
      }

  /* Calculate the matrix for CHI_2_I */

  for (int iclass=0; iclass<nclass; iclass++)
    for (int jclass=0; jclass<nclass; jclass++)
    {
      double value = 0;
      for (int ic=0; ic<nclass; ic++)
        value += mati.getValue(iclass,ic) * getDDStatProp(ic) * CHI(jclass,ic);
      chi2i.setValue(iclass, jclass, value);
    }
  return chi2i;
}

bool AnamDiscreteDD::_serialize(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret = ret && AnamDiscrete::_serialize(os, verbose);
  ret = ret && _recordWrite<double>(os, "Change of support coefficient", getSCoef());
  ret = ret && _recordWrite<double>(os, "Additional Mu coefficient", getMu());
  ret = ret && _tableWrite(os, "PCA Z2Y", getNCut() * getNCut(), getPcaZ2Fs().getValues());
  ret = ret && _tableWrite(os, "PCA Y2Z", getNCut() * getNCut(), getPcaF2Zs().getValues());
  return ret;
}

bool AnamDiscreteDD::_deserialize(std::istream& is, bool verbose)
{
  MatrixSquareGeneral pcaf2z, pcaz2f;
  double s = TEST;
  double mu = TEST;

  bool ret = true;
  ret = ret && AnamDiscrete::_deserialize(is, verbose);
  ret = ret && _recordRead<double>(is, "Anamorphosis 's' coefficient", s);
  ret = ret && _recordRead<double>(is, "Anamorphosis 'mu' coefficient", mu);

  int ncut = getNCut();
  if (ret)
  {
    VectorDouble local(ncut * ncut);
    ret = ret && _tableRead(is, "PCA Z2Y", getNCut() * getNCut(), local.data());
    pcaz2f.resetFromVD(ncut, ncut, local);
  }

  if (ret)
  {
    VectorDouble local(ncut * ncut);
    ret = ret && _tableRead(is, "PCA Y2Z", getNCut() * getNCut(), local.data());
    pcaf2z.resetFromVD(ncut, ncut, local);
  }

  if (ret)
  {
    setRCoef(s);
    setMu(mu);
    setPcaF2Z(pcaf2z);
    setPcaZ2F(pcaz2f);
  }
  return ret;
}

double AnamDiscreteDD::computeVariance(double sval) const
{
  if (! allowChangeSupport()) return TEST;
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
  if (! allowChangeSupport()) return 1;
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
 ** \param[in] selectivity Selectivity structure to be filled
 **
 *****************************************************************************/
void AnamDiscreteDD::_globalSelectivity(Selectivity* selectivity)
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
  int ncut = selloc->getNCuts();

  /* Calculate the Grade-Tonnage curves */

  for (int icut = 0; icut < ncut; icut++)
  {
    double zval = (icut == ncut - 1) ? 0. : getZCut(ncut - icut - 2);
    double tval = 0.;
    double qval = 0.;
    for (int jclass = 0; jclass <= icut; jclass++)
    {
      int ic = ncut - jclass - 1;
      tval += getDDStatProp(ic);
      qval += getDDStatProp(ic) * getDDStatZmoy(ic);
    }
    selloc->setZcut(ncut-icut-1, zval);
    selloc->setTest(ncut-icut-1, tval);
    selloc->setQest(ncut-icut-1, qval);
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
 **  Case of Discrete Diffusion
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  cols_est     Array of UIDs for factor estimation
 ** \param[in]  cols_std     Array of UIDs for factor st. dev.
 ** \param[in]  iptr0        Rank for storing the results
 **
 *****************************************************************************/
int AnamDiscreteDD::factor2Selectivity(Db *db,
                                       Selectivity* selectivity,
                                       const VectorInt& cols_est,
                                       const VectorInt& cols_std,
                                       int iptr0)
{
  int nclass   = getNClass();
  int nech     = db->getSampleNumber();
  int nb_est   = (int) cols_est.size();
  int nb_std   = (int) cols_std.size();
  int ncleff   = MAX(nb_est, nb_std);
  bool cutDefined = (selectivity->getNCuts() > 0);

  /* Preliminary checks */

  if (db == nullptr)
  {
    messerr("You must define a Db");
    return 1;
  }
  int nvar = MAX(nb_est, nb_std);

  /* Get the number of initial cutoffs */

  int nmax = getNClass();
  if (nvar >= nmax)
  {
    messerr("Number of factors (%d) must be smaller than Number of classes (%d)",
            nvar, nmax);
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

  /* Modeling the diffusion process */

  VectorDouble chi = factors_mod();
  if (chi.empty()) return 1;
  MatrixSquareGeneral ct = chi2I(chi, 1);
  MatrixSquareGeneral cq = chi2I(chi, 2);
  MatrixSquareGeneral cb = chi2I(chi, 3);

  /* Calculate the Recovery Functions from the factors */

  for (int iech = 0; iech < nech; iech++)
  {
    if (_isSampleSkipped(db, iech, cols_est, cols_std)) continue;

    /* Tonnage: Estimation */

    for (int iclass = 0; iclass < ncleff; iclass++)
    {
      double total = ct.getValue(0, iclass);
      for (int ivar = 0; ivar < nb_est; ivar++)
      {
        double value = db->getArray(iech, cols_est[ivar]);
        total += value * ct.getValue(ivar + 1, iclass);
      }
      selloc->setTest(iclass, total);
    }

    /* Correct Order relationship for Tonnage (optional) */

    selloc->correctTonnageOrder();

    /* Tonnage: Standard Deviation */

    if (selectivity->isUsedStD(ESelectivity::T))
    {
      for (int iclass = 0; iclass < ncleff; iclass++)
      {
        double total = 0.;
        for (int ivar = 0; ivar < ncleff - 1; ivar++)
        {
          double value = (ivar < nb_std) ? db->getArray(iech, cols_std[ivar]) : 1.;
          double prod = value * ct.getValue(ivar + 1, iclass);
          total += prod * prod;
        }
        selloc->setTstd(iclass, sqrt(total));
      }
    }

    /* Metal Quantity: Estimation */

    if (selectivity->isUsedEst(ESelectivity::Q))
    {
      selloc->setQest(nclass-1, getDDStatZmoy(ncleff - 1)
                     * selloc->getTest(ncleff - 1));
      for (int iclass = ncleff - 2; iclass >= 0; iclass--)
        selloc->setQest(iclass, selloc->getQest(iclass + 1)
                       + getDDStatZmoy(iclass) *
                       (selloc->getTest(iclass) - selloc->getTest(iclass + 1)));
    }

    /* Metal Quantity: Standard Deviation */

    if (selectivity->isUsedStD(ESelectivity::Q))
    {
      for (int iclass = 0; iclass < ncleff; iclass++)
      {
        double total = 0.;
        for (int ivar = 0; ivar < ncleff - 1; ivar++)
        {
          double value = (ivar < nb_std) ?
              db->getArray(iech, cols_std[ivar]) : 1.;
          double prod = value * cq.getValue(ivar + 1, iclass);
          total += prod * prod;
        }
        selloc->setQstd(iclass, sqrt(total));
      }
    }

    /* Z: Estimation */

    double zestim = 0.;
    if (selectivity->isUsedEst(ESelectivity::Z))
    {
      zestim = getDDStatZmoy(ncleff - 1) * selloc->getTest(ncleff - 1);
      for (int iclass = 0; iclass < ncleff - 1; iclass++)
        zestim += getDDStatZmoy(iclass)
            * (selloc->getTest(iclass) - selloc->getTest(iclass + 1));
    }

    /* Z: Standard Deviation */

    double zstdev = 0.;
    if (selectivity->isUsedStD(ESelectivity::Z))
    {
      double total = 0;
      for (int ivar = 0; ivar < ncleff - 1; ivar++)
      {
        double value = (ivar < nb_std) ?
            db->getArray(iech, cols_std[ivar]) : 1.;
        double prod = value * getDDStatCnorm(ivar);
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
  return 0;
}

