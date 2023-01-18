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
#include "geoslib_old_f.h"
#include "geoslib_enum.h"

#include "Enum/ECalcMember.hpp"

#include "Space/ASpaceObject.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Anamorphosis/AnamContinuous.hpp"
#include "Polynomials/Hermite.hpp"
#include "Basic/Interval.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Basic/ASerializable.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovLMC.hpp"
#include "Stats/Selectivity.hpp"

#include <math.h>

#define ANAM_YMIN -10.
#define ANAM_YMAX  10.
#define YPAS       0.1

AnamHermite::AnamHermite(int nbpoly, bool flagBound, double rCoef)
    : AnamContinuous(),
      _flagBound(flagBound),
      _rCoef(rCoef),
      _psiHn()
{
  _psiHn.resize(nbpoly);
}

AnamHermite::AnamHermite(const AnamHermite &m)
    : AnamContinuous(m),
      _flagBound(m._flagBound),
      _rCoef(m._rCoef),
      _psiHn(m._psiHn)
{
}

AnamHermite& AnamHermite::operator=(const AnamHermite &m)
{
  if (this != &m)
  {
    AnamContinuous::operator=(m);
    _flagBound = m._flagBound;
    _rCoef = m._rCoef;
    _psiHn = m._psiHn;
  }
  return *this;
}

AnamHermite::~AnamHermite()
{
}

String AnamHermite::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  int nbpoly = getNbPoly();
  if (nbpoly <= 0) return sstr.str();

  sstr << toTitle(1,"Hermitian Anamorphosis");

  sstr << AnamContinuous::toString(strfmt);

  sstr << "Number of Hermite polynomials = " << nbpoly << std::endl;
  if (isChangeSupportDefined())
    sstr << "Change of Support Coefficient = " << _rCoef << std::endl;

  if (! _isFitted()) return sstr.str();

  sstr << toVector("Normalized coefficients for Hermite polynomials (punctual variable)",
                   _psiHn);

  return sstr.str();
}

AnamHermite* AnamHermite::createFromNF(const String& neutralFilename, bool verbose)
{
  AnamHermite* anam = nullptr;
  std::ifstream is;
  anam = new AnamHermite();
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

AnamHermite* AnamHermite::create(int nbpoly, bool flagBound, double rCoef)
{
  return new AnamHermite(nbpoly, flagBound, rCoef);
}

void AnamHermite::reset(double pymin,
                        double pzmin,
                        double pymax,
                        double pzmax,
                        double aymin,
                        double azmin,
                        double aymax,
                        double azmax,
                        double r,
                        const VectorDouble &psi_hn)
{
  setPsiHns(psi_hn);
  setRCoef(r);
  calculateMeanAndVariance();
  setABounds(azmin, azmax, aymin, aymax);
  setPBounds(pzmin, pzmax, pymin, pymax);
}

double AnamHermite::RawToTransformValue(double z) const
{
  double y,y1,y2,yg,z1,z2,zg,dz,dzmax,dy,dymax;
  int i,iter;

  /* Initializations */

  if (getNbPoly() < 1) return(TEST);
  if (FFFF(z)) return TEST;

  /* Check the bounds */

  y = y1 = y2 = z1 = z2 = 0.;
  if (_flagBound)
  {
    if (_az.isOutsideBelow(z)) return(_ay.getVmin());
    if (_az.isOutsideAbove(z)) return(_ay.getVmax());
    if (_pz.isOutsideBelow(z))
    {
      if(_pz.getVmin() == _az.getVmin()) return(_py.getVmin());
      return(_ay.getVmin() + (_py.getVmin() - _ay.getVmin()) *
          (z - _az.getVmin()) / (_pz.getVmin() - _az.getVmin()));
    }

    if (_pz.isOutsideAbove(z))
    {
      if(_pz.getVmax() == _az.getVmax()) return(_py.getVmax());
      return(_ay.getVmax() + (_py.getVmax() - _ay.getVmax()) *
          (z - _az.getVmax()) / (_pz.getVmax() - _az.getVmax()));
    }
  }

  /* Calculate the precision on Z */

  z1 = TransformToRawValue(-1);

  z2 = TransformToRawValue( 1);
  dzmax = ABS((z2 - z1)/100000.);

  /* Look for a first interval in Y containing Z */

  dy = YPAS;
  y1 = 0.;
  z1 = TransformToRawValue(y1);

  if (z > z1)
  {
    for( i=0 ; i<101 ; i++ )
    {
      y2 = y1 + dy;
      z2 = TransformToRawValue(y2);
      if(z2 > z) break;
      y1 = y2;
      z1 = z2;
    }
    if (y1 > ANAM_YMAX) return(ANAM_YMAX + 1.);
  }
  else
  {
    y2 = y1;
    z2 = z1;
    for( i=0 ; i<101 ; i++ )
    {
      y1 = y2 - dy;
      z1 = TransformToRawValue(y1);
      if(z1 < z) break;
      y2 = y1;
      z2 = z1;
    }
    if (y1 < ANAM_YMIN) return(ANAM_YMIN - 1.);
  }

  dz = z2 - z1;
  dy = 1.;
  iter = 0;
  dymax = 0.0000001;

  while( iter<1000000 && dz>dzmax && dy>dymax )
  {
    yg = (y1 + y2)/2.;
    zg = TransformToRawValue(yg);

    if(zg > z)
    {
      z2 = zg;
      y2 = yg;
    }
    else
    {
      z1 = zg;
      y1 = yg;
    }

    dy = y2 - y1;
    dz = z2 - z1;
    iter++;
  }

  dz = z2 - z1;
  if(dz == 0.)
    y = (y1 + y2)/2.;
  else
    y = y1 + (z - z1)*(y2 - y1)/dz;

  if (_flagBound)
  {
    if (y < _ay.getVmin()) y = _ay.getVmin();
    if (y > _ay.getVmax()) y = _ay.getVmax();
  }
  return(y);
}

double AnamHermite::TransformToRawValue(double y) const
{
  double z;
  if (getNbPoly() < 1) return(TEST);
  if (FFFF(y)) return TEST;

  /* Check the bounds */

  if (_flagBound)
  {
    if (_ay.isOutsideBelow(y)) return(_az.getVmin());
    if (_ay.isOutsideAbove(y)) return(_az.getVmax());

    if (_py.isOutsideBelow(y))
    {
      if (_py.getVmin() == _ay.getVmin()) return(_pz.getVmin());
      return(_az.getVmin() + (_pz.getVmin() - _az.getVmin()) *
          (y - _ay.getVmin()) / (_py.getVmin() - _ay.getVmin()));
    }

    if (_py.isOutsideAbove(y))
    {
      if (_py.getVmax() == _ay.getVmax()) return(_pz.getVmax());
      return(_az.getVmax() + (_pz.getVmax() - _az.getVmax()) *
          (y - _ay.getVmax()) / (_py.getVmax() - _ay.getVmax()));
    }
  }

  /* Normal inversion */

  z = hermiteCondExpElement(y, 0., getPsiHns());

  /* Truncate within the bounds */

  if (_flagBound)
  {
    if (z < _az.getVmin()) z = _az.getVmin();
    if (z > _az.getVmax()) z = _az.getVmax();
  }

  return(z);
}

/**
 * Compute the Gaussian covariance from Raw covariance: Sum_n psi_n^2 C^n
 * @param chh
 * @return
 */
double AnamHermite::computeVariance(double chh) const
{
  int nbpoly = getNbPoly();
  double rho = 1.;
  double var = 0.;
  for (int ih = 1; ih < nbpoly; ih++)
  {
    rho *= chh;
    var += getPsiHn(ih) * getPsiHn(ih) * rho;
  }
  return var;
}

void AnamHermite::calculateMeanAndVariance()
{
  _mean = _psiHn[0];
  _variance = computeVariance(1.);
}

void AnamHermite::setRCoef(double r_coef)
{
  _rCoef = r_coef;
  calculateMeanAndVariance();
}

int AnamHermite::fitFromArray(const VectorDouble& tab, const VectorDouble& wt)
{
  int     icl,ih,ncl;
  double  Gcy1,Gcy2,Gy1,Gy2;
  VectorDouble psi, h1, h2, zs, ys;

  int nech = static_cast<int> (tab.size());
  if (nech <= 0) return 1;

  int nbpoly = getNbPoly();
  zs.resize(nech+2);
  ys.resize(nech+2);
  _psiHn.resize(nbpoly, 0.);

  /* Sort the data by classes */

  ncl = _data_sort(nech,tab,wt,zs,ys);
  if (ncl <= 0) return 1;

  /* Calculate Hermite_0 coefficients */

  Gcy1 = 0.;
  for (icl=0; icl<ncl; icl++)
  {
    Gcy2 = law_cdf_gaussian(ys[icl]);
    _psiHn[0] += zs[icl] * (Gcy2 - Gcy1);
    Gcy1 = Gcy2;
  }

  /* Calculate the Hermite coefficients */

  h1 = hermitePolynomials(ys[0],1.,nbpoly);
  Gy1 = 0.;

  for (icl=0 ; icl<ncl ; icl++)
  {
    h2 = hermitePolynomials(ys[icl],1.,nbpoly);

    Gy2 = law_df_gaussian(ys[icl]);

    for( ih=1 ; ih<nbpoly ; ih++ )
      _psiHn[ih] += zs[icl] * (h2[ih-1]*Gy2 - h1[ih-1]*Gy1) / sqrt((double) ih);

    Gy1 = Gy2;
    for( ih=0 ; ih<nbpoly ; ih++) h1[ih] = h2[ih];
  }

  /* Ultimate calculations */

  calculateMeanAndVariance();
  _defineBounds(ys[0], zs[0], ys[ncl - 2] + EPSILON5, zs[ncl - 1], _ay.getVmin(),
                _az.getVmin(), _ay.getVmax(), _az.getVmax());

  return 0;
}

double AnamHermite::getPsiHn(int ih) const
{
  if (! _isIndexValid(ih)) return TEST;
  double value = _psiHn[ih];
  if (isChangeSupportDefined())
    value *= pow(_rCoef, (double) ih);
  return value;
}

VectorDouble AnamHermite::getPsiHns() const
{
  if (isChangeSupportDefined())
  {
    VectorDouble psi = _psiHn;
    double rval = 1.;
    for (int ih = 1; ih < getNbPoly(); ih++)
    {
      rval *= _rCoef;
      psi[ih] *= rval;
    }
    return psi;
  }
  else
    return _psiHn;
}

void AnamHermite::setPsiHn(int i, double psi_hn)
{
  if (! _isIndexValid(i)) return;
  _psiHn[i] = psi_hn;
}

bool AnamHermite::_isIndexValid(int i) const
{
  int nbpoly = getNbPoly();
  if (i < 0 || i >= nbpoly)
  {
    mesArg("Hermite Polynomial Index",i,nbpoly);
    return false;
  }
  return true;
}

void AnamHermite::_defineBounds(double pymin,
                                 double pzmin,
                                 double pymax,
                                 double pzmax,
                                 double aymin,
                                 double azmin,
                                 double aymax,
                                 double azmax)
{
  int npas,ind,ind0;
  VectorDouble ym,zm;

  // Switch off the flagBound during the calculation of bounds

  bool flagBoundMemo = getFlagBound();
  setFlagBound(false);

  /* Initializations */

  npas = (int)((ANAM_YMAX - ANAM_YMIN) / YPAS) + 1;
  if (FFFF(azmin)) azmin = pzmin;
  if (FFFF(aymin)) aymin = pymin;
  if (FFFF(azmax)) azmax = pzmax;
  if (FFFF(aymax)) aymax = pymax;

  /* Core allocation */

  ym.resize(npas + 2);
  zm.resize(npas + 2);

  /* Evaluate the experimental anamorphosis */

  ind0 = (npas - 1) / 2;
  ym[ind0] = 0.;
  zm[ind0] = TransformToRawValue(ym[ind0]);
  for (ind=ind0-1; ind>=0; ind--)
  {
    ym[ind] = ym[ind+1] - YPAS;
    zm[ind] = TransformToRawValue(ym[ind]);
  }
  for (ind=ind0+1; ind<npas; ind++)
  {
    ym[ind] = ym[ind-1] + YPAS;
    zm[ind] = TransformToRawValue(ym[ind]);
  }

  /* Look for a starting search point */

  for (ind0=0; ind0<npas; ind0++)
    if (ym[ind0] > pymin) break;
  for (; ind0 < npas; ind0++)
    if (zm[ind0] > (pzmin+pzmax) / 2.) break;
  if (ind0 == npas) ind0 = (int)(npas / 2.);

  /* Look for the first non-monotonous point, starting from the median */

  for (ind=ind0; ind>0; ind--)
  {
    if (zm[ind] < azmin)
    {
      _az.setVmin(zm[ind+1]);
      _ay.setVmin(RawToTransformValue(_az.getVmin()));
      break;
    }
    else if (ind == 0)
      break;
    else if (FFFF(_pz.getVmin()) && zm[ind-1] > zm[ind])
    {
      _py.setVmin(ym[ind]);
      _pz.setVmin(zm[ind]);
    }
  }
  for (ind=ind0; ind<npas-1; ind++)
  {
    if (zm[ind] > azmax)
    {
      _az.setVmax(zm[ind]);
      _ay.setVmax(RawToTransformValue(_az.getVmax()));
      break;
    }
    else if (ind == npas-1)
      break;
    else if (FFFF(_pz.getVmax()) && zm[ind] > zm[ind+1])
    {
      _py.setVmax(ym[ind]);
      _pz.setVmax(zm[ind]);
    }
  }

  /* Final assignments */

  if (FFFF(_az.getVmin()))
  {
    _ay.setVmin(ANAM_YMIN);
    _az.setVmin(TransformToRawValue(_ay.getVmin()));
  }
  if (FFFF(_pz.getVmin()))
  {
    _py.setVmin(MAX(ANAM_YMIN,_ay.getVmin()));
    _pz.setVmin(TransformToRawValue(_py.getVmin()));
  }
  if (FFFF(_az.getVmax()))
  {
    _ay.setVmax(ANAM_YMAX);
    _az.setVmax(TransformToRawValue(_ay.getVmax()));
  }
  if (FFFF(_pz.getVmax()))
  {
    _py.setVmax(MIN(ANAM_YMAX,_ay.getVmax()));
    _pz.setVmax(TransformToRawValue(_py.getVmax()));
  }

  // Set the FlagBound to its original status
  setFlagBound(flagBoundMemo);
  return;
}

int AnamHermite::_data_sort(int nech,
                            const VectorDouble& z,
                            const VectorDouble& wt,
                            VectorDouble& zs,
                            VectorDouble& ys)
{
  double *tmp, sum, frc, eps, wgt;
  int    *ind,i,ncl,nval;

  /* Initializations */

  frc = ncl = nval = 0;
  tmp = nullptr;
  ind = nullptr;

  /* Copy the variable in arrays zs and ys eliminating undefined values */

  sum = 0.;
  for( i=nval=0 ; i<nech ; i++)
  {
    if (FFFF(z[i])) continue;
    zs[nval] = z[i];
    wgt = 1.;
    if (! wt.empty())
    {
      if (FFFF(wt[i])) continue;
      if (wt[i] <= 0.) continue;
      wgt = wt[i];
    }
    ys[nval] = wgt;
    sum += wgt;
    nval++;
  }
  if (nval <= 0) return(ncl);

  /* Sorting the data */

  if (!wt.empty())
  {
    tmp = (double *) mem_alloc(sizeof(double) * nval,0);
    if (tmp == nullptr) goto label_end;
    ind = (int    *) mem_alloc(sizeof(int)    * nval,0);
    if (ind == nullptr) goto label_end;
    for (i = 0; i < nval; i++)  ind[i] = i;
    ut_sort_double(0,nval,ind,zs.data());
    for (i = 0; i < nval; i++) tmp[i] = ys[ind[i]];
    for (i = 0; i < nval; i++) ys[i]  = tmp[i];
    tmp = (double *) mem_free((char *) tmp);
    ind = (int    *) mem_free((char *) ind);
  }
  else
  {
    ut_sort_double(0,nval,NULL,zs.data());
  }

  /* Loop on the data */

  eps = EPSILON5 * (zs[nval-1] - zs[0]);
  for( i=0 ; i<nval-1 ; i++)
  {
    frc += ys[i];
    if( zs[i] < zs[i+1] )
    {
      zs[ncl] = zs[i];
      ys[ncl] = law_invcdf_gaussian(frc/sum);
      ncl++;
    }
  }

  /* Adding the upper bound */

  zs[ncl] = zs[nval-1];
  ys[ncl] = ys[ncl-1] + .5;
  ncl++;
  zs[ncl] = zs[ncl-1] + eps;
  ys[ncl] = ANAM_YMAX + 1.;
  ncl++;

  /* Adding the lower bound */

  for( i=ncl ; i>0 ; i-- )
  {
    ys[i] = ys[i-1];
    zs[i] = zs[i-1];
  }
  ncl++;
  zs[0] = zs[1] - eps;
  ys[0] = ys[1] - .5;

label_end:
  tmp = (double *) mem_free((char *) tmp);
  ind = (int    *) mem_free((char *) ind);
  return(ncl);
}

bool AnamHermite::_serialize(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret && ret && AnamContinuous::_serialize(os, verbose);
  ret = ret && _recordWrite<int>(os, "Number of Hermite Polynomials", getNbPoly());
  ret = ret && _recordWrite<double>(os,"Change of support coefficient", getRCoef());
  ret = ret && _tableWrite(os, "Hermite Polynomial", getNbPoly(), getPsiHns());
  return ret;
}

bool AnamHermite::_deserialize(std::istream& is, bool verbose)
{
  VectorDouble hermite;
  double r = TEST;
  int nbpoly = 0;

  bool ret = true;

  ret = ret && AnamContinuous::_deserialize(is, verbose);
  ret = _recordRead<int>(is, "Number of Hermite Polynomials", nbpoly);
  if (ret) hermite.resize(nbpoly);
  ret = ret && _recordRead<double>(is, "Change of Support Coefficient", r);
  if (ret) setRCoef(r);
  ret = ret && _tableRead(is, nbpoly, hermite.data());

  return ret;
}

VectorDouble AnamHermite::z2factor(double z, const VectorInt& ifacs) const
{
  return hermitePolynomials(z, 1., ifacs);
}

int AnamHermite::updatePointToBlock(double r_coef)
{
  if (! allowChangeSupport()) return 1;
  setRCoef(r_coef);

  /* Update mean and variance */

  calculateMeanAndVariance();

  return 0;
}

/****************************************************************************/
/*!
 **  Calculate the theoretical grade tonnage value (Gaussian case)
 **
 *****************************************************************************/
void AnamHermite::_globalSelectivity(Selectivity* selectivity)
{
  int nbpoly = getNbPoly();
  setFlagBound(0);
  int ncut = selectivity->getNCuts();

  /* Loop on the cutoff values */

  for (int icut = 0; icut < ncut; icut++)
  {
    double zval = selectivity->getZcut(icut);
    double yval = RawToTransformValue(zval);
    double tval = 1. - law_cdf_gaussian(yval);
    double gval = law_df_gaussian(yval);
    VectorDouble hn = hermitePolynomials(yval, 1., nbpoly);
    double qval = getPsiHn(0) * (1. - law_cdf_gaussian(yval));
    for (int ih = 1; ih < nbpoly; ih++)
      qval -= getPsiHn(ih) * hn[ih - 1] * gval / sqrt((double) ih);
    selectivity->setTest(icut, zval);
    selectivity->setTest(icut, tval);
    selectivity->setQest(icut, qval);
  }
  selectivity->calculateBenefitAndGrade();
}

/*****************************************************************************/
/*!
 **  Calculate Experimental Grade-Tonnage curves from factors
 **  Case of Hermite Anamorphosis
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
int AnamHermite::factor2Selectivity(Db *db,
                                    Selectivity* selectivity,
                                    const VectorInt& cols_est,
                                    const VectorInt& cols_std,
                                    int iptr0)
{
  setFlagBound(1);
  int nbpoly = getNbPoly();
  bool need_T = selectivity->isNeededT();
  bool need_Q = selectivity->isNeededQ();
  int ncut = selectivity->getNCuts();
  int nb_est = (int) cols_est.size();
  int nb_std = (int) cols_std.size();

  /* Preliminary checks */

  if (db == nullptr)
  {
    messerr("You must define a Db");
    return 1;
  }

  /* Get the number of initial cutoffs */

  int nfactor = MAX(nb_est, nb_std);
  if (nfactor >= getNbPoly())
  {
    messerr("Number of Factors (%d) must be smaller than Number of Hermite polynomials (%d)",
        nfactor, getNbPoly());
    return 1;
  }

  /* Loop on the samples */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (_isSampleSkipped(db, iech, cols_est, cols_std)) continue;

    /* Z: Estimation */

    double zestim = 0.;
    if (selectivity->isUsedEst(ESelectivity::Z))
    {
      double total = getPsiHn(0);
      for (int ivar = 0; ivar < nb_est; ivar++)
      {
        double value = db->getArray(iech, cols_est[ivar]);
        double coeff = getPsiHn(ivar + 1);
        total += coeff * value;
      }
      zestim = total;
    }

    /* Z: Standard Deviation */

    double zstdev = 0.;
    if (selectivity->isUsedStD(ESelectivity::Z))
    {
      double total = 0.;
      for (int ivar = 0; ivar < nb_std; ivar++)
      {
        double value = db->getArray(iech, cols_std[ivar]);
        double coeff = getPsiHn(ivar + 1);
        total += coeff * coeff * value;
      }
      zstdev = sqrt(total);
    }

    /* Loop on the cutoffs */

    for (int icut = 0; icut < ncut; icut++)
    {
      double yc = RawToTransformValue(selectivity->getZcut(icut));
      if (need_T)
      {
        VectorDouble s_cc = hermiteCoefIndicator(yc, nbpoly);

        /* Tonnage estimation */

        if (selectivity->isUsedEst(ESelectivity::T))
        {
          double total = s_cc[0];
          for (int ivar = 0; ivar < nb_est; ivar++)
          {
            double value = db->getArray(iech, cols_est[ivar]);
            total += s_cc[ivar + 1] * value;
          }
          selectivity->setTest(icut, total);
        }

        /* Tonnage: Standard Deviation */

        if (selectivity->isUsedStD(ESelectivity::T))
        {
          double total = 0.;
          for (int ivar = 0; ivar < nb_std; ivar++)
          {
            double value = db->getArray(iech, cols_std[ivar]);
            total += s_cc[ivar + 1] * s_cc[ivar + 1] * value;
          }
          selectivity->setTstd(icut, sqrt(total));
        }
      }

      if (need_Q)
      {
        MatrixSquareGeneral TAU = hermiteIncompleteIntegral(yc, nbpoly);

        /* Metal Quantity: Estimation */

        if (selectivity->isUsedEst(ESelectivity::Q))
        {
          double total = 0.;
          for (int ivar = 0; ivar < nb_est; ivar++)
          {
            double value = db->getArray(iech, cols_est[ivar]);
            double fn = 0.;
            for (int jvar = 0; jvar < nbpoly; jvar++)
            {
              double coeff = getPsiHn(jvar);
              fn += coeff * TAU.getValue(ivar, jvar);
            }
            total += fn * value;
          }
          selectivity->setQest(icut, total);
        }

        /* Metal Quantity: Standard Deviation */

        if (selectivity->isUsedStD(ESelectivity::Q))
        {
          double total = 0.;
          for (int ivar = 0; ivar < nb_std; ivar++)
          {
            double value = db->getArray(iech, cols_std[ivar]);
            double fn = 0.;
            for (int jvar = 0; jvar < nbpoly; jvar++)
            {
              double coeff = getPsiHn(jvar);
              fn += coeff * TAU.getValue(ivar, jvar);
            }
            total += fn * fn * value;
          }
          selectivity->setQstd(icut, sqrt(total));
        }
      }
    }

    /* Storage */

    selectivity->calculateBenefitAndGrade();
    selectivity->storeInDb(db, iech, iptr0, zestim, zstdev);
  }
  return (0);
}

double AnamHermite::evalSupportCoefficient(int option,
                                           Model *model,
                                           const VectorDouble &dxs,
                                           const VectorInt &ndisc,
                                           const VectorDouble& angles,
                                           bool verbose)
{
  // Dispatch

  if (option == 1)
  {

    // DGM1 Method

    model->setActiveFactor(0); // Z variable
    double cvv = model->evalCvv(dxs, ndisc, angles);
    double r1  = sqrt(invertVariance(cvv));
    if (verbose)
      message("Change of Support coefficient (DGM-1) = %6.3lf\n", r1);
    return r1;
  }

  if (option == 2)
  {
    model->setActiveFactor(1); // Y Variable
    double cvv = model->evalCvv(dxs, ndisc, angles);
    double r2 = sqrt(cvv);
    if (verbose)
      message("Change of Support coefficient (DGM2) = %6.3lf\n",r2);
    return r2;
  }

  messerr("The argument 'option'(%d) should be 1 or 2",option);
  return TEST;
}
