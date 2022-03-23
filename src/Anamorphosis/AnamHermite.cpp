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
#include "Anamorphosis/AnamHermite.hpp"
#include "Anamorphosis/AnamContinuous.hpp"
#include "Polynomials/Hermite.hpp"
#include "Basic/Interval.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Basic/ASerializable.hpp"
#include "Db/Db.hpp"
#include "Covariances/ECalcMember.hpp"

#include <math.h>

#define ANAM_YMIN -10.
#define ANAM_YMAX  10.
#define YPAS       0.1

AnamHermite::AnamHermite(int nbpoly, bool flagBound, double rCoef)
    : AnamContinuous(),
      _nbPoly(nbpoly),
      _flagBound(flagBound),
      _rCoef(rCoef),
      _psiHn()
{
  _psiHn.resize(nbpoly);
}

AnamHermite::AnamHermite(const AnamHermite &m)
    : AnamContinuous(m),
      _nbPoly(m._nbPoly),
      _flagBound(m._flagBound),
      _rCoef(m._rCoef),
      _psiHn(m._psiHn)
{

}

AnamHermite& AnamHermite::operator=(const AnamHermite &m)
{
  if (this != &m)
  {
    _nbPoly = m._nbPoly;
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

  sstr << toTitle(1,"Hermitian Anamorphosis");

  sstr << AnamContinuous::toString(strfmt);

  sstr << "Number of Hermite polynomials = " << _nbPoly << std::endl;
  if (_rCoef > 0. && _rCoef < 1.)
    sstr << "Change of Support Coefficient = " << _rCoef << std::endl;
  sstr << toVector("Normalized coefficients for Hermite polynomials",_psiHn);

  return sstr.str();
}

int AnamHermite::dumpToNF(const String& neutralFilename, bool verbose) const
{
  std::ofstream os;
  int ret = 1;
  if (_fileOpenWrite(neutralFilename, "AnamHermite", os, verbose))
  {
    ret = _serialize(os, verbose);
    if (ret && verbose) messerr("Problem writing in the Neutral File.");
    os.close();
  }
  return ret;
}

AnamHermite* AnamHermite::createFromNF(const String& neutralFilename, bool verbose)
{
  AnamHermite* anam = nullptr;
  std::ifstream is;
  if (_fileOpenRead(neutralFilename, "AnamHermite", is, verbose))
  {
    anam = new AnamHermite();
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


AnamHermite* AnamHermite::create(int nbpoly, bool flagBound, double rCoef)
{
  return new AnamHermite(nbpoly, flagBound, rCoef);
}

double AnamHermite::RawToGaussianValue(double z) const
{
  double y,y1,y2,yg,z1,z2,zg,dz,dzmax,dy,dymax;
  int i,iter;

  /* Initializations */

  if (_nbPoly < 1) return(TEST);

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

  z1 = GaussianToRawValue(-1);

  z2 = GaussianToRawValue( 1);
  dzmax = ABS((z2 - z1)/100000.);

  /* Look for a first interval in Y containing Z */

  dy = YPAS;
  y1 = 0.;
  z1 = GaussianToRawValue(y1);

  if (z > z1)
  {
    for( i=0 ; i<101 ; i++ )
    {
      y2 = y1 + dy;
      z2 = GaussianToRawValue(y2);
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
      z1 = GaussianToRawValue(y1);
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
    zg = GaussianToRawValue(yg);

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

double AnamHermite::GaussianToRawValue(double y) const
{
  double z;
  if (_nbPoly < 1) return(TEST);

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

  z = hermiteCondExpElement(y, 0., _psiHn);

  /* Truncate within the bounds */

  if (_flagBound)
  {
    if (z < _az.getVmin()) z = _az.getVmin();
    if (z > _az.getVmax()) z = _az.getVmax();
  }

  return(z);
}

double AnamHermite::calculateVarianceFromPsi(double chh)
{
  double rho = 1.;
  double var = 0.;
  for (int ih = 1; ih < _nbPoly; ih++)
  {
    rho *= chh;
    var += _psiHn[ih] * _psiHn[ih] * rho;
  }
  return var;
}

void AnamHermite::calculateMeanAndVariance()
{
  _mean = _psiHn[0];
  _variance = calculateVarianceFromPsi(1.);
}

int AnamHermite::fit(const VectorDouble& tab, const VectorDouble& wt)
{
  int     icl,ih,ncl;
  double  Gcy1,Gcy2,Gy1,Gy2;
  VectorDouble psi, h1, h2, zs, ys;

  int nech = static_cast<int> (tab.size());
  if (nech <= 0) return 1;
  zs.resize(nech+2);
  ys.resize(nech+2);
  for (ih=0; ih<_nbPoly; ih++) _psiHn[ih] = 0.;

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

  h1 = hermitePolynomials(ys[0],1.,_nbPoly);
  Gy1 = 0.;

  for (icl=0 ; icl<ncl ; icl++)
  {
    h2 = hermitePolynomials(ys[icl],1.,_nbPoly);

    Gy2 = law_df_gaussian(ys[icl]);

    for( ih=1 ; ih<_nbPoly ; ih++ )
      _psiHn[ih] += zs[icl] * (h2[ih-1]*Gy2 - h1[ih-1]*Gy1) /sqrt((double) ih);

    Gy1 = Gy2;
    for( ih=0 ; ih<_nbPoly ; ih++) h1[ih] = h2[ih];
  }

  /* Ultimate calculations */

  calculateMeanAndVariance();
  _defineBounds(ys[0], zs[0], ys[ncl - 2] + EPSILON5, zs[ncl - 1], _ay.getVmin(),
                _az.getVmin(), _ay.getVmax(), _az.getVmax());

  return 0;
}

int AnamHermite::fit(Db *db, const ELoc& locatorType)
{
  int number = db->getLocatorNumber(locatorType);
  if (number != 1)
  {
    messerr("The number of items for locator(%d) is %d. It should be 1",
            locatorType.getValue(),number);
    return 1;
  }
  VectorDouble tab = db->getColumnByLocator(locatorType,0,true);
  VectorDouble wt;
  if (db->hasWeight())
    wt = db->getColumnByLocator(ELoc::W,0,true);

  return fit(tab, wt);
}

int AnamHermite::fit(Db *db, const String& name)
{
  VectorDouble tab = db->getColumn(name,true);
  VectorDouble wt;
  if (db->hasWeight())
    wt = db->getColumnByLocator(ELoc::W,0,true);

  return fit(tab, wt);
}

double AnamHermite::getPsiHn(int i) const
{
  if (! _isIndexValid(i)) return TEST;
  return _psiHn[i];
}

void AnamHermite::setPsiHn(VectorDouble psi_hn)
{
  _psiHn = psi_hn;
  _nbPoly = static_cast<int> (_psiHn.size());
}

void AnamHermite::setPsiHn(int i, double psi_hn)
{
  if (! _isIndexValid(i)) return;
  _psiHn[i] = psi_hn;
}

bool AnamHermite::_isIndexValid(int i) const
{
  if (i < 0 || i >= _nbPoly)
  {
    mesArg("Hermite Polynomial Index",i,_nbPoly);
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
  zm[ind0] = GaussianToRawValue(ym[ind0]);
  for (ind=ind0-1; ind>=0; ind--)
  {
    ym[ind] = ym[ind+1] - YPAS;
    zm[ind] = GaussianToRawValue(ym[ind]);
  }
  for (ind=ind0+1; ind<npas; ind++)
  {
    ym[ind] = ym[ind-1] + YPAS;
    zm[ind] = GaussianToRawValue(ym[ind]);
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
      _ay.setVmin(RawToGaussianValue(_az.getVmin()));
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
      _ay.setVmax(RawToGaussianValue(_az.getVmax()));
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
    _az.setVmin(GaussianToRawValue(_ay.getVmin()));
  }
  if (FFFF(_pz.getVmin()))
  {
    _py.setVmin(MAX(ANAM_YMIN,_ay.getVmin()));
    _pz.setVmin(GaussianToRawValue(_py.getVmin()));
  }
  if (FFFF(_az.getVmax()))
  {
    _ay.setVmax(ANAM_YMAX);
    _az.setVmax(GaussianToRawValue(_ay.getVmax()));
  }
  if (FFFF(_pz.getVmax()))
  {
    _py.setVmax(MIN(ANAM_YMAX,_ay.getVmax()));
    _pz.setVmax(GaussianToRawValue(_py.getVmax()));
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

int AnamHermite::_serialize(std::ostream& os, bool verbose) const
{
  if (AnamContinuous::_serialize(os, verbose)) return 1;

  bool ret = _recordWrite<int>(os, "Number of Hermite Polynomials", getNbPoly());
  ret = ret && _recordWrite<double>(os,"Change of support coefficient", getRCoef());
  ret = ret && _tableWrite(os, "Hermite Polynomial", getNbPoly(), getPsiHn());

  return ret ? 0 : 1;
}

int AnamHermite::_deserialize(std::istream& is, bool verbose)
{
  VectorDouble hermite;
  double r = TEST;
  int nbpoly = 0;

  if (! AnamContinuous::_deserialize(is, verbose)) return 1;

  bool ret = _recordRead<int>(is, "Number of Hermite Polynomials", nbpoly);
  ret = ret && _recordRead<double>(is, "Change of Support Coefficient", r);
  if (! ret) return 1;

  hermite.resize(nbpoly);
  if (_tableRead(is, nbpoly, hermite.data())) return 1;

  setNbPoly(nbpoly);
  setRCoef(r);
  setPsiHn(hermite);

  return 0;
}

double AnamHermite::modifyCov(const ECalcMember& member,
                              int iclass,
                              double dist,
                              double /*cov0*/,
                              double cov1,
                              double /*cov2*/) const
{
  double cov;
  double coeff = 0.;
  double rn = pow(_rCoef, (double) iclass);
  if (dist <= 0.)
  {
    switch (member.toEnum())
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
    double rhon = pow(cov1, (double) iclass);
    switch (member.toEnum())
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
  return cov;
}
