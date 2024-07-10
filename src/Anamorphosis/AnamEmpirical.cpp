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
#include "Anamorphosis/AnamEmpirical.hpp"
#include "Anamorphosis/AnamContinuous.hpp"
#include "Matrix/Table.hpp"
#include "Basic/Law.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorHelper.hpp"

#include <math.h>

#define ANAM_YMIN -10.
#define ANAM_YMAX  10.

AnamEmpirical::AnamEmpirical(int ndisc, double sigma2e, bool flagDilution, bool flagGaussian)
    : AnamContinuous(),
      _flagDilution(flagDilution),
      _flagGaussian(flagGaussian),
      _nDisc(ndisc),
      _sigma2e(sigma2e),
      _ZDisc(),
      _YDisc()
{
  _ZDisc.resize(ndisc);
  _YDisc.resize(ndisc);
}

AnamEmpirical::AnamEmpirical(const AnamEmpirical &m)
    : AnamContinuous(m),
      _flagDilution(m._flagDilution),
      _flagGaussian(m._flagGaussian),
      _nDisc(m._nDisc),
      _sigma2e(m._sigma2e),
      _ZDisc(m._ZDisc),
      _YDisc(m._YDisc)
{
}

AnamEmpirical& AnamEmpirical::operator=(const AnamEmpirical &m)
{
  if (this != &m)
  {
    AnamContinuous::operator=(m);
    _flagDilution = m._flagDilution;
    _flagGaussian = m._flagGaussian;
    _nDisc = m._nDisc;
    _sigma2e = m._sigma2e;
    _ZDisc = m._ZDisc;
    _YDisc = m._YDisc;
  }
  return *this;
}

AnamEmpirical::~AnamEmpirical()
{

}

String AnamEmpirical::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << toTitle(1,"Empirical Anamorphosis");

  if (_flagDilution)
  {
    if (_flagGaussian)
      sstr << "Using Gaussian Dilution" << std::endl;
    else
      sstr << "Using Lognormal Dilution" << std::endl;

    sstr << "Number of discretization lags = " << _nDisc << std::endl;
    sstr << "Additional variance           = " << _sigma2e << std::endl;
  }

  if (! _isFitted()) return sstr.str();

  Table ZY(_nDisc, 2, false, true);
  ZY.setTitle("Discretization Intervals");
  ZY.setColumnName(0, "Z");
  ZY.setColumnName(1, "Y");
  ZY.setColumn(0, _ZDisc);
  ZY.setColumn(1, _YDisc);
  sstr << ZY.toString() << std::endl;
  return sstr.str();
}

AnamEmpirical* AnamEmpirical::createFromNF(const String& neutralFilename, bool verbose)
{
  AnamEmpirical* anam = nullptr;
  std::ifstream is;
  anam = new AnamEmpirical();
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

void AnamEmpirical::reset(int ndisc,
                          double pymin,
                          double pzmin,
                          double pymax,
                          double pzmax,
                          double aymin,
                          double azmin,
                          double aymax,
                          double azmax,
                          double sigma2e,
                          const VectorDouble& zdisc,
                          const VectorDouble& ydisc)
{
  setNDisc(ndisc);
  setSigma2e(sigma2e);
  setDisc(zdisc, ydisc);
  setABounds(azmin, azmax, aymin, aymax);
  setPBounds(pzmin, pzmax, pymin, pymax);
}

AnamEmpirical* AnamEmpirical::create(int ndisc, double sigma2e)
{
  return new AnamEmpirical(ndisc, sigma2e);
}

void AnamEmpirical::setNDisc(int ndisc)
{
  _nDisc = ndisc;
  _ZDisc.resize(ndisc);
  _YDisc.resize(ndisc);
}

void AnamEmpirical::setDisc(const VectorDouble& zdisc,
                            const VectorDouble& ydisc)
{
  if ((int) zdisc.size() != (int) ydisc.size())
  {
    messerr("Argumznts 'zdisc' and 'ydisc' should have the same dimension");
    return;
  }
  _ZDisc = zdisc;
  _YDisc = ydisc;
  _nDisc = static_cast<int> (zdisc.size());
}

double AnamEmpirical::rawToTransformValue(double zz) const
{
  double za,zb,ya,yb;
  int    idisc,found;

  /* Initialization */

  if (zz < _ZDisc[0])        zz = _ZDisc[0];
  if (zz > _ZDisc[_nDisc-1]) zz = _ZDisc[_nDisc-1];
  ya = yb = za = zb = zz;

  for (idisc = found = 0; idisc < _nDisc && found == 0; idisc++)
  {
    if (zz > _ZDisc[idisc]) continue;
    yb = _YDisc[idisc];
    zb = _ZDisc[idisc];
    found = 1;
  }

  for (idisc = _nDisc - 1, found = 0; idisc >= 0 && found == 0; idisc--)
  {
    if (zz < _ZDisc[idisc]) continue;
    ya = _YDisc[idisc];
    za = _ZDisc[idisc];
    found = 1;
  }

  double yy;
  if (za >= zb)
    yy = ya;
  else
    yy = ((zb - zz) * ya + (zz - za) * yb) / (zb - za);

  return(yy);
}

double AnamEmpirical::transformToRawValue(double yy) const
{
  double zz,za,zb,ya,yb;
  int    idisc,found;

  /* Initialization */

  if (yy < _YDisc[0])        yy = _YDisc[0];
  if (yy > _YDisc[_nDisc-1]) yy = _YDisc[_nDisc-1];
  zz = za = zb = ya = yb = yy;

  for (idisc=found=0; idisc<_nDisc && found==0; idisc++)
  {
    if (yy > _YDisc[idisc]) continue;
    yb = _YDisc[idisc];
    zb = _ZDisc[idisc];
    found = 1;
  }

  for (idisc=_nDisc-1, found=0; idisc>=0 && found==0; idisc--)
  {
    if (yy < _YDisc[idisc]) continue;
    ya = _YDisc[idisc];
    za = _ZDisc[idisc];
    found = 1;
  }

  if (ya >= yb)
    zz = za;
  else
    zz = ((yb - yy) * za + (yy - ya) * zb) / (yb - ya);

  return(zz);
}

void AnamEmpirical::calculateMeanAndVariance()
{
  messerr("This function is not available for Empirical Anamorphosis");
}

int AnamEmpirical::_getStatistics(const VectorDouble &tab,
                                  int *count,
                                  double *mean,
                                  double *mean2,
                                  double *mini,
                                  double *maxi,
                                  double *var)
{
  int nech = static_cast<int> (tab.size());
  int number = 0;
  double dmean  = 0.;
  double dmean2 = 0.;
  double dmaxi  = -1.e30;
  double dmini  =  1.e30;
  for (int iech = 0; iech < nech; iech++)
  {
    double value = tab[iech];
    if (FFFF(value)) continue;
    number ++;
    dmean  += value;
    dmean2 += value * value;
    if (value > dmaxi) dmaxi = value;
    if (value < dmini) dmini = value;
  }
  if (number <= 0)
  {
    messerr("The number of strictly positive data is zero");
    return 1;
  }
  dmean /= number;
  dmean2 = dmean2 / number;
  double variance = dmean2 - dmean * dmean;

  // Returning arguments

  *count = number;
  *mean = dmean;
  *mean2 = dmean2;
  *mini = dmini;
  *maxi = dmaxi;
  *var = variance;
  return 0;
}

int AnamEmpirical::_fitWithDilutionGaussian(const VectorDouble &tab)
{
  int number;
  double dmean, dmean2, dmini, dmaxi, variance;

  /* Calculate the constants */

  if (_getStatistics(tab, &number, &dmean, &dmean2, &dmini, &dmaxi, &variance)) return 1;
  int nech = static_cast<int> (tab.size());

  /* Calculation parameters */

  if (FFFF(_sigma2e)) _sigma2e = variance / (2. * number);
  double sigma = sqrt(_sigma2e);
  double ecart  = dmaxi - dmini;
  double disc_val  = 3 * ecart / (_nDisc - 2);
  double disc_init = dmini - MIN(disc_val / 2., dmini / 10000.);

  /* Fill the discretized array */

  _ZDisc[0] = disc_init - disc_val;
  _ZDisc[1] = disc_init;
  for (int idisc = 2; idisc < _nDisc; idisc++)
    _ZDisc[idisc] = disc_init + (idisc - 1) * disc_val;

  for (int idisc = 0; idisc < _nDisc; idisc++)
  {
    double zval  = _ZDisc[idisc];
    double total = 0.;
    for (int iech = 0; iech < nech; iech++)
    {
      double value = tab[iech];
      if (FFFF(value) || value <= 0) continue;
      total += law_cdf_gaussian((zval - value) / sigma);
    }
    _YDisc[idisc] = total / number;
  }

  /* Re-allocate memory to the only necessary bits */

  int ndisc_util = 0;
  for (int idisc = 0; idisc < _nDisc; idisc++)
  {
    if (_YDisc[idisc] >= 1.) break;
    ndisc_util++;
  }
  _ZDisc.resize(ndisc_util);
  _YDisc.resize(ndisc_util);
  _nDisc = ndisc_util;

  for (int idisc=0; idisc<_nDisc; idisc++)
    _YDisc[idisc] = law_invcdf_gaussian(_YDisc[idisc]);

  return 0;
}

int AnamEmpirical::_fitWithDilutionLognormal(const VectorDouble &tab)
{
  int number;
  double dmean, dmean2, dmini, dmaxi, variance;

  /* Calculate the constants */

  if (_getStatistics(tab, &number, &dmean, &dmean2, &dmini, &dmaxi, &variance)) return 1;
  int nech = static_cast<int> (tab.size());
  if (dmini < 0.)
  {
    messerr("The Anamorphosis by Lognormal Dilution is not compatible with negative data values");
    return 1;
  }

  /* Calculation parameters */

  if (FFFF(_sigma2e)) _sigma2e = variance / (2. * number);
  double sigma  = sqrt(log(1. + _sigma2e / dmean2));
  double ecart  = dmaxi - dmini;
  double disc_val  = 3 * ecart / (_nDisc - 2);
  double disc_init = dmini - MIN(disc_val / 2., dmini / 10000.);

  /* Fill the discretized array */

  _ZDisc[0] = disc_init - disc_val;
  _ZDisc[1] = disc_init;
  for (int idisc = 2; idisc < _nDisc; idisc++)
    _ZDisc[idisc] = disc_init + (idisc - 1) * disc_val;

  _YDisc[0] = 1.;
  for (int idisc = 1; idisc < _nDisc; idisc++)
  {
    double zval  = _ZDisc[idisc];
    double total = 0.;
    for (int iech = 0; iech < nech; iech++)
    {
      double value = tab[iech];
      if (FFFF(value) || value <= 0) continue;
      total += 1. - law_cdf_gaussian(sigma / 2. + log(zval / value) / sigma);
    }
    _YDisc[idisc] = total / number;
  }

  /* Re-allocate memory to the only necessary bits */

  int ndisc_util = 0;
  for (int idisc = 0; idisc < _nDisc; idisc++)
    if (_YDisc[idisc] > 0) ndisc_util++;
  _ZDisc.resize(ndisc_util);
  _YDisc.resize(ndisc_util);
  _nDisc = ndisc_util;

  for (int idisc=0; idisc<_nDisc; idisc++)
    _YDisc[idisc] = law_invcdf_gaussian(1. - _YDisc[idisc]);
  return 0;
}

int AnamEmpirical::_fitNormalScore(const VectorDouble &tab)
{
  _ZDisc.clear();
  _YDisc.clear();

  for (int iech = 0, nech = (int) tab.size(); iech < nech; iech++)
  {
    if (FFFF(tab[iech])) continue;
    _ZDisc.push_back(tab[iech]);
  }

  // Sort the Z-values by ascending order
  VH::sortInPlace(_ZDisc);
  _YDisc = VH::normalScore(_ZDisc);
  _nDisc = (int) _ZDisc.size();

  return 0;
}

int AnamEmpirical::fitFromArray(const VectorDouble& tab,
                                const VectorDouble& wt)
{
  DECLARE_UNUSED(wt);

  if (_flagDilution)
  {
    if (_flagGaussian)
    {
      if (_fitWithDilutionGaussian(tab)) return 1;
    }
    else
    {
      if (_fitWithDilutionLognormal(tab)) return 1;
    }
  }
  else
  {
    if (_fitNormalScore(tab)) return 1;
  }

  for (int idisc=0; idisc<_nDisc; idisc++)
  {
    if (_YDisc[idisc] < ANAM_YMIN) _YDisc[idisc] = ANAM_YMIN;
    if (_YDisc[idisc] > ANAM_YMAX) _YDisc[idisc] = ANAM_YMAX;
  }

  /* Evaluate the bounds */

  double pzmin =  1.e30;
  double pymin =  1.e30;
  double pzmax = -1.e30;
  double pymax = -1.e30;
  for (int idisc=0; idisc<_nDisc; idisc++)
  {
    if (_ZDisc[idisc] < pzmin) pzmin = _ZDisc[idisc];
    if (_ZDisc[idisc] > pzmax) pzmax = _ZDisc[idisc];
    if (_YDisc[idisc] < pymin) pymin = _YDisc[idisc];
    if (_YDisc[idisc] > pymax) pymax = _YDisc[idisc];
  }

  /* Save the results */

  setABounds(pzmin, pzmax, pymin, pymax);
  setPBounds(pzmin, pzmax, pymin, pymax);

  return 0;
}

bool AnamEmpirical::_serialize(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret = ret && AnamContinuous::_serialize(os, verbose);
  ret = ret && _recordWrite<int>(os, "Number of Discretization lags", getNDisc());
  ret = ret && _recordWrite<double>(os, "additional variance", getSigma2e());
  ret = ret && _tableWrite(os, "Z Values", getNDisc(), getZDisc());
  ret = ret && _tableWrite(os, "Y Values", getNDisc(), getYDisc());
  return ret;
}

bool AnamEmpirical::_deserialize(std::istream& is, bool verbose)
{
  int ndisc = 0;
  double sigma2e = TEST;
  VectorDouble zdisc, ydisc;

  bool ret = true;
  ret = ret && AnamContinuous::_deserialize(is, verbose);
  ret = ret && _recordRead<int>(is, "Number of Discretization classes", ndisc);
  ret = ret && _recordRead<double>(is, "Experimental Error Variance", sigma2e);

  if (ret)
  {
    zdisc.resize(ndisc);
    ret = ret && _tableRead(is, "Z Values", ndisc, zdisc.data());
    ydisc.resize(ndisc);
    ret = ret && _tableRead(is, "Y Values", ndisc, ydisc.data());
  }

  if (ret)
  {
    setNDisc(ndisc);
    setSigma2e(sigma2e);
    setDisc(zdisc, ydisc);
  }
  return ret;
}
