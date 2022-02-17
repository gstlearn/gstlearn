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
#include "Anamorphosis/AnamEmpirical.hpp"
#include "Anamorphosis/AnamContinuous.hpp"
#include "Basic/Law.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Basic/AStringable.hpp"
#include "geoslib_f.h"
#include <math.h>

#define ANAM_YMIN -10.
#define ANAM_YMAX  10.

#define YD(i)           (_tDisc[(i)])
#define ZD(i)           (_tDisc[(i) + _nDisc])
#define ZDC(i)          (_tDisc[(i) + ndisc_util])

AnamEmpirical::AnamEmpirical(int ndisc, double sigma2e)
    : AnamContinuous(),
      _nDisc(ndisc),
      _sigma2e(sigma2e),
      _tDisc()
{
  _tDisc.resize(2 * ndisc);
}

AnamEmpirical::AnamEmpirical(const AnamEmpirical &m)
    : AnamContinuous(m),
      _nDisc(m._nDisc),
      _sigma2e(m._sigma2e),
      _tDisc(m._tDisc)
{
}

AnamEmpirical& AnamEmpirical::operator=(const AnamEmpirical &m)
{
  if (this != &m)
  {
    _nDisc = m._nDisc;
    _sigma2e = m._sigma2e;
    _tDisc = m._tDisc;
  }
  return *this;
}

AnamEmpirical::~AnamEmpirical()
{

}

String AnamEmpirical::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  sstr << toTitle(1,"Empirical Anamorphosis");

  sstr << AAnam::toString(strfmt);

  sstr << "Number of discretization lags = " << _nDisc << std::endl;
  sstr << "Additional variance           = " << _sigma2e << std::endl;
  sstr << std::endl;
  sstr << "Discretization intervals Y - Z" << std::endl;
  sstr << toMatrix(String(), VectorString(), VectorString(), true, 2, _nDisc, _tDisc);

  return sstr.str();
}

int AnamEmpirical::dumpToNF(const String& neutralFilename, bool verbose) const
{
  FILE* file = _fileOpen(neutralFilename, "AnamEmpirical", "w", verbose);
  if (file == nullptr) return 1;

  if (_serialize(file, verbose))
  {
    if (verbose) messerr("Problem writing in the Neutral File.");
    _fileClose(file, verbose);
    return 1;
  }
  _fileClose(file, verbose);
  return 0;
}

int AnamEmpirical::dumpToNF2(const String& neutralFilename, bool verbose) const
{
  std::ofstream os;
  int ret = 1;
  if (_fileOpenWrite2(neutralFilename, "AnamEmpirical", os, verbose))
  {
    ret = _serialize2(os, verbose);
    if (ret && verbose) messerr("Problem writing in the Neutral File.");
    os.close();
  }
  return ret;
}

AnamEmpirical* AnamEmpirical::createFromNF(const String& neutralFilename, bool verbose)
{
  FILE* file = _fileOpen(neutralFilename, "AnamEmpirical", "r", verbose);
  if (file == nullptr) return nullptr;

  AnamEmpirical* anam = new AnamEmpirical();
  if (anam->_deserialize(file, verbose))
  {
    if (verbose) messerr("Problem reading the Neutral File");
    delete anam;
    anam = nullptr;
  }
  _fileClose(file, verbose);
  return anam;
}

AnamEmpirical* AnamEmpirical::createFromNF2(const String& neutralFilename, bool verbose)
{
  AnamEmpirical* anam = nullptr;
  std::ifstream is;
  if (_fileOpenRead2(neutralFilename, "AnamEmpirical", is, verbose))
  {
    anam = new AnamEmpirical();
    if (anam->_deserialize2(is, verbose))
    {
      if (verbose) messerr("Problem reading the Neutral File");
      delete anam;
      anam = nullptr;
    }
    is.close();
  }
  return anam;
}


AnamEmpirical* AnamEmpirical::create(int ndisc, double sigma2e)
{
  return new AnamEmpirical(ndisc, sigma2e);
}

void AnamEmpirical::setNDisc(int ndisc)
{
  _nDisc = ndisc;
  _tDisc.resize(2 * ndisc);
}

void AnamEmpirical::setTDisc(const VectorDouble& tdisc)
{
  _tDisc = tdisc;
  _nDisc = static_cast<int> (tdisc.size()) / 2;
}

double AnamEmpirical::RawToGaussianValue(double zz) const
{
  double yy,za,zb,ya,yb;
  int    idisc,found;

  /* Initialization */

  if (zz < ZD(0))        zz = ZD(0);
  if (zz > ZD(_nDisc-1)) zz = ZD(_nDisc-1);
  yy = ya = yb = za = zb = zz;

  for (idisc=found=0; idisc<_nDisc && found==0; idisc++)
  {
    if (zz > ZD(idisc)) continue;
    yb = YD(idisc);
    zb = ZD(idisc);
    found = 1;
  }

  for (idisc=_nDisc-1, found=0; idisc>=0 && found==0; idisc--)
  {
    if (zz < ZD(idisc)) continue;
    ya = YD(idisc);
    za = ZD(idisc);
    found = 1;
  }

  if (za >= zb)
    yy = ya;
  else
    yy = ((zb - zz) * ya + (zz - za) * yb) / (zb - za);

  return(yy);
}

double AnamEmpirical::GaussianToRawValue(double yy) const
{
  double zz,za,zb,ya,yb;
  int    idisc,found;

  /* Initialization */

  if (yy < YD(0))        yy = YD(0);
  if (yy > YD(_nDisc-1)) yy = YD(_nDisc-1);
  zz = za = zb = ya = yb = yy;

  for (idisc=found=0; idisc<_nDisc && found==0; idisc++)
  {
    if (yy > YD(idisc)) continue;
    yb = YD(idisc);
    zb = ZD(idisc);
    found = 1;
  }

  for (idisc=_nDisc-1, found=0; idisc>=0 && found==0; idisc--)
  {
    if (yy < YD(idisc)) continue;
    ya = YD(idisc);
    za = ZD(idisc);
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
  my_throw("This function is not available for Empirical Anamorphosis");
}

int AnamEmpirical::fit(const VectorDouble& tab)
{
  int     iech,number,idisc,id,ndisc_util;
  double  value,dmean,dmean2,dmini,dmaxi,sigma,zval,total,variance;
  double  disc_val,disc_init,pzmin,pzmax,pymin,pymax,ecart;

  /* Calculate the constants */

  int nech = static_cast<int> (tab.size());
  number = 0;
  dmean  = dmean2 = 0.;
  dmaxi  = -1.e30;
  dmini  =  1.e30;
  for (iech=0; iech<nech; iech++)
  {
    value = tab[iech];
    if (FFFF(value)) continue;
    if (value < 0.)
    {
      messerr("The Anamorphosis by lognormal dilution");
      messerr("is not compatible with negative data values");
      return 1;
    }
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
  variance = dmean2 - dmean * dmean;
  if (FFFF(_sigma2e)) _sigma2e = variance / (2. * number);
  sigma  = sqrt(log(1. + _sigma2e / dmean2));
  ecart  = dmaxi - dmini;

  /* Calculation parameters */

  disc_val  = 3 * ecart / (_nDisc - 2);
  disc_init = dmini - MIN(disc_val / 2., dmini / 10000.);

  /* Fill the discretized array */

  idisc = id = 0;
  ZD(0) = disc_init - disc_val;
  ZD(1) = disc_init;
  for (idisc=2; idisc<_nDisc; idisc++)
    ZD(idisc) = disc_init + (idisc - 1) * disc_val;

  YD(0) = 1.;
  for (idisc=1; idisc<_nDisc; idisc++)
  {
    zval  = ZD(idisc);
    total = 0.;
    for (iech=0; iech<nech; iech++)
    {
      value = tab[iech];
      if (FFFF(value) || value <= 0) continue;
      total += 1. - law_cdf_gaussian(sigma / 2. + log(zval / value) / sigma);
    }
    YD(idisc) = total / number;
  }

  /* Re-allocate memory to the only necessary bits */

  ndisc_util = 0;
  for (idisc=0; idisc<_nDisc; idisc++)
    if (YD(idisc) > 0 && idisc > ndisc_util) ndisc_util = idisc;
  for (idisc=0; idisc<ndisc_util; idisc++) ZDC(idisc) = ZD(idisc);
  _tDisc.resize(2 * ndisc_util);
  _nDisc = ndisc_util;

  for (idisc=0; idisc<_nDisc; idisc++)
  {
    YD(idisc) = law_invcdf_gaussian(1. - YD(idisc));
    if (YD(idisc) < ANAM_YMIN) YD(idisc) = ANAM_YMIN;
    if (YD(idisc) > ANAM_YMAX) YD(idisc) = ANAM_YMAX;
  }

  /* Evaluate the bounds */

  pzmin = pymin =  1.e30;
  pzmax = pymax = -1.e30;
  for (idisc=0; idisc<_nDisc; idisc++)
  {
    if (ZD(idisc) < pzmin) pzmin = ZD(idisc);
    if (ZD(idisc) > pzmax) pzmax = ZD(idisc);
    if (YD(idisc) < pymin) pymin = YD(idisc);
    if (YD(idisc) > pymax) pymax = YD(idisc);
  }

  /* Save the results */

  setABounds(_az.getVmin(), _az.getVmax(), _ay.getVmin(), _ay.getVmax());
  setPBounds(pzmin, pzmax, pymin, pymax);

  return 0;
}

bool AnamEmpirical::isTDiscIndexValid(int i) const
{
  if (i < 0 || i >= 2 * _nDisc)
  {
    mesArg("TDisc Index",i,2 * _nDisc);
    return false;
  }
  return true;
}

int AnamEmpirical::_serialize(FILE* file, bool verbose) const
{
  AnamContinuous::_serialize(file, verbose);

  _recordWrite(file, "%d", getNDisc());
  _recordWrite(file, "#", "Number of Discretization lags");
  _recordWrite(file, "%ld", getSigma2e());
  _recordWrite(file, "#", "additional variance");
  _tableWrite(file, "Coefficients", 2 * getNDisc(), getTDisc().data());

  return 0;
}

int AnamEmpirical::_serialize2(std::ostream& os, bool verbose) const
{
  if (AnamContinuous::_serialize2(os, verbose)) return 1;

  bool ret = _recordWrite2<int>(os, "Number of Discretization lags", getNDisc());
  ret = ret && _recordWrite2<double>(os, "additional variance", getSigma2e());
  ret = ret && _tableWrite2(os, "Coefficients", 2 * getNDisc(), getTDisc());

  return ret ? 0 : 1;
}

int AnamEmpirical::_deserialize(FILE* file, bool verbose)
{
  int ndisc = 0;
  double sigma2e = TEST;
  VectorDouble tdisc;

  if (AnamContinuous::_deserialize(file, verbose)) goto label_end;

  if (_recordRead(file, "Number of Discretization classes", "%d", &ndisc))
    goto label_end;
  if (_recordRead(file, "Experimental Error Variance", "%lf", &sigma2e))
    goto label_end;
  tdisc.resize(2 * ndisc);
  if (_tableRead(file, 2 * ndisc, tdisc.data())) goto label_end;

  setNDisc(ndisc);
  setSigma2e(sigma2e);
  setTDisc(tdisc);

  label_end:
  return 0;
}

int AnamEmpirical::_deserialize2(std::istream& is, bool verbose)
{
  int ndisc = 0;
  double sigma2e = TEST;
  VectorDouble tdisc;

  if (! AnamContinuous::_deserialize2(is, verbose)) return 1;

  bool ret = _recordRead2<int>(is, "Number of Discretization classes", ndisc);
  ret = ret && _recordRead2<double>(is, "Experimental Error Variance", sigma2e);
  if (! ret) return 1;

  tdisc.resize(2 * ndisc);
  if (_tableRead2(is, 2 * ndisc, tdisc.data())) return 1;

  setNDisc(ndisc);
  setSigma2e(sigma2e);
  setTDisc(tdisc);
  return 0;
}
