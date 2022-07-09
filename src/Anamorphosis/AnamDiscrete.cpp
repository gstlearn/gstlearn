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
#include "Anamorphosis/AnamDiscrete.hpp"
#include "Stats/Selectivity.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/AMatrix.hpp"
#include "geoslib_f.h"

#include <math.h>

#define ANAM_KD_NELEM 6

AnamDiscrete::AnamDiscrete()
    : AAnam(),
      _nCut(0),
      _nElem(ANAM_KD_NELEM),
      _mean(TEST),
      _variance(TEST),
      _zCut(),
      _stats()
{
  _resize();
}

AnamDiscrete::AnamDiscrete(const AnamDiscrete &m)
    : AAnam(m),
      _nCut(m._nCut),
      _nElem(m._nElem),
      _mean(m._mean),
      _variance(m._variance),
      _zCut(m._zCut),
      _stats(m._stats)
{
}

AnamDiscrete& AnamDiscrete::operator=(const AnamDiscrete &m)
{
  if (this != &m)
  {
    AAnam::operator= (m);
    _nCut = m._nCut;
    _nElem = m._nElem;
    _zCut = m._zCut;
    _mean = m._mean;
    _variance = m._variance;
    _stats = m._stats;
  }
  return *this;
}

AnamDiscrete::~AnamDiscrete()
{

}

void AnamDiscrete::_resize()
{
  int ncut   = getNCut();
  int nclass = getNClass();
  int nelem  = getNElem();

  _zCut.resize(ncut,0.);
  _stats.reset(nclass,nelem,0.);
}

String AnamDiscrete::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  sstr << "Number of cutoffs = " << _nCut << std::endl;
  sstr << "Number of classes = " << getNClass() << std::endl;
  if (! FFFF(_mean))
    sstr << "Mean              = " << _mean << std::endl;
  if (! FFFF(_variance))
    sstr << "Variance          = " << _variance << std::endl;

  sstr << std::endl;
  sstr << toMatrix("Cutoffs", VectorString(), VectorString(), true, 1, _nCut, _zCut);

  sstr << toMatrix(String(),VectorString(),VectorString(),true,getNElem(),getNClass(),
                   getStats().getValues());

  return sstr.str();
}

void AnamDiscrete::calculateMeanAndVariance()
{
  _mean = TEST;
  _variance = TEST;
}

double AnamDiscrete::getDDStatProp  (int iclass) const
{
  if (! _isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass,0);
}
double AnamDiscrete::getDDStatZmoy  (int iclass) const
{
  if (! _isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass,1);
}
double AnamDiscrete::getDDStatCnorm (int iclass) const
{
  if (! _isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass,2);
}
double AnamDiscrete::getDDStatLambda(int iclass) const
{
  if (! _isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass,3);
}
double AnamDiscrete::getDDStatU     (int iclass) const
{
  if (! _isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass,4);
}
double AnamDiscrete::getDDStatMul   (int iclass) const
{
  if (! _isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass,5);
}
void   AnamDiscrete::setDDStatProp  (int iclass, double value)
{
  if (! _isClassValid(iclass)) return;
  _stats.setValue(iclass,0,value);
}
void   AnamDiscrete::setDDStatZmoy  (int iclass, double value)
{
  if (! _isClassValid(iclass)) return;
  _stats.setValue(iclass,1,value);
}
void   AnamDiscrete::setDDStatCnorm (int iclass, double value)
{
  if (! _isClassValid(iclass)) return;
  _stats.setValue(iclass,2,value);
}
void   AnamDiscrete::setDDStatLambda(int iclass, double value)
{
  if (! _isClassValid(iclass)) return;
  _stats.setValue(iclass,3,value);
}
void   AnamDiscrete::setDDStatU     (int iclass, double value)
{
  if (! _isClassValid(iclass)) return;
  _stats.setValue(iclass,4,value);
}
void   AnamDiscrete::setDDStatMul   (int iclass, double value)
{
  if (! _isClassValid(iclass)) return;
  _stats.setValue(iclass,5,value);
}

// Function for using Stats in IR anamorphosis
double AnamDiscrete::getIRStatT   (int iclass) const
{
  if (! _isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass,0);
}
double AnamDiscrete::getIRStatQ   (int iclass) const
{
  if (! _isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass,1);
}
double AnamDiscrete::getIRStatZ   (int iclass) const
{
  if (! _isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass,2);
}
double AnamDiscrete::getIRStatB   (int iclass) const
{
  if (! _isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass,3);
}
double AnamDiscrete::getIRStatR   (int iclass) const
{
  if (! _isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass,4);
}
double AnamDiscrete::getIRStatRV  (int iclass) const
{
  if (! _isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass,5);
}
void AnamDiscrete::setIRStatT(int iclass, double value)
{
  if (! _isClassValid(iclass)) return;
  _stats.setValue(iclass,0,value);
}
void AnamDiscrete::setIRStatQ(int iclass, double value)
{
  if (! _isClassValid(iclass)) return;
  _stats.setValue(iclass,1,value);
}
void AnamDiscrete::setIRStatZ(int iclass, double value)
{
  if (! _isClassValid(iclass)) return;
  _stats.setValue(iclass,2,value);
}
void AnamDiscrete::setIRStatB(int iclass, double value)
{
  if (! _isClassValid(iclass)) return;
  _stats.setValue(iclass,3,value);
}
void AnamDiscrete::setIRStatR(int iclass, double value)
{
  if (! _isClassValid(iclass)) return;
  _stats.setValue(iclass,4,value);
}
void AnamDiscrete::setIRStatRV(int iclass, double value)
{
  if (! _isClassValid(iclass)) return;
  _stats.setValue(iclass,5,value);
}

bool AnamDiscrete::_isClassValid(int iclass) const
{
  if (iclass < 0 || iclass >= getNClass())
  {
    mesArg("Class Index",iclass,getNClass());
    return false;
  }
  return true;
}

bool AnamDiscrete::_serialize(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;
  ret = ret && _recordWrite<int>(os, "Number of Cuttofs", getNCut());
  ret = ret && _recordWrite<int>(os, "Number of classes", getNClass());
  ret = ret && _recordWrite<int>(os, "Number of elements", getNElem());
  ret = ret && _tableWrite(os, "Cutoff value", getNCut(), getZCut());
  ret = ret && _tableWrite(os, "DD Stats", getNClass() * getNElem(), getStats().getValues());
  return ret;
}

bool AnamDiscrete::_deserialize(std::istream& is, bool /*verbose*/)
{
  VectorDouble zCut, stats;
  int nCut = 0;
  int nClass = 0;
  int nElem = 0;

  bool ret = true;
  ret = ret && _recordRead<int>(is, "Number of Cutoffs", nCut);
  ret = ret && _recordRead<int>(is, "Number of Classes", nClass);
  ret = ret && _recordRead<int>(is, "Number of Statistic Columns", nElem);

  if (ret)
  {
    zCut.resize(nCut);
    ret = ret && _tableRead(is, nCut, zCut.data());
  }

  if (ret)
  {
    stats.resize(nClass * nElem);
    ret = ret && _tableRead(is, nClass * nElem, stats.data());
  }

  if (ret)
  {
    setNCut(nCut);
    setNElem(nElem);
    setZCut(zCut);
    setStats(stats);
  }
  return ret;
}

void AnamDiscrete::setNCut(int ncut)
{
  _nCut = ncut;
  _resize();
}

void AnamDiscrete::setZCut(const VectorDouble& zcut)
{
  _nCut = (int) zcut.size();
  _resize();
  _zCut = zcut;
};

void AnamDiscrete::setNElem(int nelem)
{
  _nElem = nelem;
  _resize();
}

void AnamDiscrete::setStats(const VectorDouble& stats)
{
  int nclass = getNClass();
  int nelem = getNElem();
  if ((int) stats.size() != nclass * nelem)
  {
    messerr("Argument 'stats' incorrect. Its dimension (%d) should be %d * %d",
            stats.size(),nclass,nelem);
    return;
  }
  _stats.setValues(stats);
}

/*****************************************************************************/
/*!
 **  Interpolate the QT curves (Local estimation)
 **
 ** \param[in]  z_max    Maximum grade value (if defined)
 ** \param[in]  zcutmine Array of the requested cutoffs
 ** \param[in]  calest   Selectivity
 **
 ** \param[out] calcut   Interpolated Selectivity
 **
 *****************************************************************************/
void AnamDiscrete::_interpolateQTLocal(double z_max,
                                       const VectorDouble& zcutmine,
                                       Selectivity& calest,
                                       Selectivity& calcut) const
{
  double tval, qval;

  int nclass = calest.getNCuts();
  int ncutmine = calcut.getNCuts();
  VectorDouble zz(nclass + 2);
  VectorDouble TT(nclass + 2);
  VectorDouble QQ(nclass + 2);

  /* Load arrays */

  int ncleff = 1;
  TT[0] = QQ[0] = 0.;
  for (int iclass = 0; iclass < nclass; iclass++)
  {
    int jclass = nclass - iclass - 1;
    if (calest.getTest(jclass) <= TT[ncleff - 1]) continue;
    TT[ncleff] = calest.getTest(jclass);
    QQ[ncleff] = calest.getQest(jclass);
    ncleff++;
  }
  zz[0] = z_max;
  for (int iclass = 0; iclass < ncleff - 1; iclass++)
    zz[iclass + 1] = (QQ[iclass + 2] - QQ[iclass]) / (TT[iclass + 2] - TT[iclass]);
  zz[ncleff - 1] = 0.;
  if (FFFF(z_max)) zz[0] = 2 * zz[1];

  for (int icut = 0; icut < ncutmine; icut++)
  {
    double zval = zcutmine[icut];
    calcut.setZcut(icut, zval);

    /* Find interval [zz[iclass]; zz[iclass+1]] to which cutoffs belongs */

    int iclass = -1;
    for (int jclass = 0; jclass < ncleff && iclass < 0; jclass++)
      if ((zval - zz[jclass]) * (zval - zz[jclass + 1]) <= 0) iclass = jclass;

    /* Assuming that cutoffs belongs to the interval the class 'iclass' */

    double zi0 = zz[iclass];
    double zi1 = (iclass + 1 > ncleff - 1) ? 0. : zz[iclass + 1];
    double ti0 = TT[iclass];
    double ti1 = (iclass + 1 > ncleff - 1) ? 0. : TT[iclass + 1];
    double qi0 = QQ[iclass];
    double qi1 = QQ[iclass + 1];
    _interpolateInterval(zval, zi0, zi1, ti0, ti1, qi0, qi1, &tval, &qval);
    calcut.setTest(icut, tval);
    calcut.setQest(icut, qval);
  }
  return;
}

/*****************************************************************************/
/*!
 **  Interpolate the QT within an interval
 **
 ** \param[in]  zval     Cutoff value
 ** \param[in]  zi0      Lower cutoff of the interval
 ** \param[in]  zi1      Upper cutoff of the interval
 ** \param[in]  ti0      Lower tonnage of the interval
 ** \param[in]  ti1      Upper tonnage of the interval
 ** \param[in]  qi0      Lower metal quantity of the interval
 ** \param[in]  qi1      Upper metal quantity of the interval
 **
 ** \param[out] tval     Tonnage for the current cutoff
 ** \param[out] qval     Metal quantity for the current cutoff
 **
 *****************************************************************************/
void AnamDiscrete::_interpolateInterval(double zval,
                                        double zi0,
                                        double zi1,
                                        double ti0,
                                        double ti1,
                                        double qi0,
                                        double qi1,
                                        double *tval,
                                        double *qval) const
{
  double dzi, dti, u, aa0, zmoy;
  static double tol = 1.e-3;

  dzi = zi1 - zi0;
  dti = ti1 - ti0;
  zmoy = (qi1 - qi0) / (ti1 - ti0);
  aa0 = (zi1 - zmoy) / (zmoy - zi0);

  if (ABS(zval - zi0) < tol)
  {
    (*tval) = ti0;
    (*qval) = qi0;
    return;
  }

  if (ABS(zval - zi1) < tol)
  {
    (*tval) = ti1;
    (*qval) = qi1;
    return;
  }

  u = (zval - zi0) / dzi;
  (*tval) = (u <= 0.) ? ti0 : ti0 + dti * pow(u, 1. / aa0);
  (*qval) = (u <= 0.) ? qi0 :
      qi0 + zi0 * ((*tval) - ti0)
      + dzi * dti * pow(u, 1. + 1. / aa0) / (1. + aa0);
}

