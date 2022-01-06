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
#include "Matrix/MatrixRectangular.hpp"
#include "geoslib_f.h"

#define ANAM_KD_NELEM 6

AnamDiscrete::AnamDiscrete(const EAnam& type)
    : Anam(type),
      _nCut(0),
      _nElem(ANAM_KD_NELEM),
      _mean(TEST),
      _variance(TEST),
      _zCut(),
      _stats()
{
}

AnamDiscrete::AnamDiscrete(const AnamDiscrete &m)
    : Anam(m),
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
    _nCut = m._nCut;
    _nElem = m._nElem;;
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
  int ncut = getNCut();
  int nclass = getNClass();
  int nelem = getNElem();

  _zCut.resize(ncut,0.);
  _stats.reset(nclass,nelem,0.);
}

String AnamDiscrete::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  sstr << "Number of cutoffs = " << _nCut << std::endl;
  sstr << "Number of classes = " << getNClass() << std::endl;
  sstr << "Mean              = " << _mean << std::endl;
  sstr << "Variance          = " << _variance << std::endl;
  sstr << std::endl;
  sstr << toMatrix("Cutoffs", VectorString(), VectorString(), true, 1, _nCut, _zCut);
  return sstr.str();
}

void AnamDiscrete::setNCut(int ncut)
{
  _nCut = ncut;
  _resize();
}

void AnamDiscrete::setZCut(const VectorDouble& zcut)
{
  _zCut = zcut;
  _nCut = static_cast<int> (_zCut.size());
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
