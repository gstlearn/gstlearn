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
#include "Matrix/AMatrix.hpp"
#include "geoslib_f.h"

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

int AnamDiscrete::_serialize(FILE* file, bool verbose) const
{
  _recordWrite(file, "%d", getNCut());
  _recordWrite(file, "#", "Number of cutoffs");
  _recordWrite(file, "%d", getNClass());
  _recordWrite(file, "#", "Number of classes");
  _recordWrite(file, "%d", getNElem());
  _recordWrite(file, "#", "Number of elements");
  _tableWrite(file, "Cutoff value", getNCut(), getZCut().data());
  _tableWrite(file, "DD Stats", getNClass() * getNElem(), getStats().getValues().data());

  return 0;
}

int AnamDiscrete::_serialize2(std::ostream& os, bool verbose) const
{
  bool ret = _recordWrite2<int>(os, "Number of Cuttofs", getNCut());
  ret = ret && _recordWrite2<int>(os, "Number of classes", getNClass());
  ret = ret && _recordWrite2<int>(os, "Number of elements", getNElem());
  ret = ret && _tableWrite2(os, "Cutoff value", getNCut(), getZCut());
  ret = ret && _tableWrite2(os, "DD Stats", getNClass() * getNElem(), getStats().getValues());

  return ret ? 0 : 1;
}


int AnamDiscrete::_deserialize(FILE* file, bool verbose)
{
  VectorDouble zCut, stats;
  int nCut = 0;
  int nClass = 0;
  int nElem = 0;

  if (_recordRead(file, "Number of Cutoffs", "%d", &nCut))
    goto label_end;
  if (_recordRead(file, "Number of Classes", "%d", &nClass))
    goto label_end;
  if (_recordRead(file, "Number of Statistic Columns", "%d", &nElem))
    goto label_end;

  zCut.resize(nCut);
  if (_tableRead(file, nCut, zCut.data())) goto label_end;

  stats.resize(nClass * nElem);
  if (_tableRead(file, nClass * nElem, stats.data())) goto label_end;

  setNCut(nCut);
  setNElem(nElem);
  setZCut(zCut);
  setStats(stats);

  label_end:
  return 0;
}

int AnamDiscrete::_deserialize2(std::istream& is, bool verbose)
{
  VectorDouble zCut, stats;
  int nCut = 0;
  int nClass = 0;
  int nElem = 0;

  bool ret = _recordRead2<int>(is, "Number of Cutoffs", nCut);
  ret = ret && _recordRead2<int>(is, "Number of Classes", nClass);
  ret = ret && _recordRead2<int>(is, "Number of Statistic Columns", nElem);
  if (! ret) return 1;

  zCut.resize(nCut);
  if (_tableRead2(is, nCut, zCut.data())) goto label_end;

  stats.resize(nClass * nElem);
  if (_tableRead2(is, nClass * nElem, stats.data())) goto label_end;

  setNCut(nCut);
  setNElem(nElem);
  setZCut(zCut);
  setStats(stats);

  label_end:
  return 0;
}

void AnamDiscrete::setNCut(int ncut)
{
  _nCut = ncut;
  _resize();
}

void AnamDiscrete::setZCut(const VectorDouble& zcut)
{
  _nCut = zcut.size();
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
