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
#include "Anamorphosis/AnamDiscrete.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Stats/Selectivity.hpp"

#include <cmath>

#define ANAM_KD_NELEM 6

namespace gstlrn
{
AnamDiscrete::AnamDiscrete()
  : AAnam()
  , _nCut(0)
  , _nElem(ANAM_KD_NELEM)
  , _mean(TEST)
  , _variance(TEST)
  , _zCut()
  , _stats()
{
  _resize();
}

AnamDiscrete::AnamDiscrete(const AnamDiscrete& m)
  : AAnam(m)
  , _nCut(m._nCut)
  , _nElem(m._nElem)
  , _mean(m._mean)
  , _variance(m._variance)
  , _zCut(m._zCut)
  , _stats(m._stats)
{
}

AnamDiscrete& AnamDiscrete::operator=(const AnamDiscrete& m)
{
  if (this != &m)
  {
    AAnam::operator=(m);
    _nCut     = m._nCut;
    _nElem    = m._nElem;
    _zCut     = m._zCut;
    _mean     = m._mean;
    _variance = m._variance;
    _stats    = m._stats;
  }
  return *this;
}

AnamDiscrete::~AnamDiscrete()
{
}

void AnamDiscrete::_resize()
{
  auto ncut   = getNCut();
  auto nclass = getNClass();
  auto nelem  = getNElem();

  _zCut.resize(ncut, 0.);
  _stats.resetFromValue(nclass, nelem, 0.);
}

String AnamDiscrete::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  if (!_isFitted()) return sstr.str();

  sstr << "Number of cutoffs = " << _nCut << std::endl;
  sstr << "Number of classes = " << getNClass() << std::endl;
  if (!FFFF(_mean))
    sstr << "Mean              = " << _mean << std::endl;
  if (!FFFF(_variance))
    sstr << "Variance          = " << _variance << std::endl;

  sstr << std::endl;
  sstr << toMatrix("Cutoffs", VectorString(), VectorString(), true, _nCut, 1, _zCut);

  sstr << toMatrix(String(), VectorString(), VectorString(), true, getNClass(), getNElem(),
                   getStats().getValues());

  return sstr.str();
}

void AnamDiscrete::calculateMeanAndVariance()
{
  _mean     = TEST;
  _variance = TEST;
}

double AnamDiscrete::getDDStatProp(Id iclass) const
{
  if (!_isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass, 0);
}
double AnamDiscrete::getDDStatZmoy(Id iclass) const
{
  if (!_isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass, 1);
}
double AnamDiscrete::getDDStatCnorm(Id iclass) const
{
  if (!_isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass, 2);
}
double AnamDiscrete::getDDStatLambda(Id iclass) const
{
  if (!_isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass, 3);
}
double AnamDiscrete::getDDStatU(Id iclass) const
{
  if (!_isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass, 4);
}
double AnamDiscrete::getDDStatMul(Id iclass) const
{
  if (!_isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass, 5);
}
void AnamDiscrete::setDDStatProp(Id iclass, double value)
{
  if (!_isClassValid(iclass)) return;
  _stats.setValue(iclass, 0, value);
}
void AnamDiscrete::setDDStatZmoy(Id iclass, double value)
{
  if (!_isClassValid(iclass)) return;
  _stats.setValue(iclass, 1, value);
}
void AnamDiscrete::setDDStatCnorm(Id iclass, double value)
{
  if (!_isClassValid(iclass)) return;
  _stats.setValue(iclass, 2, value);
}
void AnamDiscrete::setDDStatLambda(Id iclass, double value)
{
  if (!_isClassValid(iclass)) return;
  _stats.setValue(iclass, 3, value);
}
void AnamDiscrete::setDDStatU(Id iclass, double value)
{
  if (!_isClassValid(iclass)) return;
  _stats.setValue(iclass, 4, value);
}
void AnamDiscrete::setDDStatMul(Id iclass, double value)
{
  if (!_isClassValid(iclass)) return;
  _stats.setValue(iclass, 5, value);
}

// Function for using Stats in IR anamorphosis
double AnamDiscrete::getIRStatT(Id iclass) const
{
  if (!_isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass, 0);
}
double AnamDiscrete::getIRStatQ(Id iclass) const
{
  if (!_isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass, 1);
}
double AnamDiscrete::getIRStatZ(Id iclass) const
{
  if (!_isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass, 2);
}
double AnamDiscrete::getIRStatB(Id iclass) const
{
  if (!_isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass, 3);
}
double AnamDiscrete::getIRStatR(Id iclass) const
{
  if (!_isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass, 4);
}
double AnamDiscrete::getIRStatRV(Id iclass) const
{
  if (!_isClassValid(iclass)) return TEST;
  return _stats.getValue(iclass, 5);
}
void AnamDiscrete::setIRStatT(Id iclass, double value)
{
  if (!_isClassValid(iclass)) return;
  _stats.setValue(iclass, 0, value);
}
void AnamDiscrete::setIRStatQ(Id iclass, double value)
{
  if (!_isClassValid(iclass)) return;
  _stats.setValue(iclass, 1, value);
}
void AnamDiscrete::setIRStatZ(Id iclass, double value)
{
  if (!_isClassValid(iclass)) return;
  _stats.setValue(iclass, 2, value);
}
void AnamDiscrete::setIRStatB(Id iclass, double value)
{
  if (!_isClassValid(iclass)) return;
  _stats.setValue(iclass, 3, value);
}
void AnamDiscrete::setIRStatR(Id iclass, double value)
{
  if (!_isClassValid(iclass)) return;
  _stats.setValue(iclass, 4, value);
}
void AnamDiscrete::setIRStatRV(Id iclass, double value)
{
  if (!_isClassValid(iclass)) return;
  _stats.setValue(iclass, 5, value);
}

bool AnamDiscrete::_isClassValid(Id iclass) const
{
  return checkArg("Class Index", iclass, getNClass());
}

bool AnamDiscrete::_serializeAscii(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;
  ret      = ret && _recordWrite<Id>(os, "Number of Cuttofs", getNCut());
  ret      = ret && _recordWrite<Id>(os, "Number of classes", getNClass());
  ret      = ret && _recordWrite<Id>(os, "Number of elements", getNElem());
  ret      = ret && _tableWrite(os, "Cutoff value", getNCut(), getZCut());
  ret      = ret && _tableWrite(os, "DD Stats", getNClass() * getNElem(), getStats().getValues());
  return ret;
}

bool AnamDiscrete::_deserializeAscii(std::istream& is, bool /*verbose*/)
{
  VectorDouble zCut, stats;
  Id nCut   = 0;
  Id nClass = 0;
  Id nElem  = 0;

  bool ret = true;
  ret      = ret && _recordRead<Id>(is, "Number of Cutoffs", nCut);
  ret      = ret && _recordRead<Id>(is, "Number of Classes", nClass);
  ret      = ret && _recordRead<Id>(is, "Number of Statistic Columns", nElem);

  if (ret)
  {
    zCut.resize(nCut);
    ret = ret && _tableRead(is, "Cutoff value", nCut, zCut.data());
  }

  if (ret)
  {
    stats.resize(nClass * nElem);
    ret = ret && _tableRead(is, "DD Stats", nClass * nElem, stats.data());
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

void AnamDiscrete::setNCut(Id ncut)
{
  _nCut = ncut;
  _resize();
}

void AnamDiscrete::setZCut(const VectorDouble& zcut)
{
  _nCut = static_cast<Id>(zcut.size());
  _resize();
  _zCut = zcut;
};

void AnamDiscrete::setNElem(Id nelem)
{
  _nElem = nelem;
  _resize();
}

void AnamDiscrete::setStats(const VectorDouble& stats)
{
  auto nclass = getNClass();
  auto nelem  = getNElem();
  if (static_cast<Id>(stats.size()) != nclass * nelem)
  {
    messerr("Argument 'stats' incorrect. Its dimension (%d) should be %d * %d",
            stats.size(), nclass, nelem);
    return;
  }
  _stats.setValues(stats);
}
#ifdef HDF5
bool AnamDiscrete::_deserializeH5(H5::Group& grp, [[maybe_unused]] bool verbose)
{
  auto anamG = SerializeHDF5::getGroup(grp, "AnamDiscrete");
  if (!anamG)
  {
    return false;
  }

  /* Read the grid characteristics */
  bool ret  = true;
  Id ncut   = 0;
  Id nclass = 0;
  Id nelem  = 0;
  VectorDouble zcuts;
  VectorDouble stats;

  ret = ret && SerializeHDF5::readValue(*anamG, "NCut", ncut);
  ret = ret && SerializeHDF5::readValue(*anamG, "NClass", nclass);
  ret = ret && SerializeHDF5::readValue(*anamG, "NElem", nelem);
  ret = ret && SerializeHDF5::readVec(*anamG, "Cuts", zcuts);
  ret = ret && SerializeHDF5::readVec(*anamG, "Stats", stats);

  if (ret)
  {
    setNCut(ncut);
    setNElem(nelem);
    setZCut(zcuts);
    setStats(stats);
  }

  return ret;
}

bool AnamDiscrete::_serializeH5(H5::Group& grp, [[maybe_unused]] bool verbose) const
{
  auto anamG = grp.createGroup("AnamDiscrete");

  bool ret = true;

  ret = ret && SerializeHDF5::writeValue(anamG, "NCut", getNCut());
  ret = ret && SerializeHDF5::writeValue(anamG, "NClass", getNClass());
  ret = ret && SerializeHDF5::writeValue(anamG, "NElem", getNElem());
  ret = ret && SerializeHDF5::writeVec(anamG, "Cuts", getZCut());
  ret = ret && SerializeHDF5::writeVec(anamG, "Stats", getStats().getValues());

  return ret;
}
#endif
} // namespace gstlrn