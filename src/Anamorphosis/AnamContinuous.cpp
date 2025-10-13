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
#include "Anamorphosis/AnamContinuous.hpp"
#include "Db/Db.hpp"
#include "Basic/SerializeHDF5.hpp"

namespace gstlrn
{
AnamContinuous::AnamContinuous()
    : AAnam(),
      _az(),
      _ay(),
      _pz(),
      _py(),
      _mean(TEST),
      _variance(TEST)
{
}

AnamContinuous::AnamContinuous(const AnamContinuous &m)
    : AAnam(m),
      _az(m._az),
      _ay(m._ay),
      _pz(m._pz),
      _py(m._py),
      _mean(m._mean),
      _variance(m._variance)
{
}

AnamContinuous& AnamContinuous::operator=(const AnamContinuous &m)
{
  if (this != &m)
  {
    AAnam::operator= (m);
    _az = m._az;
    _ay = m._ay;
    _pz = m._pz;
    _py = m._py;
    _mean = m._mean;
    _variance = m._variance;
  }
  return *this;
}

AnamContinuous::~AnamContinuous()
{

}

void AnamContinuous::setABounds(double azmin,
                                double azmax,
                                double aymin,
                                double aymax)
{
  _az.init(azmin, azmax);
  _ay.init(aymin, aymax);
}

void AnamContinuous::setPBounds(double pzmin,
                                double pzmax,
                                double pymin,
                                double pymax)
{
  _pz.init(pzmin, pzmax);
  _py.init(pymin, pymax);
}

String AnamContinuous::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);
  std::stringstream sstr;

  if (! _isFitted()) return sstr.str();

  sstr << "Minimum absolute value for Y  = " << _ay.getVmin() << std::endl;
  sstr << "Maximum absolute value for Y  = " << _ay.getVmax() << std::endl;
  sstr << "Minimum absolute value for Z  = " << _az.getVmin() << std::endl;
  sstr << "Maximum absolute value for Z  = " << _az.getVmax() << std::endl;
  sstr << "Minimum practical value for Y = " << _py.getVmin() << std::endl;
  sstr << "Maximum practical value for Y = " << _py.getVmax() << std::endl;
  sstr << "Minimum practical value for Z = " << _pz.getVmin() << std::endl;
  sstr << "Maximum practical value for Z = " << _pz.getVmax() << std::endl;
  sstr << "Mean                          = " << _mean << std::endl;
  sstr << "Variance                      = " << _variance << std::endl;

  return sstr.str();
}

void AnamContinuous::calculateMeanAndVariance()
{
  _mean = TEST;
  _variance = TEST;
}

VectorDouble AnamContinuous::rawToGaussianVector(const VectorDouble &z) const
{
  Id number = static_cast<Id>(z.size());
  VectorDouble y;
  y.resize(number);
  for (Id i = 0; i < number; i++)
    y[i] = (FFFF(z[i])) ? TEST : rawToTransformValue(z[i]);
  return y;
}

VectorDouble AnamContinuous::gaussianToRawVector(const VectorDouble &y) const
{
  Id number = static_cast<Id>(y.size());
  VectorDouble z;
  z.resize(number);
  for (Id i = 0; i < number; i++)
    z[i] = (FFFF(y[i])) ? TEST : transformToRawValue(y[i]);
  return z;
}

/**
 * Calculate Anamorphosis function for a set of Y-values
 * @param ndisc Number of discretization steps
 * @param aymin Minimum Y value
 * @param aymax Maximum Y value
 * @return AnamContinuousFit structure
 */
AnamContinuousFit AnamContinuous::sample(Id ndisc, double aymin, double aymax)
{
  VectorDouble y(ndisc);
  VectorDouble z(ndisc);
  double pas = (aymax - aymin) / ndisc;
  Id ind0 = ndisc / 2;
  y[ind0] = 0.;
  z[ind0] = transformToRawValue(y[ind0]);

  /* Calculating the values below y=0 */

  for (Id ind = ind0 - 1; ind >= 0; ind--)
  {
    y[ind] = y[ind + 1] - pas;
    z[ind] = transformToRawValue(y[ind]);
  }

  /* Calculating the values above y=0 */

  for (Id ind = ind0 + 1; ind < ndisc; ind++)
  {
    y[ind] = y[ind - 1] + pas;
    z[ind] = transformToRawValue(y[ind]);
  }

  // Preparing the returned structure
  AnamContinuousFit retfit;
  retfit.setY(y);
  retfit.setZ(z);
  retfit.setAylim(_ay.getBounds());
  retfit.setAzlim(_az.getBounds());
  retfit.setPylim(_py.getBounds());
  retfit.setPzlim(_pz.getBounds());

  return retfit;
}

bool AnamContinuous::_serializeAscii(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;
  ret = ret && _recordWrite<double>(os,"", getAzmin());
  ret = ret && _recordWrite<double>(os, "Absolute Values for Z", getAzmax());
  ret = ret && _recordWrite<double>(os, "", getAymin());
  ret = ret && _recordWrite<double>(os, "Absolute Values for Y", getAymax());
  ret = ret && _recordWrite<double>(os, "", getPzmin());
  ret = ret && _recordWrite<double>(os, "Practical Values for Z", getPzmax());
  ret = ret && _recordWrite<double>(os, "", getPymin());
  ret = ret && _recordWrite<double>(os, "Practical Values for Y", getPymax());
  ret = ret && _recordWrite<double>(os, "Calculated mean", getMean());
  ret = ret && _recordWrite<double>(os, "Calculated variance", getVariance());

  return ret;
}

bool AnamContinuous::_deserializeAscii(std::istream& is, bool /*verbose*/)
{
  double azmin = 0.;
  double azmax = 0.;
  double aymin = 0.;
  double aymax = 0.;
  double pzmin = 0.;
  double pzmax = 0.;
  double pymin = 0.;
  double pymax = 0;
  double mean = TEST;
  double variance = TEST;

  bool ret = true;
  ret = ret && _recordRead<double>(is, "Minimum absolute Z-value", azmin);
  ret = ret && _recordRead<double>(is, "Maximum absolute Z-value", azmax);
  ret = ret && _recordRead<double>(is, "Minimum absolute Y-value", aymin);
  ret = ret && _recordRead<double>(is, "Maximum absolute Y-value", aymax);
  ret = ret && _recordRead<double>(is, "Minimum Experimental Z-value", pzmin);
  ret = ret && _recordRead<double>(is, "Maximum Experimental Z-value", pzmax);
  ret = ret && _recordRead<double>(is, "Minimum Experimental Y-value", pymin);
  ret = ret && _recordRead<double>(is, "Maximum Experimental Y-value", pymax);
  ret = ret && _recordRead<double>(is, "Experimental Mean", mean);
  ret = ret && _recordRead<double>(is, "Experimental Variance", variance);
  if (! ret) return ret;

  setPymin(pymin);
  setPzmin(pzmin);
  setPymax(pymax);
  setPzmax(pzmax);
  setAymin(aymin);
  setAzmin(azmin);
  setAymax(aymax);
  setAzmax(azmax);
  setMean(mean);
  setVariance(variance);

  return ret;
}
#ifdef HDF5
bool AnamContinuous::_deserializeH5(H5::Group& grp, [[maybe_unused]] bool verbose)
{
  auto anamG = SerializeHDF5::getGroup(grp, "AnamContinuous");
  if (!anamG)
  {
    return false;
  }

  /* Read the grid characteristics */
  bool ret = true;
  double azmin = 0.;
  double azmax = 0.;
  double aymin = 0.;
  double aymax = 0.;
  double pzmin = 0.;
  double pzmax = 0.;
  double pymin = 0.;
  double pymax = 0.;
  double mean = 0.;
  double variance = 0.;

  // ret      = ret && SerializeHDF5::readVec(*anamG, "NX", nx);

  ret = ret && SerializeHDF5::readValue(*anamG, "Azmin", azmin);
  ret = ret && SerializeHDF5::readValue(*anamG, "Azmax", azmax);
  ret = ret && SerializeHDF5::readValue(*anamG, "Aymin", aymin);
  ret = ret && SerializeHDF5::readValue(*anamG, "Aymax", aymax);

  ret = ret && SerializeHDF5::readValue(*anamG, "Pzmin", pzmin);
  ret = ret && SerializeHDF5::readValue(*anamG, "Pzmax", pzmax);
  ret = ret && SerializeHDF5::readValue(*anamG, "Pymin", pymin);
  ret = ret && SerializeHDF5::readValue(*anamG, "Pymax", pymax);

  ret = ret && SerializeHDF5::readValue(*anamG, "Mean", mean);
  ret = ret && SerializeHDF5::readValue(*anamG, "Variance", variance);

  if (ret)
  {
    setPymin(pymin);
    setPzmin(pzmin);
    setPymax(pymax);
    setPzmax(pzmax);
    setAymin(aymin);
    setAzmin(azmin);
    setAymax(aymax);
    setAzmax(azmax);
    setMean(mean);
    setVariance(variance);
  }

  return ret;
}

bool AnamContinuous::_serializeH5(H5::Group& grp, [[maybe_unused]] bool verbose) const
{
  auto anamG = grp.createGroup("AnamContinuous");

  bool ret = true;

  ret = ret && SerializeHDF5::writeValue(anamG, "Azmin", getAzmin());
  ret = ret && SerializeHDF5::writeValue(anamG, "Azmax", getAzmax());
  ret = ret && SerializeHDF5::writeValue(anamG, "Aymin", getAymin());
  ret = ret && SerializeHDF5::writeValue(anamG, "Aymax", getAymax());

  ret = ret && SerializeHDF5::writeValue(anamG, "Pzmin", getPzmin());
  ret = ret && SerializeHDF5::writeValue(anamG, "Pzmax", getPzmax());
  ret = ret && SerializeHDF5::writeValue(anamG, "Pymin", getPymin());
  ret = ret && SerializeHDF5::writeValue(anamG, "Pymax", getPymax());

  ret = ret && SerializeHDF5::writeValue(anamG, "Mean", getMean());
  ret = ret && SerializeHDF5::writeValue(anamG, "Variance", getVariance());

  return ret;
}
#endif
}