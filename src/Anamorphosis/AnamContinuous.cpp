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
#include "Anamorphosis/AnamContinuous.hpp"
#include "Db/Db.hpp"
#include "Basic/String.hpp"
#include "Basic/NamingConvention.hpp"
#include "geoslib_f.h"

AnamContinuous::AnamContinuous(const EAnam& type)
    : Anam(type),
      _az(),
      _ay(),
      _pz(),
      _py(),
      _mean(TEST),
      _variance(TEST)
{
}

AnamContinuous::AnamContinuous(const AnamContinuous &m)
    : Anam(m),
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

String AnamContinuous::toString(int level) const
{
  std::stringstream sstr;
  sstr << Anam::toString(level);

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

VectorDouble AnamContinuous::RawToGaussianVector(const VectorDouble& z) const
{
  int number = static_cast<int> (z.size());
  VectorDouble y;
  y.resize(number);
  for (int i=0; i<number; i++)
    y[i] = RawToGaussianValue(z[i]);
  return y;
}

VectorDouble AnamContinuous::GaussianToRawVector(const VectorDouble& y) const
{
  int number = static_cast<int> (y.size());
  VectorDouble z;
  z.resize(number);
  for (int i=0; i<number; i++)
    z[i] = GaussianToRawValue(y[i]);
  return z;
}

int AnamContinuous::RawToGaussian(Db *db, const String& name, NamingConvention namconv)
{
  if (db == nullptr) return 1;
  VectorString exp_names = expandList(db->getNames(),name);
  if (exp_names.empty()) return 1;

  int nadd = static_cast<int> (exp_names.size());
  int iatt = db->addFields(nadd);
  if (iatt < 0) return 1;

  for (int i = 0; i < nadd; i++)
  {
    VectorDouble z = db->getField(exp_names[i],true);
    if (z.size() <= 0) return 1;
    VectorDouble y = RawToGaussianVector(z);
    db->setFieldByAttribute(y,iatt+i,true);
    namconv.setNamesAndLocators(exp_names[i],db,iatt+i,String(),1,false);
  }

  namconv.setLocators(db, iatt, nadd);
  return 0;
}

int AnamContinuous::RawToGaussian(Db *db, const ELoc& locatorType, NamingConvention namconv)
{
  if (db == nullptr) return 1;
  int number = db->getLocatorNumber(locatorType);
  if (number <= 0) return 1;

  int iatt = db->addFields(number);

  for (int item=0; item<number; item++)
  {
    VectorDouble z = db->getFieldByLocator(locatorType,item,true);
    if (z.size() <= 0) continue;
    VectorDouble y = RawToGaussianVector(z);
    db->setFieldByAttribute(y,iatt+item,true);
  }
  namconv.setNamesAndLocators(db,locatorType,-1,db,iatt,String(),1,true);
  return 0;
}

int AnamContinuous::GaussianToRaw(Db *db,const String& name, NamingConvention namconv)
{
  if (db == nullptr) return 1;
  VectorString exp_names = expandList(db->getNames(), name);
  if (exp_names.empty()) return 1;

  int nadd = static_cast<int> (exp_names.size());
  int iatt = db->addFields(nadd);
  if (iatt < 0) return 1;

  for (int i = 0; i < nadd; i++)
  {
    VectorDouble y = db->getField(exp_names[i],true);
    if (y.size() <= 0) return 1;
    VectorDouble z = GaussianToRawVector(y);
    db->setFieldByAttribute(z,iatt+i,true);
    namconv.setNamesAndLocators(exp_names[i],db,iatt+i,String(),1,false);
  }

  namconv.setLocators(db, iatt, nadd);
  return 0;
}

int AnamContinuous::GaussianToRaw(Db *db, const ELoc& locatorType, NamingConvention namconv)
{
  if (db == nullptr) return 1;
  int number = db->getLocatorNumber(locatorType);
  if (number <= 0) return 1;

  int iatt = db->addFields(number);

  for (int item=0; item<number; item++)
  {
    VectorDouble y = db->getFieldByLocator(locatorType,item,true);
    if (y.size() <= 0) continue;
    VectorDouble z = GaussianToRawVector(y);
    db->setFieldByAttribute(z,iatt+item,true);
  }

  namconv.setNamesAndLocators(db,locatorType,-1,db,iatt,String(),1,true);
  return 0;
}

/**
 * Calculate Anamorphosis function for a set of Y-values
 * @param ndisc Number of discretization steps
 * @param aymin Minimum Y value
 * @param aymax Maximum Y value
 * @return AnamContinuousFit structure
 */
AnamContinuousFit AnamContinuous::sample(int ndisc, double aymin, double aymax)
{
  VectorDouble y(ndisc);
  VectorDouble z(ndisc);
  double pas  = (aymax - aymin) / ndisc;
  int ind0 = ndisc / 2;
  y[ind0] = 0.;
  z[ind0] = GaussianToRawValue(y[ind0]);

  /* Calculating the values below y=0 */

  for (int ind = ind0-1; ind >= 0; ind--)
  {
    y[ind] = y[ind+1] - pas;
    z[ind] = GaussianToRawValue(y[ind]);
  }

/* Calculating the values above y=0 */

  for (int ind = ind0+1; ind<ndisc; ind++)
  {
    y[ind] = y[ind-1] + pas;
    z[ind] = GaussianToRawValue(y[ind]);
  }

  // Preparing the returned structure
  AnamContinuousFit retfit;
  retfit.y = y;
  retfit.z = z;
  retfit.aylim = _ay.getBounds();
  retfit.azlim = _az.getBounds();
  retfit.pylim = _py.getBounds();
  retfit.pzlim = _pz.getBounds();

  return retfit;
}
