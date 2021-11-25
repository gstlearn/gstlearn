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
#pragma once

#include "gstlearn_export.hpp"
#include "Anamorphosis/Anam.hpp"
#include "Basic/Interval.hpp"
#include "Basic/Vector.hpp"
#include "Basic/NamingConvention.hpp"
#include "Db/ELoc.hpp"

/**
 * Output structure
 */
typedef struct {
  VectorDouble y;
  VectorDouble z;
  VectorDouble aylim;
  VectorDouble azlim;
  VectorDouble pylim;
  VectorDouble pzlim;
} AnamContinuousFit;

class GSTLEARN_EXPORT AnamContinuous: public Anam
{
public:
  AnamContinuous(const EAnam& type = EAnam::UNDEFINED);
  AnamContinuous(const AnamContinuous &m);
  AnamContinuous& operator= (const AnamContinuous &m);
  virtual ~AnamContinuous();

  virtual String toString(int level) const override;

  void setABounds(double azmin = TEST,
                  double azmax = TEST,
                  double aymin = TEST,
                  double aymax = TEST);
  void setPBounds(double pzmin = TEST,
                  double pzmax = TEST,
                  double pymin = TEST,
                  double pymax = TEST);

  virtual void calculateMeanAndVariance();
  virtual double RawToGaussianValue(double z) const = 0;
  virtual double GaussianToRawValue(double y) const = 0;
  VectorDouble RawToGaussianVector(const VectorDouble& z) const;
  VectorDouble GaussianToRawVector(const VectorDouble& y) const;
  int RawToGaussian(Db *db, const String& name,NamingConvention namconv = NamingConvention("Y"));
  int GaussianToRaw(Db *db, const String& name,NamingConvention namconv = NamingConvention("Z"));
  int RawToGaussian(Db *db, const ELoc& locatorType = ELoc::Z, NamingConvention namconv = NamingConvention("Y"));
  int GaussianToRaw(Db *db, const ELoc& locatorType = ELoc::Z, NamingConvention namconv = NamingConvention("Z"));

  AnamContinuousFit sample(int ndisc = 100,
                           double aymin = -10,
                           double aymax = +10);

  double getMean()     const { return _mean; }
  double getVariance() const { return _variance; }
  double getAymax()    const { return _ay.getVmax(); }
  double getAymin()    const { return _ay.getVmin(); }
  double getAzmax()    const { return _az.getVmax(); }
  double getAzmin()    const { return _az.getVmin(); }
  double getPymax()    const { return _py.getVmax(); }
  double getPymin()    const { return _py.getVmin(); }
  double getPzmax()    const { return _pz.getVmax(); }
  double getPzmin()    const { return _pz.getVmin(); }

  void setAzmin(double azmin)       { _az.setVmin(azmin); }
  void setAzmax(double azmax)       { _az.setVmax(azmax); }
  void setAymin(double aymin)       { _ay.setVmin(aymin); }
  void setAymax(double aymax)       { _ay.setVmax(aymax); }
  void setPzmin(double pzmin)       { _pz.setVmin(pzmin); }
  void setPzmax(double pzmax)       { _pz.setVmax(pzmax); }
  void setPymin(double pymin)       { _py.setVmin(pymin); }
  void setPymax(double pymax)       { _py.setVmax(pymax); }
  void setMean(double mean)         { _mean = mean; }
  void setVariance(double variance) { _variance = variance; }

protected:
  Interval _az;
  Interval _ay;
  Interval _pz;
  Interval _py;
  double   _mean;
  double   _variance;
};
