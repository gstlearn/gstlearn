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
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"
#include "Anamorphosis/AAnam.hpp"
#include "Basic/Interval.hpp"

/**
 * Output structure
 */
class AnamContinuousFit
{
  private:
  VectorDouble _y;
  VectorDouble _z;
  VectorDouble _aylim;
  VectorDouble _azlim;
  VectorDouble _pylim;
  VectorDouble _pzlim;

  public:
  const VectorDouble& getAylim() const { return _aylim; }
  void setAylim(const VectorDouble& aylim) { _aylim = aylim; }
  const VectorDouble& getAzlim() const { return _azlim; }
  void setAzlim(const VectorDouble& azlim) { _azlim = azlim; }
  const VectorDouble& getPylim() const { return _pylim; }
  void setPylim(const VectorDouble& pylim) { _pylim = pylim; }
  const VectorDouble& getPzlim() const { return _pzlim; }
  void setPzlim(const VectorDouble& pzlim) { _pzlim = pzlim; }
  const VectorDouble& getY() const { return _y; }
  void setY(const VectorDouble& y) { _y = y; }
  const VectorDouble& getZ() const { return _z; }
  void setZ(const VectorDouble& z) { _z = z; }
};

class GSTLEARN_EXPORT AnamContinuous: public AAnam
{
public:
  AnamContinuous();
  AnamContinuous(const AnamContinuous &m);
  AnamContinuous& operator=(const AnamContinuous &m);
  virtual ~AnamContinuous();

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// AAnam interface
  bool hasGaussian() const override { return true; }

  /// Interface for AnamContinuous
  virtual void   calculateMeanAndVariance();
  virtual double getVariance() const override { return _variance; }

  void setABounds(double azmin = TEST,
                  double azmax = TEST,
                  double aymin = TEST,
                  double aymax = TEST);
  void setPBounds(double pzmin = TEST,
                  double pzmax = TEST,
                  double pymin = TEST,
                  double pymax = TEST);

  VectorDouble rawToGaussianVector(const VectorDouble &z) const;
  VectorDouble gaussianToRawVector(const VectorDouble &y) const;

  AnamContinuousFit sample(int ndisc = 100,
                           double aymin = -10,
                           double aymax = +10);

  double getMean() const { return _mean; }

  double getAymax() const { return _ay.getVmax(); }
  double getAymin() const { return _ay.getVmin(); }
  double getAzmax() const { return _az.getVmax(); }
  double getAzmin() const { return _az.getVmin(); }
  double getPymax() const { return _py.getVmax(); }
  double getPymin() const { return _py.getVmin(); }
  double getPzmax() const { return _pz.getVmax(); }
  double getPzmin() const { return _pz.getVmin(); }

  void setAzmin(double azmin) { _az.setVmin(azmin); }
  void setAzmax(double azmax) { _az.setVmax(azmax); }
  void setAymin(double aymin) { _ay.setVmin(aymin); }
  void setAymax(double aymax) { _ay.setVmax(aymax); }
  void setPzmin(double pzmin) { _pz.setVmin(pzmin); }
  void setPzmax(double pzmax) { _pz.setVmax(pzmax); }
  void setPymin(double pymin) { _py.setVmin(pymin); }
  void setPymax(double pymax) { _py.setVmax(pymax); }
  void setMean(double mean)   { _mean = mean; }
  void setVariance(double variance) { _variance = variance; }

protected:
  /// ASerializable Interface
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "AnamContinuous"; }

protected:
  Interval _az;
  Interval _ay;
  Interval _pz;
  Interval _py;
  double _mean;
  double _variance;
};
