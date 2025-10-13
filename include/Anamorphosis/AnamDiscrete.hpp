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

#include "Anamorphosis/AAnam.hpp"
#include "Matrix/MatrixDense.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class Selectivity;

class AAnam; // Forward declaration

class GSTLEARN_EXPORT AnamDiscrete: public AAnam
{
public:
  AnamDiscrete();
  AnamDiscrete(const AnamDiscrete& m);
  AnamDiscrete& operator=(const AnamDiscrete& m);
  virtual ~AnamDiscrete();

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// AAnam interface
  bool hasGaussian() const override { return false; }
  Id getNClass() const override { return _nCut + 1; }

  /// Interface for AnamDiscrete
  virtual void calculateMeanAndVariance();
  double getVariance() const override { return _variance; }

  Id getNCut() const { return _nCut; }
  Id getNElem() const { return _nElem; }
  const VectorDouble& getZCut() const { return _zCut; }
  double getZCut(Id i) const { return _zCut[i]; }
  double getMean() const { return _mean; }

  void setMean(double mean) { _mean = mean; }
  void setVariance(double variance) { _variance = variance; }
  void setNCut(Id ncut);
  void setZCut(const VectorDouble& zcut);
  void setNElem(Id nelem);
  void setStats(const VectorDouble& stats);

  // Function for using Stats in DD anamorphosis
  double getDDStatProp(Id iclass) const;
  double getDDStatZmoy(Id iclass) const;
  double getDDStatCnorm(Id iclass) const;
  double getDDStatLambda(Id iclass) const;
  double getDDStatU(Id iclass) const;
  double getDDStatMul(Id iclass) const;
  void setDDStatProp(Id iclass, double value);
  void setDDStatZmoy(Id iclass, double value);
  void setDDStatCnorm(Id iclass, double value);
  void setDDStatLambda(Id iclass, double value);
  void setDDStatU(Id iclass, double value);
  void setDDStatMul(Id iclass, double value);

  // Function for using Stats in IR anamorphosis
  double getIRStatT(Id iclass) const;
  double getIRStatQ(Id iclass) const;
  double getIRStatZ(Id iclass) const;
  double getIRStatB(Id iclass) const;
  double getIRStatR(Id iclass) const;
  double getIRStatRV(Id iclass) const;
  void setIRStatT(Id iclass, double value);
  void setIRStatQ(Id iclass, double value);
  void setIRStatZ(Id iclass, double value);
  void setIRStatB(Id iclass, double value);
  void setIRStatR(Id iclass, double value);
  void setIRStatRV(Id iclass, double value);

  const MatrixDense& getStats() const { return _stats; }

protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "AnamDiscrete"; }

  bool _isClassValid(Id iclass) const;
  void _resize();

private:
  Id _nCut;
  Id _nElem;
  double _mean;
  double _variance;
  VectorDouble _zCut;
  MatrixDense _stats;
};
} // namespace gstlrn