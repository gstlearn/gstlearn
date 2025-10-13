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

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Enum/EAnam.hpp"

#include "Anamorphosis/AnamContinuous.hpp"

/**
 * Gaussian Anamorphosis using Empirical Method
 *
 * This class is meant in order to construct the transfer function from Raw to Gaussian scale
 * directly based on the data.
 *
 * It essentially maps the cumulative function (CDF) of the raw values into the CDF of
 * the theoretical Gaussian distribution.
 *
 * This can be performed directly on the experimental CDF (normal score) or by diluting the data
 * values beforehand. In the latter solution, each (valid) datum is replaced by a small local
 * distribution. This is meant to smooth the stepwise CDF.
 *
 * The dilution function (implemented at any data point) can be either a Gaussian or a Lognormal one.
 * In the Gaussian case, the variance (width of the dilution function) is considered as constant
 * (either provided by the user or defaulted by the program)*
 * In the lognormal case, the logarithmic variance is constant (hence the width is proportional
 * to the square of the value).
 */

namespace gstlrn
{
class GSTLEARN_EXPORT AnamEmpirical: public AnamContinuous
{
public:
  AnamEmpirical(Id ndisc         = 100,
                double sigma2e    = TEST,
                bool flagDilution = false,
                bool flagGaussian = true);
  AnamEmpirical(const AnamEmpirical& m);
  AnamEmpirical& operator=(const AnamEmpirical& m);
  virtual ~AnamEmpirical();

  /// ICloneable Interface
  IMPLEMENT_CLONING(AnamEmpirical)

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ASerializable Interface
  static AnamEmpirical* createFromNF(const String& NFFilename, bool verbose = true);

  void reset(Id ndisc,
             double pymin,
             double pzmin,
             double pymax,
             double pzmax,
             double aymin,
             double azmin,
             double aymax,
             double azmax,
             double sigma2e,
             const VectorDouble& zdisc,
             const VectorDouble& ydisc);

  /// AAnam Interface
  const EAnam& getType() const override { return EAnam::fromKey("EMPIRICAL"); }
  Id getNFactor() const override { return _nDisc; }
  Id fitFromArray(const VectorDouble& tab,
                   const VectorDouble& wt = VectorDouble()) override;

  /// AnamContinuous Interface
  void calculateMeanAndVariance() override;
  double rawToTransformValue(double zz) const override;
  double transformToRawValue(double yy) const override;
  bool isChangeSupportDefined() const override { return false; }

  static AnamEmpirical* create(Id ndisc = 100, double sigma2e = TEST);
  Id getNDisc() const { return _nDisc; }
  double getSigma2e() const { return _sigma2e; }
  const VectorDouble& getZDisc() const { return _ZDisc; }
  const VectorDouble& getYDisc() const { return _YDisc; }
  bool isFlagDilution() const { return _flagDilution; }
  bool isFlagGaussian() const { return _flagGaussian; }

  void setSigma2e(double sigma2e) { _sigma2e = sigma2e; }
  void setNDisc(Id ndisc);
  void setDisc(const VectorDouble& zdisc, const VectorDouble& ydisc);
  void setFlagDilution(bool flagDilution) { _flagDilution = flagDilution; }
  void setFlagGaussian(bool flagGaussian) { _flagGaussian = flagGaussian; }

protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "AnamEmpirical"; }

private:
  static Id _getStatistics(const VectorDouble& tab,
                            Id* count,
                            double* mean,
                            double* mean2,
                            double* mini,
                            double* maxi,
                            double* var);
  Id _fitWithDilutionGaussian(const VectorDouble& tab);
  Id _fitWithDilutionLognormal(const VectorDouble& tab);
  Id _fitNormalScore(const VectorDouble& tab);

private:
  bool _flagDilution;
  bool _flagGaussian;
  Id _nDisc;
  double _sigma2e;
  VectorDouble _ZDisc;
  VectorDouble _YDisc;
};
} // namespace gstlrn