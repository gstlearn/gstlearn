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

#include "Enum/EAnam.hpp"

#include "Basic/ASerializable.hpp"
#include "Anamorphosis/AnamContinuous.hpp"

class GSTLEARN_EXPORT AnamEmpirical: public AnamContinuous
{
public:
  AnamEmpirical(int ndisc = 100,
                double sigma2e = TEST,
                bool flagDilution = false,
                bool flagGaussian = true);
  AnamEmpirical(const AnamEmpirical &m);
  AnamEmpirical& operator= (const AnamEmpirical &m);
  virtual ~AnamEmpirical();

  /// ICloneable Interface
  IMPLEMENT_CLONING(AnamEmpirical)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ASerializable Interface
  static AnamEmpirical* createFromNF(const String& neutralFilename,
                                     bool verbose = true);

  void reset(int ndisc,
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
  int getNFactor() const override { return _nDisc; }
  int fitFromArray(const VectorDouble &tab,
                   const VectorDouble &wt = VectorDouble()) override;

  /// AnamContinuous Interface
  void    calculateMeanAndVariance() override;
  double  rawToTransformValue(double zz) const override;
  double  transformToRawValue(double yy) const override;
  bool    isChangeSupportDefined() const override { return false; }

  AnamEmpirical* create(int ndisc = 100, double sigma2e = TEST);
  int    getNDisc() const { return _nDisc; }
  double getSigma2e() const { return _sigma2e; }
  const  VectorDouble& getZDisc() const { return _ZDisc; }
  const  VectorDouble& getYDisc() const { return _YDisc; }
  bool isFlagDilution() const { return _flagDilution; }
  bool isFlagGaussian() const { return _flagGaussian; }

  void setSigma2e(double sigma2e) { _sigma2e = sigma2e; }
  void setNDisc(int ndisc);
  void setDisc(const VectorDouble& zdisc, const VectorDouble& ydisc);
  void setFlagDilution(bool flagDilution) { _flagDilution = flagDilution; }
  void setFlagGaussian(bool flagGaussian) { _flagGaussian = flagGaussian; }

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "AnamEmpirical"; }

private:
  int _getStatistics(const VectorDouble &tab,
                     int *count,
                     double *mean,
                     double *mean2,
                     double *mini,
                     double *maxi,
                     double *var);
  int _fitWithDilutionGaussian(const VectorDouble &tab);
  int _fitWithDilutionLognormal(const VectorDouble &tab);
  int _fitNormalScore(const VectorDouble &tab);

private:
  bool   _flagDilution;
  bool   _flagGaussian;
  int    _nDisc;
  double _sigma2e;
  VectorDouble _ZDisc;
  VectorDouble _YDisc;
};
