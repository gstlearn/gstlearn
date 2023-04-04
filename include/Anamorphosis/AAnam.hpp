/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "geoslib_define.h"

#include "Enum/EAnam.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/NamingConvention.hpp"

class Db;

class GSTLEARN_EXPORT AAnam : public AStringable, public ASerializable, public ICloneable
{
public:
  AAnam();
  AAnam(const AAnam &m);
  AAnam& operator= (const AAnam &m);
  virtual ~AAnam();

  /// Interface for AAnam
  virtual const EAnam& getType() const = 0;
  virtual double       getVariance() const { return TEST; }
  virtual bool         hasFactor() const { return false; }
  virtual int          getNFactor() const { return 0; }
  virtual int          getNClass() const { return 0; }
  virtual bool         isChangeSupportDefined() const = 0;
  virtual VectorDouble z2factor(double z, const VectorInt& ifqcs) const;
  virtual double       computeVariance(double sval) const;
  virtual int          updatePointToBlock(double r_coef);
  virtual bool         allowChangeSupport() const { return false; }
  virtual bool         hasGaussian() const { return false; }
  virtual double       rawToTransformValue(double z) const;
  virtual double       transformToRawValue(double y) const;
  virtual int          fitFromArray(const VectorDouble &tab,
                                    const VectorDouble &wt = VectorDouble()) { DECLARE_UNUSED(tab,wt); return 0;}

  double invertVariance(double cvv) const;
  VectorDouble rawToTransformVec(const VectorDouble& z) const;
  VectorDouble transformToRawVec(const VectorDouble& z) const;

  int fitFromLocator(Db *db, const ELoc& locatorType = ELoc::fromKey("Z"));
  int fit(Db *db, const String& name);

  int rawToGaussianByLocator(Db *db,
                             const NamingConvention &namconv = NamingConvention(
                                 "Y"));
  int rawToGaussian(Db *db,
                    const String &name,
                    const NamingConvention &namconv = NamingConvention("Y"));
  int normalScore(Db *db,
                  const String& name,
                  const NamingConvention &namconv = NamingConvention("Gaussian"));
  int gaussianToRawByLocator(Db *db,
                             const NamingConvention &namconv = NamingConvention("Z"));
  int gaussianToRaw(Db *db,
                    const String &name,
                    const NamingConvention &namconv = NamingConvention("Z"));

  int rawToFactorByRanks(Db *db,
                         const VectorInt &ifacs,
                         const NamingConvention &namconv = NamingConvention(
                             "Factor"));
  int rawToFactor(Db *db,
                  int nfactor,
                  const NamingConvention &namconv = NamingConvention("Factor"));

protected:
  bool _isSampleSkipped(Db *db,
                        int iech,
                        const VectorInt& cols_est,
                        const VectorInt& cols_std);
  bool _isFitted() const { return _flagFitted; }

private:
  bool _isNcutValid(int ncut) const;
  bool _isProbaValid(double proba) const;
  void _printQTvars(const char *title, int type, int number) const;

private:
  bool _flagFitted;
};
