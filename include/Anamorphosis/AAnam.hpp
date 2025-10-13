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
#include "Basic/ICloneable.hpp"
#include "Basic/NamingConvention.hpp"


namespace gstlrn
{
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
  virtual Id          getNFactor() const { return 0; }
  virtual Id          getNClass() const { return 0; }
  virtual bool         isChangeSupportDefined() const = 0;
  virtual VectorDouble z2factor(double z, const VectorInt& ifacs) const;
  virtual double       computeVariance(double sval) const;
  virtual Id          updatePointToBlock(double r_coef);
  virtual bool         allowChangeSupport() const { return false; }
  virtual bool         hasGaussian() const { return false; }
  virtual double       rawToTransformValue(double z) const;
  virtual double       transformToRawValue(double y) const;
  virtual Id          fitFromArray(const VectorDouble &tab,
                                    const VectorDouble &wt = VectorDouble()) { DECLARE_UNUSED(tab,wt); return 0;}

  double invertVariance(double cvv) const;
  VectorDouble rawToTransformVec(const VectorDouble& z) const;
  VectorDouble transformToRawVec(const VectorDouble& y) const;

  Id fitFromLocator(Db *db, const ELoc& locatorType = ELoc::fromKey("Z"));
  Id fit(Db *db, const String& name);

  Id rawToGaussianByLocator(Db *db,
                             const NamingConvention &namconv = NamingConvention("Y"));
  Id rawToGaussian(Db *db,
                    const String &name,
                    const NamingConvention &namconv = NamingConvention("Y"));
  Id normalScore(Db *db,
                  const String& name,
                  const NamingConvention &namconv = NamingConvention("Gaussian"));
  Id gaussianToRawByLocator(Db *db,
                             const NamingConvention &namconv = NamingConvention("Z"));
  Id gaussianToRaw(Db *db,
                    const String &name,
                    const NamingConvention &namconv = NamingConvention("Z"));

  Id rawToFactorByRanks(Db *db,
                         const VectorInt &ifacs,
                         const NamingConvention &namconv = NamingConvention(
                             "Factor"));
  Id rawToFactor(Db *db,
                  Id nfactor,
                  const NamingConvention &namconv = NamingConvention("Factor"));

protected:
  static bool _isSampleSkipped(Db* db,
                               Id iech,
                               const VectorInt& cols_est,
                               const VectorInt& cols_std);
  bool _isFitted() const { return _flagFitted; }

private:
  static bool _isNcutValid(Id ncut);
  static bool _isProbaValid(double proba);
  static void _printQTvars(const char *title, Id type, Id number);

private:
  bool _flagFitted;
};
}