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
  virtual double       RawToTransformValue(double z) const;
  virtual double       TransformToRawValue(double y) const;
  virtual int          fitFromArray(const VectorDouble &tab,
                                    const VectorDouble &wt = VectorDouble()) { SYMBOL_UNUSED(tab,wt); return 0;}

  double invertVariance(double cvv) const;
  VectorDouble RawToTransformVec(const VectorDouble& z) const;
  VectorDouble TransformToRawVec(const VectorDouble& z) const;

  int fitFromLocator(Db *db, const ELoc& locatorType = ELoc::fromKey("Z"));
  int fit(Db *db, const String& name);

  int RawToGaussianByLocator(Db *db,
                             const NamingConvention &namconv = NamingConvention(
                                 "Y"));
  int RawToGaussian(Db *db,
                    const String &name,
                    const NamingConvention &namconv = NamingConvention("Y"));
  int NormalScore(Db *db,
                  const NamingConvention &namconv = NamingConvention("Gaussian"));
  int GaussianToRawByLocator(Db *db,
                             const NamingConvention &namconv = NamingConvention("Z"));
  int GaussianToRaw(Db *db,
                    const String &name,
                    const NamingConvention &namconv = NamingConvention("Z"));

  int RawToFactor(Db *db,
                  const VectorInt &ifacs,
                  const NamingConvention &namconv = NamingConvention("Factor"));
  int RawToFactor(Db *db,
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
