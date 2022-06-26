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

#include "Anamorphosis/EAnam.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/NamingConvention.hpp"

class ECalcMember;
class Db;
class Selectivity;

class GSTLEARN_EXPORT AAnam : public AStringable, public ASerializable
{
public:
  AAnam();
  AAnam(const AAnam &m);
  AAnam& operator= (const AAnam &m);
  virtual ~AAnam();

  /// Interface for AAnam
  virtual const EAnam& getType() const = 0;
  virtual bool   hasFactor() const { return false; }
  virtual int    getNFactor() const { return 0; }
  virtual bool   isChangeSupportDefined() const = 0;
  virtual VectorDouble z2factor(double /*z*/, const VectorInt& /*nfact*/) const;
  virtual double getBlockVariance(double /*sval*/, double /*power*/ = 1) const;
  virtual int    updatePointToBlock(double /*r_coef*/);
  virtual bool   allowChangeSupport() const { return false; }
  virtual bool   hasGaussian() const { return false; }
  virtual double RawToTransformValue(double z) const;
  virtual double TransformToRawValue(double y) const;

  double calculateR(double cvv, double power);
  void recoveryLocal(Db *db,
                      int iech0,
                      int iptr,
                      const VectorInt& codes,
                      const VectorInt& qt_vars,
                      double zestim,
                      double zstdev,
                      const Selectivity& calest);
  int codeAnalyze(bool verbose,
                  const VectorInt& codes,
                  int nb_est,
                  int nb_std,
                  int ncut,
                  double proba,
                  int flag_inter,
                  VectorInt& qt_vars) const;
  int DbZToFactor(Db *db,
                  const VectorInt& ifacs,
                  const NamingConvention& namconv = NamingConvention("Factor"));

protected:
  bool _isSampleSkipped(Db *db,
                        int iech,
                        const VectorInt& cols_est,
                        const VectorInt& cols_std);

private:
  bool _isNcutValid(int ncut) const;
  bool _isProbaValid(double proba) const;
  void _printQTvars(const char *title, int type, int number) const;
};
