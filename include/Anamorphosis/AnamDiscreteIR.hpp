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

#include "Anamorphosis/AnamDiscrete.hpp"
#include "Stats/Selectivity.hpp"

namespace gstlrn
{
class Db;

class GSTLEARN_EXPORT AnamDiscreteIR: public AnamDiscrete
{
public:
  AnamDiscreteIR(double rcoef = 0.);
  AnamDiscreteIR(const AnamDiscreteIR& m);
  AnamDiscreteIR& operator=(const AnamDiscreteIR& m);
  virtual ~AnamDiscreteIR();

  /// ICloneable Interface
  IMPLEMENT_CLONING(AnamDiscreteIR)

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ASerializable Interface
  static AnamDiscreteIR* createFromNF(const String& NFFilename, bool verbose = true);

  /// AAnam Interface
  const EAnam& getType() const override { return EAnam::fromKey("DISCRETE_IR"); }
  bool hasFactor() const override { return true; }
  VectorDouble z2factor(double z, const VectorInt& ifacs) const override;
  double computeVariance(double sval) const override;
  Id updatePointToBlock(double r_coef) override;
  bool allowChangeSupport() const override { return true; }
  bool isChangeSupportDefined() const override { return (_sCoef > 0.); }
  Id fitFromArray(const VectorDouble& tab,
                   const VectorDouble& wt = VectorDouble()) override;

  /// AnamDiscrete Interface
  void calculateMeanAndVariance() override;

  static AnamDiscreteIR* create(double rcoef = 0.);
  void reset(Id ncut,
             double r_coef,
             const VectorDouble& zcut,
             const VectorDouble& stats);

  double getRCoef() const { return _sCoef; }
  void setRCoef(double rcoef) { _sCoef = rcoef; }

  Id factor2Selectivity(Db* db,
                         Selectivity* selectivity,
                         const VectorInt& cols_est,
                         const VectorInt& cols_std,
                         Id iptr0);

protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "AnamDiscreteIR"; }

private:
  Id _stats_residuals(Id verbose,
                       Id nech,
                       const VectorDouble& tab,
                       Id* nsorted,
                       double* mean,
                       double* residuals,
                       double* T,
                       double* Q);
  double _getResidual(Id iclass, double z) const;
  void _globalSelectivity(Selectivity* selectivity);

private:
  double _sCoef;

  friend class Selectivity;
};
} // namespace gstlrn