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

#include "Anamorphosis/AnamDiscrete.hpp"
#include "Stats/Selectivity.hpp"

class Db;

class GSTLEARN_EXPORT AnamDiscreteIR: public AnamDiscrete
{
public:
  AnamDiscreteIR(double rcoef = 0.);
  AnamDiscreteIR(const AnamDiscreteIR &m);
  AnamDiscreteIR& operator= (const AnamDiscreteIR &m);
  virtual ~AnamDiscreteIR();

  /// ICloneable Interface
  IMPLEMENT_CLONING(AnamDiscreteIR)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ASerializable Interface
  static AnamDiscreteIR* createFromNF(const String& neutralFilename, bool verbose = true);

  /// AAnam Interface
  const EAnam&  getType() const override { return EAnam::fromKey("DISCRETE_IR"); }
  bool hasFactor() const override { return true; }
  VectorDouble z2factor(double z, const VectorInt& ifacs) const override;
  double computeVariance(double sval) const override;
  int updatePointToBlock(double r_coef) override;
  bool allowChangeSupport() const override { return true; }
  bool isChangeSupportDefined() const override { return (_sCoef > 0.); }
  int fitFromArray(const VectorDouble &tab,
                   const VectorDouble &wt = VectorDouble()) override;

  /// AnamDiscrete Interface
  void calculateMeanAndVariance() override;

  static AnamDiscreteIR* create(double rcoef = 0.);
  void reset(int ncut,
             double r_coef,
             const VectorDouble &zcut,
             const VectorDouble &stats);

  double getRCoef() const { return _sCoef; }
  void   setRCoef(double rcoef) { _sCoef = rcoef; }

  int factor2Selectivity(Db *db,
                         Selectivity* selectivity,
                         const VectorInt& cols_est,
                         const VectorInt& cols_std,
                         int iptr0);

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "AnamDiscreteIR"; }

private:
  int _stats_residuals(int verbose,
                       int nech,
                       const VectorDouble& tab,
                       int *nsorted,
                       double *mean,
                       double *residuals,
                       double *T,
                       double *Q);
  double _getResidual(int iclass, double z) const;
  void _globalSelectivity(Selectivity* selectivity);

private:
  double _sCoef;

  friend class Selectivity;
};
