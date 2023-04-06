/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Enum/EAnam.hpp"

#include "Anamorphosis/AnamDiscrete.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Stats/PCA.hpp"
#include "Stats/Selectivity.hpp"

class GSTLEARN_EXPORT AnamDiscreteDD: public AnamDiscrete
{
public:
  AnamDiscreteDD(double mu = 1., double scoef = 0.);
  AnamDiscreteDD(const AnamDiscreteDD &m);
  AnamDiscreteDD& operator= (const AnamDiscreteDD &m);
  virtual ~AnamDiscreteDD();

  /// ICloneable Interface
  IMPLEMENT_CLONING(AnamDiscreteDD)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ASerializable Interface
  static AnamDiscreteDD* createFromNF(const String& neutralFilename, bool verbose = true);

  /// AAnam Interface
  const EAnam&  getType() const override { return EAnam::fromKey("DISCRETE_DD"); }
  bool hasFactor() const override { return true; }
  int getNFactor() const override { return 0; }
  VectorDouble z2factor(double z, const VectorInt& ifacs) const override;
  double computeVariance(double sval) const override;
  int  updatePointToBlock(double r_coef) override;
  bool allowChangeSupport() const override { return true; }
  bool isChangeSupportDefined() const override { return (_sCoef > 0.); }
  int fitFromArray(const VectorDouble &tab,
                   const VectorDouble &wt = VectorDouble()) override;

  /// AnamDiscrete Interface
  void calculateMeanAndVariance() override;

  VectorDouble factors_exp(bool verbose = false);
  VectorDouble factors_maf(bool verbose = false);
  VectorDouble factors_mod();
  VectorDouble chi2I(const VectorDouble& chi, int mode);

  AnamDiscreteDD* create(double mu = 1., double scoef = 0.);
  void reset(int ncut,
             double scoef,
             double mu,
             const VectorDouble &zcut,
             const VectorDouble &pcaz2f,
             const VectorDouble &pcaf2z,
             const VectorDouble &stats);

  PCA& getMAF() { return _maf; }
  double getMu() const { return _mu; }
  double getSCoef() const { return _sCoef; }
  const VectorDouble& getI2Chi() const { return _i2Chi; }
  VectorDouble getPcaZ2F() const { return _maf.getZ2F(); }
  VectorDouble getPcaF2Z() const { return _maf.getF2Z(); }

  void setMu(double mu) { _mu = mu; }
  void setRCoef(double rcoef) { _sCoef = rcoef; }
  void setPcaZ2F(VectorDouble pcaz2f) { _maf.setPcaZ2F(pcaz2f); }
  void setPcaF2Z(VectorDouble pcaf2z) { _maf.setPcaF2Z(pcaf2z); }
  void setI2Chi(const VectorDouble& i2Chi) { _i2Chi = i2Chi; }

  int factor2Selectivity(Db *db,
                         Selectivity* selectivity,
                         const VectorInt& cols_est,
                         const VectorInt& cols_std,
                         int iptr0);

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "AnamDiscreteDD"; }

private:
  int _stats(int nech, const VectorDouble& tab);
  VectorDouble _generator(const VectorDouble& vecc,
                          const VectorDouble& veca,
                          const VectorDouble& vecb,
                          VectorDouble& eigvec,
                          VectorDouble& eigval);
  void _lambdaToMul();
  void _blockAnamorphosis(const VectorDouble& chi);
  void _globalSelectivity(Selectivity* selectivity);

private:
  double _mu;
  double _sCoef;
  PCA    _maf;
  VectorDouble _i2Chi;

  friend class Selectivity;
};
