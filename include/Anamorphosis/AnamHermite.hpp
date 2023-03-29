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

#include "Enum/EAnam.hpp"

#include "Anamorphosis/AnamContinuous.hpp"
#include "Basic/ASerializable.hpp"

class Db;
class Selectivity;
class Model;

class GSTLEARN_EXPORT AnamHermite: public AnamContinuous
{
public:
  AnamHermite(int nbpoly=0, bool flagBound=true, double rCoef=1.);
  AnamHermite(const AnamHermite &m);
  AnamHermite& operator= (const AnamHermite &m);
  virtual ~AnamHermite();

  /// ICloneable Interface
  IMPLEMENT_CLONING(AnamHermite)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface AAnam
  const EAnam&  getType() const override { return EAnam::fromKey("HERMITIAN"); }
  bool hasFactor() const override { return true; }
  int  getNFactor() const override { return getNbPoly(); }
  VectorDouble z2factor(double z, const VectorInt& ifacs) const override;
  double computeVariance(double sval) const override;
  int updatePointToBlock(double r_coef) override;
  bool allowChangeSupport() const override { return true; }
  bool isChangeSupportDefined() const override { return (_rCoef < 1.); }
  int getNClass() const override { return getNbPoly(); }

  int fitFromArray(const VectorDouble &tab,
                   const VectorDouble &wt = VectorDouble()) override;

  /// ASerializable Interface
  static AnamHermite* createFromNF(const String& neutralFilename, bool verbose = true);

  /// AnamContinuous Interface
  double rawToTransformValue(double z) const override;
  double transformToRawValue(double y) const override;
  void   calculateMeanAndVariance() override;

  static AnamHermite* create(int nbpoly=0, bool flagBound=true, double rCoef=1.);

  void reset(double pymin,
             double pzmin,
             double pymax,
             double pzmax,
             double aymin,
             double azmin,
             double aymax,
             double azmax,
             double r,
             const VectorDouble &psi_hn);

  int    getNbPoly() const { return (int) _psiHn.size(); }
  VectorDouble getPsiHns() const;
  double getPsiHn(int ih) const;
  double getRCoef() const { return _rCoef; }
  bool   getFlagBound() const { return _flagBound; }

  void   setPsiHns(const VectorDouble& psi_hn) { _psiHn = psi_hn; }
  void   setFlagBound(bool flagBound) { _flagBound = flagBound; }
  void   setPsiHn(int i, double psi_hn);
  void   setRCoef(double r_coef);

  int factor2Selectivity(Db *db,
                         Selectivity* selectivity,
                         const VectorInt& cols_est,
                         const VectorInt& cols_std,
                         int iptr0);

  double evalSupportCoefficient(int option,
                                Model* model,
                                const VectorDouble &dxs,
                                const VectorInt &ndisc,
                                const VectorDouble& angles = VectorDouble(),
                                bool verbose = true);

  VectorDouble cumulateVarianceRatio(double chh) const;

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "AnamHermite"; }

private:
  bool _isIndexValid(int i) const;
  void _defineBounds(double pymin,
                     double pzmin,
                     double pymax,
                     double pzmax,
                     double aymin,
                     double azmin,
                     double aymax,
                     double azmax);
  int _data_sort(int nech,
                 const VectorDouble& z,
                 const VectorDouble& wt,
                 VectorDouble& zs,
                 VectorDouble& ys);
  void _globalSelectivity(Selectivity* selectivity);

private:
  bool   _flagBound;
  double _rCoef;
  VectorDouble _psiHn;

  friend class Selectivity;
};
