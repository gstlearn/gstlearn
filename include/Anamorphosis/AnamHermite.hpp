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

#include "Enum/EAnam.hpp"

#include "Anamorphosis/AnamContinuous.hpp"

namespace gstlrn
{
class Db;
class Selectivity;
class Model;

class GSTLEARN_EXPORT AnamHermite: public AnamContinuous
{
public:
  AnamHermite(Id nbpoly = 3, bool flagBound = true, double rCoef = 1.);
  AnamHermite(const AnamHermite& m);
  AnamHermite& operator=(const AnamHermite& m);
  virtual ~AnamHermite();

  /// ICloneable Interface
  IMPLEMENT_CLONING(AnamHermite)

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface AAnam
  const EAnam& getType() const override { return EAnam::fromKey("HERMITIAN"); }
  bool hasFactor() const override { return true; }
  Id getNFactor() const override { return getNbPoly(); }
  VectorDouble z2factor(double z, const VectorInt& ifacs) const override;
  double computeVariance(double chh) const override;
  Id updatePointToBlock(double r_coef) override;
  bool allowChangeSupport() const override { return true; }
  bool isChangeSupportDefined() const override { return (_rCoef < 1.); }
  Id getNClass() const override { return getNbPoly(); }

  Id fitFromArray(const VectorDouble& tab,
                   const VectorDouble& wt = VectorDouble()) override;

  /// ASerializable Interface
  static AnamHermite* createFromNF(const String& NFFilename, bool verbose = true);

  /// AnamContinuous Interface
  double rawToTransformValue(double z) const override;
  double transformToRawValue(double y) const override;
  void calculateMeanAndVariance() override;

  static AnamHermite* create(Id nbpoly = 0, bool flagBound = true, double rCoef = 1.);

  void reset(double pymin,
             double pzmin,
             double pymax,
             double pzmax,
             double aymin,
             double azmin,
             double aymax,
             double azmax,
             double r,
             const VectorDouble& psi_hn);

  Id getNbPoly() const { return static_cast<Id>(_psiHn.size()); }
  VectorDouble getPsiHns() const;
  double getPsiHn(Id ih) const;
  double getRCoef() const { return _rCoef; }
  bool getFlagBound() const { return _flagBound; }

  void setPsiHns(const VectorDouble& psi_hn) { _psiHn = psi_hn; }
  void setFlagBound(bool flagBound) { _flagBound = flagBound; }
  void setPsiHn(Id i, double psi_hn);
  void setRCoef(double r_coef);

  Id factor2Selectivity(Db* db,
                         Selectivity* selectivity,
                         const VectorInt& cols_est,
                         const VectorInt& cols_std,
                         Id iptr0);

  double evalSupportCoefficient(Id option,
                                Model* model,
                                const VectorDouble& dxs,
                                const VectorInt& ndisc,
                                const VectorDouble& angles = VectorDouble(),
                                bool verbose               = true);

  VectorDouble cumulateVarianceRatio(double chh) const;

protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "AnamHermite"; }

private:
  bool _isIndexValid(Id i) const;
  void _defineBounds(double pymin,
                     double pzmin,
                     double pymax,
                     double pzmax,
                     double aymin,
                     double azmin,
                     double aymax,
                     double azmax);
  static Id _data_sort(Id nech,
                        const VectorDouble& z,
                        const VectorDouble& wt,
                        VectorDouble& zs,
                        VectorDouble& ys);
  void _globalSelectivity(Selectivity* selectivity);

private:
  bool _flagBound;
  double _rCoef;
  VectorDouble _psiHn;

  friend class Selectivity;
};
} // namespace gstlrn