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

#include "Anamorphosis/AnamContinuous.hpp"
#include "Anamorphosis/EAnam.hpp"
#include "Basic/ASerializable.hpp"
#include "Stats/Selectivity.hpp"

class Db;

class GSTLEARN_EXPORT AnamHermite: public AnamContinuous
{
public:
  AnamHermite(int nbpoly=0, bool flagBound=true, double rCoef=1.);
  AnamHermite(const AnamHermite &m);
  AnamHermite& operator= (const AnamHermite &m);
  virtual ~AnamHermite();

  /// IClonable Interface
  virtual IClonable* clone() const override { return new AnamHermite(*this); };

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface AAnam
  const EAnam&  getType() const override { return EAnam::HERMITIAN; }
  bool hasFactor() const override { return true; }
  int getNFactor() const override { return getNbPoly(); }
  VectorDouble z2factor(double z, const VectorInt& ifacs) const override;
  double computeVariance(double sval) const override;
  int updatePointToBlock(double r_coef) override;
  bool allowChangeSupport() const override { return true; }
  bool isChangeSupportDefined() const override { return (_rCoef < 1.); }

  /// ASerializable Interface
  static AnamHermite* createFromNF(const String& neutralFilename, bool verbose = true);

  /// AnamContinuous Interface
  double RawToTransformValue(double z) const override;
  double TransformToRawValue(double y) const override;
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
  const  VectorDouble& getPsiHn() const { return _psiHn; }
  double getPsiHn(int i) const;
  double getRCoef() const { return _rCoef; }
  bool   getFlagBound() const { return _flagBound; }

  void   setPsiHn(VectorDouble psi_hn) { _psiHn = psi_hn; }
  void   setFlagBound(bool flagBound) { _flagBound = flagBound; }
  void   setPsiHn(int i, double psi_hn);
  void   setRCoef(double r_coef) { _rCoef = r_coef; }

  int    fit(const VectorDouble& tab,
             const VectorDouble& wt = VectorDouble());
  int    fit(Db *db, const ELoc& locatorType = ELoc::Z);
  int    fit(Db *db, const String& name);

  Selectivity calculateSelectivity(const VectorDouble& zcut);

  int factor2QT(Db *db,
                const VectorDouble& cutmine,
                const VectorInt& cols_est,
                const VectorInt& cols_std,
                int iptr,
                const VectorInt& codes,
                VectorInt& qt_vars);

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

private:
  bool   _flagBound;
  double _rCoef;
  VectorDouble _psiHn;
};
