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

class Db;
class ECalcMember;

class GSTLEARN_EXPORT AnamHermite: public AnamContinuous
{
public:
  AnamHermite(int nbpoly=0, bool flagBound=true, double rCoef=1.);
  AnamHermite(const AnamHermite &m);
  AnamHermite& operator= (const AnamHermite &m);
  virtual ~AnamHermite();

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface AAnam
  const EAnam&  getType() const override { return EAnam:: HERMITIAN; }
  double modifyCov(const ECalcMember& member,
                   int iclass,
                   double dist,
                   double cov0,
                   double cov1,
                   double cov2) const;

  /// ASerializable Interface
  int dumpToNF(const String& neutralFilename, bool verbose = false) const;
  static AnamHermite* createFromNF(const String& neutralFilename, bool verbose = false);

  /// AnamContinuous Interface
  double RawToGaussianValue(double z) const override;
  double GaussianToRawValue(double y) const override;
  void   calculateMeanAndVariance() override;

  static AnamHermite* create(int nbpoly=0, bool flagBound=true, double rCoef=1.);

  int    getNbPoly() const { return _nbPoly; }
  const  VectorDouble& getPsiHn() const { return _psiHn; }
  double getPsiHn(int i) const;
  double getRCoef() const { return _rCoef; }
  bool   getFlagBound() const { return _flagBound; }

  void   setNbPoly(int nbPoly) { _nbPoly = nbPoly; };
  void   setPsiHn(VectorDouble psi_hn);
  void   setFlagBound(bool flagBound) { _flagBound = flagBound; }
  void   setPsiHn(int i, double psi_hn);
  void   setRCoef(double r_coef) { _rCoef = r_coef; }

  double calculateVarianceFromPsi(double chh);
  int    fit(const VectorDouble& tab,
             const VectorDouble& wt = VectorDouble());
  int    fit(Db *db, const ELoc& locatorType = ELoc::Z);
  int    fit(Db *db, const String& name);

protected:
  /// ASerializable Interface
  virtual int _deserialize(std::istream& is, bool verbose) override;
  virtual int _serialize(std::ostream& os, bool verbose = false) const override;

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
  int    _nbPoly;
  bool   _flagBound;
  double _rCoef;
  VectorDouble _psiHn;
};
