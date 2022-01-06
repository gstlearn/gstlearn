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

class Db;

class GSTLEARN_EXPORT AnamHermite: public AnamContinuous
{
public:
  AnamHermite(int nbpoly=0, bool flagBound=true, double rCoef=1.);
  AnamHermite(const AnamHermite &m);
  AnamHermite& operator= (const AnamHermite &m);
  virtual ~AnamHermite();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int    getNbPoly() const { return _nbPoly; }
  const  VectorDouble& getPsiHn() const { return _psiHn; }
  double getPsiHn(int i) const;
  double getRCoef() const { return _rCoef; }
  bool   getFlagBound() const { return _flagBound; }

  void   setPsiHn(VectorDouble psi_hn);
  void   setFlagBound(bool flagBound) { _flagBound = flagBound; }
  void   setPsiHn(int i, double psi_hn);
  void   setRCoef(double r_coef) { _rCoef = r_coef; }

  double RawToGaussianValue(double z) const override;
  double GaussianToRawValue(double y) const override;
  double calculateVarianceFromPsi(double chh);
  void   calculateMeanAndVariance() override;
  int    fit(const VectorDouble& tab,
             const VectorDouble& wt = VectorDouble());
  int    fit(Db *db, const ELoc& locatorType = ELoc::Z);
  int    fit(Db *db, const String& name);

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
