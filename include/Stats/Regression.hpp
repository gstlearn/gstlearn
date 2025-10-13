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
#include "Basic/AStringable.hpp"

#include "geoslib_define.h"

namespace gstlrn
{
class Db;
class Model;

class GSTLEARN_EXPORT Regression: public AStringable
{
public:
  Regression();
  Regression(const Regression& r);
  Regression& operator=(const Regression& r);
  virtual ~Regression();

  /// Interface for AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  void setCoeffs(const VectorDouble& coeffs) { _coeffs = coeffs; }
  void setCount(Id count) { _count = count; }
  void setFlagCst(bool flagCst) { _flagCst = flagCst; }
  void setNvar(Id nvar) { _nvar = nvar; }
  void setVariance(double variance) { _variance = variance; }
  void setVarres(double varres) { _varres = varres; }

  VectorDouble getCoeffs() const { return _coeffs; }
  double getCoeff(Id i) const { return _coeffs[i]; }
  Id getNvar() const { return _nvar; }
  Id getCount() const { return _count; }
  double getVariance() const { return _variance; }
  double getVarres() const { return _varres; }

  Id apply(Db *db1,
            Id iptr0,
            const String &nameResp,
            const VectorString &nameAux,
            Id mode = 0,
            bool flagCst = false,
            Db *db2 = nullptr,
            const Model *model = nullptr) const;

private:
  Id _count;
  Id _nvar;
  bool _flagCst;
  VectorDouble _coeffs;
  double _variance;
  double _varres;
};

GSTLEARN_EXPORT Regression regression(Db *db1,
                                      const String &nameResp,
                                      const VectorString &nameAux = VectorString(),
                                      Id mode = 0,
                                      bool flagCst = false,
                                      Db *db2 = nullptr,
                                      const Model *model = nullptr);
GSTLEARN_EXPORT VectorDouble regressionDeming(const VectorDouble &x,
                                              const VectorDouble &y,
                                              double delta = 1);

}