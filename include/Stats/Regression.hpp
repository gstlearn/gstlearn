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
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void setCoeffs(const VectorDouble& coeffs) { _coeffs = coeffs; }
  void setCount(int count) { _count = count; }
  void setFlagCst(bool flagCst) { _flagCst = flagCst; }
  void setNvar(int nvar) { _nvar = nvar; }
  void setVariance(double variance) { _variance = variance; }
  void setVarres(double varres) { _varres = varres; }

  VectorDouble getCoeffs() const { return _coeffs; }
  double getCoeff(int i) const { return _coeffs[i]; }
  int getNvar() const { return _nvar; }
  int getCount() const { return _count; }
  double getVariance() const { return _variance; }
  double getVarres() const { return _varres; }

  int apply(Db *db1,
            int iptr0,
            const String &nameResp,
            const VectorString &nameAux,
            int mode = 0,
            bool flagCst = false,
            Db *db2 = nullptr,
            const Model *model = nullptr) const;

private:
  int _count;
  int _nvar;
  bool _flagCst;
  VectorDouble _coeffs;
  double _variance;
  double _varres;
};

GSTLEARN_EXPORT Regression regression(Db *db1,
                                      const String &nameResp,
                                      const VectorString &nameAux = VectorString(),
                                      int mode = 0,
                                      bool flagCst = false,
                                      Db *db2 = nullptr,
                                      const Model *model = nullptr);
GSTLEARN_EXPORT VectorDouble regressionDeming(const VectorDouble &x,
                                              const VectorDouble &y,
                                              double delta = 1);

