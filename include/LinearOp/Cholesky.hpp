/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "LinearOp/ALinearOp.hpp"
#include "Basic/VectorNumT.hpp"

class cs; /// TODO : Dependency to csparse to be removed
class css;
class csn;
class GSTLEARN_EXPORT Cholesky: public ALinearOp
{
public:
  Cholesky(const cs* mat = nullptr, bool flagDecompose=true);
  Cholesky(const Cholesky &m);
  Cholesky& operator=(const Cholesky &m);
  virtual ~Cholesky();

  int getSize() const override;
  void evalInverse(const VectorDouble& inv, VectorDouble& outv) const override;

  int  reset(const cs* mat = nullptr, bool flagDecompose = true);
  void simulate(VectorDouble& inv, VectorDouble& outv);
  void stdev(VectorDouble& vcur, bool flagStDev = false);
  void printout(const char *title, bool verbose = false) const;
  double computeLogDet() const;

protected:
  void _evalDirect(const VectorDouble& inv, VectorDouble& outv) const override;

private:
  bool _isDefined() const { return _mat != nullptr; }
  bool _isDecomposed() const { return _matS != nullptr && _matN != nullptr; }
  void _clean();
  void _decompose(bool verbose = false) const;

private:
  const cs* _mat; // Copy of the pointer
  mutable css *_matS;
  mutable csn *_matN;
  mutable VectorDouble _work;
};
