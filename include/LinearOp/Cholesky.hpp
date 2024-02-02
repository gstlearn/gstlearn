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

#include "LinearOp/ALinearOp.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Basic/VectorNumT.hpp"

class css; /// TODO : Dependency to csparse to be removed
class csn;
class GSTLEARN_EXPORT Cholesky: public ALinearOp
{
public:
  Cholesky(const MatrixSparse* mat = nullptr, bool flagDecompose=true);
  Cholesky(const Cholesky &m);
  Cholesky& operator=(const Cholesky &m);
  virtual ~Cholesky();

  int getSize() const override;
  void evalInverse(const VectorDouble& vecin, VectorDouble& vecout) const override;

  bool isDefined() const { return _matCS != nullptr; }
  bool isCholeskyDecomposed() const { return _matS != nullptr && _matN != nullptr; }

  int  reset(const MatrixSparse* mat = nullptr, bool flagDecompose = true);
  void decompose(bool verbose = false) const;

  void simulate(VectorDouble& vecin, VectorDouble& vecout);
  void stdev(VectorDouble& vcur, bool flagStDev = false);
  void printout(const char *title, bool verbose = false) const;
  double computeLogDet() const;

protected:
  void _evalDirect(const VectorDouble& inv, VectorDouble& outv) const override;

private:
  void _clean();

private:
  const MatrixSparse* _matCS; // Copy of the pointer
  mutable css *_matS;
  mutable csn *_matN;
  mutable VectorDouble _work;
};
