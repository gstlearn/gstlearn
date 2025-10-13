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

#include "Basic/VectorNumT.hpp"
#include "Matrix/MatrixSquare.hpp"

namespace gstlrn
{
class AMatrix;
class EOperator;

/**
 * Square Symmetric matrices
 */
class GSTLEARN_EXPORT MatrixSymmetric: public MatrixSquare
{

public:
  MatrixSymmetric(Id nrow = 0);
  MatrixSymmetric(const MatrixSymmetric& m);
  MatrixSymmetric(const AMatrix& m);
  MatrixSymmetric& operator=(const MatrixSymmetric& m);
  virtual ~MatrixSymmetric();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// ICloneable interface
  IMPLEMENT_CLONING(MatrixSymmetric)

  /// Interface to AMatrix
  bool mustBeSymmetric() const final { return true; }
  bool isSymmetric(double eps = EPSILON10, bool printWhyNot = false) const final
  {
    DECLARE_UNUSED(printWhyNot);
    DECLARE_UNUSED(eps);
    return true;
  }
  void resetFromVVD(const VectorVectorDouble& tab, bool byCol = true) override;

  void normMatrix(const AMatrix& y, const MatrixSquare& x = MatrixSquare(), bool transpose = false);

  static MatrixSymmetric* createFromVVD(const VectorVectorDouble& X);
  static MatrixSymmetric* createFromVD(const VectorDouble& X);
  static MatrixSymmetric* createFromTLTU(Id neq,
                                         const VectorDouble& tl);
  static MatrixSymmetric* createFromTriangle(Id mode,
                                             Id neq,
                                             const VectorDouble& tl);
  static MatrixSymmetric* createRandomDefinitePositive(Id neq, Id seed = 13242);
  static bool sample(MatrixSymmetric& res,
                     const MatrixSymmetric& A,
                     const VectorInt& rowKeep,
                     bool flagInvert = false);

  Id computeEigen(bool optionPositive = true);
  Id computeGeneralizedEigen(const MatrixSymmetric& b, bool optionPositive = true);
  Id computeGeneralizedInverse(MatrixSymmetric& tabout,
                               double maxicond = 1.e20,
                               double eps      = EPSILON20);
  bool isDefinitePositive();
  Id minimizeWithConstraintsInPlace(const VectorDouble& gmat,
                                    const MatrixDense& aemat,
                                    const VectorDouble& bemat,
                                    const MatrixDense& aimat,
                                    const VectorDouble& bimat,
                                    VectorDouble& xmat);

  bool _isPhysicallyPresent(Id irow, Id icol) const override;
  void _setValues(const double* values, bool byCol = true) override;
  Id _invert() override;

  // Local functions (old style algebra)
  Id _matrix_qo(const VectorDouble& gmat, VectorDouble& xmat);
  Id _matrix_qoc(bool flag_invert,
                 const VectorDouble& gmat,
                 Id na,
                 const MatrixDense& amat,
                 const VectorDouble& bmat,
                 VectorDouble& xmat,
                 VectorDouble& lambda);
  Id _constraintsError(const VectorInt& active,
                       const MatrixDense& aimat,
                       const VectorDouble& bimat,
                       const VectorDouble& xmat,
                       VectorDouble& vmat,
                       VectorInt& flag);
  static Id _constraintsConcatenateMat(Id nae,
                                       Id nai,
                                       Id neq,
                                       const VectorInt& active,
                                       const MatrixDense& tabemat,
                                       const MatrixDense& tabimat,
                                       MatrixDense& tabout);
  static Id _constraintsConcatenateVD(Id nae,
                                      Id nai,
                                      const VectorInt& active,
                                      const VectorDouble& tabemat,
                                      const VectorDouble& tabimat,
                                      VectorDouble& tabout);
  static Id _constraintsCount(Id nai, VectorInt& active);
  Id _terminateEigen(const VectorDouble& eigenValues,
                     const VectorDouble& eigenVectors,
                     bool optionPositive = true,
                     bool changeOrder    = false);
  MatrixSymmetric compress0MatLC(const MatrixDense& matLC);

private:
  Id _getTriangleSize() const;
};
} // namespace gstlrn