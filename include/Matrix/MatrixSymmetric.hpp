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

class AMatrix;
class EOperator;

/**
 * Square Symmetric matrices
 */
class GSTLEARN_EXPORT MatrixSymmetric : public MatrixSquare {

public:
  MatrixSymmetric(int nrow = 0);
  MatrixSymmetric(const MatrixSymmetric &m);
  MatrixSymmetric(const AMatrix &m);
  MatrixSymmetric& operator= (const MatrixSymmetric &m);
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

  static MatrixSymmetric* createFromVVD(const VectorVectorDouble &X);
  static MatrixSymmetric* createFromVD(const VectorDouble &X);
  static MatrixSymmetric* createFromTLTU(int neq,
                                               const VectorDouble &tl);
  static MatrixSymmetric* createFromTriangle(int mode,
                                                   int neq,
                                                   const VectorDouble &tl);
  static MatrixSymmetric* createRandomDefinitePositive(int neq, int seed = 13242);
  static MatrixSymmetric* sample(const MatrixSymmetric* A,
                                       const VectorInt& rowKeep,
                                       bool flagInvert = false);

  int computeEigen(bool optionPositive = true);
  int computeGeneralizedEigen(const MatrixSymmetric& b, bool optionPositive = true);
  int computeGeneralizedInverse(MatrixSymmetric &tabout,
                                double maxicond = 1.e20,
                                double eps = EPSILON20);
  bool isDefinitePositive();
  int minimizeWithConstraintsInPlace(const VectorDouble& gmat,
                                     const MatrixDense& aemat,
                                     const VectorDouble& bemat,
                                     const MatrixDense& aimat,
                                     const VectorDouble& bimat,
                                     VectorDouble& xmat);

  virtual bool _isPhysicallyPresent(int irow, int icol) const override;
  virtual void _setValues(const double* values, bool byCol = true) override;
  virtual int _invert() override;

  // Local functions (old style algebra)
  int _matrix_qo(const VectorDouble& gmat, VectorDouble& xmat);
  int _matrix_qoc(bool flag_invert,
                  const VectorDouble& gmat,
                  int na,
                  const MatrixDense& amat,
                  const VectorDouble& bmat,
                  VectorDouble& xmat,
                  VectorDouble& lambda);
  int _constraintsError(const VectorInt& active,
                        const MatrixDense& aimat,
                        const VectorDouble& bimat,
                        const VectorDouble& xmat,
                        VectorDouble& vmat,
                        VectorInt& flag);
  static int _constraintsConcatenateMat(int nae,
                                        int nai,
                                        int neq,
                                        const VectorInt& active,
                                        const MatrixDense& tabemat,
                                        const MatrixDense& tabimat,
                                        MatrixDense& tabout);
  static int _constraintsConcatenateVD(int nae,
                                       int nai,
                                       const VectorInt& active,
                                       const VectorDouble& tabemat,
                                       const VectorDouble& tabimat,
                                       VectorDouble& tabout);
  static int _constraintsCount(int nai, VectorInt& active);
  int _terminateEigen(const VectorDouble& eigenValues,
                      const VectorDouble& eigenVectors,
                      bool optionPositive = true,
                      bool changeOrder    = false);
  MatrixSymmetric compress0MatLC(const MatrixDense& matLC);

private:
  int _getTriangleSize() const;
};
