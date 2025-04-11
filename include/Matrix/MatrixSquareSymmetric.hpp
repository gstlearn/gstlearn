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
#include "Matrix/AMatrixSquare.hpp"

class AMatrix;
class EOperator;

/**
 * Square Symmetric matrices
 */
class GSTLEARN_EXPORT MatrixSquareSymmetric : public AMatrixSquare {

public:
  MatrixSquareSymmetric(int nrow = 0);
  MatrixSquareSymmetric(const MatrixSquareSymmetric &m);
  MatrixSquareSymmetric(const AMatrix &m);
  MatrixSquareSymmetric& operator= (const MatrixSquareSymmetric &m);
	virtual ~MatrixSquareSymmetric();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// ICloneable interface
  IMPLEMENT_CLONING(MatrixSquareSymmetric)

  /// Interface to AMatrix
  bool mustBeSymmetric() const final { return true; }
  bool isSymmetric(double eps = EPSILON10, bool printWhyNot = false) const final
  {
    DECLARE_UNUSED(printWhyNot);
    DECLARE_UNUSED(eps);
    return true;
  }
  void resetFromVVD(const VectorVectorDouble& tab, bool byCol = true) override;

  void normMatrix(const AMatrix& y, const AMatrixSquare& x = AMatrixSquare(), bool transpose = false);

  static MatrixSquareSymmetric* createFromVVD(const VectorVectorDouble &X);
  static MatrixSquareSymmetric* createFromVD(const VectorDouble &X);
  static MatrixSquareSymmetric* createFromTLTU(int neq,
                                               const VectorDouble &tl);
  static MatrixSquareSymmetric* createFromTriangle(int mode,
                                                   int neq,
                                                   const VectorDouble &tl);
  static MatrixSquareSymmetric* createRandomDefinitePositive(int neq, int seed = 13242);
  static MatrixSquareSymmetric* sample(const MatrixSquareSymmetric* A,
                                       const VectorInt& rowKeep,
                                       bool flagInvert = false);

  int computeEigen(bool optionPositive = true);
  int computeGeneralizedEigen(const MatrixSquareSymmetric& b, bool optionPositive = true);
  int computeGeneralizedInverse(MatrixSquareSymmetric &tabout,
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
  MatrixSquareSymmetric compress0MatLC(const MatrixDense& matLC);

private:
  int _getTriangleSize() const;
};
