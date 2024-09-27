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
  /*! Say if the matrix must be symmetric */
  bool mustBeSymmetric() const final { return true; }
  /// Is the matrix symmetrical ?
  bool isSymmetric(bool printWhyNot = false, double eps = EPSILON10) const final
  {
    DECLARE_UNUSED(printWhyNot);
    DECLARE_UNUSED(eps);
    return true;
  }

  void normMatrix(const AMatrix& y, const AMatrixSquare& x = AMatrixSquare(), bool transpose = false);

  static MatrixSquareSymmetric* createFromVVD(const VectorVectorDouble &X);
  static MatrixSquareSymmetric* createFromVD(const VectorDouble &X,
                                             int nrow);
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
                                     const MatrixRectangular& aemat,
                                     const VectorDouble& bemat,
                                     const MatrixRectangular& aimat,
                                     const VectorDouble& bimat,
                                     VectorDouble& xmat);

  // Next methods regards the Cholesky decomposition. They also focus on the specific storage mode
  // used for symmetric matrices, i.e. the Cholesky decomposition, giving room to the upper or lower
  // triangular storage.
  // This is temporarily ensured as a VectorDouble handelde within this class. It should probably
  // become a sperate class in the future.
  int getTriangleSize() const;
  int computeCholesky();
  int invertCholesky();
  int solveCholeskyMat(const MatrixRectangular& b, MatrixRectangular& x);
  int solveCholesky(const VectorDouble& b, VectorDouble& x);
  
  VectorDouble getCholeskyTL() const;
  double getCholeskyTL(int i, int j) const;
  double getCholeskyTL(int iad) const;
  VectorDouble getCholeskyXL() const;
  double getCholeskyXL(int i, int j) const;
  static MatrixRectangular productCholeskyInPlace(int mode,
                                                  int neq,
                                                  int nrhs,
                                                  const VectorDouble &tl,
                                                  const MatrixRectangular &a);
  static MatrixSquareSymmetric normCholeskyInPlace(int mode,
                                                   int neq,
                                                   const VectorDouble &tl,
                                                   const MatrixSquareSymmetric &a);
  double computeCholeskyLogDeterminant() const;
  
  virtual bool    _isPhysicallyPresent(int irow, int icol) const override;
  virtual void    _setValues(const double* values, bool byCol = true) override;
  virtual int     _invert() override;

  void    _recopy(const MatrixSquareSymmetric& r);

  // Local functions (old style algebra)
  int _matrix_qo(const VectorDouble& gmat, VectorDouble& xmat);
  int _matrix_qoc(bool flag_invert,
                  const VectorDouble& gmat,
                  int na,
                  const MatrixRectangular& amat,
                  const VectorDouble& bmat,
                  VectorDouble& xmat,
                  VectorDouble& lambda);
  int _constraintsError(const VectorInt& active,
                        const MatrixRectangular& aimat,
                        const VectorDouble& bimat,
                        const VectorDouble& xmat,
                        VectorDouble& vmat,
                        VectorInt& flag);
  static int _constraintsConcatenateMat(int nae,
                                        int nai,
                                        int neq,
                                        const VectorInt& active,
                                        const MatrixRectangular& tabemat,
                                        const MatrixRectangular& tabimat,
                                        MatrixRectangular& tabout);
  static int _constraintsConcatenateVD(int nae,
                                       int nai,
                                       const VectorInt& active,
                                       const VectorDouble& tabemat,
                                       const VectorDouble& tabimat,
                                       VectorDouble& tabout);
  static int _constraintsCount(int nai, VectorInt& active);
  bool _checkCholeskyAlreadyPerformed(int status) const;
  int _terminateEigen(const VectorDouble &eigenValues,
                      const VectorDouble &eigenVectors,
                      bool optionPositive = true,
                      bool changeOrder = false);

private:
  bool _flagCholeskyDecompose;
  bool _flagCholeskyInverse;
  VectorDouble _tl; // Lower triangular matrix (after Cholesky decomposition)
  VectorDouble _xl; // Lower triangular matrix (inverse of _tl)
  Eigen::LLT<Eigen::MatrixXd> _factor; // Cholesky decomposition (Eigen format)
};
