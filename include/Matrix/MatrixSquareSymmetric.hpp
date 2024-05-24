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
 * Square Symmetric matrices are stored as Lower Triangular matrices stored by column
 */
class GSTLEARN_EXPORT MatrixSquareSymmetric : public AMatrixSquare {

public:
  MatrixSquareSymmetric(int nrow = 0, int opt_eigen = -1);
  MatrixSquareSymmetric(const MatrixSquareSymmetric &m);
  MatrixSquareSymmetric(const AMatrix &m);
  MatrixSquareSymmetric& operator= (const MatrixSquareSymmetric &r);
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

  /// Interface for AMatrixDense
  void    _setValue(int irow, int icol, double value) override;
  double  _getValue(int irow, int icol) const override;
  void    _updValue(int irow, int icol, const EOperator& oper, double value) override;
  void    _allocate_() override;

  void normMatrix(const AMatrix& y, const AMatrixSquare& x = AMatrixSquare(), bool transpose = false);

  static MatrixSquareSymmetric* createFromVVD(const VectorVectorDouble &X, int opt_eigen = -1);
  static MatrixSquareSymmetric* createFromVD(const VectorDouble &X,
                                             int nrow,
                                             int opt_eigen = -1);
  static MatrixSquareSymmetric* createFromTLTU(int neq,
                                               const VectorDouble &tl,
                                               int opt_eigen = -1);
  static MatrixSquareSymmetric* createFromTriangle(int mode,
                                                   int neq,
                                                   const VectorDouble &tl,
                                                   int opt_eigen = -1);

  int computeEigen(bool optionPositive = true);
  int computeGeneralizedEigen(const MatrixSquareSymmetric& b, bool optionPositive = true);
  int computeGeneralizedInverse(MatrixSquareSymmetric &tabout,
                                double maxicond = 1.e20,
                                double eps = EPSILON20);
  bool isDefinitePositive();  int minimizeWithConstraintsInPlace(const VectorDouble& gmat,
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
  VectorDouble getCholeskyXL() const;
  MatrixRectangular productCholeskyInPlace(int mode,
                                           int neq,
                                           int nrhs,
                                           const VectorDouble &tl,
                                           const MatrixRectangular &a);
  MatrixSquareSymmetric normCholeskyInPlace(int mode,
                                            int neq,
                                            const VectorDouble &tl,
                                            const MatrixSquareSymmetric &a);
  double computeCholeskyLogDeterminant() const;

private:
  /// Interface for AMatrix
  virtual bool    _isCompatible(const AMatrix& m) const override
  {
    return (isSameSize(m) && m.isSymmetric());
  }
  virtual int     _getMatrixPhysicalSize_() const override;
  virtual double& _getValueRef_(int irow, int icol) override;

  virtual int     _getIndexToRank_(int irow,int icol) const override;

  virtual double  _getValueByRank_(int irank) const override;
  virtual void    _setValueByRank_(int irank, double value) override;
  virtual void    _transposeInPlace_() override { return ; } // Nothing to do
  virtual void    _prodMatVecInPlacePtr_(const double *x,double *y, bool transpose = false) const override;
  virtual void    _prodVecMatInPlacePtr_(const double *x,double *y, bool transpose = false) const override;

  virtual bool    _isPhysicallyPresent(int irow, int icol) const override;
  virtual void    _setValues(const double* values, bool byCol = true) override;
  virtual int     _invert() override;
  virtual int     _solve(const VectorDouble& b, VectorDouble& x) const override;

private:
  // The subsequent methods rely on the specific local storage ('squareSymMatrix')
  void    _recopy(const MatrixSquareSymmetric& r);

  // Local functions (old style algebra)
  int _matrix_geigen(const double *a,
                     const double *b,
                     int neq,
                     double *value,
                     double *vector) const;
  void _matrix_triangular_product(int neq,
                                  int mode,
                                  const double *al,
                                  const double *b,
                                  double *x) const;
  int _matrix_solve(VectorDouble &at,
                    VectorDouble &b,
                    VectorDouble &x,
                    int neq,
                    int nrhs,
                    double eps = EPSILON20) const;
  int  _matrix_invert_triangle(int neq, double *tl);
  void _matrix_tri2sq(int neq, const double *tl, double *a);
  void _matrix_sq2tri(int mode, int neq, const double *a, double *tl);

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
  int _constraintsConcatenateMat(int nae,
                                 int nai,
                                 int neq,
                                 const VectorInt &active,
                                 const MatrixRectangular &tabemat,
                                 const MatrixRectangular &tabimat,
                                 MatrixRectangular &tabout);
  int _constraintsConcatenateVD(int nae,
                                int nai,
                                const VectorInt &active,
                                const VectorDouble &tabemat,
                                const VectorDouble &tabimat,
                                VectorDouble &tabout);
  int _constraintsCount(int nai, VectorInt& active);
  bool _checkCholeskyAlreadyPerformed(int status) const;
  int _terminateEigen(const VectorDouble &eigenValues,
                      const VectorDouble &eigenVectors,
                      bool optionPositive = true,
                      bool changeOrder = false);


private:
  VectorDouble _squareSymMatrix; // Classical storage
  bool _flagCholeskyDecompose;
  bool _flagCholeskyInverse;
  VectorDouble _tl; // Lower triangular matrix (after Cholesky decomposition)
  VectorDouble _xl; // Lower triangular matrix (inverse of _tl)

  Eigen::LLT<Eigen::MatrixXd> _factor; // Cholesky decomposition (Eigen format)
};
