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

#include "LinearOp/ALinearOp.hpp"
#include "gstlearn_export.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Basic/WarningMacro.hpp"
#include "Basic/VectorNumT.hpp"
#include "Matrix/AMatrix.hpp"
#include <Eigen/src/Core/Matrix.h>

#ifndef SWIG
DISABLE_WARNING_PUSH
DISABLE_WARNING_COND_EXPR_CONSTANT
DISABLE_WARNING_UNUSED_BUT_SET_VARIABLE
DISABLE_WARNING_DECLARATION_HIDE_GLOBAL
#include <Eigen/Sparse>
DISABLE_WARNING_POP
#endif

class Cholesky;
class cs;
class EOperator;

/**
 * Sparse Matrix
 * 
 * Handle a sparse matrix that can be symmetrical, square or not.
 * Storage relies either on Eigen3 Library (see opt_eigen flag) or cs code.
 * Default storage option can be set globally by using setGlobalFlagEigen
 */
class GSTLEARN_EXPORT MatrixSparse : public AMatrix, public ALinearOp {

public:
  MatrixSparse(int nrow = 0, int ncol = 0, int opt_eigen = -1);
#ifndef SWIG
  MatrixSparse(const cs* A);
#endif
  MatrixSparse(const MatrixSparse &m);
  MatrixSparse& operator= (const MatrixSparse &m);
  virtual ~MatrixSparse();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// Cloneable interface
  IMPLEMENT_CLONING(MatrixSparse)

  int getSize() const override { return getNRows();} // It assumes that the matrix is symmetric. Maybe a class MatrixSparseSymmetric would be interesting
                                                     // to inherit from ALinearOp
  bool isFlagEigen() const { return _flagEigen; }
  /// Interface for AMatrix
  /*! Returns if the current matrix is Sparse */
  bool isSparse() const override { return true; }
  /*! Returns if the matrix belongs to the MatrixSparse class (avoids dynamic_cast) */
  bool isDense() const override { return false; }

  /*! Set the value for a matrix cell */
  void setValue(int irow, int icol, double value, bool flagCheck=true) override;
  /*! Get the value from a matrix cell */
  double getValue(int row, int col, bool flagCheck=true) const override;
  /*! Modifies the contents of a matrix cell */
  void updValue(int irow,
                int icol,
                const EOperator &oper,
                double value,
                bool flagCheck = true) override;

#ifndef SWIG
  int addVecInPlace(const Eigen::Map<const Eigen::VectorXd>& xm,
                    Eigen::Map<Eigen::VectorXd>& ym) const;
  void addProdMatVecInPlaceToDest(const constvect& in,
                                  vect& out,
                                  bool transpose = false) const;
#endif
  /*! Set the contents of a Column */
  virtual void setColumn(int icol,
                         const VectorDouble &tab,
                         bool flagCheck = true) override;
  /*! Set the contents of a Row */
  virtual void setRow(int irow,
                      const VectorDouble& tab,
                      bool flagCheck=true) override;
  /*! Set the contents of the (main) Diagonal */
  virtual void setDiagonal(const VectorDouble& tab, bool flagCheck=true) override;
  /*! Set the contents of the (main) Diagonal to a constant value */
  virtual void setDiagonalToConstant(double value = 1.) override;
  /*! Transpose the matrix and return it as a copy*/
  virtual MatrixSparse* transpose() const override;
  /*! Add a value to each matrix component */
  virtual void addScalar(double v) override;
  /*! Add value to matrix diagonal */
  virtual void addScalarDiag(double v) override;
  /*! Multiply each matrix component by a value */
  virtual void prodScalar(double v) override;
  /*! Set all the values of the matrix at once */
  virtual void fill(double value) override;
  /*! Multiply the matrix row-wise */
  virtual void multiplyRow(const VectorDouble& vec) override;
  /*! Multiply the matrix column-wise */
  virtual void multiplyColumn(const VectorDouble& vec) override;
  /*! Divide the matrix row-wise */
  virtual void divideRow(const VectorDouble& vec) override;
  /*! Divide the matrix column-wise */
  virtual void divideColumn(const VectorDouble& vec) override;
  /*! Perform y = x %*% 'this' */
  virtual VectorDouble prodVecMat(const VectorDouble& x, bool transpose = false) const override;
  /*! Perform y = 'this' %*% x */
  virtual VectorDouble prodMatVec(const VectorDouble& x, bool transpose = false) const override;
  /*! Multiply matrix 'x' by matrix 'y' and store the result in 'this' */
  virtual void prodMatMatInPlace(const AMatrix *x,
                                 const AMatrix *y,
                                 bool transposeX = false,
                                 bool transposeY = false) override;

  /*! Extract the contents of the matrix */
  virtual NF_Triplet getMatrixToTriplet(int shiftRow=0, int shiftCol=0) const override;

  //// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  // Static functions
  static MatrixSparse* createFromTriplet(const NF_Triplet &NF_T,
                                         int nrow = 0,
                                         int ncol = 0,
                                         int opt_eigen = -1);
  static MatrixSparse* addMatMat(const MatrixSparse *x,
                                 const MatrixSparse *y,
                                 double cx = 1.,
                                 double cy = 1.);
  static MatrixSparse* diagVec(const VectorDouble& vec,
                               int opt_eigen = -1);
  static MatrixSparse* diagConstant(int number, double value = 1.,
                                    int opt_eigen = -1);
  static MatrixSparse* diagMat(MatrixSparse *A, int oper_choice,
                                    int opt_eigen = -1);
  static MatrixSparse* glue(const MatrixSparse *A1,
                            const MatrixSparse *A2,
                            bool flagShiftRow,
                            bool flagShiftCol);
  static void  glueInPlace(MatrixSparse *A1,
                           const MatrixSparse *A2,
                           bool flagShiftRow,
                           bool flagShiftCol);
  /// The next functions use specific definition of matrix (to avoid dynamic_cast)
  /// rather than manipulating AMatrix. They are no more generic of AMatrix
  /*! Add a matrix (multiplied by a constant) */
  virtual void addMatInPlace(const MatrixSparse& y, double cx = 1., double cy = 1.);
  /*! Product 't(A)' %*% 'M' %*% 'A' or 'A' %*% 'M' %*% 't(A)' stored in 'this'*/
  virtual void prodNormMatMatInPlace(const MatrixSparse* a,
                                     const MatrixSparse* m,
                                     bool transpose = false);
  /*! Product 't(A)' %*% ['vec'] %*% 'A' or 'A' %*% ['vec'] %*% 't(A)' stored in 'this'*/
  virtual void prodNormMatInPlace(const MatrixSparse* a,
                                  const VectorDouble& vec = VectorDouble(),
                                  bool transpose = false);

#ifndef SWIG
  /*! Returns a pointer to the Sparse storage */
  const cs* getCS() const;
  void setCS(cs* cs);
  void freeCS();
  /*! Temporary function to get the CS contents of Sparse Matrix */
  cs* getCSUnprotected() const;
#endif

  //virtual void reset(int nrows, int ncols) override; // Use base class method
  virtual void resetFromValue(int nrows, int ncols, double value) override;
  virtual void resetFromArray(int nrows, int ncols, const double* tab, bool byCol = true) override;
  virtual void resetFromVD(int nrows, int ncols, const VectorDouble &tab, bool byCol = true) override;
  virtual void resetFromVVD(const VectorVectorDouble& tab, bool byCol = true) override;

  void resetFromTriplet(const NF_Triplet& NF_T);

  /*! Dump a specific range of samples from the internal storage */
  static void dumpElements(const String& title, int ifrom, int ito);

  /*! Set all the values of the Matrix with random values */
  void fillRandom(int seed = 432432, double zeroPercent = 0);

  // Cholesky functions
  int    computeCholesky();
  int    solveCholesky(const VectorDouble& b, VectorDouble& x);

  #ifndef SWIG
    int  solveCholesky(const constvect& b, std::vector<double>& x);
    int  simulateCholesky(const constvect &b, vect &x);
    int  addVecInPlace(const constvect& x, vect& y);
  #endif
  
  int    simulateCholesky(const VectorDouble &b, VectorDouble &x);
  double computeCholeskyLogDeterminant();

  void   addValue(int row, int col, double value);

  double L1Norm() const;
  void   getStats(int *nrows, int *ncols, int *count, double *percent) const;
  int    scaleByDiag();
  int    addVecInPlace(const VectorDouble& x, VectorDouble& y);
  void   setConstant(double value);
  VectorDouble extractDiag(int oper_choice = 1) const;
  void   prodNormDiagVecInPlace(const VectorDouble &vec, int oper = 1);
  #ifndef SWIG
  const Eigen::SparseMatrix<double>& getEigenMatrix() const { return _eigenMatrix; }
  void setEigenMatrix(const Eigen::SparseMatrix<double> &eigenMatrix) { _eigenMatrix = eigenMatrix; }
  #endif
  
  MatrixSparse* extractSubmatrixByRanks(const VectorInt &rank_rows,
                                        const VectorInt &rank_cols) const;
  MatrixSparse* extractSubmatrixByColor(const VectorInt &colors,
                                        int ref_color,
                                        bool row_ok,
                                        bool col_ok);
  VectorInt colorCoding() const;
  int getNonZeros() const { return _getMatrixPhysicalSize(); }
  void gibbs(int iech, const VectorDouble& zcur, double* yk, double* sk);

#ifndef SWIG
  protected:
  virtual int _addToDest(const constvect& inv,
                          vect& outv) const override;
#endif

#ifndef SWIG
  public :
  void setDiagonal(const Eigen::Map<const Eigen::VectorXd>& tab);
#endif
protected:
  /// Interface for AMatrix
  bool _isPhysicallyPresent(int irow, int icol) const override
  {
    DECLARE_UNUSED(irow, icol);
    return true;
  }
  virtual void    _allocate() override;
  virtual void    _deallocate() override;

  virtual double& _getValueRef(int irow, int icol) override;
  virtual int     _getMatrixPhysicalSize() const override;
  virtual void    _setValueByRank(int rank, double value) override;
  virtual void    _setValues(const double* values, bool byCol) override;
  virtual double  _getValueByRank(int rank) const override;
  virtual int     _getIndexToRank(int irow,int icol) const override;
  virtual void    _transposeInPlace() override;

  virtual void    _prodMatVecInPlacePtr(const double *x, double *y, bool transpose = false) const override;
  virtual void    _prodVecMatInPlacePtr(const double *x,double *y, bool transpose = false) const override;
  virtual void    _addProdMatVecInPlaceToDestPtr(const double *x, double *y, bool transpose = false) const override;

  virtual int     _invert() override;
  virtual int     _solve(const VectorDouble& b, VectorDouble& x) const override;

  void _clear() override;
  bool _isElementPresent(int irow, int icol) const;

private:
  static bool _defineFlagEigen(int opt_eigen);
  static void _forbiddenForSparse(const String& func);
  int  _eigen_findColor(int imesh,
                        int ncolor,
                        VectorInt &colors,
                        VectorInt &temp) const;

private:
#ifndef SWIG
  cs* _csMatrix; // Classical storage for Sparse matrix
  Eigen::SparseMatrix<double> _eigenMatrix; // Eigen storage in Eigen Library (always stored Eigen::ColMajor)
#endif
  Cholesky* _factor; // Cholesky decomposition
  bool _flagEigen;
};

/*! Transform any matrix into a Sparse format */
GSTLEARN_EXPORT MatrixSparse *createFromAnyMatrix(const AMatrix* mat, int opt_eigen = -1);
GSTLEARN_EXPORT void setUpdateNonZeroValue(int status = 2);
GSTLEARN_EXPORT int getUpdateNonZeroValue();

/*! Product 't(A)' %*% 'M' %*% 'A' or 'A' %*% 'M' %*% 't(A)' */
GSTLEARN_EXPORT MatrixSparse* prodNormMatMat(const MatrixSparse* a,
                                             const MatrixSparse* m,
                                             bool transpose = false);
/*! Product 't(A)' %*% ['vec'] %*% 'A' or 'A' %*% ['vec'] %*% 't(A)' stored in 'this'*/
GSTLEARN_EXPORT MatrixSparse* prodNormMat(const MatrixSparse* a,
                                          const VectorDouble& vec = VectorDouble(),
                                          bool transpose = false);
/*! Product 'Diag(vec)' %*% 'A' %*% 'Diag(vec)' */
GSTLEARN_EXPORT MatrixSparse* prodNormDiagVec(const MatrixSparse* a,
                                              const VectorDouble& vec,
                                              int oper_choice = 1);


/// Manage global flag for EIGEN
GSTLEARN_EXPORT void setGlobalFlagEigen(bool flagEigen);
GSTLEARN_EXPORT bool isGlobalFlagEigen();


// Not exported method

#ifndef SWIG
GSTLEARN_EXPORT Eigen::SparseMatrix<double> AtMA(const Eigen::SparseMatrix<double>& A,
                                                 const Eigen::SparseMatrix<double>& M);
#endif
