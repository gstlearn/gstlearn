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
#include "Matrix/AMatrixDense.hpp"

/**
 * Rectangular matrices are stored by columns
 */
class GSTLEARN_EXPORT MatrixRectangular : public AMatrixDense {

public:
  MatrixRectangular(int nrow = 0, int ncol = 0);
  MatrixRectangular(const MatrixRectangular &r);
  MatrixRectangular(const AMatrix &m);
  MatrixRectangular& operator= (const MatrixRectangular &r);
	virtual ~MatrixRectangular();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

	/// Cloneable interface
  IMPLEMENT_CLONING(MatrixRectangular)

  /// Interface for AMatrix
  /*! Say if the matrix must be symmetric */
  bool mustBeSymmetric() const override { return false; }

  static MatrixRectangular* createFromVVD(const VectorVectorDouble& X);
  static MatrixRectangular* createFromVD(const VectorDouble &X,
                                         int nrow,
                                         int ncol,
                                         bool byCol = false,
                                         bool invertColumnOrder = false);
  static MatrixRectangular* glue(const AMatrix *A1,
                                 const AMatrix *A2,
                                 bool flagShiftRow,
                                 bool flagShiftCol);

  /*! Adding a Row or a Column (at the bottom or right of Rectangular Matrix) */
  void addRow(int nrow_added=1);
  void addColumn(int ncolumn_added = 1);
};
