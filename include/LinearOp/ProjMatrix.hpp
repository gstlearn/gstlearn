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

#include "IProjMatrix.hpp"
#include "Basic/VectorNumT.hpp"
#include "Matrix/MatrixSparse.hpp"

class AMesh;
class Db;

class GSTLEARN_EXPORT ProjMatrix: public IProjMatrix, public MatrixSparse
{
public:
  ProjMatrix();
  ProjMatrix(const Db* db,
             const AMesh* a_mesh,
             int rankZ    = -1,
             bool verbose = false);
  ProjMatrix(const ProjMatrix& m);
  ProjMatrix(const MatrixSparse* m);
  ProjMatrix& operator=(const ProjMatrix& m);
  virtual ~ProjMatrix();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// Cloneable interface
  IMPLEMENT_CLONING(ProjMatrix)

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface for IProjMatrix
  
  #ifndef SWIG
  protected:
    int _addMesh2point(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const override;
    int _addPoint2mesh(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const override;
  #endif 
  public:

  int getApexNumber() const override { return getNCols(); }
  int getPointNumber() const override { return getNRows(); }

  static ProjMatrix* create(const Db *db,
                            const AMesh *a_mesh,
                            int rankZ = -1,
                            bool verbose = false);
  void resetFromMeshAndDb(const Db* db,
                          const AMesh* a_mesh,
                          int rankZ = -1,
                          bool verbose = false);
//  int resetFromDbByNeigh(const Db *db,   // currently unused feature
//                         AMesh *amesh,
//                         double radius,
//                         int flag_exact = 0,
//                         bool verbose = false);
};
