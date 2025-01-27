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
#include "LinearOp/IProj.hpp"
#include "Matrix/MatrixSparse.hpp"

class AMesh;
class Db;

class GSTLEARN_EXPORT ProjMatrix: public IProj, public MatrixSparse
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

  /// Interface for IProj

#ifndef SWIG
  protected:
    int _addMesh2point(const constvect inv, vect outv) const override;
    int _addPoint2mesh(const constvect inv, vect outv) const override;
#endif
  public:

  int getNApex() const override { return getNCols(); }
  int getNPoint() const override { return getNRows(); }

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
