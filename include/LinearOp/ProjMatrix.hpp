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
#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Matrix/LinkMatrixSparse.hpp"

class cs; /// TODO : Dependency to csparse to be removed
class AMesh;
class Db;

class GSTLEARN_EXPORT ProjMatrix: public IProjMatrix, public AStringable
{
public:
  ProjMatrix();
  ProjMatrix(const Db* db, const AMesh *a_mesh, int rankZ = -1, bool verbose = false);
#ifndef SWIG
  ProjMatrix(int npoint, int napices, const cs *aproj);
#endif
  ProjMatrix(const ProjMatrix &m);
  ProjMatrix& operator= (const ProjMatrix &m);
  virtual ~ProjMatrix();

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface for IProjMatrix
  int point2mesh(const VectorDouble& inv, VectorDouble& outv) const override;
  int mesh2point(const VectorDouble& inv, VectorDouble& outv) const override;
  int getApexNumber() const override { return _nApices; }
  int getPointNumber() const override { return _nPoint; }

  static ProjMatrix* create(const Db *db,
                            const AMesh *a_mesh,
                            int rankZ = -1,
                            bool verbose = false);
  int resetFromDb(const Db *db,
                  const AMesh *a_mesh,
                  int rankZ = -1,
                  bool verbose = false);
#ifndef SWIG
  int resetFromPoints(int npoint, int napices, const cs *aproj);
#endif
  int resetFromDbByNeigh(const Db *db,
                         AMesh *amesh,
                         double radius,
                         int flag_exact = 0,
                         bool verbose = false);

#ifndef SWIG
  const cs* getAproj() const { return _Aproj; }
#endif
  Triplet getAprojToTriplet(bool flag_from_1 = false) const;

private:
  int  _nPoint;
  int  _nApices;
  cs*  _Aproj;
};
