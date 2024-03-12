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
#include "LinearOp/ProjMatrix.hpp"
#include "Db/Db.hpp"
#include "Mesh/AMesh.hpp"
#include "Mesh/MeshEStandard.hpp"
#include "Matrix/NF_Triplet.hpp"
#include "Matrix/LinkMatrixSparse.hpp"
#include "Basic/AException.hpp"
#include "geoslib_old_f.h"

ProjMatrix::ProjMatrix() 
    : MatrixSparse()
{
}

ProjMatrix::ProjMatrix(const Db* db, const AMesh* a_mesh, int rankZ, bool verbose)
    : MatrixSparse()
{
  resetFromMeshAndDb(db, a_mesh, rankZ, verbose);
}

ProjMatrix::ProjMatrix(const ProjMatrix& m)
    : MatrixSparse(m)
{
}

ProjMatrix::ProjMatrix(const MatrixSparse* m)
    : MatrixSparse(*m)
{
}

ProjMatrix& ProjMatrix::operator= (const ProjMatrix &m)
{
  if (this != &m)
  {
    MatrixSparse::operator =(m);
  }
  return *this;
}

ProjMatrix::~ProjMatrix() 
{
}

ProjMatrix* ProjMatrix::create(const Db* db, const AMesh* a_mesh, int rankZ, bool verbose)
{
  return new ProjMatrix(db, a_mesh, rankZ, verbose);
}

void ProjMatrix::resetFromMeshAndDb(const Db* db, const AMesh* a_mesh, int rankZ, bool verbose)
{
  if (a_mesh == nullptr)
  {
    my_throw("ProjMatrix::resetFromMeshAndDb: Mesh cannot be null");
  }
  if (db != nullptr)
  {
    a_mesh->resetProjMatrix(this, db, rankZ, verbose);
  }
  else
  {
    _setNRows(0);
    _setNCols(a_mesh->getNApices());
  }

}

/**
 * Returns the projection matrix of a set of points (contained in a Db) onto a meshing
 * @param db Db structure
 * @param amesh Meshing structure
 * @param radius Neighborhood radius
 * @param flag_exact Type of test
 * @param verbose Verbose flag
 *
 * @remarks When flag_exact is TRUE, for each active sample of Db, a vertex
 * @remarks of the mesh is active as soon as it lies within the vicinity
 * @remarks of the sample.
 * @remarks If flag_exact is FALSE, all vertices of a mesh are considered as
 * @remarks active as soon as the mesh intersects the ball around a sample.
 * @remarks The vicinity is defined as any point located at a distance
 * @remarks from the sample smaller than 'radius'. The distance is calculated
 * @remarks as the Euclidean distance over the space whose dimension is
 * @remarks if the smallest value between the Db et Mesh space dimensions.
 * @return
 */
//int ProjMatrix::resetFromDbByNeigh(const Db *db,
//                                   AMesh *amesh,
//                                   double radius,
//                                   int flag_exact,
//                                   bool verbose)
//{
//  int nactive;
//  int *ranks;
//  if (_AprojCS != nullptr) delete _AprojCS;
//  _AprojCS = db_mesh_neigh(db, amesh, radius, flag_exact, verbose, &nactive, &ranks);
//  if (_AprojCS == nullptr) return 1;
//  if (ranks != nullptr) ranks = (int *) mem_free((char *) ranks);
//  _nPoint  = nactive;
//  _nApices = amesh->getNMeshes();
//  return 0;
//}

int ProjMatrix::point2mesh(const VectorDouble& inv, VectorDouble& outv) const
{
  if ((int) inv.size() != getPointNumber())
  {
    messerr("point2mesh: Error in the dimension of argument 'inv'(%d). It should be (%d)",
            inv.size(),getPointNumber());
    return 1;
  }
  if ((int) outv.size() != getApexNumber())
  {
    messerr("point2mesh: Error in the dimension of argument 'outv'(%d). It should be (%d)",
            outv.size(),getApexNumber());
    return 1;
  }

  prodMatVecInPlace(inv, outv, true);
  return 0;
}

int ProjMatrix::mesh2point(const VectorDouble& inv, VectorDouble& outv) const
{
  if ((int) inv.size() != getApexNumber())
  {
    messerr("mesh2point: Error in the dimension of argument 'inv'(%d). It should be (%d)",
            inv.size(),getApexNumber());
    return 1;
  }
  if ((int) outv.size() != getPointNumber())
  {
    messerr("mesh2point: Error in the dimension of argument 'outv'(%d). It should be (%d)",
            outv.size(),getPointNumber());
    return 1;
  }

  prodMatVecInPlace(inv, outv, false);
  return 0;
}

String ProjMatrix::toString(const AStringFormat* strfmt) const
{
  return MatrixSparse::toString(strfmt);
}
