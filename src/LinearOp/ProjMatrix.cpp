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
#include "Mesh/AMesh.hpp"
#include "Db/Db.hpp"
#include "Mesh/MeshEStandard.hpp"
#include "geoslib_old_f.h"
#include "Matrix/LinkMatrixSparse.hpp"

// External library /// TODO : Dependency to csparse to be removed
#include "csparse_d.h"
#include "csparse_f.h"

ProjMatrix::ProjMatrix() 
  : AStringable()
  , _nPoint(0)
  , _nApices(0)
  , _AprojCS(nullptr)
{
}

ProjMatrix::ProjMatrix(const Db *db, const AMesh *a_mesh, int rankZ, bool verbose)
  : AStringable()
  , _nPoint(0)
  , _nApices(0)
  , _AprojCS(nullptr)
{
  if (resetFromDb(db,a_mesh,rankZ, verbose))
  {
    messerr("Problem in Constructor of ProjMatrix");
    return;
  }
}

#ifndef SWIG
ProjMatrix::ProjMatrix(int npoint, int napices, MatrixSparse *aproj)
  : AStringable()
  , _nPoint(0)
  , _nApices(0)
  , _AprojCS(nullptr)
{
  if (resetFromPoints(npoint, napices, aproj))
  {
    messerr("Problem in Constructor of ProjMatrix");
    return;
  }
}
#endif

ProjMatrix::ProjMatrix(const ProjMatrix &m)
    : AStringable(m),
      _nPoint(m._nPoint),
      _nApices(m._nApices),
      _AprojCS(nullptr)
{
  if (_AprojCS != nullptr) delete _AprojCS;
  _AprojCS = m._AprojCS->clone();
}

ProjMatrix& ProjMatrix::operator= (const ProjMatrix &m)
{
  if (this != &m)
   {
     AStringable::operator =(m);
     _nPoint = m._nPoint;
     _nApices = m._nApices;
     _AprojCS = m._AprojCS->clone();
   }
  return *this;
}

ProjMatrix::~ProjMatrix() 
{
  delete _AprojCS;
}

ProjMatrix* ProjMatrix::create(const Db* db, const AMesh *a_mesh, int rankZ, bool verbose)
{
  return new ProjMatrix(db,a_mesh,rankZ, verbose);
}

int ProjMatrix::resetFromDb(const Db* db, const AMesh *a_mesh, int rankZ, bool verbose)
{
  if (db != nullptr)
  {
    _AprojCS = a_mesh->getMeshToDb(db, rankZ, verbose);
    if (_AprojCS == nullptr) return 1;
    _nPoint = _AprojCS->getNRows();
    _nApices = a_mesh->getNApices();
  }
  else
  {
    _AprojCS = nullptr;
    _nPoint = 0;
    _nApices = a_mesh->getNApices();
  }
  return 0;
}

#ifndef SWIG
int ProjMatrix::resetFromPoints(int npoint, int napices, MatrixSparse *aproj)
{
  _AprojCS = aproj->clone();
  _nPoint  = npoint;
  _nApices = napices;
  return 0;
}
#endif

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
int ProjMatrix::resetFromDbByNeigh(const Db *db,
                                   AMesh *amesh,
                                   double radius,
                                   int flag_exact,
                                   bool verbose)
{
  int nactive;
  int *ranks;
  _AprojCS = db_mesh_neigh(db, amesh, radius, flag_exact, verbose, &nactive, &ranks);
  if (_AprojCS == nullptr) return 1;
  if (ranks != nullptr) ranks = (int *) mem_free((char *) ranks);
  _nPoint  = nactive;
  _nApices = amesh->getNMeshes();
  return 0;
}

int ProjMatrix::point2mesh(const VectorDouble& inv, VectorDouble& outv) const
{
  if ((int) inv.size() != _nPoint)
  {
    messerr("Point2mesh: Error in the dimension of argument 'inv'(%d). It should be (%d)",
            inv.size(),_nPoint);
    return 1;
  }
  if ((int) outv.size() != _nApices)
  {
    messerr("Point2mesh: Error in the dimension of argument 'outv'(%d). It should be (%d)",
            outv.size(),_nApices);
    return 1;
  }

  matCS_tmulvec(_AprojCS,(int) outv.size(),inv.data(),outv.data());
  return 0;
}

int ProjMatrix::mesh2point(const VectorDouble& inv, VectorDouble& outv) const
{
  if ((int) inv.size() != _nApices)
  {
    messerr("mesh2point: Error in the dimension of argument 'inv'(%d). It should be (%d)",
            inv.size(),_nApices);
    return 1;
  }
  if ((int) outv.size() != _nPoint)
  {
    messerr("mesh2point: Error in the dimension of argument 'outv'(%d). It should be (%d)",
            outv.size(),_nPoint);
    return 1;
  }

  matCS_mulvec(_AprojCS,_nPoint,inv.data(),outv.data());
  return 0;
}

String ProjMatrix::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  sstr << toStringDim("Projection Matrix",_AprojCS->getCS());

  AStringFormat sf;
  if (strfmt != nullptr) sf = *strfmt;
  if (sf.getLevel() > 0)
  {
    sstr << toStringRange(String(),_AprojCS->getCS());
    sstr << toMatrix(String(), _AprojCS->getCS());
  }
  return sstr.str();
}

Triplet ProjMatrix::getAprojToTriplet(bool flag_from_1) const
{
  return getAproj()->getSparseToTriplet(flag_from_1);
}
