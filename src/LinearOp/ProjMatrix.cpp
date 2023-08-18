/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
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
  , _Aproj(nullptr)
{
}

ProjMatrix::ProjMatrix(const Db* db,const  AMesh *a_mesh, int verbose)
  : AStringable()
  , _nPoint(0)
  , _nApices(0)
  , _Aproj(nullptr)
{
  if (resetFromDb(db,a_mesh,verbose))
  {
    messerr("Problem in Constructor of ProjMatrix");
    return;
  }
}

#ifndef SWIG
ProjMatrix::ProjMatrix(int npoint, int napices, const cs *aproj)
  : AStringable()
  , _nPoint(0)
  , _nApices(0)
  , _Aproj(nullptr)
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
      _Aproj(nullptr)
{
  if (_Aproj != nullptr) _Aproj = cs_spfree(_Aproj);
  _Aproj = cs_duplicate(m._Aproj);
}

ProjMatrix& ProjMatrix::operator= (const ProjMatrix &m)
{
  if (this != &m)
   {
     AStringable::operator =(m);
     _nPoint = m._nPoint;
     _nApices = m._nApices;
     _Aproj = cs_duplicate(m._Aproj);
   }
  return *this;
}

ProjMatrix::~ProjMatrix() 
{
  _Aproj = cs_spfree(_Aproj);
}

ProjMatrix* ProjMatrix::create(const Db* db, const AMesh *a_mesh, int verbose)
{
  return new ProjMatrix(db,a_mesh,verbose);
}

int ProjMatrix::resetFromDb(const Db* db, const AMesh *a_mesh, int verbose)
{
  if (db != nullptr)
  {
    _Aproj = a_mesh->getMeshToDb(db,verbose);
    if (_Aproj == nullptr) return 1;
    _nPoint = db->getSampleNumber(true);
    _nApices = a_mesh->getNApices();
  }
  else
  {
    _Aproj = nullptr;
    _nPoint = 0;
    _nApices = a_mesh->getNApices();
  }
  return 0;
}

#ifndef SWIG
int ProjMatrix::resetFromPoints(int npoint, int napices, const cs *aproj)
{
  _Aproj = cs_duplicate(aproj);
  if (_Aproj == nullptr) return 1;
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
                                   int verbose)
{
  int nactive;
  int *ranks;
  _Aproj = db_mesh_neigh(db, amesh, radius, flag_exact, verbose, &nactive,
                         &ranks);
  if (_Aproj == nullptr) return 1;
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

  cs_tmulvec(_Aproj,(int) outv.size(),inv.data(),outv.data());
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

  cs_mulvec(_Aproj,_nPoint,inv.data(),outv.data());
  return 0;
}

String ProjMatrix::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  sstr << toStringDim("Projection Matrix",_Aproj);

  AStringFormat sf;
  if (strfmt != nullptr) sf = *strfmt;
  if (sf.getLevel() > 0)
  {
    sstr << toStringRange(String(),_Aproj);
    sstr << toMatrix(String(), _Aproj);
  }
  return sstr.str();
}

Triplet ProjMatrix::getAprojToTriplet(bool flag_from_1) const
{
  return csToTriplet(getAproj(), flag_from_1);
}
