/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* Created on: 9 avr. 2019 by N. Desassis                                     */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_old_f.h"
#include "LinearOp/ProjMatrix.hpp"
#include "Mesh/AMesh.hpp"
#include "Basic/Vector.hpp"
#include "Db/Db.hpp"
#include "Mesh/MeshEStandard.hpp"
#include "csparse_f.h"

ProjMatrix::ProjMatrix() 
  : AStringable()
  , _nPoint(0)
  , _nApices(0)
  , _Aproj(nullptr)
{
}

ProjMatrix::ProjMatrix(const Db* db, AMesh *a_mesh, int verbose)
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

int ProjMatrix::resetFromDb(const Db* db, AMesh *a_mesh, int verbose)
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

int ProjMatrix::resetFromPoints(int npoint, int napices, const cs *aproj)
{
  _Aproj = cs_duplicate(aproj);
  if (_Aproj == nullptr) return 1;
  _nPoint  = npoint;
  _nApices = napices;
  return 0;
}

int ProjMatrix::resetFromDbOldStyle(Db* db, SPDE_Mesh* s_mesh, int verbose)
{
  MeshEStandard amesh;
  amesh.convertFromOldMesh(s_mesh, 0);
  _Aproj = db_mesh_sparse(db, &amesh, verbose);
  if (_Aproj == nullptr) return 1;
  _nPoint = db->getSampleNumber(true);
  _nApices = s_mesh->nmesh;
  return 0;
}

/**
 * Returns the projection matrix of a set of points (contained in a Db) onto a meshing
 * @param db Db structure
 * @param s_mesh Meshing structure
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
int ProjMatrix::resetFromDbByNeighOldStyle(const Db* db,
                                            SPDE_Mesh* s_mesh,
                                            double radius,
                                            int flag_exact,
                                            int verbose)
{
  int nactive;
  int *ranks;
  _Aproj   = db_mesh_neigh(db,s_mesh,radius,flag_exact,verbose,
                           &nactive,&ranks);
  if (_Aproj == nullptr) return 1;
  if (ranks != nullptr) ranks = (int *) mem_free((char *) ranks);
  _nPoint  = nactive;
  _nApices = s_mesh->nmesh;
  return 0;
}

int ProjMatrix::point2mesh(const VectorDouble& in, VectorDouble& out) const
{
  if ((int) in.size() != _nPoint)
  {
    messerr("Point2mesh: Error in the dimension of argument 'in'(%d). It should be (%d)",
            in.size(),_nPoint);
    return 1;
  }
  if ((int) out.size() != _nApices)
  {
    messerr("Point2mesh: Error in the dimension of argument 'out'(%d). It should be (%d)",
            out.size(),_nApices);
    return 1;
  }

  cs_tmulvec(_Aproj,(int) out.size(),in.data(),out.data());
  return 0;
}

int ProjMatrix::mesh2point(const VectorDouble& in, VectorDouble& out) const
{
  if ((int) in.size() != _nApices)
  {
    messerr("Mesh2Point: Error in the dimension of argument 'in'(%d). It should be (%d)",
            in.size(),_nApices);
    return 1;
  }
  if ((int) out.size() != _nPoint)
  {
    messerr("Mesh2Point: Error in the dimension of argument 'out'(%d). It should be (%d)",
            out.size(),_nPoint);
    return 1;
  }
  cs_mulvec(_Aproj,_nPoint,in.data(),out.data());
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
