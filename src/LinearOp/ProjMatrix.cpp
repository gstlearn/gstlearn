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
#include "../../include/LinearOp/ProjMatrix.hpp"
#include "Mesh/AMesh.hpp"
#include "Basic/Vector.hpp"
#include "geoslib_e.h"

ProjMatrix::ProjMatrix() 
  : IProjMatrix()
  , _nPoint(0)
  , _nApices(0)
  , _Aproj(nullptr)
{
}

ProjMatrix::ProjMatrix(const Db* db, AMesh *a_mesh, int verbose)

  : IProjMatrix()
  , _nPoint(0)
  , _nApices(0)
  , _Aproj(nullptr)
{
  if (init(db,a_mesh,verbose))
  {
    messerr("Problem in Constructor of ProjMatrix");
    return;
  }
}

ProjMatrix::ProjMatrix(int npoint, int napices, const cs *aproj)
  : IProjMatrix()
  , _nPoint(0)
  , _nApices(0)
  , _Aproj(nullptr)
{
  if (init(npoint, napices, aproj))
  {
    messerr("Problem in Constructor of ProjMatrix");
    return;
  }
}

ProjMatrix::ProjMatrix(const ProjMatrix &m)
    : IProjMatrix(m),
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
     IProjMatrix::operator =(m);
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

int ProjMatrix::init(const Db* db, AMesh *a_mesh, int verbose)
{
  _Aproj = a_mesh->getMeshToDb(db,verbose);
  if (_Aproj == (cs *) NULL) return 1;
  _nPoint = db->getActiveSampleNumber();
  _nApices = a_mesh->getNApices();
  return 0;
}

int ProjMatrix::init(int npoint, int napices, const cs *aproj)
{
  _Aproj   = cs_duplicate(aproj);
  if (_Aproj == (cs *) NULL) return 1;
  _nPoint  = npoint;
  _nApices = napices;
  return 0;
}

int ProjMatrix::init(Db* db, SPDE_Mesh* s_mesh, int verbose)
{
  MeshEStandard amesh;
  amesh.convertFromOldMesh(s_mesh, 0);
  _Aproj = db_mesh_sparse(db, &amesh, verbose);
  if (_Aproj == (cs *) NULL) return 1;
  _nPoint = db->getActiveSampleNumber();
  _nApices = s_mesh->nmesh;
  return 0;
}

int ProjMatrix::init(const Db* db,
                     SPDE_Mesh* s_mesh,
                     double radius,
                     int flag_exact,
                     int verbose)
{
  int nactive;
  int *ranks;
  _Aproj   = db_mesh_neigh(db,s_mesh,radius,flag_exact,verbose,
                           &nactive,&ranks);
  if (_Aproj == (cs *) NULL) return 1;
  if (ranks != (int *) NULL) ranks = (int *) mem_free((char *) ranks);
  _nPoint  = nactive;
  _nApices = s_mesh->nmesh;
  return 0;
}

int ProjMatrix::point2mesh(const VectorDouble& in, VectorDouble& out) const
{
  if ((int) in.size() != _nPoint)
  {
    messerr("Error in the dimension of argument 'in'(%d). It should be (%d)",
            in.size(),_nPoint);
    return 1;
  }
  if ((int) out.size() != _nApices)
  {
    messerr("Error in the dimension of argument 'out'(%d). It should be (%d)",
            out.size(),_nApices);
    return 1;
  }

  cs_tmulvec(_Aproj,in.data(),out.data());
  return 0;
}

int ProjMatrix::mesh2point(const VectorDouble& in, VectorDouble& out) const
{
  if ((int) in.size() != _nApices)
  {
    messerr("Error in the dimension of argument 'in'(%d). It should be (%d)",
            in.size(),_nApices);
    return 1;
  }
  if ((int) out.size() != _nPoint)
  {
    messerr("Error in the dimension of argument 'out'(%d). It should be (%d)",
            out.size(),_nPoint);
    return 1;
  }
  cs_mulvec(_Aproj,_nPoint,in.data(),out.data());
  return 0;
}

String ProjMatrix::toString(int level) const
{
  std::stringstream sstr;

  sstr << toStringDim("Projection Matrix",_Aproj);
  if (level < 1) return sstr.str();
  sstr << toStringRange(String(),_Aproj);
  sstr << toMatrix(String(), _Aproj);
  return sstr.str();
}
