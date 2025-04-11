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
#include "Matrix/LinkMatrixSparse.hpp"

ProjMatrix::ProjMatrix()
  : MatrixSparse()
{
}

ProjMatrix::ProjMatrix(const Db* db,
                       const AMesh* a_mesh,
                       int rankZ,
                       bool verbose)
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
    messerr("ProjMatrix::resetFromMeshAndDb: Mesh cannot be null. Nothing is done");
    return;
  }
  if (db != nullptr)
  {
    a_mesh->resetProjFromDb(this, db, rankZ, verbose);
  }
  else
  {
    _setNRows(0);
    _setNCols(a_mesh->getNApices());
  }
}

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

/* int ProjMatrix::point2mesh(const VectorDouble& inv, VectorDouble& outv) const
{
  if ((int) inv.size() != getNPoint())
  {
    messerr("point2mesh: Error in the dimension of argument 'inv'(%d). It should be (%d)",
            inv.size(),getNPoint());
    return 1;
  }
  if ((int) outv.size() != getNApex())
  {
    messerr("point2mesh: Error in the dimension of argument 'outv'(%d). It should be (%d)",
            outv.size(),getNApex());
    return 1;
  }

  prodMatVecInPlace(inv, outv, true);
  return 0;
}
 */
int ProjMatrix::_addMesh2point(const constvect inv, vect outv) const
{
  if ((int) inv.size() != getNApex())
  {
    messerr("mesh2point: Error in the dimension of argument 'inv'(%d). It should be (%d)",
            inv.size(),getNApex());
    return 1;
  }
  if ((int) outv.size() != getNPoint())
  {
    messerr("mesh2point: Error in the dimension of argument 'outv'(%d). It should be (%d)",
            outv.size(),getNPoint());
    return 1;
  }

  addProdMatVecInPlaceToDest(inv, outv, false);
  return 0;
}

int ProjMatrix::_addPoint2mesh(const constvect inv, vect outv) const
{
  if ((int) inv.size() != getNPoint())
  {
    messerr("point2mesh: Error in the dimension of argument 'inv'(%d). It should be (%d)",
            inv.size(),getNPoint());
    return 1;
  }
  if ((int) outv.size() != getNApex())
  {
    messerr("point2mesh: Error in the dimension of argument 'outv'(%d). It should be (%d)",
            outv.size(),getNApex());
    return 1;
  }

  addProdMatVecInPlaceToDest(inv, outv, true);
  return 0;
}

/* int ProjMatrix::mesh2point(const VectorDouble& inv, VectorDouble& outv) const
{
  if ((int) inv.size() != getNApex())
  {
    messerr("mesh2point: Error in the dimension of argument 'inv'(%d). It should be (%d)",
            inv.size(),getNApex());
    return 1;
  }
  if ((int) outv.size() != getNPoint())
  {
    messerr("mesh2point: Error in the dimension of argument 'outv'(%d). It should be (%d)",
            outv.size(),getNPoint());
    return 1;
  }

  prodMatVecInPlace(inv, outv, false);
  return 0;
}
 */
String ProjMatrix::toString(const AStringFormat* strfmt) const
{
  return MatrixSparse::toString(strfmt);
}

void ProjMatrix::dumpVerticesUsed(int npmax) const
{
  mestitle(1, "Vertices used in the projection matrix");
  if (npmax > 0) message("(Display is limited to %d samples)\n", npmax);
  int np = (npmax > 0) ? npmax : getNPoint();
  for (int ip = 0; ip < np; ip++)
  {
    message("Sample %3d: ", ip);
    for (int ic = 0, nc = getNApex(); ic < nc; ic++)
    {
      if (getValue(ip, ic) > 0)
        message(" %3d [%5.2lf]", ic, getValue(ip,ic));
    }
    message("\n");
  }
}