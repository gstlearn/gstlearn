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

namespace gstlrn
{
ProjMatrix::ProjMatrix()
  : MatrixSparse()
{
}

ProjMatrix::ProjMatrix(const Db* db,
                       const AMesh* a_mesh,
                       Id rankZ,
                       bool verbose)
  : MatrixSparse()
{
  resetFromMeshAndDb(db, a_mesh, rankZ, verbose);
}

ProjMatrix::ProjMatrix(const ProjMatrix& m)
  : IProj(m)
  , MatrixSparse(m)
{
}

ProjMatrix::ProjMatrix(const MatrixSparse& m)
  : MatrixSparse(m)
{
}

ProjMatrix& ProjMatrix::operator=(const ProjMatrix& m)
{
  if (this != &m)
  {
    MatrixSparse::operator=(m);
  }
  return *this;
}

ProjMatrix::~ProjMatrix()
{
}

ProjMatrix* ProjMatrix::create(const Db* db, const AMesh* a_mesh, Id rankZ, bool verbose)
{
  return new ProjMatrix(db, a_mesh, rankZ, verbose);
}

void ProjMatrix::resetFromMeshAndDb(const Db* db, const AMesh* a_mesh, Id rankZ, bool verbose)
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

Id ProjMatrix::_addMesh2point(const constvect inv, vect outv) const
{
  if (static_cast<Id>(inv.size()) != getNApex())
  {
    messerr("_addMesh2point: Error in the dimension of argument 'inv'(%d). It should be (%d)",
            inv.size(), getNApex());
    return 1;
  }
  if (static_cast<Id>(outv.size()) != getNPoint())
  {
    messerr("_addMesh2point: Error in the dimension of argument 'outv'(%d). It should be (%d)",
            outv.size(), getNPoint());
    return 1;
  }

  addProdMatVecInPlaceC(inv, outv, false);
  return 0;
}

Id ProjMatrix::_addPoint2mesh(const constvect inv, vect outv) const
{
  if (static_cast<Id>(inv.size()) != getNPoint())
  {
    messerr("_addPoint2mesh: Error in the dimension of argument 'inv'(%d). It should be (%d)",
            inv.size(), getNPoint());
    return 1;
  }
  if (static_cast<Id>(outv.size()) != getNApex())
  {
    messerr("_addPoint2mesh: Error in the dimension of argument 'outv'(%d). It should be (%d)",
            outv.size(), getNApex());
    return 1;
  }

  addProdMatVecInPlaceC(inv, outv, true);
  return 0;
}

String ProjMatrix::toString(const AStringFormat* strfmt) const
{
  return MatrixSparse::toString(strfmt);
}

void ProjMatrix::dumpVerticesUsed(Id npmax) const
{
  mestitle(1, "Vertices used in the projection matrix");
  if (npmax > 0) message("(Display is limited to %d samples)\n", npmax);
  Id np = (npmax > 0) ? npmax : getNPoint();
  for (Id ip = 0; ip < np; ip++)
  {
    message("Sample %3d: ", ip);
    for (Id ic = 0, nc = getNApex(); ic < nc; ic++)
    {
      if (getValue(ip, ic) > 0)
        message(" %3d [%5.2lf]", ic, getValue(ip, ic));
    }
    message("\n");
  }
}
} // namespace gstlrn