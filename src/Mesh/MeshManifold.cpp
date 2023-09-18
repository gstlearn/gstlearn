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
#include "geoslib_old_f.h"

#include "Mesh/MeshManifold.hpp"
#include "Mesh/AMesh.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixInt.hpp"
#include "Db/Db.hpp"
#include "Basic/ASerializable.hpp"

MeshManifold::MeshManifold()
  : AMesh()
  , _apices()
  , _meshes()
{
}

MeshManifold::MeshManifold(const MeshManifold &m)
  : AMesh(m)
{
  _recopy(m);
}

MeshManifold& MeshManifold::operator= (const MeshManifold &m)
{
  _recopy(m);
  return *this;
}

MeshManifold::~MeshManifold()
{
}

/****************************************************************************/
/*!
** Returns the number of Apices
**
** \returns Number of apices
**
*****************************************************************************/
int MeshManifold::getNApices() const
{
  return _apices.getNRows();
}

/****************************************************************************/
/*!
** Returns the number of Meshes
**
** \returns Number of meshes
**
*****************************************************************************/
int MeshManifold::getNMeshes() const
{
  return static_cast<int> (_meshes.size()) / getNApexPerMesh();
}

/****************************************************************************/
/*!
** Print the contents of the meshing
**
** \param[in] strfmt    Format for printout
**
*****************************************************************************/
String MeshManifold::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  sstr << toTitle(0,"Manifold Meshing characteristics");
  AMesh::toString(strfmt);
  return sstr.str();
}

/****************************************************************************/
/*!
** Returns the rank of the Apex 'rank' of the Mesh 'imesh'
**
** \returns The rank of the target  apex
**
** \param[in]  imesh    Rank of the Mesh (from 0 to _nMeshes-1))
** \param[in]  rank     Rank of the Apex within a Mesh (from 0 to _nApices-1)
**
*****************************************************************************/
int MeshManifold::getApex(int imesh, int rank) const
{
  return _meshes.getValue(imesh, rank) - 1;
}

/****************************************************************************/
/*!
** Returns the coordinate 'ic' of the Apex 'rank' of the Mesh 'imesh'
**
** \returns The coordinate of the target apex
**
** \param[in]  imesh    Rank of the Mesh (from 0 to _nMeshes-1))
** \param[in]  rank     Rank of the Apex within a Mesh (from 0 to _nApices-1)
** \param[in]  idim     Rank of the coordinate (from 0 to _ndimh-1)
**
*****************************************************************************/
double MeshManifold::getCoor(int imesh, int rank, int idim) const
{
  return _apices(getApex(imesh,rank),idim);
}

double MeshManifold::getApexCoor(int i, int idim) const
{
  return _apices(i,idim);
}

/**
 * Create a MeshManifold by loading the contents of a Neutral File
 *
 * @param neutralFilename Name of the Neutral File (MeshEStandard format)
 * @param verbose         Verbose
 */
MeshManifold* MeshManifold::createFromNF(const String& neutralFilename, bool verbose)
{
  MeshManifold* mesh = nullptr;
  std::ifstream is;
  if (_fileOpenRead(neutralFilename, is, verbose))
  {
    mesh = new MeshManifold;
    if (! mesh->deserialize(is, verbose))
    {
      delete mesh;
      mesh = nullptr;
    }
    is.close();
  }
  return mesh;
}
#ifndef SWIG
cs* MeshManifold::getMeshToDb(const Db *db, int rankZ, bool verbose) const
{
  /// TODO getMeshToDb
  DECLARE_UNUSED(db);
  DECLARE_UNUSED(rankZ);
  DECLARE_UNUSED(verbose);
  return nullptr;
}
#endif
void MeshManifold::_defineBoundingBox(void)
{
  VectorDouble extendmin;
  VectorDouble extendmax;
  double coor,mini,maxi;
  int ndim = getNDim();

  // Initializations
  extendmin.resize(ndim);
  extendmax.resize(ndim);

  // Loop on the Space dimensions
  for (int idim=0; idim<ndim; idim++)
  {
    mini =  1.e30;
    maxi = -1.e30;

    // Loop on the apices
    for (int i=0; i<getNApices(); i++)
    {
      coor = getApexCoor(i,idim);
      if (coor < mini) mini = coor;
      if (coor > maxi) maxi = coor;
    }
    extendmin[idim] = mini;
    extendmax[idim] = maxi;
  }

  // Store the Bounding Box extension
  (void) _setExtend(extendmin,extendmax);
}

int MeshManifold::_recopy(const MeshManifold &m)
{
  _apices = m._apices;
  _meshes = m._meshes;
  return(0);
}

int MeshManifold::_deserialize(std::istream& is, bool /*verbose*/)
{
  int ndim, napices, nmeshes, napexpermesh;

  bool ret = true;
  ret = ret && _recordRead<int>(is, "Space Dimension", ndim);
  ret = ret && _recordRead<int>(is, "Napices", napices);
  ret = ret && _recordRead<int>(is, "Number of Apices per Mesh", napexpermesh);
  ret = ret && _recordRead<int>(is, "Number of Meshes", nmeshes);

  VectorDouble local;
  ret = ret && _recordReadVec<double>(is, "Apices", local, ndim * napices);
  _apices = MatrixRectangular(napices, ndim);
  _apices.setValues(local);
  ret = ret && _recordReadVec<int>(is, "Meshes", _meshes.getValues(), nmeshes * napexpermesh);
  return 0;
}

int MeshManifold::_serialize(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;
  ret = ret && _recordWrite<int>(os, "Space Dimension", getNDim());
  ret = ret && _recordWrite<int>(os, "Napices", getNApices());
  ret = ret && _recordWrite<int>(os, "Number of Apices per Mesh", getNApexPerMesh());
  ret = ret && _recordWrite<int>(os, "Number of Meshes", getNMeshes());

  ret = ret && _recordWriteVec<double>(os, "Apices", _apices.getValues());
  ret = ret && _recordWriteVec<int>(os, "Meshes", _meshes.getValues());
  return 0;
}
