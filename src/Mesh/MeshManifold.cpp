/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "Mesh/MeshManifold.hpp"
#include "Mesh/AMesh.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Db/Db.hpp"
#include "Basic/Vector.hpp"
#include "csparse_f.h"

MeshManifold::MeshManifold()
  : AMesh()
  , _apices()
  , _meshes()
  , _units()
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
int MeshManifold::getApex(int imesh,
                           int rank) const
{
  return (_meshes[getNApexPerMesh() * imesh + rank] - 1);
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
  (void) setExtend(extendmin,extendmax);
}

int MeshManifold::_recopy(const MeshManifold &m)
{
  _apices = m._apices;
  _meshes = m._meshes;
  _units = m._units;
  return(0);
}
