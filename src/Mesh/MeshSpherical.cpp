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
#include "Mesh/MeshSpherical.hpp"
#include "Mesh/AMesh.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixInt.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "Db/Db.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "Tree/Ball.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpaceSN.hpp"

MeshSpherical::MeshSpherical(const MatrixDense &apices,
                             const MatrixInt &meshes)
    : AMesh(),
      _apices(apices),
      _meshes(meshes)
{
  int ndim = apices.getNCols();
  _setNDim(ndim);
}

MeshSpherical::MeshSpherical(const MeshSpherical &m) 
  : AMesh(m)
{
  _recopy(m);
}

MeshSpherical& MeshSpherical::operator= (const MeshSpherical &m)
{
  _recopy(m);
  return *this;
}

MeshSpherical::~MeshSpherical()
{
}

/****************************************************************************/
/*!
** Returns the number of Apices
**
** \returns Number of apices
**
*****************************************************************************/
int MeshSpherical::getNApices() const
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
int MeshSpherical::getNMeshes() const
{
  return static_cast<int> (_meshes.size()) / getNApexPerMesh();
}

/****************************************************************************/
/*!
** Returns the size of the Mesh 'imesh'
**
** \returns mesh dimension
**
** \param[in]  imesh    Rank of the Mesh (from 0 to _nMeshes-1))
**
*****************************************************************************/
double MeshSpherical::getMeshSize(int imesh) const
{
  return GH::geodeticTriangleSurface(getCoor(imesh, 0, 0), getCoor(imesh, 0, 1),
                                     getCoor(imesh, 1, 0), getCoor(imesh, 1, 1),
                                     getCoor(imesh, 2, 0), getCoor(imesh, 2, 1));
}

/****************************************************************************/
/*!
** Print the contents of the meshing
**
** \param[in] strfmt    Format for printout
**
*****************************************************************************/
String MeshSpherical::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  sstr << toTitle(0,"Spherical Meshing characteristics");
  sstr << AMesh::toString(strfmt);
  return sstr.str();
}

/**
 * Create a MeshSpherical by loading the contents of a Neutral File
 *
 * @param neutralFilename Name of the Neutral File (MeshEStandard format)
 * @param verbose         Verbose
 */
MeshSpherical* MeshSpherical::createFromNF(const String& neutralFilename, bool verbose)
{
  MeshSpherical* mesh = nullptr;
  std::ifstream is;
  mesh = new MeshSpherical;
  bool success = false;
  if (mesh->_fileOpenRead(neutralFilename, is, verbose))
  {
    success =  mesh->deserialize(is, verbose);
  }
  if (! success)
  {
    delete mesh;
    mesh = nullptr;
  }
  return mesh;
}

MeshSpherical* MeshSpherical::create(const MatrixDense &apices,
                                     const MatrixInt &meshes)
{
  return new MeshSpherical(apices, meshes);
}

/****************************************************************************/
/*!
** Create the meshing (from mesh information)
**
** \param[in]  ndim            Space Dimension
** \param[in]  napexpermesh    Number of apices per mesh
** \param[in]  apices          Vector of Apex information
** \param[in]  meshes          Vector of mesh indices
** \param[in]  byCol           true for Column major; false for Row Major
** \param[in]  verbose         Verbose flag
**
** \remark The argument 'byCol' concerns 'apices' and 'meshes'
**
*****************************************************************************/
int MeshSpherical::reset(int ndim,
                         int napexpermesh,
                         const VectorDouble &apices,
                         const VectorInt &meshes,
                         bool byCol,
                         bool verbose)
{
  _setNDim(ndim);
  int npoints = static_cast<int> (apices.size()) / ndim;
  int nmeshes = static_cast<int> (meshes.size()) / napexpermesh;

  // Core allocation

  _apices.reset(npoints,ndim);
  _apices.setValues(apices, byCol);
  _meshes.reset(nmeshes,napexpermesh);
  _meshes.setValues(meshes, byCol);

  // Check consistency

  _checkConsistency();

  // Define and store the Bounding Box extension

  _defineBoundingBox();

  // Optional printout

  if (verbose) messageFlush(toString());

  return(0);
}

/****************************************************************************/
/*!
**  Check if a point, defined by its coordinates, belongs to a Mesh
**  The geometry is expressed on the Sphere where ndim=2 and ndimEmbedded=3
**
** \return true if the point belongs to the Mesh; false otherwise
**
** \param[in]  coor      Vector of target coordinates
** \param[in]  corners   Vector of coordinates of mesh apices
** \param[in]  meshsize  Dimension of the mesh
** \param[in]  eps       Tolerance
**
** \param[out] weights   Array of barycentric weights (Dim: NApexPerMesh)
**
** \remarks The argument 'meshsize' is used to speed the calculations
**
*****************************************************************************/
bool MeshSpherical::_weightsInMesh(const VectorDouble& coor,
                                   const VectorVectorDouble& corners,
                                   double meshsize,
                                   VectorDouble& weights,
                                   double eps) const
{
  DECLARE_UNUSED(meshsize);
  DECLARE_UNUSED(eps);
  return GH::isInSphericalTriangleOptimized(coor.data(),
                                            corners[0].data(), corners[1].data(), corners[2].data(),
                                            weights.data());
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
int MeshSpherical::getApex(int imesh, int rank) const
{
  return _meshes.getValue(imesh,rank);
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
double MeshSpherical::getCoor(int imesh, int rank, int idim) const
{
  return _apices(getApex(imesh,rank),idim);
}

double MeshSpherical::getApexCoor(int i, int idim) const
{
  return _apices(i,idim);
}

void MeshSpherical::_getCoordOnSphere(double longitude,
                                      double latitude,
                                      VectorDouble& coords)
{
  double radius = EARTH_RADIUS;
  if (isDefaultSpaceSphere())
  {
    const ASpace* space    = getDefaultSpaceSh().get();
    const SpaceSN* spaceSn = dynamic_cast<const SpaceSN*>(space);
    if (spaceSn != nullptr) radius = spaceSn->getRadius();
  }
  GH::convertSph2Cart(longitude, latitude,
                      &coords.at(0), &coords.at(1), &coords.at(2), radius);
}

void MeshSpherical::getEmbeddedCoorPerMesh(int imesh, int ic, VectorDouble& coords) const
{
  _getCoordOnSphere(getCoor(imesh, ic, 0), getCoor(imesh, ic, 1), coords);
}

void MeshSpherical::getEmbeddedCoorPerApex(int iapex, VectorDouble& coords) const
{
  _getCoordOnSphere(getApexCoor(iapex, 0), getApexCoor(iapex, 1), coords);
}

/**
 * Calculate the Mesh of each Mesh (using approximated calculations)
 * @return
 */
VectorDouble MeshSpherical::_defineUnits(void) const
{
  int nmeshes = getNMeshes();
  VectorDouble units(nmeshes);
  for (int imesh=0; imesh<nmeshes; imesh++)
  {
    VectorVectorDouble corners = getCoordinatesPerMesh(imesh);
    units[imesh] = _getMeshUnit(corners);
  }
  return units;
}

void MeshSpherical::_defineBoundingBox(void)
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

double MeshSpherical::_closestValue(double ref, double coor, double period)
{
  double dref = ABS(coor - ref);
  double d1   = ABS(coor - period - ref);
  if (d1 < dref) return coor - period;
  return coor;
}

int MeshSpherical::_recopy(const MeshSpherical &m)
{
  _apices = m._apices;
  _meshes = m._meshes;
  AMesh::_recopy(m);
  return(0);
}

bool MeshSpherical::_deserialize(std::istream& is, bool /*verbose*/)
{
  int ndim = 0;
  int napices = 0;
  int nmeshes = 0;
  int napexpermesh = 0;

  bool ret = true;
  ret = ret && _recordRead<int>(is, "Space Dimension", ndim);
  ret = ret && _recordRead<int>(is, "Napices", napices);
  ret = ret && _recordRead<int>(is, "Number of Apices per Mesh", napexpermesh);
  ret = ret && _recordRead<int>(is, "Number of Meshes", nmeshes);

  if (ret)
  {
    VectorDouble apices_local;
    ret = ret && _recordReadVec<double>(is, "Apices", apices_local, ndim * napices);
    _apices = MatrixDense(napices, ndim);
    _apices.setValues(apices_local);
  }

  if (ret)
  {
    VectorInt meshes_local;
    ret = ret && _recordReadVec<int>(is, "Meshes", meshes_local, nmeshes * napexpermesh);
    _meshes = MatrixInt(nmeshes, napexpermesh);
    _meshes.setValues(meshes_local);
  }
  return ret;
}

bool MeshSpherical::_serialize(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;
  ret = ret && _recordWrite<int>(os, "Space Dimension", getNDim());
  ret = ret && _recordWrite<int>(os, "Napices", getNApices());
  ret = ret && _recordWrite<int>(os, "Number of Apices per Mesh", getNApexPerMesh());
  ret = ret && _recordWrite<int>(os, "Number of Meshes", getNMeshes());
  ret = ret && _recordWriteVec<double>(os, "Apices", _apices.getValues());
  ret = ret && _recordWriteVec<int>(os, "Meshes", _meshes.getValues());
  return ret;
}

/**
 * This function checks the consistency between the number of points
 * and the vertices indication
 */
void MeshSpherical::_checkConsistency() const
{
  for (int imesh = 0; imesh < getNMeshes(); imesh++)
    for (int ic = 0; ic < getNApexPerMesh(); ic++)
    {
      int apex = getApex(imesh, ic);
      if (apex < 0 || apex >= getNApices())
      {
        my_throw("Mesh indices are not compatible with the Points");
      }
    }
}

void MeshSpherical::getBarycenterInPlace(int imesh, VectorDouble& coord) const
{
  int ndimE   = getEmbeddedNDim();
  int ncorner = getNApexPerMesh();

  // Calculate the center of gravity (in the Embedded space)
  VectorVectorDouble coordE = getEmbeddedCoordinatesPerMesh(imesh);
  VectorDouble centerE(ndimE);

  for (int idimE = 0; idimE < ndimE; idimE++)
  {
    double local = 0.;
    for (int ic = 0; ic < ncorner; ic++)
      local += coordE[ic][idimE];
    centerE[idimE] = local / ncorner;
  }

  // Turn the gravity center from embedded to long/lat coordinates
  GH::convertCart2Sph(centerE[0], centerE[1], centerE[2],
                      &coord.at(0), &coord.at(1), TEST);
}
