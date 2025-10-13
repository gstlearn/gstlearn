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
#include "Basic/SerializeHDF5.hpp"
#include "Db/Db.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixInt.hpp"
#include "Mesh/AMesh.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpaceSN.hpp"
#include "Tree/Ball.hpp"

namespace gstlrn
{
MeshSpherical::MeshSpherical(const MatrixDense& apices,
                             const MatrixInt& meshes)
  : AMesh()
  , _apices(apices)
  , _meshes(meshes)
{
  auto ndim = apices.getNCols();
  _setNDim(ndim);
}

MeshSpherical::MeshSpherical(const MeshSpherical& m)
  : AMesh(m)
{
  _recopy(m);
}

MeshSpherical& MeshSpherical::operator=(const MeshSpherical& m)
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
Id MeshSpherical::getNApices() const
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
Id MeshSpherical::getNMeshes() const
{
  return static_cast<Id>(_meshes.size()) / getNApexPerMesh();
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
double MeshSpherical::getMeshSize(Id imesh) const
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
  sstr << toTitle(0, "Spherical Meshing characteristics");
  sstr << AMesh::toString(strfmt);
  return sstr.str();
}

/**
 * Create a MeshSpherical by loading the contents of a Neutral File
 *
 * @param NFFilename Name of the Neutral File (MeshEStandard format)
 * @param verbose    Verbose
 */
MeshSpherical* MeshSpherical::createFromNF(const String& NFFilename, bool verbose)
{
  MeshSpherical* mesh = new MeshSpherical;
  if (mesh->_fileOpenAndDeserialize(NFFilename, verbose)) return mesh;
  delete mesh;
  return nullptr;
}

MeshSpherical* MeshSpherical::create(const MatrixDense& apices,
                                     const MatrixInt& meshes)
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
Id MeshSpherical::reset(Id ndim,
                        Id napexpermesh,
                        const VectorDouble& apices,
                        const VectorInt& meshes,
                        bool byCol,
                        bool verbose)
{
  _setNDim(ndim);
  Id npoints = static_cast<Id>(apices.size()) / ndim;
  Id nmeshes = static_cast<Id>(meshes.size()) / napexpermesh;

  // Core allocation

  _apices.reset(npoints, ndim);
  _apices.setValues(apices, byCol);
  _meshes.reset(nmeshes, napexpermesh);
  _meshes.setValues(meshes, byCol);

  // Check consistency

  _checkConsistency();

  // Define and store the Bounding Box extension

  _defineBoundingBox();

  // Optional printout

  if (verbose) messageFlush(toString());

  return (0);
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
Id MeshSpherical::getApex(Id imesh, Id rank) const
{
  return _meshes.getValue(imesh, rank);
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
double MeshSpherical::getCoor(Id imesh, Id rank, Id idim) const
{
  return _apices(getApex(imesh, rank), idim);
}

double MeshSpherical::getApexCoor(Id i, Id idim) const
{
  return _apices(i, idim);
}

void MeshSpherical::_getCoordOnSphere(double longitude,
                                      double latitude,
                                      VectorDouble& coords)
{
  double radius = EARTH_RADIUS;
  if (isDefaultSpaceSphere())
  {
    const ASpace* space = getDefaultSpaceSh().get();
    const auto* spaceSn = dynamic_cast<const SpaceSN*>(space);
    if (spaceSn != nullptr) radius = spaceSn->getRadius();
  }
  GH::convertSph2Cart(longitude, latitude,
                      &coords.at(0), &coords.at(1), &coords.at(2), radius);
}

void MeshSpherical::getEmbeddedCoorPerMesh(Id imesh, Id ic, VectorDouble& coords) const
{
  _getCoordOnSphere(getCoor(imesh, ic, 0), getCoor(imesh, ic, 1), coords);
}

void MeshSpherical::getEmbeddedCoorPerApex(Id iapex, VectorDouble& coords) const
{
  _getCoordOnSphere(getApexCoor(iapex, 0), getApexCoor(iapex, 1), coords);
}

/**
 * Calculate the Mesh of each Mesh (using approximated calculations)
 * @return
 */
VectorDouble MeshSpherical::_defineUnits(void) const
{
  auto nmeshes = getNMeshes();
  VectorDouble units(nmeshes);
  for (Id imesh = 0; imesh < nmeshes; imesh++)
  {
    VectorVectorDouble corners = getCoordinatesPerMesh(imesh);
    units[imesh]               = _getMeshUnit(corners);
  }
  return units;
}

void MeshSpherical::_defineBoundingBox(void)
{
  VectorDouble extendmin;
  VectorDouble extendmax;
  double coor, mini, maxi;
  auto ndim = getNDim();

  // Initializations
  extendmin.resize(ndim);
  extendmax.resize(ndim);

  // Loop on the Space dimensions
  for (Id idim = 0; idim < ndim; idim++)
  {
    mini = MAXIMUM_BIG;
    maxi = MINIMUM_BIG;

    // Loop on the apices
    for (Id i = 0; i < getNApices(); i++)
    {
      coor = getApexCoor(i, idim);
      if (coor < mini) mini = coor;
      if (coor > maxi) maxi = coor;
    }
    extendmin[idim] = mini;
    extendmax[idim] = maxi;
  }

  // Store the Bounding Box extension
  (void)_setExtend(extendmin, extendmax);
}

double MeshSpherical::_closestValue(double ref, double coor, double period)
{
  double dref = ABS(coor - ref);
  double d1   = ABS(coor - period - ref);
  if (d1 < dref) return coor - period;
  return coor;
}

Id MeshSpherical::_recopy(const MeshSpherical& m)
{
  _apices = m._apices;
  _meshes = m._meshes;
  AMesh::_recopy(m);
  return (0);
}

bool MeshSpherical::_deserializeAscii(std::istream& is, bool /*verbose*/)
{
  Id ndim         = 0;
  Id napices      = 0;
  Id nmeshes      = 0;
  Id napexpermesh = 0;

  bool ret = true;
  ret      = ret && _recordRead<Id>(is, "Space Dimension", ndim);
  ret      = ret && _recordRead<Id>(is, "Napices", napices);
  ret      = ret && _recordRead<Id>(is, "Number of Apices per Mesh", napexpermesh);
  ret      = ret && _recordRead<Id>(is, "Number of Meshes", nmeshes);

  if (ret)
  {
    VectorDouble apices_local;
    ret     = ret && _recordReadVec<double>(is, "Apices", apices_local, ndim * napices);
    _apices = MatrixDense(napices, ndim);
    _apices.setValues(apices_local);
  }

  if (ret)
  {
    VectorInt meshes_local;
    ret     = ret && _recordReadVec<Id>(is, "Meshes", meshes_local, nmeshes * napexpermesh);
    _meshes = MatrixInt(nmeshes, napexpermesh);
    _meshes.setValues(meshes_local);
  }
  return ret;
}

bool MeshSpherical::_serializeAscii(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;
  ret      = ret && _recordWrite<Id>(os, "Space Dimension", getNDim());
  ret      = ret && _recordWrite<Id>(os, "Napices", getNApices());
  ret      = ret && _recordWrite<Id>(os, "Number of Apices per Mesh", getNApexPerMesh());
  ret      = ret && _recordWrite<Id>(os, "Number of Meshes", getNMeshes());
  ret      = ret && _recordWriteVec<double>(os, "Apices", _apices.getValues());
  ret      = ret && _recordWriteVec<Id>(os, "Meshes", _meshes.getValues());
  return ret;
}

/**
 * This function checks the consistency between the number of points
 * and the vertices indication
 */
void MeshSpherical::_checkConsistency() const
{
  for (Id imesh = 0; imesh < getNMeshes(); imesh++)
    for (Id ic = 0; ic < getNApexPerMesh(); ic++)
    {
      auto apex = getApex(imesh, ic);
      if (apex < 0 || apex >= getNApices())
      {
        my_throw("Mesh indices are not compatible with the Points");
      }
    }
}

void MeshSpherical::getBarycenterInPlace(Id imesh, vect coord) const
{
  auto ndimE   = getEmbeddedNDim();
  auto ncorner = getNApexPerMesh();

  // Calculate the center of gravity (in the Embedded space)
  VectorVectorDouble coordE = getEmbeddedCoordinatesPerMesh(imesh);
  VectorDouble centerE(ndimE);
  double rlong;
  double rlat;

  for (Id idimE = 0; idimE < ndimE; idimE++)
  {
    double local = 0.;
    for (Id ic = 0; ic < ncorner; ic++)
      local += coordE[ic][idimE];
    centerE[idimE] = local / ncorner;
  }

  // Turn the gravity center from embedded to long/lat coordinates
  GH::convertCart2Sph(centerE[0], centerE[1], centerE[2],
                      &rlong, &rlat, TEST);
  coord[0] = rlong;
  coord[1] = rlat;
}

#ifdef HDF5
bool MeshSpherical::_deserializeH5(H5::Group& grp, [[maybe_unused]] bool verbose)
{
  auto meshG = SerializeHDF5::getGroup(grp, "MeshSpherical");
  if (!meshG)
  {
    return false;
  }

  /* Read the grid characteristics */
  bool ret    = true;
  Id ndim     = 0;
  Id napices  = 0;
  Id npermesh = 0;
  Id nmeshes  = 0;
  VectorDouble apices;
  VectorInt meshes;

  ret = ret && SerializeHDF5::readValue(*meshG, "NDim", ndim);
  ret = ret && SerializeHDF5::readValue(*meshG, "NApices", napices);
  ret = ret && SerializeHDF5::readValue(*meshG, "NPerMesh", npermesh);
  ret = ret && SerializeHDF5::readValue(*meshG, "NMeshes", nmeshes);
  ret = ret && SerializeHDF5::readVec(*meshG, "Apices", apices);
  ret = ret && SerializeHDF5::readVec(*meshG, "Meshes", meshes);

  if (ret)
  {
    _apices = MatrixDense(napices, ndim);
    _apices.setValues(apices);
    _meshes = MatrixInt(nmeshes, npermesh);
    _meshes.setValues(meshes);
  }
  return ret;
}

bool MeshSpherical::_serializeH5(H5::Group& grp, [[maybe_unused]] bool verbose) const
{
  auto meshG = grp.createGroup("MeshSpherical");

  bool ret = true;

  ret = ret && SerializeHDF5::writeValue(meshG, "NDim", getNDim());
  ret = ret && SerializeHDF5::writeValue(meshG, "NApices", getNApices());
  ret = ret && SerializeHDF5::writeValue(meshG, "NPerMesh", getNApexPerMesh());
  ret = ret && SerializeHDF5::writeValue(meshG, "NMeshes", getNMeshes());
  ret = ret && SerializeHDF5::writeVec(meshG, "Apices", _apices.getValues());
  ret = ret && SerializeHDF5::writeVec(meshG, "Meshes", _meshes.getValues());

  return ret;
}
#endif
} // namespace gstlrn
