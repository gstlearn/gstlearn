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
#include "Matrix/MatrixDense.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "Mesh/MeshEStandard.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Basic/AException.hpp"
#include "Basic/SerializeHDF5.hpp"

namespace gstlrn
{ 

MeshEStandard::MeshEStandard()
  : AMesh()
  , _apices()
  , _meshes()
{
}

MeshEStandard::MeshEStandard(const MeshEStandard &m) 
  : AMesh(m)
{
  _recopy(m);
}

MeshEStandard& MeshEStandard::operator= (const MeshEStandard &m)
{
  _recopy(m);
  return *this;
}

MeshEStandard::~MeshEStandard()
{
  _deallocate();
}

/****************************************************************************/
/*!
** Returns the number of Apices
**
** \returns Number of apices
**
*****************************************************************************/
Id MeshEStandard::getNApices() const 
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
Id MeshEStandard::getNMeshes() const
{
  return (static_cast<Id> (_meshes.size()) / getNApexPerMesh());
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
double MeshEStandard::getMeshSize(Id imesh) const
{
  VectorVectorDouble corners = getCoordinatesPerMesh(imesh);
  return _getMeshUnit(corners);
}

Id MeshEStandard::resetFromTurbo(const MeshETurbo& turbo, bool verbose)
{
  auto ndim     = turbo.getNDim();
  auto napices  = turbo.getNApices();
  auto nmeshes  = turbo.getNMeshes();
  auto npermesh = turbo.getNApexPerMesh();

  // Dimension the members

  _apices = MatrixDense(napices, ndim);
  _meshes = MatrixInt(nmeshes, npermesh);

  // Load the apices;
  VectorDouble local(ndim);
  for (Id ip = 0; ip < napices; ip++)
  {
    turbo.getApexCoordinatesInPlace(ip, local);
    for (Id idim = 0; idim < ndim; idim++)
      _apices.setValue(ip, idim, local[idim]);
  }

  // Load the meshes
  for (Id imesh = 0; imesh < nmeshes; imesh++)
    for (Id rank = 0; rank < npermesh; rank++)
      _meshes.setValue(imesh,  rank, turbo.getApex(imesh, rank));

  // Check consistency

  _checkConsistency();

  // Define and store the Bounding Box extension

  _defineBoundingBox();

  // Optional printout

  if (verbose) messageFlush(toString());

  return 0;
}

/****************************************************************************/
/*!
** Create the meshing (from mesh information)
**
** \param[in]  apices          Matrix of Apices
** \param[in]  meshes          Array of mesh indices
** \param[in]  verbose         Verbose flag
**
*****************************************************************************/
Id MeshEStandard::reset(const MatrixDense& apices,
                         const MatrixInt& meshes,
                         bool verbose)
{
  auto ndim = apices.getNCols();
  _setNDim(ndim);

  // Core allocation

  _apices = apices;
  _meshes = meshes;

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
Id MeshEStandard::reset(Id ndim,
                         Id napexpermesh,
                         const VectorDouble &apices,
                         const VectorInt &meshes,
                         bool byCol,
                         bool verbose)
{
  _setNDim(ndim);
  Id npoints = static_cast<Id> (apices.size()) / ndim;
  Id nmeshes = static_cast<Id> (meshes.size()) / napexpermesh;

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
** Returns the rank of the Apex 'rank' of the Mesh 'imesh'
**
** \returns The rank of the target apex
**
** \param[in]  imesh    Rank of Mesh (from 0 to _nMeshes-1))
** \param[in]  rank     Rank of Apex within a Mesh (from 0 to _nApexPerMesh)
**
*****************************************************************************/
Id MeshEStandard::getApex(Id imesh, Id rank) const
{
  return _meshes.getValue(imesh,rank);
}

void MeshEStandard::_setApex(Id imesh, Id rank, Id value)
{
  _meshes.setValue(imesh,rank,value);
}

/****************************************************************************/
/*!
** Returns the coordinate 'idim' of the Apex 'rank' of the Mesh 'imesh'
**
** \returns The coordinate of the target apex
**
** \param[in]  imesh    Rank of the Mesh (from 0 to _nMeshes-1))
** \param[in]  rank     Rank of Apex within a Mesh (from 0 to _nApexPerMesh-1)
** \param[in]  idim     Rank of the coordinate (from 0 to _nDim-1)
**
*****************************************************************************/
double MeshEStandard::getCoor(Id imesh,
                              Id rank,
                              Id idim) const
{
  return _apices(getApex(imesh,rank),idim);
}

/****************************************************************************/
/*!
** Print the contents of the meshing
**
** \param[in] strfmt    Format for printout
**
*****************************************************************************/
String MeshEStandard::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  sstr << toTitle(0,"Standard Meshing");
  sstr << AMesh::toString(strfmt);
  return sstr.str();
}

/**
 * Create a MeshEStandard by loading the contents of a Neutral File
 *
 * @param NFFilename Name of the Neutral File (MeshEStandard format)
 * @param verbose    Verbose
 */
MeshEStandard* MeshEStandard::createFromNF(const String& NFFilename, bool verbose)
{
  MeshEStandard* mesh = new MeshEStandard;
  if (mesh->_fileOpenAndDeserialize(NFFilename, verbose)) return mesh;
  delete mesh;
  return nullptr;
}

MeshEStandard* MeshEStandard::createFromExternal(const MatrixDense &apices,
                                                 const MatrixInt &meshes,
                                                 bool verbose)
{
  MeshEStandard* mesh = new MeshEStandard;
  mesh->reset(apices, meshes, verbose);
  return mesh;
}

/**
 * Returns the list of mesh vertex information
 * This List is organized as a single Vector of Double
 * It is dimensioned to Nrows=getNApices() and Ncols=getNDim()
 * @param byCol true if the values must be sorted by column
 * @return
 */
VectorDouble MeshEStandard::getPointList(bool byCol) const
{
  VectorDouble list;

  if (byCol)
  {
    for (Id idim = 0; idim < getNDim(); idim++)
      for (Id irow = 0; irow < getNApices(); irow++)
      {
        list.push_back(getApexCoor(irow, idim));
      }
  }
  else
  {
    for (Id irow = 0; irow < getNApices(); irow++)
      for (Id idim = 0; idim < getNDim(); idim++)
      {
        list.push_back(getApexCoor(irow, idim));
      }
  }
  return list;
}

double MeshEStandard::getApexCoor(Id i, Id idim) const
{
  return _apices(i, idim);
}

/****************************************************************************/
/*!
**  Check if a point, defined by its coordinates, belong to a Mesh
**
** \return true if the point belongs to the Mesh; false otherwise
**
** \param[in]  coor      Array of target coordinates
** \param[in]  imesh     Mesh Index
** \param[in]  meshsize  Dimension of the mesh
**
** \param[out] weights   Array of barycentric weights (Dim: NApexPerMesh)
** 
** \remarks The argument 'meshsize' is used to speed the calculations
**
*****************************************************************************/
bool MeshEStandard::_coorInMesh(const VectorDouble& coor,
                                Id imesh,
                                double meshsize,
                                VectorDouble& weights) const
{
  VectorVectorDouble corners = getCoordinatesPerMesh(imesh);
  return _weightsInMesh(coor, corners, meshsize, weights);
}

void MeshEStandard::_deallocate()
{

}

Id MeshEStandard::_recopy(const MeshEStandard &m)
{
  _deallocate();

  _apices = m._apices;
  _meshes = m._meshes;
  AMesh::_recopy(m);
  return(0);
}

/**
 * This function checks the consistency between the number of points
 * and the vertices indication
 */
void MeshEStandard::_checkConsistency() const
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

bool MeshEStandard::_deserializeAscii(std::istream& is, bool verbose)
{
  DECLARE_UNUSED(verbose);
  Id ndim = 0;
  Id napices = 0;
  Id napexpermesh = 0;
  Id nmeshes = 0;

  bool ret = true;
  ret = ret && _recordRead<Id>(is, "Space Dimension", ndim);
  ret = ret && _recordRead<Id>(is, "Napices", napices);
  ret = ret && _recordRead<Id>(is, "Number of Apices per Mesh", napexpermesh);
  ret = ret && _recordRead<Id>(is, "Number of Meshes", nmeshes);

  if (ret)
  {
    VectorDouble apices_local;
    ret = ret && _recordReadVec<double>(is, "Apices", apices_local, napices * ndim);
    _apices = MatrixDense(napices, ndim);
    _apices.setValues(apices_local);
  }

  if (ret)
  {
    VectorInt meshes_local;
    ret = ret && _recordReadVec<Id>(is, "Meshes", meshes_local, nmeshes * napexpermesh);
    _meshes = MatrixInt(nmeshes, napexpermesh);
    _meshes.setValues(meshes_local);
  }
  return ret;
}

bool MeshEStandard::_serializeAscii(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;

  ret = ret && _recordWrite<Id>(os, "Space Dimension", getNDim());
  ret = ret && _recordWrite<Id>(os, "Napices", getNApices());
  ret = ret && _recordWrite<Id>(os, "Number of Apices per Mesh", getNApexPerMesh());
  ret = ret && _recordWrite<Id>(os, "Number of Meshes", getNMeshes());
  ret = ret && _recordWriteVec<double>(os, "Apices", _apices.getValues());
  ret = ret && _recordWriteVec<Id>(os, "Meshes", _meshes.getValues());
  return ret;
}

/**
 * Validate Meshing, in particular when transiting from Old Meshing to New one
 *
 * @remark For safety and considering the rank of Apices stored in 'meshes',
 * @remark the minimum value is considered
 * @remark If not equal to 0 or 1, a fatal error is issued
 * @remark If equal to 1 (old style), values are decreased by 1.
 */
void MeshEStandard::_validate()
{
  // For safety, the minimum value of the array meshes is considered
  // If equal to , the mesh numbering must be modified in order to start from 0
  // This code is temporary and there to cope with different numberings

  auto nmeshes  = getNMeshes();
  auto ncorner  = getNApexPerMesh();
  Id shift_min = 1000;
  for (Id imesh = 0; imesh < nmeshes; imesh++)
    for (Id ic = 0; ic < ncorner; ic++)
    {
      auto ipos = getApex(imesh, ic);
      if (ipos < shift_min) shift_min = ipos;
    }

  if (shift_min != 0 && shift_min != 1)
    my_throw("Wrong minimum shift: it should be 0 or 1");

  if (shift_min == 1)
  {
    for (Id imesh = 0; imesh < nmeshes; imesh++)
      for (Id ic = 0; ic < ncorner; ic++)
        _setApex(imesh, ic, getApex(imesh, ic) - shift_min);
  }
}

void MeshEStandard::_defineBoundingBox(void)
{
  VectorDouble extendmin, extendmax;
  double coor,mini,maxi;
  auto ndim = getNDim();

  // Initializations
  extendmin.resize(ndim);
  extendmax.resize(ndim);

  // Loop on the Space dimensions
  for (Id idim=0; idim<ndim; idim++)
  {
    mini = MAXIMUM_BIG;
    maxi = MINIMUM_BIG;

    // Loop on the apices
    for (Id i=0; i<getNApices(); i++)
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

#ifdef HDF5
bool MeshEStandard::_deserializeH5(H5::Group& grp, [[maybe_unused]] bool verbose)
{
  auto meshG = SerializeHDF5::getGroup(grp, "MeshEStandard");
  if (!meshG)
  {
    return false;
  }

  /* Read the grid characteristics */
  bool ret = true;
  Id ndim = 0;
  Id napices = 0;
  Id npermesh = 0;
  Id nmeshes = 0;
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

bool MeshEStandard::_serializeH5(H5::Group& grp, [[maybe_unused]] bool verbose) const
{
  auto meshG = grp.createGroup("MeshEStandard");

  bool ret = true;
  ret      = ret && SerializeHDF5::writeValue(meshG, "NDim", getNDim());
  ret      = ret && SerializeHDF5::writeValue(meshG, "NApices", getNApices());
  ret      = ret && SerializeHDF5::writeValue(meshG, "NPerMesh", getNApexPerMesh());
  ret      = ret && SerializeHDF5::writeValue(meshG, "NMeshes", getNMeshes());
  ret      = ret && SerializeHDF5::writeVec(meshG, "Apices", _apices.getValues());
  ret      = ret && SerializeHDF5::writeVec(meshG, "Meshes", _meshes.getValues());

  return ret;
}
#endif
}