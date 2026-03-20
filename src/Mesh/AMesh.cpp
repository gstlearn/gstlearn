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
#include "Mesh/AMesh.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixInt.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixSquare.hpp"
#include "Matrix/NF_Triplet.hpp"
#include "Space/SpacePoint.hpp"
#include "Tree/Ball.hpp"

#include <algorithm>

namespace gstlrn
{
AMesh::AMesh()
  : AStringable()
  , ASerializable()
  , _nDim(0)
  , _extendMin()
  , _extendMax()
  , _adjacencyMatrix(new MatrixSparse())
  , _ballTree(nullptr)
{
}

AMesh::AMesh(const AMesh& m)
  : AStringable(m)
  , ASerializable(m)
{
  _recopy(m);
}

AMesh& AMesh::operator=(const AMesh& m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    ASerializable::operator=(m);
    _recopy(m);
  }
  return *this;
}

AMesh::~AMesh()
{
  delete _adjacencyMatrix;
}

Id AMesh::_setExtend(const VectorDouble& extendmin,
                     const VectorDouble& extendmax)
{
  _extendMin = extendmin;
  _extendMax = extendmax;
  return (0);
}

String AMesh::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  if (_nDim <= 0) return sstr.str();

  if (getVariety() == 0)
    sstr << "Euclidean Geometry" << std::endl;
  else
    sstr << "Geometry defined on the Sphere" << std::endl;

  sstr << "Space Dimension           = " << _nDim << std::endl;
  sstr << "Number of Apices per Mesh = " << getNApexPerMesh() << std::endl;
  sstr << "Number of Meshes          = " << getNMeshes() << std::endl;
  sstr << "Number of Apices          = " << getNApices() << std::endl;

  if (!_extendMin.empty() && !_extendMax.empty())
  {
    sstr << toStrTitle(1, "Bounding Box Extension");
    for (Id idim = 0; idim < _nDim; idim++)
      sstr << "Dim #" << idim + 1 << " - Min:" << _extendMin[idim] << " - Max:" << _extendMax[idim] << std::endl;
  }

  AStringFormat sf;
  if (strfmt != nullptr) sf = *strfmt;
  if (sf.getLevel() > 1)
  {
    MatrixDense apices;
    MatrixInt meshes;
    getElements(apices, meshes);

    sstr << "List of Apices" << std::endl;
    sstr << apices.toString(strfmt);
    sstr << "List of Meshes" << std::endl;
    sstr << meshes.toString(strfmt);
  }
  return sstr.str();
}

void AMesh::getCoordinatesPerMeshInPlace(Id imesh, Id rank, VectorDouble& coords) const
{
  for (Id idim = 0; idim < getNDim(); idim++)
    coords[idim] = getCoor(imesh, rank, idim);
}

void AMesh::getApexCoordinatesInPlace(Id i, VectorDouble& coords) const
{
  for (Id idim = 0; idim < getNDim(); idim++)
    coords[idim] = getApexCoor(i, idim);
}

/****************************************************************************/
/*!
** Checks that the Db is compatible with the Meshing
**
** \return 1 if Db and Meshing are incompatible; 0 otherwise
**
** \param[in]  db        Db structure
**
*****************************************************************************/
Id AMesh::isCompatibleDb(const Db* db) const
{
  if (getNDim() == db->getNDim()) return 0;

  messerr("Inconsistent Space dimension between Meshing (%d) and Db (%d)",
          getNDim(), db->getNDim());
  return 1;
}

VectorDouble AMesh::getMeshSizes() const
{
  VectorDouble units;
  for (Id imesh = 0; imesh < getNMeshes(); imesh++)
    units.push_back(getMeshSize(imesh));
  return units;
}

void AMesh::printMesh(Id imesh0) const
{
  mestitle(0, "Mesh Information");
  message("- Number of Meshes = %d\n", getNMeshes());
  message("- Number of Apices = %d\n", getNApices());
  if (imesh0 < 0) return;

  Id ideb = (imesh0 >= 0) ? imesh0 : 0;
  Id ifin = (imesh0 >= 0) ? imesh0 + 1 : getNMeshes();
  for (Id imesh = ideb; imesh < ifin; imesh++)
  {
    message("Mesh #%d\n", imesh + 1);
    for (Id icorn = 0; icorn < getNApexPerMesh(); icorn++)
    {
      message("Point #%d", getApex(imesh, icorn));
      for (Id idim = 0; idim < getNDim(); idim++)
      {
        message(" %lf", getCoor(imesh, icorn, idim));
      }
      message("\n");
    }
  }
}

void AMesh::printMeshes(Id level, Id nline_max) const
{
  mestitle(0, "Mesh Information");
  message("- Number of Meshes = %d\n", getNMeshes());
  message("- Number of Apices = %d\n", getNApices());

  if (level == 0) return;

  if (level == 1)
    _printMeshListByIndices(nline_max);

  if (level == 2)
    _printMeshListByCoordinates(nline_max);
}

void AMesh::_recopy(const AMesh& m)
{
  _nDim      = m._nDim;
  _extendMin = m._extendMin;
  _extendMax = m._extendMax;
  delete _adjacencyMatrix;
  _adjacencyMatrix = new MatrixSparse(*m._adjacencyMatrix);
}

VectorDouble AMesh::getCoordinatesPerApex(Id idim) const
{
  if (!_isSpaceDimensionValid(idim)) return VectorDouble();
  auto np = getNApices();
  VectorDouble coor(np);
  for (Id ip = 0; ip < np; ip++)
    coor[ip] = getApexCoor(ip, idim);
  return coor;
}

/**
 * Returns the coordinates of all meshes:
 * - the first dimension if the space dimension
 * - the second dimension is the number of apices
 * @return
 */
VectorVectorDouble AMesh::getAllCoordinates() const
{
  auto napices = getNApices();
  VectorDouble local(_nDim);
  VectorVectorDouble coords(_nDim);
  for (Id idim = 0; idim < _nDim; idim++)
    coords[idim].resize(napices);

  for (Id ip = 0; ip < napices; ip++)
  {
    getApexCoordinatesInPlace(ip, local);
    for (Id idim = 0; idim < _nDim; idim++)
      coords[idim][ip] = local[idim];
  }
  return coords;
}

/**
 * Returns the information of all meshes:
 * - the first dimension is the number of apices (nrow)
 * - the second dimension if the space dimension (ncol)
 * @return
 */
MatrixInt AMesh::getAllMeshes() const
{
  auto nper    = getNApexPerMesh();
  auto nmeshes = getNMeshes();
  MatrixInt meshes(nmeshes, nper);

  for (Id imesh = 0; imesh < nmeshes; imesh++)
    for (Id iper = 0; iper < nper; iper++)
      meshes.setValue(imesh, iper, getApex(imesh, iper));
  return meshes;
}

/**
 * Returns the coordinates of the Center of Gravity of a Mesh
 * @param imesh Rank of the Mesh
 * @param idim Index of the space dimension
 * @return
 */
double AMesh::getCenterCoordinate(Id imesh, Id idim) const
{
  double coor  = 0.;
  auto ncorner = getNApexPerMesh();
  for (Id icorner = 0; icorner < ncorner; icorner++)
    coor += getCoor(imesh, icorner, idim);
  return (coor / static_cast<double>(ncorner));
}

VectorVectorDouble AMesh::getAllCenterCoordinates() const
{
  auto ncorner = getNApexPerMesh();
  auto nmeshes = getNMeshes();
  VectorVectorDouble coords(_nDim);
  for (Id idim = 0; idim < _nDim; idim++)
    coords[idim].resize(nmeshes);

  for (Id imesh = 0; imesh < nmeshes; imesh++)
  {
    for (Id idim = 0; idim < _nDim; idim++)
    {
      double total = 0.;
      for (Id ic = 0; ic < ncorner; ic++)
        total += getCoor(imesh, ic, idim);
      coords[idim][imesh] = total / ncorner;
    }
  }
  return coords;
}

/**
 * Returns the information about all apices:
 * - the first dimension is the number of meshes (nrow)
 * - the second dimension if the space dimension (ncol)
 * @return
 */
MatrixDense AMesh::getAllApices() const
{
  auto napices = getNApices();
  MatrixDense apices(napices, _nDim);
  for (Id ip = 0; ip < napices; ip++)
    for (Id idim = 0; idim < _nDim; idim++)
      apices.setValue(ip, idim, getApexCoor(ip, idim));
  return apices;
}

VectorInt AMesh::getMeshByApexPair(Id apex1, Id apex2) const
{
  VectorInt list;
  auto ncorner = getNApexPerMesh();
  Id found, apex0;

  for (Id imesh = 0; imesh < getNMeshes(); imesh++)
  {
    found = 0;
    for (Id ic = 0; ic < ncorner; ic++)
    {
      apex0 = getApex(imesh, ic);
      if (apex0 == apex1)
        found++;
      if (apex0 == apex2)
        found++;
      if (found == 2)
      {
        list.push_back(imesh);
        break;
      }
    }
  }
  return list;
}

/****************************************************************************/
/*!
** Extract the elements of the meshing
**
** \param[out]  apices  Pointer on the array of Apices
** \param[out]  meshes  Pointer on the array of Meshes
**
*****************************************************************************/
void AMesh::getElements(MatrixDense& apices, MatrixInt& meshes) const
{
  auto nmeshes = getNMeshes();
  auto ndim    = getNDim();
  auto napices = getNApices();
  auto ncorner = getNApexPerMesh();

  // Dimension the returned containers

  apices.reset(napices, ndim);
  meshes.reset(nmeshes, ncorner);

  // Load the Apices

  VectorDouble local(ndim);
  for (Id i = 0; i < napices; i++)
  {
    getApexCoordinatesInPlace(i, local);
    for (Id idim = 0; idim < ndim; idim++)
      apices.setValue(i, idim, local[idim]);
  }

  // Load the Meshes

  for (Id imesh = 0; imesh < nmeshes; imesh++)
    for (Id icorner = 0; icorner < ncorner; icorner++)
      meshes.setValue(imesh, icorner, getApex(imesh, icorner));
}

bool AMesh::_isSpaceDimensionValid(Id idim) const
{
  return checkArg("SPace Dimension Index", idim, _nDim);
}

VectorDouble AMesh::getExtrema(Id idim) const
{
  VectorDouble vec(2);
  vec[0] = getExtendMin(idim);
  vec[1] = getExtendMax(idim);
  return vec;
}

VectorDouble AMesh::getCoordinatesPerMesh(Id imesh, Id idim, bool flagClose) const
{
  VectorDouble vec;
  auto ncorner = getNApexPerMesh();

  if (flagClose)
    vec.resize(ncorner + 1);
  else
    vec.resize(ncorner);

  for (Id ic = 0; ic < ncorner; ic++)
    vec[ic] = getCoor(imesh, ic, idim);
  if (flagClose) vec[ncorner] = getCoor(imesh, 0, idim);

  return vec;
}

/**
 * Fill the coordinates of a corner of a mesh in embedded space
 * @param imesh  Mesh rank
 * @param ic     Corner index
 * @param coords Array of coordinates
 */
void AMesh::getEmbeddedCoorPerMesh(Id imesh, Id ic, VectorDouble& coords) const
{
  getCoordinatesPerMeshInPlace(imesh, ic, coords);
}

/**
 * Fill the coordinates of an apex in embedded space
 * @param iapex  Apex index
 * @param coords Array of coordinates
 */
void AMesh::getEmbeddedCoorPerApex(Id iapex, VectorDouble& coords) const
{
  getApexCoordinatesInPlace(iapex, coords);
}

/*! Returns the Sparse Matrix for projecting the Mesh to a Db */
ProjMatrix* AMesh::createProjMatrix(const Db* db, Id rankZ, bool verbose) const
{
  auto* m = new ProjMatrix();
  resetProjFromDb(m, db, rankZ, verbose);
  return m;
}

/**
 * Fill the array of coordinates of all apices of a mesh in embedded space
 * Storage [ndim, ncorner]
 * @param imesh Mesh rank
 * @param vec   Returned array
 */
void AMesh::getEmbeddedCoordinatesPerMeshInPlace(Id imesh, VectorVectorDouble& vec) const
{
  auto ncorner = getNApexPerMesh();

  for (Id ic = 0; ic < ncorner; ic++)
    getEmbeddedCoorPerMesh(imesh, ic, vec[ic]);
}

/**
 * Returns the array of coordinates of all apices of any mesh in embedded space
 * Its dimensions are: ncorner * ndim
 * @param imesh Mesh rank
 * @return
 */
VectorVectorDouble AMesh::getEmbeddedCoordinatesPerMesh(Id imesh) const
{
  auto ndim    = getEmbeddedNDim();
  auto ncorner = getNApexPerMesh();
  VectorVectorDouble vec(ncorner);
  for (auto& e: vec)
    e = VectorDouble(ndim);

  getEmbeddedCoordinatesPerMeshInPlace(imesh, vec);
  return vec;
}

VectorVectorDouble AMesh::getCoordinatesPerMesh(Id imesh) const
{
  auto ndim    = getNDim();
  auto ncorner = getNApexPerMesh();
  VectorVectorDouble vec(ncorner);
  for (auto& e: vec)
    e = VectorDouble(ndim);

  for (Id ic = 0; ic < ncorner; ic++)
    for (Id idim = 0; idim < ndim; idim++)
      vec[ic][idim] = getCoor(imesh, ic, idim);
  return vec;
}

/**
 * Returns the coordinates of the Mesh apices expressed in the embedded space
 * The returned vector is organized by coordinate
 * @return
 */
VectorVectorDouble AMesh::getEmbeddedCoordinatesPerApex() const
{
  auto ndim    = getEmbeddedNDim();
  auto napices = getNApices();
  VectorVectorDouble vec(ndim);
  for (auto& e: vec)
    e = VectorDouble(napices);

  VectorDouble local(ndim);
  for (Id ip = 0; ip < napices; ip++)
  {
    getEmbeddedCoorPerApex(ip, local);
    for (Id idim = 0; idim < ndim; idim++)
      vec[idim][ip] = local[idim];
  }
  return vec;
}

VectorDouble AMesh::getApexCoordinates(Id iapex) const
{
  VectorDouble vec(_nDim);
  for (Id idim = 0; idim < _nDim; idim++)
    vec[idim] = getApexCoor(iapex, idim);
  return vec;
}

VectorInt AMesh::getMeshApicesFromCoordinates(const VectorDouble& coords) const
{
  auto ncorner = getNApexPerMesh();
  VectorInt meshApices(ncorner, -1);

  Id found  = getMeshFromCoordinates(coords);
  Id ip_max = 0;

  if (found >= 0)
  {
    for (Id icorn = 0; icorn < ncorner; icorn++)
    {
      auto ip = getApex(found, icorn);
      if (ip > ip_max) ip_max = ip;
      meshApices[icorn] = ip;
    }
  }

  return meshApices;
}

Id AMesh::getMeshFromCoordinates(const VectorDouble& coords) const
{
  auto ncorner = getNApexPerMesh();
  VectorDouble weights(ncorner, 0);

  return getMeshAndInPlaceWeightsFromCoordinates(coords, weights);
}

Id AMesh::getMeshAndInPlaceWeightsFromCoordinates(const VectorDouble& coords, VectorDouble& weights) const
{
  auto nmeshes = getNMeshes();
  VectorDouble units(nmeshes, 0.);
  if (getVariety() != 1)
    units = _defineUnits();

  /* Instantiate a Ball Tree for quick search */
  // Note: this Ball tree is defined in 3D despite the space dimension of mesh
  if (_ballTree == nullptr)
    _ballTree = std::make_unique<Ball>(this, 10, true);
  // if (verbose) ball.display(1);

  /* Loop on the samples */
  VectorInt neighs;
  VectorDouble distances;

  /* Loop on the elligible meshes */
  Id nb_neigh = MIN(5, getNMeshes());
  ;
  (void)_ballTree->queryOneInPlace(coords, nb_neigh, neighs, distances);
  Id found = _findBarycenter(coords, units, nb_neigh, neighs, weights);

  // If search has failed with a small number of neighbors, try with a larger one
  if (found < 0)
  {
    nb_neigh = MIN(50, getNMeshes());

    (void)_ballTree->queryOneInPlace(coords, nb_neigh, neighs, distances);
    found = _findBarycenter(coords, units, nb_neigh, neighs, weights);
  }
  if (found < 0)
  {
    Id ndim = getNDim();

    /* Printout if a point does not belong to any mesh */
    if (ndim == 3)
    {
      messerr("Point (%lf %lf %lf) does not belong to any mesh (nb_neigh=%d)",
              coords[0], coords[1], coords[2], nb_neigh);
    }
    else
    {
      messerr("Point (%lf %lf) does not belong to any mesh (nb_neigh=%d)",
              coords[0], coords[1], nb_neigh);
    }
  }

  return found;
}

VectorDouble AMesh::getDistances(Id iapex0, const VectorInt& japices) const
{
  VectorInt jlocal = japices;
  if (jlocal.empty()) jlocal = VH::sequence(getNApices());
  Id number = static_cast<Id>(jlocal.size());
  VectorDouble vec(number, 0.);
  //TODO[space]: Use default space
  SpacePoint P1(getApexCoordinates(iapex0), -1);

  for (Id iapex = 0; iapex < number; iapex++)
  {
    //TODO[space]: Use default space
    SpacePoint P2(getApexCoordinates(jlocal[iapex]), -1);
    vec[iapex] = P1.getDistance(P2);
  }
  return vec;
}

/**
 * Returns the list of neighboring meshes
 * This is a complex structure which stands as a vector of vectors of integers
 * - the first dimension is the number of apices
 * - for each apex, the second vector gives the indices of the neighboring meshes
 * @return
 */
VectorVectorInt AMesh::getNeighborhoodPerMesh() const
{
  auto napices  = getNApices();
  auto nmeshes  = getNMeshes();
  auto npermesh = getNApexPerMesh();

  VectorVectorInt Vmesh;
  Vmesh.resize(napices);

  // Loop on the meshes

  for (Id imesh = 0; imesh < nmeshes; imesh++)
  {
    // Loop on the apices for each mesh

    for (Id rank = 0; rank < npermesh; rank++)
    {
      auto ip = getApex(imesh, rank);
      Vmesh[ip].push_back(imesh);
    }
  }
  return Vmesh;
}

/**
 * Returns the list of neighboring apices
 * This is a complex structure which stands as a vector of vectors of integers
 * - the first dimension is the number of apices
 * - for each apex, the second vector gives the indices of the neighboring apices
 * @return
 */
VectorVectorInt AMesh::getNeighborhoodPerApex() const
{
  auto napices  = getNApices();
  auto npermesh = getNApexPerMesh();

  VectorVectorInt Vapex;
  Vapex.resize(napices);

  // Elaborate the Meshing Neighborhood by Mesh first
  VectorVectorInt Vmesh = getNeighborhoodPerMesh();

  for (Id ip = 0; ip < napices; ip++)
  {
    VectorInt vec;

    // Loop on the meshes neighboring the target apex 'ip'
    Id nmesh = static_cast<Id>(Vmesh[ip].size());
    for (Id i = 0; i < nmesh; i++)
    {
      // Index of the neighboring mesh
      Id imesh = Vmesh[ip][i];

      // Loop on the apices of the neighboring mesh
      for (Id rank = 0; rank < npermesh; rank++)
      {
        auto jp = getApex(imesh, rank);

        // Skip the target apex itself
        if (jp == ip) continue;

        // Add the new apex to the list
        vec.push_back(jp);
      }
    }

    // Sort the list and suppress the duplicates
    std::sort(vec.begin(), vec.end());
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());

    // Store the list
    Vapex[ip] = vec;
  }
  return Vapex;
}

VectorInt AMesh::getAdjacentApices(Id iapex) const
{
  if (_adjacencyMatrix->empty())
  {
    buildAdjacencyMatrix();
  }
  return _adjacencyMatrix->getNonZeroCols(iapex);
}

thread_local VectorBool _visitedWorkspace;

VectorInt AMesh::getNRingsAdjacentApices( Id start_apex, Id n_rings) const
{
    if (n_rings <= 0) return {start_apex};
    
    VectorInt all_neighbors;
    VectorInt current_layer;
    size_t nApices = getNApices();

    if (_visitedWorkspace.size() != nApices) {
        _visitedWorkspace.fill(false, nApices);
    }

    _visitedWorkspace[start_apex] = true;
    current_layer.push_back(start_apex);
    all_neighbors.push_back(start_apex);

    for (Id r = 0; r < n_rings; ++r) {
        VectorInt next_layer;
        for (Id apex : current_layer) {
            VectorInt neighbors = getAdjacentApices(apex);
            
            for (Id v : neighbors) {
                if (!_visitedWorkspace[v]) {
                    _visitedWorkspace[v] = true;
                    next_layer.push_back(v);
                    all_neighbors.push_back(v);
                }
            }
        }
        if (next_layer.empty()) break;
        current_layer = next_layer; // TODO use std::move
    }

    for (Id idx : all_neighbors) {
        _visitedWorkspace[idx] = false;
    }

    return all_neighbors;
}

const MatrixSparse* AMesh::getAdjacentApicesMatrix() const
{
  if (_adjacencyMatrix->empty())
  {
    buildAdjacencyMatrix();
  }
  return _adjacencyMatrix;
}

void AMesh::buildAdjacencyMatrix() const
{
  delete _adjacencyMatrix;
  Id nmeshes         = getNMeshes();
  Id ncorner         = getNApexPerMesh();

  // Define Sl as the sparse matrix giving the clutter of apices among vertices
  NF_Triplet NF_T;
  for (Id imesh = 0; imesh < nmeshes; imesh++)
  {
    for (Id ic = 0; ic < ncorner; ic++)
    {
      Id iapex = getApex(imesh, ic);
      NF_T.add(iapex, imesh, 1.);
    }
  }
  auto* Sl = MatrixSparse::createFromTriplet(NF_T);

  // Operate the product Sl * t(Sl) to get the final matrix Sret
  _adjacencyMatrix = prodNormMat(Sl, false);
  delete Sl;

  // Blank out the contents of the sparse matrix
  _adjacencyMatrix->setConstant(1.);

}

void AMesh::dumpNeighborhood(std::vector<VectorInt>& Vmesh, Id nline_max)
{
  mestitle(1, "List of Meshing Neighborhood");
  Id nmax = static_cast<Id>(Vmesh.size());
  if (nline_max > 0) nmax = MIN(nmax, nline_max);
  for (Id irow = 0; irow < nmax; irow++)
  {
    printVector(Vmesh[irow], String(), true, true);
  }
}

bool AMesh::_deserializeAscii(std::istream& is)
{
  bool ret = true;
  ret      = ret && _recordRead<Id>(is, "Space Dimension", _nDim);
  ret      = ret && _recordReadVec<double>(is, "Minimum Extension", _extendMin, _nDim);
  ret      = ret && _recordReadVec<double>(is, "Maximum Extension", _extendMax, _nDim);
  return ret;
}

bool AMesh::_serializeAscii(std::ostream& os) const
{
  bool ret = true;
  ret      = ret && _recordWrite<Id>(os, "Space Dimension", getNDim());
  ret      = ret && _recordWriteVec<double>(os, "Minimum Extension", _extendMin);
  ret      = ret && _recordWriteVec<double>(os, "Maximum Extension", _extendMax);
  return ret;
}

/****************************************************************************/
/*!
**  Check if a point, defined by its coordinates, belongs to a Mesh
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
bool AMesh::_weightsInMesh(const VectorDouble& coor,
                           const VectorVectorDouble& corners,
                           double meshsize,
                           VectorDouble& weights,
                           double eps) const
{
  static double FACDIM[] = {0., 1., 2., 6.};

  // Initializations
  auto ncorner = getNApexPerMesh();
  auto ndim    = getNDim();

  // Loop on the vertices
  double total = 0.;
  for (Id icorn = 0; icorn < ncorner; icorn++)
  {

    // Build the determinant
    MatrixSquare mat(ndim);
    Id kcorn = 0;
    for (Id jcorn = 0; jcorn < ncorner; jcorn++)
    {
      if (icorn == jcorn) continue;
      for (Id idim = 0; idim < ndim; idim++)
        mat.setValue(idim, kcorn, corners[jcorn][idim] - coor[idim]);
      kcorn++;
    }
    double ratio = ABS(mat.determinant()) / meshsize / FACDIM[ndim];
    if (ratio < -eps || ratio > 1 + eps) return false;
    weights[icorn] = ratio;
    total += ratio;
  }
  return (ABS(total - 1) <= eps);
}

/****************************************************************************/
/*!
** Returns the size of the Mesh 'imesh'
**
** \returns mesh dimension
**
** \param[in]  corners   Vector of coordinates of mesh apices
**
*****************************************************************************/
double AMesh::_getMeshUnit(const VectorVectorDouble& corners) const
{
  double unit;
  static double facdim[] = {0., 1., 2., 6.};

  // Initializations
  auto ndim    = getNDim();
  auto ncorner = getNApexPerMesh();

  // Calculate the mesh size
  MatrixSquare mat;
  mat.reset(ndim, ndim);
  for (Id icorn = 1; icorn < ncorner; icorn++)
    for (Id idim = 0; idim < ndim; idim++)
      mat.setValue(icorn - 1, idim, corners[icorn][idim] - corners[0][idim]);
  unit = ABS(mat.determinant()) / facdim[ndim];

  return unit;
}

void AMesh::_printMeshListByCoordinates(Id nline_max) const
{
  auto ndim    = getNDim();
  auto nmesh   = getNMeshes();
  auto ncorner = getNApexPerMesh();

  Id iline = 0;
  for (Id imesh = 0; imesh < nmesh; imesh++)
  {
    message("Mesh #%5d/%5d\n", imesh + 1, nmesh);
    for (Id icorn = 0; icorn < ncorner; icorn++)
    {
      message(" Apex %4d: ", getApex(imesh, icorn));
      for (Id idim = 0; idim < ndim; idim++)
        message(" %lf", getCoor(imesh, icorn, idim));
      message("\n");
    }

    iline++;
    if (nline_max > 0 && iline >= nline_max) return;
  }
}

void AMesh::_printMeshListByIndices(Id nline_max) const
{
  auto nmesh   = getNMeshes();
  auto ncorner = getNApexPerMesh();

  Id iline = 0;
  for (Id imesh = 0; imesh < nmesh; imesh++)
  {
    message("Mesh #%d/%d: ", imesh + 1, nmesh);
    for (Id icorn = 0; icorn < ncorner; icorn++)
      message(" %d", getApex(imesh, icorn));
    message("\n");

    iline++;
    if (nline_max > 0 && iline >= nline_max) return;
  }
}

void AMesh::getBarycenterInPlace(Id imesh, vect coord) const
{
  auto ndim    = getNDim();
  auto ncorner = getNApexPerMesh();

  VectorVectorDouble coords = getCoordinatesPerMesh(imesh);

  for (Id idim = 0; idim < ndim; idim++)
  {
    double local = 0.;
    for (Id ic = 0; ic < ncorner; ic++)
      local += coords[ic][idim];
    coord[idim] = local / ncorner;
  }
}

/****************************************************************************/
/*!
** Returns the Sparse Matrix used to project a Db onto the Meshing
**
** \param[out] m         Projection matrix to be initialized
** \param[in]  db        Db structure
** \param[in]  rankZ     Rank of the Z-locator to be tested (see remarks)
** \param[in]  verbose   Verbose flag
**
** \remarks If rankZ>=0, a sample is only considered if the value
** \remarks of the corresponding variable is defined
**
*****************************************************************************/
void AMesh::resetProjFromDb(ProjMatrix* m,
                            const Db* db,
                            Id rankZ,
                            bool verbose) const
{
  auto ndim    = getNDim();
  auto nvertex = getNApices();
  auto ncorner = getNApexPerMesh();
  Id nech      = db->getNSample();
  auto nmeshes = getNMeshes();
  VectorDouble units(nmeshes, 0.);
  if (getVariety() != 1)
    units = _defineUnits();

  // Preliminary checks
  if (isCompatibleDb(db)) return;

  /* Instantiate a Sparse matrix structure (Triplets) */
  NF_Triplet NF_T;

  /* Optional title */
  if (verbose) mestitle(0, "Mesh Barycenter");

  /* Loop on the samples */
  Id ip_max = 0;
  Id iech   = 0;
  Id nout   = 0;
  Id nvalid = 0;
  VectorInt neighs;
  VectorDouble distances;
  VectorDouble target(ndim);
  VectorDouble weight(ncorner, 0);
  for (Id jech = 0; jech < nech; jech++)
  {
    if (!db->isActive(jech)) continue;
    if (rankZ >= 0)
    {
      double testval = db->getFromLocator(ELoc::Z, jech, rankZ);
      if (FFFF(testval)) continue;
    }
    nvalid++;

    // Identification of the target point
    db->getCoordinatesInPlace(target, jech);

    Id found = getMeshAndInPlaceWeightsFromCoordinates(target, weight);

    if (found >= 0)
    {
      /* Store the items in the sparse matrix */

      if (verbose) message("Sample %4d in Mesh %4d :", jech + 1, found + 1);
      for (Id icorn = 0; icorn < ncorner; icorn++)
      {
        auto ip = getApex(found, icorn);
        ip_max = MAX(ip_max, ip);
        if (verbose) message(" %4d (%4.2lf)", ip, weight[icorn]);
        NF_T.add(iech, ip, weight[icorn]);
      }
      if (verbose) message("\n");
    }
    else
    {
      nout++;
    }
    iech++;
  }

  /* Add the extreme value to force dimension */

  if (ip_max < nvertex - 1)
  {
    NF_T.force(nvalid, nvertex);
  }

  /* Convert the triplet into a sparse matrix */

  if (verbose && nout > 0)
    messerr("%d / %d samples which do not belong to the Meshing",
            nout, db->getNSample(true));

  return m->resetFromTriplet(NF_T);
}

Id AMesh::_findBarycenter(const VectorDouble& target,
                          const VectorDouble& units,
                          Id nb_neigh,
                          VectorInt& neighs,
                          VectorDouble& weight) const
{
  for (Id jm = 0; jm < nb_neigh; jm++)
  {
    Id im                      = neighs[jm];
    VectorVectorDouble corners = getCoordinatesPerMesh(im);
    if (!_weightsInMesh(target, corners, units[im], weight)) continue;
    return im;
  }
  return -1;
}

VectorDouble AMesh::_defineUnits(void) const
{
  auto nmeshes = getNMeshes();
  VectorDouble units(nmeshes);
  for (Id imesh = 0; imesh < nmeshes; imesh++)
    units[imesh] = getMeshSize(imesh);
  return units;
}
#ifdef HDF5
bool AMesh::deserializeH5(H5::Group& grp)
{
  auto ameshG = SerializeHDF5::getGroup(grp, "AMesh");
  if (!ameshG)
  {
    return false;
  }

  /* Read the grid characteristics */
  bool ret = true;
  ret      = ret && SerializeHDF5::readValue(*ameshG, "NDim", _nDim);
  ret      = ret && SerializeHDF5::readVec(*ameshG, "ExtendMin", _extendMin);
  ret      = ret && SerializeHDF5::readVec(*ameshG, "ExtendMax", _extendMax);

  return ret;
}

bool AMesh::serializeH5(H5::Group& grp) const
{
  auto ameshG = grp.createGroup("AMesh");

  bool ret = true;
  ret      = ret && SerializeHDF5::writeValue(ameshG, "NDim", getNDim());
  ret      = ret && SerializeHDF5::writeVec(ameshG, "ExtendMin", _extendMin);
  ret      = ret && SerializeHDF5::writeVec(ameshG, "ExtendMax", _extendMax);

  return ret;
}
#endif

} // namespace gstlrn
