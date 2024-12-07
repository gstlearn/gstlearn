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
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/NF_Triplet.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "Mesh/AMesh.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Mesh/Delaunay.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/Grid.hpp"
#include "Geometry/Rotation.hpp"

#include <math.h>

MeshETurbo::MeshETurbo(int mode)
    : AMesh(),
      _grid(),
      _nPerCell(0),
      _isPolarized(false),
      _meshIndirect(mode),
      _gridIndirect(mode)
{
}

MeshETurbo::MeshETurbo(const VectorInt& nx,
                       const VectorDouble& dx,
                       const VectorDouble& x0,
                       const VectorDouble& angles,
                       bool flag_polarized,
                       bool verbose,
                       int mode)
    : AMesh(),
      _grid(),
      _nPerCell(0),
      _isPolarized(flag_polarized),
      _meshIndirect(mode),
      _gridIndirect(mode)
{
  (void) initFromGridByAngles(nx, dx, x0, angles, VectorDouble(), flag_polarized, verbose);
}

MeshETurbo::MeshETurbo(const DbGrid *dbgrid,
                       bool flag_polarized,
                       bool verbose,
                       int mode)
    : AMesh(),
      _grid(),
      _nPerCell(0),
      _isPolarized(flag_polarized),
      _meshIndirect(mode),
      _gridIndirect(mode)
{
  if (!dbgrid->isGrid()) return;
  VectorDouble sel = dbgrid->getSelections();
  (void)initFromGridByMatrix(dbgrid->getNXs(), dbgrid->getDXs(),
                             dbgrid->getX0s(), dbgrid->getRotMat(), sel,
                             flag_polarized, verbose);
}

MeshETurbo::MeshETurbo(const MeshETurbo &r)
    : AMesh(r),
      _grid(),
      _nPerCell(0),
      _isPolarized(false),
      _meshIndirect(r._meshIndirect),
      _gridIndirect(r._gridIndirect)
{
  _grid = r._grid;
}

MeshETurbo& MeshETurbo::operator= (const MeshETurbo &r)
{
  _grid = r._grid;
  _nPerCell = r._nPerCell;
  _isPolarized = r._isPolarized;
  _meshIndirect = r._meshIndirect;
  _gridIndirect = r._gridIndirect;

  return *this;
}

MeshETurbo::~MeshETurbo()
{
}

/**
 * Returns the total number of apices of the whole grid
 * (not accounting for possible mask on meshes)
 * @return
 */
int MeshETurbo::getNApices() const
{
  if (_gridIndirect.isDefined()) return _gridIndirect.getRelSize();
  return _grid.getNTotal();
}

/**
 * Returns the total number of possible meshes built using the whole grid
 * (not accounting for possible mask on triangles)
 * @return
 */
int MeshETurbo::_nmeshInCompleteGrid() const
{
  int nmesh = 1;
  for (int idim = 0, ndim = getNDim(); idim < ndim; idim++)
    nmesh *= (_grid.getNX(idim) - 1);
  nmesh *= _nPerCell;
  // _meshIndirect.getRelSize();
  return nmesh;
}

/**
 * Actual number of (active) meshes
 * @return
 */
int MeshETurbo::getNMeshes() const
{
  if (_meshIndirect.isDefined()) return _meshIndirect.getRelSize();
  return _nmeshInCompleteGrid();
}

double MeshETurbo::getMeshSize(int /*imesh*/) const
{
  double size = 1.;
  for (int idim=0, ndim=getNDim(); idim<ndim; idim++)
    size *= _grid.getDX(idim);
  size /= _nPerCell;
  return size;
}

static std::vector<int> indg;

/****************************************************************************/
/*!
** Returns the Apex 'rank' of the Mesh 'imesh'
**
** \returns The rank of the target apex
**
** \param[in]  imesh    Rank of active Mesh (starting from 0)
** \param[in]  rank     Rank of Apex within a Mesh (from 0 to _nApexPerMesh-1)
**
*****************************************************************************/
int MeshETurbo::getApex(int imesh, int rank) const
{
  int node,icas;
  int ndim = getNDim();
  indg.resize(ndim);

  int jmesh = _meshIndirect.getRToA(imesh);
  _getGridFromMesh(jmesh,&node,&icas);
  _grid.rankToIndice(node,indg);
  int ipol = _getPolarized(indg);

  for (int idim = 0; idim < ndim; idim++)
    indg[idim] += MSS(ndim, ipol, icas, rank, idim);

  int igrid = _grid.indiceToRank(indg);

  int irel = _gridIndirect.getAToR(igrid);
  if (irel < 0)
  {
    messerr("Problem for mesh=%d rank=%d grid=%d -> Mesh relative rank is negative",
            imesh+1,rank+1,igrid+1);
  }
  return irel;
}

double MeshETurbo::getCoor(int imesh, int rank, int idim) const
{
  indg.resize(getNDim());

  int irel = getApex(imesh,rank);
  int iabs = _gridIndirect.getRToA(irel);
  _grid.rankToIndice(iabs, indg);
  return _grid.indiceToCoordinate(idim, indg);
}

void MeshETurbo::getCoordinatesInPlace(int imesh, int rank, VectorDouble& coords) const
{
  indg.resize(getNDim());

  int irel = getApex(imesh,rank);
  int iabs = _gridIndirect.getRToA(irel);
  _grid.rankToCoordinatesInPlace(iabs, coords);
}

double MeshETurbo::getApexCoor(int i, int idim) const
{
  // _meshIndirect.getRelSize();
  indg.resize(getNDim());

  int iabs = _gridIndirect.getRToA(i);
  _grid.rankToIndice(iabs, indg);
  return _grid.indiceToCoordinate(idim, indg);
}

void MeshETurbo::getApexIndicesInPlace(int i, VectorInt& indg) const
{
  // _meshIndirect.getRelSize();
  indg.resize(getNDim());

  int iabs = _gridIndirect.getRToA(i);
  _grid.rankToIndice(iabs, indg);
}

void MeshETurbo::getApexCoordinatesInPlace(int i, VectorDouble& coords) const
{
  indg.resize(getNDim());

  int iabs = _gridIndirect.getRToA(i);
  _grid.rankToIndice(iabs, indg);

  for (int idim = 0; idim < getNDim(); idim++)
    coords[idim] = _grid.indiceToCoordinate(idim, indg);
}

int MeshETurbo::initFromGridByAngles(const VectorInt& nx,
                                     const VectorDouble& dx,
                                     const VectorDouble& x0,
                                     const VectorDouble& angles,
                                     const VectorDouble& sel,
                                     bool flag_polarized,
                                     bool verbose)
{
  int ndim = static_cast<int>(nx.size());
  _setNDim(ndim);

  /* Create the internal (rotated) grid by Angles */

  if (_grid.resetFromVector(nx, dx, x0)) return 1;
  _grid.setRotationByAngles(angles);

  return _initFromGridInternal(sel, flag_polarized, verbose);
}

int MeshETurbo::initFromGridByMatrix(const VectorInt& nx,
                                     const VectorDouble& dx,
                                     const VectorDouble& x0,
                                     const VectorDouble& rotmat,
                                     const VectorDouble& sel,
                                     bool flag_polarized,
                                     bool verbose)
{
  int ndim = static_cast<int> (nx.size());
  _setNDim(ndim);

  /* Create the internal (rotated) grid by Rotation Matrix */

  if (_grid.resetFromVector(nx, dx, x0)) return 1;
  _grid.setRotationByVector(rotmat);

  return _initFromGridInternal(sel, flag_polarized, verbose);
}

int MeshETurbo::_initFromGridInternal(const VectorDouble& sel,
                                      bool flag_polarized,
                                      bool verbose)
{
  int ndim = getNDim();

  // Get grid extension
  // TODO: the grid extension should be calculated in Grid and take
  // case of a possible rotation
  VectorDouble extendmin(ndim);
  VectorDouble extendmax(ndim);
  for (int idim = 0; idim < ndim; idim++)
  {
    extendmin[idim] = _grid.getX0(idim);
    extendmax[idim] =
      _grid.getX0(idim) + (_grid.getNX(idim) - 1) * _grid.getDX(idim);
  }
  if (_setExtend(extendmin, extendmax)) return (1);

  // Define the number of Elements per Cell

  _setNumberElementPerCell();

  // Set polarization

  _isPolarized = flag_polarized;

  // Convert optional selection array into mask on meshes

  _buildMaskInMeshing(sel);

  // Optional printout

  if (verbose) messageFlush(toString());

  return 0;
}

void MeshETurbo::_buildMaskInMeshing(const VectorDouble& sel)
{
  int node, icas;
  std::map<int,int> map;

  // If no selection is defined on the grid, the vector of Meshing Mask is cancelled
  if (sel.empty()) return;

  // Creating the Masking information for Meshing 'meshActiveToAbsolute'
  // which gives the Absolute meshing index from its Active index

  int ndim = getNDim();
  int nmesh = _nmeshInCompleteGrid();
  int ncorner = getNApexPerMesh();
  VectorInt indg0(ndim);
  indg.resize(ndim);

  // Loop on all possible meshes
  int meshNactive = 0;
  for (int imesh = 0; imesh < nmesh; imesh++)
  {
    _getGridFromMesh(imesh,&node,&icas);
    _grid.rankToIndice(node,indg0);
    int ipol = _getPolarized(indg0);

    // Loop on the corners of the mesh (polarization is taken into account)
    bool flagMasked = false;
    for (int icorner=0; icorner<ncorner && !flagMasked; icorner++)
    {

      // Generate the indices of the mesh apex
      for (int idim = 0; idim < ndim; idim++)
        indg[idim] = indg0[idim] + MSS(ndim, ipol, icas, icorner, idim);
      int iad = _grid.indiceToRank(indg);
      if (sel[iad] == 0.) flagMasked = true;
    }

    // The triangle is not masked, store its index
    if (flagMasked) continue;
    map[imesh] = meshNactive;
    meshNactive++;
  }

  // Creating the Indirection information for Meshing

  _meshIndirect.buildFromMap(map, nmesh);

  // Creating the selection of the active grid nodes
  // It is at least equal to 'sel'. In addition, non active grid nodes are also discarded

  VectorDouble selbis = VectorDouble(sel.size(),0.);
  for (int imesh = 0; imesh < meshNactive; imesh++)
  {
    int jmesh = _meshIndirect.getRToA(imesh);
    _getGridFromMesh(jmesh,&node,&icas);
    _grid.rankToIndice(node,indg0);
    int ipol = _getPolarized(indg0);
    for (int icorner=0; icorner<ncorner; icorner++)
    {
      for (int idim = 0; idim < ndim; idim++)
        indg[idim] = indg0[idim] + MSS(ndim, ipol, icas, icorner, idim);
      int iad = _grid.indiceToRank(indg);
      selbis[iad] = 1;
    }
  }

  _gridIndirect.buildFromSel(selbis);
}

/**
 * Create a MeshETurbo by loading the contents of a Neutral File
 *
 * @param neutralFilename Name of the Neutral File (MeshEStandard format)
 * @param verbose         Verbose
 */
MeshETurbo* MeshETurbo::createFromNF(const String& neutralFilename, bool verbose)
{
  MeshETurbo* mesh = nullptr;
  std::ifstream is;
  mesh = new MeshETurbo;
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

MeshETurbo* MeshETurbo::create(const VectorInt& nx,
                               const VectorDouble& dx,
                               const VectorDouble& x0,
                               const VectorDouble& angles,
                               bool flag_polarized,
                               bool verbose)
{
  MeshETurbo* mesh = new MeshETurbo(nx, dx, x0, angles, flag_polarized, verbose);
  return mesh;
}

MeshETurbo* MeshETurbo::createFromGrid(const DbGrid *dbgrid,
                                       bool flag_polarized,
                                       bool verbose,
                                       int mode)
{
  MeshETurbo* mesh = new MeshETurbo(dbgrid, flag_polarized, verbose, mode);
  return mesh;
}

MeshETurbo* MeshETurbo::createFromGridInfo(const Grid *grid,
                                           bool flag_polarized,
                                           bool verbose,
                                           int mode)
{
  MeshETurbo* mesh = new MeshETurbo(mode);
  if (mesh->initFromGridByMatrix(grid->getNXs(), grid->getDXs(), grid->getX0s(),
                         grid->getRotMat(), VectorDouble(), flag_polarized,
                         verbose))
    return nullptr;
  return mesh;
}

MeshETurbo* MeshETurbo::createFromCova(const CovAniso &cova,
                                       const Db *field,
                                       double ratio,
                                       int nbExt,
                                       bool useSel,
                                       bool flagNoStatRot,
                                       bool verbose)
{
  MeshETurbo* mesh = new MeshETurbo();
  if (mesh->initFromCova(cova, field, ratio, nbExt, useSel, flagNoStatRot, verbose))
    return nullptr;
  return mesh;
}

/****************************************************************************/
/*!
** Create the meshing
**
** \param[in]  extendmin       Minimum of the dilated rotated bounding box
** \param[in]  extendmax       Minimum of the dilated rotated bounding box
** \param[in]  cellsize        Array giving the cell size (see details)
** \param[in]  rotmat          Rotation matrix (optional)
** \param[in]  flag_polarized  Switching ON/OFF the polarization
** \param[in]  verbose         Verbose flag
**
*****************************************************************************/
int MeshETurbo::initFromExtend(const VectorDouble &extendmin,
                               const VectorDouble &extendmax,
                               const VectorDouble &cellsize,
                               const VectorDouble &rotmat,
                               bool flag_polarized,
                               bool verbose)
{
  int ndim = static_cast<int> (extendmin.size());
  _setNDim(ndim);
  if (_setExtend(extendmin,extendmax)) return(1);

  /* Create the internal (rotated) grid */

  if (_defineGrid(cellsize)) return(1);
  _grid.setRotationByVector(rotmat);

  // Define the number of Elements per Cell

  _setNumberElementPerCell();

  // Set polarization

  _isPolarized = flag_polarized;

  // Optional printout

  if (verbose) messageFlush(toString());

  return 0;
}

static std::vector<int> indices;
static std::vector<double> lambda;

bool MeshETurbo::_addElementToTriplet(NF_Triplet& NF_T,
                                      int iech,
                                      const VectorDouble &coor,
                                      const VectorInt &indg0,
                                      bool verbose) const
{
  int ncorner = getNApexPerMesh();
  indices.resize(ncorner);
  lambda.resize(ncorner);

  for (int icas = 0; icas < _nPerCell; icas++)
  {
    if (_addWeights(icas, indg0, coor, indices, lambda, verbose) == 0)
    {
      for (int icorner = 0; icorner < ncorner; icorner++)
        NF_T.add(iech, indices[icorner], lambda[icorner]);
      return true;
    }
  }
  return false;
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
void MeshETurbo::resetProjMatrix(ProjMatrix* m, const Db *db, int rankZ, bool verbose) const
{
  int ndim = getNDim();
  VectorInt indg0(ndim);
  VectorDouble coor(ndim);

  // Preliminary checks

  if (isCompatibleDb(db)) return;

  // Core allocation

  NF_Triplet NF_T;

  /* Optional title */

  if (verbose) mestitle(0,"Mesh Barycenter");

  /* Loop on the samples */

  int iech = 0;
  int nout = 0;
  int nvalid = 0;
  for (int jech=0; jech<db->getSampleNumber(); jech++)
  {
    if (! db->isActive(jech)) continue;
    if (rankZ >= 0)
    {
      double testval = db->getFromLocator(ELoc::Z, jech, rankZ);
      if (FFFF(testval)) continue;
    }
    nvalid++;

    /* Identification */

    for (int idim=0; idim<ndim; idim++)
      coor[idim] = db->getCoordinate(jech,idim);

    /* Calculate the grid indices */

    if (_grid.coordinateToIndicesInPlace(coor,indg0) != 0)
    {
      messerr("Sample #%d does not belong to the grid",jech+1);
      continue;
    }

    /* Optional printout */

    if (verbose)
      message("Sample %4d assigned to Grid Node %4d :",
              jech+1,_grid.indiceToRank(indg0)+1);

    // Finding the active mesh to which the sample belongs

    bool found = _addElementToTriplet(NF_T, iech, coor, indg0, verbose);

    // In the case the target coordinate is on the edge of the grid
    // try to shift the point down by one node
    if (! found)
    {
      bool flag_correct = false;
      for (int idim = 0; idim < ndim; idim++)
      {
        if (indg0[idim] != _grid.getNX(idim)-1) continue;
        indg0[idim] -= 1;
        flag_correct = true;
      }
      if (flag_correct)
        found = _addElementToTriplet(NF_T, iech, coor, indg0, verbose);
    }

    // The point does not belong to any active mesh, issue a message (optional)
    if (! found)
    {
      nout++;
      if (verbose)
        messerr("Sample #%d does not belong to the meshing",jech+1);
    }
    iech++;
  }

  /* Add the extreme value to force dimension */

  NF_T.force(nvalid, getNApices());

  /* Convert the triplet into a sparse matrix */

  if (verbose && nout > 0)
    messerr("%d / %d samples which do not belong to the Meshing",
            nout, db->getSampleNumber(true));

  return m->resetFromTriplet(NF_T);
}

/****************************************************************************/
/*!
** Print the contents of the meshing
**
** \param[in] strfmt    Format for printout
**
*****************************************************************************/
String MeshETurbo::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  sstr << toTitle(0,"Turbo Meshing");
  if (_isPolarized) sstr << "Diamond construction is activated" << std::endl;
  sstr << _grid.toString(strfmt);
  sstr << AMesh::toString(strfmt);

  if (_meshIndirect.isDefined())
  {
    sstr << toTitle(2, "Mask Information");
    sstr << "Mesh Masking Indexing" << std::endl;
    sstr << _meshIndirect.toString(strfmt) << std::endl;
    sstr << "Grid Masking Indexing" << std::endl;
    sstr << _gridIndirect.toString(strfmt) << std::endl;
  }

  return sstr.str();
}

/****************************************************************************/
/*!
** Define the internal grid
**
** \param[in]  cellsize  Array giving the cell size (see details)
**
*****************************************************************************/
int MeshETurbo::_defineGrid(const VectorDouble& cellsize)

{
  int ndim;

  // Initializations

  if (cellsize.empty())
  {
    messerr("The argument 'cellsize' must be provided");
    return(1);
  }
  ndim = getNDim();

  // Create the grid internal structure

  _grid.resetFromSpaceDimension(ndim);

  // Copy the grid main characteristics

  for (int idim=0; idim<ndim; idim++)
  {
    _grid.setX0(idim, getExtendMin(idim));
    _grid.setDX(idim, cellsize[idim]);
    _grid.setNX(idim, static_cast<int> (ceil((getExtendMax(idim) - getExtendMin(idim)) /
                           cellsize[idim]) + 1));
  }

  return 0;
}

void MeshETurbo::_setNumberElementPerCell()
{
  int ndim = getNDim();

  if (ndim == 1)
    _nPerCell = 1;
  else if (ndim == 2)
    _nPerCell = 2;
  else if (ndim == 3)
    _nPerCell = 6;
}

static std::vector<double> rhs;
static std::vector<int> indgg;

/**
 * Return the weights assigned to the corners
 * @param icas   Corner indication
 * @param indg0  Indices of the starting grid node
 * @param coor   Coordinates of the targte point
 * @param indices Grid indices of the target (in active ranks)
 * @param lambda  Weights
 * @param verbose Verbose flag
 * @return
 *
 * @remark The function returns 1 if:
 * @remark - the grid node corresponding to a mesh apex is outside the grid
 * @remark - the grid node corresponding to a mesh apex is not active
 */
int MeshETurbo::_addWeights(
  int icas,
  const constvectint indg0,
  const constvect coor,
  const vectint indices, // Returned indices (active grid nodes)
  const vect lambda,
  bool verbose) const
{
  int ndim    =  getNDim();
  int ncorner =  getNApexPerMesh();
  int ipol    = _getPolarized(indg0);
  MatrixSquareGeneral lhs;
  rhs.resize(ncorner);
  indgg.resize(ndim);

  // Build the LHS matrix

  lhs.reset(ncorner,ncorner);
  for (int icorner=0; icorner<ncorner; icorner++)
  {
    // Generate the indices of the mesh apex
    for (int idim=0; idim<ndim; idim++)
      indgg[idim] = indg0[idim] + MSS(ndim,ipol,icas,icorner,idim);
    int igrid = _grid.indiceToRank(indgg);
    if (igrid < 0) return 1; // grid node outside grid

    indices[icorner] = _gridIndirect.getAToR(igrid);
    if (indices[icorner] < 0) return 1; // grid node not active

    // Update the LHS matrix
    for (int idim=0; idim<ndim; idim++)
      lhs.setValue(idim,icorner,_grid.indiceToCoordinate(idim,indgg));
    lhs.setValue(ndim,icorner,1.);
  }

  // Generate the right-hand side vector
  for (int idim=0; idim<ndim; idim++)
    rhs[idim] = coor[idim];
  rhs[ndim] = 1;

  // Invert the matrix
  if (lhs.invert()) return 1;

  // Calculate the weights
  lhs.prodMatVecInPlace(rhs, lambda);

  // Check that all weights are positive
  for (int icorner=0; icorner<ncorner; icorner++)
  {
    if (lambda[icorner] < -EPSILON6) return 1;
    if (lambda[icorner] < 0) lambda[icorner] = 0.;
    if (lambda[icorner] > 1 + EPSILON6) return 1;
    if (lambda[icorner] > 1) lambda[icorner] = 1.;
  }

  // Optional printout
  if (verbose)
  {
    for (int icorner=0; icorner<ncorner; icorner++)
      message(" %4d (%4.2lf)",indices[icorner],lambda[icorner]);
    message("\n");
  }

  return 0;
}

int MeshETurbo::_getPolarized(const constvectint indg) const
{
  int ndim = getNDim();
  if (! _isPolarized) return(0);

  // Polarization has only been coded for the 2-D case
  if (ndim != 2) return(0);
  if ((indg[0] + indg[1]) % 2 == 1)
    return(0);
  return(1);
}

/**
 * Returns the (starting) grid node, given the absolute rank of the mesh
 * @param imesh Absolute Rank of the  mesh
 * @param node  Starting grid node
 * @param icas  Sorting used for reviewing grid meshes (takes polarization into account)
 */
void MeshETurbo::_getGridFromMesh(int imesh, int *node, int *icas) const
{
  indg.resize(getNDim());
  int ncas = _nPerCell;
  int rank = imesh / ncas;
  *icas    = imesh - rank * ncas;
  _grid.rankToIndice(rank,indg,true);
  *node = _grid.indiceToRank(indg);
}

int MeshETurbo::initFromCova(const CovAniso& cova,
                             const Db* field,
                             double ratio,
                             int nbExt,
                             bool useSel,
                             bool flagNoStatRot,
                             bool verbose)
{
  // Initializations
  int ndim = cova.getNDim();
  int nval = (int) pow(2., ndim);

  // Get the rotation linked to the covariance
  const Rotation& rot = cova.getAnisoRotation();

  // Project the corners of the grid
  VectorDouble extendMinRot(ndim, TEST);
  VectorDouble extendMaxRot(ndim, TEST);
  VectorDouble cornerRot(ndim);
  VectorInt ic(ndim,0);
  VectorVectorDouble extremesData;
  VectorDouble cornerRef(ndim,0.);
  const DbGrid* fieldGrid = nullptr;
  if (field->isGrid())
  {
    fieldGrid = dynamic_cast<const DbGrid*>(field);
    cornerRef = fieldGrid->getGrid().getCoordinatesByCorner(ic);
  }
  else
  {
    cornerRef = field->getCoorMinimum(useSel);
    extremesData = field->getExtremas(useSel);
  }

  for (int icorner = 0; icorner < nval; icorner++)
  {

    // Get one corner
    int jcorner = icorner;
    for (int idim = 0 ; idim < ndim; idim++)
    {
      ic[idim] = jcorner % 2;
      jcorner /= 2;
    }
    VectorDouble corner1(ndim,0.);
    if (field->isGrid())
    {
      corner1 = fieldGrid->getGrid().getCoordinatesByCorner(ic);
    }
    else
    {
      for (int idim = 0; idim < ndim; idim++)
        corner1[idim] = (ic[idim] == 0) ? extremesData[idim][0] : extremesData[idim][1];
    }
    VH::subtractInPlace(corner1, cornerRef);

    // Rotate this corner in the Covariance Rotation system
    rot.rotateInverse(corner1, cornerRot);

    // Calculate the minimum and maximum in the Covariance rotated system
    for (int idim = 0; idim < ndim; idim++)
    {
      if (FFFF(extendMinRot[idim]) || cornerRot[idim] < extendMinRot[idim])
        extendMinRot[idim] = cornerRot[idim];
      if (FFFF(extendMaxRot[idim]) || cornerRot[idim] > extendMaxRot[idim])
        extendMaxRot[idim] = cornerRot[idim];
    }
  }

  // Calculating the Mesh of the Grid
  VectorInt    nx(ndim);
  VectorDouble dx(ndim);
  VectorDouble x0(ndim);

  double dxmin = 1.e30;
  int nxmax = 300;
  for (int idim = 0; idim < ndim; idim++)
  {
    double delta = extendMaxRot[idim] - extendMinRot[idim];
    dx[idim] = cova.getRange(idim) / ratio;
    nx[idim] = (int) ceil(delta / dx[idim]) + 2 * nbExt + 1;
    if (nx[idim] > nxmax)
    {
      nx[idim] = nxmax;
      dx[idim] = delta / (nxmax - 2 * nbExt);
    }
    if (dx[idim] < dxmin) dxmin = dx[idim];
  }

  if (flagNoStatRot)
  {
    // In case of non-staionarity on the anisotropy rotation angle
    // use the minimum mesh (for internal grid)
    for (int idim = 0; idim < ndim; idim++)
      dx[idim] = dxmin;

    for (int idim = 0; idim < ndim; idim++)
    {
      double delta = extendMaxRot[idim] - extendMinRot[idim];
      nx[idim] = (int) ceil(delta / dx[idim]) + 2 * nbExt + 1;
    }
  }

  for (int idim = 0; idim < ndim; idim++)
  {
    extendMinRot[idim] -= nbExt * dx[idim];
    extendMaxRot[idim] += nbExt * dx[idim];
  }

  // Get the rotated Bounding Box in the initial system
  rot.rotateDirect(extendMinRot, x0);
  VH::addInPlace(x0, cornerRef);

  initFromGridByMatrix(nx,dx,x0,rot.getMatrixDirectVec(),VectorDouble(),true,verbose);
  return 0;
}

bool MeshETurbo::_deserialize(std::istream& is, bool /*verbose*/)
{
  int ndim = 0;
  VectorInt nx;
  VectorDouble dx;
  VectorDouble x0;
  VectorDouble rotmat;
  VectorInt relranks;
  int flag_polarized = 0;
  int nmesh_active = 0;
  int ngrid_active = 0;
  int nmesh_mask = 0;
  int ngrid_mask = 0;
  int mode = 0;

  bool ret = true;
  ret = ret && _recordRead<int>(is, "Space Dimension", ndim);
  ret = ret && _recordReadVec<int>(is, "NX", nx, ndim);
  ret = ret && _recordReadVec<double>(is, "DX", dx, ndim);
  ret = ret && _recordReadVec<double>(is, "X0", x0, ndim);
  ret = ret && _recordReadVec<double>(is, "Rotation", rotmat, ndim * ndim);
  ret = ret && _recordRead<int>(is, "Polarization", flag_polarized);
  ret = ret && _recordRead<int>(is, "Storing Mode", mode);

  if (ret)
  {
    _meshIndirect.setMode(mode);
    _gridIndirect.setMode(mode);
  }

  if (ret)
    (void) initFromGridByMatrix(nx, dx, x0, rotmat, VectorDouble(), (bool) flag_polarized, 0);

  ret = ret && _recordRead<int>(is, "Mesh Active Count", nmesh_active);
  ret = ret && _recordRead<int>(is, "Mesh Masking Count", nmesh_mask);
  if (ret && nmesh_mask > 0)
  {
    ret = ret && _recordReadVec<int>(is, "Mesh Masking", relranks, nmesh_active);
    if (ret) _meshIndirect.buildFromRankRInA(relranks, _nmeshInCompleteGrid());
  }

  ret = ret && _recordRead<int>(is, "Grid Active Count", ngrid_active);
  ret = ret && _recordRead<int>(is, "Mesh Masking Count", ngrid_mask);
  if (ret && ngrid_mask > 0)
  {
    ret = ret && _recordReadVec<int>(is, "Grid Masking", relranks, ngrid_active);
    if (ret) _gridIndirect.buildFromRankRInA(relranks, _grid.getNTotal());
  }
  return ret;
}

bool MeshETurbo::_serialize(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;
  ret = ret && _recordWrite<int>(os, "Space Dimension", getNDim());
  ret = ret && _recordWriteVec<int>(os, "NX", _grid.getNXs());
  ret = ret && _recordWriteVec<double>(os, "DX", _grid.getDXs());
  ret = ret && _recordWriteVec<double>(os, "X0", _grid.getX0s());
  ret = ret && _recordWriteVec<double>(os, "Rotation", _grid.getRotMat());
  ret = ret && _recordWrite<int>(os, "Polarization", _isPolarized);
  ret = ret && _recordWrite<int>(os, "Storing Mode", _meshIndirect.getMode());

  // Dumping the Mesh Masking map
  int nmesh_mask = (int) _meshIndirect.getRelRanks().size();
  ret = ret && _recordWrite<int>(os, "Mesh Active Count", getNMeshes());
  ret = ret && _recordWrite<int>(os, "Mesh Masking Count", nmesh_mask);
  if (nmesh_mask > 0)
    ret = ret && _recordWriteVec<int>(os, "Mesh Masking", _meshIndirect.getRelRanks());

  // Dumping the Grid Masking map
  int ngrid_mask = (int) _gridIndirect.getRelRanks().size();
  ret = ret && _recordWrite<int>(os, "Grid Active Count", getNApices());
  ret = ret && _recordWrite<int>(os, "Grid Masking Count", ngrid_mask);
  if (ngrid_mask > 0)
    ret = ret && _recordWriteVec<int>(os, "Grid Masking", _gridIndirect.getRelRanks());

  // Dumping the Grid Masking map
  return ret;
}
