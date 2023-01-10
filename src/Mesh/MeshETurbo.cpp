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
#include "geoslib_old_f.h"

#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Mesh/AMesh.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/AException.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/Grid.hpp"
#include "Geometry/Rotation.hpp"
#include "csparse_f.h"

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
                       const VectorDouble& rotmat,
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
  (void) initFromGrid(nx, dx, x0, rotmat, VectorDouble(), flag_polarized, verbose);
}

MeshETurbo::MeshETurbo(const DbGrid* dbgrid, bool verbose, int mode)
    : AMesh(),
      _grid(),
      _nPerCell(0),
      _isPolarized(false),
      _meshIndirect(mode),
      _gridIndirect(mode)
{
  if (!dbgrid->isGrid()) return;
  VectorDouble sel = dbgrid->getSelection();
  (void) initFromGrid(dbgrid->getNXs(), dbgrid->getDXs(), dbgrid->getX0s(),
                      dbgrid->getRotMat(), sel, true, verbose);
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
  if (_gridIndirect.isDefined())
  {
    return _gridIndirect.getRelSize();
  }
  else
  {
    return _grid.getNTotal();
  }
}

/**
 * Returns the total number of possible meshes built using the whole grid
 * (not accounting for possible mask on triangles)
 * @return
 */
int MeshETurbo::_nmeshInCompleteGrid() const
{
  int ndim  = getNDim();

  int nmesh = 1;
  for (int idim=0; idim<ndim; idim++)
    nmesh *= (_grid.getNX(idim) - 1);
  nmesh *= _nPerCell; _meshIndirect.getRelSize();
  return nmesh;
}

/**
 * Actual number of (active) meshes
 * @return
 */
int MeshETurbo::getNMeshes() const
{
  if (_meshIndirect.isDefined())
  {
    return _meshIndirect.getRelSize();
  }
  else
  {
    return _nmeshInCompleteGrid();
  }
}

double MeshETurbo::getMeshSize(int /*imesh*/) const
{
  int ndim = getNDim();

  double size = 1.;
  for (int idim=0; idim<ndim; idim++)
    size *= _grid.getDX(idim);
  size /= _nPerCell;
  return size;
}

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
  VectorInt indg(ndim);

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
  VectorInt indg(getNDim());

  int irel = getApex(imesh,rank);
  int iabs = _gridIndirect.getRToA(irel);
  _grid.rankToIndice(iabs, indg);
  return _grid.indiceToCoordinate(idim, indg, VectorDouble());
}

void MeshETurbo::getCoordinatesInPlace(int imesh, int rank, VectorDouble& coords) const
{
  VectorInt indg(getNDim());

  int irel = getApex(imesh,rank);
  int iabs = _gridIndirect.getRToA(irel);
  _grid.rankToIndice(iabs, indg);

  for (int idim = 0; idim < getNDim(); idim++)
    coords[idim] = _grid.indiceToCoordinate(idim, indg, VectorDouble());
}

double MeshETurbo::getApexCoor(int i, int idim) const
{ _meshIndirect.getRelSize();
  VectorInt indg(getNDim());

  int iabs = _gridIndirect.getRToA(i);
  _grid.rankToIndice(iabs, indg);
  return _grid.indiceToCoordinate(idim, indg, VectorDouble());
}

void MeshETurbo::getApexCoordinatesInPlace(int i, VectorDouble& coords) const
{
  VectorInt indg(getNDim());

  int iabs = _gridIndirect.getRToA(i);
  _grid.rankToIndice(iabs, indg);

  for (int idim = 0; idim < getNDim(); idim++)
    coords[idim] = _grid.indiceToCoordinate(idim, indg, VectorDouble());
}

int MeshETurbo::initFromGrid(const VectorInt&    nx,
                             const VectorDouble& dx,
                             const VectorDouble& x0,
                             const VectorDouble& rotmat,
                             const VectorDouble& sel,
                             bool  flag_polarized,
                             bool  verbose)
{
  int ndim = static_cast<int> (nx.size());
  _setNDim(ndim);

  /* Create the internal (rotated) grid */

  if (_grid.resetFromVector(nx, dx, x0)) return 1;
  _grid.setRotationByVector(rotmat);

  // Get grid extension
  // TODO: the grid extension should be calculated in Grid and take
  // case of a possible rotation
  VectorDouble extendmin(ndim);
  VectorDouble extendmax(ndim);
  for (int idim = 0; idim < ndim; idim++)
  {
    extendmin[idim] = _grid.getX0(idim);
    extendmax[idim] = _grid.getX0(idim) + (_grid.getNX(idim) - 1) * _grid.getDX(idim);
  }
  if (_setExtend(extendmin,extendmax)) return(1);

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
  VectorInt indg(ndim);

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

  return;
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

MeshETurbo* MeshETurbo::create(const VectorInt &nx,
                               const VectorDouble &dx,
                               const VectorDouble &x0,
                               const VectorDouble &rotmat,
                               bool flag_polarized,
                               bool verbose)
{
  MeshETurbo* mesh = new MeshETurbo(nx, dx, x0, rotmat, flag_polarized, verbose);
  return mesh;
}

MeshETurbo* MeshETurbo::createFromGrid(const DbGrid* dbgrid, bool verbose, int mode)
{
  MeshETurbo* mesh = new MeshETurbo(dbgrid, verbose, mode);
  return mesh;
}

MeshETurbo* MeshETurbo::createFromGridInfo(const Grid* grid, bool verbose, int mode)
{
  MeshETurbo* mesh = new MeshETurbo(mode);
  if (mesh->initFromGrid(grid->getNXs(), grid->getDXs(), grid->getX0s(),
                         grid->getRotMat(), VectorDouble(), true, verbose))
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

bool MeshETurbo::_addElementToCS(cs *Atriplet,
                                 int iech,
                                 const VectorDouble &coor,
                                 const VectorInt &indg0,
                                 bool verbose) const
{
  int ncorner = getNApexPerMesh();
  VectorInt indices(ncorner);
  VectorDouble lambda(ncorner);

  for (int icas = 0; icas < _nPerCell; icas++)
  {
    if (_addWeights(icas, indg0, coor, indices, lambda, verbose) == 0)
    {
      for (int icorner = 0; icorner < ncorner; icorner++)
      {
        if (!cs_entry(Atriplet, iech, indices[icorner], lambda[icorner]))
          return false;
      }
      return true;
    }
  }
  return false;
}

/****************************************************************************/
/*!
** Returns the Sparse Matrix used to project a Db onto the Meshing
**
** \return A Sparse matrix (cs structure)
**
** \param[in]  db        Db structure
** \param[in]  verbose   Verbose flag
**
*****************************************************************************/
cs* MeshETurbo::getMeshToDb(const Db *db, bool verbose) const
{
  cs* Atriplet = nullptr;
  cs* A        = nullptr;
  int ndim     = getNDim();
  VectorInt indg0(ndim);
  VectorDouble coor(ndim);

  // Preliminary checks

  if (isCompatibleDb(db)) return NULL;

  // Core allocation

  Atriplet = cs_spalloc(0, 0, 1, 1, 1);
  if (Atriplet == nullptr) return NULL;

  /* Optional title */

  if (verbose) mestitle(0,"Mesh Barycenter");

  /* Loop on the samples */

  int iech = 0;
  int nout = 0;
  for (int jech=0; jech<db->getSampleNumber(); jech++)
  {
    if (! db->isActive(jech)) continue;

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

    bool found = _addElementToCS(Atriplet, iech, coor, indg0, verbose);

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
        found = _addElementToCS(Atriplet, iech, coor, indg0, verbose);
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

  cs_force_dimension(Atriplet, db->getSampleNumber(true), getNApices());

  /* Convert the triplet into a sparse matrix */

  A = cs_triplet(Atriplet);

  // Set the error return code

  if (verbose && nout > 0)
    messerr("%d / %d samples which do not belong to the Meshing",
            nout, db->getSampleNumber(true));

  Atriplet  = cs_spfree(Atriplet);
  return(A);
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
  _grid.display();
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
int MeshETurbo::_addWeights(int icas,
                            const VectorInt& indg0,
                            const VectorDouble& coor,
                            VectorInt& indices, // Returned indices (active grid nodes)
                            VectorDouble& lambda,
                            bool verbose) const
{
  int ndim    =  getNDim();
  int ncorner =  getNApexPerMesh();
  int ipol    = _getPolarized(indg0);
  MatrixSquareGeneral lhs;
  VectorDouble rhs(ncorner);
  VectorInt indgg(ndim);

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
      lhs.setValue(icorner,idim,_grid.indiceToCoordinate(idim,indgg));
    lhs.setValue(icorner,ndim,1.);
  }

  // Generate the right-hand side vector
  for (int idim=0; idim<ndim; idim++)
    rhs[idim] = coor[idim];
  rhs[ndim] = 1;

  // Invert the matrix
  if (lhs.invert()) return 1;

  // Calculate the weights
  lhs.prodVector(rhs,lambda);

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

int MeshETurbo::_getPolarized(VectorInt indg) const
{
  int ndim = getNDim();
  if (! _isPolarized) return(0);

  // Polarization has only been coded for the 2-D case
  if (ndim != 2) return(0);
  if ((indg[0] + indg[1]) % 2 == 1)
    return(0);
  else
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
  VectorInt indg(getNDim());
  int ncas = _nPerCell;
  int rank = imesh / ncas;
  *icas    = imesh - rank * ncas;
  _grid.rankToIndice(rank,indg,true);
  *node = _grid.indiceToRank(indg);
}

int MeshETurbo::initFromCova(const CovAniso& cova,
                             const DbGrid* field,
                             double ratio,
                             int nbExt,
                             bool /*useSel*/,
                             bool verbose)
{
  // Preliminary checks
  if (! field->isGrid())
  {
    messerr("This function is limited to 'field' defined as a Grid");
    return 1;
  }

  // Initializations
  int ndim = cova.getNDim();
  int nval = (int) pow(2., ndim);

  // Get the rotation linked to the covariance
  const Rotation& rot = cova.getAnisoRotation();

  // Project the corners of the grid
  VectorDouble extendMinNoRot(ndim, TEST);
  VectorDouble extendMaxNoRot(ndim, TEST);
  VectorDouble cornerRot(ndim);
  VectorInt ic(ndim);
  for (int icorner = 0; icorner < nval; icorner++)
  {

    // Get one corner
    int jcorner = icorner;
    for (int idim = 0 ; idim < ndim; idim++)
    {
      ic[idim] = jcorner % 2;
      jcorner /= 2;
    }
    VectorDouble corner1 = field->getGrid().getCoordinatesByCorner(ic);

    // Rotate this corner in the Covariance Rotation system
    rot.rotateDirect(corner1, cornerRot);

    // Calculate the minimum and maximum in the Covariance rotated system
    for (int idim = 0; idim < ndim; idim++)
    {
      if (FFFF(extendMinNoRot[idim]) || cornerRot[idim] < extendMinNoRot[idim])
        extendMinNoRot[idim] = cornerRot[idim];
      if (FFFF(extendMaxNoRot[idim]) || cornerRot[idim] > extendMaxNoRot[idim])
        extendMaxNoRot[idim] = cornerRot[idim];
    }
  }

  // Calculating the Mesh of the Grid
  VectorInt    nx(ndim);
  VectorDouble dx(ndim);
  for (int idim = 0; idim < ndim; idim++)
  {
    double delta = extendMaxNoRot[idim] - extendMinNoRot[idim];
    dx[idim] = cova.getRange(idim) / ratio;
    nx[idim] = (int) ceil(delta / dx[idim]) + 2 * nbExt;
    if( nx[idim] > 400 )
    {
      nx[idim] = 400;
      dx[idim] = delta / (nx[idim] - 2 * nbExt);
    }
    extendMinNoRot[idim] -= nbExt * dx[idim];
    extendMaxNoRot[idim] += nbExt * dx[idim];
  }

  // Get the rotated Bounding Box in the initial system
  VectorDouble x0(ndim);
  rot.rotateInverse(extendMinNoRot, x0);

  initFromGrid(nx,dx,x0,rot.getMatrixDirectByVector(),VectorDouble(),true,verbose);
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
    (void) initFromGrid(nx, dx, x0, rotmat, VectorDouble(), (bool) flag_polarized, 0);

  ret = ret && _recordRead<int>(is, "Mesh Active Count", nmesh_active);
  ret = ret && _recordReadVec<int>(is, "Mesh Masking", relranks, nmesh_active);
  if (ret) _meshIndirect.buildFromRankRInA(relranks, _nmeshInCompleteGrid());

  ret = ret && _recordRead<int>(is, "Grid Active Count", ngrid_active);
  ret = ret && _recordReadVec<int>(is, "Grid Masking", relranks, ngrid_active);
  if (ret) _gridIndirect.buildFromRankRInA(relranks, _grid.getNTotal());

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
  ret = ret && _recordWrite<int>(os, "Mesh Active Count", getNMeshes());
  ret = ret && _recordWriteVec<int>(os, "Mesh Masking", _meshIndirect.getRelRanks());

  // Dumping the Grid Masking map
  ret = ret && _recordWrite<int>(os, "Grid Active Count", getNApices());
  ret = ret && _recordWriteVec<int>(os, "Grid Masking", _gridIndirect.getRelRanks());

  // Dumping the Grid Masking map
  return ret;
}
