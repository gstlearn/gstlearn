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
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Mesh/AMesh.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/Vector.hpp"
#include "Basic/Rotation.hpp"
#include "Basic/AException.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/Grid.hpp"
#include "csparse_f.h"

#include <math.h>

MeshETurbo::MeshETurbo()
    : AMesh(),
      _grid(),
      _nPerCell(0),
      _isPolarized(true),
      _isMaskDefined(false),
      _maskGrid(nullptr)
{
}

MeshETurbo::MeshETurbo(const VectorInt& nx,
                       const VectorDouble& dx,
                       const VectorDouble& x0,
                       const VectorDouble& rotmat,
                       bool flag_polarized,
                       int verbose)
    : AMesh(),
      _grid(),
      _nPerCell(0),
      _isPolarized(flag_polarized),
      _isMaskDefined(false),
      _maskGrid(nullptr)
{
  (void) initFromGrid(nx, dx, x0, rotmat, flag_polarized, verbose);
}

MeshETurbo::MeshETurbo(const DbGrid* db, int verbose)
    : AMesh(),
      _grid(),
      _nPerCell(0),
      _isPolarized(true),
      _isMaskDefined(false),
      _maskGrid(nullptr)
{
  if (!db->isGrid()) return;
  (void) initFromGrid(db->getNXs(), db->getDXs(), db->getX0s(), db->getRotMat(), true,
                      verbose);
}

MeshETurbo::MeshETurbo(const MeshETurbo &r)
    : AMesh(r),
      _grid(),
      _nPerCell(0),
      _isPolarized(false),
      _isMaskDefined(false),
      _maskGrid(nullptr)
{
  _grid = r._grid;
  _maskGrid = (bool *) mem_free((char *) _maskGrid);
  if (r._isMaskDefined)
    _maskGrid = (bool *) 
      mem_copy((char *) r._maskGrid,sizeof(bool) * _grid.getNTotal(),1);
}

MeshETurbo& MeshETurbo::operator= (const MeshETurbo &r)
{
  _grid = r._grid;
  _maskGrid = (bool *) mem_free((char *) _maskGrid);
  if (r._isMaskDefined)
    _maskGrid = (bool *) 
      mem_copy((char *) r._maskGrid,sizeof(bool) * _grid.getNTotal(),1);
  return *this;
}

MeshETurbo::~MeshETurbo()
{
  _maskGrid = (bool *) mem_free((char *) _maskGrid);
}

int MeshETurbo::getNApices() const
{
  int ndim  = getNDim();
  int napex = 1;
  for (int idim=0; idim<ndim; idim++)
    napex *= _grid.getNX(idim);
  return napex;
}

int MeshETurbo::getNMeshes() const
{
  int ndim  = getNDim();

  int nmesh = 1;
  for (int idim=0; idim<ndim; idim++)
    nmesh *= (_grid.getNX(idim) - 1);
  nmesh *= _nPerCell;
  return nmesh;
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
** \param[in]  imesh    Rank of Mesh (from 0 to _nMeshes-1))
** \param[in]  rank     Rank of Apex within a Mesh (from 0 to _nApexPerMesh-1)
**
*****************************************************************************/
int MeshETurbo::getApex(int imesh, int rank) const
{
  int node,icas;
  int ndim = getNDim();
  VectorInt indg(ndim);

  _fromMeshToIndex(imesh,&node,&icas);
  _grid.rankToIndice(node,indg);
  int ipol = _getPolarized(indg);

  for (int idim = 0; idim < ndim; idim++)
    indg[idim] += MSS(ndim, ipol, icas, rank, idim);

  return _grid.indiceToRank(indg);
}

double MeshETurbo::getCoor(int imesh,
                           int rank,
                           int idim) const
{
  VectorInt indg(getNDim());

  int node = getApex(imesh,rank);
  _grid.rankToIndice(node, indg);
  return _grid.indiceToCoordinate(idim, indg, VectorDouble());
}

double MeshETurbo::getApexCoor(int i, int idim) const
{
  VectorInt indg(getNDim());

  _grid.rankToIndice(i, indg);
  return _grid.indiceToCoordinate(idim, indg, VectorDouble());
}

void MeshETurbo::setMaskArrayFromDouble(double *array)
{
  int ntotal  = _grid.getNTotal();
  int nactive = 0;

  for (int i=0; i<ntotal; i++)
  {
    double value = array[i];
    _maskGrid[i] = (value != 0);
    nactive     += (value != 0);
  }    
  _isMaskDefined = (nactive > 0);
}

void MeshETurbo::setMaskArrayFromInt(int *array)
{
  int ntotal  = _grid.getNTotal();
  int nactive = 0;

  for (int i=0; i<ntotal; i++)
  {
    double value = array[i];
    _maskGrid[i] = (value != 0);
    nactive     += (value != 0);
  }    
  _isMaskDefined = (nactive > 0);
}

int MeshETurbo::initFromGrid(const VectorInt&    nx,
                             const VectorDouble& dx,
                             const VectorDouble& x0,
                             const VectorDouble& rotmat,
                             bool flag_polarized,
                             int verbose)
{
  int ndim = static_cast<int> (nx.size());
  setNDim(ndim);

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
  if (setExtend(extendmin,extendmax)) return(1);

  // Define the number of Elements per Cell

  _setNumberElementPerCell();

  // Set polarization

  _isPolarized = flag_polarized;

  // Optional printout

  if (verbose) messageFlush(toString());

  return 0;
}

int MeshETurbo::dumpToNF(const String& neutralFilename, bool verbose) const
{
  std::ofstream os;
  int ret = 1;
  if (_fileOpenWrite(neutralFilename, "MeshETurbo", os, verbose))
  {
    ret = _serialize(os, verbose);
    if (ret && verbose) messerr("Problem writing in the Neutral File.");
    os.close();
  }
  return ret;
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
  if (_fileOpenRead(neutralFilename, "MeshETurbo", is, verbose))
  {
    mesh = new MeshETurbo;
    if (mesh->_deserialize(is, verbose))
    {
      if (verbose) messerr("Problem reading the Neutral File.");
      delete mesh;
      mesh = nullptr;
    }
    is.close();
  }
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
int MeshETurbo::initFromExtend(const VectorDouble& extendmin,
                                 const VectorDouble& extendmax,
                                 const VectorDouble& cellsize,
                                 const VectorDouble& rotmat,
                                 bool flag_polarized,
                                 int verbose)
{
  int ndim = static_cast<int> (extendmin.size());
  setNDim(ndim);
  if (setExtend(extendmin,extendmax)) return(1);

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


bool MeshETurbo::isNodeMasked(int iabs) const
{
  if (! _isMaskDefined) return(false);
  if (_maskGrid == (bool *) NULL) return(false);
  return _maskGrid[iabs];
}

/****************************************************************************/
/*!
** Returns the Sparse Matrix used to project a Db onto the Meshing
**
** \return A Sparse matrix (cs structure)
**
** \param[in]  db        Db structure
** \param[in]  fatal     Error type when point does not belong to Meshing
** \param[in]  verbose   Verbose flag
**
*****************************************************************************/
cs* MeshETurbo::getMeshToDb(const Db *db, bool fatal, bool verbose) const
{
  int ip_max,iech;
  
  // Initializations

  int error    = 1;
  cs* Atriplet = nullptr;
  cs* A = nullptr;
  double* rhs      = nullptr;
  double* lambda   = nullptr;
  int ndim     = getNDim();
  int ncorner  = getNApexPerMesh();
  VectorInt indg0(ndim);
  VectorInt indgg(ndim);
  VectorInt indices(ncorner);
  VectorDouble coor(ndim);

  // Preliminary checks

  if (isCompatibleDb(db)) goto label_end;

  // Core allocation

  Atriplet = cs_spalloc(0, 0, 1, 1, 1);
  if (Atriplet == nullptr) goto label_end;
  rhs    = (double *) mem_alloc(sizeof(double) * ncorner,0);
  if (rhs    == nullptr) goto label_end;
  lambda = (double *) mem_alloc(sizeof(double) * ncorner,0);
  if (lambda == nullptr) goto label_end;

  /* Optional title */

  if (verbose) mestitle(0,"Mesh Barycenter");

  /* Loop on the samples */

  ip_max = iech = 0;
  for (int jech=0; jech<db->getSampleNumber(); jech++)
  {
    if (! db->isActive(jech)) continue;

    /* Identification */

    for (int idim=0; idim<ndim; idim++)
      coor[idim] = db->getCoordinate(jech,idim);
    
    /* Calculate the grid indices */

    if (_grid.coordinateToIndice(coor,indg0) != 0) 
    {
      messerr("Sample #%d does not belong to the meshing",jech+1);
      if (fatal) goto label_end;
      continue;
    }

    /* Optional printout */

    if (verbose)
      message("Sample %4d in Mesh %4d :",jech+1,
              _grid.indiceToRank(indg0)+1);

    /* Loop on the different meshes constituting the cell */

    int found = -1;
    for (int icas=0; icas<_nPerCell && found<0; icas++)
    {
      if (_addWeights(verbose,icas,indg0,indgg,coor,indices,rhs,lambda) == 0)
      {
        for (int icorner=0; icorner<ncorner; icorner++)
        {
          if (! cs_entry(Atriplet,iech,indices[icorner],lambda[icorner])) 
            goto label_end;
          if (indices[icorner] > ip_max) ip_max = indices[icorner];
        }
        found = icas;
      }
    }
    if (found < 0)
    {
      messerr("Sample #%d does not belong to the meshing",jech+1);
      if (fatal) goto label_end;
      continue;
    }
    iech++;
  }

  /* Add the extreme value to force dimension */

  if (ip_max < getNApices() - 1 || iech < db->getSampleNumber(true) - 1)
  {
    if (!cs_entry(Atriplet, db->getSampleNumber(true) - 1,
                  getNApices() - 1, 0.)) goto label_end;
  }
  
  /* Convert the triplet into a sparse matrix */
  
  A = cs_triplet(Atriplet);

  // Set the error return code

  error = 0;

label_end:
  Atriplet  = cs_spfree(Atriplet);
  rhs       = (double *) mem_free((char *) rhs);
  lambda    = (double *) mem_free((char *) lambda);
  if (error) A = cs_spfree(A);
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
  sstr << toTitle(1,"Turbo Meshing");
  if (_isPolarized) sstr << "Diamond construction is activated" << std::endl;
  _grid.display();
  sstr << AMesh::toString(strfmt);
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

int MeshETurbo::_addWeights(const int verbose,
                            const int icas,
                            const VectorInt& indg0,
                            VectorInt& indgg, // work array
                            const VectorDouble& coor,
                            VectorInt& indices,
                            double *rhs,
                            double *lambda) const
{
  int ndim    =  getNDim();
  int ncorner =  getNApexPerMesh();
  int ipol    = _getPolarized(indg0);
  MatrixSquareGeneral lhs;

  // Build the LHS matrix

  lhs.reset(ncorner,ncorner);
  for (int icorner=0; icorner<ncorner; icorner++)
  {
    // Generate the indices of the mesh apex
    for (int idim=0; idim<ndim; idim++) 
      indgg[idim] = indg0[idim] + MSS(ndim,ipol,icas,icorner,idim);
    indices[icorner] = _grid.indiceToRank(indgg);

    if(indices[icorner]<0) return 1; // outside grid

    // Update the LHS matrix
    for (int idim=0; idim<ndim; idim++)
      lhs.setValue(icorner,idim,_grid.indiceToCoordinate(idim,indgg,VectorDouble()));
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
    if (lambda[icorner] < -EPSILON8) return 1;
    if (lambda[icorner] < 0) lambda[icorner] = 0.;
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

void MeshETurbo::_fromMeshToIndex(int imesh,
                                  int *node,
                                  int *icas) const
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
                             int verbose)
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

  initFromGrid(nx,dx,x0,rot.getMatrixDirectByVector(),true,verbose);
  return 0;
}

int MeshETurbo::_deserialize(std::istream& is, bool /*verbose*/)
{
  int ndim;
  VectorInt nx;
  VectorDouble dx;
  VectorDouble x0;
  VectorDouble rotmat;
  int flag_polarized;

  bool ret = true;
  ret = ret && _recordRead<int>(is, "Space Dimension", ndim);
  ret = ret && _recordReadVec<int>(is, "Nx", nx, ndim);
  ret = ret && _recordReadVec<double>(is, "Dx", dx, ndim);
  ret = ret && _recordReadVec<double>(is, "X0", x0, ndim);
  ret = ret && _recordReadVec<double>(is, "Rotation", rotmat, ndim * ndim);
  ret = ret && _recordRead<int>(is, "Polarization", flag_polarized);

  (void) initFromGrid(nx, dx, x0, rotmat, (bool) flag_polarized);

  return 0;
}

int MeshETurbo::_serialize(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;
  ret = ret && _recordWrite<int>(os, "Space Dimension", getNDim());
  ret = ret && _recordWriteVec<int>(os, "Nx", _grid.getNXs());
  ret = ret && _recordWriteVec<double>(os, "Dx", _grid.getDXs());
  ret = ret && _recordWriteVec<double>(os, "X0", _grid.getX0s());
  ret = ret && _recordWriteVec<double>(os, "Rotation", _grid.getRotMat());
  ret = ret && _recordWrite<int>(os, "Polarization", _isPolarized);
  return 0;
}
