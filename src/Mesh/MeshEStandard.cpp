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
#include "Mesh/MeshEStandard.hpp"
#include "Basic/AException.hpp"
#include "Db/Db.hpp"
#include "Mesh/tetgen.h"
#include "csparse_f.h"


MeshEStandard::MeshEStandard()
  : AMesh()
  , _apices()
  , _meshes()
  , _units()
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
int MeshEStandard::getNApices() const 
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
int MeshEStandard::getNMeshes() const
{
  return (static_cast<int> (_meshes.size()) / getNApexPerMesh());
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
double MeshEStandard::getMeshSize(int imesh) const
{
  double unit;
  static double facdim[] = {0., 1., 2., 6.};

  // Initializations
  int ndim    = getNDim();
  int ncorner = getNApexPerMesh();

  // Calculate the mesh size
  MatrixSquareGeneral mat;
  mat.reset(ndim,ndim);
  for (int icorn=1; icorn<ncorner; icorn++)
    for (int idim=0; idim<ndim; idim++)
      mat.setValue(icorn-1,idim,
                   getCoor(imesh,icorn,idim) - getCoor(imesh,0,idim));
  unit = ABS(mat.determinant()) / facdim[ndim];

  return unit;
}

/****************************************************************************/
/*!
** Create the meshing
**
** \param[in]  dbin            Pointer to the input Db (optional)
** \param[in]  dbout           Pointer to the output Db (optional)
** \param[in]  dilate          Dilation of the Bounding box (optional)
** \param[in]  triswitch       Construction switch
** \param[in]  verbose         Verbose flag
**
*****************************************************************************/
int MeshEStandard::resetFromDb(Db*                 dbin,
                          Db*                 dbout,
                          const VectorDouble& dilate,
                          const String&       triswitch,
                          bool                verbose)
{
  int error = 1;

  // Discover the space dimension
  int ndim = 0;
  if (dbin != NULL)
    ndim = dbin->getNDim();
  if (dbout != NULL)
  {
    if (ndim > 0 && ndim != dbout->getNDim())
      my_throw("'dbin' and 'dbout' are both defined but with different 'ndim'");
    ndim = dbout->getNDim();
  }
  setNDim(ndim);

  // Preliminary checks 

  if (! dilate.empty())
  {
    for (int idim=0; idim<ndim; idim++)
      if (dilate[idim] < 0)
      {
        messerr("The dilation[%d] is negative",idim+1);
        return(1);
      }
  }

  // Dispatch according to SPace Dimension

  if (ndim == 1)
    error = _create1D(1,verbose,dbin,dbout,dilate);
  else if (ndim == 2)
    error = _create2D(2,verbose,dbin,dbout,dilate,triswitch.c_str());
  else if (ndim == 3)
    error = _create3D(3,verbose,dbin,dbout,dilate,triswitch.c_str());
  else
  {
    messerr("Meshing is only provided for Space Dimension 1, 2 or 3");
  }
  if (error) return(error);
  
  // Define and store the Bounding Box extension

  _defineBoundingBox();

  // Calculate and store the units

  _defineUnits();

  // Optional printout
  
  if (verbose) messageFlush(toString());

  return(error);
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
int MeshEStandard::reset(const MatrixRectangular& apices,
                          const VectorInt&          meshes,
                          bool                      verbose)
{
  int ndim = apices.getNCols();
  setNDim(ndim);

  // Core allocation

  _apices = apices;
  _meshes = meshes;

  // Check consistency

  _checkConsistency();

  // Define and store the Bounding Box extension

  _defineBoundingBox();

  // Calculate and store the units

  _defineUnits();

  // Optional printout
  
  if (verbose) messageFlush(toString());

  return(0);
}

/****************************************************************************/
/*!
** Create the meshing (from mesh information)
**
** \param[in]  ndim            Space Dimension
** \param[in]  apices          Vector of Apex information
** \param[in]  meshes          Vector of mesh indices
** \param[in]  verbose         Verbose flag
**
*****************************************************************************/
int MeshEStandard::resetOldStyle(int                 ndim,
                          const VectorDouble& apices,
                          const VectorInt&    meshes,
                          bool                verbose)
{
  setNDim(ndim);
  int npoints = static_cast<int> (apices.size()) / ndim;

  // Core allocation

  _apices.reset(npoints,ndim);
  _apices.setValues(apices);
  _meshes = meshes;

  // Check consistency

  _checkConsistency();

  // Define and store the Bounding Box extension

  _defineBoundingBox();

  // Calculate and store the units

  _defineUnits();

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
int MeshEStandard::getApex(int imesh,
                           int rank) const
{
  return (_meshes[getNApexPerMesh() * imesh + rank]);
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
double MeshEStandard::getCoor(int imesh,
                              int rank,
                              int idim) const
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
  sstr << toTitle(0,"Standard Meshing characteristics");

  sstr << toTitle(1,"Standard Meshing");
  sstr << AMesh::toString(strfmt);
  return sstr.str();
}

/****************************************************************************/
/*!
** Returns the duplicate information (if any) 
**
** \param[in]  verbose         Verbose flag
** \param[in]  dbin            Pointer to the input Db (optional)
** \param[in]  dbout           Pointer to the output Db (optional)
**
** \param[out]  nbdupl  Number of duplicates
** \param[out]  dupl1   Array giving ranks of duplicate samples in Db1
** \param[out]  dupl2   Array giving ranks of duplicate samples in Db2
**
*****************************************************************************/
void MeshEStandard::getDuplicates(int   verbose,
                                  Db   *dbin,
                                  Db   *dbout,
                                  int  *nbdupl,
                                  int **dupl1,
                                  int **dupl2) const
{
  Vercoloc *vercoloc;

  // Initializations

  *nbdupl = 0;
  *dupl1 = nullptr;
  *dupl2 = nullptr;
  if (dbin == nullptr) return;
  
  // Look for duplicates

  vercoloc = vercoloc_manage(verbose,1,dbin,dbout,1,NULL);
  *nbdupl = vercoloc->ndupl;
  *dupl1  = (int *) mem_copy((char *) vercoloc->dupl_dabs,
                             sizeof(int) * vercoloc->ndupl,1);
  *dupl2  = (int *) mem_copy((char *) vercoloc->dupl_grid,
                             sizeof(int) * vercoloc->ndupl,1);
  vercoloc = vercoloc_manage(verbose,-1,dbin,dbout,1,vercoloc);
}

/****************************************************************************/
/*!
** Returns the Sparse Matrix used to project a Db onto the Meshing
**
** \return A Sparse matrix (cs structure)
**
** \param[in]  db        Db structure
** \param[in]  verbose   Verbosity flag
**
*****************************************************************************/
cs* MeshEStandard::getMeshToDb(const Db *db,
                               int verbose) const
{
  double *coor,*container,*weight;
  int     error,imesh,imesh0,ncorner,ip,found,iech,ip_max,ndim,nmeshes;
  cs     *A,*Atriplet;
 
  /* Initializations */

  error     = 1;
  coor      = nullptr;
  container = nullptr;
  weight    = nullptr;
  Atriplet  = A = nullptr;
  nmeshes   = getNMeshes();
  ncorner   = getNApexPerMesh();
  ndim      = getNDim();

  // Preliminary checks 

  if (isCompatibleDb(db)) goto label_end;

  /* Core allocation */

  Atriplet = cs_spalloc(0, 0, 1, 1, 1);
  if (Atriplet == nullptr) goto label_end;

  coor   = db_sample_alloc(db,ELoc::X);
  if (coor   == nullptr) goto label_end;
  weight = (double *) mem_alloc(sizeof(double) * ncorner,0);
  if (weight == nullptr) goto label_end;
  container = _defineContainers();
  if (container == nullptr) goto label_end;
  
  /* Optional title */

  if (verbose) mestitle(0,"Mesh Barycenter");

  /* Loop on the samples */

  imesh0 = ip_max = iech = 0;
  for (int jech=0; jech<db->getSampleNumber(); jech++)
  {
    if (! db->isActive(jech)) continue;

    /* Identification */

    for (int idim=0; idim<ndim; idim++)
      coor[idim] = db->getCoordinate(jech,idim);
    
    /* Loop on the meshes */
    
    found = -1;
    for (int jmesh=0; jmesh<nmeshes; jmesh++)
    {
      imesh = imesh0 + jmesh;
      if (imesh >= nmeshes) imesh -= nmeshes;
      if (! _coorInMeshContainer(coor,imesh,container)) continue;
      if (! _coorInMesh(coor,imesh,_units[imesh],weight)) continue;

      /* Store the items in the sparse matrix */

      if (verbose) message("Sample %4d in Mesh %4d :",jech+1,imesh+1);
      for (int icorn=0; icorn<ncorner; icorn++)
      {
        ip = getApex(imesh,icorn);
        if (ip > ip_max) ip_max = ip;
        if (verbose) message(" %4d (%4.2lf)",ip,weight[icorn]);
        if (! cs_entry(Atriplet,iech,ip,weight[icorn])) goto label_end;
      }
      if (verbose) message("\n");
      imesh0 = found = imesh;
      break;
    }

    /* Printout if a point does not belong to any mesh */

    if (found < 0)
    {
      messerr("Point %d does not belong to any mesh",jech+1);
      for (int idim=0; idim<ndim; idim++)
        messerr(" Coordinate #%d = %lf",idim+1,coor[idim]);
    }
    iech++;
  }
  
  /* Add the extreme value to force dimension */

  if (ip_max < getNApices() - 1)
  {
    if (!cs_entry(Atriplet, db->getSampleNumber(true) - 1,
                  getNApices() - 1, 0.)) goto label_end;
  }
  
  /* Convert the triplet into a sparse matrix */

  A = cs_triplet(Atriplet);
  
  /* Set the error return code */

  error = 0;

label_end:
  Atriplet  = cs_spfree(Atriplet);
  container = (double *) mem_free((char *) container);
  weight    = (double *) mem_free((char *) weight);
  coor      = db_sample_free(coor);
  if (error) A = cs_spfree(A);
  return(A);
}

/****************************************************************************/
/*!
** Interpolates an array defined on the Meshes to the Db locations
**
** \return The array of interpolated values
**
** \param[in]  db        Db structure
** \param[in]  mtab      Array of values defined on the meshes
**
** \remarks The newly allocated array is dimensioned to the total number
** \remarks of samples in the Db (masked included).
** \remarks It must be freed by the calling function
**
*****************************************************************************/
double* MeshEStandard::interpolateMeshToDb(Db *db,
                                           double* mtab) const
{
  double *coor,*container,*weight,*dtab;
  int     error,imesh,imesh0,ncorner,ip,found,nech,ndim,nmeshes,iech;
 
  /* Initializations */
  
  error     = 1;
  coor      = nullptr;
  container = nullptr;
  weight    = nullptr;
  dtab      = nullptr;
  nmeshes   = getNMeshes();
  ncorner   = getNApexPerMesh();
  ndim      = getNDim();
  nech      = db->getSampleNumber(true);

  // Preliminary checks

  if (isCompatibleDb(db)) goto label_end;

  /* Core allocation */

  dtab   = (double *) mem_alloc(sizeof(double) * nech,0);
  if (dtab == nullptr) goto label_end;
  coor   = db_sample_alloc(db,ELoc::X);
  if (coor   == nullptr) goto label_end;
  weight = (double *) mem_alloc(sizeof(double) * ncorner,0);
  if (weight == nullptr) goto label_end;
  container = _defineContainers();
  if (container == nullptr) goto label_end;
  
  /* Loop on the samples */

  imesh0 = iech = 0;
  for (int jech=0; jech<db->getSampleNumber(); jech++)
  {
    if (! db->isActive(jech)) continue;

    /* Identification */

    for (int idim=0; idim<ndim; idim++)
      coor[idim] = db->getCoordinate(jech,idim);
    
    /* Loop on the meshes */
    
    found = -1;
    for (int jmesh=0; jmesh<nmeshes; jmesh++)
    {
      imesh = imesh0 + jmesh;
      if (imesh >= nmeshes) imesh -= nmeshes;
      if (! _coorInMeshContainer(coor,imesh,container)) continue;
      if (! _coorInMesh(coor,imesh,_units[imesh],weight)) continue;

      // Perform the interpolation

      double total = 0.;
      for (int icorn=0; icorn<ncorner; icorn++)
      {
        ip = getApex(imesh,icorn);
        total += mtab[ip] * weight[icorn];
      }
      dtab[iech] = total;
      imesh0 = found = imesh;
      break;
    }

    /* Printout if a point does not belong to any mesh */

    if (found < 0)
    {
      messerr("Point %d does not belong to any mesh",jech+1);
      for (int idim=0; idim<ndim; idim++)
        messerr(" Coordinate #%d = %lf",idim+1,coor[idim]);
    }
    iech++;
  }
  
  /* Set the error return code */

  error = 0;

label_end:
  container = (double *) mem_free((char *) container);
  weight    = (double *) mem_free((char *) weight);
  coor      = db_sample_free(coor);
  if (error) dtab = (double *) mem_free((char *) dtab);
  return(dtab);
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
    for (int icol = 0; icol < getNDim(); icol++)
      for (int irow = 0; irow < getNApices(); irow++)
      {
        list.push_back(getApexCoor(irow, icol));
      }
  }
  else
  {
    for (int irow = 0; irow < getNApices(); irow++)
      for (int icol = 0; icol < getNDim(); icol++)
      {
        list.push_back(getApexCoor(irow, icol));
      }
  }
  return list;
}

/****************************************************************************/
/*!
** Create the meshing in 1D
**
** \param[in]  ndim_ref        Space dimension
** \param[in]  verbose         Verbose flag
** \param[in]  dbin            Pointer to the input Db (optional)
** \param[in]  dbout           Pointer to the output Db (optional)
** \param[in]  dilate          Dilation of the Bounding box (optional)
**
*****************************************************************************/
int MeshEStandard::_create1D(int                 ndim_ref,
                             int                 verbose,
                             Db*                 dbin,
                             Db*                 dbout,
                             const VectorDouble& dilate)
{
  int       error,ndim,ncorner,napices,nmeshes,flag_defined;
  segmentio in,out;
  Vercoloc *vercoloc;
  double   *apices;
  int      *meshes;

  /* Initializations */

  error    = 1;
  vercoloc = nullptr;
  apices = nullptr;
  meshes = nullptr;

  /* Initialize the Meshing structure */
  
  meshes_1D_init(&in);
  meshes_1D_init(&out);
  
  /* Set the control points for the triangulation */
  
  flag_defined = 0;
  if (dbout != nullptr)
  {
    if (meshes_1D_from_db(dbout,0,NULL,&in)) goto label_end;
    flag_defined = 1;
  }
  if (dbin != nullptr)
  {
    vercoloc = vercoloc_manage(verbose,1,dbin,dbout,1,NULL);
    if (meshes_1D_from_db(dbin,vercoloc->ndupl,vercoloc->dupl_dabs,
                          &in)) goto label_end;
    vercoloc = vercoloc_manage(verbose,-1,dbin,dbout,1,vercoloc);
    flag_defined = 1;
  }
  if (! flag_defined)
  {
    meshes_1D_default(dbin,dbout,&in);
    flag_defined = 1;
  }
  
  /* Extend the domain if 'dilate' is specified */
  
  if (dbout != nullptr && ! dilate.empty())
    meshes_1D_extended_domain(dbout,dilate.data(),&in);
  
  /* Perform the triangulation */
  
  meshes_1D_create(verbose,&in,&out);
  
  /* Coordinates of the triangle vertices */
  
  meshes_1D_load_vertices(&out,"Points",&napices,&ndim,(void **) apices);
  if (ndim != ndim_ref)
  {
    messerr("The space dimension (%d) is not correct. It should be (%d)",
            ndim,ndim_ref);
    goto label_end;
  }
  _apices.reset(napices, ndim);
  _apices.setValues(apices, false);
  
  meshes_1D_load_vertices(&out,"Segments",&nmeshes,&ncorner,(void **) &meshes);
  if (ncorner != getNApexPerMesh())
  {
    messerr("The number of Apices per Mesh (%d) is not correct (%d)",
            ncorner,getNApexPerMesh());
    goto label_end;
  }
  _meshes.resize(nmeshes * ncorner);
  _meshes.assign(meshes, meshes + ncorner * nmeshes);
  ut_ivector_addval(_meshes,-1);

  /* Set the error return code */

  error = 0;

label_end:
  vercoloc = vercoloc_manage(verbose,-1,dbin,dbout,1,vercoloc);
  apices = (double *) mem_free((char *) apices);
  meshes = (int    *) mem_free((char *) meshes);
  meshes_1D_free(&in,1);
  meshes_1D_free(&out,0);
  return(error);
}

/****************************************************************************/
/*!
** Create the meshing in 2D
**
** \param[in]  ndim_ref        Space dimension
** \param[in]  verbose         Verbose flag
** \param[in]  dbin            Pointer to the input Db (optional)
** \param[in]  dbout           Pointer to the output Db (optional)
** \param[in]  dilate          Dilation of the Bounding box (optional)
** \param[in]  triswitch       Construction switch
**
*****************************************************************************/
int MeshEStandard::_create2D(int                 ndim_ref,
                             int                 verbose,
                             Db*                 dbin,
                             Db*                 dbout,
                             const VectorDouble& dilate,
                             const char*         triswitch)
{
  int           error,napices,ndim,ncorner,nmeshes,flag_defined;
  triangulateio in,out,vorout;
  Vercoloc     *vercoloc;
  double *apices;
  int    *meshes;

  /* Initializations */

  error    = 1;
  vercoloc = nullptr;
  apices   = nullptr;
  meshes   = nullptr;

  /* Initialize the Meshing structure */
  
  meshes_2D_init(&in);
  meshes_2D_init(&out);
  meshes_2D_init(&vorout);
  
  /* Set the control points for the triangulation */
  
  flag_defined = 0;
  if (dbout != nullptr)
  {
    if (meshes_2D_from_db(dbout,1,0,NULL,&in)) goto label_end;
    flag_defined = 1;
  }
  if (dbin != nullptr)
  {
    vercoloc = vercoloc_manage(verbose,1,dbin,dbout,1,NULL);
    if (meshes_2D_from_db(dbin,1,vercoloc->ndupl,vercoloc->dupl_dabs,
                          &in)) goto label_end;
    vercoloc = vercoloc_manage(verbose,-1,dbin,dbout,1,vercoloc);
    flag_defined = 1;
  }
  if (! flag_defined)
  {
    meshes_2D_default(dbin,dbout,&in);
    flag_defined = 1;
  }
  
  /* Extend the domain if 'dilate' is specified */
  
  if (dbout != nullptr && ! dilate.empty())
    meshes_2D_extended_domain(dbout,dilate.data(),&in);
  
  /* Perform the triangulation */
  
  meshes_2D_create(verbose,triswitch,&in,&out,&vorout);
  
  /* Coordinates of the triangle vertices */
  
  meshes_2D_load_vertices(&out,"Points",&napices,&ndim,(void **) &apices);
  if (ndim != ndim_ref)
  {
    messerr("Space dimension (%d) is not correct (%d)",
            ndim,ndim_ref);
    goto label_end;
  }
  _apices.reset(napices, ndim);
  _apices.setValues(apices, false);
  
  meshes_2D_load_vertices(&out,"Triangles",&nmeshes,&ncorner,(void **) &meshes);
  if (ncorner != getNApexPerMesh())
  {
    messerr("Number of Apices per Mesh (%d) is not correct. It should be (%d)",
            ncorner,getNApexPerMesh());
    goto label_end;
  }
  _meshes.resize(nmeshes * ncorner);
  _meshes.assign(meshes, meshes + ncorner * nmeshes);
  ut_ivector_addval(_meshes,-1);

  /* Set the error return code */

  error = 0;

label_end:
  vercoloc = vercoloc_manage(verbose,-1,dbin,dbout,1,vercoloc);
  apices = (double *) mem_free((char *) apices);
  meshes = (int    *) mem_free((char *) meshes);
  meshes_2D_free(&in,1);
  meshes_2D_free(&out,0);
  meshes_2D_free(&vorout,0);
  return(error);
}

/****************************************************************************/
/*!
** Create the meshing in 3D
**
** \param[in]  ndim_ref        Space dimension
** \param[in]  verbose         Verbose flag
** \param[in]  dbin            Pointer to the input Db (optional)
** \param[in]  dbout           Pointer to the output Db (optional)
** \param[in]  dilate          Dilation of the Bounding box (optional)
** \param[in]  triswitch       Construction switch
**
*****************************************************************************/
int MeshEStandard::_create3D(int                 ndim_ref,
                             int                 verbose,
                             Db*                 dbin,
                             Db*                 dbout,
                             const VectorDouble& dilate,
                             const char*         triswitch)
{
  int      *meshes,error,ndim,ncorner,napices,nmeshes,flag_defined;
  double   *apices;
  tetgenio  in,out;
  Vercoloc *vercoloc;

  /* Initializations */

  error    = 1;
  vercoloc = nullptr;
  apices   = nullptr;
  meshes   = nullptr;

  /* Set the control points for the tetrahedralization */

  flag_defined = 0;
  if (dbout != nullptr)
  {
    if (meshes_3D_from_db(dbout,0,NULL,&in)) goto label_end;
    flag_defined = 1;
  }
  if (dbin != nullptr)
  {
    vercoloc = vercoloc_manage(verbose,1,dbin,dbout,1,NULL);
    if (meshes_3D_from_db(dbin,vercoloc->ndupl,vercoloc->dupl_dabs,
                          &in)) goto label_end;
    vercoloc = vercoloc_manage(verbose,-1,dbin,dbout,1,vercoloc);
    flag_defined = 1;
  }
  if (! flag_defined)
    meshes_3D_default(dbin,dbout,&in);

  /* Extend the domain if 'dilate' is specified */

  if (dbout != nullptr && ! dilate.empty())
    meshes_3D_extended_domain(dbout,dilate.data(),&in);

  /* Perform the tetrahedralization */

  meshes_3D_create(verbose,triswitch,&in,&out);

  /* Coordinates of the vertices */

  meshes_3D_load_vertices(&out,"Points",&napices,&ndim,(void **) &apices);
  if (ndim != ndim_ref)
  {
    messerr("The space dimension (%d) is not correct (%d)",
            ndim,ndim_ref);
    goto label_end;
  }
  _apices.reset(napices, ndim);
  _apices.setValues(apices, false);
  
  meshes_3D_load_vertices(&out,"Tetrahedra",&nmeshes,&ncorner,(void **) &meshes);
  if (ncorner != getNApexPerMesh())
  {
    messerr("The number of Apices per Mesh (%d) is not correct (%d)",
            ncorner,getNApexPerMesh());
    goto label_end;
  }
  _meshes.resize(nmeshes * ncorner);
  _meshes.assign(meshes, meshes + ncorner * nmeshes);
  ut_ivector_addval(_meshes,-1);

  /* Set the error return code */

  error = 0;

label_end:
  vercoloc = vercoloc_manage(verbose,-1,dbin,dbout,1,vercoloc);
  apices   = (double *) mem_free((char *) apices);
  meshes   = (int    *) mem_free((char *) meshes);
  meshes_3D_free(&in);
  meshes_3D_free(&out);
  return(error);
}

double MeshEStandard::getApexCoor(int i, int idim) const
{
  return _apices(i, idim);
}

void MeshEStandard::_defineUnits(void)
{
  int nmeshes = getNMeshes();

  // Core allocation

  _units.resize(nmeshes);

  // Loop on the meshes 

  for (int imesh=0; imesh<nmeshes; imesh++)
    _units[imesh ] = getMeshSize(imesh);
}

void MeshEStandard::_defineBoundingBox(void)
{
  VectorDouble extendmin, extendmax;
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

/****************************************************************************/
/*!
**  Define the container for each mesh
**
** \return Pointer to the array containing the containers
**
*****************************************************************************/
double* MeshEStandard::_defineContainers() const

{
  double *container,value,vmin,vmax;
  int ncorner,ndim,nmesh;

  /* Initializations */

  ndim    = getNDim();
  nmesh   = getNMeshes();
  ncorner = getNApexPerMesh();
  
  /* Allocation */

  container = (double *) mem_alloc(sizeof(double) * 2 * ndim * nmesh,0);
  if (container == nullptr) return(container);
  
  /* Loop on the meshes */

  for (int imesh=0; imesh<nmesh; imesh++)
  {

    /* Loop on the dimensions */
    
    for (int idim=0; idim<ndim; idim++)
    {
      vmin =  1.e30;
      vmax = -1.e30;

      /* Loop on the corners */
      
      for (int icorn=0; icorn<ncorner; icorn++)
      {
        value = getCoor(imesh,icorn,idim);
        if (value < vmin) vmin = value;
        if (value > vmax) vmax = value;
      }
      _setContainer(container,imesh,idim,vmin,vmax);
    }
  }
  return(container);  
}

/****************************************************************************/
/*!
**  Check if a point, defined by its coordinates, belong to the container
**  of a Mesh
**
** \return true if the point belongs to the container; false otherwise
**
** \param[in]  coor      Array of target coordinates
** \param[in]  imesh     Mesh Index
** \param[in]  container Array of containers
**
*****************************************************************************/
bool MeshEStandard::_coorInMeshContainer(double *coor,
                                         int     imesh,
                                         double *container) const
{
  double center,vmin,vmax;
  
  /* Initializations */
  
  if (container == nullptr) return true;
  
  /* Loop on the space dimension */

  for (int idim=0; idim<getNDim(); idim++)
  {
    center = coor[idim];
    _getContainer(container,imesh,idim,&vmin,&vmax);
    if ((center - vmin) * (center - vmax) > 0) return false;
  }
  return true;  
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
bool MeshEStandard::_coorInMesh(double    *coor,
                                int        imesh,
                                double     meshsize,
                                double    *weights) const
{
  double ratio,total;
  int    ncorner,ndim,kcorn;
  static double FACDIM[] = {0., 1., 2., 6.};

  // Initialiations
  ncorner = getNApexPerMesh();
  ndim    = getNDim();
  
  // Loop on the vertices
  total = 0.;
  for (int icorn=0; icorn<ncorner; icorn++)
  {
    
    // Build the determinant
    MatrixSquareGeneral mat;
    mat.reset(ndim,ndim);
    kcorn = 0;
    for (int jcorn=0; jcorn<ncorner; jcorn++)
    {
      if (icorn == jcorn) continue;
      for (int idim=0; idim<ndim; idim++) 
        mat.setValue(idim,kcorn,getCoor(imesh,jcorn,idim) - coor[idim]);
      //        mat[ecr++] = getCoor(imesh,jcorn,idim) - coor[idim];
      kcorn++;
    }
    ratio = ABS(mat.determinant()) / meshsize / FACDIM[ndim];
    if (ratio < 0 || ratio > 1) return false;
    weights[icorn] = ratio;
    total += ratio;
  }
  if (ABS(total - 1) > 1.e-3) return false;
  return true;
}

void MeshEStandard::_setContainer(double *container,
                                  int     imesh,
                                  int     idim,
                                  double  vmin,
                                  double  vmax) const
{
  int ndim = getNDim();
  container[(0) + 2*((idim) + ndim * (imesh))] = vmin;
  container[(1) + 2*((idim) + ndim * (imesh))] = vmax;
}

void MeshEStandard::_getContainer(double *container,
                                  int     imesh,
                                  int     idim,
                                  double *vmin,
                                  double *vmax) const
{
  int ndim = getNDim();
  *vmin = container[(0) + 2*((idim) + ndim * (imesh))];
  *vmax = container[(1) + 2*((idim) + ndim * (imesh))];
}

void MeshEStandard::_printContainers(double *container) const
{
  int    ndim,nmesh,ncorner;
  double vmin,vmax;

  /* Initializations */

  ndim    = getNDim();
  nmesh   = getNMeshes();
  ncorner = getNApexPerMesh();
  
  /* Loop on the meshes */

  for (int imesh=0; imesh<nmesh; imesh++)
  {
    
    /* Printing the mesh */
    
    message("Mesh #%d/%d\n",imesh+1,nmesh);
    for (int icorn=0; icorn<ncorner; icorn++)
    {
      for (int idim=0; idim<ndim; idim++)
        message(" %lf",getCoor(imesh,icorn,idim));
      message("\n");
    }
          
    /* Printing the container */
    
    message(" Container\n");
    for (int idim=0; idim<ndim; idim++)
    {
      _getContainer(container,imesh,idim,&vmin,&vmax);
      message(" %lf - %lf\n",vmin,vmax);
    }
  }
}

void MeshEStandard::_deallocate()
{

}

int MeshEStandard::_recopy(const MeshEStandard &m)
{
  _deallocate();

  _apices = m._apices;
  _meshes = m._meshes;
  _units = m._units;

  return(0);
}

/**
 * This function checks the consistency between the number of points
 * and the vertices indication
 */
void MeshEStandard::_checkConsistency() const
{
  for (int imesh = 0; imesh < getNMeshes(); imesh++)
    for (int ic = 0; ic < getNApexPerMesh(); ic++)
    {
      int apex = getApex(imesh, ic);
      if (apex < 0 || apex >= getNApices())
        throw("Mesh indices are not compatible with the Points");
    }
}

/**
 * Convert the Meshing (Old class) into a meshing (New class)
 * @param s_mesh Meshing characteristics (Old class)
 * @param verbose Verbosity flag
 */
int MeshEStandard::convertFromOldMesh(SPDE_Mesh* s_mesh, int verbose)
{
  int lec;
  int napices = s_mesh->nvertex;
  int ndim    = s_mesh->ndim;
  MatrixRectangular points(napices, ndim);
  lec = 0;
  for (int ip = 0; ip < napices; ip++)
    for (int idim = 0 ; idim < ndim; idim++, lec++)
    {
      points.setValue(ip, idim, s_mesh->points[lec]);
    }

  // For safety, the minimum value of the array meshes is considered
  // If equal to , the mesh numbering must be modified in order to start from 0
  // This code is temporary and there to cope with different numberings

  int nmeshes = s_mesh->nmesh;
  int ncorner = s_mesh->ncorner;
  int shift_min = 1000;
  for (int i = 0; i < nmeshes * ncorner; i++)
    if (s_mesh->meshes[i] < shift_min) shift_min = s_mesh->meshes[i];
  if (shift_min != 0 && shift_min != 1)
    throw("Wrong minimum shift: it should be 0 or 1");

  VectorInt meshes = VectorInt(nmeshes * ncorner);
  lec = 0;
  for (int imesh = 0; imesh < nmeshes; imesh++)
    for (int ic = 0; ic < ncorner; ic++, lec++)
    {
      meshes[lec] = s_mesh->meshes[lec] - shift_min;
    }

  int error = reset(points, meshes, verbose);
  return error;
}
