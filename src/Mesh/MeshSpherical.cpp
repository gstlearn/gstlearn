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
#include "Mesh/MeshSpherical.hpp"
#include "Mesh/AMesh.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Db/Db.hpp"
#include "Basic/Vector.hpp"
#include "csparse_f.h"

MeshSpherical::MeshSpherical()
  : AMesh()
  , _apices()
  , _meshes()
  , _units()
{
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
  double unit;

  unit = ut_geodetic_triangle_surface(getCoor(imesh,0,0),
                                      getCoor(imesh,0,1),
                                      getCoor(imesh,1,0),
                                      getCoor(imesh,1,1),
                                      getCoor(imesh,2,0),
                                      getCoor(imesh,2,1));
  return unit;
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
  AMesh::toString(strfmt);
  return sstr.str();
}

/****************************************************************************/
/*!
** Create the meshing
**
** \param[in]  dbin            Pointer to the input Db (optional)
** \param[in]  dbout           Pointer to the output Db (optional)
** \param[in]  triswitch       Construction switch
** \param[in]  verbose         Verbose flag
**
*****************************************************************************/
int MeshSpherical::reset(Db* dbin,
                         Db *dbout,
                         const String& triswitch,
                         int verbose)
{
  int *meshes, ndim, ncorner, napices, nmeshes, error;
  double *apices;
  SphTriangle in;
  Vercoloc *vercoloc;

  /* Initializations */

  error = 1;
  int ndim_ref = 2;
  vercoloc = nullptr;
  meshes = nullptr;
  apices = nullptr;

  /* Define the environment */

  setNDim(ndim_ref);

  /* Initialize the Meshing output structure */
  
  meshes_2D_sph_init(&in);

  /* Set the control points for the triangulation */

  if (dbout != nullptr)
  {
    if (meshes_2D_sph_from_db(dbout,0,NULL,&in)) goto label_end;
  }
  if (dbin != nullptr)
  {
    vercoloc = vercoloc_manage(verbose,1,dbin,dbout,1,NULL);
    if (meshes_2D_sph_from_db(dbin,vercoloc->ndupl,vercoloc->dupl_data,
                              &in)) goto label_end;
    vercoloc = vercoloc_manage(verbose,-1,dbin,dbout,1,vercoloc);
  }

  /* Add auxiliary random points */

  if (meshes_2D_sph_from_auxiliary(triswitch.c_str(),&in)) goto label_end;
  
  /* Perform the triangulation */

  if (meshes_2D_sph_create(verbose,&in)) goto label_end;

  /* Coordinates of the triangle vertices */

  meshes_2D_sph_load_vertices(&in,"Points",
                              &napices,&ndim,(void **) &apices);
  if (ndim != ndim_ref)
  {
    messerr("The space dimension (%d) is not correct. It should be (%d)",
            ndim,ndim_ref);
    goto label_end;
  }
  else
  {
    _apices.reset(napices, ndim);
    _apices.setValues(apices,false);
  }

  meshes_2D_sph_load_vertices(&in,"Triangles",
                              &nmeshes,&ncorner,(void **) &meshes);
  if (ncorner != getNApexPerMesh())
  {
    messerr("Number of Apices per Mesh (%d) is not correct. It should be (%d)",
            ncorner,getNApexPerMesh());
    goto label_end;
  }
  else
  {
    _meshes.resize(nmeshes * ncorner);
    _meshes.assign(meshes, meshes + ncorner * nmeshes);
  }
  
  /* Define and store the Bounding Box extension */

  _defineBoundingBox();

  // Calculate and store the units

  _defineUnits();

  // Optional printout
  
  if (verbose) messageFlush(toString());

  /* Set the error return code */

  error = 0;

label_end:
  vercoloc = vercoloc_manage(verbose,-1,dbin,dbout,1,vercoloc);
  meshes_2D_sph_free(&in,0);
  return error;
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
void MeshSpherical::getDuplicates(Db   *dbin,
                                  Db   *dbout,
                                  int  *nbdupl,
                                  int **dupl1,
                                  int **dupl2,
                                  int   verbose) const
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
** \param[in]  fatal     Error type when point does not belong to Meshing
** \param[in]  verbose   Verbosity flag
**
*****************************************************************************/
cs* MeshSpherical::getMeshToDb(const Db *db, bool fatal, bool /*verbose*/) const
{
  double *coor,*weight;
  int     error,imesh,imesh0,ncorner,ip,found;
  int     iech_max,ip_max,nmeshes,ndim,nvertex;
  cs     *A,*Atriplet;
 
  /* Initializations */
  
  error     = 1;
  coor      = nullptr;
  weight    = nullptr;
  Atriplet  = A = nullptr;
  nmeshes   = getNMeshes();
  nvertex   = getNApices();
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
  
  /* Loop on the samples */

  imesh0 = ip_max = iech_max = 0;
  for (int iech=0; iech<db->getSampleNumber(); iech++)
  {
    if (! db->isActive(iech)) continue;
    if (iech > iech_max) iech_max = iech;

    /* Identification */

    for (int idim=0; idim<ndim; idim++)
      coor[idim] = db->getCoordinate(iech,idim);
    
    /* Loop on the meshes */
    
    found = -1;
    for (int jmesh=0; jmesh<nmeshes; jmesh++)
    {
      imesh = imesh0 + jmesh;
      if (imesh >= nmeshes) imesh -= nmeshes;
      if (! _coorInMesh(coor,imesh,weight)) continue;

      /* Store the items in the sparse matrix */

      for (int icorn=0; icorn<ncorner; icorn++)
      {
        ip = getApex(imesh,icorn);
        if (ip > ip_max) ip_max = ip;
        if (! cs_entry(Atriplet,iech,ip,weight[icorn])) goto label_end;
      }
      imesh0 = found = imesh;
      break;
    }

    /* Printout if a point does not belong to any mesh */

    if (found < 0)
    {
      messerr("Point %d does not belong to any mesh",iech+1);
      for (int idim=0; idim<ndim; idim++)
        messerr(" Coordinate #%d = %lf",idim+1,coor[idim]);
      if (fatal) goto label_end;
    }
  }
  
  /* Add the extreme value to force dimension */

  if (ip_max < nvertex - 1)
  {
    if (! cs_entry(Atriplet,iech_max,nvertex-1,0.)) goto label_end;
  }
  
  /* Convert the triplet into a sparse matrix */

  A = cs_triplet(Atriplet);
  
  /* Set the error return code */

  error = 0;

label_end:
  Atriplet = cs_spfree(Atriplet);
  weight   = (double *) mem_free((char *) weight);
  coor     = db_sample_free(coor);
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
** \remarks The newly allocated array is dimensionned to the total number 
** \remarks of samples in the Db (masked included).
** \remarks It must be freed by the calling function
*,nech;
*****************************************************************************/
double* MeshSpherical::interpolateMeshToDb(Db *db,
                                           double* mtab) const
{
  double *coor,*weight,*dtab,total;
  int     error,imesh,imesh0,ncorner,ip,found;
  int     iech_max,ip_max,nmeshes,ndim,nech;
 
  /* Initializations */
  
  error     = 1;
  dtab      = nullptr;
  coor      = nullptr;
  weight    = nullptr;
  nmeshes   = getNMeshes();
  ncorner   = getNApexPerMesh();
  ndim      = getNDim();
  nech      = db->getSampleNumber();

  // Preliminary checks 

  if (isCompatibleDb(db)) goto label_end;

  /* Core allocation */

  dtab   = (double *) mem_alloc(sizeof(double) * nech,0);
  if (dtab == nullptr) goto label_end;
  coor   = db_sample_alloc(db,ELoc::X);
  if (coor   == nullptr) goto label_end;
  weight = (double *) mem_alloc(sizeof(double) * ncorner,0);
  if (weight == nullptr) goto label_end;
  
  /* Loop on the samples */

  imesh0 = ip_max = iech_max = 0;
  for (int iech=0; iech<db->getSampleNumber(); iech++)
  {
    if (! db->isActive(iech)) continue;
    if (iech > iech_max) iech_max = iech;

    /* Identification */

    for (int idim=0; idim<ndim; idim++)
      coor[idim] = db->getCoordinate(iech,idim);
    
    /* Loop on the meshes */
    
    found = -1;
    for (int jmesh=0; jmesh<nmeshes; jmesh++)
    {
      imesh = imesh0 + jmesh;
      if (imesh >= nmeshes) imesh -= nmeshes;
      if (! _coorInMesh(coor,imesh,weight)) continue;

      /* Store the items in the sparse matrix */

      total = 0.;
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
      messerr("Point %d does not belong to any mesh",iech+1);
      for (int idim=0; idim<ndim; idim++)
        messerr(" Coordinate #%d = %lf",idim+1,coor[idim]);
    }
  }
  
  /* Set the error return code */

  error = 0;

label_end:
  weight   = (double *) mem_free((char *) weight);
  coor     = db_sample_free(coor);
  if (error) dtab = (double *) mem_free((char *) dtab);
  return(dtab);
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
int MeshSpherical::getApex(int imesh,
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
double MeshSpherical::getCoor(int imesh, int rank, int idim) const
{
  return _apices(getApex(imesh,rank),idim);
}

double MeshSpherical::getApexCoor(int i, int idim) const
{
  return _apices(i,idim);
}

void MeshSpherical::getEmbeddedCoor(int imesh, int ic, VectorDouble& coords) const
{
  /* The Variety is defined in the Global Environment */
  /* The required radius is set to the radius of Earth (6371m) */

  int variety_sphere;
  double r;
  variety_query(&variety_sphere);
  if(variety_sphere==1)
  {
    variety_get_characteristics(&r);
  }
  else
  {
    r = 6371.;
  }

  util_convert_sph2cart(getCoor(imesh, ic, 0)-180.,
                        getCoor(imesh, ic, 1),
                        &coords[0], &coords[1], &coords[2],
                        r);

}

void MeshSpherical::_defineUnits(void)
{
  int nmeshes = getNMeshes();

  // Core allocation

  _units.resize(nmeshes);

  // Loop on the meshes 

  for (int imesh=0; imesh<nmeshes; imesh++)
    _units[imesh ] = getMeshSize(imesh);
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
  (void) setExtend(extendmin,extendmax);
}

/****************************************************************************/
/*!
**  Check if a point, defined by its coordinates, belong to a Mesh
**
** \return true if the point belongs to the Mesh; false otherwise
**
** \param[in]  coor      Array of target coordinates
** \param[in]  imesh     Mesh Index
**
** \param[out] weights   Array of barycentric weights (Dim: NApexPerMesh)
** 
** \remarks The argument 'meshsize' is used to speed the calculations
**
*****************************************************************************/
bool MeshSpherical::_coorInMesh(double    *coor,
                                int        imesh,
                                double    *weights) const
{
  double pts[3][2];

  for (int i=0; i<3; i++)
    for (int j=0; j<2; j++)
      pts[i][j] = getCoor(imesh,i,j);
  return (is_in_spherical_triangle_optimized(coor, pts[0], pts[1], pts[2],
                                             weights));
}

int MeshSpherical::_recopy(const MeshSpherical &m)
{
  _apices = m._apices;
  _meshes = m._meshes;
  _units = m._units;
  return(0);
}
