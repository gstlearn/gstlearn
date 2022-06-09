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
#include "Matrix/MatrixInt.hpp"
#include "Db/Db.hpp"
#include "Basic/Vector.hpp"
#include "Basic/Geometry.hpp"
#include "Space/SpaceSN.hpp"
#include "csparse_f.h"

MeshSpherical::MeshSpherical()
  : AMesh()
  , _apices()
  , _meshes()
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
  return ut_geodetic_triangle_surface(getCoor(imesh, 0, 0),
                                      getCoor(imesh, 0, 1),
                                      getCoor(imesh, 1, 0),
                                      getCoor(imesh, 1, 1),
                                      getCoor(imesh, 2, 0),
                                      getCoor(imesh, 2, 1));
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
    _meshes.reset(nmeshes,ncorner);
    _meshes.setValues(meshes, false);
  }
  
  /* Define and store the Bounding Box extension */

  _defineBoundingBox();

  // Optional printout
  
  if (verbose) messageFlush(toString());

  /* Set the error return code */

  error = 0;

label_end:
  vercoloc = vercoloc_manage(verbose,-1,dbin,dbout,1,vercoloc);
  meshes_2D_sph_free(&in,0);
  return error;
}

int MeshSpherical::dumpToNF(const String& neutralFilename, bool verbose) const
{
  std::ofstream os;
  int ret = 1;
  if (_fileOpenWrite(neutralFilename, "MeshSpherical", os, verbose))
  {
    ret = _serialize(os, verbose);
    if (ret && verbose) messerr("Problem writing in the Neutral File.");
    os.close();
  }
  return ret;
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
  if (_fileOpenRead(neutralFilename, "MeshSpherical", is, verbose))
  {
    mesh = new MeshSpherical;
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
  bool flag_approx = true;
 
  /* Initializations */
  
  cs* Atriplet   = nullptr;
  cs* A          = nullptr;
  int nmeshes    = getNMeshes();
  int nvertex    = getNApices();
  int ncorner    = getNApexPerMesh();
  int ndim       = getNDim();

  // Preliminary checks 

  if (isCompatibleDb(db)) return nullptr;

  /* Core allocation */

  Atriplet = cs_spalloc(0, 0, 1, 1, 1);
  if (Atriplet == nullptr) return nullptr;
  VectorDouble weight(ncorner,0);
  VectorDouble units = _defineUnits();
  
  /* Loop on the samples */

  int ip_max = 0;
  int iech_max = 0;
  for (int iech=0; iech<db->getSampleNumber(); iech++)
  {
    if (! db->isActive(iech)) continue;
    if (iech > iech_max) iech_max = iech;
    VectorDouble coorLongLat = db->getSampleCoordinates(iech);
    
    /* Loop on the meshes */
    
    int found = -1;
    for (int imesh=0; imesh<nmeshes && found < 0; imesh++)
    {
      if (! _coorInMesh(coorLongLat,imesh,units[imesh],weight,flag_approx)) continue;

      /* Store the items in the sparse matrix */

      for (int icorn=0; icorn<ncorner; icorn++)
      {
        int ip = getApex(imesh,icorn);
        if (ip > ip_max) ip_max = ip;
        if (! cs_entry(Atriplet,iech,ip,weight[icorn]))
        {
          Atriplet = cs_spfree(Atriplet);
          return nullptr;
        }
      }
      found = imesh;
    }

    /* Printout if a point does not belong to any mesh */

    if (found < 0)
    {
      messerr("Point %d does not belong to any mesh",iech+1);
      for (int idim=0; idim<ndim; idim++)
        messerr(" Coordinate #%d = %lf",idim+1,coorLongLat[idim]);
      if (fatal)
      {
        Atriplet = cs_spfree(Atriplet);
        return nullptr;
      }
    }
  }
  
  /* Add the extreme value to force dimension */

  if (ip_max < nvertex - 1)
  {
    if (! cs_entry(Atriplet,iech_max,nvertex-1,0.))
    {
      Atriplet = cs_spfree(Atriplet);
      return nullptr;
    };
  }
  
  /* Convert the triplet into a sparse matrix */

  A = cs_triplet(Atriplet);
  Atriplet = cs_spfree(Atriplet);
  return(A);
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
  return _meshes.getValue(imesh,rank) - 1;
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

void MeshSpherical::getEmbeddedCoorPerMesh(int imesh, int ic, VectorDouble& coords) const
{
  /* The Variety is defined in the Global Environment */
  /* The required radius is set to the radius of Earth (6371m) */

  int variety_sphere;
  double r;
  variety_query(&variety_sphere);
  if (variety_sphere == 1)
  {
    variety_get_characteristics(&r);
  }
  else
  {
    r = EARTH_RADIUS;
  }
  util_convert_sph2cart(getCoor(imesh, ic, 0)-180., getCoor(imesh, ic, 1),
                        &coords[0], &coords[1], &coords[2],
                        r);

}

void MeshSpherical::getEmbeddedCoorPerApex(int iapex, VectorDouble& coords) const
{
  int variety_sphere;
  double r;
  variety_query(&variety_sphere);
  if (variety_sphere == 1)
  {
    variety_get_characteristics(&r);
  }
  else
  {
    r = EARTH_RADIUS;
  }
  util_convert_sph2cart(getApexCoor(iapex, 0)-180., getApexCoor(iapex, 1),
                        &coords[0], &coords[1], &coords[2],
                        r);
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
** \param[in]  meshsize  Dimension of the mesh (approx)
** \param[in]  flag_approx Approcimation flag
**
** \param[out] weights   Array of barycentric weights (Dim: NApexPerMesh)
** 
** \remark If flag_approx is False, calculations are performed with geodetic
** \remark distances
** \remark Otherwise calculations are performed by simply interpolating
** \remark in the longitude-latitude space
**
*****************************************************************************/
bool MeshSpherical::_coorInMesh(const VectorDouble& coor,
                                int imesh,
                                double meshsize,
                                VectorDouble& weights,
                                bool flag_approx) const
{
  VectorVectorDouble corners = getCoordinatesPerMesh(imesh);

  if (! flag_approx)
  {
    return (is_in_spherical_triangle_optimized(coor.data(),
                                               corners[0].data(),
                                               corners[1].data(),
                                               corners[2].data(),
                                               weights.data()));
  }
  else
  {

    // Round the angles with respect to the target coordinates

    for (int i = 0; i < 3; i++)
      corners[i][0] = _closestValue(coor[0], corners[i][0], 360.);
    return _weightsInMesh(coor, corners, meshsize, weights);
  }
}

double MeshSpherical::_closestValue(double ref, double coor, double period) const
{
  double dref = ABS(coor - ref);
  double d1 = ABS(coor - period - ref);
  if (d1 < dref)
    return coor - period;
  else
    return coor;
}

int MeshSpherical::_recopy(const MeshSpherical &m)
{
  _apices = m._apices;
  _meshes = m._meshes;
  return(0);
}

int MeshSpherical::_deserialize(std::istream& is, bool /*verbose*/)
{
  int ndim, napices, nmeshes, napexpermesh;

  bool ret = true;
  ret = ret && _recordRead<int>(is, "Space Dimension", ndim);
  ret = ret && _recordRead<int>(is, "Napices", napices);
  ret = ret && _recordRead<int>(is, "Number of Apices per Mesh", napexpermesh);
  ret = ret && _recordRead<int>(is, "Number of Meshes", nmeshes);

  VectorDouble apices_local;
  ret = ret && _recordReadVec<double>(is, "Apices", apices_local, ndim * napices);
  _apices = MatrixRectangular(napices, ndim);
  _apices.setValues(apices_local);

  VectorInt meshes_local;
  ret = ret && _recordReadVec<int>(is, "Meshes", meshes_local, nmeshes * napexpermesh);
  _meshes = MatrixInt(nmeshes, napexpermesh);
  _meshes.setValues(meshes_local);

  return 0;
}

int MeshSpherical::_serialize(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;

  ret = ret && _recordWrite<int>(os, "Space Dimension", getNDim());
  ret = ret && _recordWrite<int>(os, "Napices", getNApices());
  ret = ret && _recordWrite<int>(os, "Number of Apices per Mesh", getNApexPerMesh());
  ret = ret && _recordWrite<int>(os, "Number of Meshes", getNMeshes());

  ret = ret && _recordWriteVec<double>(os, "Apices", _apices.getValues());
  ret = ret && _recordWriteVec<int>(os, "Meshes", _meshes.getValues());
  return 0;
}
