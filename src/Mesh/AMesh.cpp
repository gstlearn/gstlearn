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
#include "Mesh/AMesh.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Db/Db.hpp"
#include "Basic/Vector.hpp"
#include "Space/SpacePoint.hpp"

#include <algorithm>

AMesh::AMesh()
  : AStringable()
  , _nDim(0)
  , _extendMin()
  , _extendMax()
{
}

AMesh::AMesh(const AMesh &m)
  : AStringable(m)
{
  _recopy(m);
}

AMesh& AMesh::operator= (const AMesh &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _recopy(m);
  }
  return *this;
}

AMesh::~AMesh()
{

}

int AMesh::setExtend(const VectorDouble extendmin,
                     const VectorDouble extendmax)
{
  _extendMin = extendmin;
  _extendMax = extendmax;
  return(0);
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
** \remarks This function only makes sens when the Meshing is built
** \remarks based on two Dbs where samples may coincide
**
*****************************************************************************/
void AMesh::getDuplicates(int   /*verbose*/,
                          Db* /*dbin*/,
                          Db* /*dbout*/,
                          int  *nbdupl,
                          int **dupl1,
                          int **dupl2) const
{
  *nbdupl = 0;
  *dupl1 = nullptr;
  *dupl2 = nullptr;
}

String AMesh::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

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
    sstr << toTitle(1, "Bounding Box Extension");
    for (int idim = 0; idim < _nDim; idim++)
      sstr << "Dim #" << idim+1 << " - Min:" << _extendMin[idim] << " - Max:" <<
              _extendMax[idim] << std::endl;
  }

  AStringFormat sf;
  if (strfmt != nullptr) sf = *strfmt;
  if (sf.getLevel() > 0)
  {
    MatrixRectangular apices;
    VectorInt meshes;
    getElements(apices, meshes);

    sstr << "List of Apices" << std::endl;
    sstr << apices.toString(strfmt);

    sstr << toMatrix("List of Meshes", VectorString(), VectorString(), false,
                    getNApexPerMesh(), getNMeshes(), meshes);
  }
  return sstr.str();
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
int AMesh::isCompatibleDb(const Db *db) const
{
  if (getNDim() == db->getNDim()) return 0;

  messerr("Inconsistent Space dimension between Meshing (%d) and Db (%d)",
          getNDim(),db->getNDim());
  return 1;
}

VectorDouble AMesh::getMeshSizes() const
{
  VectorDouble units;
  for (int imesh = 0; imesh < getNMeshes(); imesh++)
    units.push_back(getMeshSize(imesh));
  return units;
}

void AMesh::printMeshes(int imesh0) const
{
  mestitle(0,"Mesh Information");
  message("- Number of Meshes = %d\n",getNMeshes());
  message("- Number of Apices = %d\n",getNApices());

  int ideb = (imesh0 >= 0) ? imesh0 : 0;
  int ifin = (imesh0 >= 0) ? imesh0 + 1: getNMeshes();
  for (int imesh=ideb; imesh<ifin; imesh++)
  {
    message("Mesh #%d\n",imesh+1);
    for (int icorn=0; icorn<getNApexPerMesh(); icorn++)
    {
      message("Point #%d");
      for (int idim=0; idim<getNDim(); idim++)
      {
        message(" %lf",getCoor(imesh,icorn,idim));
      }
      message("\n");
    }
  }
}

void AMesh::_recopy(const AMesh &m)
{
  _nDim    = m._nDim;
  _extendMin = m._extendMin;
  _extendMax = m._extendMax;
}

/**
 * Convert the Meshing (New class) into a meshing (Old class)
 * @param a_mesh Meshing characteristics (New class)
 * @return Pointer to the Meshing (Old class)
 */
SPDE_Mesh* AMesh::_convertToOldMesh(AMesh* a_mesh) const
{
  SPDE_Mesh *s_mesh;
  int number;

  s_mesh = spde_mesh_manage(1, NULL);
  if (s_mesh == nullptr) return (s_mesh);
  s_mesh->ndim = a_mesh->getNDim();
  s_mesh->ncorner = a_mesh->getNApexPerMesh();
  s_mesh->nmesh = a_mesh->getNMeshes();
  s_mesh->nvertex = a_mesh->getNApices();

  // Retrieve the elements
  MatrixRectangular points;
  VectorInt meshes;
  a_mesh->getElements(points, meshes);

  number = points.getNTotal();
  VectorDouble vect = points.getValues();
  s_mesh->points = (double *) mem_alloc(sizeof(double) * number,1);
  for (int i=0; i<number; i++) s_mesh->points[i] = vect[i];

  number = static_cast<int> (meshes.size());
  s_mesh->meshes = (int *) mem_alloc(sizeof(int) * number,1);
  for (int i=0; i<number; i++) s_mesh->meshes[i] = meshes[i];

  return s_mesh;
}

VectorDouble AMesh::getCoordinates(int idim) const
{
  if (! _isSpaceDimensionValid(idim)) return VectorDouble();
  int np = getNApices();
  VectorDouble coor(np);
  for (int ip = 0; ip < np; ip++)
    coor[ip] = getApexCoor(ip, idim);
  return coor;
}

VectorInt AMesh::getMeshByApexPair(int apex1, int apex2) const
{
  VectorInt list;
  int ncorner = getNApexPerMesh();
  int found, apex0;

  for (int imesh = 0; imesh < getNMeshes(); imesh++)
  {
    found = 0;
    for (int ic = 0; ic < ncorner; ic++)
    {
      apex0 = getApex(imesh,ic);
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
void AMesh::getElements(MatrixRectangular& apices, VectorInt& meshes) const
{
  int nmeshes = getNMeshes();
  int ndim    = getNDim();
  int napices = getNApices();
  int ncorner = getNApexPerMesh();

  // Dimension the returned containers

  apices.reset(napices,ndim);
  meshes.clear();

  // Load the Apices

  for (int i = 0; i < napices; i++)
    for (int idim = 0; idim < ndim; idim++)
      apices.setValue(i,idim,getApexCoor(i, idim));

  // Load the Meshes

  for (int imesh = 0; imesh < nmeshes; imesh++)
    for (int icorner= 0; icorner < ncorner; icorner++)
      meshes.push_back(getApex(imesh, icorner));
 }

bool AMesh::_isSpaceDimensionValid(int idim) const
{
  if (idim < 0 || idim >= _nDim)
  {
    mesArg("SPace Dimension Index",idim,_nDim);
    return false;
  }
  return true;
}

VectorDouble AMesh::getExtrema(int idim) const
{
  VectorDouble vec(2);
  vec[0] = getExtendMin(idim);
  vec[1] = getExtendMax(idim);
  return vec;
}

VectorDouble AMesh::getCoordinatesPerMesh(int imesh, int idim, bool flagClose) const
{
  VectorDouble vec;
  int ncorner = getNApexPerMesh();
  if (flagClose)
    vec.resize(ncorner+1);
  else
    vec.resize(ncorner);

  for (int ic = 0; ic < ncorner; ic++)
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
void AMesh::getEmbeddedCoorPerMesh(int imesh, int ic, VectorDouble& coords) const
{
  for (int idim = 0; idim < _nDim; idim++)
    coords[idim] = getCoor(imesh, ic, idim);
}

/**
 * Fill the coordinates of an apex in embedded space
 * @param iapex  Apex index
 * @param coords Array of coordinates
 */
void AMesh::getEmbeddedCoorPerApex(int iapex, VectorDouble& coords) const
{
  for (int idim = 0; idim < _nDim; idim++)
    coords[idim] = getApexCoor(iapex, idim);
}


/**
 * Fill the array of coordinates of all apices of a mesh in embedded space
 * Storage [ndim, ncorner]
 * @param imesh Mesh rank
 * @param vec   Returned array
 */
void AMesh::getEmbeddedCoordinatesPerMesh(int imesh, VectorVectorDouble& vec) const
{
  int ncorner = getNApexPerMesh();

  for (int ic = 0; ic < ncorner; ic++)
    getEmbeddedCoorPerMesh(imesh, ic, vec[ic]);
}

/**
 * Returns the array of coordinates of all apices of a mesh in embedded space
 * @param imesh Mesh rank
 * @return
 */
VectorVectorDouble AMesh::getEmbeddedCoordinatesPerMesh(int imesh) const
{
  int ndim = getEmbeddedNDim();
  int ncorner = getNApexPerMesh();
  VectorVectorDouble vec(ncorner);
  for (auto &e: vec)
    e = VectorDouble(ndim);

  getEmbeddedCoordinatesPerMesh(imesh, vec);
  return vec;
}

/**
 * Returns the coordinates of the Mesh apices expressed in the embedded space
 * The returned vector is organized by coordinate
 * @return
 */
VectorVectorDouble AMesh::getEmbeddedApexCoordinates() const
{
  int ndim = getEmbeddedNDim();
  int napices = getNApices();
  VectorVectorDouble vec(ndim);
  for (auto &e: vec)
    e = VectorDouble(napices);

  VectorDouble local(ndim);
  for (int ip = 0; ip < napices; ip++)
  {
    getEmbeddedCoorPerApex(ip, local);
    for (int idim = 0; idim < ndim; idim++)
      vec[idim][ip] = local[idim];
  }
  return vec;
}

VectorDouble AMesh::getApexCoordinates(int iapex) const
{
  VectorDouble vec(_nDim);
  for (int idim = 0; idim < _nDim; idim++)
    vec[idim] = getApexCoor(iapex, idim);
  return vec;
}

VectorDouble AMesh::getDistance(int imesh0, const VectorInt& jmeshs)
{
  VectorDouble vec;
  VectorInt jlocal = jmeshs;
  if (jlocal.empty()) jlocal = ut_ivector_sequence(getNMeshes());
  int number = (int) jlocal.size();

  SpacePoint P1 = SpacePoint(getApexCoordinates(imesh0));

  for (int jmesh = 0; jmesh < number; jmesh++)
  {
    SpacePoint P2 = SpacePoint(getApexCoordinates(jmesh));
    vec[jmesh] = P1.getDistance(P2);
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
std::vector<VectorInt> AMesh::getNeighborhoodPerMesh() const
{

  int napices  = getNApices();
  int nmeshes  = getNMeshes();
  int npermesh = getNApexPerMesh();

  std::vector<VectorInt> Vmesh;
  Vmesh.resize(napices);

  // Loop on the meshes

  for (int imesh = 0; imesh < nmeshes; imesh++)
  {
    // Loop on the apices for each mesh

    for (int rank = 0; rank < npermesh; rank++)
    {
      int ip = getApex(imesh, rank);
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
std::vector<VectorInt> AMesh::getNeighborhoodPerApex() const
{
  int napices  = getNApices();
  int npermesh = getNApexPerMesh();

  std::vector<VectorInt> Vapex;
  Vapex.resize(napices);

  // Elaborate the Meshing Neighborhood by Mesh first
  std::vector<VectorInt> Vmesh = getNeighborhoodPerMesh();

  for (int ip = 0; ip < napices; ip++)
  {
    VectorInt vec;

    // Loop on the meshes neighboring the target apex 'ip'
    int nmesh = (int) Vmesh[ip].size();
    for (int i = 0; i < nmesh; i++)
    {
      // Index of the neighboring mesh
      int imesh = Vmesh[ip][i];

      // Loop on the apices of the neighboring mesh
      for (int rank = 0; rank < npermesh; rank++)
      {
        int jp = getApex(imesh, rank);

        // Skip the target apex itself
        if (jp == ip) continue;

        // Add the new apex to the list
        vec.push_back(jp);
      }
    }

    // Sort the list and suppress the duplicates
    std::sort( vec.begin(), vec.end() );
    vec.erase( std::unique( vec.begin(), vec.end() ), vec.end() );

    // Store the list
    Vapex[ip] = vec;
  }
  return Vapex;
}

void AMesh::dumpNeighborhood(std::vector<VectorInt>& Vmesh)
{
  mestitle(1,"List of Meshing Neighborhood");
  int nrow = (int) Vmesh.size();
  for (int irow = 0; irow < nrow; irow++)
  {
    ut_ivector_display(String(), Vmesh[irow]);
  }
}
