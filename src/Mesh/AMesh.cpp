/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "geoslib_old_f.h"
#include "Mesh/AMesh.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixInt.hpp"
#include "Db/Db.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/AStringable.hpp"
#include "Space/SpacePoint.hpp"

#include <algorithm>

AMesh::AMesh()
    : AStringable(),
      ASerializable(),
      _nDim(0),
      _extendMin(),
      _extendMax()
{
}

AMesh::AMesh(const AMesh &m)
  : AStringable(m),
    ASerializable(m)
{
  _recopy(m);
}

AMesh& AMesh::operator= (const AMesh &m)
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

}

int AMesh::_setExtend(const VectorDouble extendmin,
                     const VectorDouble extendmax)
{
  _extendMin = extendmin;
  _extendMax = extendmax;
  return(0);
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
    sstr << toTitle(1, "Bounding Box Extension");
    for (int idim = 0; idim < _nDim; idim++)
      sstr << "Dim #" << idim+1 << " - Min:" << _extendMin[idim] << " - Max:" <<
              _extendMax[idim] << std::endl;
  }

  AStringFormat sf;
  if (strfmt != nullptr) sf = *strfmt;
  if (sf.getLevel() > 1)
  {
    MatrixRectangular apices;
    MatrixInt meshes;
    getElements(apices, meshes);

    sstr << "List of Apices" << std::endl;
    sstr << apices.toString(strfmt);
    sstr << "List of Meshes" << std::endl;
    sstr << meshes.toString(strfmt);
  }
  return sstr.str();
}

void AMesh::getCoordinatesInPlace(int imesh, int rank, VectorDouble& coords) const
{
  for (int idim = 0; idim < getNDim(); idim++)
    coords[idim] = getCoor(imesh, rank, idim);
}

void AMesh::getApexCoordinatesInPlace(int i, VectorDouble& coords) const
{
  for (int idim = 0; idim < getNDim(); idim++)
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

void AMesh::printMesh(int imesh0) const
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
      message("Point #%d",getApex(imesh, icorn));
      for (int idim=0; idim<getNDim(); idim++)
      {
        message(" %lf",getCoor(imesh,icorn,idim));
      }
      message("\n");
    }
  }
}

void AMesh::printMeshes(int level, int nline_max) const
{
  mestitle(0,"Mesh Information");
  message("- Number of Meshes = %d\n",getNMeshes());
  message("- Number of Apices = %d\n",getNApices());

  if (level == 0) return;

  if (level == 1)
    _printMeshListByIndices(nline_max);

  if (level == 2)
    _printMeshListByCoordinates(nline_max);
}

void AMesh::_recopy(const AMesh &m)
{
  _nDim    = m._nDim;
  _extendMin = m._extendMin;
  _extendMax = m._extendMax;
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

VectorVectorDouble AMesh::getAllCoordinates() const
{
  int napices = getNApices();
  VectorDouble local(_nDim);
  VectorVectorDouble coords(_nDim);
  for (int idim = 0; idim < _nDim; idim++)
    coords[idim].resize(napices);

  for (int ip = 0; ip < napices; ip++)
  {
    getApexCoordinatesInPlace(ip, local);
    for (int idim = 0; idim < _nDim; idim++)
      coords[idim][ip] = local[idim];
  }
  return coords;
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
void AMesh::getElements(MatrixRectangular& apices, MatrixInt& meshes) const
{
  int nmeshes = getNMeshes();
  int ndim    = getNDim();
  int napices = getNApices();
  int ncorner = getNApexPerMesh();

  // Dimension the returned containers

  apices.reset(napices,ndim);
  meshes.reset(nmeshes, ncorner);

  // Load the Apices

  VectorDouble local(ndim);
  for (int i = 0; i < napices; i++)
  {
    getApexCoordinatesInPlace(i, local);
    for (int idim = 0; idim < ndim; idim++)
      apices.setValue(i,idim,local[idim]);
  }

  // Load the Meshes

  for (int imesh = 0; imesh < nmeshes; imesh++)
    for (int icorner= 0; icorner < ncorner; icorner++)
      meshes.setValue(imesh, icorner, getApex(imesh, icorner));

  return;
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
  getCoordinatesInPlace(imesh, ic, coords);
}

/**
 * Fill the coordinates of an apex in embedded space
 * @param iapex  Apex index
 * @param coords Array of coordinates
 */
void AMesh::getEmbeddedCoorPerApex(int iapex, VectorDouble& coords) const
{
  getApexCoordinatesInPlace(iapex, coords);
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

VectorVectorDouble AMesh::getCoordinatesPerMesh(int imesh) const
{
  int ndim = getNDim();
  int ncorner = getNApexPerMesh();
  VectorVectorDouble vec(ncorner);
  for (auto &e: vec)
    e = VectorDouble(ndim);

  for (int ic = 0; ic < ncorner; ic++)
    for (int idim = 0; idim < ndim; idim++)
      vec[ic][idim] = getCoor(imesh, ic, idim);
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

VectorDouble AMesh::getDistances(int iapex0, const VectorInt& japices)
{
  VectorInt jlocal = japices;
  if (jlocal.empty()) jlocal = VH::sequence(getNApices());
  int number = (int) jlocal.size();
  VectorDouble vec(number,0.);

  SpacePoint P1(getApexCoordinates(iapex0));

  for (int iapex = 0; iapex < number; iapex++)
  {
    SpacePoint P2(getApexCoordinates(jlocal[iapex]));
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
    VH::display(String(), Vmesh[irow]);
  }
}

bool AMesh::_deserialize(std::istream& is, bool /*verbose*/)
{
  bool ret = true;
  ret = ret && _recordRead<int>(is, "Space Dimension", _nDim);
  ret = ret && _recordReadVec<double>(is, "Minimum Extension", _extendMin, _nDim);
  ret = ret && _recordReadVec<double>(is, "Maximum Extension", _extendMax, _nDim);
  return ret;
}

bool AMesh::_serialize(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;
  ret = ret && _recordWrite<int>(os, "Space Dimension", getNDim());
  ret = ret && _recordWriteVec<double>(os, "Minimum Extension", _extendMin);
  ret = ret && _recordWriteVec<double>(os, "Maximum Extension", _extendMax);
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
  int ncorner = getNApexPerMesh();
  int ndim    = getNDim();

  // Loop on the vertices
  double total = 0.;
  for (int icorn=0; icorn<ncorner; icorn++)
  {

    // Build the determinant
    MatrixSquareGeneral mat(ndim);
    int kcorn = 0;
    for (int jcorn=0; jcorn<ncorner; jcorn++)
    {
      if (icorn == jcorn) continue;
      for (int idim=0; idim<ndim; idim++)
        mat.setValue(idim,kcorn,corners[jcorn][idim] - coor[idim]);
      kcorn++;
    }
    double ratio = ABS(mat.determinant()) / meshsize / FACDIM[ndim];
    if (ratio < -eps || ratio > 1 + eps) return false;
    weights[icorn] = ratio;
    total += ratio;
  }
  if (ABS(total - 1) > eps) return false;
  return true;
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
  int ndim    = getNDim();
  int ncorner = getNApexPerMesh();

  // Calculate the mesh size
  MatrixSquareGeneral mat;
  mat.reset(ndim,ndim);
  for (int icorn=1; icorn<ncorner; icorn++)
    for (int idim=0; idim<ndim; idim++)
      mat.setValue(icorn-1,idim, corners[icorn][idim] - corners[0][idim]);
  unit = ABS(mat.determinant()) / facdim[ndim];

  return unit;
}

void AMesh::_printMeshListByCoordinates(int nline_max) const
{
  int ndim    = getNDim();
  int nmesh   = getNMeshes();
  int ncorner = getNApexPerMesh();

  int iline = 0;
  for (int imesh=0; imesh<nmesh; imesh++)
  {
    message("Mesh #%5d/%5d\n",imesh+1,nmesh);
    for (int icorn=0; icorn<ncorner; icorn++)
    {
      message(" Apex %4d: ",getApex(imesh,icorn));
      for (int idim=0; idim<ndim; idim++)
        message(" %lf",getCoor(imesh,icorn,idim));
      message("\n");
    }

    iline++;
    if (nline_max > 0 && iline >= nline_max) return;
  }
}

void AMesh::_printMeshListByIndices(int nline_max) const
{
  int nmesh   = getNMeshes();
  int ncorner = getNApexPerMesh();

  int iline = 0;
  for (int imesh=0; imesh<nmesh; imesh++)
  {
    message("Mesh #%d/%d: ",imesh+1,nmesh);
    for (int icorn=0; icorn<ncorner; icorn++)
      message(" %d",getApex(imesh,icorn));
    message("\n");

    iline++;
    if (nline_max > 0 && iline >= nline_max) return;
  }
}

