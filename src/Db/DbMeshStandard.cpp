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
#include "Db/Db.hpp"
#include "Db/DbMeshStandard.hpp"
#include "Db/DbStringFormat.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"

DbMeshStandard::DbMeshStandard(int ndim,
                               int napexpermesh,
                               const VectorDouble& apices,
                               const VectorInt& meshes,
                               const ELoadBy& order,
                               const VectorDouble& tab,
                               const VectorString& names,
                               const VectorString& locatorNames,
                               bool verbose)
  : Db()
  , _mesh()
{
  _mesh.reset(ndim, napexpermesh, apices, meshes, true, verbose);
  int napices = _mesh.getNApices();

  addColumns(apices, "x", ELoc::X, 0, false, 0., ndim);

  if (!tab.empty())
    (void)resetFromSamples(napices, order, tab, names, locatorNames);
}

DbMeshStandard::DbMeshStandard(const DbMeshStandard& r)
  : Db(r)
  , _mesh(r._mesh)
{
}

DbMeshStandard& DbMeshStandard::operator=(const DbMeshStandard& r)
{
  if (this != &r)
  {
    Db::operator=(r);
    _mesh = r._mesh;
  }
  return *this;
}

DbMeshStandard::~DbMeshStandard()
{
}

String DbMeshStandard::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  const DbStringFormat* dbfmt = dynamic_cast<const DbStringFormat*>(strfmt);
  DbStringFormat dsf;
  if (dbfmt != nullptr) dsf = *dbfmt;

  sstr << toTitle(0, "Data Base for Standard Meshing");

  sstr << _toStringCommon(&dsf);

  if (dsf.matchResume())
  {
    sstr << _summaryString();
  }

  sstr << _mesh.toString(strfmt);

  return sstr.str();
}

DbMeshStandard* DbMeshStandard::create(int ndim,
                                       int napexpermesh,
                                       const VectorDouble& apices,
                                       const VectorInt& meshes,
                                       const ELoadBy& order,
                                       const VectorDouble& tab,
                                       const VectorString& names,
                                       const VectorString& locatorNames,
                                       bool verbose)
{
  DbMeshStandard* dbmesh =
    new DbMeshStandard(ndim, napexpermesh, apices, meshes, order, tab, names,
                       locatorNames, verbose);
  if (dbmesh == nullptr)
  {
    messerr("Error when creating DbMeshStandard from Samples");
    delete dbmesh;
    return nullptr;
  }
  return dbmesh;
}

DbMeshStandard*
DbMeshStandard::createFromExternal(const MatrixRectangular& apices,
                                   const MatrixInt& meshes,
                                   const ELoadBy& order,
                                   const VectorDouble& tab,
                                   const VectorString& names,
                                   const VectorString& locatorNames,
                                   bool verbose)
{
  DbMeshStandard* dbmesh = new DbMeshStandard;
  dbmesh->_mesh.reset(apices, meshes, verbose);
  if (dbmesh == nullptr)
  {
    messerr("Error when creating DbMeshStandard from Samples");
    delete dbmesh;
    return nullptr;
  }
  if (!tab.empty())
  {
    int nech = apices.getNCols();
    (void)dbmesh->resetFromSamples(nech, order, tab, names, locatorNames);
  }
  return dbmesh;
}

bool DbMeshStandard::_deserialize(std::istream& is, bool verbose)
{
  int ndim = 0;
  bool ret = true;

  // Reading the header

  ret      = ret && _recordRead<int>(is, "Space Dimension", ndim);

  // Reading the meshing information

  ret      = ret && _mesh.deserialize(is);

  // Reading the Db information

  ret      = ret && Db::_deserialize(is, verbose);

  return ret;
}

bool DbMeshStandard::_serialize(std::ostream& os, bool verbose) const
{
  bool ret = true;

  /* Writing the header */

  ret      = ret && _recordWrite<int>(os, "Space Dimension", getNDim());

  // Writing the Meshing information

  ret      = ret && _mesh.serialize(os);

  /* Writing the tail of the file */

  ret      = ret && Db::_serialize(os, verbose);

  return ret;
}

/**
 * Create a DbMesh by loading the contents of a Neutral File
 *
 * @param neutralFilename Name of the Neutral File (DbMesh format)
 * @param verbose         Verbose
 *
 * @remarks The name does not need to be completed in particular when defined by absolute path
 * @remarks or read from the Data Directory (in the gstlearn distribution)
 */
DbMeshStandard* DbMeshStandard::createFromNF(const String& neutralFilename, bool verbose)
{
  DbMeshStandard* dbmesh = new DbMeshStandard;
  std::ifstream is;
  bool success = false;
  if (dbmesh->_fileOpenRead(neutralFilename, is, verbose))
  {
    success = dbmesh->deserialize(is, verbose);
  }
  if (! success)
  {
    delete dbmesh;
    dbmesh = nullptr;
  }
  return dbmesh;
}

/**
 * @brief Check if the contents of private member of this class is compatible
 * with the number of samples stored in the DbGrid
 * @return true if everything is OK; false if a problem occurs
 */
bool DbMeshStandard::isConsistent() const
{
  // Check on the count of addresses
  int nech = getNSample();
  if (_mesh.getNApices() > nech)
  {
    messerr("Number of meshes (%d)", _mesh.getNApices());
    messerr("must not be larger than Sample Number (%d)", nech);
    return false;
  }
  return true;
}

double DbMeshStandard::getCoor(int imesh, int rank, int idim) const
{
  return getCoordinate(_mesh.getApex(imesh, rank), idim);
}
void DbMeshStandard::getCoordinatesPerMeshInPlace(int imesh, int rank, VectorDouble& coords) const
{
  for (int idim = 0; idim < getNDim(); idim++)
    coords[idim] = getCoor(imesh, rank, idim);
}
double DbMeshStandard::getApexCoor(int i, int idim) const
{
  return getCoordinate(i, idim);
}
void DbMeshStandard::getApexCoordinatesInPlace(int i, VectorDouble& coords) const
{
  for (int idim = 0; idim < getNDim(); idim++)
    coords[idim] = getApexCoor(i, idim);
}
VectorDouble DbMeshStandard::getCoordinatesPerMesh(int imesh, int idim, bool flagClose) const
{
  VectorDouble vec;
  int ncorner = _mesh.getNApexPerMesh();

  if (flagClose)
    vec.resize(ncorner + 1);
  else
    vec.resize(ncorner);

  for (int ic = 0; ic < ncorner; ic++) vec[ic] = getCoor(imesh, ic, idim);
  if (flagClose) vec[ncorner] = getCoor(imesh, 0, idim);

  return vec;
}