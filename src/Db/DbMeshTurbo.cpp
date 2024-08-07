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
#include "Db/DbMeshTurbo.hpp"
#include "Db/DbStringFormat.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"

DbMeshTurbo::DbMeshTurbo(const VectorInt& nx,
                         const VectorDouble& dx,
                         const VectorDouble& x0,
                         const VectorDouble& angles,
                         const ELoadBy& order,
                         const VectorDouble& tab,
                         const VectorString& names,
                         const VectorString& locatorNames,
                         bool flag_polarized,
                         bool verbose,
                         int mode)
  : DbGrid()
  , _mesh(nx, dx, x0, angles, flag_polarized, verbose, mode)
{
  (void)reset(nx, dx, x0, angles, order, tab, names, locatorNames, true, false);
}

DbMeshTurbo::DbMeshTurbo(const DbMeshTurbo& r)
  : DbGrid(r)
  , _mesh(r._mesh)
{
}

DbMeshTurbo& DbMeshTurbo::operator=(const DbMeshTurbo& r)
{
  if (this != &r)
  {
    DbGrid::operator=(r);
    _mesh = r._mesh;
  }
  return *this;
}

DbMeshTurbo::~DbMeshTurbo()
{
}

String DbMeshTurbo::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  const DbStringFormat* dbfmt = dynamic_cast<const DbStringFormat*>(strfmt);
  DbStringFormat dsf;
  if (dbfmt != nullptr) dsf = *dbfmt;

  sstr << toTitle(0, "Data Base for Turbo Meshing");

  sstr << _toStringCommon(&dsf);

  if (dsf.matchResume())
  {
    sstr << _summaryString();
  }

  sstr << _mesh.toString(strfmt);

  return sstr.str();
}

DbMeshTurbo* DbMeshTurbo::create(const VectorInt& nx,
                                 const VectorDouble& dx,
                                 const VectorDouble& x0,
                                 const VectorDouble& angles,
                                 const ELoadBy& order,
                                 const VectorDouble& tab,
                                 const VectorString& names,
                                 const VectorString& locatorNames,
                                 bool flag_polarized,
                                 bool verbose)
{
  // Creating the MeshETurbo internal storage
  DbMeshTurbo* dbmesh = new DbMeshTurbo(nx, dx, x0, angles, order, tab, names,
                                        locatorNames, flag_polarized, verbose);
  if (dbmesh == nullptr)
  {
    messerr("Error when creating DbMeshTurbo from Samples");
    delete dbmesh;
    return nullptr;
  }
  return dbmesh;
}

bool DbMeshTurbo::_deserialize(std::istream& is, bool verbose)
{
  int ndim = 0;
  bool ret = true;

  // Reading the header

  ret      = ret && _recordRead<int>(is, "Space Dimension", ndim);

  // Reading the meshing information

  ret      = ret && _mesh.deserialize(is);

  // Reading the Db information

  ret      = ret && DbGrid::_deserialize(is, verbose);

  return ret;
}

bool DbMeshTurbo::_serialize(std::ostream& os, bool verbose) const
{
  bool ret = true;

  /* Writing the header */

  ret      = ret && _recordWrite<int>(os, "Space Dimension", getNDim());

  // Writing the Meshing information

  ret      = ret && _mesh.serialize(os);

  /* Writing the tail of the file */

  ret      = ret && DbGrid::_serialize(os, verbose);

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
DbMeshTurbo* DbMeshTurbo::createFromNF(const String& neutralFilename, bool verbose)
{
  DbMeshTurbo* dbmesh = new DbMeshTurbo;
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
bool DbMeshTurbo::isConsistent() const
{
  // Check on the count of addresses
  int nech = getSampleNumber();
  if (_mesh.getNApices() > nech)
  {
    messerr("Number of meshes (%d)", _mesh.getNApices());
    messerr("must not be larger than Sample Number (%d)", nech);
    return false;
  }
  return true;
}
