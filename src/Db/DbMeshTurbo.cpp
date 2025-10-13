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
#include "Db/DbMeshTurbo.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Basic/VectorNumT.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"

namespace gstlrn
{
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
                         Id mode)
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

  const auto* dbfmt = dynamic_cast<const DbStringFormat*>(strfmt);
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
  auto* dbmesh = new DbMeshTurbo(nx, dx, x0, angles, order, tab, names,
                                        locatorNames, flag_polarized, verbose);
  if (dbmesh == nullptr)
  {
    messerr("Error when creating DbMeshTurbo from Samples");
    delete dbmesh;
    return nullptr;
  }
  return dbmesh;
}

bool DbMeshTurbo::_deserializeAscii(std::istream& is, bool verbose)
{
  Id ndim = 0;
  bool ret = true;

  // Reading the header

  ret = ret && _recordRead<Id>(is, "Space Dimension", ndim);

  // Reading the meshing information

  ret = ret && _mesh._deserializeAscii(is);

  // Reading the Db information

  ret = ret && DbGrid::_deserializeAscii(is, verbose);

  return ret;
}

bool DbMeshTurbo::_serializeAscii(std::ostream& os, bool verbose) const
{
  bool ret = true;

  /* Writing the header */

  ret = ret && _recordWrite<Id>(os, "Space Dimension", getNDim());

  // Writing the Meshing information

  ret = ret && _mesh._serializeAscii(os);

  /* Writing the tail of the file */

  ret = ret && DbGrid::_serializeAscii(os, verbose);

  return ret;
}

/**
 * Create a DbMesh by loading the contents of a Neutral File
 *
 * @param NFFilename Name of the Neutral File (DbMesh format)
 * @param verbose    Verbose
 *
 * @remarks The name does not need to be completed in particular when defined by absolute path
 * @remarks or read from the Data Directory (in the gstlearn distribution)
 */
DbMeshTurbo* DbMeshTurbo::createFromNF(const String& NFFilename, bool verbose)
{
  DbMeshTurbo* dbmesh = new DbMeshTurbo;
  if (dbmesh->_fileOpenAndDeserialize(NFFilename, verbose)) return dbmesh;
  delete dbmesh;
  return nullptr;
}

/**
 * @brief Check if the contents of private member of this class is compatible
 * with the number of samples stored in the DbGrid
 * @return true if everything is OK; false if a problem occurs
 */
bool DbMeshTurbo::isConsistent() const
{
  // Check on the count of addresses
  auto nech = getNSample();
  if (_mesh.getNApices() > nech)
  {
    messerr("Number of meshes (%d)", _mesh.getNApices());
    messerr("must not be larger than Sample Number (%d)", nech);
    return false;
  }
  return true;
}
#ifdef HDF5
bool DbMeshTurbo::_deserializeH5(H5::Group& grp, [[maybe_unused]] bool verbose)
{
  auto dbg = SerializeHDF5::getGroup(grp, "DbMeshTurbo");
  if (!dbg)
  {
    return false;
  }

  /* Read the grid characteristics */
  bool ret = true;
  Id ndim = 0;

  ret = ret && SerializeHDF5::readValue(*dbg, "NDim", ndim);

  // Writing the Meshing information
  ret = ret && _mesh._deserializeH5(*dbg, verbose);

  /* Writing the tail of the file */
  ret = ret && DbGrid::_deserializeH5(*dbg, verbose);

  return ret;
}

bool DbMeshTurbo::_serializeH5(H5::Group& grp, [[maybe_unused]] bool verbose) const
{
  auto dbG = grp.createGroup("DbMeshTurbo");

  bool ret = true;

  ret = ret && SerializeHDF5::writeValue(dbG, "NDim", getNDim());

  // Writing the Meshing information
  ret = ret && _mesh._serializeH5(dbG, verbose);

  /* Writing the tail of the file */
  ret = ret && DbGrid::_serializeH5(dbG, verbose);

  return ret;
}
#endif
}