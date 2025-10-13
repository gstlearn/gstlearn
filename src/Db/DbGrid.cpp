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
#include "Db/DbGrid.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Grid.hpp"
#include "Basic/Law.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/VectorNumT.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Estimation/CalcImage.hpp"
#include "Polygon/Polygons.hpp"
#include "Space/SpaceTarget.hpp"
#include "Stats/Classical.hpp"
#include "geoslib_define.h"

#include <algorithm>
#include <cmath>

namespace gstlrn
{

DbGrid::DbGrid()
  : Db()
  , _grid()
{
  _clear();
}

DbGrid::DbGrid(const DbGrid& r)
  : Db(r)
  , _grid(r._grid)
{
}

DbGrid& DbGrid::operator=(const DbGrid& r)
{
  if (this != &r)
  {
    Db::operator=(r);
    _grid = r._grid;
  }
  return *this;
}

DbGrid::~DbGrid()
{
}

String DbGrid::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  const auto* dbfmt = dynamic_cast<const DbStringFormat*>(strfmt);
  DbStringFormat dsf;
  if (dbfmt != nullptr) dsf = *dbfmt;

  sstr << toTitle(0, "Data Base Grid Characteristics");

  if (dsf.matchResume())
  {
    sstr << _summaryString();
    sstr << _grid.toString();
  }

  sstr << _toStringCommon(&dsf);

  return sstr.str();
}

/**
 * Creating a Db regular grid of any dimension
 *
 * @param nx            A vector of the number of grid meshes.
 *                      The number of items in this argument gives the dimension of the space.
 *                      (size = ndim)
 * @param dx            Vector cell meshes size in each direction (size = ndim) (by default, use 1)
 * @param x0            Vector of origin coordinates (size = ndim) (by default, use 0)
 * @param angles        Array giving the rotation angles (only for dimension 2 or 3).
 *                      The first angle corresponds to the rotation around OZ axis,
 *                      the second to a rotation around OY'and the third one around Ox.
 *                      The dimension of this array cannot exceed the space dimension.
 * @param order         Flag for values order in 'tab' (defined ELoadBy.hpp)
 * @param tab           Variable values array (size = nvar * nsamples)
 * @param names         Names of the Variables of 'tab' (size = nvar)
 * @param locatorNames  Locators for each variable of array 'tab' (size = nvar)
 * @param flagAddSampleRank If true, add an automatic rank variable
 * @param flagAddCoordinates If TRUE, add the grid coordinates
 */
Id DbGrid::reset(const VectorInt& nx,
                 const VectorDouble& dx,
                 const VectorDouble& x0,
                 const VectorDouble& angles,
                 const ELoadBy& order,
                 const VectorDouble& tab,
                 const VectorString& names,
                 const VectorString& locatorNames,
                 bool flagAddSampleRank,
                 bool flagAddCoordinates)
{
  _clear();

  Id ndim = static_cast<Id>(nx.size());
  Id nech = 1;
  for (Id idim = 0; idim < ndim; idim++)
    nech *= nx[idim];
  Id ntab   = (tab.empty()) ? 0 : static_cast<Id>(tab.size() / nech);
  Id number = 0;
  if (flagAddSampleRank) number += 1;
  if (flagAddCoordinates) number += ndim;
  Id ncol = number + ntab;

  // Create the grid

  if (gridDefine(nx, dx, x0, angles)) return 1;
  resetDims(ncol, nech);

  // Load the data

  _loadData(tab, names, locatorNames, order, number);

  // Additional fields

  if (flagAddSampleRank) _createRank(0);

  if (flagAddCoordinates) _createGridCoordinates(flagAddSampleRank);

  // Create the names (for the remaining variables)

  _defineDefaultNames(number, names);

  // Create the locators

  if (flagAddCoordinates)
  {
    Id jcol = 0;
    if (flagAddSampleRank) jcol++;
    setLocatorsByUID(ndim, jcol, ELoc::X, 0);
    _defineDefaultLocators(number, locatorNames);
  }
  initThread();
  return 0;
}

/**
 * Creating a Grid Db which covers the extension of the input 'Db'
 *
 * @param db       Input Db from which the newly created Db is constructed
 * @param nx       Vector of the expected number of grid nodes (default = 10)
 * @param dx       Vector of the expected sizes for the grid meshes (in distance)
 * @param x0       Vector of the expected origin of the grid (in coordinate)
 * @param margin   Vector of the expected margins of the grid (in distance)
 *
 * @remarks Arguments 'nodes' and 'dcell' are disjunctive. If both defined, 'dcell' prevails
 */
Id DbGrid::resetCoveringDb(const Db* db,
                           const VectorInt& nx,
                           const VectorDouble& dx,
                           const VectorDouble& x0,
                           const VectorDouble& margin)
{
  _clear();
  Id ndim = db->getNDim();

  // Derive the Grid parameters

  VectorInt nx_new(ndim);
  VectorDouble x0_new(ndim);
  VectorDouble dx_new(ndim);
  Id nech = 1;
  for (Id idim = 0; idim < ndim; idim++)
  {
    VectorDouble coor = db->getExtrema(idim, true);

    double marge = 0.;
    if (ndim == static_cast<Id>(margin.size())) marge = margin[idim];

    double x0loc = coor[0];
    if (ndim == static_cast<Id>(x0.size())) x0loc = x0[idim];
    x0loc -= marge;

    double ext = coor[1] - x0loc + marge;

    // Constraints specified by the number of nodes
    Id nxloc = 10;
    if (ndim == static_cast<Id>(nx.size()))
      nxloc = nx[idim];
    double dxloc = ext / (static_cast<double>(nxloc) - 1.);

    // Constraints specified by the cell sizes
    if (ndim == static_cast<Id>(dx.size()))
    {
      dxloc = dx[idim];
      nxloc = ceil((ext - dxloc / 2.) / dxloc) + 1;
    }

    nx_new[idim] = nxloc;
    dx_new[idim] = dxloc;
    x0_new[idim] = x0loc;
    nech *= nxloc;
  }

  // Create the grid

  if (gridDefine(nx_new, dx_new, x0_new)) return 1;
  resetDims(ndim, nech);

  /// Load the data

  _createGridCoordinates(0);

  // Create the locators

  Id jcol = 0;
  setLocatorsByUID(ndim, jcol, ELoc::X, 0);

  return 0;
}

/**
 * Creating a regular grid Db which covers the input Polygon
 *
 * @param polygon    Pointer to the input Polygon
 * @param nodes      Vector of the expected number of nodes
 * @param dcell      Vector of the expected dimensions for the grid cells
 * @param flagAddSampleRank true if the sample rank must be generated
 */
Id DbGrid::resetFromPolygon(Polygons* polygon,
                            const VectorInt& nodes,
                            const VectorDouble& dcell,
                            bool flagAddSampleRank)
{
  _clear();
  double xmin, xmax, ymin, ymax;
  Id ndim = 2;

  polygon->getExtension(&xmin, &xmax, &ymin, &ymax);

  // Derive the Grid parameters

  VectorInt nx_tab;
  VectorDouble x0_tab;
  VectorDouble dx_tab;
  Id nech = 1;
  for (Id idim = 0; idim < ndim; idim++)
  {
    double x0  = (idim == 0) ? xmin : ymin;
    double ext = (idim == 0) ? xmax - xmin : ymax - ymin;

    Id nx     = 10;
    double dx = ext / static_cast<double>(nx);
    if (ndim == static_cast<Id>(nodes.size()))
    {
      nx = nodes[idim];
      dx = ext / static_cast<double>(nx);
    }
    if (ndim == static_cast<Id>(dcell.size()))
    {
      dx = dcell[idim];
      nx = static_cast<Id>(ext / dx);
    }

    nx_tab.push_back(nx);
    x0_tab.push_back(x0);
    dx_tab.push_back(dx);
    nech *= nx;
  }
  Id ncol = (flagAddSampleRank) ? ndim + 1 : ndim;

  // Create the grid

  if (gridDefine(nx_tab, dx_tab, x0_tab)) return 1;
  resetDims(ncol, nech);

  /// Load the data

  if (flagAddSampleRank) _createRank(0);
  _createGridCoordinates(flagAddSampleRank);

  // Create the locators

  Id jcol = 0;
  if (flagAddSampleRank) jcol++;
  setLocatorsByUID(ndim, jcol, ELoc::X, 0);

  return 0;
}

DbGrid* DbGrid::create(const VectorInt& nx,
                       const VectorDouble& dx,
                       const VectorDouble& x0,
                       const VectorDouble& angles,
                       const ELoadBy& order,
                       const VectorDouble& tab,
                       const VectorString& names,
                       const VectorString& locatorNames,
                       bool flagAddSampleRank,
                       bool flagAddCoordinates)
{
  auto* dbgrid = new DbGrid;
  if (dbgrid->reset(nx, dx, x0, angles, order, tab, names, locatorNames,
                    flagAddSampleRank, flagAddCoordinates))
  {
    messerr("Error when creating DbGrid from Grid");
    delete dbgrid;
    return nullptr;
  }
  return dbgrid;
}

DbGrid* DbGrid::createCoveringDb(const Db* db,
                                 const VectorInt& nx,
                                 const VectorDouble& dx,
                                 const VectorDouble& x0,
                                 const VectorDouble& margin)
{
  DbGrid* dbgrid = new DbGrid;
  if (dbgrid->resetCoveringDb(db, nx, dx, x0, margin))
  {
    messerr("Error when creating DbGrid covering another Db");
    delete dbgrid;
    return nullptr;
  }
  return dbgrid;
}

DbGrid* DbGrid::createFromPolygon(Polygons* polygon,
                                  const VectorInt& nodes,
                                  const VectorDouble& dcell,
                                  bool flagAddSampleRank)
{
  DbGrid* dbgrid = new DbGrid;
  if (dbgrid->resetFromPolygon(polygon, nodes, dcell, flagAddSampleRank))
  {
    messerr("Error when creating DbGrid from Polygon");
    delete dbgrid;
    return nullptr;
  }
  return dbgrid;
}

DbGrid* DbGrid::coarsify(const VectorInt& nmult)
{
  return createCoarse(this, nmult, 1);
}

DbGrid* DbGrid::createCoarse(DbGrid* dbin,
                             const VectorInt& nmult,
                             bool flagCell,
                             bool flagAddSampleRank)
{
  DbGrid* dbgrid;
  Id ndim = dbin->getNDim();

  // Get the new grid characteristics
  VectorInt nx(ndim);
  VectorDouble dx(ndim);
  VectorDouble x0(ndim);
  dbin->getGrid().multiple(nmult, flagCell, nx, dx, x0);

  // Create the new grid
  dbgrid = create(nx, dx, x0, dbin->getAngles(), ELoadBy::SAMPLE,
                  VectorDouble(), VectorString(), VectorString(), flagAddSampleRank);

  // Migrate all variables (except 'rank' and coordinates
  (void)dbgrid->migrateAllVariables(dbin, true, true, false, flagAddSampleRank);

  return dbgrid;
}

/**
 * Create a new Grid, starting from an initial Grid, and extending its space dimensions
 * A set of Top and Bottom variables is provided which serve in designing the
 * Top and Bottom of the new coordinates.
 * @param gridIn Initial
 * @param tops   Vector of Variable names which define the Tops
 * @param bots   Vector of Variable names which define the Bottoms
 * @param nxnew  Vector giving the number of meshes for each additional space dimension
 * @param verbose Verbose flag
 * @param eps    Each new coordinate is calculated from the top to bottom extension
 *               inflated by eps
 * @return
 */
DbGrid* DbGrid::createFromGridExtend(const DbGrid& gridIn,
                                     const VectorString& tops,
                                     const VectorString& bots,
                                     const VectorInt& nxnew,
                                     bool verbose,
                                     double eps)
{
  DbGrid* gridnew = new DbGrid;

  Id ncoor = static_cast<Id>(nxnew.size());
  if (ncoor <= 0)
  {
    messerr("You must provide a non-empty vector of meshing dimensions");
    return gridnew;
  }
  if (ncoor != static_cast<Id>(tops.size()))
  {
    messerr("Arguments 'tops' and 'nxnew' should have the same dimension");
    return gridnew;
  }
  if (ncoor != static_cast<Id>(bots.size()))
  {
    messerr("Arguments 'bots' and 'nxnew' should have the same dimension");
    return gridnew;
  }

  // Calculate extremes on new coordinates variables

  VectorDouble mini(ncoor);
  VectorDouble maxi(ncoor);
  double coteB, coteT;
  for (Id icoor = 0; icoor < ncoor; icoor++)
  {
    coteB = gridIn.getMinimum(bots[icoor]);
    coteT = gridIn.getMinimum(tops[icoor]);
    if (FFFF(coteB) || FFFF(coteT))
    {
      messerr("The grid extension along variable (%d) is not possible", icoor + 1);
      messerr("The variable has no valid value available or all values are equal");
      return gridnew;
    }
    mini[icoor] = MIN(coteB, coteT);

    coteB = gridIn.getMaximum(bots[icoor]);
    coteT = gridIn.getMaximum(tops[icoor]);
    if (FFFF(coteB) || FFFF(coteT))
    {
      messerr("The grid extension along variable (%d) is not possible", icoor + 1);
      messerr("The variable has no valid value available or all values are equal");
      return gridnew;
    }
    maxi[icoor] = MAX(coteB, coteT);

    if (maxi[icoor] <= mini[icoor])
    {
      messerr("The grid extension along variable (%d) is not possible", icoor + 1);
      messerr("The variable has no valid value available or all values are equal");
      return gridnew;
    }
    if (nxnew[icoor] < 2)
    {
      messerr("The number of meshes along new direction5%d) should be larger than 1", icoor + 1);
      return gridnew;
    }

    if (verbose)
      message("Additional coordinate %d: Minimum = %lf - Maximum = %lf - Nstep = %d\n",
              icoor + 1, mini[icoor], maxi[icoor], nxnew[icoor]);
  }

  // Get the characteristics of Input Grid
  auto ndim           = gridIn.getNDim();
  VectorInt nx        = gridIn.getNXs();
  VectorDouble x0     = gridIn.getX0s();
  VectorDouble dx     = gridIn.getDXs();
  VectorDouble angles = gridIn.getAngles();

  // Extend the characteristics for the new file dimension
  Id ndimnew = ndim + ncoor;
  nx.resize(ndimnew);
  dx.resize(ndimnew);
  x0.resize(ndimnew);
  angles.resize(ndimnew);

  for (Id icoor = 0; icoor < ncoor; icoor++)
  {
    double delta         = maxi[icoor] - mini[icoor];
    nx[ndim + icoor]     = nxnew[icoor];
    x0[ndim + icoor]     = mini[icoor] - delta * eps / 2.;
    dx[ndim + icoor]     = delta * (1. + eps) / nxnew[icoor];
    angles[ndim + icoor] = 0.;
  }

  // Creating the new grid
  gridnew = create(nx, dx, x0, angles);

  return gridnew;
}

/**
 * Create a new grid, from an Initial Grid, by suppressing a set of space dimensions
 * @param gridIn       Initial grid
 * @param deletedRanks Vector of indices of space dimensions to be suppressed
 * @return
 */
DbGrid* DbGrid::createFromGridShrink(const DbGrid& gridIn,
                                     const VectorInt& deletedRanks)
{
  auto* gridnew = new DbGrid();
  auto ndim     = gridIn.getNDim();

  for (Id i = 0; i < static_cast<Id>(deletedRanks.size()); i++)
  {
    if (i < 0 || i >= ndim)
    {
      messerr("The dimension to be removed (%d) should lie within [0,%d[",
              i + 1, ndim);
      return gridnew;
    }
  }
  VectorInt ranks = deletedRanks;
  auto last       = std::unique(ranks.begin(), ranks.end());
  ranks.erase(last, ranks.end());
  std::sort(ranks.begin(), ranks.end());
  std::reverse(ranks.begin(), ranks.end());

  // Get the characteristics of Input Grid
  VectorInt nx        = gridIn.getNXs();
  VectorDouble x0     = gridIn.getX0s();
  VectorDouble dx     = gridIn.getDXs();
  VectorDouble angles = gridIn.getAngles();

  // Suppress the dimensions of the grid
  for (Id i = 0; i < static_cast<Id>(ranks.size()); i++)
  {
    nx.erase(nx.begin() + ranks[i]);
    dx.erase(dx.begin() + ranks[i]);
    x0.erase(x0.begin() + ranks[i]);
    angles.erase(angles.begin() + ranks[i]);
  }

  // Creating the new grid
  gridnew = create(nx, dx, x0, angles);

  return gridnew;
}

VectorInt DbGrid::getNXsExt(Id ndimMax) const
{
  VectorInt nxs = getNXs();
  nxs.resize(ndimMax, 1);
  return nxs;
}

DbGrid* DbGrid::refine(const VectorInt& nmult)
{
  return createRefine(this, nmult, 0);
}

DbGrid* DbGrid::createRefine(DbGrid* dbin,
                             const VectorInt& nmult,
                             bool flagCell,
                             bool flagAddSampleRank)
{
  DbGrid* dbgrid;
  Id ndim = dbin->getNDim();

  // Get the new grid characteristics
  VectorInt nx(ndim);
  VectorDouble dx(ndim);
  VectorDouble x0(ndim);
  dbin->getGrid().divider(nmult, flagCell, nx, dx, x0);

  // Create the new grid
  dbgrid = create(nx, dx, x0, dbin->getAngles(), ELoadBy::SAMPLE,
                  VectorDouble(), VectorString(), VectorString(), flagAddSampleRank);

  // Migrate all variables (except 'rank'  and coordinates
  (void)dbgrid->migrateAllVariables(dbin, true, true, false, flagAddSampleRank);

  return dbgrid;
}

/**
 * Migrate all the variables (Z_locator) from 'dbin' on the nodes of 'this'
(grid)
 * @param dbin  Input Db
 * @param flagAddSampleRank true if the rank of the samples must be aaded
 * @param flag_fill  Filling option
 * @param flag_inter Interpolation
 * @param flag_ball  Use BallTree sorting algorithm when available
 * @return
 */
bool DbGrid::migrateAllVariables(Db* dbin, bool flag_fill, bool flag_inter, bool flag_ball, bool flagAddSampleRank)
{
  ELoc locatorType;
  Id locatorIndex;

  // Constitute the list of Variables to be migrated

  VectorInt icols;
  for (Id icol = 0; icol < dbin->getNColumn(); icol++)
  {
    // Skip the rank
    if (flagAddSampleRank && icol == 0) continue;

    // Skip the coordinates
    String name = dbin->getNameByColIdx(icol);
    if (dbin->getLocatorByColIdx(icol, &locatorType, &locatorIndex))
    {
      if (locatorType == ELoc::X) continue;
    }
    icols.push_back(icol);
  }
  Id ncol = static_cast<Id>(icols.size());
  if (ncol <= 0) return true;

  // Migrate the variables
  auto icolOut = getNColumn();
  if (migrateByAttribute(dbin, this, icols, 2, VectorDouble(), flag_fill,
                         flag_inter, flag_ball, NamingConvention(String())))
    return false;

  // Duplicate the locators
  for (Id icol = 0; icol < ncol; icol++)
  {
    if (dbin->getLocatorByColIdx(icols[icol], &locatorType, &locatorIndex))
      setLocatorByColIdx(icolOut + icol, locatorType, locatorIndex);
    else
      setLocatorByColIdx(icolOut + icol, ELoc::UNKNOWN, 0);
  }
  return true;
}

/**
 * Paint the ndim columns starting from 'icol0' with grid coordinates
 * @param icol0 Starting column
 */
void DbGrid::_createGridCoordinates(Id icol0)
{
  auto ndim = getNDim();

  // Set the Names

  for (Id idim = 0; idim < ndim; idim++)
    _setNameByColIdx(icol0 + idim, getLocatorName(ELoc::X, idim));

  // Set the locators

  setLocatorsByUID(getNDim(), icol0, ELoc::X, 0);

  // Generate the vector of coordinates

  std::vector<double> coors(ndim);
  std::vector<Id> indices;
  _grid.iteratorInit();
  for (Id iech = 0; iech < getNSample(); iech++)
  {
    _grid.iteratorNext(indices);
    _grid.indicesToCoordinateInPlace(indices, coors);
    for (Id idim = 0; idim < ndim; idim++)
      setArray(iech, icol0 + idim, coors[idim]);
  }
}

bool DbGrid::isSameGrid(const Grid& grid) const
{
  if (grid.empty())
  {
    messerr("Both files should be organized as grids");
    return false;
  }
  return _grid.isSame(grid);
}

void DbGrid::gridCopyParams(Id mode, const Grid& gridaux)
{
  _grid.copyParams(mode, gridaux);
}

bool DbGrid::isSameGridMesh(const DbGrid& dbaux) const
{
  if (!dbaux.isGrid())
  {
    messerr("Both files should be organized as grids");
    return false;
  }
  return _grid.isSameMesh(dbaux.getGrid());
}

bool DbGrid::isSameGridRotation(const DbGrid& dbaux) const
{
  if (!dbaux.isGrid())
  {
    messerr("Both files should be organized as grids");
    return false;
  }
  if (!isGridRotated() && !dbaux.isGridRotated()) return true;
  return _grid.isSameRotation(dbaux.getGrid());
}

bool DbGrid::isGridRotated() const
{
  return (_grid.isRotated());
}

/**
 * Return the coordinate of a sample along one Space Dimension
 * @param iech Rank of the sample
 * @param idim Rank of the Space Dimension
 * @param flag_rotate Use the rotation (only for Grid)
 * @return
 */
double DbGrid::getCoordinate(Id iech, Id idim, bool flag_rotate) const
{
  if (idim >= getNDim()) return TEST;
  return _grid.getCoordinate(iech, idim, flag_rotate);
}

void DbGrid::initThread() const
{
  _grid.initThread();
}
void DbGrid::getCoordinatesInPlace(VectorDouble& coor, Id iech, bool flag_rotate) const
{
  const VectorDouble& vec = _grid.getCoordinatesByRank(iech, flag_rotate);
  std::copy(vec.begin(), vec.begin() + getNDim(), coor.begin());
}

Id DbGrid::getNDim() const
{
  return (_grid.getNDim());
}

/**
 * Set dimension
 * @param ncol Number of columns (= variables)
 * @param nech Number of samples (ignore in case of Grid)
 */
void DbGrid::resetDims(Id ncol, Id /*nech*/)
{
  Id nech = _grid.getNTotal();
  Db::resetDims(ncol, nech);
}

bool DbGrid::isConsistent() const
{
  return _grid.getNTotal() == getNSample();
}

bool DbGrid::_deserializeAscii(std::istream& is, bool verbose)
{
  bool ret = true;

  ret = ret && _grid._deserializeAscii(is, verbose);

  ret = ret && Db::_deserializeAscii(is, verbose);

  return ret;
}

bool DbGrid::_serializeAscii(std::ostream& os, bool verbose) const
{
  bool ret = true;

  /* Writing the grid characteristics */

  ret = ret && _grid._serializeAscii(os, verbose);

  /* Writing the tail of the file */

  ret = ret && Db::_serializeAscii(os, verbose);

  return ret;
}

#ifdef HDF5
bool DbGrid::_deserializeH5(H5::Group& grp, bool verbose)
{
  auto dbgridG = SerializeHDF5::getGroup(grp, "DbGrid");
  if (!dbgridG)
  {
    return false;
  }

  bool ret = true;

  // call _deserializeAscii on each member with the current class Group
  ret = ret && _grid._deserializeH5(*dbgridG, verbose);

  // call _deserializeAscii on the parent class with the current class Group
  ret = ret && Db::_deserializeH5(*dbgridG, verbose);

  return ret;
}

bool DbGrid::_serializeH5(H5::Group& grp, bool verbose) const
{
  auto dbgridG = grp.createGroup("DbGrid");

  bool ret = true;

  // call _serializeAscii on each member with the current class Group
  ret = ret && _grid._serializeH5(dbgridG, verbose);

  // call _serializeAscii on the parent class with the current class Group
  ret = ret && Db::_serializeH5(dbgridG, verbose);

  return ret;
}
#endif

double DbGrid::getUnit(Id idim) const
{
  return _grid.getDX(idim);
}

Id DbGrid::gridDefine(const VectorInt& nx,
                      const VectorDouble& dx,
                      const VectorDouble& x0,
                      const VectorDouble& angles)
{
  return _grid.resetFromVector(nx, dx, x0, angles);
}

/**
 * Create a Db by loading the contents of a Neutral File
 *
 * @param NFFilename Name of the Neutral File (Db format)
 * @param verbose    Verbose
 *
 * @remarks The name does not need to be completed in particular when defined by absolute path
 * @remarks or read from the Data Directory (in the gstlearn distribution)
 */
DbGrid* DbGrid::createFromNF(const String& NFFilename, bool verbose)
{
  auto* dbgrid = new DbGrid;
  if (dbgrid->_fileOpenAndDeserialize(NFFilename, verbose)) return dbgrid;
  delete dbgrid;
  return nullptr;
}

VectorDouble DbGrid::getColumnSubGrid(const String& name,
                                      Id idim0,
                                      Id rank,
                                      bool useSel)
{
  VectorDouble vec;
  if (!isGrid())
  {
    messerr("This method is only available for Grid Db");
    return vec;
  }

  // Define optional selection

  VectorDouble sel;
  if (useSel) sel = getSelections();

  // Loop on the samples

  _grid.iteratorInit();
  for (Id iech = 0; iech < getNSample(); iech++)
  {
    VectorInt indices = _grid.iteratorNext();
    if (indices[idim0] != rank) continue;
    Id iabs = _grid.indiceToRank(indices);

    double value = getValue(name, iabs);
    if (useSel && !sel.empty() && sel[iech] == 0) value = TEST;
    vec.push_back(value);
  }
  return vec;
}

void DbGrid::getGridPileInPlace(Id iuid,
                                const VectorInt& indg,
                                Id idim0,
                                VectorDouble& vec) const
{
  auto nz = getNX(idim0);
  if (nz != static_cast<Id>(vec.size())) vec.resize(nz);

  // Loop on the samples

  VectorInt indices = indg;
  VectorInt iechs(nz);
  for (Id iz = 0; iz < nz; iz++)
  {
    indices[idim0] = iz;
    iechs[iz]      = _grid.indiceToRank(indices);
  }
  getArrayVec(iechs, iuid, vec);
}

void DbGrid::setGridPileInPlace(Id iuid,
                                const VectorInt& indg,
                                Id idim0,
                                const VectorDouble& vec)
{
  auto nz = getNX(idim0);
  if (static_cast<Id>(vec.size()) != nz) return;

  // Loop on the samples

  VectorInt indices = indg;
  VectorInt iechs(nz);
  for (Id iz = 0; iz < nz; iz++)
  {
    indices[idim0] = iz;
    iechs[iz]      = _grid.indiceToRank(indices);
  }
  setArrayVec(iechs, iuid, vec);
}

void DbGrid::generateCoordinates(const String& radix)
{
  if (!isGrid())
  {
    messerr("This method is only available in the case of Grid. Nothing done");
    return;
  }
  auto ndim = getNDim();
  VectorDouble coors(ndim);
  (void)addColumnsByConstant(ndim, 0., radix, ELoc::X);
  for (Id iech = 0; iech < getNSample(); iech++)
  {
    _grid.rankToCoordinatesInPlace(iech, coors);
    for (Id idim = 0; idim < ndim; idim++)
      setCoordinate(iech, idim, coors[idim]);
  }
}

/**
 * Returns the contents of one slice extracted from a DbGrid
 * @param name Name of the targte variable
 * @param posx Rank of the first extracted coordinate (in [0, ndim[)
 * @param posy Rank of the second extracted coordinate (in [0, ndim[)
 * @param corner  Vector giving a reference node that belongs to the extracted section
 * @param useSel Use of the current Selection
 * @return
 *
 * @remark The argument 'corner' gives the indices of a node that belongs to the
 * @remarks extracted section. Obviously corner[posx] and corner[posy] are not used
 */
VectorDouble DbGrid::getOneSlice(const String& name,
                                 Id posx,
                                 Id posy,
                                 const VectorInt& corner,
                                 bool useSel) const
{
  VectorDouble tab;
  auto ndim = getNDim();
  if (getNDim() < 2)
  {
    messerr("This method is limited to Grid with space dimension >= 2");
    return tab;
  }
  if (posx < 0 || posx >= ndim)
  {
    messerr("Argument 'posx'(%d) should lie in [0,%d[", posx, ndim);
    return tab;
  }
  if (posy < 0 || posy >= ndim)
  {
    messerr("Argument 'posy'(%d) should lie in [0,%d[", posy, ndim);
    return tab;
  }
  if (posx == posy)
  {
    messerr("Arguments 'posx' and 'posy' should not be similar");
    return tab;
  }
  VectorInt cornloc = corner;
  if (cornloc.empty())
    cornloc.resize(ndim, 0);
  if (ndim != static_cast<Id>(cornloc.size()))
  {
    messerr("The dimension of 'corner' should be equal to 'ndim'");
    return tab;
  }
  auto iuid = getUID(name);
  if (iuid < 0)
  {
    messerr("The Variable %s is not found", name.c_str());
    return tab;
  }

  auto n1 = getNX(posx);
  auto n2 = getNX(posy);
  tab.resize(n1 * n2, TEST);

  VectorInt indices = cornloc;

  Id ecr = 0;
  for (Id i2 = 0; i2 < n2; i2++)
    for (Id i1 = 0; i1 < n1; i1++, ecr++)
    {
      indices[posx] = i1;
      indices[posy] = i2;
      Id iech       = indiceToRank(indices);
      if (!useSel || isActive(iech))
        tab[ecr] = getArray(iech, iuid);
      else
        tab[ecr] = TEST;
    }
  return tab;
}

/**
 * Returns the contents of one slice extracted from a DbGrid
 * @param idim Rank of the target coordinate
 * @param posx Rank of the first extracted coordinate (in [0, ndim[)
 * @param posy Rank of the second extracted coordinate (in [0, ndim[)
 * @param corner  Vector giving a reference node that belongs to the extracted section
 * @param useSel Use of the current Selection
 * @return
 *
 * @remark If idim does not match the Space dimension of the DbGrid, empty vector if returned
 * @remark If the variable exists physically, this variable is read
 * @remark Otherwise, the coordinate is generated on the fly
 *
 * @remark The argument 'corner' gives the indices of a node that belongs to the
 * @remark extracted section. Obviously corner[posx] and corner[posy] are not used
 *
 */
VectorDouble DbGrid::getOneSliceForCoordinate(Id idim,
                                              Id posx,
                                              Id posy,
                                              const VectorInt& corner,
                                              bool useSel) const
{
  VectorDouble tab;
  auto ndim = getNDim();
  if (getNDim() < 2)
  {
    messerr("This method is limited to Grid with space dimension >= 2");
    return tab;
  }
  if (posx < 0 || posx >= ndim)
  {
    messerr("Argument 'posx'(%d) should lie in [0,%d[", posx, ndim);
    return tab;
  }
  if (posy < 0 || posy >= ndim)
  {
    messerr("Argument 'posy'(%d) should lie in [0,%d[", posy, ndim);
    return tab;
  }
  if (posx == posy)
  {
    messerr("Arguments 'posx' and 'posy' should not be similar");
    return tab;
  }
  VectorInt cornloc = corner;
  if (cornloc.empty())
    cornloc.resize(ndim, 0);
  if (ndim != static_cast<Id>(cornloc.size()))
  {
    messerr("The dimension of 'corner' should be equal to 'ndim'");
    return tab;
  }
  if (idim < 0 || idim >= ndim)
  {
    messerr("Argument 'idim'(%d) should lie in [0,%d[", idim, ndim);
    return tab;
  }
  // Check if the variable name already exists
  if (getNLoc(ELoc::X) > 0)
  {
    String name = getNameByLocator(ELoc::X, idim);
    return getOneSlice(name, posx, posy, corner, useSel);
  }
  // The variable does not exist, it must be generated on the fly
  auto n1 = getNX(posx);
  auto n2 = getNX(posy);
  tab.resize(n1 * n2, TEST);

  VectorInt indices = cornloc;
  VectorDouble coord(ndim);

  Id ecr = 0;
  for (Id i2 = 0; i2 < n2; i2++)
    for (Id i1 = 0; i1 < n1; i1++, ecr++)
    {
      indices[posx] = i1;
      indices[posy] = i2;

      indicesToCoordinateInPlace(indices, coord);
      Id iech = indiceToRank(indices);
      if (!useSel || isActive(iech))
        tab[ecr] = coord[idim];
      else
        tab[ecr] = TEST;
    }
  return tab;
}

/**
 * Set all elements of a column (1-D) along a given space dimension
 * to a constant value
 * @param name   Name of the target variable
 * @param idim   Rank of the Space dimension
 * @param rank   Rank of the target Column
 * @param value  Assigned value
 * @param useSel Use the selection
 */
Id DbGrid::assignGridColumn(const String& name,
                            Id idim,
                            Id rank,
                            double value,
                            bool useSel)
{
  if (idim < 0 || idim >= getNDim())
  {
    messerr("Argument 'idim'(%d) is incompatible with Grid dimension(%d)",
            idim, getNDim());
    return 1;
  }
  if (rank < 0 || rank >= getNX(idim))
  {
    messerr("Argument 'rank'(%d) is incompatible with number of cells(%d)",
            rank, getNX(idim));
    return 1;
  }

  _grid.iteratorInit();
  for (Id iech = 0; iech < getNSample(); iech++)
  {
    VectorInt indices = _grid.iteratorNext();
    if (indices[idim] != rank) continue;
    if (useSel && !isActive(iech)) continue;
    setValue(name, iech, value);
  }
  return 0;
}

Id DbGrid::coordinateToRank(const VectorDouble& coor,
                            bool centered,
                            double eps) const
{
  return _grid.coordinateToRank(coor, centered, eps);
}

VectorInt DbGrid::coordinateToIndices(const VectorDouble& coor,
                                      bool centered,
                                      double eps) const
{
  return _grid.coordinateToIndices(coor, centered, eps);
}

Id DbGrid::coordinateToIndicesInPlace(const VectorDouble& coor,
                                      VectorInt& indices,
                                      bool centered,
                                      double eps) const
{
  return _grid.coordinateToIndicesInPlace(coor, indices, centered, eps);
}

Id DbGrid::centerCoordinateInPlace(VectorDouble& coor, bool centered, bool stopIfOut, double eps) const
{
  Id ndim = static_cast<Id>(coor.size());
  VectorInt indice(ndim);
  Id err = coordinateToIndicesInPlace(coor, indice, centered, eps);
  if (stopIfOut && err > 0) return -1;
  indicesToCoordinateInPlace(indice, coor);
  return 0;
}

/**
 * Extracts a slice from a 3-D Grid
 * @param name   Name of the target variable
 * @param pos    Type of section: 0 for YoZ; 1 for XoZ and 2 for XoY
 * @param indice Rank of the section
 * @param useSel Use the active selection
 * @return A VectorVectorDouble with 4 columns, i.e: X, Y, Z, Var
 *
 * @remark In presence of a selection and if useSel is TRUE,
 * @remarks values are returned but set to TEST
 */
VectorVectorDouble DbGrid::getSlice(const String& name,
                                    Id pos,
                                    Id indice,
                                    bool useSel) const
{
  VectorVectorDouble tab;
  Id nvect = 4;
  if (getNDim() != 3)
  {
    messerr("This method is limited to 3-D Grid data base");
    return tab;
  }
  if (!checkArg("Argument 'pos'", pos, 3)) return tab;
  auto iuid = getUID(name);
  if (iuid < 0)
  {
    messerr("The Variable %s is not found", name.c_str());
    return tab;
  }

  tab.resize(nvect);
  VectorInt indices(3);
  VectorDouble coor(3);

  if (pos == 0)
  {
    // Section YoZ
    auto n1 = getNX(1);
    auto n2 = getNX(2);
    auto n3 = getNX(0);
    Id nech = n1 * n2;
    for (Id i = 0; i < nvect; i++) tab[i].resize(nech, TEST);
    if (!checkArg("Error in argument 'indice'", indice, n3)) return VectorVectorDouble();
    indices[0] = indice;

    Id ecr = 0;
    for (Id i1 = 0; i1 < n1; i1++)
      for (Id i2 = 0; i2 < n2; i2++, ecr++)
      {
        indices[1] = i1;
        indices[2] = i2;
        Id iech    = indiceToRank(indices);
        getCoordinatesInPlace(coor, iech);
        tab[0][ecr] = coor[0];
        tab[1][ecr] = coor[1];
        tab[2][ecr] = coor[2];
        if (!useSel || isActive(iech))
          tab[3][ecr] = getArray(iech, iuid);
        else
          tab[3][ecr] = TEST;
      }
  }
  else if (pos == 1)
  {
    // Section XoZ
    auto n1 = getNX(0);
    auto n2 = getNX(2);
    auto n3 = getNX(1);
    Id nech = n1 * n2;
    for (Id i = 0; i < nvect; i++)
      tab[i].resize(nech, TEST);
    if (!checkArg("Error in argument 'indice'", indice, n3)) return VectorVectorDouble();
    indices[1] = indice;

    Id ecr = 0;
    for (Id i1 = 0; i1 < n1; i1++)
      for (Id i2 = 0; i2 < n2; i2++, ecr++)
      {
        indices[0] = i1;
        indices[2] = i2;
        Id iech    = indiceToRank(indices);
        getCoordinatesInPlace(coor, iech);
        tab[0][ecr] = coor[0];
        tab[1][ecr] = coor[1];
        tab[2][ecr] = coor[2];
        if (!useSel || isActive(iech))
          tab[3][ecr] = getArray(iech, iuid);
        else
          tab[3][ecr] = TEST;
      }
  }
  else
  {
    // Section XoY
    auto n1 = getNX(0);
    auto n2 = getNX(1);
    auto n3 = getNX(2);
    Id nech = n1 * n2;
    for (Id i = 0; i < nvect; i++) tab[i].resize(nech, TEST);
    if (!checkArg("Error in argument 'indice'", indice, n3)) return VectorVectorDouble();
    indices[2] = indice;

    Id ecr = 0;
    for (Id i1 = 0; i1 < n1; i1++)
      for (Id i2 = 0; i2 < n2; i2++, ecr++)
      {
        indices[0] = i1;
        indices[1] = i2;
        Id iech    = indiceToRank(indices);
        getCoordinatesInPlace(coor, iech);
        tab[0][ecr] = coor[0];
        tab[1][ecr] = coor[1];
        tab[2][ecr] = coor[2];
        if (!useSel || isActive(iech))
          tab[3][ecr] = getArray(iech, iuid);
        else
          tab[3][ecr] = TEST;
      }
  }
  return tab;
}

/**
 * Return the VectorVectorDouble containing the borders of a cell
 * @param node Target cell
 * @param forceGridMesh When TRUE, returns the edges of the standard grid mesh
 *                      even if a variable block extension is defined
 * @return
 */
VectorVectorDouble DbGrid::getCellEdges(Id node, bool forceGridMesh) const
{
  VectorVectorDouble coords(2);
  coords[0].resize(5);
  coords[1].resize(5);

  auto ndim = getNDim();
  VectorInt icorner(ndim, 0);
  VectorDouble local;

  // Get the extension of the target cell (possibly variable)
  VectorDouble dxsPerCell;
  if (forceGridMesh)
    dxsPerCell = getDXs();
  else
    dxsPerCell = getBlockExtensions(node);

  icorner[0]   = -1;
  icorner[1]   = -1;
  local        = getGrid().getCellCoordinatesByCorner(node, icorner, dxsPerCell);
  coords[0][0] = local[0];
  coords[1][0] = local[1];

  icorner[0]   = -1;
  icorner[1]   = 1;
  local        = getGrid().getCellCoordinatesByCorner(node, icorner, dxsPerCell);
  coords[0][1] = local[0];
  coords[1][1] = local[1];

  icorner[0]   = 1;
  icorner[1]   = 1;
  local        = getGrid().getCellCoordinatesByCorner(node, icorner, dxsPerCell);
  coords[0][2] = local[0];
  coords[1][2] = local[1];

  icorner[0]   = 1;
  icorner[1]   = -1;
  local        = getGrid().getCellCoordinatesByCorner(node, icorner, dxsPerCell);
  coords[0][3] = local[0];
  coords[1][3] = local[1];

  icorner[0]   = -1;
  icorner[1]   = -1;
  local        = getGrid().getCellCoordinatesByCorner(node, icorner, dxsPerCell);
  coords[0][4] = local[0];
  coords[1][4] = local[1];

  return coords;
}

VectorVectorDouble DbGrid::getAllCellsEdges(bool forceGridMesh) const
{
  VectorVectorDouble coords(2);
  auto ndim = getNDim();
  VectorInt icorner(ndim, 0);
  VectorDouble local;

  // Get the extension of the target cell (possibly variable)
  VectorDouble dxsPerCell;
  if (forceGridMesh) dxsPerCell = getDXs();

  for (Id node = 0; node < getNTotal(); node++)
  {
    if (!forceGridMesh) dxsPerCell = getBlockExtensions(node);

    icorner[0] = -1;
    icorner[1] = -1;
    local      = getGrid().getCellCoordinatesByCorner(node, icorner, dxsPerCell);
    coords[0].push_back(local[0]);
    coords[1].push_back(local[1]);

    icorner[0] = -1;
    icorner[1] = 1;
    local      = getGrid().getCellCoordinatesByCorner(node, icorner, dxsPerCell);
    coords[0].push_back(local[0]);
    coords[1].push_back(local[1]);

    icorner[0] = 1;
    icorner[1] = 1;
    local      = getGrid().getCellCoordinatesByCorner(node, icorner, dxsPerCell);
    coords[0].push_back(local[0]);
    coords[1].push_back(local[1]);

    icorner[0] = 1;
    icorner[1] = -1;
    local      = getGrid().getCellCoordinatesByCorner(node, icorner, dxsPerCell);
    coords[0].push_back(local[0]);
    coords[1].push_back(local[1]);
  }
  return coords;
}

/**
 * Return the VectorVectorDouble containing the borders of the grid
 * @return
 */
VectorVectorDouble DbGrid::getGridEdges() const
{
  VectorVectorDouble coords(2);
  coords[0].resize(5);
  coords[1].resize(5);

  auto ndim = getNDim();
  VectorInt icorner(ndim, 0);
  VectorDouble local;

  icorner[0]   = 0;
  icorner[1]   = 0;
  local        = getGrid().getCoordinatesByCorner(icorner);
  coords[0][0] = local[0];
  coords[1][0] = local[1];

  icorner[0]   = 0;
  icorner[1]   = 1;
  local        = getGrid().getCoordinatesByCorner(icorner);
  coords[0][1] = local[0];
  coords[1][1] = local[1];

  icorner[0]   = 1;
  icorner[1]   = 1;
  local        = getGrid().getCoordinatesByCorner(icorner);
  coords[0][2] = local[0];
  coords[1][2] = local[1];

  icorner[0]   = 1;
  icorner[1]   = 0;
  local        = getGrid().getCoordinatesByCorner(icorner);
  coords[0][3] = local[0];
  coords[1][3] = local[1];

  icorner[0]   = 0;
  icorner[1]   = 0;
  local        = getGrid().getCoordinatesByCorner(icorner);
  coords[0][4] = local[0];
  coords[1][4] = local[1];

  return coords;
}

VectorDouble DbGrid::getCodir(const VectorInt& grincr) const
{
  VectorDouble codir = getGrid().indicesToCoordinate(grincr);
  VH::subtractInPlace(codir, getGrid().getX0s());
  VH::normalize(codir);
  return codir;
}

VectorDouble DbGrid::getBlockExtensions(Id node) const
{
  VectorDouble dxsPerCell;
  if (!hasLocVariable(ELoc::BLEX))
    dxsPerCell = getDXs();
  else
  {
    auto ndim  = getNDim();
    dxsPerCell = getLocVariables(ELoc::BLEX, node, ndim);
  }
  return dxsPerCell;
}

Id DbGrid::morpho(const EMorpho& oper,
                  double vmin,
                  double vmax,
                  Id option,
                  const VectorInt& radius,
                  bool flagDistErode,
                  bool verbose,
                  const NamingConvention& namconv)
{
  return dbMorpho(this, oper, vmin, vmax, option, radius, flagDistErode, verbose, namconv);
}

Id DbGrid::smooth(ANeigh* neigh,
                  Id type,
                  double range,
                  const NamingConvention& namconv)
{
  return dbSmoother(this, neigh, type, range, namconv);
}

/****************************************************************************/
/*!
 **  Create a 2-D Db structure
 **
 ** \return  Pointer to the newly created 2-D Db grid structure
 **
 ** \param[in]  order     Manner in which values in tab are ordered
 **                       (ELoadBy)
 ** \param[in]  flagAddSampleRank true to add 'rank' as a supplementary field
 **
 ** \param[in]  nx        Number of grid nodes along X
 ** \param[in]  ny        Number of grid nodes along Y
 ** \param[in]  x0        Grid origin along X
 ** \param[in]  y0        Grid origin along Y
 ** \param[in]  dx        Grid mesh along X
 ** \param[in]  dy        Grid mesh along Y
 ** \param[in]  angle     Rotation angle
 ** \param[in]  tab       Array containing the data
 **
 *****************************************************************************/
DbGrid* DbGrid::createGrid2D(const ELoadBy& order,
                             Id nx,
                             Id ny,
                             double x0,
                             double y0,
                             double dx,
                             double dy,
                             double angle,
                             bool flagAddSampleRank,
                             const VectorDouble& tab)
{
  VectorInt nn(2);
  VectorDouble xx(2);
  VectorDouble dd(2);
  VectorDouble angles(2);

  nn[0]     = nx;
  nn[1]     = ny;
  dd[0]     = dx;
  dd[1]     = dy;
  xx[0]     = x0;
  xx[1]     = y0;
  angles[0] = angle;
  angles[1] = 0.;

  DbGrid* db = DbGrid::create(nn, dd, xx, angles, order, tab, VectorString(),
                              VectorString(), flagAddSampleRank);

  return db;
}

/**
 * Creating a new Db loaded with random values
 * @param nx Vector of mesh indices
 * @param nvar Number of variables
 * @param nfex Number of external drift functions
 * @param ncode Number of codes (no code when 0)
 * @param varmax Maximum value for the measurement error
 * @param selRatio Percentage of samples that must be masked off (between 0 and 1)
 * @param heteroRatio Vector of proportions of NA to be generated per variable
 * @param means Vector of means per variable (optional)
 * @param x0 Vector of coordinates of the origin of the grid (optional)
 * @param seed Value for the Random Generator seed
 * @return A pointer to the newly created DbGrid
 *
 * @remarks
 * The generated grid lies within a [0,1] hypercube.
 * The variance of measurement error is created only if 'varmax' is
 * positive. Then a field is created for each variable. this field is filled
 * with random values uniformly generated in [0, varmax] The external drift
 * values are generated according to Normal distribution.
 */
DbGrid* DbGrid::createFillRandom(const VectorInt& nx,
                                 Id nvar,
                                 Id nfex,
                                 Id ncode,
                                 double varmax,
                                 double selRatio,
                                 const VectorDouble& heteroRatio,
                                 const VectorDouble& means,
                                 const VectorDouble& x0,
                                 Id seed)
{
  // Set the seed
  law_set_random_seed(seed);

  // Create the Db
  Id ndim = static_cast<Id>(nx.size());
  VectorDouble dx(ndim);
  for (Id idim = 0; idim < ndim; idim++) dx[idim] = 1. / nx[idim];
  DbGrid* dbgrid = DbGrid::create(nx, dx, x0);
  Id ndat        = VH::product(nx);

  // Generate the Vectors of Variance of measurement error (optional)
  if (varmax > 0.)
  {
    VectorVectorDouble varm(nvar);
    for (Id ivar = 0; ivar < nvar; ivar++)
      varm[ivar] = VH::simulateUniform(ndat, 0., varmax);
    dbgrid->addColumnsByVVD(varm, "v", ELoc::V);
  }

  // Generate the External Drift functions (optional)
  if (nfex > 0)
  {
    VectorVectorDouble fex(nfex);
    for (Id ifex = 0; ifex < nfex; ifex++) fex[ifex] = VH::simulateGaussian(ndat);
    dbgrid->addColumnsByVVD(fex, "f", ELoc::F);
  }

  // Generate the selection (optional)
  if (selRatio > 0)
  {
    VectorDouble sel(ndat);
    VectorDouble rnd = VH::simulateUniform(ndat);
    for (Id idat = 0; idat < ndat; idat++) sel[idat] = (rnd[idat] > selRatio) ? 1. : 0.;
    dbgrid->addColumns(sel, "sel", ELoc::SEL);
  }

  // Generate the variables
  bool flag_hetero = (static_cast<Id>(heteroRatio.size()) == nvar);
  VectorVectorDouble vars(nvar);
  for (Id ivar = 0; ivar < nvar; ivar++)
  {
    double mean = (means.empty()) ? 0. : means[ivar];
    vars[ivar]  = VH::simulateGaussian(ndat, mean);
    if (flag_hetero)
    {
      VectorDouble rnd = VH::simulateUniform(ndat);
      for (Id idat = 0; idat < ndat; idat++)
        if (rnd[idat] <= heteroRatio[ivar]) vars[ivar][idat] = TEST;
    }
  }
  dbgrid->addColumnsByVVD(vars, "z", ELoc::Z);

  // Generate the code (optional)
  if (ncode > 0)
  {
    VectorDouble codes = VH::simulateUniform(ndat);
    for (Id idat = 0; idat < ndat; idat++) codes[idat] = floor(ncode * codes[idat]);
    dbgrid->addColumns(codes, "code", ELoc::C);
  }

  return dbgrid;
}

/**
 * @brief Create an empty DbGrid structure from the characteristics of a Grid
 *
 * @param grid  Grid structure containing the grid architecture definition
 * @return DbGrid* A Pointer to the newly created data base
 */
DbGrid* DbGrid::createFromGrid(const Grid& grid)
{
  DbGrid* dbgrid = DbGrid::create(grid.getNXs(), grid.getDXs(), grid.getX0s());
  dbgrid->setRotation(grid.getRotation());
  return dbgrid;
}

void DbGrid::_interpolate(const DbGrid* grid3D,
                          Id idim0,
                          double top,
                          double bot,
                          const VectorDouble& vecin,
                          VectorDouble& vecout) const
{
  Id nzin      = grid3D->getNX(idim0);
  double z0out = getX0(idim0);
  double dzout = getDX(idim0);
  auto nzout   = getNX(idim0);

  // Blank out the output vector
  vecout.fill(TEST);

  // Get the top and bottom indices in the output vector
  Id indtop = ceil((top - z0out) / dzout);
  Id indbot = floor((bot - z0out) / dzout);

  for (Id iz = indbot; iz <= indtop; iz++)
  {
    if (iz < 0 || iz >= nzout) continue;
    double zz = z0out + iz * dzout;

    // Find the index in the input vector
    Id izin = static_cast<Id>(static_cast<double>(nzin) * (zz - bot) / (top - bot));
    if (izin < 0 || izin <= nzin) continue;

    // Assign the value
    vecout[iz] = vecin[izin];
  }
}

/**
 * Create the sub-grid, extracted from 'gridIn' and reduced to the vector of limits
 * @param gridIn Input grid
 * @param limits A vector of Min and Max per space dimension (Dimension: [ndim][2])
 * @param flagAddCoordinates True if the grid coordinates must be included in the output file
 * @return
 */
DbGrid* DbGrid::createSubGrid(const DbGrid* gridIn, VectorVectorInt limits, bool flagAddCoordinates)
{
  DbGrid* gridOut = nullptr;
  if (gridIn == nullptr) return gridOut;
  Id ndim = gridIn->getNDim();

  // Preliminary checks
  if (ndim != static_cast<Id>(limits.size()))
  {
    messerr("The argument 'limits' should have dimension ndim x 2");
    return gridOut;
  }

  // Get the list of variables to be copied (rank and coordinates excluded)
  VectorString names = gridIn->getAllNames(true);
  VectorInt iuidIn   = gridIn->getUIDs(names);
  Id nvar            = static_cast<Id>(names.size());

  // Create the characteristics of the new grid
  VectorInt NXs       = gridIn->getNXs();
  VectorDouble DXs    = gridIn->getDXs();
  VectorDouble X0s    = gridIn->getX0s();
  VectorDouble angles = gridIn->getAngles();

  for (Id idim = 0; idim < ndim; idim++)
  {
    NXs[idim] = limits[idim][1] - limits[idim][0];
    X0s[idim] += limits[idim][0] * DXs[idim];
  }

  // Create the new grid
  gridOut = DbGrid::create(NXs, DXs, X0s, angles, ELoadBy::fromKey("SAMPLE"),
                           VectorDouble(), VectorString(), VectorString(), 1,
                           flagAddCoordinates);

  // Add the variables of interest
  VectorInt iuidOut(nvar);
  for (Id ivar = 0; ivar < nvar; ivar++)
    iuidOut[ivar] = gridOut->addColumnsByConstant(1, TEST, names[ivar]);

  // Loop on the nodes of the output sub-grid
  VectorInt indg(ndim);
  double value;
  Id igin;
  for (Id igout = 0, nout = gridOut->getNSample(); igout < nout; igout++)
  {
    // Get the indices in the output grid
    gridOut->rankToIndice(igout, indg);

    // Convert them to the indices in the input grid
    for (Id idim = 0; idim < ndim; idim++)
      indg[idim] += limits[idim][0];

    // Convert in rank in the input grid
    igin = gridIn->indiceToRank(indg);

    // Loop on the variables
    for (Id ivar = 0; ivar < nvar; ivar++)
    {
      value = gridIn->getArray(igin, iuidIn[ivar]);
      gridOut->setArray(igout, iuidOut[ivar], value);
    }
  }
  return gridOut;
}

/**
 * Creating a 3D grid by squeeze-and-stretch forwards, i.e. from structural to sugar box, from:
 * - a 3D grid containing the relevant information
 * - a 2D grid containing the top and bottom information
 * @param grid3Din 3D input grid
 * @param surf2D 2D grid of surfaces
 * @param nameTop Name of the variable in 'surf2D' containing the top information
 * @param nameBot Name of the variable in 'surf2D' containing the bottom information
 * @param names   Vector of names in 'grid3D' to be exported in output 3D grid
 * @param nzout   Number of Vertical meshes in the output 3D grid
 * @param thickmin The algorithm is not applied if
 * @return The output 3D grid (or nullptr in case of error)
 *
 * @remarks:
 * - the grid files 'surf2D' and 'grid3Din' should match (in 2-D)
 * - the grid 'surf2D' contains the top and bottom (identified by the corresponding locators)
 * - the grid 'surf2D' contains a selection which designates the only pixels where
 *   'top' and 'bot' are defined and correctly ordered (o save time)
 */
DbGrid* DbGrid::createSqueezeAndStretchForward(const DbGrid* grid3Din,
                                               const DbGrid* surf2D,
                                               const String& nameTop,
                                               const String& nameBot,
                                               const VectorString& names,
                                               Id nzout,
                                               double thickmin)
{
  DbGrid* grid3Dout = nullptr;

  // Preliminary checks

  if (surf2D == nullptr) return grid3Dout;

  if (surf2D->getNDim() != 2)
  {
    messerr("The grid 'surf2D' must be defined in the 2-D space");
    return grid3Dout;
  }
  if (grid3Din->getNDim() != 3)
  {
    messerr("The grid 'grid3Din' must be defined in the 3-D space");
    return grid3Dout;
  }
  if (!grid3Din->isSameGrid(surf2D->getGrid()))
  {
    messerr("The grid files 'grid3Din' and 'surf2D' should match (in 2D)");
    return grid3Dout;
  }
  if (nzout <= 0.)
  {
    messerr("The number of vertical grid meshes 'nzout' must be strictly positive");
    return grid3Dout;
  }
  if (names.empty())
  {
    messerr("You must designate variable(s) to be copied from input to output 3D grid");
    return grid3Dout;
  }
  Id ndim  = 3;
  Id idim0 = ndim - 1;
  Id nvar  = static_cast<Id>(names.size());

  // Getting relevant information from the top and bottom surfaces (using the selection)
  VectorDouble botArray = surf2D->getColumn(nameBot, true);
  double botmin         = VH::minimum(botArray);
  VectorDouble topArray = surf2D->getColumn(nameTop, true);
  double topmax         = VH::maximum(topArray);
  if (topmax <= botmin)
  {
    messerr("The thickness of the target Layer seems too small for a Squeeze-and-Stretch");
    return grid3Dout;
  }

  // Retrieve information from the 3D input grid
  VectorDouble X0s    = grid3Din->getX0s();
  VectorDouble DXs    = grid3Din->getDXs();
  VectorInt NXs       = grid3Din->getNXs();
  VectorDouble angles = grid3Din->getAngles();

  // Modify these characteristics for the output 3D Grid

  Id nzin   = NXs[idim0];
  double z0 = X0s[idim0];
  double dz = DXs[idim0];

  NXs[idim0] = nzout;
  DXs[idim0] = (topmax - botmin) / static_cast<double>(nzout);
  X0s[idim0] = 0.;

  // Create the output 3D grid (with no coordinate)
  grid3Dout = create(NXs, DXs, X0s, angles, ELoadBy::fromKey("SAMPLE"),
                     VectorDouble(), VectorString(), VectorString(), 1, false);

  // Create the variables in the 3D grid and identify their UIDs
  for (Id ivar = 0; ivar < nvar; ivar++)
    grid3Dout->addColumnsByConstant(1, TEST, names[ivar]);
  VectorInt iuids = grid3Dout->getUIDs(names);

  // Define local variables
  Id iuidTop = surf2D->getUID(nameTop);
  Id iuidBot = surf2D->getUID(nameBot);
  VectorDouble vecin(nzin);
  VectorDouble vecout(nzout);
  VectorInt indg(ndim, 0);
  Id ig2D;
  double top;
  double bot;
  double thick;

  // Loop on the 3-D vertical columns of the 3-D grid
  for (Id ix = 0, nx = grid3Dout->getNX(0); ix < nx; ix++)
    for (Id iy = 0, ny = grid3Dout->getNX(1); iy < ny; iy++)
    {
      indg[0] = ix;
      indg[1] = iy;

      // Identify the corresponding node in 'surf2D'
      ig2D = surf2D->indiceToRank(indg);

      if (!surf2D->isActive(ig2D)) continue;

      // Read the Top and Bottom
      top   = surf2D->getArray(ig2D, iuidTop);
      bot   = surf2D->getArray(ig2D, iuidBot);
      thick = top - bot;
      if (thick < thickmin) continue;

      // Loop on the variables to be transformed
      for (Id ivar = 0; ivar < nvar; ivar++)
      {
        // Read the pile from the 3D input grid
        grid3Din->getGridPileInPlace(iuids[ivar], indg, idim0, vecin);

        // Perform the squeeze-and-stretch forward operation
        VH::squeezeAndStretchInPlaceForward(vecin, vecout, z0, dz, top, bot);

        // Write the resulting pile in the output 3D grid
        grid3Dout->setGridPileInPlace(iuids[ivar], indg, idim0, vecout);
      }
    }
  return grid3Dout;
}

/**
 * Creating a 3D grid by squeeze-and-stretch backwards, i.e. from sugar box to structural, from:
 * - a 3D grid containing the relevant information
 * - a 2D grid containing the top and bottom information
 * @param grid3Din 3D input grid
 * @param surf2D 2D grid of surfaces
 * @param nameTop Name of the variable in 'surf2D' containing the top information
 * @param nameBot Name of the variable in 'surf2D' containing the bottom information
 * @param names   Vector of names in 'grid3D' to be exported in output 3D grid
 * @param nzout, z0out, dzout Specification along third dimension of the output 3D Grid
 * @return The output 3D grid (or nullptr in case of error)
 *
 * @remarks:
 * - the grid files 'surf2D' and 'grid3Din' should match (in 2-D)
 * - the grid 'surf2D' contains the top and bottom (identified by the corresponding locators)
 * - the grid 'surf2D' contains a selection which designates the only pixels where
 *   'top' and 'bot' are defined and correctly ordered (o save time)
 */
DbGrid* DbGrid::createSqueezeAndStretchBackward(const DbGrid* grid3Din,
                                                const DbGrid* surf2D,
                                                const String& nameTop,
                                                const String& nameBot,
                                                const VectorString& names,
                                                Id nzout,
                                                double z0out,
                                                double dzout)
{
  DbGrid* grid3Dout = nullptr;

  // Preliminary checks

  if (surf2D == nullptr) return grid3Dout;

  if (surf2D->getNDim() != 2)
  {
    messerr("The grid 'surf2D' must be defined in the 2-D space");
    return grid3Dout;
  }
  if (grid3Din->getNDim() != 3)
  {
    messerr("The grid 'grid3Din' must be defined in the 3-D space");
    return grid3Dout;
  }
  if (!grid3Din->isSameGrid(surf2D->getGrid()))
  {
    messerr("The grid files 'grid3Din' and 'surf2D' should match (in 2D)");
    return grid3Dout;
  }
  if (names.empty())
  {
    messerr("You must designate variable(s) to be copied from input to output 3D grid");
    return grid3Dout;
  }
  Id ndim  = 3;
  Id idim0 = ndim - 1;
  Id nvar  = static_cast<Id>(names.size());

  // Getting relevant information from the top and bottom surfaces (using the selection)
  VectorDouble botArray = surf2D->getColumn(nameBot, true);
  double botmin         = VH::minimum(botArray);
  VectorDouble topArray = surf2D->getColumn(nameTop, true);
  double topmax         = VH::maximum(topArray);
  if (topmax <= botmin)
  {
    messerr("The thickness of the target Layer seems too small for a Squeeze-and-Stretch");
    return grid3Dout;
  }

  // Retrieve information from the 3D input grid
  VectorDouble X0s    = grid3Din->getX0s();
  VectorDouble DXs    = grid3Din->getDXs();
  VectorInt NXs       = grid3Din->getNXs();
  VectorDouble angles = grid3Din->getAngles();

  // Modify these characteristics for the output 3D Grid
  Id nzin    = NXs[idim0];
  NXs[idim0] = nzout;
  DXs[idim0] = dzout;
  X0s[idim0] = z0out;

  // Create the output 3D grid
  grid3Dout = create(NXs, DXs, X0s, angles, ELoadBy::fromKey("SAMPLE"),
                     VectorDouble(), VectorString(), VectorString(), 1, false);

  // Create the variables in the 3D grid and identify their UIDs
  for (Id ivar = 0; ivar < nvar; ivar++)
    grid3Dout->addColumnsByConstant(1, TEST, names[ivar]);
  VectorInt iuids = grid3Dout->getUIDs(names);

  // Define local variables
  Id iuidTop = surf2D->getUID(nameTop);
  Id iuidBot = surf2D->getUID(nameBot);
  VectorDouble vecin(nzin);
  VectorDouble vecout(nzout);
  VectorInt indg(ndim, 0);
  Id ig2D;
  double top;
  double bot;

  // Loop on the 3-D vertical columns of the 3-D grid
  for (Id ix = 0, nx = grid3Dout->getNX(0); ix < nx; ix++)
    for (Id iy = 0, ny = grid3Dout->getNX(1); iy < ny; iy++)
    {
      indg[0] = ix;
      indg[1] = iy;

      // Identify the corresponding node in 'surf2D'
      ig2D = surf2D->indiceToRank(indg);

      if (!surf2D->isActive(ig2D)) continue;

      // Read the Top and Bottom
      top = surf2D->getArray(ig2D, iuidTop);
      bot = surf2D->getArray(ig2D, iuidBot);

      // Loop on the variables to be transformed
      for (Id ivar = 0; ivar < nvar; ivar++)
      {
        // Read the pile from the 3D input grid
        grid3Din->getGridPileInPlace(iuids[ivar], indg, idim0, vecin);

        // Perform the squeeze-and-stretch operation backwards
        VH::squeezeAndStretchInPlaceBackward(vecin, vecout, z0out, dzout, top, bot);

        // Write the resulting pile in the output 3D grid
        grid3Dout->setGridPileInPlace(iuids[ivar], indg, idim0, vecout);
      }
    }
  return grid3Dout;
}

/**
 * Returns the minimum and maximum indices of the subgrid
 * where variables 'nameTop' and 'nameBot' are both defined
 * @param nameTop Name of the Top variable
 * @param nameBot Name of the Bottom variable
 * @param dimExclude Array giving excluding dimension (see details)
 * @return A vector of Min and Max per space dimension (Dimension: [ndim][2])
 *
 * @details: When a dimension is 'excluded', the reduction of the output grid
 * should not be applied to this dimension
 */
VectorVectorInt DbGrid::getLimitsFromVariableExtend(const String& nameTop,
                                                    const String& nameBot,
                                                    const VectorInt& dimExclude) const
{
  auto ndim = getNDim();
  VectorVectorInt vec(ndim);
  for (Id idim = 0; idim < ndim; idim++)
  {
    vec[idim].resize(2);
    vec[idim][0] = 0;
    vec[idim][1] = getNX(idim);
  }
  if (nameTop.empty() || nameBot.empty()) return vec;

  // Find the set of Min and Max indices of the subgrid

  auto nech = getNSample(true);
  VectorInt indmin(ndim, 10000000);
  VectorInt indmax(ndim, -10000000);
  VectorInt indg(ndim);
  Id iuid_top = getUID(nameTop);
  Id iuid_bot = getUID(nameBot);

  for (Id iech = 0; iech < nech; iech++)
  {
    // Discard not relevant pixels
    if (!isActive(iech)) continue;
    double top = getArray(iech, iuid_top);
    double bot = getArray(iech, iuid_bot);
    if (FFFF(top) || FFFF(bot) || bot > top) continue;

    rankToIndice(iech, indg);
    for (Id idim = 0; idim < ndim; idim++)
    {
      Id indloc = indg[idim];
      if (indloc < indmin[idim]) indmin[idim] = indloc;
      if (indloc > indmax[idim]) indmax[idim] = indloc;
    }
  }

  // Discard the case where the sub-grid does not exist

  bool flag_exist = true;
  for (Id idim = 0; idim < ndim && flag_exist; idim++)
  {
    if (indmin[idim] > indmax[idim]) flag_exist = false;
  }
  if (!flag_exist) return vec;

  // Get the sub-grid characteristics

  vec.resize(ndim);
  for (Id idim = 0; idim < ndim; idim++)
  {
    vec[idim].resize(2);
    vec[idim][0] = indmin[idim];
    vec[idim][1] = indmax[idim];
  }

  for (Id iexc = 0, nexc = static_cast<Id>(dimExclude.size()); iexc < nexc; iexc++)
  {
    vec[iexc][0] = 0;
    vec[iexc][1] = getNX(iexc);
  }
  return vec;
}

/**
 * Defines a selection in the current grid where variables 'nameTop' and 'nameBot'
 * are both defined and ordered properly
 * @param nameTop Name of the Top variable
 * @param nameBot Name of the Bottom variable
 *
 * @details: This method also adds a selection in the current grid
 * which masks off the pixels where 'nameTop' and 'nameBot' are defined
 * but not correctly ordered.
 * This is the reason why this method cannot be 'const'
 */
Id DbGrid::setSelectionFromVariableExtend(const String& nameTop, const String& nameBot)
{
  // Create the selection new variable
  Id iuidSel = addColumnsByConstant(1, 1, "SelLayer", ELoc::SEL);

  if (nameTop.empty() || nameBot.empty()) return -1;

  // Find the set of Min and Max indices of the subgrid

  auto nech   = getNSample(true);
  Id iuid_top = getUID(nameTop);
  Id iuid_bot = getUID(nameBot);

  for (Id iech = 0; iech < nech; iech++)
  {
    // Discard not relevant pixels
    if (!isActive(iech)) continue;
    double top = getArray(iech, iuid_top);
    double bot = getArray(iech, iuid_bot);
    if (FFFF(top) || FFFF(bot) || bot > top)
    {
      setArray(iech, iuidSel, 0);
      continue;
    }
  }
  return iuidSel;
}

/**
 * Clean the contents of a 3D file by using surfaces extracted from the 2D file
 * @param names   Vector of variable of the current grid which must be processed
 * @param surf2D  Name of the 2-D grid file containing the surfaces
 * @param nameTop Name of the Top surface (optional)
 * @param nameBot Name of the Bottom surface (optional)
 * @param verbose Verbose flag
 *
 * @remark The input file 'surf2D' and the current grid should match (in 2-D)
 */
void DbGrid::clean3DFromSurfaces(const VectorString& names,
                                 const DbGrid* surf2D,
                                 const String& nameTop,
                                 const String& nameBot,
                                 bool verbose)
{
  if (surf2D == nullptr) return;
  if (surf2D->getNDim() != 2)
  {
    messerr("The grid 'surf2D' must be defined in the 2-D space");
    return;
  }
  if (getNDim() != 3)
  {
    messerr("The current grid must be defined in the 3-D space");
    return;
  }
  if (!isSameGrid(surf2D->getGrid()))
  {
    messerr("The grid 'surf2D' and the current one should coincide horizontally");
    return;
  }
  if (names.empty())
  {
    messerr("You must define some variable to be processed");
    return;
  }

  bool limitsDefined = !nameTop.empty() && !nameBot.empty();
  Id nvar            = static_cast<Id>(names.size());

  // Loop on the vertical columns of the 3-D grid

  Id ndim    = 3;
  Id idim0   = ndim - 1;
  double top = MAXIMUM_BIG;
  double bot = MINIMUM_BIG;
  double z0  = getX0(idim0);
  double dz  = getDX(idim0);
  auto nz    = getNX(idim0);
  Id indzmin = 0; // included
  Id indzmax;
  VectorInt indg(ndim, 0);
  VectorDouble vec(nz);
  VectorDouble vecempty(nz, TEST);
  VectorInt iuids = getUIDs(names);

  Id nmodif3D   = 0;
  Id nmodif2D   = 0;
  Id rank2D     = 0;
  double thick  = 0;
  double thickA = 0;
  for (Id ix = 0, nx = getNX(0); ix < nx; ix++)
    for (Id iy = 0, ny = getNX(1); iy < ny; iy++)
    {
      indg[0] = ix;
      indg[1] = iy;
      indg[2] = 0;

      // Identify the corresponding node in 'surf2D'
      rank2D = surf2D->indiceToRank(indg);

      // Discard if the coordinate does not belong to 'surf2D' extension
      bool flagRead = rank2D >= 0;

      // Avoid processing the column if masked in the 2-D file
      if (flagRead)
      {
        if (!surf2D->isActive(rank2D)) flagRead = false;
      }

      // Get the Top and bottom information from the surf2D grid
      // Converts them into minimum (included) and maximum (excluded) layer indices
      // Note: No use to test correct surface ordering as it is already captured in the selection
      indzmin = 0;
      indzmax = nz;
      if (flagRead && limitsDefined)
      {
        top = surf2D->getValue(nameTop, rank2D);
        if (FFFF(top)) flagRead = false;

        if (flagRead)
        {
          bot = surf2D->getValue(nameBot, rank2D);
          if (FFFF(bot)) flagRead = false;
        }

        if (flagRead)
        {
          indzmin = floor((bot - z0) / dz);
          if (indzmin < 0) indzmin = 0;
          indzmax = ceil((top - z0) / dz);
          if (indzmax > nz) indzmax = nz;
          thick = dz * (indzmax - indzmin + 1);
          if (thick > thickA) thickA = thick;
        }
      }

      // Loop on the variables
      if (flagRead)
      {
        for (Id ivar = 0; ivar < nvar; ivar++)
        {
          // Partial update
          getGridPileInPlace(iuids[ivar], indg, idim0, vec);

          // Blank out
          for (Id iz = 0; iz < indzmin; iz++, nmodif3D++)
            vec[iz] = TEST;
          for (Id iz = indzmax; iz < nz; iz++, nmodif3D++)
            vec[iz] = TEST;

          setGridPileInPlace(iuids[ivar], indg, idim0, vec);
        }
      }
      else
      {
        // Complete update
        nmodif2D++;
        for (Id ivar = 0; ivar < nvar; ivar++, nmodif3D++)
        {
          setGridPileInPlace(iuids[ivar], indg, idim0, vecempty);
        }
      }
    }

  // Optional printout
  if (verbose)
  {
    message("Blanking out the 3-D grid file:\n");
    message("- Number of nodes              = %d %d %d\n", getNX(0), getNX(1), getNX(2));
    message("- Number of variables          = %d\n", nvar);
    message("- Total number of piles        = %d\n", getNX(0) * getNX(1));
    message("- Number of piles blanked out  = %d\n", nmodif2D);
    message("- Number of values blanked out = %d\n", nmodif3D);
    message("- Maximum Layer thickness      = %lf\n", thickA);
  }
}

VectorInt DbGrid::locateDataInGrid(const Db* data,
                                   const VectorInt& rankIn,
                                   bool centered,
                                   bool useSel) const
{
  VectorInt rankOut;

  if (data == nullptr) return VectorInt();

  if (!rankIn.empty())
  {

    // Locate the samples defined by their ranks stored in 'rankIn'

    for (Id ip = 0; ip < static_cast<Id>(rankIn.size()); ip++)
    {
      VectorDouble coor = data->getSampleCoordinates(rankIn[ip]);
      rankOut.push_back(coordinateToRank(coor, centered));
    }
  }
  else
  {

    // Locate all samples (using useSel criterion)

    for (Id ip = 0, np = data->getNSample(useSel); ip < np; ip++)
    {
      if (data->isActive(ip) || !useSel)
      {
        VectorDouble coor = data->getSampleCoordinates(ip);
        rankOut.push_back(coordinateToRank(coor, centered));
      }
    }
  }
  return rankOut;
}

bool DbGrid::hasSingleBlock() const
{
  for (Id idim = 0; idim < getNDim(); idim++)
    if (getNX(idim) == 1) return true;
  return false;
}

/**
 * Create a selection based on the count of active samples of 'db'
 * @param db  Db used for statistics
 * @param nmin Minimum number of samples
 * @param radius Radius of the cell neighborhood used when counting the samples
 * @param option Type of structuring element: 0 for Cross and 1 for Block
 * @param dilation Vector giving the radius extension for Dilation operation
 * @param verbose Verbose flag
 * @param namconv Naming convention
 * @return
 */
Id DbGrid::addSelectionFromDbByMorpho(Db* db,
                                      Id nmin,
                                      Id radius,
                                      Id option,
                                      const VectorInt& dilation,
                                      bool verbose,
                                      const NamingConvention& namconv)
{
  if (db == nullptr)
  {
    messerr("You must define a valid Db");
    return 1;
  }

  auto nech = getNSample();

  VectorString names = db->getNamesByColIdx({0});
  Id iuid            = addColumnsByConstant(1);
  if (dbStatisticsInGridTool(db, this, names, EStatOption::NUM, radius, iuid)) return 1;
  VectorDouble stats = getColumnByUID(iuid, false, false);
  for (Id iech = 0; iech < nech; iech++)
    stats[iech] = (stats[iech] <= nmin) ? 0. : 1.;
  setColumnByUID(stats, iuid, false);
  setLocatorByUID(iuid, ELoc::Z, 0, true);

  Id err = morpho(EMorpho::DILATION, 0.5, 1.5, option, dilation, false, verbose, namconv);

  deleteColumnByUID(iuid);
  return err;
}

void DbGrid::getSampleAsSTInPlace(Id iech, SpaceTarget& P) const
{
  // Initiate the SpacePoint (performed in Db class)
  Db::getSampleAsSTInPlace(iech, P);

  // Load the extension
  if (P.checkExtend())
    P.setExtend(getBlockExtensions(iech));
}

/**
 * Generate a set of discretization locations, relative to the block center
 * Dimension: number of discretization locations, space dimension
 * @param ndiscs Number of discretization (per space dimension)
 * @param iech   Rank of the target sample (used if flagPerCell = true)
 * @param flagPerCell TRUE when the cell dimension are read from the Db (BLEX)
 * @param flagRandom TRUE if the discretization location must be randomized
 * @param seed Seed for random number generator
 * @return
 *
 * @remark: Although randomization can be performed, this process does not consume
 * random numbers.
 */
VectorVectorDouble DbGrid::getDiscretizedBlock(const VectorInt& ndiscs,
                                               Id iech,
                                               bool flagPerCell,
                                               bool flagRandom,
                                               Id seed) const
{
  auto ndim = getNDim();
  Id ntot   = VH::product(ndiscs);
  auto memo = law_get_random_seed();
  law_set_random_seed(seed);
  VectorVectorDouble discs(ntot);
  for (Id i = 0; i < ntot; i++) discs[i].resize(ndim);

  /* Loop on the discretization points */

  for (Id i = 0; i < ntot; i++)
  {
    Id jech = i;
    Id nval = ntot;
    for (Id idim = ndim - 1; idim >= 0; idim--)
    {
      double taille = (!flagPerCell) ? getDX(idim) : getLocVariable(ELoc::BLEX, iech, idim);
      Id nd         = ndiscs[idim];
      nval /= nd;
      Id j = jech / nval;
      jech -= j * nval;
      double local = taille * ((j + 0.5) / nd - 0.5);
      if (!flagRandom)
        discs[i][idim] = local;
      else
        discs[i][idim] = local + taille * law_uniform(-0.5, 0.5) / static_cast<double>(nd);
    }
  }
  law_set_random_seed(memo);

  return discs;
}

/****************************************************************************/
/*!
 **  Create a Grid Db as a multiple of another Grid Db
 **
 ** \return  Pointer to the newly created Db grid structure
 **
 ** \param[in]  dbin      Initial Db Grid
 ** \param[in]  nmult     Array of multiplicity coefficients
 ** \param[in]  flagAddSampleRank true to add the 'rank' as first column
 **  **
 *****************************************************************************/
DbGrid* DbGrid::createMultiple(DbGrid* dbin,
                               const VectorInt& nmult,
                               bool flagAddSampleRank)
{
  DbGrid* dbout = nullptr;
  if (dbin == nullptr) return (dbin);
  Id ndim = dbin->getNDim();

  /* Core allocation */

  VectorInt nx(ndim);
  VectorDouble dx(ndim);
  VectorDouble x0(ndim);

  /* Get the new grid characteristics */

  dbin->getGrid().multiple(nmult, 1, nx, dx, x0);

  /* Create the new grid */

  dbout = DbGrid::create(nx, dx, x0, dbin->getAngles(), ELoadBy::COLUMN,
                         VectorDouble(), VectorString(), VectorString(),
                         flagAddSampleRank);

  return dbout;
}

/****************************************************************************/
/*!
 **  Create a Grid Db as a divider of another Grid Db
 **
 ** \return  Pointer to the newly created Db grid structure
 **
 ** \param[in]  dbin      Initial Db Grid
 ** \param[in]  nmult     Array of subdivision coefficients
 ** \param[in]  flagAddSampleRank true to add the 'rank' as first column
 **
 *****************************************************************************/
DbGrid* DbGrid::createDivider(DbGrid* dbin,
                              const VectorInt& nmult,
                              bool flagAddSampleRank)
{
  DbGrid* dbout = nullptr;
  if (dbin == nullptr) return dbin;

  Id ndim = dbin->getNDim();
  VectorInt nx(ndim);
  VectorDouble dx(ndim);
  VectorDouble x0(ndim);

  /* Get the new grid characteristics */

  dbin->getGrid().divider(nmult, 1, nx, dx, x0);

  /* Create the new grid */

  dbout = DbGrid::create(nx, dx, x0, dbin->getAngles(), ELoadBy::COLUMN,
                         VectorDouble(), VectorString(), VectorString(),
                         flagAddSampleRank);

  return dbout;
}

VectorDouble DbGrid::getDistanceToOrigin(const VectorInt& origin,
                                         const VectorDouble& radius)
{
  auto ndim           = getNDim();
  auto nech           = getNSample();
  VectorDouble coor0  = getCoordinatesByIndice(origin);
  VectorDouble radloc = radius;
  if (ndim != static_cast<Id>(radloc.size())) radloc = VectorDouble(ndim, 1.);

  VectorDouble coor(ndim);
  VectorDouble distvec(nech, 0.);
  for (Id iech = 0; iech < nech; iech++)
  {
    getCoordinatesInPlace(coor, iech);
    double dist = 0.;
    for (Id idim = 0; idim < ndim; idim++)
    {
      double delta = (coor[idim] - coor0[idim]) / radloc[idim];
      dist += delta * delta;
    }
    distvec[iech] = sqrt(dist);
  }
  return distvec;
}
} // namespace gstlrn