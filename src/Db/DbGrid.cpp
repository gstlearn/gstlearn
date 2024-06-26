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
#include "geoslib_f_private.h"
#include "geoslib_old_f.h"
#include "geoslib_define.h"

#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Polygon/Polygons.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"
#include "Basic/CSVformat.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Limits.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Law.hpp"
#include "Basic/AException.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/GlobalEnvironment.hpp"
#include "Stats/Classical.hpp"
#include "Estimation/CalcImage.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Morpho/Morpho.hpp"
#include "Space/SpaceTarget.hpp"

#include <algorithm>
#include <functional>
#include <math.h>

DbGrid::DbGrid()
    : Db(),
      _grid(0)
{
  _clear();
}

DbGrid::DbGrid(const DbGrid& r)
    : Db(r),
      _grid(r._grid)
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

  const DbStringFormat* dbfmt = dynamic_cast<const DbStringFormat*>(strfmt);
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
 * @param names         Variable names (size = nvar)
 * @param locatorNames  Locators for each variable (size = nvar)
 * @param flagAddSampleRank If true, add an automatic rank variable
 * @param flagAddCoordinates If TRUE, add the grid coordinates
 */
int DbGrid::reset(const VectorInt& nx,
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

  int ndim = static_cast<int> (nx.size());
  int nech = 1;
  for (int idim = 0; idim < ndim; idim++)
    nech *= nx[idim];
  int ntab = (tab.empty()) ? 0 : (int) (tab.size() / nech);
  int number = 0;
  if (flagAddSampleRank) number += 1;
  if (flagAddCoordinates) number += ndim;
  int ncol = number + ntab;

  // Create the grid

  if (gridDefine(nx, dx, x0, angles)) return 1;
  resetDims(ncol, nech);

  // Load the data

  _loadData(tab, names, locatorNames, order, number);

  // Additional fields

  if (flagAddSampleRank) _createRank(0);

  if (flagAddCoordinates) _createCoordinatesGrid(flagAddSampleRank);

  // Create the names (for the remaining variables)

  _defineDefaultNames(number, names);

  // Create the locators

  if (flagAddCoordinates)
  {
    int jcol = 0;
    if (flagAddSampleRank) jcol++;
    setLocatorsByUID(ndim, jcol, ELoc::X);
    _defineDefaultLocators(number, locatorNames);
  }

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
int DbGrid::resetCoveringDb(const Db* db,
                            const VectorInt& nx,
                            const VectorDouble& dx,
                            const VectorDouble& x0,
                            const VectorDouble& margin)
{
  _clear();
  int ndim = db->getNDim();

  // Derive the Grid parameters

  VectorInt    nx_new(ndim);
  VectorDouble x0_new(ndim);
  VectorDouble dx_new(ndim);
  int nech = 1;
  for (int idim = 0; idim < ndim; idim++)
  {
    VectorDouble coor = db->getExtrema(idim, true);

    double marge = 0.;
    if (ndim == (int) margin.size()) marge = margin[idim];

    double x0loc = coor[0];
    if (ndim == (int) x0.size()) x0loc = x0[idim];
    x0loc -= marge;

    double ext = coor[1] - x0loc + marge;

    // Constraints specified by the number of nodes
    int nxloc = 10;
    if (ndim == (int) nx.size())
      nxloc = nx[idim];
    double dxloc = ext / ((double) nxloc - 1.);

    // Constraints specified by the cell sizes
    if (ndim == (int) dx.size())
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
  resetDims(ndim,nech);

  /// Load the data

  _createCoordinatesGrid(0);

  // Create the locators

  int jcol = 0;
  setLocatorsByUID(ndim, jcol, ELoc::X);

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
int DbGrid::resetFromPolygon(Polygons* polygon,
                             const VectorInt& nodes,
                             const VectorDouble& dcell,
                             bool flagAddSampleRank)
{
  _clear();
  double xmin, xmax, ymin, ymax;
  int ndim = 2;

  polygon->getExtension(&xmin, &xmax, &ymin, &ymax);

  // Derive the Grid parameters

  VectorInt    nx_tab;
  VectorDouble x0_tab;
  VectorDouble dx_tab;
  int nech = 1;
  for (int idim = 0; idim < ndim; idim++)
  {
    double x0  = (idim == 0) ? xmin : ymin;
    double ext = (idim == 0) ? xmax - xmin : ymax - ymin;

    int nx = 10;
    double dx = ext / (double) nx;
    if (ndim == (int) nodes.size())
    {
      nx = nodes[idim];
      dx = ext / (double) nx;
    }
    if (ndim == (int) dcell.size())
    {
      dx = dcell[idim];
      nx = static_cast<int> (ext / dx);
    }

    nx_tab.push_back(nx);
    x0_tab.push_back(x0);
    dx_tab.push_back(dx);
    nech *= nx;
  }
  int ncol = (flagAddSampleRank) ? ndim + 1 : ndim;

  // Create the grid

  if (gridDefine(nx_tab, dx_tab, x0_tab)) return 1;
  resetDims(ncol, nech);

  /// Load the data

  if (flagAddSampleRank) _createRank(0);
  _createCoordinatesGrid(flagAddSampleRank);

  // Create the locators

  int jcol = 0;
  if (flagAddSampleRank) jcol++;
  setLocatorsByUID(ndim, jcol, ELoc::X);

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
  DbGrid* dbgrid = new DbGrid;
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

DbGrid* DbGrid::coarsify(const VectorInt &nmult)
{
  return createCoarse(this,nmult,1);
}

DbGrid* DbGrid::createCoarse(DbGrid *dbin,
                             const VectorInt &nmult,
                             int flag_cell,
                             bool flagAddSampleRank)
{
  DbGrid *dbgrid;
  int ndim = dbin->getNDim();

  // Get the new grid characteristics
  VectorInt nx(ndim);
  VectorDouble dx(ndim);
  VectorDouble x0(ndim);
  dbin->getGrid().multiple(nmult, flag_cell, nx, dx, x0);

  // Create the new grid
  dbgrid = create(nx, dx, x0, dbin->getAngles(), ELoadBy::SAMPLE,
                  VectorDouble(), VectorString(), VectorString(), flagAddSampleRank);

  // Migrate all variables (except 'rank' and coordinates
  (void) migrateAllVariables(dbin, dbgrid, flagAddSampleRank);

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
DbGrid* DbGrid::createFromGridExtend(const DbGrid &gridIn,
                                     const VectorString& tops,
                                     const VectorString& bots,
                                     const VectorInt &nxnew,
                                     bool verbose,
                                     double eps)
{
  DbGrid *gridnew = new DbGrid;

  int ncoor = (int) nxnew.size();
  if (ncoor <= 0)
  {
    messerr("You must provide a non-empty vector of meshing dimensions");
    return gridnew;
  }
  if (ncoor != (int) tops.size())
  {
    messerr("Arguments 'tops' and 'nxnew' should have the same dimension");
    return gridnew;
  }
  if (ncoor != (int) bots.size())
  {
    messerr("Arguments 'bots' and 'nxnew' should have the same dimension");
    return gridnew;
  }

  // Calculate extremes on new coordinates variables

  VectorDouble mini(ncoor);
  VectorDouble maxi(ncoor);
  double coteB, coteT;
  for (int icoor = 0; icoor < ncoor; icoor++)
  {
    coteB = gridIn.getMinimum(bots[icoor]);
    coteT = gridIn.getMinimum(tops[icoor]);
    if (FFFF(coteB) || FFFF(coteT))
    {
      messerr("The grid extension along variable (%d) is not possible",icoor+1);
      messerr("The variable has no valid value available or all values are equal");
      return gridnew;
    }
    mini[icoor] = MIN(coteB, coteT);

    coteB = gridIn.getMaximum(bots[icoor]);
    coteT = gridIn.getMaximum(tops[icoor]);
    if (FFFF(coteB) || FFFF(coteT))
    {
      messerr("The grid extension along variable (%d) is not possible",icoor+1);
      messerr("The variable has no valid value available or all values are equal");
      return gridnew;
    }
    maxi[icoor] = MAX(coteB, coteT);

    if (maxi[icoor] <= mini[icoor])
    {
      messerr("The grid extension along variable (%d) is not possible",icoor+1);
      messerr("The variable has no valid value available or all values are equal");
      return gridnew;
    }
    if (nxnew[icoor] < 2)
    {
      messerr("The number of meshes along new direction5%d) should be larger than 1",icoor+1);
      return gridnew;
    }

    if (verbose)
      message("Additional coordinate %d: Minimum = %lf - Maximum = %lf - Nstep = %d\n",
              icoor+1, mini[icoor], maxi[icoor], nxnew[icoor]);

  }

  // Get the characteristics of Input Grid
  int ndim = gridIn.getNDim();
  VectorInt nx = gridIn.getNXs();
  VectorDouble x0 = gridIn.getX0s();
  VectorDouble dx = gridIn.getDXs();
  VectorDouble angles = gridIn.getAngles();

  // Extend the characteristics for the new file dimension
  int ndimnew = ndim + ncoor;
  nx.resize(ndimnew);
  dx.resize(ndimnew);
  x0.resize(ndimnew);
  angles.resize(ndimnew);

  for (int icoor = 0; icoor < ncoor; icoor++)
  {
    double delta = maxi[icoor] - mini[icoor];
    nx[ndim + icoor] = nxnew[icoor];
    x0[ndim + icoor] = mini[icoor] - delta * eps / 2.;
    dx[ndim + icoor] = delta * (1. + eps) / nxnew[icoor];
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
DbGrid* DbGrid::createFromGridShrink(const DbGrid &gridIn,
                                     const VectorInt& deletedRanks)
{
  DbGrid* gridnew = new DbGrid();
  int ndim = gridIn.getNDim();

  for (int i = 0; i < (int) deletedRanks.size(); i++)
  {
    if (i < 0 || i >= ndim)
    {
      messerr("The dimension to be removed (%d) should lie within [0,%d[",
              i+1, ndim);
      return gridnew;
    }
  }
  VectorInt ranks = deletedRanks;
  (void) std::unique(ranks.begin(), ranks.end());
  std::sort(ranks.begin(), ranks.end());
  std::reverse(ranks.begin(), ranks.end());

  // Get the characteristics of Input Grid
  VectorInt nx = gridIn.getNXs();
  VectorDouble x0 = gridIn.getX0s();
  VectorDouble dx = gridIn.getDXs();
  VectorDouble angles = gridIn.getAngles();

  // Suppress the dimensions of the grid
  for (int i = 0; i < (int) ranks.size(); i++)
  {
    nx.erase(nx.begin()+ranks[i]);
    dx.erase(dx.begin()+ranks[i]);
    x0.erase(x0.begin()+ranks[i]);
    angles.erase(angles.begin()+ranks[i]);
  }

  // Creating the new grid
  gridnew = create(nx, dx, x0, angles);

  return gridnew;
}

VectorInt DbGrid::getNXsExt(int ndimMax) const
{
  VectorInt nxs = getNXs();
  nxs.resize(ndimMax, 1);
  return nxs;
}

DbGrid* DbGrid::refine(const VectorInt &nmult)
{
  return createRefine(this,nmult,0);
}

DbGrid* DbGrid::createRefine(DbGrid *dbin,
                             const VectorInt &nmult,
                             int flag_cell,
                             bool flagAddSampleRank)
{
  DbGrid *dbgrid;
  int ndim = dbin->getNDim();

  // Get the new grid characteristics
  VectorInt nx(ndim);
  VectorDouble dx(ndim);
  VectorDouble x0(ndim);
  dbin->getGrid().divider(nmult, flag_cell, nx, dx, x0);

  // Create the new grid
  dbgrid = create(nx, dx, x0, dbin->getAngles(), ELoadBy::SAMPLE,
                  VectorDouble(), VectorString(), VectorString(), flagAddSampleRank);

  // Migrate all variables (except 'rank'  and coordinates
  (void) migrateAllVariables(dbin, dbgrid, flagAddSampleRank);

  return dbgrid;
}

/**
 * Migrate all the variables (Z_locator) from 'dbin' on the nodes of 'dbout' (grid)
 * @param dbin  Input Db
 * @param dbout Output db
 * @param flagAddSampleRank true if the rank of the samples must be aaded
 * @return
 */
bool DbGrid::migrateAllVariables(Db *dbin, Db *dbout, bool flagAddSampleRank)
{
  ELoc locatorType;
  int  locatorIndex;

  // Constitute the list of Variables to be migrated

  VectorInt icols;
  for (int icol = 0; icol < dbin->getColumnNumber(); icol++)
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
  int ncol = (int) icols.size();
  if (ncol <= 0) return true;

  // Migrate the variables
  int icolOut = dbout->getColumnNumber();
  if (migrateByAttribute(dbin, dbout, icols, 2, VectorDouble(), true, true, false,
                         NamingConvention(String()))) return false;

  // Duplicate the locators
  for (int icol = 0; icol < ncol; icol++)
  {
    if (dbin->getLocatorByColIdx(icols[icol], &locatorType, &locatorIndex))
      dbout->setLocatorByColIdx(icolOut + icol, locatorType, locatorIndex);
    else
      dbout->setLocatorByColIdx(icolOut + icol, ELoc::UNKNOWN, 0);
  }
  return true;
}

/**
 * Paint the ndim columns starting from 'icol0' with grid coordinates
 * @param icol0 Starting column
 */
void DbGrid::_createCoordinatesGrid(int icol0)
{
  int ndim = getNDim();

  // Set the Names

  for (int idim = 0; idim < ndim; idim++)
    _setNameByColIdx(icol0 + idim, getLocatorName(ELoc::X, idim));

  // Set the locators

  setLocatorsByUID(getNDim(), icol0, ELoc::X);

  // Generate the vector of coordinates

  VectorDouble coors(ndim);
  _grid.iteratorInit();
  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    VectorInt indices = _grid.iteratorNext();
    _grid.indicesToCoordinateInPlace(indices, coors);
    for (int idim = 0; idim < ndim; idim++)
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

void DbGrid::gridCopyParams(int mode, const Grid& gridaux)
{
  _grid.copyParams(mode, gridaux);
}

bool DbGrid::isSameGridMesh(const DbGrid& dbaux) const
{
  if (! dbaux.isGrid())
  {
    messerr("Both files should be organized as grids");
    return false;
  }
  return _grid.isSameMesh(dbaux.getGrid());
}

bool DbGrid::isSameGridRotation(const DbGrid& dbaux) const
{
  if (! dbaux.isGrid())
  {
    messerr("Both files should be organized as grids");
    return false;
  }
  if (! isGridRotated() && ! dbaux.isGridRotated()) return true;
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
double DbGrid::getCoordinate(int iech, int idim, bool flag_rotate) const
{
  if (idim >= getNDim()) return TEST;
  return _grid.getCoordinate(iech, idim, flag_rotate);
}

void DbGrid::getCoordinatesPerSampleInPlace(int iech, VectorDouble& coor, bool flag_rotate) const
{
  VectorDouble vec = _grid.getCoordinatesByRank(iech, flag_rotate);
  coor = vec;
}

int DbGrid::getNDim() const
{
  return (_grid.getNDim());
}

/**
 * Set dimension
 * @param ncol Number of columns (= variables)
 * @param nech Number of samples (ignore in case of Grid)
 */
void DbGrid::resetDims(int ncol, int /*nech*/)
{
  int nech = _grid.getNTotal();
  Db::resetDims(ncol, nech);
}

bool DbGrid::_deserialize(std::istream& is, bool verbose)
{
  int ndim = 0;
  VectorInt nx;
  VectorString locators;
  VectorString names;
  VectorDouble x0;
  VectorDouble dx;
  VectorDouble angles;
  VectorDouble values;
  VectorDouble allvalues;

  /* Initializations */

  bool ret = true;
  ret = ret && _recordRead<int>(is, "Space Dimension", ndim);

  /* Core allocation */

  nx.resize(ndim);
  dx.resize(ndim);
  x0.resize(ndim);
  angles.resize(ndim);

  /* Read the grid characteristics */

  for (int idim = 0; ret && idim < ndim; idim++)
  {
    ret = ret && _recordRead<int>(is, "Grid Number of Nodes", nx[idim]);
    ret = ret && _recordRead<double>(is, "Grid Origin", x0[idim]);
    ret = ret && _recordRead<double>(is, "Grid Mesh", dx[idim]);
    ret = ret && _recordRead<double>(is, "Grid Angles", angles[idim]);
  }

  // Create the Grid characteristics
  (void) gridDefine(nx, dx, x0, angles);

  ret && Db::_deserialize(is, verbose);

  return ret;
}

bool DbGrid::_serialize(std::ostream& os, bool verbose) const
{
  bool ret = true;

  /* Writing the header */

  ret = ret && _recordWrite<int>(os, "Space Dimension", getNDim());

  /* Writing the grid characteristics */

  ret = ret && _commentWrite(os, "Grid characteristics (NX,X0,DX,ANGLE)");
  for (int idim = 0; ret && idim < getNDim(); idim++)
  {
    ret = ret && _recordWrite<int>(os, "",  getNX(idim));
    ret = ret && _recordWrite<double>(os, "", getX0(idim));
    ret = ret && _recordWrite<double>(os, "", getDX(idim));
    ret = ret && _recordWrite<double>(os, "", getAngle(idim));
    ret = ret && _commentWrite(os, "");
  }

  /* Writing the tail of the file */

  ret && Db::_serialize(os, verbose);

  return ret;
}

double DbGrid::getUnit(int idim) const
{
  return _grid.getDX(idim);
}

int DbGrid::gridDefine(const VectorInt& nx,
                       const VectorDouble& dx,
                       const VectorDouble& x0,
                       const VectorDouble& angles)
{
  return (_grid.resetFromVector(nx, dx, x0, angles));
}

/**
 * Create a Db by loading the contents of a Neutral File
 *
 * @param neutralFilename Name of the Neutral File (Db format)
 * @param verbose         Verbose
 *
 * @remarks The name does not need to be completed in particular when defined by absolute path
 * @remarks or read from the Data Directory (in the gstlearn distribution)
 */
DbGrid* DbGrid::createFromNF(const String& neutralFilename, bool verbose)
{
  DbGrid* dbgrid = nullptr;
  std::ifstream is;
  dbgrid = new DbGrid;
  bool success = false;
  if (dbgrid->_fileOpenRead(neutralFilename, is, verbose))
  {
    success = dbgrid->deserialize(is, verbose);
  }
  if (! success)
  {
    delete dbgrid;
    dbgrid = nullptr;
  }
  return dbgrid;
}

VectorDouble DbGrid::getColumnSubGrid(const String& name,
                                     int idim0,
                                     int rank,
                                     bool useSel)
{
  VectorDouble vec;
  if (! isGrid())
  {
    messerr("This method is only available for Grid Db");
    return vec;
  }

  // Define optional selection

  VectorDouble sel;
  if (useSel) sel = getSelections();

  // Loop on the samples

  _grid.iteratorInit();
  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    VectorInt indices = _grid.iteratorNext();
    if (indices[idim0] != rank) continue;
    int iabs = _grid.indiceToRank(indices);

    double value = getValue(name, iabs);
    if (useSel && !sel.empty() && sel[iech] == 0) value = TEST;
    vec.push_back(value);
  }
  return vec;
}

void DbGrid::getGridPileInPlace(int iuid,
                                const VectorInt &indg,
                                int idim0,
                                VectorDouble &vec) const
{
  int nz = getNX(idim0);
  if (nz != (int) vec.size()) vec.resize(nz);

  // Loop on the samples

  VectorInt indices = indg;
  VectorInt iechs(nz);
  for (int iz = 0; iz < nz; iz++)
  {
    indices[idim0] = iz;
    iechs[iz] = _grid.indiceToRank(indices);
  }
  getArrayVec(iechs, iuid, vec);
}

void DbGrid::setGridPileInPlace(int iuid,
                                const VectorInt &indg,
                                int idim0,
                                const VectorDouble &vec)
{
  int nz = getNX(idim0);
  if ((int) vec.size() != nz) return;

  // Loop on the samples

  VectorInt indices = indg;
  VectorInt iechs(nz);
  for (int iz = 0; iz < nz; iz++)
  {
    indices[idim0] = iz;
    iechs[iz] = _grid.indiceToRank(indices);
  }
  setArrayVec(iechs, iuid, vec);
}

void DbGrid::generateCoordinates(const String& radix)
{
  if (! isGrid())
  {
    messerr("This method is only available in the case of Grid. Nothing done");
    return;
  }
  int ndim = getNDim();
  VectorDouble coors(ndim);
  (void) addColumnsByConstant(ndim, 0., radix, ELoc::X);
  display();
  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    _grid.rankToCoordinatesInPlace(iech, coors);
    for (int idim = 0; idim < ndim; idim++)
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
                                 int posx,
                                 int posy,
                                 const VectorInt& corner,
                                 bool useSel) const
{
  VectorDouble tab;
  int ndim = getNDim();
  if (getNDim() < 2)
  {
    messerr("This method is limited to Grid with space dimension >= 1");
    return tab;
  }
  if (posx < 0 || posx >= ndim)
  {
    messerr("Argument 'posx'(%d) should lie in [0,%d[",posx,ndim);
    return tab;
  }
  if (posy < 0 || posy >= ndim)
  {
    messerr("Argument 'posy'(%d) should lie in [0,%d[",posy,ndim);
    return tab;
  }
  if (posx == posy)
  {
    messerr("Arguments 'posx' and 'posy' should not be similar");
    return tab;
  }
  VectorInt cornloc = corner;
  if (cornloc.empty())
    cornloc.resize(ndim,0);
  if (ndim != (int) cornloc.size())
  {
    messerr("The dimension of 'corner' should be equal to 'ndim'");
    return tab;
  }
  int iuid = getUID(name);
  if (iuid < 0)
  {
    messerr("The Variable %s is not found",name.c_str());
    return tab;
  }

  int n1 = getNX(posx);
  int n2 = getNX(posy);
  tab.resize(n1 * n2, TEST);

  VectorInt indices = cornloc;

  int ecr = 0;
  for (int i2 = 0; i2 < n2; i2++)
    for (int i1 = 0; i1 < n1; i1++, ecr++)
    {
      indices[posx] = i1;
      indices[posy] = i2;
      int iech = indiceToRank(indices);
      if (! useSel || isActive(iech))
        tab[ecr] = getArray(iech, iuid);
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
int DbGrid::assignGridColumn(const String& name,
                             int idim,
                             int rank,
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
  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    VectorInt indices = _grid.iteratorNext();
    if (indices[idim] != rank) continue;
    if (useSel && ! isActive(iech)) continue;
    setValue(name, iech, value);
  }
  return 0;
}

int DbGrid::coordinateToRank(const VectorDouble &coor,
                             bool centered,
                             double eps) const
{
  return _grid.coordinateToRank(coor,centered,eps);
}

VectorInt DbGrid::coordinateToIndices(const VectorDouble &coor,
                                      bool centered,
                                      double eps) const
{
  return _grid.coordinateToIndices(coor, centered, eps);
}

int DbGrid::coordinateToIndicesInPlace(const VectorDouble &coor,
                                       VectorInt &indices,
                                       bool centered,
                                       double eps) const
{
  return _grid.coordinateToIndicesInPlace(coor, indices, centered, eps);
}

int DbGrid::centerCoordinateInPlace(VectorDouble& coor, bool centered, bool stopIfOut, double eps) const
{
  int ndim = (int) coor.size();
  VectorInt indice(ndim);
  int err = coordinateToIndicesInPlace(coor,indice,centered,eps);
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
                                    int pos,
                                    int indice,
                                    bool useSel) const
{
  VectorVectorDouble tab;
  int nvect = 4;
  if (getNDim() != 3)
  {
    messerr("This method is limited to 3-D Grid data base");
    return tab;
  }
  if (pos < 0 || pos > 2)
  {
    mesArg("Argument 'pos'", pos, 3);
    return tab;
  }
  int iuid = getUID(name);
  if (iuid < 0)
  {
    messerr("The Variable %s is not found",name.c_str());
    return tab;
  }

  tab.resize(nvect);
  VectorInt indices(3);
  VectorDouble coor(3);

  if (pos == 0)
  {
    // Section YoZ
    int n1 = getNX(1);
    int n2 = getNX(2);
    int n3 = getNX(0);
    int nech = n1 * n2;
    for (int i = 0; i < nvect; i++) tab[i].resize(nech,TEST);
    if (indice < 0 || indice >= n3)
    {
      mesArg("Error in argument 'indice'",indice,n3);
      return VectorVectorDouble();
    }
    indices[0] = indice;

    int ecr = 0;
    for (int i1 = 0; i1 < n1; i1++)
      for (int i2 = 0; i2 < n2; i2++, ecr++)
      {
        indices[1] = i1;
        indices[2] = i2;
        int iech = indiceToRank(indices);
        getCoordinatesPerSampleInPlace(iech, coor);
        tab[0][ecr] = coor[0];
        tab[1][ecr] = coor[1];
        tab[2][ecr] = coor[2];
        if (! useSel || isActive(iech))
          tab[3][ecr] = getArray(iech, iuid);
        else
          tab[3][ecr] = TEST;
      }
  }
  else if (pos == 1)
  {
    // Section XoZ
    int n1 = getNX(0);
    int n2 = getNX(2);
    int n3 = getNX(1);
    int nech = n1 * n2;
    for (int i = 0; i < nvect; i++)
      tab[i].resize(nech, TEST);
    if (indice < 0 || indice >= n3)
    {
      mesArg("Error in argument 'indice'",indice,n3);
      return VectorVectorDouble();
    }
    indices[1] = indice;

    int ecr = 0;
    for (int i1 = 0; i1 < n1; i1++)
      for (int i2 = 0; i2 < n2; i2++, ecr++)
      {
        indices[0] = i1;
        indices[2] = i2;
        int iech = indiceToRank(indices);
        getCoordinatesPerSampleInPlace(iech, coor);
        tab[0][ecr] = coor[0];
        tab[1][ecr] = coor[1];
        tab[2][ecr] = coor[2];
        if (! useSel || isActive(iech))
          tab[3][ecr] = getArray(iech, iuid);
        else
          tab[3][ecr] = TEST;
      }
  }
  else
  {
    // Section XoY
    int n1 = getNX(0);
    int n2 = getNX(1);
    int n3 = getNX(2);
    int nech = n1 * n2;
    for (int i = 0; i < nvect; i++)
      tab[i].resize(nech, TEST);
    if (indice < 0 || indice >= n3)
    {
      mesArg("Error in argument 'indice'",indice,n3);
      return VectorVectorDouble();
    }
    indices[2] = indice;

    int ecr = 0;
    for (int i1 = 0; i1 < n1; i1++)
      for (int i2 = 0; i2 < n2; i2++, ecr++)
      {
        indices[0] = i1;
        indices[1] = i2;
        int iech = indiceToRank(indices);
        getCoordinatesPerSampleInPlace(iech, coor);
        tab[0][ecr] = coor[0];
        tab[1][ecr] = coor[1];
        tab[2][ecr] = coor[2];
        if (! useSel || isActive(iech))
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
VectorVectorDouble DbGrid::getCellEdges(int node, bool forceGridMesh) const
{
  VectorVectorDouble coords(2);
  coords[0].resize(5);
  coords[1].resize(5);

  int ndim = getNDim();
  VectorInt icorner(ndim,0);
  VectorDouble local;

  // Get the extension of the target cell (possibly variable)
  VectorDouble dxsPerCell;
  if (forceGridMesh)
    dxsPerCell = getDXs();
  else
    dxsPerCell = getBlockExtensions(node);

  icorner[0] = -1;
  icorner[1] = -1;
  local = getGrid().getCellCoordinatesByCorner(node, icorner, dxsPerCell);
  coords[0][0] = local[0];
  coords[1][0] = local[1];

  icorner[0] = -1;
  icorner[1] =  1;
  local = getGrid().getCellCoordinatesByCorner(node, icorner, dxsPerCell);
  coords[0][1] = local[0];
  coords[1][1] = local[1];

  icorner[0] = 1;
  icorner[1] = 1;
  local = getGrid().getCellCoordinatesByCorner(node, icorner, dxsPerCell);
  coords[0][2] = local[0];
  coords[1][2] = local[1];

  icorner[0] = 1;
  icorner[1] = -1;
  local = getGrid().getCellCoordinatesByCorner(node, icorner, dxsPerCell);
  coords[0][3] = local[0];
  coords[1][3] = local[1];

  icorner[0] = -1;
  icorner[1] = -1;
  local = getGrid().getCellCoordinatesByCorner(node, icorner, dxsPerCell);
  coords[0][4] = local[0];
  coords[1][4] = local[1];

  return coords;
}

VectorVectorDouble DbGrid::getAllCellsEdges(bool forceGridMesh) const
{
  VectorVectorDouble coords(2);
  int ndim = getNDim();
  VectorInt icorner(ndim,0);
  VectorDouble local;

  // Get the extension of the target cell (possibly variable)
  VectorDouble dxsPerCell;
  if (forceGridMesh) dxsPerCell = getDXs();

  for (int node = 0; node < getNTotal(); node++)
  {
    if (! forceGridMesh) dxsPerCell = getBlockExtensions(node);

    icorner[0] = -1;
    icorner[1] = -1;
    local = getGrid().getCellCoordinatesByCorner(node, icorner, dxsPerCell);
    coords[0].push_back(local[0]);
    coords[1].push_back(local[1]);

    icorner[0] = -1;
    icorner[1] = 1;
    local = getGrid().getCellCoordinatesByCorner(node, icorner, dxsPerCell);
    coords[0].push_back(local[0]);
    coords[1].push_back(local[1]);

    icorner[0] = 1;
    icorner[1] = 1;
    local = getGrid().getCellCoordinatesByCorner(node, icorner, dxsPerCell);
    coords[0].push_back(local[0]);
    coords[1].push_back(local[1]);

    icorner[0] = 1;
    icorner[1] = -1;
    local = getGrid().getCellCoordinatesByCorner(node, icorner, dxsPerCell);
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

  int ndim = getNDim();
  VectorInt icorner(ndim,0);
  VectorDouble local;

  icorner[0] = 0;
  icorner[1] = 0;
  local = getGrid().getCoordinatesByCorner(icorner);
  coords[0][0] = local[0];
  coords[1][0] = local[1];

  icorner[0] = 0;
  icorner[1] = 1;
  local = getGrid().getCoordinatesByCorner(icorner);
  coords[0][1] = local[0];
  coords[1][1] = local[1];

  icorner[0] = 1;
  icorner[1] = 1;
  local = getGrid().getCoordinatesByCorner(icorner);
  coords[0][2] = local[0];
  coords[1][2] = local[1];

  icorner[0] = 1;
  icorner[1] = 0;
  local = getGrid().getCoordinatesByCorner(icorner);
  coords[0][3] = local[0];
  coords[1][3] = local[1];

  icorner[0] = 0;
  icorner[1] = 0;
  local = getGrid().getCoordinatesByCorner(icorner);
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

VectorDouble DbGrid::getBlockExtensions(int node) const
{
  int ndim = getNDim();

  VectorDouble dxsPerCell = getDXs();
  if (hasLocVariable(ELoc::BLEX))
  {
    for (int idim = 0; idim < ndim; idim++)
    {
      double value = getLocVariable(ELoc::BLEX,node, idim);
      if (! FFFF(value)) dxsPerCell[idim] = value;
    }
  }
  return dxsPerCell;
}

int DbGrid::morpho(const EMorpho &oper,
                   double vmin,
                   double vmax,
                   int option,
                   const VectorInt &radius,
                   bool flagDistErode,
                   bool verbose,
                   const NamingConvention &namconv)
{
  return dbMorpho(this, oper, vmin, vmax, option, radius, flagDistErode, verbose, namconv);
}

int DbGrid::smooth(ANeigh *neigh,
                   int type,
                   double range,
                   const NamingConvention &namconv)
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
DbGrid* DbGrid::createGrid2D(const ELoadBy &order,
                             int nx,
                             int ny,
                             double x0,
                             double y0,
                             double dx,
                             double dy,
                             double angle,
                             bool flagAddSampleRank,
                             const VectorDouble &tab)
{
  VectorInt nn(2);
  VectorDouble xx(2);
  VectorDouble dd(2);
  VectorDouble angles(2);

  nn[0] = nx;
  nn[1] = ny;
  dd[0] = dx;
  dd[1] = dy;
  xx[0] = x0;
  xx[1] = y0;
  angles[0] = angle;
  angles[1] = 0.;

  DbGrid *db = DbGrid::create(nn, dd, xx, angles, order, tab, VectorString(),
                              VectorString(), flagAddSampleRank);

  return db;
}

void DbGrid::_interpolate(const DbGrid *grid3D,
                          int idim0,
                          double top,
                          double bot,
                          const VectorDouble &vecin,
                          VectorDouble &vecout)
{
  int    nzin  = grid3D->getNX(idim0);
  double z0out = getX0(idim0);
  double dzout = getDX(idim0);
  int    nzout = getNX(idim0);

  // Blank out the output vector
  vecout.fill(TEST);

  // Get the top and bottom indices in the output vector
  int indtop = ceil((top - z0out)  / dzout);
  int indbot = floor((bot - z0out) / dzout);

  for (int iz = indbot; iz <= indtop; iz++)
  {
    if (iz < 0 || iz >= nzout) continue;
    double zz = z0out + iz * dzout;

    // Find the index in the input vector
    int izin = (int) (double(nzin) * (zz - bot) / (top - bot));
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
  int ndim = gridIn->getNDim();

  // Preliminary checks
  if (ndim != (int) limits.size())
  {
    messerr("The argument 'limits' should have dimension ndim x 2");
    return gridOut;
  }

  // Get the list of variables to be copied (rank and coordinates excluded)
  VectorString names = gridIn->getAllNames(true);
  VectorInt iuidIn = gridIn->getUIDs(names);
  int nvar = (int) names.size();

  // Create the characteristics of the new grid
  VectorInt NXs       = gridIn->getNXs();
  VectorDouble DXs    = gridIn->getDXs();
  VectorDouble X0s    = gridIn->getX0s();
  VectorDouble angles = gridIn->getAngles();

  for (int idim = 0; idim < ndim; idim++)
  {
    NXs[idim]  = limits[idim][1] - limits[idim][0];
    X0s[idim] += limits[idim][0] * DXs[idim];
  }

  // Create the new grid
  gridOut = DbGrid::create(NXs, DXs, X0s, angles, ELoadBy::fromKey("SAMPLE"),
                           VectorDouble(), VectorString(), VectorString(), 1,
                           flagAddCoordinates);

  // Add the variables of interest
  VectorInt iuidOut(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
    iuidOut[ivar] = gridOut->addColumnsByConstant(1, TEST, names[ivar]);

  // Loop on the nodes of the output sub-grid
  VectorInt indg(ndim);
  double value;
  int igin;
  for (int igout = 0, nout = gridOut->getSampleNumber(); igout < nout; igout++)
  {
    // Get the indices in the output grid
    gridOut->rankToIndice(igout, indg);

    // Convert them to the indices in the input grid
    for (int idim = 0; idim < ndim; idim++)
      indg[idim] += limits[idim][0];

    // Convert in rank in the input grid
    igin = gridIn->indiceToRank(indg);

    // Loop on the variables
    for (int ivar = 0; ivar < nvar; ivar++)
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
                                               const DbGrid *surf2D,
                                               const String &nameTop,
                                               const String &nameBot,
                                               const VectorString &names,
                                               int nzout,
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
  if (! grid3Din->isSameGrid(surf2D->getGrid()))
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
  int ndim = 3;
  int idim0 = ndim - 1;
  int nvar = (int) names.size();

  // Getting relevant information from the top and bottom surfaces (using the selection)
  VectorDouble botArray = surf2D->getColumn(nameBot, true);
  double botmin = VH::minimum(botArray);
  VectorDouble topArray = surf2D->getColumn(nameTop, true);
  double topmax = VH::maximum(topArray);
  if (topmax <= botmin)
  {
    messerr("The thickness of the target Layer seems too small for a Squeeze-and-Stretch");
    return grid3Dout;
  }

  // Retrieve information from the 3D input grid
  VectorDouble X0s = grid3Din->getX0s();
  VectorDouble DXs = grid3Din->getDXs();
  VectorInt    NXs = grid3Din->getNXs();
  VectorDouble angles = grid3Din->getAngles();

  // Modify these characteristics for the output 3D Grid

  int nzin  = NXs[idim0];
  double z0 = X0s[idim0];
  double dz = DXs[idim0];

  NXs[idim0] = nzout;
  DXs[idim0] = (topmax - botmin) / (double) nzout;
  X0s[idim0] = 0.;

  // Create the output 3D grid (with no coordinate)
  grid3Dout = create(NXs, DXs, X0s, angles, ELoadBy::fromKey("SAMPLE"),
                     VectorDouble(), VectorString(), VectorString(), 1, false);

  // Create the variables in the 3D grid and identify their UIDs
  for (int ivar = 0; ivar < nvar; ivar++)
    grid3Dout->addColumnsByConstant(1, TEST, names[ivar]);
  VectorInt iuids = grid3Dout->getUIDs(names);

  // Define local variables
  int iuidTop = surf2D->getUID(nameTop);
  int iuidBot = surf2D->getUID(nameBot);
  VectorDouble vecin(nzin);
  VectorDouble vecout(nzout);
  VectorInt indg(ndim, 0);
  int ig2D;
  double top;
  double bot;
  double thick;

  // Loop on the 3-D vertical columns of the 3-D grid
  for (int ix = 0, nx = grid3Dout->getNX(0); ix < nx; ix++)
    for (int iy = 0, ny = grid3Dout->getNX(1); iy < ny; iy++)
    {
      indg[0] = ix;
      indg[1] = iy;

      // Identify the corresponding node in 'surf2D'
      ig2D = surf2D->indiceToRank(indg);

      if (! surf2D->isActive(ig2D)) continue;

      // Read the Top and Bottom
      top = surf2D->getArray(ig2D, iuidTop);
      bot = surf2D->getArray(ig2D, iuidBot);
      thick = top - bot;
      if (thick < thickmin) continue;

      // Loop on the variables to be transformed
      for (int ivar = 0; ivar < nvar; ivar++)
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
DbGrid* DbGrid::createSqueezeAndStretchBackward(const DbGrid *grid3Din,
                                                const DbGrid *surf2D,
                                                const String &nameTop,
                                                const String &nameBot,
                                                const VectorString &names,
                                                int nzout,
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
  if (! grid3Din->isSameGrid(surf2D->getGrid()))
  {
    messerr("The grid files 'grid3Din' and 'surf2D' should match (in 2D)");
    return grid3Dout;
  }
  if (names.empty())
  {
    messerr("You must designate variable(s) to be copied from input to output 3D grid");
    return grid3Dout;
  }
  int ndim = 3;
  int idim0 = ndim - 1;
  int nvar = (int) names.size();

  // Getting relevant information from the top and bottom surfaces (using the selection)
  VectorDouble botArray = surf2D->getColumn(nameBot, true);
  double botmin = VH::minimum(botArray);
  VectorDouble topArray = surf2D->getColumn(nameTop, true);
  double topmax = VH::maximum(topArray);
  if (topmax <= botmin)
  {
    messerr("The thickness of the target Layer seems too small for a Squeeze-and-Stretch");
    return grid3Dout;
  }

  // Retrieve information from the 3D input grid
  VectorDouble X0s = grid3Din->getX0s();
  VectorDouble DXs = grid3Din->getDXs();
  VectorInt    NXs = grid3Din->getNXs();
  VectorDouble angles = grid3Din->getAngles();

  // Modify these characteristics for the output 3D Grid
  int nzin  = NXs[idim0];
  NXs[idim0] = nzout;
  DXs[idim0] = dzout;
  X0s[idim0] = z0out;

  // Create the output 3D grid
  grid3Dout = create(NXs, DXs, X0s, angles, ELoadBy::fromKey("SAMPLE"),
                     VectorDouble(), VectorString(), VectorString(), 1, false);

  // Create the variables in the 3D grid and identify their UIDs
  for (int ivar = 0; ivar < nvar; ivar++)
    grid3Dout->addColumnsByConstant(1, TEST, names[ivar]);
  VectorInt iuids = grid3Dout->getUIDs(names);

  // Define local variables
  int iuidTop = surf2D->getUID(nameTop);
  int iuidBot = surf2D->getUID(nameBot);
  VectorDouble vecin(nzin);
  VectorDouble vecout(nzout);
  VectorInt indg(ndim, 0);
  int ig2D;
  double top;
  double bot;

  // Loop on the 3-D vertical columns of the 3-D grid
  for (int ix = 0, nx = grid3Dout->getNX(0); ix < nx; ix++)
    for (int iy = 0, ny = grid3Dout->getNX(1); iy < ny; iy++)
    {
      indg[0] = ix;
      indg[1] = iy;

      // Identify the corresponding node in 'surf2D'
      ig2D = surf2D->indiceToRank(indg);

      if (! surf2D->isActive(ig2D)) continue;

      // Read the Top and Bottom
      top = surf2D->getArray(ig2D, iuidTop);
      bot = surf2D->getArray(ig2D, iuidBot);

      // Loop on the variables to be transformed
      for (int ivar = 0; ivar < nvar; ivar++)
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
VectorVectorInt DbGrid::getLimitsFromVariableExtend(const String &nameTop,
                                                    const String &nameBot,
                                                    const VectorInt &dimExclude) const
{
  int ndim = getNDim();
  VectorVectorInt vec(ndim);
  for (int idim = 0; idim < ndim; idim++)
  {
    vec[idim].resize(2);
    vec[idim][0] = 0;
    vec[idim][1] = getNX(idim);
  }
  if (nameTop.empty() || nameBot.empty()) return vec;

  // Find the set of Min and Max indices of the subgrid

  int nech = getSampleNumber(true);
  VectorInt indmin(ndim,  10000000);
  VectorInt indmax(ndim, -10000000);
  VectorInt indg(ndim);
  int iuid_top = getUID(nameTop);
  int iuid_bot = getUID(nameBot);

  for (int iech = 0; iech < nech; iech++)
  {
    // Discard not relevant pixels
    if (! isActive(iech)) continue;
    double top = getArray(iech, iuid_top);
    double bot = getArray(iech, iuid_bot);
    if (FFFF(top) || FFFF(bot) || bot > top) continue;

    rankToIndice(iech, indg);
    for (int idim = 0; idim < ndim; idim++)
    {
      int indloc = indg[idim];
      if (indloc < indmin[idim]) indmin[idim] = indloc;
      if (indloc > indmax[idim]) indmax[idim] = indloc;
    }
  }

  // Discard the case where the sub-grid does not exist

  bool flag_exist = true;
  for (int idim = 0; idim < ndim && flag_exist; idim++)
  {
    if (indmin[idim] > indmax[idim]) flag_exist = false;
  }
  if (! flag_exist) return vec;

  // Get the sub-grid characteristics

  vec.resize(ndim);
  for (int idim = 0; idim < ndim; idim++)
  {
    vec[idim].resize(2);
    vec[idim][0] = indmin[idim];
    vec[idim][1] = indmax[idim];
  }

  for (int iexc = 0, nexc = (int) dimExclude.size(); iexc < nexc; iexc++)
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
int DbGrid::setSelectionFromVariableExtend(const String &nameTop, const String &nameBot)
{
  // Create the selection new variable
  int iuidSel = addColumnsByConstant(1, 1, "SelLayer", ELoc::SEL);

  if (nameTop.empty() || nameBot.empty()) return -1;

  // Find the set of Min and Max indices of the subgrid

  int nech = getSampleNumber(true);
  int iuid_top = getUID(nameTop);
  int iuid_bot = getUID(nameBot);

  for (int iech = 0; iech < nech; iech++)
  {
    // Discard not relevant pixels
    if (! isActive(iech)) continue;
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
                                 const DbGrid *surf2D,
                                 const String &nameTop,
                                 const String &nameBot,
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
  if (! isSameGrid(surf2D->getGrid()))
  {
    messerr("The grid 'surf2D' and the current one should coincide horizontally");
    return;
  }
  if (names.empty())
  {
    messerr("You must define some variable to be processed");
    return;
  }

  bool limitsDefined = ! nameTop.empty() && ! nameBot.empty();
  int nvar = (int) names.size();

  // Loop on the vertical columns of the 3-D grid

  int ndim  = 3;
  int idim0 = ndim - 1;
  double top = +1.e30;
  double bot = -1.e30;
  double z0 = getX0(idim0);
  double dz = getDX(idim0);
  int nz = getNX(idim0);
  int indzmin = 0; // included
  int indzmax;
  VectorInt indg(ndim, 0);
  VectorDouble vec(nz);
  VectorDouble vecempty(nz, TEST);
  VectorInt iuids = getUIDs(names);

  int nmodif3D = 0;
  int nmodif2D = 0;
  int rank2D   = 0;
  double thick = 0;
  double thickA = 0;
  for (int ix = 0, nx = getNX(0); ix < nx; ix++)
    for (int iy = 0, ny = getNX(1); iy < ny; iy++)
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
        if (! surf2D->isActive(rank2D)) flagRead = false;
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
          thick   = dz * (indzmax - indzmin + 1);
          if (thick > thickA) thickA = thick;
        }
      }

      // Loop on the variables
      if (flagRead)
      {
        for (int ivar = 0; ivar < nvar; ivar++)
        {
          // Partial update
          getGridPileInPlace(iuids[ivar], indg, idim0, vec);

          // Blank out
          for (int iz = 0; iz < indzmin; iz++, nmodif3D++)
            vec[iz] = TEST;
          for (int iz = indzmax; iz < nz; iz++, nmodif3D++)
            vec[iz] = TEST;

          setGridPileInPlace(iuids[ivar], indg, idim0, vec);
        }
      }
      else
      {
        // Complete update
        nmodif2D ++;
        for (int ivar = 0; ivar < nvar; ivar++, nmodif3D++)
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
  return;
}

VectorInt DbGrid::locateDataInGrid(const Db *data,
                                   const VectorInt &rankIn,
                                   bool centered,
                                   bool useSel) const
{
  VectorInt rankOut;

  if (data == nullptr) return VectorInt();

  if (! rankIn.empty())
  {

    // Locate the samples defined by their ranks stored in 'rankIn'

    for (int ip = 0; ip < (int) rankIn.size(); ip++)
    {
      VectorDouble coor = data->getSampleCoordinates(rankIn[ip]);
      rankOut.push_back(coordinateToRank(coor, centered));
    }
  }
  else
  {

    // Locate all samples (using useSel criterion)

    for (int ip = 0, np = data->getSampleNumber(useSel); ip < np; ip++)
    {
      if (data->isActive(ip) || ! useSel)
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
  for (int idim = 0; idim < getNDim(); idim++)
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
int DbGrid::addSelectionFromDbByMorpho(Db *db,
                                       int nmin,
                                       int radius,
                                       int option,
                                       const VectorInt &dilation,
                                       bool verbose,
                                       const NamingConvention &namconv)
{
  if (db == nullptr)
  {
    messerr("You must define a valid Db");
    return 1;
  }

  int nech = getSampleNumber();

  VectorString names = db->getNamesByColIdx({0});
  int iuid = addColumnsByConstant(1);
  if (dbStatisticsInGridTool(db, this, names, EStatOption::NUM, radius, iuid)) return 1;
  VectorDouble stats = getColumnByUID(iuid, false, false);
  for (int iech = 0; iech < nech; iech++)
    stats[iech] = (stats[iech] <= nmin) ? 0. : 1.;
  setColumnByUID(stats, iuid, false);
  setLocatorByUID(iuid, ELoc::Z, 0, true);

  int err = morpho(EMorpho::DILATION, 0.5, 1.5, option, dilation, false, verbose, namconv);

  deleteColumnByUID(iuid);
  return err;
}

void DbGrid::getSampleAsST(int iech, SpaceTarget& P) const
{
  Db::getSampleAsST(iech, P);

  // Load the target extension
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
VectorVectorDouble DbGrid::getDiscretizedBlock(const VectorInt &ndiscs,
                                               int iech,
                                               bool flagPerCell,
                                               bool flagRandom,
                                               int seed) const
{
  int ndim = getNDim();
  int ntot = VH::product(ndiscs);
  int memo = law_get_random_seed();
  law_set_random_seed(seed);
  VectorVectorDouble discs(ntot);
  for (int i = 0; i < ntot; i++) discs[i].resize(ndim);

  /* Loop on the discretization points */

  for (int i = 0; i < ntot; i++)
  {
    int jech = i;
    int nval = ntot;
    for (int idim = ndim - 1; idim >= 0; idim--)
    {
      double taille = (!flagPerCell) ? getDX(idim) : getLocVariable(ELoc::BLEX, iech, idim);
      int nd = ndiscs[idim];
      nval /= nd;
      int j = jech / nval;
      jech -= j * nval;
      double local = taille * ((j + 0.5) / nd - 0.5);
      if (!flagRandom)
        discs[i][idim] = local;
      else
        discs[i][idim] = local + taille * law_uniform(-0.5, 0.5) / (double) nd;
    }
  }
  law_set_random_seed(memo);

  return discs;
}

