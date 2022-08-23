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
#include "Basic/Vector.hpp"
#include "Basic/Law.hpp"
#include "Basic/AException.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/GlobalEnvironment.hpp"
#include "Stats/Classical.hpp"
#include "Estimation/CalcImage.hpp"

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
    DbGrid::operator=(r);
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
 * @param flag_add_rank If 1, add an automatic rank variable
 */
int DbGrid::reset(const VectorInt& nx,
                  const VectorDouble& dx,
                  const VectorDouble& x0,
                  const VectorDouble& angles,
                  const ELoadBy& order,
                  const VectorDouble& tab,
                  const VectorString& names,
                  const VectorString& locatorNames,
                  int flag_add_rank)
{
  _clear();

  int ndim = static_cast<int> (nx.size());
  int nech = 1;
  for (int idim = 0; idim < ndim; idim++)
    nech *= nx[idim];
  int ntab = (tab.empty()) ? 0 : (int) (tab.size() / nech);
  int ncol = ndim + ntab + flag_add_rank;
  resetDims(ncol, nech);

  // Create the grid

  if (gridDefine(nx, dx, x0, angles)) return 1;

  // Load the data

  _loadData(tab, names, locatorNames, order, ndim + flag_add_rank);

  // Additional fields

  if (flag_add_rank) _createRank(0);
  _createCoordinatesGrid(flag_add_rank);

  // Create the names (for the remaining variables)

  _defineDefaultNames(flag_add_rank + ndim, names);

  // Create the locators

  int jcol = 0;
  if (flag_add_rank) jcol++;
  setLocatorsByUID(ndim, jcol, ELoc::X);
  _defineDefaultLocators(flag_add_rank + ndim, locatorNames);

  return 0;
}

/**
 * Creating a Grid Db which covers the extension of the input 'Db'
 *
 * @param db       Input Db from which the newly created Db is constructed
 * @param nodes    Vector of the expected number of grid nodes (default = 10)
 * @param dcell    Vector of the expected sizes for the grid meshes
 * @param origin   Vector of the expected origin of the grid
 * @param margin   Vector of the expected margins of the grid
 *
 * @remarks Arguments 'nodes' and 'dcell' are disjunctive. If both defined, 'dcell' prevails
 */
int DbGrid::resetCoveringDb(Db* db,
                            const VectorInt& nodes,
                            const VectorDouble& dcell,
                            const VectorDouble& origin,
                            const VectorDouble& margin)
{
  _clear();
  int ndim = db->getNDim();

  // Derive the Grid parameters

  VectorInt    nx(ndim);
  VectorDouble x0(ndim);
  VectorDouble dx(ndim);
  int nech = 1;
  for (int idim = 0; idim < ndim; idim++)
  {
    VectorDouble coor = db->getExtrema(idim);

    double marge = 0.;
    if (ndim == (int) margin.size()) marge = margin[idim];

    double x0loc = coor[0];
    if (ndim == (int) origin.size()) x0loc = origin[idim];
    x0loc -= marge;

    double ext = coor[1] - x0loc + 2. * marge;

    int nxloc = 10;
    double dxloc = ext / (double) nxloc;

    // Constraints specified by the number of nodes
    if (ndim == (int) nodes.size())
    {
      nxloc = nodes[idim];
      dxloc = ext / (double) nxloc;
    }

    // Constraints specified by the cell sizes
    if (ndim == (int) dcell.size())
    {
      dxloc = dcell[idim];
      nxloc = static_cast<int> (ext / dxloc) + 1;
      ++nxloc; // one more node than intervals
    }

    nx[idim] = nxloc;
    dx[idim] = dxloc;
    x0[idim] = x0loc;
    nech *= nxloc;
  }

  resetDims(ndim,nech);

  // Create the grid

  if (gridDefine(nx, dx, x0)) return 1;

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
 * @param flag_add_rank 1 if the sample rank must be generated
 */
int DbGrid::resetFromPolygon(Polygons* polygon,
                             const VectorInt& nodes,
                             const VectorDouble& dcell,
                             int flag_add_rank)
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

  int ncol = ndim + flag_add_rank;
  resetDims(ncol, nech);

  // Create the grid

  if (gridDefine(nx_tab, dx_tab, x0_tab)) return 1;

  /// Load the data

  if (flag_add_rank) _createRank(0);
  _createCoordinatesGrid(flag_add_rank);

  // Create the locators

  int jcol = 0;
  if (flag_add_rank) jcol++;
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
                       int flag_add_rank)
{
  DbGrid* dbgrid = new DbGrid;
  if (dbgrid->reset(nx, dx, x0, angles, order, tab, names, locatorNames,
                    flag_add_rank))
  {
    messerr("Error when creating DbGrid from Grid");
    delete dbgrid;
    return nullptr;
  }
  return dbgrid;
}

DbGrid* DbGrid::createCoveringDb(Db* db,
                                 const VectorInt& nodes,
                                 const VectorDouble& dcell,
                                 const VectorDouble& origin,
                                 const VectorDouble& margin)
{
  DbGrid* dbgrid = new DbGrid;
  if (dbgrid->resetCoveringDb(db, nodes, dcell, origin, margin))
  {
    messerr("Error when creating DbGrid covaring another Db");
    delete dbgrid;
    return nullptr;
  }
  return dbgrid;

}
DbGrid* DbGrid::createFromPolygon(Polygons* polygon,
                                  const VectorInt& nodes,
                                  const VectorDouble& dcell,
                                  int flag_add_rank)
{
  DbGrid* dbgrid = new DbGrid;
  if (dbgrid->resetFromPolygon(polygon, nodes, dcell, flag_add_rank))
  {
    messerr("Error when creating DbGrid from Polygon");
    delete dbgrid;
    return nullptr;
  }
  return dbgrid;
}

DbGrid* DbGrid::coarsify(const VectorInt &nmult)
{
  return createCoarse(this,nmult,0);
}

DbGrid* DbGrid::createCoarse(DbGrid *dbin,
                             const VectorInt &nmult,
                             int flag_add_rank)
{
  DbGrid *dbgrid = new DbGrid;
  int ndim = dbin->getNDim();

  // Get the new grid characteristics
  VectorInt nx(ndim);
  VectorDouble dx(ndim);
  VectorDouble x0(ndim);
  dbin->getGrid().multiple(nmult, 1, nx, dx, x0);

  // Create the new grid
  dbgrid = create(nx, dx, x0, dbin->getAngles(), ELoadBy::SAMPLE,
                  VectorDouble(), VectorString(), VectorString(), flag_add_rank);

  // Migrate all variables (except 'rank'  and coordinates
  (void) migrateAllVariables(dbin, dbgrid, flag_add_rank);

  return dbgrid;
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
                             int flag_add_rank)
{
  DbGrid *dbgrid = new DbGrid;
  int ndim = dbin->getNDim();

  // Get the new grid characteristics
  VectorInt nx(ndim);
  VectorDouble dx(ndim);
  VectorDouble x0(ndim);
  dbin->getGrid().divider(nmult, 1, nx, dx, x0);

  // Create the new grid
  dbgrid = create(nx, dx, x0, dbin->getAngles(), ELoadBy::SAMPLE,
                  VectorDouble(), VectorString(), VectorString(), flag_add_rank);

  // Migrate all variables (except 'rank'  and coordinates
  (void) migrateAllVariables(dbin, dbgrid, flag_add_rank);

  return dbgrid;
}

bool DbGrid::migrateAllVariables(Db *dbin, Db *dbout, int flag_add_rank)
{
  ELoc locatorType;
  int  locatorIndex;

  // Constitute the list of Variables to be migrated

  VectorInt icols;
  for (int icol = 0; icol < dbin->getColumnNumber(); icol++)
  {
    // Skip the rank
    if (flag_add_rank && icol == 0) continue;

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
  if (migrateByAttribute(dbin, dbout, icols,
                         0, VectorDouble(), true, true,
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
  // Set the Names

  for (int idim = 0; idim < getNDim(); idim++)
    _setNameByColIdx(icol0 + idim, getLocatorName(ELoc::X, idim));

  // Set the locators

  setLocatorsByUID(getNDim(), icol0, ELoc::X);

  // Generate the vector of coordinates

  _grid.iteratorInit();
  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    VectorInt indices = _grid.iteratorNext();
    VectorDouble coors = _grid.indicesToCoordinate(indices);
    for (int idim = 0; idim < getNDim(); idim++)
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

bool DbGrid::isSameGridMeshOldStyle(const DbGrid* dbaux) const
{
  if (! dbaux->isGrid())
  {
    messerr("Both files should be organized as grids");
    return false;
  }
  return _grid.isSameMesh(dbaux->getGrid());
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

bool DbGrid::isSameGridRotationOldStyle(const DbGrid* dbaux) const
{
  if (! dbaux->isGrid())
  {
    messerr("Both files should be organized as grids");
    return false;
  }
  if (! isGridRotated() && ! dbaux->isGridRotated()) return true;
  return _grid.isSameRotation(dbaux->getGrid());
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

int DbGrid::getNDim() const
{
  return (_grid.getNDim());
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

  ret && Db::_deserialize(is, verbose);

  // Create the Grid coordinates
  (void) gridDefine(nx, dx, x0, angles);

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
    success =  dbgrid->deserialize(is, verbose);
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
  if (useSel) sel = getSelection();

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

/**
 * Extracts a slice from a 3-D Grid
 * @param name   Name of the target variable
 * @param pos    Type of section: 0 for YoZ; 1 for XoZ and 2 for XoY
 * @param indice Rank of the section
 * @param useSel Use the active selection
 * @return A VectorVectorDouble with 4 columns, i.e: X, Y, Z, Var
 *
 * @remark In presence of a selection, values are returned but set to TEST
 * @remark (if useSel)
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
        getCoordinatesInPlace(iech, coor);
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
        getCoordinatesInPlace(iech, coor);
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
        getCoordinatesInPlace(iech, coor);
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
 * @return
 */
VectorVectorDouble DbGrid::getCellEdges(int node) const
{
  VectorVectorDouble coords(2);
  coords[0].resize(5);
  coords[1].resize(5);

  int ndim = getNDim();
  VectorInt icorner(ndim,0);
  VectorDouble local;

  // Get the extension of the target cell (possibly variable)
  VectorDouble dxsPerCell = getBlockExtensions(node);

  icorner[0] = -1;
  icorner[1] = -1;
  local = getGrid().getCellCoordinatesByCorner(node, icorner, dxsPerCell);
  coords[0][0] = local[0];
  coords[1][0] = local[1];

  icorner[0] = -1;
  icorner[1] = 1;
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

VectorDouble DbGrid::getBlockExtensions(int node) const
{
  int ndim = getNDim();

  VectorDouble dxsPerCell = getDXs();
  if (hasBlockExtension())
  {
    for (int idim = 0; idim < ndim; idim++)
    {
      double value = getBlockExtension(node, idim);
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
                   bool dist_erode,
                   bool verbose,
                   const NamingConvention &namconv)
{
  return dbMorpho(this, oper, vmin, vmax, option, radius, dist_erode, verbose, namconv);
}

