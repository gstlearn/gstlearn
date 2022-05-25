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

int DbGrid::_deserialize(std::istream& is, bool /*verbose*/)
{
  int ndim, ntot, nech, i, ncol;
  VectorInt nx;
  VectorString locators;
  VectorString names;
  VectorDouble x0;
  VectorDouble dx;
  VectorDouble angles;
  VectorDouble values;
  VectorDouble allvalues;

  /* Initializations */

  ndim = nech = ntot = ncol = 0;

  /* Decoding the header */

  bool ret = _recordRead<int>(is, "Space Dimension", ndim);

  /* Core allocation */

  nx.resize(ndim);
  dx.resize(ndim);
  x0.resize(ndim);
  angles.resize(ndim);

  /* Read the grid characteristics */

  for (int idim = 0; idim < ndim; idim++)
  {
    ret = ret && _recordRead<int>(is, "Grid Number of Nodes", nx[idim]);
    ret = ret && _recordRead<double>(is, "Grid Origin", x0[idim]);
    ret = ret && _recordRead<double>(is, "Grid Mesh", dx[idim]);
    ret = ret && _recordRead<double>(is, "Grid Angles", angles[idim]);
  }
  ntot = ut_vector_prod(nx);

  ret = ret && _recordRead<int>(is, "Number of variables", ncol);
  if (ncol > 0)
  {
    ret = ret && _recordReadVec<String>(is, "Locators", locators);
    if (!ret || (int) locators.size() != ncol) return 1;
    ret = ret && _recordReadVec<String>(is, "Names", names);
    if (!ret || (int) names.size() != ncol) return 1;
  }

  /* Reading the tail of the file */

  while (ret)
  {
    ret = _recordReadVec<double>(is, "", values);
    if (ret)
    {
      if ((int)values.size() != ncol) return 1;
      // Concatenate values by samples
      allvalues.insert(allvalues.end(), std::make_move_iterator(values.begin()),
                                        std::make_move_iterator(values.end()));
      nech++;
    } // else "end of file"
  }

  // Decode the locators
  std::vector<ELoc> tabloc;
  VectorInt tabnum;
  int  inum = 0, mult = 0;
  ELoc iloc;
  for (auto loc : locators)
  {
    if (locatorIdentify(loc, &iloc, &inum, &mult)) return 1;
    tabloc.push_back(iloc);
    tabnum.push_back(inum);
  }

  /* Creating the Db */

  if (ncol > 0 && nech != ntot)
  {
    messerr("The number of lines read from the Grid file (%d)", nech);
    messerr("is not a multiple of the number of samples (%d)", ntot);
    messerr("The Grid Db is created with no sample attached");
    ncol = 0;
  }

  resetDims(ncol, ut_vector_prod(nx));
  (void) gridDefine(nx, dx, x0, angles);

  // Load the values
  _loadData(ELoadBy::SAMPLE, 0, allvalues);
  // Update the column names and locators
  if (ncol > 0)
    for (i = 0; i < ncol; i++)
    {
      setNameByUID(i, names[i]);
      setLocatorByUID(i, tabloc[i], tabnum[i]);
    }

  return 0;
}

int DbGrid::_serialize(std::ostream& os, bool verbose) const
{

  /* Writing the header */

  bool ret = _recordWrite<int>(os, "Space Dimension", getNDim());

  /* Writing the grid characteristics */

  ret = ret && _commentWrite(os, "Grid characteristics (NX,X0,DX,ANGLE)");
  for (int idim = 0; idim < getNDim(); idim++)
  {
    ret = ret && _recordWrite<int>(os, "",  getNX(idim));
    ret = ret && _recordWrite<double>(os, "", getX0(idim));
    ret = ret && _recordWrite<double>(os, "", getDX(idim));
    ret = ret && _recordWrite<double>(os, "", getAngle(idim));
    ret = ret && _commentWrite(os, "");
  }

  /* Writing the tail of the file */

  if (Db::_serialize(os, verbose)) return 1;

  return 0;
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

int DbGrid::dumpToNF(const String& neutralFilename, bool verbose) const
{
  std::ofstream os;
  int ret = 1;
  if (_fileOpenWrite(neutralFilename, "DbGrid", os, verbose))
  {
    ret = _serialize(os, verbose);
    if (ret && verbose) messerr("Problem writing in the Neutral File.");
    os.close();
  }
  return ret;
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
  DbGrid* db = nullptr;
  std::ifstream is;
  if (_fileOpenRead(neutralFilename, "DbGrid", is, verbose))
  {
    db = new DbGrid;
    if (db->_deserialize(is, verbose))
    {
      if (verbose) messerr("Problem reading the Neutral File.");
      delete db;
      db = nullptr;
    }
    is.close();
  }
  return db;
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
