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
#include "Db/Dbgrid.hpp"
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
#include "Basic/GlobalEnvironment.hpp"
#include "Stats/Classical.hpp"

#include <algorithm>
#include <functional>
#include <math.h>

Dbgrid::Dbgrid()
    : Db(),
      _grid(0)
{
  _clear();
}

Dbgrid::Dbgrid(const Dbgrid& r)
    : Db(r),
      _grid(r._grid)
{
}

Dbgrid& Dbgrid::operator=(const Dbgrid& r)
{
  if (this != &r)
  {
    Dbgrid::operator=(r);
    _grid = r._grid;
  }
  return *this;
}

Dbgrid::~Dbgrid()
{
}

String Dbgrid::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  const DbStringFormat* dbfmt = dynamic_cast<const DbStringFormat*>(strfmt);
  DbStringFormat dsf;
  if (dbfmt != nullptr) dsf = *dbfmt;

  sstr << toTitle(0, "Data Base Characteristics");

  if (dsf.matchResume())
  {
    sstr << _summaryString();
    sstr << _grid.toString();
  }

  sstr << _toStringCommon(strfmt);

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
int Dbgrid::reset(const VectorInt& nx,
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
  int natt = (tab.empty()) ? 0 : (int) (tab.size() / nech);
  int ncol = ndim + natt + flag_add_rank;
  resetDims(ncol, nech);

  // Create the grid

  if (gridDefine(nx, dx, x0, angles)) return 1;

  // Load the data

  if (flag_add_rank) _createRank(0);
  _createCoordinatesGrid(flag_add_rank);
  _loadData(tab, names, locatorNames, order, flag_add_rank);

  // Create the coordinate names (for the remaining variables)

  _defineDefaultNames(flag_add_rank + ndim, names);

  // Create the locators

  int jcol = 0;
  if (flag_add_rank) jcol++;
  setLocatorsByAttribute(ndim, jcol, ELoc::X);
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
int Dbgrid::resetCoveringDb(Db* db,
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
  setLocatorsByAttribute(ndim, jcol, ELoc::X);

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
int Dbgrid::resetFromPolygon(Polygons* polygon,
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
  setLocatorsByAttribute(ndim, jcol, ELoc::X);

  return 0;
}

Dbgrid* Dbgrid::create(const VectorInt& nx,
                       const VectorDouble& dx,
                       const VectorDouble& x0,
                       const VectorDouble& angles,
                       const ELoadBy& order,
                       const VectorDouble& tab,
                       const VectorString& names,
                       const VectorString& locatorNames,
                       int flag_add_rank)
{
  Dbgrid* dbgrid = new Dbgrid;
  if (dbgrid->reset(nx, dx, x0, angles, order, tab, names, locatorNames,
                    flag_add_rank))
  {
    messerr("Error when creating Dbgrid from Grid");
    delete dbgrid;
    return nullptr;
  }
  return dbgrid;
}

Dbgrid* Dbgrid::createCoveringDb(Db* db,
                                 const VectorInt& nodes,
                                 const VectorDouble& dcell,
                                 const VectorDouble& origin,
                                 const VectorDouble& margin)
{
  Dbgrid* dbgrid = new Dbgrid;
  if (dbgrid->resetCoveringDb(db, nodes, dcell, origin, margin))
  {
    messerr("Error when creating Dbgrid covaring another Db");
    delete dbgrid;
    return nullptr;
  }
  return dbgrid;

}
Dbgrid* Dbgrid::createFromPolygon(Polygons* polygon,
                                  const VectorInt& nodes,
                                  const VectorDouble& dcell,
                                  int flag_add_rank)
{
  Dbgrid* dbgrid = new Dbgrid;
  if (dbgrid->resetFromPolygon(polygon, nodes, dcell, flag_add_rank))
  {
    messerr("Error when creating Dbgrid from Polygon");
    delete dbgrid;
    return nullptr;
  }
  return dbgrid;
}

/**
 * Paint the ndim columns starting from 'icol0' with grid coordinates
 * @param icol0 Starting column
 */
void Dbgrid::_createCoordinatesGrid(int icol0)
{
  // Set the Names

  for (int idim = 0; idim < getNDim(); idim++)
    _setNameByColumn(icol0 + idim, getLocatorName(ELoc::X, idim));

  // Set the locators

  setLocatorsByAttribute(getNDim(), icol0, ELoc::X);

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

bool Dbgrid::isSameGrid(const Grid& grid) const
{
  if (grid.empty())
  {
    messerr("Both files should be organized as grids");
    return false;
  }
  return _grid.isSame(grid);
}

void Dbgrid::gridCopyParams(int mode, const Grid& gridaux)
{
  _grid.copyParams(mode, gridaux);
}

bool Dbgrid::isSameGridMesh(const Dbgrid& dbaux) const
{
  if (! dbaux.isGrid())
  {
    messerr("Both files should be organized as grids");
    return false;
  }
  return _grid.isSameMesh(dbaux.getGrid());
}

bool Dbgrid::isSameGridMeshOldStyle(const Dbgrid* dbaux) const
{
  if (! dbaux->isGrid())
  {
    messerr("Both files should be organized as grids");
    return false;
  }
  return _grid.isSameMesh(dbaux->getGrid());
}


bool Dbgrid::isSameGridRotation(const Dbgrid& dbaux) const
{
  if (! dbaux.isGrid())
  {
    messerr("Both files should be organized as grids");
    return false;
  }
  if (! isGridRotated() && ! dbaux.isGridRotated()) return true;
  return _grid.isSameRotation(dbaux.getGrid());
}

bool Dbgrid::isSameGridRotationOldStyle(const Dbgrid* dbaux) const
{
  if (! dbaux->isGrid())
  {
    messerr("Both files should be organized as grids");
    return false;
  }
  if (! isGridRotated() && ! dbaux->isGridRotated()) return true;
  return _grid.isSameRotation(dbaux->getGrid());
}

bool Dbgrid::isGridRotated() const
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
double Dbgrid::getCoordinate(int iech, int idim, bool flag_rotate) const
{
  if (idim >= getNDim()) return TEST;
  return _grid.getCoordinate(iech, idim, flag_rotate);
}

int Dbgrid::getNDim() const
{
  return (_grid.getNDim());
}

int Dbgrid::_deserialize(FILE* file, bool /*verbose*/)
{
  int ndim, ndim2, ntot, natt, nech, i, flag_grid;
  VectorInt tabnum;
  std::vector<ELoc> tabatt;
  VectorInt nx;
  VectorString tabnam;
  VectorDouble x0;
  VectorDouble dx;
  VectorDouble angles;
  VectorDouble tab;
  static int flag_add_rank = 1;

  /* Initializations */

  natt = ndim = nech = ntot = 0;

  /* Decoding the header */

  if (_recordRead(file, "Space Dimension", "%d", &ndim)) goto label_end;

  /* Core allocation */

  nx.resize(ndim);
  dx.resize(ndim);
  x0.resize(ndim);
  angles.resize(ndim);

  /* Read the grid characteristics */

  for (int idim = 0; idim < ndim; idim++)
  {
    if (_recordRead(file, "Grid Number of Nodes", "%d", &nx[idim]))
      goto label_end;
    if (_recordRead(file, "Grid Origin", "%lf", &x0[idim])) goto label_end;
    if (_recordRead(file, "Grid Mesh", "%lf", &dx[idim])) goto label_end;
    if (_recordRead(file, "Grid Angles", "%lf", &angles[idim])) goto label_end;
  }
  ntot = ut_ivector_prod(nx);

  /* Reading the tail of the file */

  _variableRead(file, &natt, &ndim2, &nech, tabatt, tabnum, tabnam, tab);

  /* Creating the Db */

  if (natt > 0 && nech != ntot)
  {
    messerr("The number of lines read from the Grid file (%d)", nech);
    messerr("is not a multiple of the number of samples (%d)", ntot);
    messerr("The Grid Db is created with no sample attached");
    natt = 0;
  }
  resetDims(natt + flag_add_rank, ut_ivector_prod(nx));
  (void) gridDefine(nx, dx, x0, angles);
  _loadData(ELoadBy::SAMPLE, flag_add_rank, tab);

  /* Loading the names */

  if (natt > 0) for (i = 0; i < natt; i++)
    setNameByAttribute(i + flag_add_rank, tabnam[i]);

  /* Create the locators */

  if (natt > 0) for (i = 0; i < natt; i++)
    setLocatorByAttribute(i + flag_add_rank, tabatt[i], tabnum[i]);

  /* Core deallocation */

  label_end: return 0;
}

int Dbgrid::_serialize(FILE* file, bool /*verbose*/) const
{
  bool onlyLocator = false;
  bool writeCoorForGrid = true;

  /* Writing the header */

  _recordWrite(file, "%d", getNDim());
  _recordWrite(file, "#", "Space Dimension");

  /* Writing the grid characteristics */

  _recordWrite(file, "#", "Grid characteristics (NX,X0,DX,ANGLE)");
  for (int idim = 0; idim < getNDim(); idim++)
  {
    _recordWrite(file, "%d",  getNX(idim));
    _recordWrite(file, "%lf", getX0(idim));
    _recordWrite(file, "%lf", getDX(idim));
    _recordWrite(file, "%lf", getAngle(idim));
    _recordWrite(file, "\n");
  }

  /* Writing the tail of the file */

  if (_variableWrite(file, true, onlyLocator, writeCoorForGrid)) return 1;

  return 0;
}

double Dbgrid::getUnit(int idim) const
{
  return _grid.getDX(idim);
}

int Dbgrid::gridDefine(const VectorInt& nx,
                       const VectorDouble& dx,
                       const VectorDouble& x0,
                       const VectorDouble& angles)
{
  return (_grid.resetFromVector(nx, dx, x0, angles));
}

int Dbgrid::dumpToNF(const String& neutralFilename, bool verbose) const
{
  FILE* file = _fileOpen(neutralFilename, "Dbgrid", "w", verbose);
  if (file == nullptr) return 1;

  if (_serialize(file, verbose))
  {
    if (verbose) messerr("Problem writing in the Neutral File.");
    _fileClose(file, verbose);
    return 1;
  }
  _fileClose(file, verbose);
  return 0;
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
Dbgrid* Dbgrid::createFromNF(const String& neutralFilename, bool verbose)
{
  FILE* file = _fileOpen(neutralFilename, "Dbgrid", "r", verbose);
  if (file == nullptr) return nullptr;

  Dbgrid* db = new Dbgrid;
  if (db->_deserialize(file, verbose))
  {
    if (verbose) messerr("Problem reading the Neutral File.");
    delete db;
    db = nullptr;
  }
  _fileClose(file, verbose);
  return db;
}

VectorDouble Dbgrid::getFieldSubGrid(const String& name,
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

void Dbgrid::generateCoordinates(const String& radix)
{
  if (! isGrid())
  {
    messerr("This method is only available in the case of Grid. Nothing done");
    return;
  }
  int ndim = getNDim();
  VectorDouble coors(ndim);
  (void) addFieldsByConstant(ndim, 0., radix, ELoc::X);
  display();
  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    _grid.rankToCoordinatesInPlace(iech, coors);
    for (int idim = 0; idim < ndim; idim++)
      setCoordinate(iech, idim, coors[idim]);
  }
}

