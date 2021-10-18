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
#include "Db/Db.hpp"
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
#include "geoslib_f_private.h"
#include "geoslib_old_f.h"
#include "geoslib_f.h"
#include <algorithm>
#include <functional>

Db::Db()
    : AStringable(),
      ASerializable(),
      _ncol(0),
      _nech(0),
      _isGrid(0),
      _array(),
      _attcol(),
      _colNames(),
      _p(),
      _grid(0)
{
  _initP();
}

Db::Db(int nech,
       const ELoadBy& order,
       const VectorDouble& tab,
       const VectorString& names,
       const VectorString& locatorNames,
       int flag_add_rank)
    : AStringable(),
      ASerializable(),
      _ncol(0),
      _nech(0),
      _isGrid(false),
      _array(),
      _attcol(),
      _colNames(),
      _p(),
      _grid(0)
{
  _initP();
  int natt = (tab.empty()) ? 0 : (int) (tab.size() / nech);
  _ncol = natt + flag_add_rank;
  _nech = nech;
  reset(_ncol, _nech);

  // Load data (if defined)

  if (flag_add_rank) _createRank(0);
  _loadData(tab, names, locatorNames, order, flag_add_rank);

}

Db::Db(const VectorInt& nx,
       const VectorDouble& dx,
       const VectorDouble& x0,
       const VectorDouble& angles,
       const ELoadBy& order,
       const VectorDouble& tab,
       const VectorString& names,
       const VectorString& locatorNames,
       int flag_add_rank)
    : AStringable(),
      ASerializable(),
      _ncol(0),
      _nech(0),
      _isGrid(true),
      _array(),
      _attcol(),
      _colNames(),
      _p(),
      _grid(0)
{
  _initP();
  int ndim = static_cast<int> (nx.size());
  int nech = 1;
  for (int idim = 0; idim < ndim; idim++)
    nech *= nx[idim];
  int natt = (tab.empty()) ? 0 : (int) (tab.size() / _nech);

  _nech = nech;
  _ncol = ndim + natt + flag_add_rank;
  reset(_ncol, _nech);

  // Create the grid

  if (gridDefine(nx, dx, x0, angles)) return;

  /// Load the data

  if (flag_add_rank) _createRank(0);
  _createGridCoordinates(flag_add_rank);
  _loadData(tab, names, locatorNames, order, flag_add_rank);

  // Create the coordinate names (for the remaining variables)

  _defineDefaultNames(flag_add_rank + ndim, names);

  // Create the locators

  int jcol = 0;
  if (flag_add_rank) jcol++;
  setLocatorsByAttribute(ndim, jcol, ELoc::X);
  _defineDefaultLocators(flag_add_rank + ndim, locatorNames);
}

/**
 * Creating a Db by reading a CSV file
 *
 * @param filename   Name of the CSV file
 * @param verbose    Verbose flag
 * @param csv        Description of the CSV format
 * @param ncol_max   Maximum number of columns
 * @param nrow_max   Maximum number of rows
 * @param flag_add_rank 1 if the sample rank must be generated
 */
Db::Db(const String& filename,
       bool verbose,
       const CSVformat& csv,
       int ncol_max,
       int nrow_max,
       int flag_add_rank)
    : AStringable(),
      ASerializable(),
      _ncol(0),
      _nech(0),
      _isGrid(false),
      _array(),
      _attcol(),
      _colNames(),
      _p(),
      _grid(0)
{
  _initP();
  VectorString names;
  VectorDouble tab;
  int ncol, nrow;

  /* Reading the CSV file */

  if (csv_table_read2(filename, (int) verbose,
                      csv.getFlagHeader(), csv.getNSkip(),
                      csv.getCharSep(), csv.getCharDec(),csv.getNaString(),
                      ncol_max, nrow_max, &ncol, &nrow, names, tab))
  {
    messerr("Problem when reading CSV file");
    return;
  }

  int natt = (tab.empty()) ? 0 : (int) (tab.size() / nrow);
  _ncol = natt + flag_add_rank;
  _nech = nrow;
  reset(_ncol, _nech);

  // Load data (if defined)

  if (flag_add_rank) _createRank(0);
  _loadData(tab, names, VectorString(), ELoadBy::SAMPLE, flag_add_rank);

  // Set the names

  _defineDefaultNames(flag_add_rank, names);

  // Locators: Try to guess them from the Names
  _defineDefaultLocatorsByNames(flag_add_rank, names);
}

/**
 * Creating a grid Db which covers the extension of the input 'Db'
 *
 * @param db       Input Db from which the newly created Db is constructed
 * @param nodes    Vector of the expected number of grid nodes (default = 10)
 * @param dcell    Vector of the expected sizes for the grid meshes
 * @param origin   Vector of the expected origin of the grid
 * @param margin   Vector of the expected margins of the grid
 * @param flag_add_rank 1 if the sample rank must be generated
 */
Db::Db(Db* db,
       const VectorInt&    nodes,
       const VectorDouble& dcell,
       const VectorDouble& origin,
       const VectorDouble& margin,
       int flag_add_rank)
    : AStringable(),
      ASerializable(),
      _ncol(0),
      _nech(0),
      _isGrid(false),
      _array(),
      _attcol(),
      _colNames(),
      _p(),
      _grid(0)
{
  _initP();
  int ndim = db->getNDim();

  // Derive the Grid parameters

  VectorInt nx_tab;
  VectorDouble x0_tab;
  VectorDouble dx_tab;
  int nech = 1;
  for (int idim = 0; idim < ndim; idim++)
  {
    VectorDouble coor = db->getExtrema(idim);

    double x0 = coor[0];
    double ext = coor[1] - coor[0];
    double marge = 0.;

    if (ndim == (int) margin.size()) marge = margin[idim];
    if (ndim == (int) origin.size()) x0 = origin[idim];
    x0 -= marge;
    ext += 2. * marge;

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

  _nech = nech;
  _ncol = ndim + flag_add_rank;
  reset(_ncol, _nech);

  // Create the grid

  if (gridDefine(nx_tab, dx_tab, x0_tab)) return;

  /// Load the data

  if (flag_add_rank) _createRank(0);
  _createGridCoordinates(flag_add_rank);

  // Create the locators

  int jcol = 0;
  if (flag_add_rank) jcol++;
  setLocatorsByAttribute(ndim, jcol, ELoc::X);
}

/**
 * Create a Db by loading the contents of a Neutral File
 *
 * @param neutralFileName Name of the Neutral File (Db format)
 * @param verbose         Verbosity flag
 */
Db::Db(const String& neutralFileName, bool verbose)
    : AStringable(),
      ASerializable(),
      _ncol(0),
      _nech(0),
      _isGrid(false),
      _array(),
      _attcol(),
      _colNames(),
      _p(),
      _grid(0)
{
  if (deSerialize(neutralFileName, verbose))
  {
    messerr("Problem reading the Neutral File.");
    messerr("The Db is not entirely created");
    reset(0,0);
  }
}

/**
 * Creating a regular grid Db which covers the input Polygon
 *
 * @param polygon    Pointer to the input Polygon
 * @param nodes      Vector of the expected number of nodes
 * @param dcell      Vector of the expected dimensions for the grid cells
 * @param flag_add_rank 1 if the sample rank must be generated
 */
Db::Db(Polygons* polygon,
       const VectorInt& nodes,
       const VectorDouble& dcell,
       int flag_add_rank)
    : AStringable(),
      ASerializable(),
      _ncol(0),
      _nech(0),
      _isGrid(false),
      _array(),
      _attcol(),
      _colNames(),
      _p(),
      _grid(0)
{
  _initP();
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

  _nech = nech;
  _ncol = ndim + flag_add_rank;
  reset(_ncol, _nech);

  // Create the grid

  if (gridDefine(nx_tab, dx_tab, x0_tab)) return;

  /// Load the data

  if (flag_add_rank) _createRank(0);
  _createGridCoordinates(flag_add_rank);

  // Create the locators

  int jcol = 0;
  if (flag_add_rank) jcol++;
  setLocatorsByAttribute(ndim, jcol, ELoc::X);
}

/**
 * Sampling an input Db to create the output Db by selecting a subset of samples
 *
 * @param dbin       Pointer to the input Db
 * @param proportion Proportion of samples to be retained
 * @param names      Vector of Names to be copied (empty: all names)
 * @param seed       Seed used for the random number generator
 * @param verbose    Verbose flag
 *
 * @remark A possible selection in 'dbin' will not be taken into account
 */
Db::Db(const Db* dbin,
       double proportion,
       const VectorString& names,
       int seed,
       bool verbose)
    : AStringable(),
      ASerializable(),
      _ncol(0),
      _nech(0),
      _isGrid(false),
      _array(),
      _attcol(),
      _colNames(),
      _p(),
      _grid(0)
{
  _initP();

  // Creating the vector of selected samples

  int nfrom = dbin->getSampleNumber();
  VectorInt ranks = ut_vector_sample(nfrom, proportion);
  _nech = static_cast<int> (ranks.size());
  if (verbose)
    message("From %d samples, the extraction concerns %d samples\n", nfrom,
            _nech);

  // Create the new data base

  VectorString namloc = names;
  if (namloc.empty())
    namloc = dbin->getNames();
  _ncol = static_cast<int> (namloc.size());
  reset(_ncol, _nech);

  // Create Variables and Locators

  ELoc locatorType;
  int locatorIndex;
  for (int icol = 0; icol < _ncol; icol++)
  {
    setName(icol, namloc[icol]);
    if (dbin->getLocator(namloc[icol],&locatorType,&locatorIndex))
      setLocator(namloc[icol],locatorType,locatorIndex);
  }

  // Load samples

  VectorDouble values(_nech);
  for (int icol = 0; icol < _ncol; icol++)
  {
    int jcol = dbin->getColumn(namloc[icol]);
    for (int iech = 0; iech < _nech; iech++)
      values[iech] = dbin->getByColumn(ranks[iech],jcol);
    setColumnByRank(values, icol);
  }
}

/**
 * Create a Db generating samples randomly
 *
 * @param nech    Number of samples to be generated
 * @param coormin Vector giving the smallest values of the coordinates
 * @param coormax Vector giving the largest values for the coordinates
 * @param ndim    Space dimension (used if 'coormin' and 'coormax' are empty)
 * @param seed    Seed for the random number generator
 * @param flag_add_rank 1 if the Sample ranks must be generated
 */
Db::Db(int nech,
       const VectorDouble& coormin,
       const VectorDouble& coormax,
       int ndim,
       int seed,
       int flag_add_rank)
    : AStringable(),
      ASerializable(),
      _ncol(0),
      _nech(0),
      _isGrid(false),
      _array(),
      _attcol(),
      _colNames(),
      _p(),
      _grid(0)
{
  _initP();
  if (! coormin.empty()) ndim = (int) coormin.size();
  if (! coormax.empty()) ndim = MIN(ndim, (int) coormax.size());
  _ncol = ndim + flag_add_rank;
  _nech = nech;
  reset(_ncol, _nech);

  // Generate the sample number
  if (flag_add_rank) _createRank(0);

  // Generate the coordinates
  law_set_random_seed(seed);
  VectorDouble tab(ndim * nech);
  int ecr = 0;
  for (int idim = 0; idim < ndim; idim++)
  {
    double mini = (coormin.empty()) ? 0. : coormin[idim];
    double maxi = (coormax.empty()) ? 1. : coormax[idim];
    for (int iech = 0; iech < nech; iech++)
      tab[ecr++] = law_uniform(mini,maxi);
  }

  // Load the coordinates
  VectorString names = generateMultipleNames("x", ndim);
  _loadData(tab, names, VectorString(), ELoadBy::SAMPLE, flag_add_rank);

  int jcol = 0;
  if (flag_add_rank) jcol++;
  setLocatorsByAttribute(ndim, jcol, ELoc::X);
}

/**
 * Create a Db from a single sample whose coordinates are provided in 'tab'
 * @param tab Array containing the coordinates of the single sample
 * @param flag_add_rank 1 if the Sample ranks must be generated
 */
Db::Db(const VectorDouble& tab, int flag_add_rank)
    : AStringable(),
      ASerializable(),
      _ncol(0),
      _nech(0),
      _isGrid(false),
      _array(),
      _attcol(),
      _colNames(),
      _p(),
      _grid(0)
{
  _initP();

  int ndim = static_cast<int> (tab.size());
  _ncol = ndim + flag_add_rank;
  _nech = 1;
  reset(_ncol, _nech);

  // Generate the sample number
  if (flag_add_rank) _createRank(0);

  // Load the coordinates
  VectorString names = generateMultipleNames("x", ndim);
  _loadData(tab, names, VectorString(), ELoadBy::SAMPLE, flag_add_rank);

  int jcol = 0;
  if (flag_add_rank) jcol++;
  setLocatorsByAttribute(ndim, jcol, ELoc::X);
}

Db::Db(const Db& r)
    : _ncol(r._ncol),
      _nech(r._nech),
      _isGrid(r._isGrid),
      _array(r._array),
      _colNames(r._colNames),
      _p(r._p),
      _grid(r._grid)

{
}

Db& Db::operator=(const Db& r)
{
  if (this != &r)
  {
    _ncol = r._ncol;
    _nech = r._nech;
    _isGrid = r._isGrid;
    _array = r._array;
    _colNames = r._colNames;
    _p = r._p;
    _grid = r._grid;
  }
  return *this;
}

Db::~Db()
{
}

bool Db::isDimensionIndexValid(int idim) const
{
  if (idim < 0 || idim >= getNDim())
  {
    mesArg("Space Dimension", idim, getNDim());
    return false;
  }
  return true;
}

bool Db::isAttributeIndexValid(int iatt) const
{
  if (iatt < 0 || iatt >= getAttributeMaxNumber())
  {
    mesArg("Attribute Index", iatt, getAttributeMaxNumber());
    return false;
  }
  return true;
}

bool Db::isColumnIndexValid(int icol) const
{
  if (icol < 0 || icol >= _ncol)
  {
    mesArg("Column Index", icol, _ncol);
    return false;
  }
  return true;
}

bool Db::isSampleIndexValid(int iech) const
{
  if (iech < 0 || iech >= _nech)
  {
    mesArg("Sample Index", iech, _nech);
    return false;
  }
  return true;
}

bool Db::isLocatorIndexValid(const ELoc& locatorType, int locatorIndex) const
{
  if (!isLocatorTypeValid(locatorType)) return false;
  bool ok = _p.at(locatorType).isLocatorIndexValid(locatorIndex);
  if (! ok)
    messerr("Problem in the identification of Locator %d",locatorType);
  return ok;
}

int Db::getColumnByAttribute(int iatt) const
{
  if (!isAttributeIndexValid(iatt)) return -1;
  int icol = _attcol[iatt];
  return icol;
}

VectorInt Db::getColumnByAttribute(const VectorInt iatts) const
{
  VectorInt cols(iatts.size());
  for (unsigned int i = 0; i < iatts.size(); i++)
    cols[i] = getColumnByAttribute(iatts[i]);
  return cols;
}

int Db::_getAttributeByColumn(int icol) const
{
  if (!isColumnIndexValid(icol)) return -1;
  for (int iatt = 0; iatt < getAttributeMaxNumber(); iatt++)
    if (_attcol[iatt] == icol) return iatt;
  return -1;
}

int Db::getAttribute(const ELoc& locatorType, int locatorIndex) const
{
  if (!isLocatorIndexValid(locatorType, locatorIndex)) return -1;
  const PtrGeos& p = _p.at(locatorType);
  return p.getLocatorByIndex(locatorIndex);
}

/**
 * Find Column for a given Locator characteristics
 * @param locatorType Locator type
 * @param locatorIndex   Locator index (starting from 0)
 * @return
 */
int Db::getColumnByLocator(const ELoc& locatorType, int locatorIndex) const
{
  if (!isLocatorTypeValid(locatorType)) return -1;
  if (!isLocatorIndexValid(locatorType,locatorIndex)) return -1;
  const PtrGeos& p = _p.at(locatorType);
  int icol = getColumnByAttribute(p.getLocatorByIndex(locatorIndex));
  return (icol);
}

int Db::getLocatorNumber(const ELoc& locatorType) const
{
  if (!isLocatorTypeValid(locatorType)) return -1;
  const PtrGeos& p = _p.at(locatorType);
  return p.getLocatorNumber();
}

int Db::_findAttributeInLocator(const ELoc& locatorType, int iatt) const
{
  if (!isLocatorTypeValid(locatorType)) return -1;
  const PtrGeos& p = _p.at(locatorType);
  if (!isAttributeIndexValid(iatt)) return -1;
  for (int locatorIndex = 0; locatorIndex < p.getLocatorNumber(); locatorIndex++)
    if (p.getLocatorByIndex(locatorIndex) == iatt) return (locatorIndex);
  return -1;
}

int Db::_findColumnInLocator(const ELoc& locatorType, int icol) const
{
  if (!isLocatorTypeValid(locatorType)) return -1;
  int iatt = _getAttributeByColumn(icol);
  return _findAttributeInLocator(locatorType, iatt);
}

/**
 * Find the locator characteristics of a given Column
 * @param icol       Index of the target column
 * @param ret_locatorType Locator type
 * @param ret_locatorIndex Locator index (starting from 0)
 * @return
 */
int Db::getLocatorByColumn(int icol,
                           ELoc* ret_locatorType,
                           int* ret_locatorIndex) const
{
  auto it = ELoc::getIterator();
  while (it.hasNext())
  {
    if (*it != ELoc::UNKNOWN)
    {
      const PtrGeos& p = _p.at(*it);
      for (int i = 0; i < p.getLocatorNumber(); i++)
      {
        int jcol = getColumnByAttribute(p.getLocatorByIndex(i));
        if (icol == jcol)
        {
          *ret_locatorType = *it;
          *ret_locatorIndex = i;
          return true;
        }
      }
    }
    it.toNext();
  }
  *ret_locatorType = ELoc::UNKNOWN;
  *ret_locatorIndex = -1;
  return false;
}

int Db::getLocator(int iatt,
                   ELoc* ret_locatorType,
                   int* ret_locatorIndex) const
{
  if (!isAttributeIndexValid(iatt)) return -1;
  int icol = getColumnByAttribute(iatt);
  return getLocatorByColumn(icol, ret_locatorType, ret_locatorIndex);
}

/**
 * Return the locator information corresponding to the input variable
 * @param name Input variable name (unique)
 * @param ret_locatorType Locator Type
 * @param ret_locatorIndex Locator Index (starting from 0)
 * @return
 */
int Db::getLocator(const String& name,
                   ELoc *ret_locatorType,
                   int *ret_locatorIndex) const
{
  VectorInt iatts = ids(name, true);
  if (iatts.empty()) return -1;
  return getLocator(iatts[0], ret_locatorType, ret_locatorIndex);
}

VectorString Db::getLocators(bool anyLocator, const ELoc& locatorType) const
{
  VectorString retval;
  ELoc type;
  int item;

  for (int icol = 0; icol < _ncol; icol++)
  {
    if (! anyLocator)
    {
      (void) getLocatorByColumn(icol, &type, &item);
      if (type != locatorType) continue;
    }
    String string = _getLocatorNameByColumn(icol);
    retval.push_back(string);
  }
  return retval;
}

bool Db::isAttributeDefined(int iatt) const
{
  if (!isAttributeIndexValid(iatt)) return false;
  int icol = getColumnByAttribute(iatt);
  if (!isColumnIndexValid(icol)) return false;
  return (_attcol[icol] >= 0);
}

VectorString Db::expandNameList(const VectorString& names) const
{
  return expandList(_colNames, names);
}

VectorString Db::expandNameList(const String& names) const
{
  return expandList(_colNames, names);
}

VectorInt Db::ids(const String& name, bool flagOne) const
{
  VectorString exp_names = expandNameList(name);
  VectorInt iatts = getAttributesBasic(exp_names);
  if (! _isCountValid(iatts, flagOne)) return VectorInt();
  return iatts;
}

VectorInt Db::ids(const VectorString& names, bool flagOne) const
{
  VectorString exp_names = expandNameList(names);
  VectorInt iatts = getAttributesBasic(exp_names);
  if (! _isCountValid(iatts, flagOne)) return VectorInt();
  return iatts;
}

VectorInt Db::ids(const ELoc& locatorType, bool flagOne) const
{
  VectorString exp_names = getNames(locatorType);
  VectorInt iatts = getAttributesBasic(exp_names);
  if (! _isCountValid(iatts, flagOne)) return VectorInt();
  return iatts;
}

VectorInt Db::ids(const VectorInt& iatts, bool flagOne) const
{
  VectorString exp_names = getNames(iatts);
  if (! _isCountValid(iatts, flagOne)) return VectorInt();
  return iatts;
}

void Db::reset(int ncol, int nech)
{
  _ncol = ncol;
  _nech = nech;
  _isGrid = 0;

  /* The attribute pointers */

  _attcol.resize(ncol);
  for (int i = 0; i < ncol; i++)
    _attcol[i] = i;

  /* The variable names */

  _colNames = generateMultipleNames("New", ncol);

  /* The variable pointers */

  auto it = ELoc::getIterator();
  while (it.hasNext())
  {
    if (*it != ELoc::UNKNOWN)
      clearLocators(*it);
    it.toNext();
  }

  /* Main array */

  if (nech * ncol > 0) _array.resize(ncol * nech, 0);
}

/**
 * Set the value by Sample and Attribute
 * @param iech  Index of the Sample
 * @param iatt  Index of the Attribute
 * @param value Value to be assigned
 */
void Db::setArray(int iech, int iatt, double value)
{
  if (!isSampleIndexValid(iech)) return;
  int icol = getColumnByAttribute(iatt);
  if (!isColumnIndexValid(icol)) return;
  _array[_getAddress(iech, icol)] = value;
}

/**
 * Returns the value of the 'iech' sample of the variable 'name'
 *
 * This function does not use 'ids' mechanism in order to allow
 * referring to a non-existing variable
 */
double Db::getValue(const String& name, int iech) const
{
  int iatt = getAttribute(name);
  if (iatt < 0) return TEST;
  return getArray(iech, iatt);
}

/**
 * Sets the value of the 'iech' sample of the variable 'name'
 *
 * This function does not use 'ids' mechanism in order to allow
 * referring to a non-existing variable
 */
void Db::setValue(const String& name, int iech, double value)
{
  int iatt  = getAttribute(name);
  if (iatt < 0) return;
  setArray(iech, iatt, value);
}

/**
 * Return the value defined by Sample and Attribute
 * @param iech Sample Index
 * @param iatt Attribute Index
 * @return
 */
double Db::getArray(int iech, int iatt) const
{
  if (!isSampleIndexValid(iech)) return (TEST);
  int icol = getColumnByAttribute(iatt);
  if (!isColumnIndexValid(icol)) return (TEST);
  return (_array[_getAddress(iech, icol)]);
}

VectorDouble Db::getArray(int iatt, bool useSel) const
{
  int nech = getSampleNumber();
  VectorDouble sel, tab;
  if (!isAttributeIndexValid(iatt)) return tab;

  tab.resize(nech);
  if (useSel) sel = getSelection();

  int ecr = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    if (useSel && !sel.empty() && sel[iech] == 0) continue;
    tab[ecr] = getArray(iech, iatt);
    ecr++;
  }
  tab.resize(ecr);
  return tab;
}

void Db::updArray(int iech, int iatt, int oper, double value)
{
  if (!isSampleIndexValid(iech)) return;
  int icol = getColumnByAttribute(iatt);
  if (!isColumnIndexValid(icol)) return;
  double oldval = getArray(iech, icol);
  setArray(iech, icol, _updateValue(oper, oldval, value));
}

VectorDouble Db::getSampleCoordinates(int iech) const
{
  VectorDouble coor(getNDim());
  getSampleCoordinates(iech, coor);
  return coor;
}

void Db::getSampleCoordinates(int iech, VectorDouble& coor) const
{
  int ndim = getNDim();
  for (int idim = 0; idim < ndim; idim++)
    coor[idim] = getCoordinate(iech, idim);
}

/**
 * Return the coordinate of a sample along one Space Dimension
 * @param iech Rank of the sample
 * @param idim Rank of the Space Dimension
 * @param flag_rotate Use the rotation (only for Grid)
 * @return
 */
double Db::getCoordinate(int iech, int idim, bool flag_rotate) const
{
  if (idim >= getNDim()) return TEST;
  if (isGrid())
  {
    return _grid.getCoordinate(iech, idim, flag_rotate);
  }
  else
  {
    return getFromLocator(ELoc::X, iech, idim);
  }
}

void Db::getCoordinate(int iech, VectorDouble& coor, bool flag_rotate) const
{
  for (int idim = 0; idim < getNDim(); idim++)
    coor[idim] = getCoordinate(iech, idim, flag_rotate);
}

double Db::getDistance1D(int iech, int jech, int idim, bool flagAbs) const
{
  double v1 = getCoordinate(iech, idim);
  if (FFFF(v1)) return TEST;
  double v2 = getCoordinate(jech, idim);
  if (FFFF(v2)) return TEST;
  double delta = v1 - v2;
  if (flagAbs) delta = ABS(delta);
  return delta;
}

/**
 * Constitute a Vector of Vector of coordinates at a given sample, for all (active) samples
 * @param useSel
 * @return
 */
VectorVectorDouble Db::getCoordinates(bool useSel) const
{
  VectorVectorDouble result;
  VectorDouble local(getNDim());
  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    if (useSel && ! isActive(iech)) continue;
    getSampleCoordinates(iech, local);
    result.push_back(local);
  }
  return result;
}

void Db::setCoordinate(int iech, int idim, double value)
{
  if (!isSampleIndexValid(iech)) return;
  int icol = getColumnByAttribute(idim);
  if (!isColumnIndexValid(icol)) return;
  _array[_getAddress(iech, icol)] = value;
}

void Db::setFromLocator(const ELoc& locatorType,
                        int iech,
                        int locatorIndex,
                        double value)
{
  if (!isSampleIndexValid(iech)) return;
  int icol = getColumnByLocator(locatorType,locatorIndex);
  if (!isColumnIndexValid(icol)) return;
  _array[_getAddress(iech, icol)] = value;
}

double Db::getFromLocator(const ELoc& locatorType,
                          int iech,
                          int locatorIndex) const
{
  if (!isSampleIndexValid(iech)) return (TEST);
  int icol = getColumnByLocator(locatorType, locatorIndex);
  if (!isColumnIndexValid(icol)) return (TEST);
  return (_array[_getAddress(iech, icol)]);
}

int Db::getFromLocatorNumber(const ELoc& locatorType) const
{
  const PtrGeos& p = _p.at(locatorType);
  return p.getLocatorNumber();
}

void Db::_initP(void)
{
  _p.clear();
  auto it = ELoc::getIterator();
  while (it.hasNext())
  {
    if (*it != ELoc::UNKNOWN)
      _p[*it].resize(0);
    it.toNext();
  }
}

int Db::_getAddress(int iech, int icol) const
{
  return ((iech) + _nech * icol);
}

void Db::printLocators(void) const
{
  /* Loop on the pointers */

  mestitle(1, "List of locators");
  int rank = 0;
  auto it = ELoc::getIterator();
  while (it.hasNext())
  {
    if (*it != ELoc::UNKNOWN)
    {
      const PtrGeos& p = _p.at(*it);
      if (p.getLocatorNumber() > 0)
      {
        p.print(rank, *it);
        message("- Columns    = ");
        for (int locatorIndex = 0; locatorIndex < p.getLocatorNumber(); locatorIndex++)
          message("%2d ", getColumnByAttribute(p.getLocatorByIndex(locatorIndex)));
        message("\n");
        rank++;
      }
    }
    it.toNext();
  }
}

void Db::printAttributes(void) const
{
  mestitle(1, "List of attributes");
  message("Maximum number of positions = %d\n", getAttributeMaxNumber());
  message("Number of Columns           = %d\n", getFieldNumber());

  /* Loop on the attributes */

  if (getAttributeMaxNumber() <= 0) return;
  message("Attribute = ");
  for (int iatt = 0; iatt < getAttributeMaxNumber(); iatt++)
    message("%2d ", _attcol[iatt]);
  message("\n");
  return;
}

void Db::clearLocators(const ELoc& locatorType)
{
  PtrGeos& p = _p[locatorType];
  p.clear();
}

/**
 * Setting the locator for a set of variables designated by their names
 * @param names        Vector if variable names
 * @param locatorType  Locator type (include ELoc::UNKNOWN)
 * @param locatorIndex Starting locator rank (starting from 0)
 * @return
 */
void Db::setLocator(const VectorString& names,
                    const ELoc& locatorType,
                    int locatorIndex)
{
  if (!isLocatorTypeValid(locatorType, true)) return;
  VectorInt iatts = ids(names, false);
  if (iatts.empty()) return;
  for (unsigned int i = 0; i < iatts.size(); i++)
    setLocatorByAttribute(iatts[i], locatorType, locatorIndex + i);
}

/**
 * Define the Locator(s) for the given variable(s)
 * @param names Set of variable names
 * @param locatorType Locator Type
 * @param locatorIndex Locator Index (for the first variable) (starting from 0)
 */
void Db::setLocator(const String& names, const ELoc& locatorType, int locatorIndex)
{
  if (!isLocatorTypeValid(locatorType, true)) return;
  VectorInt iatts = ids(names, false);
  if (iatts.empty()) return;
  for (unsigned int i = 0; i < iatts.size(); i++)
    setLocatorByAttribute(iatts[i], locatorType, locatorIndex + i);
}

/**
 * Setting the locator for a variable designated by its attribute
 * @param iatt          Index of the Attribute
 * @param locatorType   Type of locator (include ELoc::UNKNOWN)
 * @param locatorIndex  Rank in the Locator (starting from 0)
 * @return Error return code
 * @remark: At this stage, no check is performed to see if items
 * @remark: are consecutive and all defined
 * @remark: This allow using this function in any order.
 */
void Db::setLocatorByAttribute(int iatt,
                               const ELoc& locatorType,
                               int locatorIndex)
{
  if (!isAttributeIndexValid(iatt)) return;
  if (!isLocatorTypeValid(locatorType, true)) return;
  if (locatorIndex < 0) return;

  /* Cancel any locator referring to this column */

  auto it = ELoc::getIterator();
  while (it.hasNext())
  {
    if (*it != ELoc::UNKNOWN)
    {
      PtrGeos& p = _p[*it];
      int found = _findAttributeInLocator(*it, iatt);
      if (found >= 0)
        p.erase(found);
    }
    it.toNext();
  }

  // Check if this locator already exists for the current pointer type

  if (locatorType != ELoc::UNKNOWN)
  {
    PtrGeos& p = _p[locatorType];
    int nitem = p.getLocatorNumber();
    if (locatorIndex >= nitem) p.resize(locatorIndex + 1);
    p.setLocatorByIndex(locatorIndex, iatt);
  }
}

String Db::_getLocatorNameByColumn(int icol) const
{
  ELoc locatorType;
  int locatorIndex;
  (void) getLocatorByColumn(icol, &locatorType, &locatorIndex);
  return getLocatorName(locatorType, locatorIndex);
}

/**
 * Set the Locators for a set of variables identified by their attribute rank
 * @param number        Number of variables to be set
 * @param iatt          Rank of the first attribute
 * @param locatorType   Type of the Locator (include ELoc::UNKNOWN)
 * @param locatorIndex  Rank of the first Locator index (starting from 0)
 */
void Db::setLocatorsByAttribute(int number,
                                int iatt,
                                const ELoc& locatorType,
                                int locatorIndex)
{
  if (!isLocatorTypeValid(locatorType, true)) return;
  for (int i = 0; i < number; i++)
    setLocatorByAttribute(iatt+i, locatorType, locatorIndex + i);
}

/**
 * Create a set of new variables in an already existing Db and initialize
 * their contents to a constant value
 * @param nadd     Number of variables to be added
 * @param valinit  Value to be used for variable initialization
 * @param radix    Generic radix given to the newly created variables
 * @param locatorType Generic locator assigned to new variables
 * @param locatorIndex   Locator index (starting from 0)
 * @param nechInit Number of samples (used only if the Db is initially empty)
 * @return Rank of the first attribute
 */
int Db::addFields(int nadd,
                  double valinit,
                  const String& radix,
                  const ELoc& locatorType,
                  int locatorIndex,
                  int nechInit)
{
  int ncol = _ncol;
  int nmax = getAttributeMaxNumber();
  int nnew = ncol + nadd;
  if (nadd <= 0) return (-1);

  /* Case of an empty Db, define the number of samples using 'nechInit' */

  if (_nech <= 0) _nech = nechInit;

  /* Dimension the array */

  _array.resize(_nech * nnew);

  /* Dimension the attribute pointer */

  _attcol.resize(nmax + nadd);
  for (int i = 0; i < nadd; i++)
    _attcol[nmax + i] = ncol + i;

  // Set the name
  _colNames.resize(nnew);
  if (nadd == 1)
    _colNames[ncol] = radix;
  else
  {
    VectorString names = generateMultipleNames(radix,nadd);
    for (int i = 0; i < nadd; i++)
      _colNames[ncol + i] = names[i];
  }
  (void) correctNamesForDuplicates(_colNames);

  // Initialize the variables with a given value
  _columnInit(nadd, ncol, valinit);

  // Set the locator (if defined)
  if (locatorType != ELoc::UNKNOWN)
    setLocatorsByAttribute(nadd, nmax, locatorType, locatorIndex);

  _ncol += nadd;

  return (nmax);
}

/**
 * Add one or several fields to an already existing Db. This is performed
 * by providing an array of values 'tab'. Its dimension must be equal to the
 * number of samples (or active samples if 'useSel' is true, times the number
 * of variables 'nvar'.
 * @param tab    Array to be loaded
 * @param radix  Generic name for the newly created variables
 * @param locatorType Generic locator assigned to new variables
 * @param locatorIndex   Locator index (starting from 0)
 * @param useSel true if values for active samples are provided
 * @param valinit initial value (for unselected samples)
 * @param nvar   Number of variables loaded
 * @return Rank of the first attribute
 */
int Db::addFields(const VectorDouble& tab,
                  const String& radix,
                  const ELoc& locatorType,
                  int locatorIndex,
                  bool useSel,
                  double valinit,
                  int nvar)
{
  // Particular case where the Db is empty.
  // Set its dimension to the number of samples of the input array 'tab'
  if (_nech <= 0) _nech = static_cast<int> (tab.size()) / nvar;

  // Check dimensions
  int nech = (useSel) ? getActiveSampleNumber() : getSampleNumber();
  if ((int) tab.size() != nvar * nech)
  {
    messerr("Db::addFields : Incompatibility between dimension of 'tab' (%d)", tab.size());
    messerr("and 'nvar'(%d) * 'nech'(%d)", nvar, nech);
    return 1;
  }

  // Adding the new Fields
  int iatt = addFields(nvar, valinit, radix, locatorType, locatorIndex);
  if (iatt < 0) return 1;

  setFieldByAttribute(tab, iatt, useSel);

  const double* local = tab.data();
  for (int ivar = 0; ivar < nvar; ivar++)
    setFieldByAttribute(&local[ivar * nech], iatt + ivar, useSel);

  return iatt;
}

void Db::setColumnByRank(const double* tab, int icol, bool useSel)
{
  if (!isColumnIndexValid(icol)) return;
  VectorDouble sel;

  if (useSel) sel = getSelection();

  int lec = 0;
  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    bool defined = true;
    if (useSel && !sel.empty()) defined = (sel[iech] == 1);

    double value = TEST;
    if (defined)
      value = tab[lec++];
    else
    {
      value = TEST;
      if (!useSel) lec++;
    }
    setByColumn(iech, icol, value);
  }
}

void Db::setColumnByRank(const VectorDouble& tab, int icol, bool useSel)
{
  setColumnByRank(tab.data(), icol, useSel);
}

void Db::setFieldByAttribute(const double* tab, int iatt, bool useSel)
{
  if (!isAttributeIndexValid(iatt)) return;
  VectorDouble sel;

  if (useSel) sel = getSelection();

  int lec = 0;
  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    bool defined = true;
    if (useSel && !sel.empty()) defined = (sel[iech] == 1);

    double value = TEST;
    if (defined)
      value = tab[lec++];
    else
    {
      value = TEST;
      if (!useSel) lec++;
    }
    setArray(iech, iatt, value);
  }
}

void Db::setFieldByAttribute(const VectorDouble& tab, int iatt, bool useSel)
{
  setFieldByAttribute(tab.data(), iatt, useSel);
}

void Db::setField(const VectorDouble& tab, const String& name, bool useSel)
{
  VectorInt iatts = ids(name, true);
  if (iatts.empty()) return;
  setFieldByAttribute(tab.data(), iatts[0], useSel);
}

void Db::duplicateColumnByAttribute(int iatt_in, int iatt_out)
{
  if (!isAttributeIndexValid(iatt_in)) return;
  if (!isAttributeIndexValid(iatt_out)) return;
  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    double value = getArray(iech, iatt_in);
    setArray(iech, iatt_out, value);
  }
}

void Db::deleteField(const String& names)
{
  VectorInt iatts = ids(names, false);
  if (iatts.empty()) return;

  for (unsigned int i = 0; i < iatts.size(); i++)
    deleteFieldByAttribute(iatts[i]);
}

void Db::deleteField(const VectorString& names)
{
  VectorInt iatts = ids(names, false);
  if (iatts.empty()) return;

  for (unsigned int i = 0; i < iatts.size(); i++)
    deleteFieldByAttribute(iatts[i]);
}

int Db::addSelection(const VectorDouble& tab, const String& name)
{
  int iatt = addFields(tab, name, ELoc::SEL);
  return iatt;
}

/**
 * Create a selection around the only defined values of the target variable
 * @param testvar Name of the target variable
 * @param limits  Limits defining the Definition Domain to be tested (optional)
 * @param name    Name of the newly created selection
 * @return
 */
int Db::addSelection(const String& testvar,
                     const Limits& limits,
                     const String& name)
{
  int iatt = addFields(1,0.,name,ELoc::SEL);
  if (iatt < 0) return 1;

  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    double value = getValue(testvar,iech);
    double answer = 1;
    if (FFFF(value))
    {
      answer = 0;
    }
    else if (! limits.empty())
    {
      if (! limits.isInside(value)) answer = 0;
    }
    setArray(iech, iatt, answer);
  }
  return 0;
}

int Db::addSamples(int nadd, double valinit)
{
  int nech = _nech;
  int nnew = nech + nadd;
  if (nadd <= 0) return (-1);

  /* Core allocation */

  VectorDouble new_array(_ncol * nnew);
  for (int i = 0; i < _ncol * nnew; i++)
    new_array[i] = valinit;

  /* Copy the array */

  for (int icol = 0; icol < _ncol; icol++)
    for (int iech = 0; iech < nech; iech++)
    {
      int iad1 = iech + nnew * icol;
      new_array[iad1] = _array[_getAddress(iech, icol)];
    }

  /* Core deallocation */

  _array = new_array;
  _nech = nnew;
  return (nech);
}

void Db::deleteSample(int e_del)
{
  int nech = _nech;
  int nnew = nech - 1;
  if (!isSampleIndexValid(e_del)) return;

  /* Core allocation */

  VectorDouble new_array(_ncol * nnew);

  /* Copy the array */

  for (int icol = 0; icol < _ncol; icol++)
    for (int iech = 0; iech < nech; iech++)
    {
      if (iech == e_del) continue;
      int jech = (iech < e_del) ? iech :
                                  iech - 1;
      int iad1 = jech + nnew * icol;
      new_array[iad1] = _array[_getAddress(iech, icol)];
    }

  /* Core deallocation */

  _array = new_array;
  _nech = nnew;
}

void Db::deleteFieldByRank(int rank_del)
{
  deleteFieldByAttribute(rank_del - 1);
}

/**
 * Delete an Attribute
 * @param iatt_del Rank of the attribute to be deleted
 */
void Db::deleteFieldByAttribute(int iatt_del)
{
  int ncol = _ncol;
  int nech = _nech;
  int nmax = getAttributeMaxNumber();
  int nnew = ncol - 1;
  if (!isAttributeIndexValid(iatt_del)) return;

  /* Identify the column to be deleted */

  int c_del = getColumnByAttribute(iatt_del);
  if (!isColumnIndexValid(c_del)) return;
  _attcol[iatt_del] = -1;
  for (int iatt = 0; iatt < nmax; iatt++)
  {
    if (_attcol[iatt] < c_del) continue;
    _attcol[iatt]--;
  }

  /* Dimension the array */

  for (int icol = c_del + 1; icol < ncol; icol++)
    for (int iech = 0; iech < nech; iech++)
      _array[_getAddress(iech, icol - 1)] = _array[_getAddress(iech, icol)];
  _array.resize(nech * nnew);

  /* Resize the variable pointers */

  auto it = ELoc::getIterator();
  while (it.hasNext())
  {
    if (*it != ELoc::UNKNOWN)
    {
      PtrGeos& p = _p[*it];
      int found = _findAttributeInLocator(*it, iatt_del);
      if (found >= 0) p.erase(found);
    }
    it.toNext();
  }

  /* Resize the variables names */

  _colNames.erase(_colNames.begin() + c_del);

  /* Set the error return code */

  _ncol = nnew;
}

void Db::deleteFieldByLocator(const ELoc& locatorType)
{
  if (!isLocatorTypeValid(locatorType)) return;
  PtrGeos& p = _p[locatorType];
  int nitem = p.getLocatorNumber();
  // Loop is performed downwards as PtrGeos is modified by called routine
  for (int locatorIndex = nitem - 1; locatorIndex >= 0; locatorIndex--)
  {
    deleteFieldByAttribute(p.getLocatorByIndex(locatorIndex));
  }
}

/**
 * Returns a Vector containing the minimum and maximum along a Space dimension
 * @param idim   Rank of the Space dimension
 * @param useSel true if the possible selection must be taken intao account
 * @return
 */
VectorDouble Db::getExtrema(int idim, bool useSel) const
{
  VectorDouble ext;
  if (!isDimensionIndexValid(idim)) return ext;
  VectorDouble coor = getCoordinate(idim, useSel);
  ext.push_back(ut_vector_min(coor));
  ext.push_back(ut_vector_max(coor));
  return ext;
}

double Db::getExtension(int idim, bool useSel) const
{
  if (!isDimensionIndexValid(idim)) return 0.;
  VectorDouble coor = getCoordinate(idim, useSel);
  double mini = ut_vector_min(coor);
  double maxi = ut_vector_max(coor);
  return maxi - mini;
}

double Db::getFieldSize(bool useSel) const
{
  double diag = 0.;
  for (int idim = 0; idim < getNDim(); idim++)
  {
    VectorDouble ext = getExtrema(idim);
    double delta = ext[1] - ext[0];
    diag += delta * delta;
  }
  return sqrt(diag);
}

double Db::getMinimum(const String& name, bool useSel) const
{
  VectorInt iatts = ids(name, true);
  if (iatts.empty()) return TEST;
  VectorDouble tab = getFieldByAttribute(iatts[0], useSel);
  return ut_vector_min(tab);
}

double Db::getMaximum(const String& name, bool useSel) const
{
  VectorInt iatts = ids(name, true);
  if (iatts.empty()) return TEST;
  VectorDouble tab = getFieldByAttribute(iatts[0], useSel);
  return ut_vector_max(tab);
}

double Db::getMean(const String& name, bool useSel) const
{
  VectorInt iatts = ids(name, true);
  if (iatts.empty()) return TEST;
  VectorDouble tab = getFieldByAttribute(iatts[0], useSel);
  return ut_vector_mean(tab);
}

double Db::getVariance(const String& name, bool useSel) const
{
  VectorInt iatts = ids(name, true);
  if (iatts.empty()) return TEST;
  VectorDouble tab = getFieldByAttribute(iatts[0], useSel);
  return ut_vector_var(tab);
}

double Db::getStdv(const String& name, bool useSel) const
{
  VectorInt iatts = ids(name, true);
  if (iatts.empty()) return TEST;
  VectorDouble tab = getFieldByAttribute(iatts[0], useSel);
  return ut_vector_stdv(tab);
}

int Db::gridDefine(const VectorInt& nx,
                   const VectorDouble& dx,
                   const VectorDouble& x0,
                   const VectorDouble& angles)
{
  _isGrid = 1;
  return (_grid.init(nx, dx, x0, angles));
}

void Db::gridCopyParams(int mode, const GridC& gridaux)
{
  _grid.copyParams(mode, gridaux);
}

bool Db::isSameGrid(const GridC& grid) const
{
  if (! isGrid() || grid.empty())
  {
    messerr("Both files should be organized as grids");
    return false;
  }
  return _grid.isSame(grid);
}

bool Db::isSameGridMesh(const Db& dbaux) const
{
  if (! isGrid() || ! dbaux.isGrid())
  {
    messerr("Both files should be organized as grids");
    return false;
  }
  return _grid.isSameMesh(dbaux.getGrid());
}

bool Db::isSameGridRotation(const Db& dbaux) const
{
  if (! isGrid() || ! dbaux.isGrid())
  {
    messerr("Both files should be organized as grids");
    return false;
  }
  if (! isGridRotated() && ! dbaux.isGridRotated()) return true;
  return _grid.isSameRotation(dbaux.getGrid());
}

bool Db::isGridRotated() const
{
  return (_isGrid && _grid.isRotated());
}

int Db::getNDim() const
{
  if (isGrid())
  {
    return (_grid.getNDim());
  }
  else
  {
    return (_p.at(ELoc::X).getLocatorNumber());
  }
}

bool Db::hasSameDimension(const Db* dbaux, bool verbose) const
{
  bool retOK = dbaux->getNDim() == getNDim();
  if (!retOK)
    messerr("The two Data bases should have the same Space Dimension");
  return retOK;
}

/**
 * Check if the Space Dimension of 'dbaux' is larger (or equal) than the one of 'this'
 * @param dbaux    Second Db
 * @param verbose  Verbose flag
 * @return
 */
bool Db::hasLargerDimension(const Db* dbaux, bool verbose) const
{
  bool retOK = dbaux->getNDim() >= getNDim();
  if (!retOK)
  {
    messerr("The Space Dimension of the Secondary Data base (%d)",
            dbaux->getNDim());
    messerr(
        "should be larger than the Space Dimension of the Current Data Base (%d)",
        getNDim());
  }
  return retOK;
}

int Db::getNX(int idim) const
{
  if (!isGrid()) return (-1);
  if (!isDimensionIndexValid(idim)) return (1);
  return (_grid.getNX(idim));
}

double Db::getDX(int idim) const
{
  if (!isGrid()) return (-1);
  if (!isDimensionIndexValid(idim)) return (TEST);
  return (_grid.getDX(idim));
}

double Db::getX0(int idim) const
{
  if (!isGrid()) return (-1);
  if (!isDimensionIndexValid(idim)) return (0.);
  return (_grid.getX0(idim));
}

double Db::getAngles(int idim) const
{
  if (!isGrid()) return (-1);
  if (!isDimensionIndexValid(idim)) return (0.);
  return _grid.getRotAngles(idim);
}

void Db::_columnInit(int ncol, int icol0, double valinit)
{
  for (int jcol = 0; jcol < ncol; jcol++)
  {
    int icol = jcol + icol0;

    if (! GlobalEnvironment::getEnv()->isDomainReference() ||
        getFromLocatorNumber(ELoc::DOM) == 0)
      for (int iech = 0; iech < _nech; iech++)
        _array[_getAddress(iech, icol)] = valinit;
    else
      for (int iech = 0; iech < _nech; iech++)
      {
        double value = getFromLocator(ELoc::DOM, iech, 0);
        int iad = _getAddress(iech, icol);
        if (GlobalEnvironment::getEnv()->matchDomainReference(value))
          _array[iad] = valinit;
        else
          _array[iad] = TEST;
      }
  }
  return;
}

void Db::switchLocator(const ELoc& locatorType_in, const ELoc& locatorType_out)
{
  PtrGeos& p_in  = _p[locatorType_in];
  PtrGeos& p_out = _p[locatorType_out];
  int n_in  = getFromLocatorNumber(locatorType_in);
  int n_out = getFromLocatorNumber(locatorType_out);

  /* Move the gradient components into additional variables */
  p_out.resize(n_in + n_out);
  for (int i_in = 0; i_in < n_in; i_in++)
    p_out.setLocatorByIndex(n_out + i_in, p_in.getLocatorByIndex(i_in));
  p_in.clear();
}

double Db::getByColumn(int iech, int icol) const
{
  if (!isColumnIndexValid(icol)) return TEST;
  return (_array[_getAddress(iech, icol)]);
}

void Db::setByColumn(int iech, int icol, double value)
{
  if (!isColumnIndexValid(icol)) return;
  if (!isSampleIndexValid(iech)) return;
  _array[_getAddress(iech, icol)] = value;
}

int Db::getVariableNumber() const
{
  return getFromLocatorNumber(ELoc::Z);
}

/**
 * Checks the number of variables in 'this' compared to the required 'nvar'
 * - compare=0: they should be equal
 * - compare<0: 'this' should contain less (or equal) than 'nvar'
 * - compare>0: 'this' should contain more (or equal) than 'nvar'
 */
bool Db::isVariableNumberComparedTo(int nvar, int compare) const
{
  if (compare == 0)
  {
    if (! (getVariableNumber() == nvar))
    {
      messerr("This function requires %d variables (locator 'Z'). The 'Db' contains %d variables",
              nvar,getVariableNumber());
      return false;
    }
  }
  else if (compare < 0)
  {
    if (! (getVariableNumber() <= nvar))
    {
      messerr("This function requires nvar <= %d variables (locator 'Z'). The 'Db' contains %d variables",
              nvar,getVariableNumber());
      return false;
    }
  }
  else
  {
    if (! (getVariableNumber() > nvar))
    {
      messerr("This function requires nvar >= %d variables (locator 'Z'). The 'Db' contains %d variables",
              nvar,getVariableNumber());
      return false;
    }
  }
  return true;
}

bool Db::hasVariable() const
{
  return (getVariableNumber() > 0);
}

double Db::getVariable(int iech, int item) const
{
  if (!hasVariable()) return (TEST);
  return getFromLocator(ELoc::Z, iech, item);
}

void Db::setVariable(int iech, int item, double value)
{
  setFromLocator(ELoc::Z, iech, item, value);
}

/****************************************************************************/
/*!
 **  Update the value for a variable
 **
 ** \param[in]  iech    Rank of the sample
 ** \param[in]  item    Rank of the variable
 ** \param[in]  oper    Type of operation
 ** \li                 0 : New = New + Old
 ** \li                 1 : New = New * Old
 ** \li                 2 : New = New - Old
 ** \li                 3 : New = Old / New
 ** \li                 4 : New = New (only if old is defined)
 ** \li                 5 : New = MAX(New, Old)
 ** \li                 6 : New = MIN(New, Old)
 ** \param[in]  value   Update value
 **
 ** \remark  For efficiency reason, argument validity is not tested
 **
 *****************************************************************************/
void Db::updVariable(int iech, int item, int oper, double value)
{
  double oldval = getFromLocator(ELoc::Z, iech, item);
  double newval = _updateValue(oper, oldval, value);
  setFromLocator(ELoc::Z, iech, item, newval);
}

int Db::getIntervalNumber() const
{
  return MAX(getLowerIntervalNumber(), getUpperIntervalNumber());
}

int Db::getLowerIntervalNumber() const
{
  return getFromLocatorNumber(ELoc::RKLOW);
}

bool Db::hasLowerInterval() const
{
  return (getLowerIntervalNumber() > 0);
}

int Db::getUpperIntervalNumber() const
{
  return getFromLocatorNumber(ELoc::RKUP);
}

bool Db::hasUpperInterval() const
{
  return (getUpperIntervalNumber() > 0);
}

double Db::getLowerInterval(int iech, int item) const
{
  if (!hasLowerInterval()) return TEST;
  return getFromLocator(ELoc::RKLOW, iech, item);
}

double Db::getUpperInterval(int iech, int item) const
{
  if (!hasUpperInterval()) return TEST;
  return getFromLocator(ELoc::RKUP, iech, item);
}

void Db::setLowerInterval(int iech, int item, double rklow)
{
  setFromLocator(ELoc::RKLOW, iech, item, rklow);
}

void Db::setUpperInterval(int iech, int item, double rkup)
{
  setFromLocator(ELoc::RKUP, iech, item, rkup);
}

void Db::setIntervals(int iech, int item, double rklow, double rkup)
{
  if (rklow > rkup)
  {
    messerr("Setting Intervals: Lower (%lf) cannot be larger than upper (%lf)",
            rklow,rkup);
    return;
  }
  setFromLocator(ELoc::RKLOW, iech, item, rklow);
  setFromLocator(ELoc::RKUP,  iech, item, rkup);
}

int Db::getLowerBoundNumber() const
{
  return getFromLocatorNumber(ELoc::L);
}

int Db::getUpperBoundNumber() const
{
  return getFromLocatorNumber(ELoc::U);
}

bool Db::hasLowerBound() const
{
  return (getLowerBoundNumber() > 0);
}

bool Db::hasUpperBound() const
{
  return (getUpperBoundNumber() > 0);
}

double Db::getLowerBound(int iech, int item) const
{
  if (!hasLowerBound()) return TEST;
  return getFromLocator(ELoc::L, iech, item);
}

double Db::getUpperBound(int iech, int item) const
{
  if (!hasUpperBound()) return TEST;
  return getFromLocator(ELoc::U, iech, item);
}

void Db::setLowerBound(int iech, int item, double lower)
{
  setFromLocator(ELoc::L, iech, item, lower);
}

void Db::setUpperBound(int iech, int item, double upper)
{
  setFromLocator(ELoc::U, iech, item, upper);
}
void Db::setBounds(int iech, int item, double lower, double upper)
{
  if (lower > upper)
  {
    messerr("Setting bounds: Lower (%lf) cannot be larger than upper (%lf)",
            lower,upper);
    return;
  }
  setLowerBound(iech,item,lower);
  setUpperBound(iech,item,upper);
}

VectorDouble Db::getWithinBounds(int item, bool useSel) const
{
  int nech = (useSel) ? getActiveSampleNumber() : getSampleNumber();
  VectorDouble vec(nech);
  VectorDouble vecl = getFieldByLocator(ELoc::L, item, useSel);
  VectorDouble vecu = getFieldByLocator(ELoc::U, item, useSel);

  for (int iech = 0; iech < nech; iech++)
  {
    double vall = (vecl.empty()) ? TEST : vecl[iech];
    double valu = (vecu.empty()) ? TEST : vecu[iech];
    if (FFFF(vall))
    {
      if (FFFF(valu))
        vec[iech] = TEST;
      else
        vec[iech] = valu;
    }
    else
    {
      if (FFFF(valu))
        vec[iech] = vall;
      else
        vec[iech] = (vall + valu) / 2.;
    }
  }
  return vec;
}

int Db::getGradientNumber() const
{
  return getFromLocatorNumber(ELoc::G);
}

bool Db::hasGradient() const
{
  return (getGradientNumber() > 0);
}

double Db::getGradient(int iech, int item) const
{
  if (!hasGradient()) return TEST;
  return getFromLocator(ELoc::G, iech, item);
}

void Db::setGradient(int iech, int item, double value)
{
  setFromLocator(ELoc::G, iech, item, value);
}

int Db::getTangentNumber() const
{
  return getFromLocatorNumber(ELoc::TGTE);
}

bool Db::hasTangent() const
{
  return (getTangentNumber() > 0);
}

double Db::getTangent(int iech, int item) const
{
  if (!hasTangent()) return TEST;
  return getFromLocator(ELoc::TGTE, iech, item);
}

void Db::setTangent(int iech, int item, double value)
{
  setFromLocator(ELoc::TGTE, iech, item, value);
}

int Db::getProportionNumber() const
{
  return getFromLocatorNumber(ELoc::P);
}

bool Db::hasProportion() const
{
  return (getProportionNumber() > 0);
}

double Db::getProportion(int iech, int item) const
{
  if (!hasProportion()) return TEST;
  return getFromLocator(ELoc::P, iech, item);
}

void Db::setProportion(int iech, int item, double value)
{
  setFromLocator(ELoc::P, iech, item, value);
}

int Db::getSelection(int iech) const
{
  if (!hasSelection()) return 1;
  double value = getFromLocator(ELoc::SEL, iech, 0);
  if (FFFF(value)) return 1;
  int sel = (value != 0) ? 1 :
                           0;
  return (sel);
}

void Db::setSelection(int iech, int value)
{
  setFromLocator(ELoc::SEL, iech, 0, (value == 0) ? 0. :
                                                  1.);
}

bool Db::hasSelection() const
{
  return (getFromLocatorNumber(ELoc::SEL) > 0);
}

int Db::getActiveSampleNumber() const
{
  if (!hasSelection()) return (getSampleNumber());

  /* Case when a selection is present */

  int count = 0;
  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    double value = getFromLocator(ELoc::SEL, iech, 0);
    if (value != 0) count++;
  }
  return count;
}

double Db::getWeight(int iech) const
{
  if (!hasWeight()) return 1.;
  double w = getFromLocator(ELoc::W, iech, 0);
  if (FFFF(w)) w = 1.;
  if (w < 0) w = 0.;
  return (w);
}

VectorDouble Db::getWeight(bool useSel) const
{
  int icol = -1;
  int nech = getSampleNumber();
  VectorDouble sel;
  VectorDouble tab(nech);

  if (useSel) sel = getSelection();
  if (hasWeight()) icol = getColumnByLocator(ELoc::W, 0);

  int ecr = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    if (useSel && !sel.empty() && sel[iech] == 0) continue;
    if (icol >= 0)
      tab[ecr] = getByColumn(iech, icol);
    else
      tab[ecr] = 1.;
    ecr++;
  }
  tab.resize(ecr);
  return tab;
}

void Db::setWeight(int iech, double value)
{
  if (value < 0) value = 0.;
  setFromLocator(ELoc::W, iech, 0, value);
}

bool Db::hasWeight() const
{
  return (getFromLocatorNumber(ELoc::W) > 0);
}

int Db::getExternalDriftNumber() const
{
  return getFromLocatorNumber(ELoc::F);
}

bool Db::hasExternalDrift() const
{
  return (getExternalDriftNumber() > 0);
}

double Db::getExternalDrift(int iech, int item) const
{
  if (!hasExternalDrift()) return TEST;
  return getFromLocator(ELoc::F, iech, item);
}

void Db::setExternalDrift(int iech, int item, double value)
{
  setFromLocator(ELoc::F, iech, item, value);
}

int Db::getBlockExtensionNumber() const
{
  return getFromLocatorNumber(ELoc::BLEX);
}

bool Db::hasBlockExtension() const
{
  return (getBlockExtensionNumber() > 0);
}

double Db::getBlockExtension(int iech, int item) const
{
  if (!hasBlockExtension()) return TEST;
  return getFromLocator(ELoc::BLEX, iech, item);
}

void Db::setBlockExtension(int iech, int item, double value)
{
  setFromLocator(ELoc::BLEX, iech, item, value);
}

int Db::getCodeNumber() const
{
  return getFromLocatorNumber(ELoc::C);
}

bool Db::hasCode() const
{
  return (getCodeNumber() > 0);
}

double Db::getCode(int iech) const
{
  if (!hasCode()) return 0;
  double code = getFromLocator(ELoc::C, iech, 0);
  if (FFFF(code)) code = 0.;
  return code;
}

void Db::setCode(int iech, double value)
{
  setFromLocator(ELoc::C, iech, 0, value);
}

/****************************************************************************/
/*!
 **  Returns the list of Unique codes
 **
 ** \return  Pointer to the array containing a single occurence of each code
 **
 *****************************************************************************/
VectorDouble Db::getCodeList(void)
{
  int nred;
  VectorDouble tab(_nech);

  /* Load all the codes */

  int number = 0;
  for (int iech = 0; iech < _nech; iech++)
  {
    if (isActive(iech)) tab[number++] = getCode(iech);
  }

  /* Get the Unique occurrence */

  ut_tab_unique(number, tab.data(), &nred);
  if (nred < number) tab.resize(nred);
  return (tab);
}

int Db::getVarianceErrorNumber() const
{
  return getFromLocatorNumber(ELoc::V);
}

bool Db::hasVarianceError() const
{
  return (getVarianceErrorNumber() > 0);
}

double Db::getVarianceError(int iech, int item) const
{
  if (!hasVarianceError()) return 0.;
  return getFromLocator(ELoc::V, iech, item);
}

void Db::setVarianceError(int iech, int item, double value)
{
  setFromLocator(ELoc::V, iech, item, value);
}

bool Db::hasDomain() const
{
  return (getFromLocatorNumber(ELoc::DOM) > 0);
}

/****************************************************************************/
/*!
 **  Reads the domain flag of a sample
 **
 ** \return  1 if the Domain variable is not defined
 ** \return    or if the Domain variable is defined and equal to the
 ** \return    the Reference Domain value
 **
 ** \param[in]  iech Rank of the sample
 **
 *****************************************************************************/
int Db::getDomain(int iech) const
{
  if (! GlobalEnvironment::getEnv()->isDomainReference()) return 1;
  if (! hasDomain()) return 1;
  double value = getFromLocator(ELoc::DOM, iech, 0);
  if (FFFF(value)) return (0);
  if (! GlobalEnvironment::getEnv()->matchDomainReference(value)) return 1;
  return 0;
}

void Db::setDomain(int iech, int value)
{
  setFromLocator(ELoc::DOM, iech, 0, (value == 0) ? 0. : 1.);
}

int Db::getDipDirectionNumber() const
{
  return getFromLocatorNumber(ELoc::ADIR);
}

bool Db::hasDipDirection() const
{
  return (getDipDirectionNumber() > 0);
}

double Db::getDipDirection(int iech) const
{
  if (!hasDipDirection()) return TEST;
  return getFromLocator(ELoc::ADIR, iech, 0);
}

void Db::setDipDirection(int iech, double value)
{
  setFromLocator(ELoc::ADIR, iech, 0, value);
}

int Db::getDipAngleNumber() const
{
  return getFromLocatorNumber(ELoc::ADIP);
}

bool Db::hasDipAngle() const
{
  return (getDipDirectionNumber() > 0);
}

double Db::getDipAngle(int iech) const
{
  return getFromLocator(ELoc::ADIP, iech, 0);
}

void Db::setDipAngle(int iech, double value)
{
  setFromLocator(ELoc::ADIP, iech, 0, value);
}

int Db::getObjectSizeNumber() const
{
  return getFromLocatorNumber(ELoc::SIZE);
}

bool Db::hasObjectSize() const
{
  return (getObjectSizeNumber() > 0);
}

double Db::getObjectSize(int iech) const
{
  if (!hasObjectSize()) return TEST;
  return getFromLocator(ELoc::SIZE, iech, 0);
}

void Db::setObjectSize(int iech, double value)
{
  setFromLocator(ELoc::SIZE, iech, 0, value);
}

int Db::getBorderUpNumber() const
{
  return getFromLocatorNumber(ELoc::BU);
}

bool Db::hasBorderUp() const
{
  return (getBorderUpNumber() > 0);
}

double Db::getBorderUp(int iech) const
{
  if (!hasBorderUp()) return TEST;
  return getFromLocator(ELoc::BU, iech, 0);
}

void Db::setBorderUp(int iech, double value)
{
  setFromLocator(ELoc::BU, iech, 0, value);
}

int Db::getBorderDownNumber() const
{
  return getFromLocatorNumber(ELoc::BD);
}

bool Db::hasBorderDown() const
{
  return (getBorderDownNumber() > 0);
}

double Db::getBorderDown(int iech) const
{
  if (!hasBorderDown()) return TEST;
  return getFromLocator(ELoc::BD, iech, 0);
}

void Db::setBorderDown(int iech, double value)
{
  setFromLocator(ELoc::BD, iech, 0, value);
}

int Db::getDateNumber() const
{
  return getFromLocatorNumber(ELoc::DATE);
}

bool Db::hasDate() const
{
  return (getDateNumber() > 0);
}

double Db::getDate(int iech) const
{
  if (!hasDate()) return 0.;
  return getFromLocator(ELoc::DATE, iech, 0);
}

void Db::setDate(int iech, double value)
{
  setFromLocator(ELoc::DATE, iech, 0, value);
}

int Db::getSimvarRank(int isimu, int ivar, int icase, int nbsimu, int nvar)
{
  return (_getSimrank(isimu, ivar, icase, nbsimu, nvar));
}

double Db::getSimvar(const ELoc& locatorType,
                     int iech,
                     int isimu,
                     int ivar,
                     int icase,
                     int nbsimu,
                     int nvar) const
{
  int item = _getSimrank(isimu, ivar, icase, nbsimu, nvar);

  return getFromLocator(locatorType, iech, item);
}

void Db::setSimvar(const ELoc& locatorType,
                   int iech,
                   int isimu,
                   int ivar,
                   int icase,
                   int nbsimu,
                   int nvar,
                   double value)
{
  int item = _getSimrank(isimu, ivar, icase, nbsimu, nvar);
  setFromLocator(locatorType, iech, item, value);
}

void Db::updSimvar(const ELoc& locatorType,
                   int iech,
                   int isimu,
                   int ivar,
                   int icase,
                   int nbsimu,
                   int nvar,
                   int oper,
                   double value)
{
  int item = _getSimrank(isimu, ivar, icase, nbsimu, nvar);
  double oldval = getFromLocator(locatorType, iech, item);
  double newval = _updateValue(oper, oldval, value);
  setFromLocator(locatorType, iech, item, newval);
}

bool Db::isActive(int iech) const
{
  return (getSelection(iech) && getDomain(iech));
}

bool Db::isActiveAndDefined(int iech, int item) const
{
  if (!isActive(iech)) return false;;
  if (FFFF(getVariable(iech, item))) return false;
  return true;
}

/**
 * Returns the number of active samples for which the target variable (ELoc::Z)
 * is defined
 * @param item Rank of the ELoc::Z variable
 * @return Number of samples
 */
int Db::getActiveAndDefinedNumber(int item) const
{
  int nech = 0;
  for (int iech = 0; iech < _nech; iech++)
  {
    if (!isActive(iech)) continue;
    if (FFFF(getVariable(iech, item))) continue;
    nech++;
  }
  return (nech);
}

/**
 * Returns the number of active samples for which the variable 'name'
 * is defined
 * @param name Name of the Target variable
 * @return Number of samples
 */
int Db::getActiveAndDefinedNumber(const String& name) const
{
  VectorInt iatts = ids(name, true);
  if (iatts.empty()) return 0;
  VectorDouble tab = getFieldByAttribute(iatts[0], false);

  int nech = 0;
  for (int iech = 0; iech < (int) tab.size(); iech++)
  {
    if (! FFFF(tab[iech])) nech++;
  }
  return (nech);
}


/**
 * Update an Old by a New value according to 'oper'
 * @param oper   0 : New = New + old
 *               1 : New = New * Old
 *               2 : New = New - Old
 *               3 : New = Old / New
 *               4 : New = New (only if 'Old' is defined)
 *               5 : New = MAX(New, Old)
 *               6 : New = MIN(New, Old)
 * @param oldval Old value
 * @param value  New value
 */
double Db::_updateValue(int oper, double oldval, double value)
{
  double newval;

  newval = 0.;
  switch (oper)
  {
    case 0:
      if (FFFF(value) || FFFF(oldval)) return (TEST);
      newval = value + oldval;
      break;

    case 1:
      if (FFFF(value) || FFFF(oldval)) return (TEST);
      newval = value * oldval;
      break;

    case 2:
      if (FFFF(value) || FFFF(oldval)) return (TEST);
      newval = value - oldval;
      break;

    case 3:
      if (FFFF(value) || FFFF(oldval)) return (TEST);
      newval = (value == 0.) ? TEST :
                               oldval / value;
      break;

    case 4:
      newval = value;
      break;

    case 5:
      if (FFFF(value)) return (oldval);
      if (FFFF(oldval)) return (value);
      newval = MAX(newval, value);
      break;

    case 6:
      if (FFFF(value)) return (oldval);
      if (FFFF(oldval)) return (value);
      newval = MIN(newval, value);
      break;
  }
  return (newval);
}

/**
 * Returns the rank of (one of) the lastly added attribute in the Db
 * @param number 0 designates the last, 1 the one before last...
 * @return
 */
int Db::getLastAttribute(int number) const
{
  VectorInt ranks;
  for (int i = 0; i < (int) _attcol.size(); i++)
    if (_attcol[i] >= 0) ranks.push_back(i);
  int size = static_cast<int> (ranks.size());
  if (number > size)
    return -1;
  else
    return ranks[size - number - 1];
}

String Db::getLastName(int number) const
{
  int iatt = getLastAttribute(number);
  String name = getName(iatt);
  return name;
}

int Db::_getLastColumn(int number) const
{
  if (number > _ncol)
    return -1;
  else
    return (_ncol - number);
}

String Db::getName(const ELoc& locatorType, int locatorIndex) const
{
  int icol = getColumnByLocator(locatorType, locatorIndex);
  if (icol < 0) return String();
  return _colNames[icol];
}

String Db::getName(int iatt) const
{
  int icol = getColumnByAttribute(iatt);
  if (icol < 0) return ("");
  return getNameByColumn(icol);
}

VectorString Db::getNames(const ELoc& locatorType) const
{
  VectorString namelist;
  if (!isLocatorTypeValid(locatorType)) return namelist;
  int count = getFromLocatorNumber(locatorType);
  for (int i = 0; i < count; i++)
  {
    int icol = getColumnByLocator(locatorType, i);
    namelist.push_back(getNameByColumn(icol));
  }
  return namelist;
}

VectorString Db::getNames(const VectorInt& iatts) const
{
  VectorString namelist;
  int count = static_cast<int> (iatts.size());
  for (int i = 0; i < count; i++)
  {
    int icol = getColumnByAttribute(iatts[i]);
    namelist.push_back(getNameByColumn(icol));
  }
  return namelist;
}

VectorString Db::getNames(const String& name) const
{
  return expandNameList(name);
}

VectorString Db::getNames(const VectorString& names) const
{
  return expandNameList(names);
}

VectorString Db::getNames() const
{
  VectorString names = _colNames;
  return names;
}

void Db::_setNameByColumn(int icol, const String& name)
{
  if (!isColumnIndexValid(icol)) return;
  _colNames[icol] = name;
  correctNewNameForDuplicates(_colNames, icol);
}

void Db::setName(int iatt, const String& name)
{
  int icol = getColumnByAttribute(iatt);
  if (icol < 0) return;
  _colNames[icol] = name;
  correctNewNameForDuplicates(_colNames, icol);
}

void Db::setName(const String& old_name, const String& name)
{
  int icol = getColumn(old_name);
  if (icol < 0) return;
  _colNames[icol] = name;
  correctNewNameForDuplicates(_colNames, icol);
}

void Db::setName(const VectorString list, const String& name)
{
  for (int i = 0; i < (int) list.size(); i++)
  {
    int icol = getColumn(list[i]);
    if (icol < 0) continue;
    _colNames[icol] = incrementStringVersion(name, i + 1);
  }
  correctNamesForDuplicates(_colNames);
}

void Db::setName(const ELoc& locatorType, const String& name)
{
  VectorString namelist;
  if (!isLocatorTypeValid(locatorType)) return;
  int count = getFromLocatorNumber(locatorType);
  for (int i = 0; i < count; i++)
  {
    int icol = getColumnByLocator(locatorType, i);
    if (icol < 0) continue;
    _colNames[icol] = incrementStringVersion(name, i+ 1);
  }
  correctNamesForDuplicates(_colNames);
  return;
}

String Db::_summaryString(void) const
{
  std::stringstream sstr;

  sstr << toTitle(1, "Data Base Summary");

  if (isGrid())
    sstr << "File is organized as a regular grid" << std::endl;
  else
    sstr << "File is organized as a set of isolated points" << std::endl;

  sstr << "Space dimension              = " << getNDim() << std::endl;
  sstr << "Number of fields             = " << getFieldNumber() << std::endl;
  sstr << "Maximum Number of attributes = " << getAttributeMaxNumber()
       << std::endl;
  sstr << "Total number of samples      = " << getSampleNumber() << std::endl;
  if (hasSelection())
    sstr << "Number of active samples     = " << getActiveSampleNumber()
         << std::endl;

  if (isGrid()) sstr << _grid.toString(0);

  return sstr.str();
}

String Db::_summaryExtensionString(void) const
{
  std::stringstream sstr;
  int ndim = getNDim();
  if (ndim <= 0) return sstr.str();

  /* Printout */

  sstr << toTitle(1, "Data Base Extension");
  for (int idim = 0; idim < ndim; idim++)
  {
    VectorDouble coor = getCoordinate(idim, true);
    double vmin = ut_vector_min(coor);
    double vmax = ut_vector_max(coor);

    sstr << "Coor #" << idim + 1;
    sstr << " - Min = " << toDouble(vmin);
    sstr << " - Max = " << toDouble(vmax);
    sstr << " - Ext = " << vmax - vmin;
    sstr << std::endl;
  }

  return sstr.str();
}

String Db::_summaryVariableString(void) const
{
  std::stringstream sstr;

  if (getFieldNumber() <= 0) return sstr.str();
  sstr << toTitle(1, "Variables");

  for (int icol = 0; icol < getFieldNumber(); icol++)
  {
    sstr << "Field = " << icol + 1;
    sstr << " - Name = " << getNameByColumn(icol);
    sstr << " - Locator = " << _getLocatorNameByColumn(icol);
    sstr << std::endl;
  }
  return sstr.str();
}

/**
 * Print statistics on the variable
 * @param cols List of Columns of target variable (all if empty)
 * @param mode 1 for basic statistics; 2 for clas statistics
 * @param maxNClass Maximum number of printed classes
 * @return
 */
String Db::_summaryVariableStat(VectorInt cols, int mode, int maxNClass) const
{
  std::stringstream sstr;

  sstr << toTitle(1, "Data Base Statistics");

  int nval, nmask, ntest, nout;
  double vmin, vmax, delta, mean, stdv;
  int nech = getSampleNumber();
  VectorDouble tab, wgt;

  // Loop on the columns

  int ncol = (cols.empty()) ? getFieldNumber() : static_cast<int> (cols.size());
  for (int jcol = 0; jcol < ncol; jcol++)
  {
    int icol = (cols.empty()) ? jcol :
                                cols[jcol];
    if (!isColumnIndexValid(icol)) continue;

    tab = getColumnByRank(icol, true);
    wgt = getWeight(true);

    ut_statistics(static_cast<int> (tab.size()), tab.data(), NULL,
                  wgt.data(), &nval, &vmin, &vmax,
                  &delta, &mean, &stdv);

    sstr << icol + 1 << " - Name " << getNameByColumn(icol) << " - Locator "
         << _getLocatorNameByColumn(icol) << std::endl;
    sstr << " Nb of data          = " << toInt(nech) << std::endl;
    sstr << " Nb of active values = " << toInt(nval) << std::endl;
    if (nval <= 0) continue;

    /* Dispatch */

    if (mode == 1)
    {
      sstr << " Minimum value       = " << toDouble(vmin) << std::endl;
      sstr << " Maximum value       = " << toDouble(vmax) << std::endl;
      sstr << " Mean value          = " << toDouble(mean) << std::endl;
      sstr << " Standard Deviation  = " << toDouble(stdv) << std::endl;
      sstr << " Variance            = " << toDouble(stdv * stdv) << std::endl;
    }
    else
    {
      vmin = floor(vmin - 0.5);
      vmax = ceil(vmax + 0.5);
      int nclass = (int) (vmax - vmin) + 1;
      if (nclass > maxNClass)
        sstr << " Number of classes is truncated to " << maxNClass << std::endl;
      nclass = MIN(maxNClass, nclass);
      VectorInt classe(nclass);
      ut_classify(static_cast<int> (tab.size()), tab.data(), NULL,
                  nclass, vmin, 1., &nmask,
                  &ntest, &nout, classe.data());
      if (ntest > 0)
        sstr << " Unknown values      = " << toInt(ntest) << std::endl;
      if (nout > 0)
        sstr << " Outside classes     = " << toInt(nout) << std::endl;

      for (int iclass = 0; iclass < nclass; iclass++)
      {
        if (classe[iclass] <= 0) continue;
        sstr << " Class" << toInt((int) vmin + iclass);
        sstr << " = " << toInt(classe[iclass]);
        sstr << " (" << toDouble(100. * classe[iclass] / nval) << "%)";
        sstr << std::endl;
      }
    }
  }
  return sstr.str();
}

String Db::_summaryArrayString(VectorInt cols, bool flagSel) const
{
  std::stringstream sstr;

  sstr << toTitle(1, "Data Base Contents");

  int ncol = (cols.empty()) ? getFieldNumber() : static_cast<int> (cols.size());
  int number = (flagSel) ? getActiveSampleNumber() : getSampleNumber();

  VectorDouble tab;
  VectorString colnames;
  for (int jcol = 0; jcol < ncol; jcol++)
  {
    int icol = (cols.empty()) ? jcol :
                                cols[jcol];
    if (!isColumnIndexValid(icol)) continue;
    VectorDouble local = getColumnByRank(icol, flagSel);
    tab.insert(tab.end(), local.begin(), local.end());
    colnames.push_back(getNameByColumn(icol));
  }

  sstr << toMatrix(String(), colnames, VectorString(), 1, ncol, number, tab);

  return sstr.str();
}

void Db::displayMore(unsigned char params,
                     const VectorInt& cols,
                     bool flagSel,
                     int mode) const
{
  messageFlush(_display(params, cols, flagSel, mode));
}

void Db::displayMore(unsigned char params,
                     const VectorString& names,
                     bool flagSel,
                     int mode) const
{
  VectorInt iatts = ids(names, false);
  if (iatts.empty()) return;
  VectorInt cols = getColumnByAttribute(iatts);
  messageFlush(_display(params, cols, flagSel, mode));
}

/**
 * Print the characteristics of the Db structure
 * @param params      Mask specifying the printing options
 * @param cols        Optional list of selected variables
 * @param flagSel     Take the selection into account
 * @param mode        1 for basic statistics; 2 for class statistics
 * @return
 */
String Db::_display(unsigned char params,
                    const VectorInt& cols,
                    bool flagSel,
                    int mode) const
{
  static int MAX_NCLASS = 50;
  std::stringstream sstr;

  sstr << toTitle(0, "Data Base Characteristics");

  /* Print the Summary of the Db */

  if (params & FLAG_RESUME) sstr << _summaryString();

  /* Print the Extension */

  if (params & FLAG_EXTEND) sstr << _summaryExtensionString();

  /* Print the statistics */

  if (params & FLAG_STATS) sstr << _summaryVariableStat(cols, mode, MAX_NCLASS);

  /* Print the contents of the Data Base */

  if (params & FLAG_ARRAY) sstr << _summaryArrayString(cols, flagSel);

  /* Print the list of variables */

  if (params & FLAG_VARS) sstr << _summaryVariableString();

  return sstr.str();
}

String Db::toString(int level) const
{
  std::stringstream sstr;
  sstr << _display(FLAG_RESUME | FLAG_VARS);
  return sstr.str();
}

VectorDouble Db::getSelection(void) const
{
  int nech = getSampleNumber();
  VectorDouble tab;

  if (!hasSelection()) return tab;
  int icol = getColumnByLocator(ELoc::SEL,0);
  if (!isColumnIndexValid(icol)) return tab;

  tab.resize(nech);
  for (int iech = 0; iech < nech; iech++)
    tab[iech] = getByColumn(iech, icol);
  return tab;
}

VectorDouble Db::getColumnByRank(int icol, bool useSel) const
{
  int nech = getSampleNumber();
  VectorDouble tab, sel;
  if (!isColumnIndexValid(icol)) return tab;

  tab.resize(nech, TEST);
  if (useSel) sel = getSelection();

  int ecr = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    bool defined = true;
    if (useSel && !sel.empty()) defined = (sel[iech] == 1);
    if (!defined) continue;
    tab[ecr] = getByColumn(iech, icol);
    ecr++;
  }
  tab.resize(ecr);
  return tab;
}

VectorDouble Db::getFieldByAttribute(int iatt, bool useSel) const
{
  int icol = getColumnByAttribute(iatt);
  return getColumnByRank(icol, useSel);
}

VectorDouble Db::getFieldByLocator(const ELoc& locatorType,
                                   int locatorIndex,
                                   bool useSel) const
{
  int icol = getColumnByLocator(locatorType, locatorIndex);
  if (icol < 0) return VectorDouble();
  return getColumnByRank(icol, useSel);
}

VectorDouble Db::getField(const String& name, bool useSel) const
{
  VectorInt iatts = ids(name, true);
  if (iatts.empty()) return VectorDouble();
  int icol = getColumnByAttribute(iatts[0]);
  if (icol < 0) return VectorDouble();
  return getColumnByRank(icol, useSel);
}

VectorDouble Db::getFieldsByLocator(const ELoc& locatorType, bool useSel) const
{
  VectorString names = getNames(locatorType);
  return getFields(names, useSel);
}

VectorDouble Db::getFieldsByAttribute(const VectorInt& iatts, bool useSel) const
{
  int nech = (useSel) ? getActiveSampleNumber() :
                        getSampleNumber();
  int nvar = static_cast<int> (iatts.size());
  VectorDouble retval(nvar * nech);

  /* Loop on the variables to be retrieved */

  int ecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    VectorDouble local = getFieldByAttribute(iatts[ivar], useSel);
    for (int iech = 0; iech < nech; iech++)
      retval[ecr++] = local[iech];
  }
  return retval;
}

VectorDouble Db::getColumnsByRank(const VectorInt& icols, bool useSel) const
{
  int nech = getSampleNumber();
  int nvar = static_cast<int> (icols.size());
  VectorDouble retval(nvar * nech);

  /* Loop on the variables to be retrieved */

  int ecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    VectorDouble local = getColumnByRank(icols[ivar], useSel);
    for (int iech = 0; iech < nech; iech++)
      retval[ecr++] = local[iech];
  }
  return retval;
}

VectorDouble Db::getColumnsByRank(int icol_beg, int icol_end, bool useSel) const
{
  VectorInt icols;
  for (int icol = icol_beg; icol < icol_end; icol++)
    icols.push_back(icol);
  return getColumnsByRank(icols, useSel);
}

VectorDouble Db::getFieldsByAttribute(int iatt_beg,
                                      int iatt_end,
                                      bool useSel) const
{
  VectorInt iatts;
  for (int iatt = iatt_beg; iatt < iatt_end; iatt++)
    iatts.push_back(iatt);
  return getFieldsByAttribute(iatts, useSel);
}

VectorDouble Db::getFields(const VectorString& names, bool useSel) const
{
  VectorInt iatts;
  if (names.empty())
    iatts = getAttributes();
  else
    iatts = ids(names, false);
  return getFieldsByAttribute(iatts, useSel);
}

/**
 * Returns the vector of coordinates along a given Space Dimension
 * @param idim    Rank of the Space dimension
 * @param useSel  Use the Data Selection
 * @param flag_rotate Flag for rotation (only for Grid)
 * @return
 */
VectorDouble Db::getCoordinate(int idim, bool useSel, bool flag_rotate) const
{
  int nech = getSampleNumber();
  VectorDouble tab, sel;

  tab.resize(nech, TEST);
  if (useSel) sel = getSelection();

  int ecr = 0;
  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    bool defined = true;
    if (useSel && !sel.empty()) defined = (sel[iech] == 1);
    if (!defined) continue;
    tab[ecr] = getCoordinate(iech, idim, flag_rotate);
    ecr++;
  }

  tab.resize(ecr);
  return tab;
}

/**
 * Returns the rank of the Single field corresponding to 'name'
 * @param name Named for the searched field
 * @return The rank of the Single field or -1
 */
int Db::getColumn(const String& name) const
{
  VectorString exp_name = expandNameList(name);
  if (exp_name.empty()) return -1;
  return getRankInList(_colNames, exp_name[0]);
}

VectorInt Db::getColumns(const VectorString& names) const
{
  VectorString exp_names = expandNameList(names);
  if (exp_names.size() <= 0) return VectorInt();
  int number = static_cast<int> (exp_names.size());
  VectorInt icols(number);
  for (int i = 0; i < number; i++)
    icols[i] = getColumn(exp_names[i]);
  return icols;
}

VectorInt Db::getColumnsByAttribute(const ELoc& locatorType) const
{
  VectorInt icols;
  if (!isLocatorTypeValid(locatorType)) return icols;
  int number = getLocatorNumber(locatorType);
  if (number <= 0) return icols;

  icols.resize(number);
  for (int i = 0; i < number; i++)
    icols[i] = getColumnByLocator(locatorType, i);
  return icols;
}

/**
 * Returns the Single Attribute which corresponds to the searched name
 * @param name Name to be searched for
 * @return Rank of the Attribute or -1
 */
int Db::getAttribute(const String& name) const
{
  VectorInt iatts = ids(name, true);
  if (iatts.empty()) return -1;
  int icol = getColumnByAttribute(iatts[0]);
  return _getAttributeByColumn(icol);
}

/**
 * This is a BASIC function returning the vector of ranks of the Attribute
 * which corresponds to a set of existing names
 */
VectorInt Db::getAttributesBasic(const VectorString& names) const
{
  if (names.empty()) return VectorInt();

  VectorInt iatts(names.size());
  for (unsigned int i = 0; i < names.size(); i++)
  {
    int icol = getRankInList(_colNames, names[i]);
    iatts[i] = _getAttributeByColumn(icol);
  }
  return iatts;
}

VectorInt Db::getAttributes(const VectorString& names) const
{
  if (names.empty()) return VectorInt();

  VectorInt iatts(names.size());
  for (unsigned int i = 0; i < names.size(); i++)
    iatts[i] = getAttribute(names[i]);
  return iatts;
}

VectorInt Db::getAttributes(const ELoc& locatorType) const
{
  VectorInt iatts;
  if (!isLocatorTypeValid(locatorType)) return iatts;
  int number = getLocatorNumber(locatorType);
  if (number <= 0) return iatts;
  iatts.resize(number);
  for (int i = 0; i < number; i++)
    iatts[i] = getAttribute(locatorType, i);
  return iatts;
}

VectorInt Db::getAttributes() const
{
  VectorInt iatts;
  for (int i = 0; i < (int) _attcol.size(); i++)
    if (_attcol[i] >= 0) iatts.push_back(i);
  return iatts;
}

void Db::_loadData(const VectorDouble& tab,
                   const VectorString& names,
                   const VectorString& locatorNames,
                   const ELoadBy& order,
                   int shift)
{
  // Preliminary check

  if (_ncol <= 0) return;
  if (tab.empty()) return;
  if (!isMultiple(static_cast<int> (tab.size()), _nech))
  {
    messerr("The Dimension of the array (%d) is inconsistent", tab.size());
    messerr("It should be a multiple of the number of samples (%d)", _nech);
    return;
  }
  int ntab = static_cast<int> (tab.size()) / _nech;
  int ecr = 0;
  for (int icol = 0; icol < ntab; icol++)
  {
    int jcol = icol + shift;
    for (int iech = 0; iech < _nech; iech++, ecr++)
    {
      if (order == ELoadBy::SAMPLE)
        setByColumn(iech, jcol, tab[icol + ntab * iech]);
      else
        setByColumn(iech, jcol, tab[ecr]);
    }
  }

  // Set the names
  _defineDefaultNames(shift, names);

  // Set the locators
  _defineDefaultLocators(shift, locatorNames);

  return;
}

void Db::_createRank(int shift)
{
  int nech = getSampleNumber();
  for (int iech = 0; iech < nech; iech++)
    setArray(iech, shift, iech + 1);

  // Set the name

  _setNameByColumn(shift, "rank");

  // Set the locators (No particular action for the Rank)
}

void Db::_createGridCoordinates(int shift)
{
  // Set the Names

  for (int idim = 0; idim < getNDim(); idim++)
    _setNameByColumn(shift + idim, getLocatorName(ELoc::X, idim));

  // Set the locators

  setLocatorsByAttribute(getNDim(), shift, ELoc::X);

  // Generate the vector of coordinates

  _grid.iteratorInit();
  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    VectorInt indices = _grid.iteratorNext();
    VectorDouble coors = _grid.indiceToCoordinate(indices);
    for (int idim = 0; idim < getNDim(); idim++)
      setCoordinate(iech, shift + idim, coors[idim]);
  }
}

void Db::_defineDefaultNames(int shift, const VectorString& names)
{
  int ncol = getFieldNumber() - shift;
  if (!names.empty())
  {
    if ((int) names.size() != ncol) throw("Error in the dimension of 'names'");
  }

  for (int icol = 0; icol < ncol; icol++)
  {
    if (!names.empty())
      _setNameByColumn(icol + shift, names[icol]);
    else
      _setNameByColumn(icol + shift, incrementStringVersion("New", icol + 1));
  }
}

void Db::_defineDefaultLocators(int shift, const VectorString& locatorNames)
{
  if (locatorNames.empty()) return;

  int ncol = getFieldNumber() - shift;
  if ((int) locatorNames.size() != ncol)
    throw("Error in the dimension of 'locatorNames'");

  ELoc locatorType;
  int locatorIndex, mult;
  for (int icol = 0; icol < ncol; icol++)
  {
    if (!locatorIdentify(locatorNames[icol], &locatorType, &locatorIndex, &mult))
      setLocatorByAttribute(icol + shift, locatorType, locatorIndex);
  }
}

void Db::_defineDefaultLocatorsByNames(int shift, const VectorString& names)
{
  if (names.empty()) return;

  int ncol = getFieldNumber() - shift;
  if ((int) names.size() != ncol) throw("Error in the dimension of 'names'");

  ELoc locatorType;
  int locatorIndex, mult;
  for (int icol = 0; icol < ncol; icol++)
  {
    if (!locatorIdentify(names[icol], &locatorType, &locatorIndex, &mult))
      setLocatorByAttribute(icol + shift, locatorType, locatorIndex);
  }
}

VectorDouble Db::_statistics(const VectorInt& iatts,
                             const VectorString& opers,
                             bool flagIso,
                             bool flagVariableWise,
                             bool flagPrint,
                             double proba,
                             double vmin,
                             double vmax,
                             const String& title,
                             NamingConvention namconv)
{
  VectorDouble stats;

  int natt = static_cast<int> (iatts.size());
  if (natt <= 0) return stats;

  VectorInt iopers = statsList(opers);
  int noper = static_cast<int> (iopers.size());
  if (noper <= 0) return stats;

  // Add the variables for PointWise statistics
  if (!flagVariableWise)
  {
    int iattn = addFields(noper);
    if (iattn < 0) return VectorDouble();

    dbStatisticsVariables(this, iatts, iopers, iattn, vmin, vmax, proba);

    namconv.setNamesAndLocators(this, iattn);
    for (int i = 0; i < noper; i++)
      namconv.setNamesAndLocators(this, iattn + i, opers[i]);
    return VectorDouble();
  }
  else
  {
    stats = dbStatisticsMono(this, iatts, iopers, flagIso, proba, vmin, vmax);

    if (flagPrint)
    {
      VectorString varnames = getNames(iatts);
      messageFlush(statisticsMonoPrint(stats, iopers, varnames, title));
      return VectorDouble();
    }
  }
  return stats;
}

VectorDouble Db::statistics(const VectorString& names,
                            const VectorString& opers,
                            bool flagIso,
                            bool flagVariableWise,
                            bool flagPrint,
                            double proba,
                            double vmin,
                            double vmax,
                            const String& title,
                            NamingConvention namconv)
{
  VectorInt iatts = ids(names, false);
  if (iatts.empty()) return VectorDouble();
  return _statistics(iatts, opers, flagIso, flagVariableWise, flagPrint, proba,
                     vmin, vmax, title, namconv);
}

VectorDouble Db::_statisticsMulti(const VectorInt& iatts,
                                  bool flagIso,
                                  bool flagPrint,
                                  const String& title)
{
  VectorDouble stats;

  int natt = static_cast<int> (iatts.size());
  if (natt <= 0) return stats;

  stats = dbStatisticsMulti(this, iatts, flagIso);

  if (flagPrint)
  {
    VectorString varnames = getNames(iatts);
    messageFlush(statisticsMultiPrint(stats, varnames, title));
  }
  return stats;
}

VectorDouble Db::statisticsMulti(const VectorString& names,
                                 bool flagIso,
                                 bool flagPrint,
                                 const String& title)
{
  VectorInt iatts = ids(names, false);
  if (iatts.empty()) return VectorDouble();

  return _statisticsMulti(iatts, flagIso, flagPrint, title);
}

/****************************************************************************/
/*!
 **  Returns the address of the combination
 **
 ** \return  Returned value
 **
 ** \param[in]  isimu     Rank of the simulation
 ** \param[in]  ivar      Rank of the variable
 ** \param[in]  icase     Rank of the GRF / PGS
 ** \param[in]  nbsimu    Number of simulations
 ** \param[in]  nvar      Number of variables
 **
 *****************************************************************************/
int Db::_getSimrank(int isimu, int ivar, int icase, int nbsimu, int nvar) const
{
  return (isimu + nbsimu * (ivar + nvar * icase));
}

int Db::deSerialize(const String& filename, bool verbose)
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

  /* Opening the Data file */

  if (_fileOpen(filename, "Db", "r", verbose)) return 1;

  /* Check the grid organization */

  if (_recordRead("Grid Flag", "%d", &flag_grid)) goto label_end;

  /* Grid case; read the grid header */

  if (flag_grid)
  {

    /* Decoding the header */

    if (_recordRead("Space Dimension", "%d", &ndim)) goto label_end;

    /* Core allocation */

    nx.resize(ndim);
    dx.resize(ndim);
    x0.resize(ndim);
    angles.resize(ndim);

    /* Read the grid characteristics */

    for (int idim = 0; idim < ndim; idim++)
    {
      if (_recordRead("Grid Number of Nodes", "%d", &nx[idim])) goto label_end;
      if (_recordRead("Grid Origin", "%lf", &x0[idim])) goto label_end;
      if (_recordRead("Grid Mesh", "%lf", &dx[idim])) goto label_end;
      if (_recordRead("Grid Angles", "%lf", &angles[idim])) goto label_end;
    }
    ntot = ut_ivector_prod(nx);
  }

  /* Reading the tail of the file */

  _variableRead(&natt, &ndim2, &nech, tabatt, tabnum, tabnam, tab);

  /* Creating the Db */

  if (flag_grid)
  {
    if (natt > 0 && nech != ntot)
    {
      messerr("The number of lines read from the Grid file (%d)", nech);
      messerr("is not a multiple of the number of samples (%d)", ntot);
      messerr("The Grid Db is created with no sample attached");
      natt = 0;
    }
    reset(natt + flag_add_rank, ut_ivector_prod(nx));
    (void) gridDefine(nx, dx, x0, angles);
    _loadData(ELoadBy::SAMPLE, flag_add_rank, tab);
  }
  else
  {
    reset(natt + flag_add_rank, nech);
    _loadData(ELoadBy::SAMPLE, flag_add_rank, tab);
  }

  /* Loading the names */

  if (natt > 0)
    for (i = 0; i < natt; i++)
      setName(i + flag_add_rank, tabnam[i]);

  /* Create the locators */

  if (natt > 0) for (i = 0; i < natt; i++)
    setLocatorByAttribute(i + flag_add_rank, tabatt[i], tabnum[i]);

  /* Core deallocation */

  label_end:
  _fileClose(verbose);
  return 0;
}

int Db::serialize(const String& filename, bool verbose) const
{
  bool onlyLocator = false;
  bool flag_grid = isGrid();

  /* Opening the Data file */

  if (_fileOpen(filename, "Db", "w", verbose)) return 1;

  /* Writing the file organization */

  _recordWrite("%d", flag_grid);
  _recordWrite("#", "File organization (0:Points; 1:Grid)");

  if (flag_grid)
  {

    /* Writing the header */

    _recordWrite("%d", getNDim());
    _recordWrite("#", "Space Dimension");

    /* Writing the grid characteristics */

    _recordWrite("#", "Grid characteristics (NX,X0,DX,ANGLE)");
    for (int idim = 0; idim < getNDim(); idim++)
    {
      _recordWrite("%d",  getNX(idim));
      _recordWrite("%lf", getX0(idim));
      _recordWrite("%lf", getDX(idim));
      _recordWrite("%lf", getAngles(idim));
      _recordWrite("\n");
    }
  }

  /* Writing the tail of the file */

  if (_variableWrite(flag_grid, onlyLocator)) return 1;

  // Close the Neutral file

  _fileClose(verbose);

  return 0;
}

int Db::_variableWrite(bool flag_grid, bool onlyLocator) const
{
  int ecr, item, rankZ;
  ELoc locatorType = ELoc::UNKNOWN;

  /* Preliminary check */

  if (getFieldNumber() <= 0 || getSampleNumber() <= 0) return 0;

  /* Count the number of variables to be written */

  int ncol = 0;
  for (int icol = 0; icol < getFieldNumber(); icol++)
  {
    if (!getLocatorByColumn(icol, &locatorType, &item))
    {
      if (onlyLocator) continue;
      locatorType = ELoc::Z;
    }
    if (flag_grid && locatorType == ELoc::X) continue;
    ncol++;
  }
  _recordWrite("%d", ncol);
  _recordWrite("#", "Number of variables");

  /* Print the locators */

  _recordWrite("#", "Locators");
  rankZ = getLocatorNumber(ELoc::Z);
  ecr = 0;
  for (int icol =  0; icol < getFieldNumber(); icol++)
  {
    if (! getLocatorByColumn(icol, &locatorType, &item))
    {
      if (onlyLocator) continue;
      locatorType = ELoc::Z;
      item = rankZ++;
    }
    if (flag_grid && locatorType == ELoc::X) continue;
    if (ecr >= ncol) break;
    String string = getLocatorName(locatorType, item);
    _recordWrite("%s", string.c_str());
    ecr++;
  }
  _recordWrite("\n");

  /* Print the variable names */

  _recordWrite("#", "Names");
  VectorInt iatts;
  ecr = 0;
  for (int icol = 0; icol < getFieldNumber(); icol++)
  {
    if (! getLocatorByColumn(icol, &locatorType, &item))
    {
      if (onlyLocator) continue;
      locatorType = ELoc::Z;
    }
    if (flag_grid && locatorType == ELoc::X) continue;
    if (ecr >= ncol) break;
    _recordWrite("%s", getNameByColumn(icol).c_str());
    iatts.push_back(getAttribute(getNameByColumn(icol)));
    ecr++;
  }
  _recordWrite("\n");

  /* Print the array of values */

  _recordWrite("#", "Array of values");
  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    if (!flag_grid && !getSelection(iech)) continue;
    for (int icol = 0; icol < ncol; icol++)
      _recordWrite("%lf", getArray(iech, iatts[icol]));
    _recordWrite("\n");
  }
  return (0);
}

void Db::_variableRead(int *natt_r,
                       int *ndim_r,
                       int *nech_r,
                       std::vector<ELoc>& tabatt,
                       VectorInt& tabnum,
                       VectorString& tabnam,
                       VectorDouble& tab)
{
  char line[LONG_SIZE];
  int  inum, natt, ndim, nval, ecr, mult;
  ELoc iatt;
  double value;

  /* Initializations */

  natt = nval = ndim = 0;

  /* Read the number of variables */

  if (_recordRead("Number of Variables", "%d", &natt)) goto label_end;

  /* Decoding the locators */

  ecr = 0;
  while (1)
  {
    if (ecr >= natt) break;
    if (_recordRead("Locator Name", "%s", line)) goto label_end;
    if (locatorIdentify(line, &iatt, &inum, &mult)) break;
    tabatt.push_back(iatt);
    tabnum.push_back(inum);
    if (iatt == ELoc::X) ndim++;
    ecr++;
  }

  /* Decoding the names */

  ecr = 0;
  while (1)
  {
    if (ecr >= natt) break;
    if (_recordRead("Variable Name", "%s", line)) goto label_end;
    tabnam.push_back(line);
    ecr++;
  }

  /* Read the numeric values */

  while (1)
  {
    if (_recordRead("Numerical value", "%lf", &value)) goto label_end;
    tab.push_back(value);
    nval++;
  }

  label_end:

  /* Returning arguments */

  *natt_r = natt;
  *nech_r = (natt > 0) ? nval / natt : 0;
  *ndim_r = ndim;
  return;
}

void Db::_loadData(const ELoadBy& order, int flag_add_rank, const VectorDouble& tab)
{
  // Preliminary check

  if (getFieldNumber() <= 0) return;
  int jcol = 0;

  // Add the rank (optional)

  if (flag_add_rank)
  {
    for (int iech = 0; iech < getSampleNumber(); iech++)
      setByColumn(iech, jcol, iech + 1);
    setName(jcol, "rank");
    jcol++;
  }

  // Add the input array 'tab' (if provided)

  if (tab.empty()) return;
  int ntab = (flag_add_rank) ? getFieldNumber() - 1 : getFieldNumber();
  int ecr = 0;
  for (int icol = 0; icol < ntab; icol++)
  {
    for (int iech = 0; iech < getSampleNumber(); iech++, ecr++)
    {
      if (order == ELoadBy::SAMPLE)
        setByColumn(iech, jcol, tab[icol + ntab * iech]);
      else
        setByColumn(iech, jcol, tab[ecr]);
    }
    jcol++;
  }
  return;
}

bool Db::_isCountValid(const VectorInt iatts, bool flagOne) const
{
  if (iatts.empty() && flagOne)
  {
    messerr("No variable name corresponding to your criterion");
    return false;
  }
  else
  {
    if (iatts.size() > 1 && flagOne)
    {
      messerr("You wanted to designate a SINGLE variable.");
      messerr("There are several variables matching your criterion:");
      for (unsigned int i = 0; i < iatts.size(); i++)
        messerr("- %s", getName(iatts[i]).c_str());
      return false;
    }
  }
  return true;
}

/**
 * Returns the Number of different facies (labelling starts at 1)
 * The facies variable must be locatorized as ELoc::Z and be unique
 */
int Db::getFaciesNumber(void) const
{
  if (getLocatorNumber(ELoc::Z) != 1)
  {
    messerr("This function requires the number of variables (%d) to be equal to 1",
            getLocatorNumber(ELoc::Z));
    return ITEST;
  }
  int nech = getSampleNumber();

  // Find the number of Facies (labelled starting from 1)

  int nfac = 0;
  for (int iech=0; iech<nech; iech++)
  {
    if (! isActiveAndDefined(iech,0)) continue;
    int ifac = (int) getVariable(iech,0);
    if (ifac <= 0) continue;
    if (ifac > nfac) nfac = ifac;
  }
  return nfac;
}

/****************************************************************************/
/*!
**  Return the vector of ordered samples by increasing coordinate along X
**
** \return    Array containing the increasing order
**
** \remarks  The returned array must be desallocated
**
*****************************************************************************/
VectorInt Db::getSortArray() const
{
  VectorInt rindex;
  VectorDouble xval;

  /* Initializations */

  int nech = getSampleNumber();

  /* Core allocation */

  xval.resize(nech);
  rindex.resize(nech);

  /* Load the arrays */

  for (int iech=0; iech<nech; iech++)
  {
    rindex[iech] = iech;
    xval[iech]   = getCoordinate(iech,0);
  }

  /* Sorting */

  ut_sort_double(0,nech,rindex.data(),xval.data());

  return(rindex);
}

/****************************************************************************/
/*!
 **  Calculates the cosine of the angle between a reference direction
 **  and the increment between two points in the same Db
 **
 ** \return  Cosine of the angle
 **
 ** \param[in]  iech1  rank of the first sample
 ** \param[in]  iech2  rank of the second sample
 ** \param[in]  codir  Direction coefficient
 **
 *****************************************************************************/
double Db::getCosineToDirection(int iech1,
                                int iech2,
                                const VectorDouble& codir) const
{
  double cosdir = 0.;
  double dn1 = 0.;
  double dn2 = 0.;
  for (int idim = 0; idim < getNDim(); idim++)
  {
    double delta = getDistance1D(iech1, iech2, idim);
    if (FFFF(delta)) return TEST;
    cosdir += delta * codir[idim];
    dn1 += delta * delta;
    dn2 += codir[idim] * codir[idim];
  }
  double prod = dn1 * dn2;
  if (prod <= 0.) return (1.);
  return (cosdir / sqrt(prod));
}

