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
#include "geoslib_f_private.h"
#include "geoslib_old_f.h"
#include "geoslib_define.h"

#include "Db/Db.hpp"
#include "Db/PtrGeos.hpp"
#include "Db/DbStringFormat.hpp"
#include "Db/DbGrid.hpp"
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
#include "Basic/Utilities.hpp"
#include "Stats/Classical.hpp"

#include <algorithm>
#include <functional>
#include <math.h>
#include <stdio.h>

Db::Db()
    : AStringable(),
      ASerializable(),
      _ncol(0),
      _nech(0),
      _array(),
      _uidcol(),
      _colNames(),
      _p()
{
  _clear();
}

Db::Db(const Db& r)
    : AStringable(r),
      ASerializable(r),
      _ncol(r._ncol),
      _nech(r._nech),
      _array(r._array),
      _uidcol(r._uidcol),
      _colNames(r._colNames),
      _p(r._p)
{
}

Db& Db::operator=(const Db& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    ASerializable::operator=(r);
    _ncol = r._ncol;
    _nech = r._nech;
    _array = r._array;
    _uidcol = r._uidcol;
    _colNames = r._colNames;
    _p = r._p;
  }
  return *this;
}

Db::~Db()
{
}

int Db::resetFromSamples(int nech,
                         const ELoadBy& order,
                         const VectorDouble& tab,
                         const VectorString& names,
                         const VectorString& locatorNames,
                         int flag_add_rank)
{
  _clear();
  int ncol = (tab.empty()) ? 0 : (int) (tab.size() / nech);
  _ncol = ncol + flag_add_rank;
  _nech = nech;
  resetDims(_ncol, _nech);

  // Load data (if defined)

  if (flag_add_rank) _createRank(0);
  _loadData(tab, names, locatorNames, order, flag_add_rank);

  return 0;
}


/**
 * Creating a Db by reading a CSV file
 *
 * @param filename   Name of the CSV file
 * @param verbose    Verbose flag
 * @param csvfmt     Description of the CSV format
 * @param ncol_max   Maximum number of columns
 * @param nrow_max   Maximum number of rows
 * @param flag_add_rank 1 if the sample rank must be generated
 */
int Db::resetFromCSV(const String& filename,
                     bool verbose,
                     const CSVformat& csvfmt,
                     int ncol_max,
                     int nrow_max,
                     int flag_add_rank)
{
  _clear();
  VectorString names;
  VectorDouble tab;
  int ncol, nrow;

  /* Reading the CSV file */

  if (csv_table_read(filename, csvfmt, (int) verbose, ncol_max, nrow_max,
                     &ncol, &nrow, names, tab))
  {
    messerr("Problem when reading CSV file");
    return 1;
  }

  ncol = (tab.empty()) ? 0 : (int) (tab.size() / nrow);
  _ncol = ncol + flag_add_rank;
  _nech = nrow;
  resetDims(_ncol, _nech);

  // Load data (if defined)

  if (flag_add_rank) _createRank(0);
  _loadData(tab, names, VectorString(), ELoadBy::SAMPLE, flag_add_rank);

  // Set the names
  _defineDefaultNames(flag_add_rank, names);

  // Locators: Try to guess them from the Names
  _defineDefaultLocatorsByNames(flag_add_rank, names);

  return 0;
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
int Db::resetFromBox(int nech,
                     const VectorDouble& coormin,
                     const VectorDouble& coormax,
                     int ndim,
                     int seed,
                     int flag_add_rank)
{
  _clear();
  if (! coormin.empty()) ndim = (int) coormin.size();
  if (! coormax.empty()) ndim = MIN(ndim, (int) coormax.size());
  _ncol = ndim + flag_add_rank;
  _nech = nech;
  resetDims(_ncol, _nech);

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
  setLocatorsByUID(ndim, jcol, ELoc::X);

  return 0;
}

/**
 * Create a Db from a single sample whose coordinates are provided in 'tab'
 * @param tab Array containing the coordinates of the single sample
 * @param flag_add_rank 1 if the Sample ranks must be generated
 */
int Db::resetFromOnePoint(const VectorDouble& tab, int flag_add_rank)
{
  _clear();

  int ndim = static_cast<int> (tab.size());
  _ncol = ndim + flag_add_rank;
  _nech = 1;
  resetDims(_ncol, _nech);

  // Generate the sample number
  if (flag_add_rank) _createRank(0);

  // Load the coordinates
  VectorString names = generateMultipleNames("x", ndim);
  VectorDouble tabloc = tab;
  if (tabloc.empty()) tabloc.resize(ndim,0.);
  _loadData(tabloc, names, VectorString(), ELoadBy::SAMPLE, flag_add_rank);

  int jcol = 0;
  if (flag_add_rank) jcol++;
  setLocatorsByUID(ndim, jcol, ELoc::X);

  return 0;
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

bool Db::isUIDValid(int iuid) const
{
  if (iuid < 0 || iuid >= getUIDMaxNumber())
  {
    mesArg("UID Index", iuid, getUIDMaxNumber());
    return false;
  }
  return true;
}

bool Db::isColIdxValid(int icol) const
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

bool Db::isSampleIndicesValid(const VectorInt& iechs, bool useSel) const
{
  for (int i = 0; i < (int) iechs.size(); i++)
  {
    int iech = iechs[i];
    if (iech < 0 || iech >= getSampleNumber(useSel))
    {
      mesArg("Sample Index", iech, getSampleNumber(useSel));
      return false;
    }
  }
  return true;
}

bool Db::isLocatorIndexValid(const ELoc& locatorType, int locatorIndex) const
{
  if (!isLocatorTypeValid(locatorType)) return false;
  bool ok = _p.at(locatorType).isLocatorIndexValid(locatorIndex);
  if (! ok)
    messerr("Problem in the identification of Locator %d", locatorType.getValue());
  return ok;
}

int Db::getColIdxByUID(int iuid) const
{
  if (!isUIDValid(iuid)) return -1;
  int icol = _uidcol[iuid];
  return icol;
}

VectorInt Db::getColIdxsByUID(const VectorInt iuids) const
{
  VectorInt cols(iuids.size());
  for (unsigned int i = 0; i < iuids.size(); i++)
    cols[i] = getColIdxByUID(iuids[i]);
  return cols;
}

int Db::getUIDByColIdx(int icol) const
{
  if (!isColIdxValid(icol)) return -1;
  for (int iuid = 0; iuid < getUIDMaxNumber(); iuid++)
    if (_uidcol[iuid] == icol) return iuid;
  return -1;
}

int Db::getUIDByLocator(const ELoc& locatorType, int locatorIndex) const
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
int Db::getColIdxByLocator(const ELoc& locatorType, int locatorIndex) const
{
  if (!isLocatorTypeValid(locatorType)) return -1;
  if (!isLocatorIndexValid(locatorType,locatorIndex)) return -1;
  const PtrGeos& p = _p.at(locatorType);
  int icol = getColIdxByUID(p.getLocatorByIndex(locatorIndex));
  return (icol);
}

int Db::getLocatorNumber(const ELoc& locatorType) const
{
  if (!isLocatorTypeValid(locatorType)) return -1;
  const PtrGeos& p = _p.at(locatorType);
  return p.getLocatorNumber();
}

int Db::_findUIDInLocator(const ELoc& locatorType, int iuid) const
{
  if (!isLocatorTypeValid(locatorType)) return -1;
  const PtrGeos& p = _p.at(locatorType);
  if (!isUIDValid(iuid)) return -1;
  for (int locatorIndex = 0; locatorIndex < p.getLocatorNumber(); locatorIndex++)
    if (p.getLocatorByIndex(locatorIndex) == iuid) return (locatorIndex);
  return -1;
}

int Db::_findColumnInLocator(const ELoc& locatorType, int icol) const
{
  if (!isLocatorTypeValid(locatorType)) return -1;
  int iuid = getUIDByColIdx(icol);
  return _findUIDInLocator(locatorType, iuid);
}

/**
 * Find the locator characteristics of a given Column
 * @param icol       Index of the target column
 * @param ret_locatorType Locator type
 * @param ret_locatorIndex Locator index (starting from 0)
 * @return true if the target variable has a locator assigned and false otherwise
 */
bool Db::getLocatorByColIdx(int icol,
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
        int jcol = getColIdxByUID(p.getLocatorByIndex(i));
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

bool Db::getLocatorByUID(int iuid,
                         ELoc* ret_locatorType,
                         int* ret_locatorIndex) const
{
  if (!isUIDValid(iuid)) return false;
  int icol = getColIdxByUID(iuid);
  return getLocatorByColIdx(icol, ret_locatorType, ret_locatorIndex);
}

/**
 * Return the locator information corresponding to the input variable
 * @param name Input variable name (unique)
 * @param ret_locatorType Locator Type
 * @param ret_locatorIndex Locator Index (starting from 0)
 * @return
 */
bool Db::getLocator(const String& name,
                    ELoc *ret_locatorType,
                    int *ret_locatorIndex) const
{
  VectorInt iuids = _ids(name, true);
  if (iuids.empty()) return false;
  return getLocatorByUID(iuids[0], ret_locatorType, ret_locatorIndex);
}

/**
 * Check if a variable (specified by its name) matches the required locator
 * @param name         Name of the target Variable
 * @param locatorType  Characteristics of the required Locator Type
 * @param locatorIndex Index of the required Locator (or -1)
 * @return
 */
bool Db::hasLocatorDefined(const String& name,
                           const ELoc& locatorType,
                           int locatorIndex) const
{
  VectorInt iuids = _ids(name, true);
  if (iuids.empty()) return false;
  if (!isUIDValid(iuids[0])) return false;
  int icol = getColIdxByUID(iuids[0]);
  ELoc ret_locatorType;
  int ret_locatorIndex;
  getLocatorByColIdx(icol, &ret_locatorType, &ret_locatorIndex);
  if (ret_locatorType != locatorType) return false;
  if (locatorIndex >= 0 && ret_locatorIndex != locatorIndex) return false;
  return true;
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
      (void) getLocatorByColIdx(icol, &type, &item);
      if (type != locatorType) continue;
    }
    String string = _getLocatorNameByColIdx(icol);
    retval.push_back(string);
  }
  return retval;
}

bool Db::isUIDDefined(int iuid) const
{
  if (!isUIDValid(iuid)) return false;
  int icol = getColIdxByUID(iuid);
  if (!isColIdxValid(icol)) return false;
  return (_uidcol[icol] >= 0);
}

VectorString Db::expandNameList(const VectorString& names) const
{
  return expandList(_colNames, names);
}

VectorString Db::expandNameList(const String& names) const
{
  return expandList(_colNames, names);
}

VectorInt Db::_ids(const String& name, bool flagOne, bool verbose) const
{
  VectorString exp_names = expandNameList(name);
  VectorInt iuids = _getUIDsBasic(exp_names);
  if (! _isCountValid(iuids, flagOne, verbose)) return VectorInt();
  return iuids;
}

VectorInt Db::_ids(const VectorString& names, bool flagOne, bool verbose) const
{
  VectorString exp_names = expandNameList(names);
  VectorInt iuids = _getUIDsBasic(exp_names);
  if (! _isCountValid(iuids, flagOne, verbose)) return VectorInt();
  return iuids;
}

VectorInt Db::_ids(const ELoc& locatorType, bool flagOne, bool verbose) const
{
  VectorString exp_names = getNamesByLocator(locatorType);
  VectorInt iuids = _getUIDsBasic(exp_names);
  if (! _isCountValid(iuids, flagOne, verbose)) return VectorInt();
  return iuids;
}

VectorInt Db::_ids(const VectorInt& iuids, bool flagOne, bool verbose) const
{
  VectorString exp_names = getNamesByUID(iuids);
  if (! _isCountValid(iuids, flagOne, verbose)) return VectorInt();
  return iuids;
}

void Db::resetDims(int ncol, int nech)
{
  _ncol = ncol;
  _nech = nech;

  /* The UID pointers */

  _uidcol.resize(ncol);
  for (int i = 0; i < ncol; i++)
    _uidcol[i] = i;

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
 * Set the value by Sample and UID
 * @param iech  Index of the Sample
 * @param iuid  Index of the UID
 * @param value Value to be assigned
 */
void Db::setArray(int iech, int iuid, double value)
{
  if (!isSampleIndexValid(iech)) return;
  int icol = getColIdxByUID(iuid);
  if (!isColIdxValid(icol)) return;
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
  int iuid = getUID(name);
  if (iuid < 0) return TEST;
  return getArray(iech, iuid);
}

/**
 * Sets the value of the 'iech' sample of the variable 'name'
 *
 * This function does not use 'ids' mechanism in order to allow
 * referring to a non-existing variable
 */
void Db::setValue(const String& name, int iech, double value)
{
  int iuid  = getUID(name);
  if (iuid < 0) return;
  setArray(iech, iuid, value);
}

/**
 * Return the value defined by Sample and UID
 * @param iech Sample Index
 * @param iuid UID Index
 * @return
 */
double Db::getArray(int iech, int iuid) const
{
  if (!isSampleIndexValid(iech)) return (TEST);
  int icol = getColIdxByUID(iuid);
  if (!isColIdxValid(icol)) return (TEST);
  return (_array[_getAddress(iech, icol)]);
}

VectorDouble Db::getArray(int iuid, bool useSel) const
{
  int nech = getSampleNumber();
  VectorDouble sel, tab;
  if (!isUIDValid(iuid)) return tab;

  tab.resize(nech);
  if (useSel) sel = getSelection();

  int ecr = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    if (useSel && !sel.empty() && sel[iech] == 0) continue;
    tab[ecr] = getArray(iech, iuid);
    ecr++;
  }
  tab.resize(ecr);
  return tab;
}

VectorDouble Db::getArrayBySample(int iech) const
{
  VectorInt uids = getAllUIDs();
  VectorDouble vals;
  for (int iuid = 0; iuid < (int) uids.size(); iuid++)
    vals.push_back(getArray(iech, uids[iuid]));
  return vals;
}

void Db::setArrayBySample(int iech, const VectorDouble& vec)
{
  VectorInt uids = getAllUIDs();
  if ((int) uids.size() != (int) vec.size())
  {
    messerr("Dimension of 'vec'(%d) does not match number of columns(%)",
            (int) vec.size(),(int) uids.size());
    return;
  }
  for (int iuid = 0; iuid < (int) uids.size(); iuid++)
    setArray(iech, uids[iuid], vec[iuid]);
}

void Db::updArray(int iech, int iuid, int oper, double value)
{
  if (!isSampleIndexValid(iech)) return;
  double oldval = getArray(iech, iuid);
  setArray(iech, iuid, _updateValue(oper, oldval, value));
}

VectorDouble Db::getSampleCoordinates(int iech) const
{
  VectorDouble coor(getNDim());
  getSampleCoordinates(iech, coor);
  return coor;
}

VectorDouble Db::getSampleLocators(const ELoc& locatorType, int iech) const
{
  VectorDouble vec;
  if (!isLocatorTypeValid(locatorType)) return vec;
  int number = getLocatorNumber(locatorType);
  if (number <= 0) return vec;
  vec.resize(number);
  for (int i = 0; i < number; i++)
    vec[i] = getFromLocator(locatorType, iech, i);
  return vec;
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
double Db::getCoordinate(int iech, int idim, bool /*flag_rotate*/) const
{
  if (idim >= getNDim()) return TEST;
  return getFromLocator(ELoc::X, iech, idim);
}

void Db::getCoordinatesInPlace(int iech, VectorDouble& coor, bool flag_rotate) const
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

double Db::getDistance(int iech, int jech) const
{
  int ndim = getNDim();
  double dist = 0.;
  for (int idim = 0; idim < ndim; idim++)
  {
    double v1 = getCoordinate(iech, idim);
    if (FFFF(v1)) return TEST;
    double v2 = getCoordinate(jech, idim);
    if (FFFF(v2)) return TEST;
    double delta = v1 - v2;
    dist += delta * delta;
  }
  return sqrt(dist);
}

/**
 * Calculate the distance vector in place
 * @param iech Rank of the first sample
 * @param jech Rank of the second sample (from db2 if db2 provided)
 * @param dd   Vector for distances (It must be dimensioned to getNDim())
 * @param db2  Second Db if different from current one (or nullptr)
 * @return
 */
int Db::getDistanceVec(int iech, int jech, VectorDouble& dd, const Db* db2) const
{
  int ndim = getNDim();
  double v1;
  double v2;
  for (int idim = 0; idim < ndim; idim++)
  {
    v1 = getCoordinate(iech, idim);
    if (FFFF(v1)) return 1;
    if (db2 == nullptr)
      v2 = getCoordinate(jech, idim);
    else
      v2 = db2->getCoordinate(jech, idim);
    if (FFFF(v2)) return 1;
    dd[idim] = v1 - v2;
  }
  return 0;
}

/**
 * Constitute a Vector of Vector of coordinates at a given sample, for all (active) samples
 * @param useSel
 * @return
 */
VectorVectorDouble Db::getAllCoordinates(bool useSel) const
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
  int icol = getColIdxByLocator(ELoc::X, idim);
  if (!isColIdxValid(icol)) return;
  _array[_getAddress(iech, icol)] = value;
}

void Db::setFromLocator(const ELoc& locatorType,
                        int iech,
                        int locatorIndex,
                        double value)
{
  if (!isSampleIndexValid(iech)) return;
  int icol = getColIdxByLocator(locatorType,locatorIndex);
  if (!isColIdxValid(icol)) return;
  _array[_getAddress(iech, icol)] = value;
}

double Db::getFromLocator(const ELoc& locatorType,
                          int iech,
                          int locatorIndex) const
{
  if (!isSampleIndexValid(iech)) return (TEST);
  int icol = getColIdxByLocator(locatorType, locatorIndex);
  if (!isColIdxValid(icol)) return (TEST);
  return (_array[_getAddress(iech, icol)]);
}

int Db::getFromLocatorNumber(const ELoc& locatorType) const
{
  const PtrGeos& p = _p.at(locatorType);
  return p.getLocatorNumber();
}

void Db::_clear(void)
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

int Db::_getUIDcol(int iuid) const
{
  if (! isUIDValid(iuid)) return ITEST;
  return _uidcol[iuid];
}

int Db::_getAddress(int iech, int icol) const
{
  return ((iech) + _nech * icol);
}

String Db::_summaryLocators(void) const
{
  std::stringstream sstr;

  /* Loop on the pointers */

  sstr << toTitle(1, "List of locators");
  int rank = 0;
  auto it = ELoc::getIterator();
  while (it.hasNext())
  {
    if (*it != ELoc::UNKNOWN)
    {
      const PtrGeos& p = _p.at(*it);
      if (p.getLocatorNumber() > 0)
      {
        sstr << p.dumpLocator(rank, *it);
        sstr << "- Columns    = ";
        for (int locatorIndex = 0; locatorIndex < p.getLocatorNumber(); locatorIndex++)
          sstr << getColIdxByUID(p.getLocatorByIndex(locatorIndex)) << " ";
        sstr << std::endl;
        rank++;
      }
    }
    it.toNext();
  }
  return sstr.str();
}

String Db::_summaryUIDs(void) const
{
  std::stringstream sstr;

  sstr << toTitle(1, "List of unsorted UIDs");
  sstr << "Maximum number of positions = " << getUIDMaxNumber() << std::endl;
  sstr << "Number of Columns           = " << getColumnNumber() << std::endl;

  /* Loop on the UIDs */

  if (getUIDMaxNumber() <= 0) return sstr.str();

  sstr << "UID = ";
  for (int iuid = 0; iuid < getUIDMaxNumber(); iuid++)
    sstr << _uidcol[iuid] << " ";
  sstr << std::endl;
  return sstr.str();
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
 * @param cleanSameLocator When TRUE, clean variables with same locator beforehand
 *
 */
void Db::setLocators(const VectorString &names,
                     const ELoc& locatorType,
                     int locatorIndex,
                     bool cleanSameLocator)
{
  if (!isLocatorTypeValid(locatorType, true)) return;
  VectorInt iuids = _ids(names, false);
  if (iuids.empty()) return;

  if (cleanSameLocator) clearLocators(locatorType);

  for (unsigned int i = 0; i < iuids.size(); i++)
    setLocatorByUID(iuids[i], locatorType, locatorIndex + i);
}

/**
 * Define the Locator(s) for the given variable(s)
 * @param names Set of variable names
 * @param locatorType Locator Type
 * @param locatorIndex Locator Index (for the first variable) (starting from 0)
 * @param cleanSameLocator When TRUE, clean variables with same locator beforehand
 */
void Db::setLocator(const String &names,
                    const ELoc& locatorType,
                    int locatorIndex,
                    bool cleanSameLocator)
{
  if (!isLocatorTypeValid(locatorType, true)) return;
  VectorInt iuids = _ids(names, false);
  if (iuids.empty()) return;

  if (cleanSameLocator) clearLocators(locatorType);

  for (unsigned int i = 0; i < iuids.size(); i++)
    setLocatorByUID(iuids[i], locatorType, locatorIndex + i);
}

/**
 * Setting the locator for a variable designated by its UID
 * @param iuid          Index of the UID
 * @param locatorType   Type of locator (include ELoc::UNKNOWN)
 * @param locatorIndex  Rank in the Locator (starting from 0)
 * @param cleanSameLocator When TRUE, clean variables with same locator beforehand
 * @remark: At this stage, no check is performed to see if items
 * @remark: are consecutive and all defined
 * @remark: This allow using this function in any order.
 */
void Db::setLocatorByUID(int iuid,
                         const ELoc& locatorType,
                         int locatorIndex,
                         bool cleanSameLocator)
{
  if (!isUIDValid(iuid)) return;
  if (!isLocatorTypeValid(locatorType, true)) return;
  if (locatorIndex < 0) return;

  // Optional clean

  if (cleanSameLocator) clearLocators(locatorType);

  /* Cancel any locator referring to this column */

  auto it = ELoc::getIterator();
  while (it.hasNext())
  {
    if (*it != ELoc::UNKNOWN)
    {
      PtrGeos& p = _p[*it];
      int found = _findUIDInLocator(*it, iuid);
      if (found >= 0)
        p.erase(found);
    }
    it.toNext();
  }

  // Check if this locator already exists for the current pointer type
  // Warning: the following code does not forbid declaring locatorIndex
  // in incorrect order. This must be kept as long as the Demonstration files
  // use the db.locate() of unsorted ranks

  if (locatorType != ELoc::UNKNOWN)
  {
    PtrGeos& p = _p[locatorType];
    int nitem = p.getLocatorNumber();
    if (locatorIndex >= nitem)
    {
      p.resize(locatorIndex + 1);
    }
    p.setLocatorByIndex(locatorIndex, iuid);
  }
}

void Db::setLocatorByColIdx(int icol,
                            const ELoc& locatorType,
                            int locatorIndex,
                            bool cleanSameLocator)
{
  if (!isColIdxValid(icol)) return;

  int iuid = getUIDByColIdx(icol);
  setLocatorByUID(iuid, locatorType, locatorIndex, cleanSameLocator);
}

String Db::_getLocatorNameByColIdx(int icol) const
{
  ELoc locatorType;
  int locatorIndex;
  (void) getLocatorByColIdx(icol, &locatorType, &locatorIndex);
  return getLocatorName(locatorType, locatorIndex);
}

/**
 * Set the Locators for a set of variables identified by their UID
 * @param number        Number of variables to be set
 * @param iuid          Index of the first UID
 * @param locatorType   Type of the Locator (include ELoc::UNKNOWN)
 * @param locatorIndex  Rank of the first Locator index (starting from 0)
 * @param cleanSameLocator When TRUE, clean variables with same locator beforehand
 */
void Db::setLocatorsByUID(int number,
                          int iuid,
                          const ELoc& locatorType,
                          int locatorIndex,
                          bool cleanSameLocator)
{
  if (!isLocatorTypeValid(locatorType, true)) return;

  if (cleanSameLocator) clearLocators(locatorType);

  for (int i = 0; i < number; i++)
    setLocatorByUID(iuid+i, locatorType, locatorIndex + i);
}

void Db::setLocatorsByUID(const VectorInt& iuids,
                          const ELoc& locatorType,
                          int locatorIndex,
                          bool cleanSameLocator)
{
  if (!isLocatorTypeValid(locatorType, true)) return;

  if (cleanSameLocator) clearLocators(locatorType);

  int number = (int) iuids.size();
  for (int i = 0; i < number; i++)
    setLocatorByUID(iuids[i], locatorType, locatorIndex + i);
}

void Db::setLocatorsByColIdx(const VectorInt& icols,
                              const ELoc& locatorType,
                              int locatorIndex,
                              bool cleanSameLocator)
{
  if (cleanSameLocator) clearLocators(locatorType);

  for (int icol = 0; icol < (int) icols.size(); icol++)
  {
    int iuid = getUIDByColIdx(icol);
    setLocatorByUID(iuid, locatorType, locatorIndex + icol);
  }
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
 * @return Rank of the first UID
 */
int Db::addColumnsByConstant(int nadd,
                             double valinit,
                             const String &radix,
                             const ELoc &locatorType,
                             int locatorIndex,
                             int nechInit)
{
  int ncol = _ncol;
  int nmax = getUIDMaxNumber();
  int nnew = ncol + nadd;
  if (nadd <= 0) return (-1);

  /* Case of an empty Db, define the number of samples using 'nechInit' */

  if (_nech <= 0) _nech = nechInit;

  /* Dimension the array */

  _array.resize(_nech * nnew);

  /* Dimension the UID pointer */

  _uidcol.resize(nmax + nadd);
  for (int i = 0; i < nadd; i++)
    _uidcol[nmax + i] = ncol + i;

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
    setLocatorsByUID(nadd, nmax, locatorType, locatorIndex);

  _ncol += nadd;

  return (nmax);
}


void Db::addColumnsByVVD(const VectorVectorDouble tab,
                    const String& radix,
                    const ELoc& locatorType,
                    int locatorIndex,
                    bool useSel,
                    double valinit,
                    int nvar)
{
  VectorDouble tabv;
  for(auto &e : tab)
    for(auto &f : e)
      tabv.push_back(f);
  addColumns(tabv,radix,locatorType,locatorIndex,useSel,valinit,nvar);
}


/**
 * Add one or several columns to an already existing Db. This is performed
 * by providing an array of values 'tab'. Its dimension must be equal to the
 * number of samples (or active samples if 'useSel' is true, times the number
 * of variables 'nvar'.
 * @param tab    Array to be loaded
 * @param radix  Generic name for the newly created variables
 * @param locatorType Generic locator assigned to new variables
 * @param locatorIndex   Locator index (starting from 0)
 * @param useSel true if the Selection must be taken into account
 * @param valinit initial value (for unselected samples)
 * @param nvar   Number of variables loaded
 *
 * @return Rank of the first UID
 *
 * @remark When 'useSel' is used, you must have a Selection already defined. Then the number
 * @remark of samples provided in 'tab' must match the number of active samples
 * @remark When a vector 'tab' is provided, the number of variables 'nvar'
 * @remark is calculated as its size divided by the number of samples in the grid.
 */
int Db::addColumns(const VectorDouble &tab,
                   const String &radix,
                   const ELoc &locatorType,
                   int locatorIndex,
                   bool useSel,
                   double valinit,
                   int nvar)
{
  // Particular case where the Db is empty.
  // Set its dimension to the number of samples of the input array 'tab'
  if (_nech <= 0) _nech = static_cast<int> (tab.size()) / nvar;

  // Check dimensions
  int nech = getSampleNumber(useSel);
  nvar = (int) tab.size() / nech;
  if ((int) tab.size() != nvar * nech)
  {
    messerr("Db::addColumns : Incompatibility between dimension of 'tab' (%d)", tab.size());
    messerr("and 'nvar'(%d) * 'nech'(%d)", nvar, nech);
    return 1;
  }

  // Adding the new Columns
  int iuid = addColumnsByConstant(nvar, valinit, radix, locatorType, locatorIndex);
  if (iuid < 0) return 1;

  setColumnByUID(tab, iuid, useSel);

  const double* local = tab.data();
  for (int ivar = 0; ivar < nvar; ivar++)
    setColumnByUIDOldStyle(&local[ivar * nech], iuid + ivar, useSel);

  return iuid;
}

void Db::setColumnByColIdxOldStyle(const double* tab, int icol, bool useSel)
{
  if (!isColIdxValid(icol)) return;
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
    setValueByColIdx(iech, icol, value);
  }
}

void Db::setColumnByColIdx(const VectorDouble& tab, int icol, bool useSel)
{
  setColumnByColIdxOldStyle(tab.data(), icol, useSel);
}

void Db::setColumnsByColIdx(const VectorDouble& tabs, const VectorInt& icols, bool useSel)
{
  int nech = getSampleNumber(useSel);
  if ((int) icols.size() * nech != (int) tabs.size())
  {
    messerr("Dimensions of 'icols'(%d), 'nech'(%d) and 'tabs'(%d) are inconsistent",
            (int) icols.size(), nech, (int) tabs.size());
    return;
  }
  int lec = 0;
  VectorDouble tabloc(nech);
  for (int i = 0; i < (int) icols.size(); i++)
  {
    int icol = icols[i];
    for (int j = 0; j < getSampleNumber(useSel); j++) tabloc[j] = tabs[lec++];
    setColumnByColIdx(tabloc, icol, useSel);
  }
}

void Db::setColumnByUIDOldStyle(const double* tab, int iuid, bool useSel)
{
  if (!isUIDValid(iuid)) return;
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
    setArray(iech, iuid, value);
  }
}

void Db::setColumnByUID(const VectorDouble& tab, int iuid, bool useSel)
{
  setColumnByUIDOldStyle(tab.data(), iuid, useSel);
}

/**
 * Set the values for an already existing Column.
 * Note that, if the Column does not exist, this Column is added beforehand
 * @param tab    Array of values to be stored in the target Column
 * @param name   Name of the Column
 * @param locatorType Locator type
 * @param locatorIndex   Locator index (starting from 0)
 * @param useSel Should an already existing Selection be taken into account
 *
 * @remark: Arguments 'locatorType'  and 'locatorIndex' are only used
 * @remark: for newly added variables
 */
void Db::setColumn(const VectorDouble& tab, const String& name,
                   const ELoc& locatorType, int locatorIndex, bool useSel)
{
  VectorInt iuids = _ids(name, true, false);
  if (iuids.empty())
  {
    (void) addColumns(tab, name, locatorType, locatorIndex, useSel);
  }
  else
  {
    setColumnByUIDOldStyle(tab.data(), iuids[0], useSel);
  }
}

void Db::duplicateColumnByUID(int iuid_in, int iuid_out)
{
  if (!isUIDValid(iuid_in)) return;
  if (!isUIDValid(iuid_out)) return;
  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    double value = getArray(iech, iuid_in);
    setArray(iech, iuid_out, value);
  }
}

void Db::deleteColumn(const String& name)
{
  VectorInt iuids = _ids(name, false);
  if (iuids.empty()) return;

  for (unsigned int i = 0; i < iuids.size(); i++)
    deleteColumnByUID(iuids[i]);
}

void Db::deleteColumns(const VectorString& names)
{
  VectorInt iuids = _ids(names, false);
  if (iuids.empty()) return;

  for (unsigned int i = 0; i < iuids.size(); i++)
    deleteColumnByUID(iuids[i]);
}

void Db::deleteColumnsByColIdx(const VectorInt& icols)
{
  if (icols.empty()) return;

  VectorInt v = ut_ivector_sort(icols, false);

  for (unsigned int i = 0; i < v.size(); i++)
    deleteColumnByColIdx(v[i]);
}

DbGrid* Db::coveringDb(const VectorInt& nodes,
                             const VectorDouble& dcell,
                             const VectorDouble& origin,
                             const VectorDouble& margin) const
{
   return DbGrid::createCoveringDb(this,nodes,dcell,origin,margin);
}

void Db::deleteColumnsByUID(const VectorInt& iuids)
{
  if (iuids.empty()) return;

  for (unsigned int i = 0; i < iuids.size(); i++)
    deleteColumnByUID(iuids[i]);
}

/**
 * Add the contents of the 'tab' as a Selection
 * @param tab Input array
 * @param name Name given to the newly created Selection variable
 * @param combine How to combine with an already existing selection (see combineSelection() for details)
 * @return Rank of the newly created Column within the Data Base
 * @remark The Selection is set to True if tab is not zero and to False otherwise.
 * @remark If the dimension of 'tab' does not match the number of samples in the Db
 * @remark the action is cancelled (a message is issued)
 */
int Db::addSelection(const VectorDouble& tab, const String& name, const String& combine)
{
  int nech = getSampleNumber();
  VectorDouble sel(nech);

  if (tab.empty())
  {
    for (int i = 0; i < nech; i++)
      sel[i] = 1.;
  }
  else
  {
    if (nech != (int) tab.size())
    {
      messerr("Dimension of 'tab' (%d) does not match the number of samples (%d)",
              (int) tab.size(),nech);
      messerr("Action is cancelled");
      return -1;
    }

    for (int iech = 0; iech < nech; iech++)
    {
      sel[iech] = (tab[iech] != 0.) ? 1. : 0.;
    }
  }

  // Convert the input array into a selection (0 or 1)

  combineSelection(sel, combine);
  int iuid = addColumns(sel, name, ELoc::SEL);
  return iuid;
}

/**
 * Add a Selection by considering the input 'ranks' vector which give the ranks
 * of the active samples (starting from 0)
 * @param ranks   Vector of ranks of active samples
 * @param name Name given to the newly created Selection variable
 * @param combine How to combine with an already existing selection (see combineSelection() for details)
 * @return
 */
int Db::addSelectionByRanks(const VectorInt &ranks,
                            const String &name,
                            const String &combine)
{
  int nech = getSampleNumber();
  VectorDouble sel(nech, 0.);

  for (int i = 0; i < (int) ranks.size(); i++)
    sel[ranks[i]] = 1.;

  // Convert the input array into a selection (0 or 1)

  combineSelection(sel, combine);
  int iuid = addColumns(sel, name, ELoc::SEL);
  return iuid;
}

/**
 * Create a selection around the only defined values of the target variable
 * @param testvar Name of the target variable
 * @param limits  Limits defining the Definition Domain to be tested (optional)
 * @param name    Name of the newly created selection
 * @param combine How to combine with an already existing selection (see combineSelection() for details)
 * @return The rank of the newly created selection variable within the Db
 */
int Db::addSelectionByLimit(const String& testvar,
                            const Limits& limits,
                            const String& name,
                            const String& combine)
{
  int nech = getSampleNumber();
  VectorDouble sel(nech);

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
    sel[iech] = answer;
  }
  combineSelection(sel, combine);
  int iuid = addColumns(sel, name, ELoc::SEL);

  return iuid;
}

/**
 * Add samples to the Data Base
 * @param nadd    Number of samples to be added
 * @param valinit Default value given to the added samples
 * @return Index of the first newly added sample (or -1 if adding samples is not authorized)
 */
int Db::addSamples(int nadd, double valinit)
{
  if (! mayChangeSampleNumber())
  {
    messerr("This type of Data Base does not allow modifying the Count of Samples");
    return -1;
  }
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

/**
 * Deleting a sample
 * @param e_del Index of the sample to be deleted
 * @return 0 if successfull or -1 if sample deletion is not authorized
 */
int Db::deleteSample(int e_del)
{
  if (! mayChangeSampleNumber())
  {
    messerr("This type of Data Base does not allow modifying the Count of Samples");
    return 1;
  }
  int nech = _nech;
  int nnew = nech - 1;
  if (!isSampleIndexValid(e_del)) return 1;

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
  return 0;
}

void Db::deleteColumnByColIdx(int icol_del)
{
  if (! isColIdxValid(icol_del)) return;
  VectorInt iuids = _ids(_colNames[icol_del],true);
  if (iuids.empty()) return;
  deleteColumnByUID(iuids[0]);
}

/**
 * Delete an UID
 * @param iuid_del Rank of the UID to be deleted
 */
void Db::deleteColumnByUID(int iuid_del)
{
  int ncol = _ncol;
  int nech = _nech;
  int nmax = getUIDMaxNumber();
  int nnew = ncol - 1;
  if (!isUIDValid(iuid_del)) return;

  /* Identify the column to be deleted */

  int c_del = getColIdxByUID(iuid_del);
  if (!isColIdxValid(c_del)) return;
  _uidcol[iuid_del] = -1;
  for (int iuid = 0; iuid < nmax; iuid++)
  {
    if (_uidcol[iuid] < c_del) continue;
    _uidcol[iuid]--;
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
      int found = _findUIDInLocator(*it, iuid_del);
      if (found >= 0) p.erase(found);
    }
    it.toNext();
  }

  /* Resize the variables names */

  _colNames.erase(_colNames.begin() + c_del);

  /* Set the error return code */

  _ncol = nnew;
}

void Db::deleteColumnsByLocator(const ELoc& locatorType)
{
  if (!isLocatorTypeValid(locatorType)) return;
  PtrGeos& p = _p[locatorType];
  int nitem = p.getLocatorNumber();
  // Loop is performed downwards as PtrGeos is modified by called routine
  for (int locatorIndex = nitem - 1; locatorIndex >= 0; locatorIndex--)
  {
    deleteColumnByUID(p.getLocatorByIndex(locatorIndex));
  }
}

/**
 * Returns a Vector containing the minimum and maximum along a Space dimension
 * @param idim   Rank of the Space dimension
 * @param useSel True if the possible selection must be taken intao account
 * @return
 */
VectorDouble Db::getExtrema(int idim, bool useSel) const
{
  VectorDouble ext;
  if (!isDimensionIndexValid(idim)) return ext;
  VectorDouble coor = getCoordinates(idim, useSel);
  ext.push_back(ut_vector_min(coor));
  ext.push_back(ut_vector_max(coor));
  return ext;
}

VectorVectorDouble Db::getExtremas(bool useSel) const
{
  VectorVectorDouble exts;
  for (int idim = 0; idim < getNDim(); idim++)
    exts.push_back(getExtrema(idim, useSel));
  return exts;
}

VectorDouble Db::getCoorMinimum(bool useSel) const
{
  VectorDouble ext;
  for (int idim = 0; idim < getNDim(); idim++)
  {
    VectorDouble coor = getCoordinates(idim, useSel);
    ext.push_back(ut_vector_min(coor));
  }
  return ext;
}

VectorDouble Db::getCoorMaximum(bool useSel) const
{
  VectorDouble ext;
  for (int idim = 0; idim < getNDim(); idim++)
  {
    VectorDouble coor = getCoordinates(idim, useSel);
    ext.push_back(ut_vector_max(coor));
  }
  return ext;
}

VectorDouble Db::getCenter(bool useSel) const
{
  int ndim = getNDim();
  VectorDouble center(ndim);
  for (int idim = 0; idim < ndim; idim++)
    center[idim] = getCenter(idim, useSel);
  return center;
}

double Db::getCenter(int idim, bool useSel) const
{
  if (!isDimensionIndexValid(idim)) return TEST;
  VectorDouble coor = getCoordinates(idim, useSel);
  double mini = ut_vector_min(coor);
  double maxi = ut_vector_max(coor);
  return ((mini + maxi) / 2.);
}

double Db::getExtension(int idim, bool useSel) const
{
  if (!isDimensionIndexValid(idim)) return 0.;
  VectorDouble coor = getCoordinates(idim, useSel);
  double mini = ut_vector_min(coor);
  double maxi = ut_vector_max(coor);
  return maxi - mini;
}

double Db::getExtensionDiagonal(bool useSel) const
{
  int ndim = getNDim();
  double total = 0.;
  for (int idim = 0; idim < ndim; idim++)
  {
    double delta = getExtension(idim, useSel);
    total += delta * delta;
  }
  return sqrt(total);
}

void Db::getExtensionInPlace(VectorDouble &mini, VectorDouble &maxi)
{
  int ndim = getNDim();
  if (ndim != (int) mini.size()) mini.resize(ndim,TEST);
  if (ndim != (int) maxi.size()) maxi.resize(ndim,TEST);

  /* Loop on the space dimension */

  for (int idim = 0; idim < getNDim(); idim++)
  {
    VectorDouble coor = getCoordinates(idim, true);
    double vmin = ut_vector_min(coor);
    double vmax = ut_vector_max(coor);
    if (FFFF(mini[idim]) || vmin < mini[idim]) mini[idim] = vmin;
    if (FFFF(maxi[idim]) || vmax > maxi[idim]) maxi[idim] = vmax;
  }
}

/**
 * Return a Unit calculated for a Db (in a given Space dimension)
 * @param idim Rank of the Space dimension
 * @return
 *
 * @remarks This unit is defined as 1/1000 of the extension in the given space dimension
 */
double Db::getUnit(int idim) const
{
  double delta = getExtension(idim);
  return delta / 1000.;
}

double Db::getColumnSize(bool useSel) const
{
  double diag = 0.;
  for (int idim = 0; idim < getNDim(); idim++)
  {
    VectorDouble ext = getExtrema(idim, useSel);
    double delta = ext[1] - ext[0];
    diag += delta * delta;
  }
  return sqrt(diag);
}

double Db::getMinimum(const String& name, bool useSel) const
{
  VectorInt iuids = _ids(name, true);
  if (iuids.empty()) return TEST;
  VectorDouble tab = getColumnByUID(iuids[0], useSel);
  return ut_vector_min(tab);
}

double Db::getMaximum(const String& name, bool useSel) const
{
  VectorInt iuids = _ids(name, true);
  if (iuids.empty()) return TEST;
  VectorDouble tab = getColumnByUID(iuids[0], useSel);
  return ut_vector_max(tab);
}

double Db::getMean(const String& name, bool useSel) const
{
  VectorInt iuids = _ids(name, true);
  if (iuids.empty()) return TEST;
  VectorDouble tab = getColumnByUID(iuids[0], useSel);
  return ut_vector_mean(tab);
}

double Db::getVariance(const String& name, bool useSel) const
{
  VectorInt iuids = _ids(name, true);
  if (iuids.empty()) return TEST;
  VectorDouble tab = getColumnByUID(iuids[0], useSel);
  return ut_vector_var(tab);
}

double Db::getStdv(const String& name, bool useSel) const
{
  VectorInt iuids = _ids(name, true);
  if (iuids.empty()) return TEST;
  VectorDouble tab = getColumnByUID(iuids[0], useSel);
  return ut_vector_stdv(tab);
}

double Db::getCorrelation(const String& name1, const String& name2, bool useSel) const
{
  VectorInt iuids;
  iuids = _ids(name1, true);
  if (iuids.empty()) return TEST;
  VectorDouble tab1 = getColumnByUID(iuids[0], useSel);
  iuids = _ids(name2, true);
  if (iuids.empty()) return TEST;
  VectorDouble tab2 = getColumnByUID(iuids[0], useSel);
  return ut_vector_correlation(tab1, tab2);
}

int Db::getNDim() const
{
  return (_p.at(ELoc::X).getLocatorNumber());
}

bool Db::hasSameDimension(const Db* dbaux) const
{
  bool retOK = dbaux->getNDim() == getNDim();
  if (!retOK)
    messerr("The two Data bases should have the same Space Dimension");
  return retOK;
}

/**
 * Check if the Space Dimension of 'dbaux' is larger (or equal) than the one of 'this'
 * @param dbaux    Second Db
 * @return
 */
bool Db::hasLargerDimension(const Db* dbaux) const
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

double Db::getValueByColIdx(int iech, int icol) const
{
  if (!isColIdxValid(icol)) return TEST;
  return (_array[_getAddress(iech, icol)]);
}

VectorDouble Db::getValuesByColIdx(const VectorInt &iechs,
                                   const VectorInt &icols,
                                   bool bySample) const
{
  VectorDouble vec;

  if (bySample)
  {
    for (int j = 0; j < (int) iechs.size(); j++)
      for (int i = 0; i < (int) icols.size(); i++)
      {
        int iech = iechs[j];
        int icol = icols[i];
        if (!isColIdxValid(icol)) return VectorDouble();
        if (!isSampleIndexValid(iech)) return VectorDouble();
        vec.push_back(getValueByColIdx(iech, icol));
      }
  }
  else
  {
    for (int i = 0; i < (int) icols.size(); i++)
      for (int j = 0; j < (int) iechs.size(); j++)
      {
        int iech = iechs[j];
        int icol = icols[i];
        if (!isColIdxValid(icol)) return VectorDouble();
        if (!isSampleIndexValid(iech)) return VectorDouble();
        vec.push_back(getValueByColIdx(iech, icol));
      }
  }
  return vec;
}

void Db::setValueByColIdx(int iech, int icol, double value)
{
  if (!isColIdxValid(icol)) return;
  if (!isSampleIndexValid(iech)) return;
  _array[_getAddress(iech, icol)] = value;
}

void Db::setValuesByColIdx(const VectorInt &iechs,
                           const VectorInt &icols,
                           const VectorDouble &values,
                           bool bySample)
{
  if ((int) icols.size() * (int) iechs.size() != (int) values.size())
  {
    messerr("Dimensions of 'icols'(%d), 'iechs'(%d) and 'values'(%d) are inconsistent",
            (int) icols.size(), (int) iechs.size(), (int) values.size());
    return ;
  }

  int lec = 0;
  if (bySample)
  {
    for (int j = 0; j < (int) iechs.size(); j++)
      for (int i = 0; i < (int) icols.size(); i++)
      {
        int icol = icols[i];
        int iech = iechs[j];
        if (!isColIdxValid(icol)) return;
        if (!isSampleIndexValid(iech)) return;
        _array[_getAddress(iech, icol)] = values[lec++];
      }
  }
  else
  {
    for (int i = 0; i < (int) icols.size(); i++)
      for (int j = 0; j < (int) iechs.size(); j++)
      {
        int icol = icols[i];
        int iech = iechs[j];
        if (!isColIdxValid(icol)) return;
        if (!isSampleIndexValid(iech)) return;
        _array[_getAddress(iech, icol)] = values[lec++];
      }
  }
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
    if (getVariableNumber() != nvar)
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

/**
 * Check if the information (ELOC.Z) for a sample is isotopic or not
 * Isotopic means that all variables (for this sample) are defined
 * @param iech Rank of the sample
 * @param nvar_max Maximum number of variables to be checked (or -1)
 * @return
 *
 * @remark The returned answer is false is there is no variable defined
 * @remark or if the sample rank is not valid.
 * @remark If 'nvar-max' is defined, the test is performed on the 'nvar_max'
 * @remark first variables. Otherwise, it is performed on all ELOC.Z variables
 */
bool Db::isIsotopic(int iech, int nvar_max) const
{
  int nvar = getVariableNumber();
  if (nvar_max > 0) nvar = MIN(nvar, nvar_max);
  if (nvar <= 0) return false;
  if (!isSampleIndexValid(iech)) return false;

  for (int ivar = 0; ivar < nvar; ivar++)
    if (FFFF(getVariable(iech, ivar))) return false;
  return true;
}

bool Db::isAllUndefined(int iech) const
{
  int nvar = getVariableNumber();
  if (nvar <= 0) return false;
  if (!isSampleIndexValid(iech)) return false;

  for (int ivar = 0; ivar < nvar; ivar++)
    if (! FFFF(getVariable(iech, ivar))) return true;
  return false;
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
  int nech = getSampleNumber(useSel);
  VectorDouble vec(nech);
  VectorDouble vecl = getColumnByLocator(ELoc::L, item, useSel);
  VectorDouble vecu = getColumnByLocator(ELoc::U, item, useSel);

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

VectorDouble Db::getGradients(int item, bool useSel) const
{
  if (!hasGradient()) return VectorDouble();
  VectorDouble tab;

  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    if (useSel && ! isActive(iech)) continue;
    tab.push_back(getGradient(iech,item));
  }
  return tab;
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

VectorDouble Db::getTangents(int item, bool useSel) const
{
  if (!hasTangent()) return VectorDouble();
  VectorDouble tab;

  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    if (useSel && ! isActive(iech)) continue;
    tab.push_back(getTangent(iech,item));
  }
  return tab;
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

/**
 * Return the Selection value at Sample 'iech'
 * @param iech Sample number
 * @return
 * @remark If the selection value if TEST, the sample is considered as masked off.
 */
int Db::getSelection(int iech) const
{
  if (!hasSelection()) return 1;
  double value = getFromLocator(ELoc::SEL, iech, 0);
  if (FFFF(value)) return 0;
  int sel = (value != 0) ? 1 :  0;
  return (sel);
}

void Db::setSelection(int iech, int value)
{
  setFromLocator(ELoc::SEL, iech, 0, (value == 0) ? 0. : 1.);
}

bool Db::hasSelection() const
{
  return (getFromLocatorNumber(ELoc::SEL) > 0);
}

/**
 * Returns the number of active samples if a Selection is already defined.
 *
 * If no Selection is currently defined, it returns the total number of samples (see getSampleNumber())
 * @return Number of active samples
 *
 * @remark This method is deprecated and should be replaced by a call to
 * getSampleNumber()
 */
GSTLEARN_DEPRECATED int Db::getActiveSampleNumber() const
{
  if (!hasSelection()) return (getSampleNumber());

  /* Case when a selection is present */

  int count = 0;
  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    if (getFromLocator(ELoc::SEL, iech, 0) != 0) count++;
  }
  return count;
}

/**
 * Return the absolute rank of a sample from its relative rank
 * @param iech Relative rank
 * @return
 */
int Db::getActiveSampleRank(int iech) const
{
  int nech = getSampleNumber(false);
  int jech = 0;
  for (int i = 0; i < nech; i++)
  {
    if (! isActive(i)) continue;
    if (iech == jech) return i;
    jech++;
  }
  return -1;
}

/**
 * Return the absolute rank of a sample where the variable 'item' is defined
 * from its relative rank
 * @param iech Relative rank
 * @param item Rank of the item
 * @return
 */
int Db::getActiveAndDefinedSampleRank(int iech, int item) const
{
  int nech = getSampleNumber(false);
  int jech = 0;
  for (int i = 0; i < nech; i++)
  {
    if (! isActive(i)) continue;
    if (FFFF(getVariable(i, item))) continue;
    if (iech == jech) return i;
    jech++;
  }
  return -1;
}

/**
 * Returns the Number of samples
 * @param useSel When FALSE returns the total sample number.
 * When TRUE returns the number of active samples
 * @return
 */
int Db::getSampleNumber(bool useSel) const
{
  if (! hasSelection()) return _nech;

  if (! useSel)
    return _nech;
  else
  {
    int count = 0;
    for (int iech = 0; iech < getSampleNumber(); iech++)
    {
      if (getFromLocator(ELoc::SEL, iech, 0) != 0) count++;
    }
    return count;
  }
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
  if (hasWeight()) icol = getColIdxByLocator(ELoc::W, 0);

  int ecr = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    if (useSel && !sel.empty() && sel[iech] == 0) continue;
    if (icol >= 0)
      tab[ecr] = getValueByColIdx(iech, icol);
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
  VectorInt iuids = _ids(name, true);
  if (iuids.empty()) return 0;
  VectorDouble tab = getColumnByUID(iuids[0], false);

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
 * Returns the rank of (one of) the lastly added UID in the Db
 * @param number 0 designates the last, 1 the one before last...
 * @return
 */
int Db::getLastUID(int number) const
{
  VectorInt ranks;
  for (int i = 0; i < (int) _uidcol.size(); i++)
    if (_uidcol[i] >= 0) ranks.push_back(i);
  int size = static_cast<int> (ranks.size());
  if (number > size)
    return -1;
  else
    return ranks[size - number - 1];
}

String Db::getLastName(int number) const
{
  int iuid = getLastUID(number);
  String name = getNameByUID(iuid);
  return name;
}

int Db::_getLastColumn(int number) const
{
  if (number > _ncol)
    return -1;
  else
    return (_ncol - number);
}

String Db::getNameByLocator(const ELoc& locatorType, int locatorIndex) const
{
  int icol = getColIdxByLocator(locatorType, locatorIndex);
  if (icol < 0) return String();
  return _colNames[icol];
}

String Db::getNameByColIdx(int icol) const
{
  if (! isColIdxValid(icol)) return String();
  return _colNames[icol];
}

String Db::getNameByUID(int iuid) const
{
  int icol = getColIdxByUID(iuid);
  if (icol < 0) return ("");
  return getNameByColIdx(icol);
}

VectorString Db::getNamesByLocator(const ELoc& locatorType) const
{
  VectorString namelist;
  if (!isLocatorTypeValid(locatorType)) return namelist;
  int count = getFromLocatorNumber(locatorType);
  for (int i = 0; i < count; i++)
  {
    int icol = getColIdxByLocator(locatorType, i);
    namelist.push_back(getNameByColIdx(icol));
  }
  return namelist;
}

VectorString Db::getNamesByColIdx(const VectorInt& icols) const
{
  VectorString namelist;
  for (int i = 0; i < (int) icols.size(); i++)
  {
    int icol = icols[i];
    namelist.push_back(_colNames[icol]);
  }
  return namelist;
}

VectorString Db::getNamesByUID(const VectorInt& iuids) const
{
  VectorString namelist;
  if (iuids.empty()) return namelist;
  int count = static_cast<int> (iuids.size());
  for (int i = 0; i < count; i++)
  {
    int icol = getColIdxByUID(iuids[i]);
    namelist.push_back(getNameByColIdx(icol));
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

VectorString Db::getAllNames() const
{
  VectorString names = _colNames;
  return names;
}

void Db::_setNameByColIdx(int icol, const String& name)
{
  if (!isColIdxValid(icol)) return;
  _colNames[icol] = name;
  correctNewNameForDuplicates(_colNames, icol);
}

void Db::setNameByUID(int iuid, const String& name)
{
  int icol = getColIdxByUID(iuid);
  if (icol < 0) return;
  _colNames[icol] = name;
  correctNewNameForDuplicates(_colNames, icol);
}

void Db::setNameByColIdx(int icol, const String& name)
{
  if (! isColIdxValid(icol)) return;
  _colNames[icol] = name;
}

void Db::setName(const String& old_name, const String& name)
{
  int icol = getColIdx(old_name);
  if (icol < 0) return;
  _colNames[icol] = name;
  correctNewNameForDuplicates(_colNames, icol);
}

void Db::setName(const VectorString list, const String& name)
{
  for (int i = 0; i < (int) list.size(); i++)
  {
    int icol = getColIdx(list[i]);
    if (icol < 0) continue;
    _colNames[icol] = incrementStringVersion(name, i + 1);
  }
  correctNamesForDuplicates(_colNames);
}

void Db::setNameByLocator(const ELoc& locatorType, const String& name)
{
  VectorString namelist;
  if (!isLocatorTypeValid(locatorType)) return;
  int count = getFromLocatorNumber(locatorType);
  for (int i = 0; i < count; i++)
  {
    int icol = getColIdxByLocator(locatorType, i);
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
  sstr << "Number of Columns            = " << getColumnNumber() << std::endl;
  sstr << "Maximum Number of UIDs       = " << getUIDMaxNumber()
       << std::endl;
  sstr << "Total number of samples      = " << getSampleNumber() << std::endl;
  if (hasSelection())
    sstr << "Number of active samples     = " << getSampleNumber(true)
         << std::endl;
  return sstr.str();
}

String Db::_summaryExtensions(void) const
{
  std::stringstream sstr;
  int ndim = getNDim();
  if (ndim <= 0) return sstr.str();

  /* Printout */

  sstr << toTitle(1, "Data Base Extension");
  for (int idim = 0; idim < ndim; idim++)
  {
    VectorDouble coor = getCoordinates(idim, true);
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

String Db::_summaryVariables(void) const
{
  std::stringstream sstr;

  if (getColumnNumber() <= 0) return sstr.str();
  sstr << toTitle(1, "Variables");

  for (int icol = 0; icol < getColumnNumber(); icol++)
  {
    sstr << "Column = " << icol;
    sstr << " - Name = " << getNameByColIdx(icol);
    sstr << " - Locator = " << _getLocatorNameByColIdx(icol);
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
String Db::_summaryStats(VectorInt cols, int mode, int maxNClass) const
{
  std::stringstream sstr;

  int ncol = (cols.empty()) ? getColumnNumber() : static_cast<int> (cols.size());
  if (ncol <= 0) return sstr.str();

  sstr << toTitle(1, "Data Base Statistics");

  int nval, nmask, ntest, nout;
  double vmin, vmax, delta, mean, stdv;
  int nech = getSampleNumber();
  VectorDouble tab, wgt;

  // Loop on the columns

  for (int jcol = 0; jcol < ncol; jcol++)
  {
    int icol = (cols.empty()) ? jcol :
                                cols[jcol];
    if (!isColIdxValid(icol)) continue;

    tab = getColumnByColIdx(icol, true);
    wgt = getWeight(true);

    ut_statistics(static_cast<int> (tab.size()), tab.data(), NULL,
                  wgt.data(), &nval, &vmin, &vmax,
                  &delta, &mean, &stdv);

    sstr << icol + 1 << " - Name " << getNameByColIdx(icol) << " - Locator "
         << _getLocatorNameByColIdx(icol) << std::endl;
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

String Db::_summaryArrays(VectorInt cols, bool useSel) const
{
  std::stringstream sstr;

  int ncol = (cols.empty()) ? getColumnNumber() : static_cast<int> (cols.size());
  if (ncol <= 0) return sstr.str();

  sstr << toTitle(1, "Data Base Contents");

  int number = getSampleNumber(useSel);

  VectorDouble tab;
  VectorString colnames;
  for (int jcol = 0; jcol < ncol; jcol++)
  {
    int icol = (cols.empty()) ? jcol : cols[jcol];
    if (!isColIdxValid(icol)) continue;
    VectorDouble local = getColumnByColIdx(icol, useSel);
    tab.insert(tab.end(), local.begin(), local.end());
    colnames.push_back(getNameByColIdx(icol));
  }

  sstr << toMatrix(String(), colnames, VectorString(), 1, ncol, number, tab);

  return sstr.str();
}

String Db::_toStringCommon(const AStringFormat *strfmt) const
{
  std::stringstream sstr;
  static int MAX_NCLASS = 50;

  const DbStringFormat* dbfmt = dynamic_cast<const DbStringFormat*>(strfmt);
  DbStringFormat dsf;
  if (dbfmt != nullptr) dsf = *dbfmt;

  // Possibly convert 'names' into 'cols'

  VectorInt cols = dsf.getCols();
  if (cols.empty())
  {
    VectorInt iuids = _ids(dsf.getNames(), false);
    if (! iuids.empty()) cols = getColIdxsByUID(iuids);
  }

  /* Print the Extension */

  if (dsf.matchExtend())
    sstr << _summaryExtensions();

  /* Print the statistics */

  if (dsf.matchStats())
    sstr << _summaryStats(cols, dsf.getMode(), MAX_NCLASS);

  /* Print the contents of the Data Base */

  if (dsf.matchArray())
    sstr << _summaryArrays(cols, dsf.getUseSel());

  /* Print the list of variables */

  if (dsf.matchVars())
    sstr << _summaryVariables();

  /* Print the locators */

  if (dsf.matchLocator())
  {
    sstr << _summaryUIDs() << std::endl;
    sstr << _summaryLocators() << std::endl;
  }
  return sstr.str();
}

String Db::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  const DbStringFormat* dbfmt = dynamic_cast<const DbStringFormat*>(strfmt);
  DbStringFormat dsf;
  if (dbfmt != nullptr) dsf = *dbfmt;

  sstr << toTitle(0, "Data Base Characteristics");

  if (dsf.matchResume())
    sstr << _summaryString();

  sstr << _toStringCommon(&dsf);

  return sstr.str();
}

VectorDouble Db::getSelection(void) const
{
  int nech = getSampleNumber();
  VectorDouble tab;

  if (!hasSelection()) return tab;
  int icol = getColIdxByLocator(ELoc::SEL,0);
  if (!isColIdxValid(icol)) return tab;

  tab.resize(nech);
  for (int iech = 0; iech < nech; iech++)
    tab[iech] = getValueByColIdx(iech, icol);
  return tab;
}

/**
 * Returns the contents of a Column of the Db refered by its column index
 * @param icol Column index (from [0,n[)
 * @param useSel Is the selection taken into account
 * @return The vector of values
 */
VectorDouble Db::getColumnByColIdx(int icol, bool useSel) const
{
  int nech = getSampleNumber();
  if (!isColIdxValid(icol)) return VectorDouble();

  VectorDouble tab(nech, TEST);
  VectorDouble sel;
  if (useSel) sel = getSelection();

  int ecr = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    bool defined = true;
    if (useSel && !sel.empty()) defined = (sel[iech] == 1);
    if (!defined) continue;
    tab[ecr] = getValueByColIdx(iech, icol);
    ecr++;
  }
  tab.resize(ecr);
  return tab;
}

VectorDouble Db::getColumnByUID(int iuid, bool useSel) const
{
  int icol = getColIdxByUID(iuid);
  if (icol < 0) return VectorDouble();
  return getColumnByColIdx(icol, useSel);
}

VectorDouble Db::getColumnByLocator(const ELoc& locatorType,
                                   int locatorIndex,
                                   bool useSel) const
{
  int icol = getColIdxByLocator(locatorType, locatorIndex);
  if (icol < 0) return VectorDouble();
  return getColumnByColIdx(icol, useSel);
}

VectorDouble Db::getColumn(const String& name, bool useSel) const
{
  VectorInt iuids = _ids(name, true);
  if (iuids.empty()) return VectorDouble();
  int icol = getColIdxByUID(iuids[0]);
  if (icol < 0) return VectorDouble();
  return getColumnByColIdx(icol, useSel);
}

VectorDouble Db::getColumnsByLocator(const ELoc& locatorType, bool useSel) const
{
  VectorString names = getNamesByLocator(locatorType);
  return getColumns(names, useSel);
}

VectorDouble Db::getColumnsByUID(const VectorInt& iuids, bool useSel) const
{
  if (iuids.empty()) return VectorDouble();
  int nech = getSampleNumber(useSel);
  int nvar = static_cast<int> (iuids.size());
  VectorDouble retval(nvar * nech);

  /* Loop on the variables to be retrieved */

  int ecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    VectorDouble local = getColumnByUID(iuids[ivar], useSel);
    for (int iech = 0; iech < nech; iech++)
      retval[ecr++] = local[iech];
  }
  return retval;
}

VectorDouble Db::getColumnsByColIdx(const VectorInt& icols, bool useSel) const
{
  int nech = getSampleNumber();
  int nvar = static_cast<int> (icols.size());
  VectorDouble retval(nvar * nech);

  /* Loop on the variables to be retrieved */

  int ecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    VectorDouble local = getColumnByColIdx(icols[ivar], useSel);
    if (local.empty()) continue;
    for (int iech = 0; iech < nech; iech++)
      retval[ecr++] = local[iech];
  }
  return retval;
}

VectorDouble Db::getColumnsByColIdxInterval(int icol_beg, int icol_end, bool useSel) const
{
  VectorInt icols;
  for (int icol = icol_beg; icol < icol_end; icol++)
    icols.push_back(icol);
  return getColumnsByColIdx(icols, useSel);
}

VectorDouble Db::getColumnsByUIDRange(int iuid_beg, int iuid_end, bool useSel) const
{
  VectorInt iuids;
  for (int iuid = iuid_beg; iuid < iuid_end; iuid++)
    iuids.push_back(iuid);
  return getColumnsByUID(iuids, useSel);
}

VectorDouble Db::_getItem(const String& exp_name,
                          bool useSel,
                          const VectorInt& rows) const
{
  int nrows = (int) rows.size();
  VectorDouble local(nrows);

  // Read the whole column of values through possible selection
  VectorDouble allvec = getColumn(exp_name, useSel);

  // Shrink the values for the retained rows only
  for (int irow = 0; irow < nrows; irow++)
    local[irow] = allvec[rows[irow]];

  return local;
}

void Db::_setItem(const String& name,
                  const VectorInt& rows,
                  const VectorDouble& values)
{
  int icol = getUID(name);
  for (int jjrow = 0; jjrow < (int) rows.size(); jjrow++)
  {
    int jrow = rows[jjrow];
    setArray(jrow, icol, values[jjrow]);
  }
}

void Db::_setItem(const String& name,
                  bool useSel,
                  const VectorDouble& values)
{
  int icol = getUID(name);
  int nrows  = getSampleNumber();
  int jjrow = 0;
  for (int jrow = 0; jrow < nrows; jrow++)
  {
    if (useSel && ! isActive(jrow)) continue;
    setArray(jjrow, icol, values[jrow]);
    jjrow++;
  }
}

bool Db::_isValidCountRows(const VectorInt& rows,
                           bool useSel,
                           const VectorDouble& values) const
{
  if (rows.empty()) return false;
  if (! isSampleIndicesValid(rows, useSel)) return false;
  if (rows.size() != values.size())
  {
    messerr("Mismatch in dimensions:");
    messerr("- From 'values' = %d",(int) values.size());
    messerr("- From 'rows' = %d",(int) rows.size());
    return false;
  }
  return true;
}

bool Db::_isValidCountRows(bool useSel, const VectorDouble& values) const
{
  int nrows = getSampleNumber(useSel);
  if (nrows != (int) values.size())
  {
    messerr("Mismatch in dimensions:");
    messerr("- From 'values' = %d",(int) values.size());
    messerr("- From 'rows' = %d",nrows);
    return false;
  }
  return true;
}

VectorString Db::_getVarNames(const VectorString& colnames,
                              int expectedVarCount)
{
  VectorString exp_names;

  if (colnames.empty()) return exp_names;
  int number = (int) colnames.size();

  // Constitute the output list by gluing expanded parts

  for (int i = 0; i < number; i++)
  {
    VectorString sublist = expandNameList(colnames[i]);

    if (sublist.empty())
    {
      // sublist is empty: the variable must be created

      (void) addColumnsByConstant(1, TEST, colnames[i]);
      exp_names.push_back(colnames[i]);
    }
    else
    {

      // Glue the expanded variable names to the output list

      exp_names.insert(exp_names.end(),
                       sublist.begin(), sublist.end());
    }
  }

  int current = (int) exp_names.size();
  if (current > expectedVarCount)
  {
    messerr("Mismatch between dimension of 'values'(%d) and variable list");
    for (int i =0; i < current; i++)
      messerr("- %s",exp_names[i].c_str());
    return VectorString();
  }
  if (current < expectedVarCount)
  {
    // Complete the variable list by duplicating the last variable

    int missing = expectedVarCount - current;
    VectorString sublist = generateMultipleNames(colnames[number-1],
                                                 missing);
    exp_names.insert(exp_names.end(),
                     sublist.begin(), sublist.end());
  }

  return exp_names;
}

VectorVectorDouble Db::getItem(const VectorInt& rows,
                               const VectorString& colnames,
                               bool useSel) const
{
  VectorVectorDouble values;

  if (! isSampleIndicesValid(rows, useSel)) return values;
  if (rows.empty()) return values;
  VectorString exp_names = expandNameList(colnames);
  if (exp_names.empty()) return values;

  for (int icol = 0; icol < (int) exp_names.size(); icol++)
  {
    VectorDouble local = _getItem(exp_names[icol],useSel,rows);
    values.push_back(local);
  }
  return values;
}

VectorVectorDouble Db::getItem(const VectorInt& rows,
                               const String& colname,
                               bool useSel) const
{
  VectorVectorDouble values;

  if (! isSampleIndicesValid(rows, useSel)) return values;
  if (rows.empty()) return values;
  VectorString exp_names = expandNameList(colname);
  if (exp_names.empty()) return values;

  for (int icol = 0; icol < (int) exp_names.size(); icol++)
  {
    VectorDouble local = _getItem(exp_names[icol],useSel,rows);
    values.push_back(local);
  }
  return values;
}

VectorVectorDouble Db::getItem(const VectorInt& rows,
                               const ELoc& locatorType,
                               bool useSel) const
{
  VectorVectorDouble values;

  if (! isSampleIndicesValid(rows, useSel)) return values;
  if (rows.empty()) return values;
  VectorString exp_names = getNamesByLocator(locatorType);
  if (exp_names.empty()) return values;

  for (int icol = 0; icol < (int) exp_names.size(); icol++)
  {
    VectorDouble local = _getItem(exp_names[icol],useSel,rows);
    values.push_back(local);
  }
  return values;
}

VectorVectorDouble Db::getItem(const VectorString& colnames, bool useSel) const
{
  VectorVectorDouble values;

  VectorString exp_names = expandNameList(colnames);
  if (exp_names.empty()) return values;

  for (int icol = 0; icol < (int) exp_names.size(); icol++)
  {
    VectorDouble local = getColumn(exp_names[icol], useSel);
    values.push_back(local);
  }
  return values;
}

VectorVectorDouble Db::getItem(const String& colname, bool useSel) const
{
  VectorVectorDouble values;

  VectorString exp_names = expandNameList(colname);
  if (exp_names.empty()) return values;

  for (int icol = 0; icol < (int) exp_names.size(); icol++)
  {
    VectorDouble local = getColumn(exp_names[icol], useSel);
    values.push_back(local);
  }
  return values;
}

VectorVectorDouble Db::getItem(const ELoc& locatorType, bool useSel) const
{
  VectorVectorDouble values;

  VectorString exp_names = getNamesByLocator(locatorType);
  if (exp_names.empty()) return values;

  for (int icol = 0; icol < (int) exp_names.size(); icol++)
  {
    VectorDouble local = getColumn(exp_names[icol], useSel);
    values.push_back(local);
  }
  return values;
}

VectorString Db::getItemNames(const VectorString& colnames)
{
  return expandNameList(colnames);
}
VectorString Db::getItemNames(const String& colname)
{
  return expandNameList(colname);
}
VectorString Db::getItemNames(const ELoc& locatorType)
{
  return getNamesByLocator(locatorType);
}

int Db::setItem(const VectorInt& rows,
                const VectorString& colnames,
                const VectorVectorDouble& values,
                bool useSel)
{
  if (! _isValidCountRows(rows, useSel, values[0])) return 1;
  int expectedVarCount = (int) values.size();
  VectorString exp_names = _getVarNames(colnames, expectedVarCount);
  if (exp_names.empty()) return 1;

  for (int icol = 0; icol < (int) exp_names.size(); icol++)
    _setItem(exp_names[icol], rows, values[icol]);
  return 0;
}

int Db::setItem(const VectorInt& rows,
                const ELoc& locatorType,
                const VectorVectorDouble& values,
                bool useSel)
{
  if (! _isValidCountRows(rows, useSel, values[0])) return 1;
  VectorString colnames = getNamesByLocator(locatorType);
  int expectedVarCount = (int) values.size();
  VectorString exp_names = _getVarNames(colnames, expectedVarCount);
  if (exp_names.empty()) return 1;

  for (int icol = 0; icol < (int) exp_names.size(); icol++)
    _setItem(exp_names[icol], rows, values[icol]);
  return 0;
}

int Db::setItem(const VectorInt& rows,
                const String& colname,
                const VectorDouble& values,
                bool useSel)
{
  if (! _isValidCountRows(rows, useSel, values)) return 1;
  VectorString colnames(1);
  colnames[0] = colname;
  VectorString exp_names = _getVarNames(colnames, 1);
  if (exp_names.empty()) return 1;

  _setItem(exp_names[0], rows, values);
  return 0;
}

int Db::setItem(const VectorString& colnames,
                const VectorVectorDouble& values,
                bool useSel)
{
  if (! _isValidCountRows(useSel, values[0])) return 1;
  int expectedVarCount = (int) values.size();
  VectorString exp_names = _getVarNames(colnames, expectedVarCount);
  if (exp_names.empty()) return 1;

  for (int icol = 0; icol < (int) exp_names.size(); icol++)
    _setItem(exp_names[icol], useSel, values[icol]);
  return 0;
}

int Db::setItem(const ELoc& locatorType,
                const VectorVectorDouble& values,
                bool useSel)
{
  if (! _isValidCountRows(useSel, values[0])) return 1;
  VectorString colnames = getNamesByLocator(locatorType);
  int expectedVarCount = (int) values.size();
  VectorString exp_names = _getVarNames(colnames, expectedVarCount);
  if (exp_names.empty()) return 1;

  for (int icol = 0; icol < (int) exp_names.size(); icol++)
    _setItem(exp_names[icol], useSel, values[icol]);
  return 0;
}

int Db::setItem(const String& colname,
                const VectorDouble& values,
                bool useSel)
{
  if (! _isValidCountRows(useSel, values)) return 1;
  VectorString colnames(1);
  colnames[0] = colname;
  VectorString exp_names = _getVarNames(colnames, 1);
  if (exp_names.empty()) return 1;

  _setItem(colname, useSel, values);
  return 0;
}

VectorDouble Db::getAllColumns(bool useSel) const
{
  VectorInt iuids = getAllUIDs();
  return getColumnsByUID(iuids, useSel);
}

void Db::setAllColumns(const VectorVectorDouble& tabs,bool useSel)
{
  VectorInt iuids = getAllUIDs();
  for (int iuid = 0; iuid < (int) iuids.size(); iuid++)
    setColumnByUID(tabs[iuid], iuids[iuid], useSel);
}

VectorDouble Db::getColumns(const VectorString& names, bool useSel) const
{
  if (names.empty()) return VectorDouble();
  VectorInt iuids =  _ids(names, false);
  return getColumnsByUID(iuids, useSel);
}

VectorDouble Db::getFFFFs(const VectorString& names, bool useSel) const
{
  if (names.empty()) return VectorDouble();
  VectorInt iuids =  _ids(names, false);
  return getColumnsByUID(iuids, useSel);
}

/**
 * Returns the vector of coordinates along a given Space Dimension
 * @param idim    Rank of the Space dimension
 * @param useSel  Use the Data Selection
 * @param flag_rotate Flag for rotation (only for Grid)
 * @return
 */
VectorDouble Db::getCoordinates(int idim, bool useSel, bool flag_rotate) const
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
 * Returns the rank of the Single Column corresponding to 'name'
 * @param name Named for the searched column
 * @return The rank of the Single column or -1
 */
int Db::getColIdx(const String& name) const
{
  VectorString exp_name = expandNameList(name);
  if (exp_name.empty()) return -1;
  return getRankInList(_colNames, exp_name[0]);
}

VectorInt Db::getColIdxs(const String& name) const
{
  VectorString exp_names = expandNameList(name);
  return getColIdxs(exp_names);
}

VectorInt Db::getColIdxs(const VectorString& names) const
{
  VectorString exp_names = expandNameList(names);
  if (exp_names.size() <= 0) return VectorInt();
  int number = static_cast<int> (exp_names.size());
  VectorInt icols(number);
  for (int i = 0; i < number; i++)
    icols[i] = getColIdx(exp_names[i]);
  return icols;
}

VectorInt Db::getColIdxsByLocator(const ELoc& locatorType) const
{
  VectorInt icols;
  if (!isLocatorTypeValid(locatorType)) return icols;
  int number = getLocatorNumber(locatorType);
  if (number <= 0) return icols;

  icols.resize(number);
  for (int i = 0; i < number; i++)
    icols[i] = getColIdxByLocator(locatorType, i);
  return icols;
}

/**
 * Returns the Single UID which corresponds to the searched name
 * @param name Name to be searched for
 * @return Rank of the UID or -1
 */
int Db::getUID(const String& name) const
{
  VectorInt iuids = _ids(name, true);
  if (iuids.empty()) return -1;
  int icol = getColIdxByUID(iuids[0]);
  return getUIDByColIdx(icol);
}

/**
 * This is a BASIC function returning the vector of ranks of the UID
 * which corresponds to a set of existing names
 */
VectorInt Db::_getUIDsBasic(const VectorString& names) const
{
  if (names.empty()) return VectorInt();

  VectorInt iuids(names.size());
  for (unsigned int i = 0; i < names.size(); i++)
  {
    int icol = getRankInList(_colNames, names[i]);
    iuids[i] = getUIDByColIdx(icol);
  }
  return iuids;
}

VectorInt Db::getUIDs(const VectorString& names) const
{
  if (names.empty()) return VectorInt();

  VectorInt iuids(names.size());
  for (unsigned int i = 0; i < names.size(); i++)
    iuids[i] = getUID(names[i]);
  return iuids;
}

VectorInt Db::getUIDsByLocator(const ELoc& locatorType) const
{
  VectorInt iuids;
  if (!isLocatorTypeValid(locatorType)) return iuids;
  int number = getLocatorNumber(locatorType);
  if (number <= 0) return iuids;
  iuids.resize(number);
  for (int i = 0; i < number; i++)
    iuids[i] = getUIDByLocator(locatorType, i);
  return iuids;
}

VectorInt Db::getUIDsByColIdx(const VectorInt& icols) const
{
  VectorInt iuids;
  for (int i = 0; i < (int) icols.size(); i++)
    iuids.push_back(getUIDByColIdx(icols[i]));
  return iuids;
}

VectorInt Db::getAllUIDs() const
{
  VectorInt iuids;
  for (int i = 0; i < (int) _uidcol.size(); i++)
    if (_uidcol[i] >= 0) iuids.push_back(i);
  return iuids;
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
        setValueByColIdx(iech, jcol, tab[icol + ntab * iech]);
      else
        setValueByColIdx(iech, jcol, tab[ecr]);
    }
  }

  // Set the names
  _defineDefaultNames(shift, names);

  // Set the locators
  _defineDefaultLocators(shift, locatorNames);

  return;
}

void Db::generateRank(const String& radix)
{
  int nech = getSampleNumber();
  VectorDouble vec(nech);
  for (int iech = 0; iech < nech; iech++)
    vec[iech] = iech + 1;

  (void) addColumns(vec, radix);
}

/**
 * Paint the column 'icol' with sample rank
 * @param icol Rank of the column to be painted
 */
void Db::_createRank(int icol)
{
  int nech = getSampleNumber();
  for (int iech = 0; iech < nech; iech++)
    setArray(iech, icol, iech + 1);

  // Set the name

  _setNameByColIdx(icol, "rank");

  // Set the locators (No particular action for the Rank)
}

void Db::_defineDefaultNames(int shift, const VectorString& names)
{
  int ncol = getColumnNumber() - shift;
  if (!names.empty())
  {
    if ((int) names.size() != ncol) throw("Error in the dimension of 'names'");
  }

  for (int icol = 0; icol < ncol; icol++)
  {
    if (!names.empty())
      _setNameByColIdx(icol + shift, names[icol]);
    else
      _setNameByColIdx(icol + shift, incrementStringVersion("New", icol + 1));
  }
}

void Db::_defineDefaultLocators(int shift, const VectorString& locatorNames)
{
  if (locatorNames.empty()) return;

  int ncol = getColumnNumber() - shift;
  if ((int) locatorNames.size() != ncol)
    throw("Error in the dimension of 'locatorNames'");

  ELoc locatorType;
  int locatorIndex, mult;
  for (int icol = 0; icol < ncol; icol++)
  {
    if (!locatorIdentify(locatorNames[icol], &locatorType, &locatorIndex, &mult))
      setLocatorByUID(icol + shift, locatorType, locatorIndex);
  }
}

void Db::_defineDefaultLocatorsByNames(int shift, const VectorString& names)
{
  if (names.empty()) return;

  int ncol = getColumnNumber() - shift;
  if ((int) names.size() != ncol) throw("Error in the dimension of 'names'");

  ELoc locatorType;
  int locatorIndex, mult;
  for (int icol = 0; icol < ncol; icol++)
  {
    if (!locatorIdentify(names[icol], &locatorType, &locatorIndex, &mult))
      setLocatorByUID(icol + shift, locatorType, locatorIndex);
  }
}

/**
 * Return the monovariate statistics on variables given their UID
 * @param iuids              List of UID of the variables
 * @param opers              List of operations to be performed on the variables at each point.
 * @param flagIso  
 * @param flagVariableWise   If False, the new variable is added to Db
 * @param flagPrint          If True (default), the statistics are printed. If False, statistics are returned in a vector.
 * @param proba              For 'quant': the quantile for this probability is calculated 
 * @param vmin               For 'prop', 'T', 'Q', 'M', 'B': defines the lower bound of the interval to work in
 * @param vmax               For 'prop', 'T', 'Q', 'M', 'B': defines the upper bound of the interval to work in
 * @param title              If flagPrint is True, the title of the printed statistics.
 * @param namconv            Naming Convention
 * @return If flagPrint is False, returns a vector containing the statistics. If there is more than one operator and more than one variable,
 * the statistics are ordered first by variables (all the statistics of the first variable, then all the stats of the second variable...).
 * 
 * @see statisticsMulti (multivariate statistics)
 */
VectorDouble Db::statisticsByUID(const VectorInt& iuids,
                                 const std::vector<EStatOption>& opers,
                                 bool flagIso,
                                 bool flagVariableWise,
                                 bool flagPrint,
                                 double proba,
                                 double vmin,
                                 double vmax,
                                 const String& title,
                                 const NamingConvention& namconv)
{
  VectorDouble stats;

  if (iuids.empty()) return stats;

  int noper = static_cast<int> (opers.size());
  if (noper <= 0) return stats;

  // Add the variables for PointWise statistics
  if (!flagVariableWise)
  {
    int iuidn = addColumnsByConstant(noper);
    if (iuidn < 0) return VectorDouble();

    dbStatisticsVariables(this, iuids, opers, iuidn, proba, vmin, vmax);

    namconv.setNamesAndLocators(this, iuidn);
    for (int i = 0; i < noper; i++)
    {
      EStatOption oper = opers[i];
      namconv.setNamesAndLocators(this, iuidn + i, oper.getKey());
    }
    return VectorDouble();
  }
  else
  {
    stats = dbStatisticsMono(this, iuids, opers, flagIso, proba, vmin, vmax);

    if (flagPrint)
    {
      VectorString varnames = getNamesByUID(iuids);
      messageFlush(statisticsMonoPrint(stats, opers, varnames, title));
      return VectorDouble();
    }
  }
  return stats;
}

/**
 * Return the monovariate statistics on variables given their names.
 * @param names              List of names of the variables
 * @param opers              List of operations to be performed on the variables at each point.
 * @param flagIso            If True, perform statistics on isotopic samples only
 * @param flagVariableWise   If False, results are added to Db; otherwise they are produced in output vector
 * @param flagPrint          If True (default), the statistics are printed. If False, statistics are returned in a vector.
 * @param proba              For 'quant': the quantile for this probability is calculated
 * @param vmin               For 'prop', 'T', 'Q', 'M', 'B': defines the lower bound of the interval to work in
 * @param vmax               For 'prop', 'T', 'Q', 'M', 'B': defines the upper bound of the interval to work in
 * @param title              If flagPrint is True, the title of the printed statistics. 
 * @param namconv            Naming Convention
 * @return If flagPrint is False, returns a vector containing the statistics. If there is more than one operator and more than one variable,
 * the statistics are ordered first by variables (all the statistics of the first variable, then all the stats of the second variable...).
 * 
 * @see statisticsMulti (multivariate statistics)
 */
VectorDouble Db::statistics(const VectorString& names,
                            const std::vector<EStatOption>& opers,
                            bool flagIso,
                            bool flagVariableWise,
                            bool flagPrint,
                            double proba,
                            double vmin,
                            double vmax,
                            const String& title,
                            const NamingConvention& namconv)
{
  VectorInt iuids = _ids(names, false);
  if (iuids.empty()) return VectorDouble();
  return statisticsByUID(iuids, opers, flagIso, flagVariableWise, flagPrint,
                         proba, vmin, vmax, title, namconv);
}

VectorDouble Db::statisticsByLocator(const ELoc& locatorType,
                                     const std::vector<EStatOption>& opers,
                                     bool flagIso,
                                     bool flagVariableWise,
                                     bool flagPrint,
                                     double proba,
                                     double vmin,
                                     double vmax,
                                     const String& title,
                                     const NamingConvention& namconv)
{
  VectorInt iuids = getUIDsByLocator(locatorType);
  if (iuids.empty()) return VectorDouble();
  return statisticsByUID(iuids, opers, flagIso, flagVariableWise, flagPrint,
                         proba, vmin, vmax, title, namconv);
}

VectorDouble Db::statisticsMultiByUID(const VectorInt& iuids,
                                      bool flagIso,
                                      bool flagPrint,
                                      const String& title)
{
  VectorDouble stats;

  if (iuids.empty()) return stats;

  stats = dbStatisticsMulti(this, iuids, flagIso);

  if (flagPrint)
  {
    VectorString varnames = getNamesByUID(iuids);
    messageFlush(statisticsMultiPrint(stats, varnames, title));
  }
  return stats;
}

VectorDouble Db::statisticsMulti(const VectorString& names,
                                 bool flagIso,
                                 bool flagPrint,
                                 const String& title)
{
  VectorInt iuids = _ids(names, false);
  if (iuids.empty()) return VectorDouble();

  return statisticsMultiByUID(iuids, flagIso, flagPrint, title);
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

Db* Db::createFromNF(const String& neutralFilename, bool verbose)
{
  Db* db = nullptr;
  std::ifstream is;
  db = new Db;
  if (db->_fileOpenRead(neutralFilename, is, verbose))
  {
    if (! db->deserialize(is, verbose))
    {
      delete db;
      db = nullptr;
    }
  }
  else
  {
    delete db;
    db = nullptr;
  }
  return db;
}

bool Db::_serialize(std::ostream& os,bool /*verbose*/) const
{
  int ncol = getColumnNumber();
  VectorString locators = getLocators(1);
  VectorString names = getNames("*");
  bool ret = true;

  ret = ret && _recordWrite<int>(os, "Number of variables", ncol);
  ret = ret && _recordWrite<int>(os, "Number of samples", getSampleNumber());
  ret = ret && _recordWriteVec<String>(os, "Locators", locators);
  ret = ret && _recordWriteVec<String>(os, "Names", names);
  ret = ret && _commentWrite(os, "Array of values");
  for (int iech = 0; ret && iech < getSampleNumber(); iech++)
  {
    VectorDouble vals = getArrayBySample(iech);
    ret = ret && _recordWriteVec(os, "", vals);
  }
  return ret;
}

bool Db::_deserialize(std::istream& is, bool /*verbose*/)
{
  int ncol = 0, nrow = 0, nech = 0;
  VectorString locators;
  VectorString names;
  VectorDouble values;
  VectorDouble allvalues;

  // Read the file
  bool ret = true;
  ret = ret && _recordRead<int>(is, "Number of variables", ncol);
  ret = ret && _recordRead<int>(is, "Number of samples", nrow);
  if (! ret) return ret;
  if (ncol > 0)
  {
    ret = ret && _recordReadVec<String>(is, "Locators", locators, ncol);
    ret = ret && _recordReadVec<String>(is, "Names", names, ncol);
  }
  if (! ret) return ret;

  for (int iech = 0; iech < nrow; iech++)
  {
    ret = ret && _recordReadVec<double>(is, "Array of values", values, ncol);
    if (ret)
    {
      // Concatenate values by samples
      allvalues.insert(allvalues.end(), values.begin(), values.end());
      nech++;
    }
  }
  ret = (nech == nrow);

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

  // Initialize the Db
  resetDims(ncol, nech);

  // Load the values
  _loadData(ELoadBy::SAMPLE, 0, allvalues);

  // Update the column names and locators
  for (int i = 0; i < ncol; i++)
  {
    setNameByUID(i, names[i]);
    setLocatorByUID(i, tabloc[i], tabnum[i]);
  }
  return ret;
}

void Db::_loadData(const ELoadBy& order, int flag_add_rank, const VectorDouble& tab)
{
  // Preliminary check

  if (getColumnNumber() <= 0) return;
  int jcol = 0;

  // Add the rank (optional)

  if (flag_add_rank)
  {
    for (int iech = 0; iech < getSampleNumber(); iech++)
      setValueByColIdx(iech, jcol, iech + 1);
    setNameByUID(jcol, "rank");
    jcol++;
  }

  // Add the input array 'tab' (if provided)

  if (tab.empty()) return;
  int ntab = (flag_add_rank) ? getColumnNumber() - 1 : getColumnNumber();
  int ecr = 0;
  for (int icol = 0; icol < ntab; icol++)
  {
    for (int iech = 0; iech < getSampleNumber(); iech++, ecr++)
    {
      if (order == ELoadBy::SAMPLE)
        setValueByColIdx(iech, jcol, tab[icol + ntab * iech]);
      else
        setValueByColIdx(iech, jcol, tab[ecr]);
    }
    jcol++;
  }
  return;
}

bool Db::_isCountValid(const VectorInt iuids, bool flagOne, bool verbose) const
{
  if (iuids.empty() && flagOne)
  {
    if (verbose) messerr("No variable name corresponding to your criterion");
    return false;
  }
  else
  {
    if (iuids.size() > 1 && flagOne)
    {
      if (verbose)
      {
        messerr("You wanted to designate a SINGLE variable.");
        messerr("There are several variables matching your criterion:");
        for (unsigned int i = 0; i < iuids.size(); i++)
          messerr("- %s", getNameByUID(iuids[i]).c_str());
      }
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

VectorBool Db::getMaskArray() const
{
  int nech = getSampleNumber();
  VectorBool status(nech);
  for (int iech = 0; iech < nech; iech++)
    status[iech] = isActive(iech);
  return status;
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

/**
 * Sampling an input Db to create the output Db by selecting a subset of samples
 *
 * @param dbin       Pointer to the input Db
 * @param proportion Proportion of samples to be retained
 * @param number     Number of samples to be retained
 * @param names      Vector of Names to be copied (empty: all names)
 * @param seed       Seed used for the random number generator
 * @param verbose    Verbose flag
 *
 * @remark A possible selection in 'dbin' will not be taken into account
 * @remark You can use either 'proportion' or 'number'
 */
int Db::resetSamplingDb(const Db* dbin,
                        double proportion,
                        int number,
                        const VectorString& names,
                        int seed,
                        bool verbose)
{
  if (proportion <= 0. && number <= 0)
  {
    messerr("You must specify either 'proportion' or 'number'");
    return 1;
  }
  _clear();

  // Creating the vector of selected samples

  int nfrom = dbin->getSampleNumber();
  VectorInt ranks = ut_vector_sample(nfrom, proportion, number, seed);
  _nech = static_cast<int> (ranks.size());
  if (verbose)
    message("From %d samples, the extraction concerns %d samples\n", nfrom,_nech);

  // Create the new data base

  VectorString namloc = names;
  if (namloc.empty())
    namloc = dbin->getAllNames();
  _ncol = static_cast<int> (namloc.size());
  resetDims(_ncol, _nech);

  // Create Variables and Locators

  ELoc locatorType;
  int locatorIndex;
  for (int icol = 0; icol < _ncol; icol++)
  {
    setNameByUID(icol, namloc[icol]);
    if (dbin->getLocator(namloc[icol],&locatorType,&locatorIndex))
      setLocator(namloc[icol],locatorType,locatorIndex);
  }

  // Load samples

  VectorDouble values(_nech);
  for (int icol = 0; icol < _ncol; icol++)
  {
    int jcol = dbin->getColIdx(namloc[icol]);
    for (int iech = 0; iech < _nech; iech++)
      values[iech] = dbin->getValueByColIdx(ranks[iech],jcol);
    setColumnByColIdx(values, icol);
  }

  return 0;
}

/*****************************************************************************/
/*!
 **  Create a new Data Base with points generated at random
 **
 ** \return  Pointer for the new Db structure
 **
 ** \param[in]  nech        Expected number of samples
 ** \param[in]  dbgrid      Descriptor of the Db grid parameters
 ** \param[in]  seed        Seed for the random number generator
 ** \param[in]  flag_exact  True if the number of samples must not be drawn
 ** \param[in]  flag_repulsion True if repulsion is processed
 ** \param[in]  range       Repulsion range
 ** \param[in]  beta        Bending coefficient
 ** \param[in]  flag_add_rank 1 if the Rank must be generated in the output Db
 **
 *****************************************************************************/
Db* Db::createFromDbGrid(int nech,
                         DbGrid* dbgrid,
                         int seed,
                         bool flag_exact,
                         bool flag_repulsion,
                         double range,
                         double beta,
                         int flag_add_rank)
{
  Db* db = db_point_init(nech, VectorDouble(), VectorDouble(), dbgrid,
                         flag_exact, flag_repulsion, range, beta, seed,
                         flag_add_rank);
  return db;
}

/**
 * Combine 'sel' input argument with an already existing selection (if any)
 * @param sel Input selection (only 0 and 1)
 * @param combine Type of combination: "set", "not", "or", "and", "xor"
 * @remark Argument 'sel' may be modified by this procedure
 */
void Db::combineSelection(VectorDouble& sel, const String& combine) const
{
  int nech = (int) sel.size();
  if (nech <= 0) return;

  if (combine == "set")
    return;

  else if (combine == "not")
  {
    for (int iech = 0; iech < nech; iech++)
      sel[iech] = 1. - sel[iech];
    return;
  }

  else
  {
    // Read an already existing selection
    VectorDouble oldsel = getColumnByLocator(ELoc::SEL, 0);
    if (oldsel.empty()) return;

    if (combine == "or")
    {
      for (int iech = 0; iech < nech; iech++)
        sel[iech] = sel[iech] || oldsel[iech];
      return;
    }
    else if (combine == "and")
    {
      for (int iech = 0; iech < nech; iech++)
        sel[iech] = sel[iech] && oldsel[iech];
      return;
    }
    else if (combine == "xor")
    {
      for (int iech = 0; iech < nech; iech++)
        sel[iech] = sel[iech] != oldsel[iech];
      return;
    }
  }

  // The 'combine' argument is not valid

  messerr("Error in 'combine' argument. It should be one of the following ones:");
  messerr("('sel' is the current selection and 'oldsel' the already existing one)");
  messerr("'set': Do not combine with previous selection");
  messerr("'not': sel = 1 - sel");
  messerr("'or' : sel = sel || oldsel");
  messerr("'and': sel = sel && oldsel");
  messerr("'xor': sel = sel != oldsel");
}

Db* Db::create()
{
  return new Db();
}

Db* Db::createFromSamples(int nech,
                          const ELoadBy& order,
                          const VectorDouble& tab,
                          const VectorString& names,
                          const VectorString& locatorNames,
                          int flag_add_rank)
{
  Db* db = new Db;
  if (db->resetFromSamples(nech, order, tab, names, locatorNames,
                           flag_add_rank))
  {
    messerr("Error when creating Db from Samples");
    delete db;
    return nullptr;
  }
  return db;
}
Db* Db::createFromCSV(const String& filename,
                      const CSVformat& csv,
                      bool verbose,
                      int ncol_max,
                      int nrow_max,
                      int flag_add_rank)
{
  Db* db = new Db;
  if (db->resetFromCSV(filename, verbose, csv, ncol_max, nrow_max,
                       flag_add_rank))
  {
    messerr("Error when creating Db from Grid");
    delete db;
    return nullptr;
  }
  return db;
}

Db* Db::createFromBox(int nech,
                      const VectorDouble& coormin,
                      const VectorDouble& coormax,
                      int seed,
                      bool flag_exact,
                      bool flag_repulsion,
                      double range,
                      double beta,
                      int flag_add_rank)
{
  Db* db = db_point_init(nech, coormin, coormax, nullptr,
                         flag_exact, flag_repulsion, range, beta, seed,
                         flag_add_rank);
  return db;
}

Db* Db::createFromOnePoint(const VectorDouble& tab, int flag_add_rank)
{
  Db* db = new Db;
  if (db->resetFromOnePoint(tab, flag_add_rank))
  {
    messerr("Error when creating Db from One Point");
    delete db;
    return nullptr;
  }
  return db;
}
Db* Db::createSamplingDb(const Db* dbin,
                         double proportion,
                         int number,
                         const VectorString& names,
                         int seed,
                         bool verbose)
{
  Db* db = new Db;
  if (db->resetSamplingDb(dbin, proportion, number, names, seed, verbose))
  {
    messerr("Error when clearing Db by Sampling another Db");
    delete db;
    return nullptr;
  }
  return db;
}

bool Db::areSame(const String& name1,
                 const String& name2,
                 double eps,
                 bool useSel,
                 bool verbose)
{
  int ndiff = 0;
  VectorDouble tab1 = getColumn(name1, useSel);
  VectorDouble tab2 = getColumn(name2, useSel);
  if (tab1.empty() || tab2.empty()) return true;

  int nech = (int) tab1.size();
  for (int iech = 0; iech < nech; iech++)
  {
    int ntest = 0;
    if (FFFF(tab1[iech])) ntest++;
    if (FFFF(tab2[iech])) ntest++;
    if (ntest == 1) return false;
    if (ntest == 2) continue;
    double dist = tab1[iech] - tab2[iech];
    if (ABS(dist) > eps)
    {
      ndiff++;
      if (verbose)
        message("Sample #%d: V1=%lf V2=%lf\n",iech+1,tab1[iech],tab2[iech]);
    }
  }

  if (ndiff > 0)
    message("Differences between %s and %s (eps = %lf) = %d / %d\n",
            name1.c_str(),name2.c_str(),eps,ndiff,nech);
  else
    message("Variables %s and %s are similar (eps=%lf)\n",
            name1.c_str(),name2.c_str(),eps);
  return (ndiff > 0);
}

/**
 * Find the occurrence of a given range of values for a given variable
 *
 * @param name     Name of the Target variable
 * @param interval Interval definition
 * @param belowRow If specified, search must be performed below this row
 * @param aboveRow If specified, search must be performed above this row
 * @return
 */
VectorInt Db::filter(const String& name,
                     const Interval& interval,
                     int belowRow,
                     int aboveRow) const
{
  VectorInt rows;

  VectorDouble tab = getColumn(name, false);
  if (tab.empty()) return rows;

  int rankFrom = 0;
  if (! IFFFF(belowRow)) rankFrom = belowRow;
  int rankTo = getSampleNumber() - 1;
  if (! IFFFF(aboveRow)) rankTo = aboveRow;

  for (int irow = rankFrom; irow <= rankTo; irow++)
  {
    if (interval.isInside(tab[irow])) rows.push_back(irow);
  }
  return rows;
}
