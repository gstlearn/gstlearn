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
#include "Basic/VectorNumT.hpp"
#include "Basic/Law.hpp"
#include "Basic/AException.hpp"
#include "Basic/GlobalEnvironment.hpp"
#include "Basic/VectorHelper.hpp"
#include "Stats/Classical.hpp"
#include "Matrix/Table.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Space/SpacePoint.hpp"
#include "Space/SpaceTarget.hpp"

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
                         bool flagAddSampleRank)
{
  _clear();
  int ncol = (tab.empty()) ? 0 : (int) (tab.size() / nech);
  _ncol = (flagAddSampleRank) ? ncol + 1 : ncol;
  _nech = nech;
  resetDims(_ncol, _nech);

  // Load data (if defined)

  if (flagAddSampleRank) _createRank(0);
  _loadData(tab, names, locatorNames, order, (int) flagAddSampleRank);

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
 * @param flagAddSampleRank true if the sample rank must be generated
 */
int Db::resetFromCSV(const String& filename,
                     bool verbose,
                     const CSVformat& csvfmt,
                     int ncol_max,
                     int nrow_max,
                     bool flagAddSampleRank)
{
  _clear();
  VectorString names;
  VectorDouble tab;
  int ncol, nrow;

  /* Reading the CSV file */

  if (csv_table_read(filename, csvfmt, (int) verbose, ncol_max, nrow_max,
                     &ncol, &nrow, names, tab) != 0)
  {
    messerr("Problem when reading CSV file");
    return 1;
  }

  ncol = (tab.empty()) ? 0 : (int) (tab.size() / nrow);
  _ncol = (flagAddSampleRank) ? ncol + 1 : ncol;
  _nech = nrow;
  resetDims(_ncol, _nech);

  // Load data (if defined)

  if (flagAddSampleRank) _createRank(0);
  _loadData(tab, names, VectorString(), ELoadBy::SAMPLE, (int) flagAddSampleRank);

  // Set the names
  _defineDefaultNames((int) flagAddSampleRank, names);

  // Locators: Try to guess them from the Names
  _defineDefaultLocatorsByNames((int) flagAddSampleRank, names);

  return 0;
}

/**
 * Create a Db generating samples randomly
 *
 * @param nech    Number of samples to be generated
 * @param coormin Vector giving the smallest values of the coordinates
 * @param coormax Vector giving the largest values for the coordinates
 * @param ndim    Space dimension (used if 'coormin' and 'coormax' are empty)
 * @param extend  Extension of the bounding box (if positive)
 * @param seed    Seed for the random number generator
 * @param flagAddSampleRank true if the Sample ranks must be generated
 */
int Db::resetFromBox(int nech,
                     const VectorDouble& coormin,
                     const VectorDouble& coormax,
                     int ndim,
                     double extend,
                     int seed,
                     bool flagAddSampleRank)
{
  _clear();
  if (! coormin.empty()) ndim = (int) coormin.size();
  if (! coormax.empty()) ndim = MIN(ndim, (int) coormax.size());
  _ncol = (flagAddSampleRank) ? ndim + 1 : ndim;
  _nech = nech;
  resetDims(_ncol, _nech);

  // Generate the sample number
  if (flagAddSampleRank) _createRank(0);

  // Generate the coordinates
  law_set_random_seed(seed);
  VectorDouble tab(ndim * nech);
  int ecr = 0;
  for (int idim = 0; idim < ndim; idim++)
  {
    double mini = (coormin.empty()) ? 0. : coormin[idim];
    if (extend > 0.) mini -= extend;
    double maxi = (coormax.empty()) ? 1. : coormax[idim];
    if (extend > 0.) maxi += extend;
    message("idim=%d coormin=%lf mini=%lf coormax=%lf maxi=%lf\n",
            idim, coormin[idim], mini, coormax[idim], maxi);
    for (int iech = 0; iech < nech; iech++)
      tab[ecr++] = law_uniform(mini,maxi);
  }

  // Load the coordinates
  VectorString names = generateMultipleNames("x", ndim);
  _loadData(tab, names, VectorString(), ELoadBy::COLUMN, (int) flagAddSampleRank);

  int jcol = 0;
  if (flagAddSampleRank) jcol++;
  setLocatorsByUID(ndim, jcol, ELoc::X);

  return 0;
}

/**
 * Create a Db from a single sample whose coordinates are provided in 'tab'
 * @param tab Array containing the coordinates of the single sample
 * @param flagAddSampleRank true if the Sample ranks must be generated
 */
int Db::resetFromOnePoint(const VectorDouble& tab, bool flagAddSampleRank)
{
  _clear();

  int ndim = static_cast<int> (tab.size());
  _ncol = (flagAddSampleRank) ? ndim + 1 : ndim;
  _nech = 1;
  resetDims(_ncol, _nech);

  // Generate the sample number
  if (flagAddSampleRank) _createRank(0);

  // Load the coordinates
  VectorString names = generateMultipleNames("x", ndim);
  VectorDouble tabloc = tab;
  if (tabloc.empty()) tabloc.resize(ndim,0.);
  _loadData(tabloc, names, VectorString(), ELoadBy::SAMPLE, (int) flagAddSampleRank);

  int jcol = 0;
  if (flagAddSampleRank) jcol++;
  setLocatorsByUID(ndim, jcol, ELoc::X);

  return 0;
}

/**
 * Check if the argument 'idim' is a valid Space rank (0-based)
 */
bool Db::isDimensionIndexValid(int idim) const
{
  return checkArg("Space Dimension", idim, getNDim());
}

/**
 * Check if the argument 'iuid' is a valid user-designated rank
 */
bool Db::isUIDValid(int iuid) const
{
  return checkArg("UID Index", iuid, getUIDMaxNumber());
}

/**
 * Check if the argument 'icol' is a valid Column rank (0-based)
 */
bool Db::isColIdxValid(int icol) const
{
  return checkArg("Column Index", icol, _ncol);
}

/**
 * Check if the argument 'iech' is a valid Sample rank (0-based)
 */
bool Db::isSampleIndexValid(int iech) const
{
  return checkArg("Sample Index", iech, _nech);
}

/**
 * Check if the argument 'iechs' are valid Sample ranks (0-based)
 */
bool Db::isSampleIndicesValid(const VectorInt& iechs, bool useSel) const
{
  for (int i = 0; i < (int)iechs.size(); i++)
  {
    int iech = iechs[i];
    if (!checkArg("Sample Index", iech, getSampleNumber(useSel))) return false;
  }
  return true;
}

/**
 * Check if the arguments 'locatorType' and 'locatorIndex' are valid
 */
bool Db::isLocatorIndexValid(const ELoc& locatorType, int locatorIndex) const
{
  const PtrGeos& p = _p[locatorType.getValue()];
  return p.isLocatorIndexValid(locatorIndex);
}

int Db::getColIdxByUID(int iuid) const
{
  if (!isUIDValid(iuid)) return -1;
  int icol = _uidcol[iuid];
  return icol;
}

VectorInt Db::getColIdxsByUID(const VectorInt& iuids) const
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
  const PtrGeos& p = _p[locatorType.getValue()];
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
  const PtrGeos& p = _p[locatorType.getValue()];
  int number = p.getLocatorNumber();
  if (number <= 0 || locatorIndex >= number)
    return -1;
  int icol = getColIdxByUID(p.getLocatorByIndex(locatorIndex));
  return (icol);
}

int Db::getLocatorNumber(const ELoc& locatorType) const
{
  const PtrGeos& p = _p[locatorType.getValue()];
  return p.getLocatorNumber();
}

int Db::_findUIDInLocator(const ELoc& locatorType, int iuid) const
{
  const PtrGeos& p = _p[locatorType.getValue()];
  if (!isUIDValid(iuid)) return -1;
  for (int locatorIndex = 0; locatorIndex < p.getLocatorNumber(); locatorIndex++)
    if (p.getLocatorByIndex(locatorIndex) == iuid) return (locatorIndex);
  return -1;
}

int Db::_findColumnInLocator(const ELoc& locatorType, int icol) const
{
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
  int number = getNEloc();
  for (int iloc = 0; iloc < number; iloc++)
  {
    const PtrGeos& p = _p[iloc];
    for (int i = 0; i < p.getLocatorNumber(); i++)
    {
      int jcol = getColIdxByUID(p.getLocatorByIndex(i));
      if (icol == jcol)
      {
        *ret_locatorType = ELoc::fromValue(iloc);
        *ret_locatorIndex = i;
        return true;
      }
    }
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
 * Set the values of a series of samples
 * @param iechs List of sample indices
 * @param iuid Index of the UID
 * @param values List of values to be written
 * @remarks: for efficiency purpose, no check is performed on the sample ranks
 */
void Db::setArrayVec(const VectorInt& iechs, int iuid, const VectorDouble& values)
{
  int icol = getColIdxByUID(iuid);
  if (!isColIdxValid(icol)) return;
  for (int i = 0, n = (int) iechs.size(); i < n; i++)
    _array[_getAddress(iechs[i], icol)] = values[i];
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

/**
 * Get the values of a series of samples
 * @param iechs List of sample indices
 * @param iuid Index of the UID
 * @param values List of values to be written
 * @remarks: for efficiency purpose, no check is performed on the sample ranks
 */
void Db::getArrayVec(const VectorInt& iechs, int iuid, VectorDouble& values) const
{
  int icol = getColIdxByUID(iuid);
  if (!isColIdxValid(icol)) return;
  for (int i = 0, n = (int) iechs.size(); i < n; i++)
  {
    values[i] = _array[_getAddress(iechs[i], icol)];
  }
}

VectorDouble Db::getArrayByUID(int iuid, bool useSel) const
{
  int nech = getSampleNumber();
  VectorDouble sel, tab;
  if (!isUIDValid(iuid)) return tab;

  tab.resize(nech);
  if (useSel) sel = getSelections();

  int ecr = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    if (useSel && !sel.empty() && isZero(sel[iech])) continue;
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

void Db::updArray(int iech, int iuid, const EOperator& oper, double value)
{
  if (!isSampleIndexValid(iech)) return;

  int icol = getColIdxByUID(iuid);
  if (!isColIdxValid(icol)) return;

  int internalAddress = _getAddress(iech, icol);
  double oldval = _array[internalAddress];
  double newval = modifyOperator(oper, oldval, value);
  _array[internalAddress] = newval;
}

void Db::updArrayVec(const VectorInt& iechs, int iuid, const EOperator& oper, VectorDouble& values)
{
  int icol = getColIdxByUID(iuid);
  if (!isColIdxValid(icol)) return;

  int iad;
  double oldval;
  double newval;
  for (int i = 0, n = (int) iechs.size(); i < n; i++)
  {
    iad = _getAddress(iechs[i], icol);
    oldval = _array[iad];
    newval = modifyOperator(oper, oldval, values[i]);
    _array[iad] = newval;
  }
}

VectorDouble Db::getSampleCoordinates(int iech) const
{
  VectorDouble coor(getNDim());
  getSampleCoordinatesInPlace(iech, coor);
  return coor;
}

void Db::getSampleAsSPInPlace(int iech, SpacePoint& P) const
{
  getCoordinatesPerSampleInPlace(iech, P.getCoordRef());
}

VectorVectorDouble Db::getIncrements(const VectorInt& iechs, const VectorInt& jechs) const
{
  VectorVectorDouble tab;
  int ndim = getNDim();
  SpacePoint P1(ndim);
  SpacePoint P2(ndim);

  int number = (int) iechs.size();
  if ((int) jechs.size() != number)
  {
    messerr("Arguments 'iechs'(%d) and 'jechs'(%d) should share the same dimension",
            (int) iechs.size(), (int) jechs.size());
    return tab;
  }

  // Dimension the output vector
  tab.resize(ndim);
  for (int idim = 0; idim < ndim; idim++) tab[idim].resize(number);

  for (int ip = 0; ip < number; ip++)
  {
    getSampleAsSPInPlace(iechs[ip], P1);
    getSampleAsSPInPlace(jechs[ip], P2);
    VectorDouble vect = P2.getIncrement(P1);

    for (int idim = 0; idim < ndim; idim++)
      tab[idim][ip] = vect[idim];
  }
  return tab;
}

/**
 * Load a Space Target with all possible contents gathered from Db
 * @param iech Rank of the target sample
 * @param P    Space Target (used to store information)
 */
void Db::getSampleAsSTInPlace(int iech, SpaceTarget& P) const
{
  // Load the coordinates
  getSampleAsSPInPlace(iech, P);

  // Load the code (optional)
  if (P.checkCode())
  {
    if (hasLocVariable(ELoc::C)) P.setCode(getLocVariable(ELoc::C, iech, 0));
  }

  // Load the Date (optional)
  if (P.checkDate())
  {
    if (hasLocVariable(ELoc::DATE))
      P.setCode(getLocVariable(ELoc::DATE, iech, 0));
  }
}

std::vector<SpacePoint> Db::getSamplesAsSP(bool useSel) const
{
  std::vector<SpacePoint> pvec;
  VectorDouble coord(getNDim());
  SpacePoint p;
  for (int iech = 0, nech = getSampleNumber(); iech < nech; iech++)
  {
    if (isActive(iech))
    {
      getSampleAsSPInPlace(iech, p);
      pvec.push_back(p);
    }
    else
    {
      if (useSel) continue;
      p.setFFFF();
      pvec.push_back(p);
    }
  }
  return pvec;
}

VectorDouble Db::getSampleLocators(const ELoc& locatorType, int iech) const
{
  VectorDouble vec;
  int number = getLocatorNumber(locatorType);
  if (number <= 0) return vec;
  vec.resize(number);
  for (int i = 0; i < number; i++)
    vec[i] = getFromLocator(locatorType, iech, i);
  return vec;
}

void Db::getSampleCoordinatesInPlace(int iech, VectorDouble& coor) const
{
  for (int idim = 0, ndim=getNDim(); idim < ndim; idim++)
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

void Db::getCoordinatesPerSampleInPlace(int iech, VectorDouble& coor, bool flag_rotate) const
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
  VectorDouble dd(ndim);
  if (getDistanceVecInPlace(iech, jech, dd) != 0) return TEST;
  double dist = 0.;
  for (int idim = 0; idim < ndim; idim++)
  {
    double delta = dd[idim];
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
// TODO to be corrected to use SpaceDistance
int Db::getDistanceVecInPlace(int iech, int jech, VectorDouble& dd, const Db* db2) const
{
  int ndim = getNDim();
  VectorDouble v1(ndim);
  VectorDouble v2(ndim);

  getCoordinatesPerSampleInPlace(iech, v1);
  if (db2 == nullptr)
    getCoordinatesPerSampleInPlace(jech, v2);
  else
    db2->getCoordinatesPerSampleInPlace(jech, v2);
  for (int idim = 0; idim < ndim; idim++)
    dd[idim] = v1[idim] - v2[idim];
  return 0;
}

/**
 * Constitute a Vector of Vector of coordinates for all (active) samples
 * - the first dimension is the space dimension
 * - the second dimension is the number of (active) samples
 * @param useSel
 * @return
 */
VectorVectorDouble Db::getAllCoordinates(bool useSel) const
{
  VectorVectorDouble result;
  for (int idim = 0, ndim = getNDim(); idim < ndim; idim++)
  {
    VectorDouble local = getCoordinates(idim, useSel);
    result.push_back(local);
  }
  return result;
}

/**
 * Constitute a Matrix of coordinates for all (active) samples
 * - one row per sample
 * - one column by Space Dimension
 * @return
 */
MatrixRectangular Db::getAllCoordinatesMat() const
{
  int nech = getSampleNumber(true);
  int ndim = getNDim();

  MatrixRectangular mat(nech, ndim);

  VectorInt ranks = getRanksActive();
  for (int jech = 0; jech < nech; jech++)
  {
    int iech = ranks[jech];
    VectorDouble coors = getSampleCoordinates(iech);
    mat.setRow(iech, coors);
  }
  return mat;
}

void Db::setCoordinate(int iech, int idim, double value)
{
  if (!isSampleIndexValid(iech)) return;
  int icol = getColIdxByLocator(ELoc::X, idim);
  if (!isColIdxValid(icol)) return;
  _array[_getAddress(iech, icol)] = value;
}

void Db::setCoordinates(int idim, const VectorDouble& coor, bool useSel)
{
  int icol = getColIdxByLocator(ELoc::X, idim);
  if (!isColIdxValid(icol)) return;
  setColumnByColIdx(coor, icol, useSel);
}

void Db::setSampleCoordinates(int iech, const VectorDouble& coor)
{
  int ndim = getNDim();
  int size = (int)coor.size();
  if (ndim != size)
  {
    messerr("Argument 'coor' (%d) should have dimension ndim (%d)", size, ndim);
    messerr("Nothing is done");
    return;
  }
  for (int idim = 0, ndim = getNDim(); idim < ndim; idim++)
    setCoordinate(iech, idim, coor[idim]);
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

bool Db::hasLocator(const ELoc& locatorType) const
{
  const PtrGeos& p = _p[locatorType.getValue()];
  return p.hasLocator();
}

int Db::getFromLocatorNumber(const ELoc& locatorType) const
{
  const PtrGeos& p = _p[locatorType.getValue()];
  return p.getLocatorNumber();
}

int Db::getNEloc()
{
  int number = 0;
  auto it = ELoc::getIterator();
  while (it.hasNext())
  {
    if (*it != ELoc::UNKNOWN) number++;
    it.toNext();
  }
  return number;
}

void Db::_clear(void)
{
  _p.clear();
  int number = getNEloc();
  _p.resize(number);
  for (int iloc = 0; iloc < number; iloc++)
    _p[iloc] = PtrGeos();
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
  int number = getNEloc();
  for (int iloc = 0; iloc < number; iloc++)
  {
    const PtrGeos& p = _p[iloc];
    if (p.getLocatorNumber() > 0)
    {
      sstr << p.dumpLocator(rank, ELoc::fromValue(iloc));
      sstr << "- Columns    = ";
      for (int locatorIndex = 0; locatorIndex < p.getLocatorNumber(); locatorIndex++)
        sstr << getColIdxByUID(p.getLocatorByIndex(locatorIndex)) << " ";
      sstr << std::endl;
      rank++;
    }
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
  PtrGeos& p = _p[locatorType.getValue()];
  p.clear();
}

/**
 * Setting the locator for a set of variables designated by their names
 * @param names        Vector of variable names
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
  VectorInt iuids = _ids(names, false);
  if (iuids.empty()) return;

  if (cleanSameLocator) clearLocators(locatorType);

  for (unsigned int i = 0; i < iuids.size(); i++)
    setLocatorByUID(iuids[i], locatorType, locatorIndex + i);
}

/**
 * Define the Locator(s) for the given variable(s)
 * @param name Variable name
 * @param locatorType Locator Type
 * @param locatorIndex Locator Index (for the first variable) (starting from 0)
 * @param cleanSameLocator When TRUE, clean variables with same locator beforehand
 */
void Db::setLocator(const String &name,
                    const ELoc& locatorType,
                    int locatorIndex,
                    bool cleanSameLocator)
{
  VectorInt iuids = _ids(name, false);
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
  if (locatorIndex < 0) return;

  // Optional clean

  if (cleanSameLocator) clearLocators(locatorType);

  /* Cancel any locator referring to this column */

  int number = getNEloc();
  for (int iloc = 0; iloc < number; iloc++)
  {
    PtrGeos& p = _p[iloc];
    int found = p.findUIDInLocator(iuid);
    if (found >= 0)
        p.erase(found);
  }

  // Check if this locator already exists for the current pointer type
  // Warning: the following code does not forbid declaring locatorIndex
  // in incorrect order. This must be kept as long as the Demonstration files
  // use the db.locate() of unsorted ranks

  if (locatorType != ELoc::UNKNOWN)
  {
    PtrGeos& p = _p[locatorType.getValue()];
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
  if (cleanSameLocator) clearLocators(locatorType);

  for (int i = 0; i < number; i++)
    setLocatorByUID(iuid+i, locatorType, locatorIndex + i);
}

void Db::setLocatorsByUID(const VectorInt& iuids,
                          const ELoc& locatorType,
                          int locatorIndex,
                          bool cleanSameLocator)
{
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
  _columnInit(nadd, ncol, true, valinit);

  // Set the locator (if defined)
  if (locatorType != ELoc::UNKNOWN)
    setLocatorsByUID(nadd, nmax, locatorType, locatorIndex);

  _ncol += nadd;

  return (nmax);
}

/**
 * Create a set of new variables in an already existing Db and initialize
 * their contents as a ranom value (from Normal distribution)
 * @param nadd     Number of variables to be added
 * @param radix    Generic radix given to the newly created variables
 * @param locatorType Generic locator assigned to new variables
 * @param locatorIndex   Locator index (starting from 0)
 * @param seed     Seed value
 * @param nechInit Number of samples (used only if the Db is initially empty)
 * @return Rank of the first UID
 */
int Db::addColumnsRandom(int nadd,
                         const String &radix,
                         const ELoc &locatorType,
                         int locatorIndex,
                         int seed,
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

  // Initialize the variables with a random value
  law_set_random_seed(seed);
  _columnInit(nadd, ncol, false);

  // Set the locator (if defined)
  if (locatorType != ELoc::UNKNOWN)
    setLocatorsByUID(nadd, nmax, locatorType, locatorIndex);

  _ncol += nadd;

  return (nmax);
}

void Db::addColumnsByVVD(const VectorVectorDouble& tab,
                         const String &radix,
                         const ELoc &locatorType,
                         int locatorIndex,
                         bool useSel)
{
  VectorDouble tabv;
  int nvar = (int) tab.size();
  for(const auto &e : tab)
    for(const auto &f : e)
      tabv.push_back(f);
  addColumns(tabv,radix,locatorType,locatorIndex,useSel,TEST,nvar);
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
  // If the input array 'tab' is empty, nothing is done
  if (tab.empty()) return 0;

  // Particular case where the Db is empty.
  // Set its dimension to the number of samples of the input array 'tab'
  if (_nech <= 0) _nech = static_cast<int> (tab.size()) / nvar;

  // Check dimensions
  int nech = getSampleNumber(useSel);
  nvar = (int) tab.size() / nech;
  if ((int) tab.size() != nvar * nech)
  {
    messerr("Db::addColumns : Incompatibility between 'tab'(%d) and 'nvar'(%d) * 'nech'(%d)",
            tab.size(), nvar, nech);
    return 1;
  }

  // Adding the new Columns
  int iuid = addColumnsByConstant(nvar, valinit, radix, locatorType, locatorIndex);
  if (iuid < 0) return 1;

  const double* local = tab.data();
  for (int ivar = 0; ivar < nvar; ivar++)
    setColumnByUIDOldStyle(&local[ivar * nech], iuid + ivar, useSel);

  return iuid;
}

void Db::setColumnByColIdxOldStyle(const double* tab, int icol, bool useSel)
{
  if (!isColIdxValid(icol)) return;
  VectorDouble sel;

  if (useSel) sel = getSelections();

  int lec = 0;
  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    bool defined = true;
    if (useSel && !sel.empty()) defined = (isOne(sel[iech]));

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

/**
 * Update the contents of an already existing variable in a Db
 * @param tab    Vector containing the values to be written
 * @param iuid   UID of the already existing variable to be written
 * @param useSel When TRUE, take the Selection into account (seed remarks)
 *
 * @remarks When useSel=TRUE, the input vector should be dimensioned to
 * @remarks the number of active samples. Only the active samples of the Db
 * @remarks are updated using the contents of the input 'tab' vector.
 */
void Db::setColumnByUIDOldStyle(const double* tab, int iuid, bool useSel)
{
  if (!isUIDValid(iuid)) return;
  VectorDouble sel;

  if (useSel) sel = getSelections();

  int lec = 0;
  bool defined = true;
  for (int iech = 0, nech = getSampleNumber(); iech < nech; iech++)
  {
    defined = true;
    if (!sel.empty()) defined = (isOne(sel[iech]));

    if (defined)
      setArray(iech, iuid, tab[lec++]);
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

/**
 * Delete one Column specified by its name
 *
 */
void Db::deleteColumn(const String& name)
{
  VectorInt iuids = _ids(name, false);
  if (iuids.empty()) return;

  for (unsigned int i = 0; i < iuids.size(); i++)
    deleteColumnByUID(iuids[i]);
}

/**
 * Delete a set of variables specified by their names
 *
 */
void Db::deleteColumns(const VectorString& names)
{
  VectorInt iuids = _ids(names, false);
  if (iuids.empty()) return;

  for (unsigned int i = 0; i < iuids.size(); i++)
    deleteColumnByUID(iuids[i]);
}

/**
 * Delete a set of variables specified by their column ranks (0 based)
 *
 */
void Db::deleteColumnsByColIdx(const VectorInt& icols)
{
  if (icols.empty()) return;

  // Reverse order of the columns in order to start by the furthest one.
  VectorInt v = VH::sort(icols, false);

  for (unsigned int i = 0; i < v.size(); i++)
    deleteColumnByColIdx(v[i]);
}

/**
 * Delete a set of variables specified by their user-identification ranks (0 based)
 *
 */
void Db::deleteColumnsByUID(const VectorInt& iuids)
{
  if (iuids.empty()) return;

  for (unsigned int i = 0; i < iuids.size(); i++)
    deleteColumnByUID(iuids[i]);
}

void Db::deleteColumnsByUIDRange(int i_del, int n_del)
{
  if (i_del <= 0) return;
  for (int i = n_del - 1; i >= 0; i--)
    deleteColumnByUID(i_del + i);
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
int Db::addSelection(const VectorDouble &tab,
                     const String &name,
                     const String &combine)
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
      sel[iech] = (! isZero(tab[iech])) ? 1. : 0.;
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
 * Create a Selection based on the Convex Hull of the active samples of 'Db'
 * @param db       Data Base providing the (active) samples
 * @param dilate   The convex hull can be dilated: this gives the radius
 * @param verbose  Verbose option
 * @param namconv  Naming Convention
 * @return
 */
int Db::addSelectionFromDbByConvexHull(Db *db,
                                       double dilate,
                                       bool verbose,
                                       const NamingConvention &namconv)
{
  if (db == nullptr)
  {
    messerr("You must define a valid Db");
    return 1;
  }

  return db_selhull(db, this, dilate, verbose, namconv);
}

/**
 * Create a Selection based on a proportion of active samples
 * @param prop   Proportion of active samples (between 0 and 1)
 * @param seed   Seed for the random number generator
 * @param name   Name of the newly created selection
 * @param combine How to combine with an already existing selection (see combineSelection() for details)
 * @return
 */
int Db::addSelectionRandom(double prop,
                           int seed,
                           const String &name,
                           const String &combine)
{
  VectorInt ranks = VH::sampleRanks(getSampleNumber(false),prop,-1,seed,1);
  return addSelectionByRanks(ranks, name, combine);
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

int Db::deleteSamples(const VectorInt& e_dels)
{
  if (e_dels.empty()) return 0;

  // Reverse order of the samples in order to start by the furthest one.
  VectorInt v = VH::sort(e_dels, false);

  for (unsigned int i = 0; i < v.size(); i++)
    if (deleteSample(v[i]) != 0) return 1;

  return 0;
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
      int jech = (iech < e_del) ? iech : iech - 1;
      int iad1 = jech + nnew * icol;
      new_array[iad1] = _array[_getAddress(iech, icol)];
    }

  /* Core deallocation */

  _array = new_array;
  _nech = nnew;
  return 0;
}

/**
 * Delete a variablesspecified by its column number (0 based)
 *
 */
void Db::deleteColumnByColIdx(int icol_del)
{
  if (! isColIdxValid(icol_del)) return;
  VectorInt iuids = _ids(_colNames[icol_del],true);
  if (iuids.empty()) return;
  deleteColumnByUID(iuids[0]);
}

/**
 * Delete a variable specified by its user-identification rank
 *
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

  int number = getNEloc();
  for (int iloc = 0; iloc < number; iloc++)
  {
    PtrGeos& p = _p[iloc];
    int found = p.findUIDInLocator(iuid_del);
    if (found >= 0) p.erase(found);
  }

  /* Resize the variables names */

  _colNames.erase(_colNames.begin() + c_del);

  /* Set the error return code */

  _ncol = nnew;
}

/**
 * Delete a set of variables specified by their locator type
 *
 */
void Db::deleteColumnsByLocator(const ELoc& locatorType)
{
  const PtrGeos& p = _p[locatorType.getValue()];
  int nitem = p.getLocatorNumber();
  // Loop is performed downwards as PtrGeos is modified by called routine
  for (int locatorIndex = nitem - 1; locatorIndex >= 0; locatorIndex--)
  {
    deleteColumnByUID(p.getLocatorByIndex(locatorIndex));
  }
}

/**
 * Returns the extreme coordinates for the target space dimension
 *
 */
VectorDouble Db::getExtrema(int idim, bool useSel) const
{
  VectorDouble ext;
  if (!isDimensionIndexValid(idim)) return ext;
  VectorDouble coor = getCoordinates(idim, useSel);
  ext.push_back(VH::minimum(coor));
  ext.push_back(VH::maximum(coor));
  return ext;
}

/**
 * Returns the extreme coordinates for all space dimensions
 *
 */
VectorVectorDouble Db::getExtremas(bool useSel) const
{
  VectorVectorDouble exts;
  for (int idim = 0; idim < getNDim(); idim++)
    exts.push_back(getExtrema(idim, useSel));
  return exts;
}

/**
 * Returns the minimum coordinates for all space dimensions
 *
 */
VectorDouble Db::getCoorMinimum(bool useSel) const
{
  VectorDouble ext;
  for (int idim = 0; idim < getNDim(); idim++)
  {
    VectorDouble coor = getCoordinates(idim, useSel);
    ext.push_back(VH::minimum(coor));
  }
  return ext;
}

/**
 * Returns the maximum coordinates for all space dimensions
 *
 */
VectorDouble Db::getCoorMaximum(bool useSel) const
{
  VectorDouble ext;
  for (int idim = 0; idim < getNDim(); idim++)
  {
    VectorDouble coor = getCoordinates(idim, useSel);
    ext.push_back(VH::maximum(coor));
  }
  return ext;
}

/**
 * Returns the coordinates of the center of the (active) samples
 *
 */
VectorDouble Db::getCenters(bool useSel) const
{
  int ndim = getNDim();
  VectorDouble center(ndim);
  for (int idim = 0; idim < ndim; idim++)
    center[idim] = getCenter(idim, useSel);
  return center;
}

/**
 * Returns the center of the (active) samples for the target space dimension
 *
 */
double Db::getCenter(int idim, bool useSel) const
{
  if (!isDimensionIndexValid(idim)) return TEST;
  VectorDouble coor = getCoordinates(idim, useSel);
  double mini = VH::minimum(coor);
  double maxi = VH::maximum(coor);
  return ((mini + maxi) / 2.);
}

/**
 * Returns the extension (distance between minimum and maximum) for the target space dimension
 *
 */
double Db::getExtension(int idim, bool useSel) const
{
  if (!isDimensionIndexValid(idim)) return 0.;
  VectorDouble coor = getCoordinates(idim, useSel);
  double mini = VH::minimum(coor);
  double maxi = VH::maximum(coor);
  return maxi - mini;
}

/**
 * Returns the diagonal of the rectangle containing all (active) samples and parallel to main axes
 *
 */
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

/**
 * Returns the extensions (distance between minimum and maximum) for all space dimensions
 *
 */
void Db::getExtensionInPlace(VectorDouble& mini,
                             VectorDouble& maxi,
                             bool flagPreserve,
                             bool useSel) const
{
  int ndim = getNDim();
  if (ndim != (int) mini.size()) mini.resize(ndim,TEST);
  if (ndim != (int)maxi.size()) maxi.resize(ndim, TEST);

  // If flagPreserve is false, the output arguments are reset beforehand
  if (!flagPreserve)
  {
    for (int idim = 0; idim < ndim; idim++)
    {
      mini[idim] = maxi[idim] = TEST;
    }
  }

  /* Loop on the space dimension */

  for (int idim = 0; idim < getNDim(); idim++)
  {
    VectorDouble coor = getCoordinates(idim, useSel);
    double vmin = VH::minimum(coor);
    double vmax = VH::maximum(coor);
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

/**
 * Identify the list of names. These names are searched in the following order:
 * - within the list of input variable names (possibly expanded)
 * - within the names of the locators
 * @param names Names to be be identified
 * @return List of variable names
 */
VectorString Db::identifyNames(const VectorString& names) const
{
  VectorString ret_names;
  VectorString namloc;
  ELoc locatorType;
  int locatorIndex, mult;

  // Constitute the list of the locator names
  VectorString locnames;
  for (int j = 0; j < getColumnNumber(); j++)
  {
    if (! getLocatorByColIdx(j, &locatorType, &locatorIndex)) continue;
    String local = getLocatorName(locatorType, locatorIndex);
    locnames.push_back(local);
  }

  for (int i = 0; i < (int) names.size(); i++)
  {
    // Look within the list of names
    namloc = getName(names[i]);
    if (! namloc.empty())
    {
      for (int j = 0; j < (int) namloc.size(); j++)
        ret_names.push_back(namloc[j]);
      continue;
    }

    // Look within the list of locators
    namloc = expandList(locnames, names[i]);
    if (! namloc.empty())
    {
      for (int j = 0; j < (int) namloc.size(); j++)
      {
        if (locatorIdentify(namloc[j], &locatorType, &locatorIndex, &mult) != 0) continue;
        String local = getNameByLocator(locatorType, locatorIndex);
        ret_names.push_back(local);
      }
      continue;
    }
  }
  return ret_names;
}

/**
 * Returns the minimum of the target variable
 */
double Db::getMinimum(const String& name, bool useSel) const
{
  VectorInt iuids = _ids(name, true);
  if (iuids.empty()) return TEST;
  VectorDouble tab = getColumnByUID(iuids[0], useSel);
  return VH::minimum(tab);
}

/**
 * Returns the maximum of the target variable
 */
double Db::getMaximum(const String& name, bool useSel) const
{
  VectorInt iuids = _ids(name, true);
  if (iuids.empty()) return TEST;
  VectorDouble tab = getColumnByUID(iuids[0], useSel);
  return VH::maximum(tab);
}

/**
 * Returns a vector containing the minimum and maximum of the target variable
 */
VectorDouble Db::getRange(const String& name, bool useSel) const
{
  VectorDouble range(2);
  range[0] = getMinimum(name, useSel);
  range[1] = getMaximum(name, useSel);
  return range;
}

/**
 * Returns the mean of the target variable
 */
double Db::getMean(const String& name, bool useSel) const
{
  VectorInt iuids = _ids(name, true);
  if (iuids.empty()) return TEST;
  VectorDouble tab = getColumnByUID(iuids[0], useSel);
  return VH::mean(tab);
}

/**
 * Returns the variance of the target variable
 */
double Db::getVariance(const String& name, bool useSel) const
{
  VectorInt iuids = _ids(name, true);
  if (iuids.empty()) return TEST;
  VectorDouble tab = getColumnByUID(iuids[0], useSel);
  return VH::variance(tab);
}

/**
 * Returns the standard deviation (square root of the variance) of the target variable
 */
double Db::getStdv(const String& name, bool useSel) const
{
  VectorInt iuids = _ids(name, true);
  if (iuids.empty()) return TEST;
  VectorDouble tab = getColumnByUID(iuids[0], useSel);
  return VH::stdv(tab);
}

/**
 * Returns the correlation coefficient between two target variables
 */
double Db::getCorrelation(const String& name1, const String& name2, bool useSel) const
{
  VectorInt iuids;
  iuids = _ids(name1, true);
  if (iuids.empty()) return TEST;
  VectorDouble tab1 = getColumnByUID(iuids[0], useSel);
  iuids = _ids(name2, true);
  if (iuids.empty()) return TEST;
  VectorDouble tab2 = getColumnByUID(iuids[0], useSel);
  return VH::correlation(tab1, tab2);
}

int Db::getNDim() const
{
  return _p[ELoc::X.getValue()].getLocatorNumber();
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


void Db::_columnInit(int ncol, int icol0, bool flagCst, double valinit)
{
  double value;
  if (flagCst)
    value = valinit;
  else
    value = law_gaussian();
  for (int jcol = 0; jcol < ncol; jcol++)
  {
    int icol = jcol + icol0;

    if (!GlobalEnvironment::getEnv()->isDomainReference() || !hasLocator(ELoc::DOM))
    {
      for (int iech = 0; iech < _nech; iech++)
        _array[_getAddress(iech, icol)] = value;
    }
    else
    {
      for (int iech = 0; iech < _nech; iech++)
      {
        value   = getFromLocator(ELoc::DOM, iech, 0);
        int iad = _getAddress(iech, icol);
        if (GlobalEnvironment::getEnv()->matchDomainReference(value))
          _array[iad] = value;
        else
          _array[iad] = TEST;
      }
    }
  }
}

void Db::switchLocator(const ELoc& locatorType_in, const ELoc& locatorType_out)
{
  PtrGeos& p_in  = _p[locatorType_in.getValue()];
  PtrGeos& p_out = _p[locatorType_out.getValue()];
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

VectorDouble Db::getValuesByNames(const VectorInt &iechs,
                                  const VectorString &names,
                                  bool bySample) const
{
  VectorInt icols = getColIdxs(names);
  return getValuesByColIdx(iechs, icols, bySample);
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

void Db::setValuesByNames(const VectorInt &iechs,
                          const VectorString &names,
                          const VectorDouble &values,
                          bool bySample)
{
  VectorInt icols = getColIdxs(names);
  setValuesByColIdx(iechs, icols, values, bySample);
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

/**
 * Returns the number of fields corresponding to the target locator present in the Db
 *
 * @return Number of fields
 */
int Db::getLocNumber(const ELoc& loctype) const
{
  if (loctype == ELoc::UNKNOWN) return 0;
  const PtrGeos& p = _p[loctype.getValue()];
  return p.getLocatorNumber();
}
int Db::getZNumber() const
{
  const PtrGeos& p = _p[ELoc::Z.getValue()];
  return p.getLocatorNumber();
}

/**
 * Check if there is at least one field corresponding to the target locator
 *
 * @return TRUE if at least one field corresponds to 'loctype' locator; FALSE otherwise
 */
bool Db::hasLocVariable(const ELoc& loctype) const
{
  if (loctype == ELoc::UNKNOWN) return false;
  return (int) hasLocator(loctype);
}
bool Db::hasZVariable() const
{
  return (int) hasLocator(ELoc::Z);
}

/**
 * Get the value of the field corresponding to the target locator (and its target item) at the target sample
 *
 * @return Returned value
 */
double Db::getLocVariable(const ELoc& loctype, int iech, int item) const
{
  if (!hasLocVariable(loctype)) return (TEST);
  return getFromLocator(loctype, iech, item);
}
double Db::getZVariable(int iech, int item) const
{
  return getFromLocator(ELoc::Z, iech, item);
}
VectorDouble Db::getLocVariables(const ELoc& loctype, int iech, int nitemax) const
{
  VectorDouble vec;
  int number = getFromLocatorNumber(loctype);
  if (number <= 0) return vec;
  int nitem = (nitemax > 0) ? MIN(nitemax, number) : number;

  vec.resize(nitem, TEST);
  for (int item = 0; item < nitem; item++)
    vec[item] = getLocVariable(loctype, iech, item);
  return vec;
}

/**
 *  Set the value of the field corresponding to the target locator (and its target item) at the target sample
 *
 */
void Db::setLocVariable(const ELoc& loctype, int iech, int item, double value)
{
  if (loctype == ELoc::UNKNOWN) return;
  setFromLocator(loctype, iech, item, value);
}
void Db::setZVariable(int iech, int item, double value)
{
  setFromLocator(ELoc::Z, iech, item, value);
}
void Db::setLocVariables(const ELoc& loctype,
                         int iech,
                         const VectorDouble& values)
{
  int number = getFromLocatorNumber(loctype);
  int size = (int) values.size();
  if (number != size)
  {
    messerr("Dimension of 'values' (%d) does not match number of elements in "
            "locator (%d)",
            size, number);
    messerr("Nothing is done");
    return;
  }
  for (int i = 0; i < number; i++) setFromLocator(loctype, iech, i, values[i]);
}

/**
 *  Update the value of the field corresponding to the target locator (and its target item) at the target sample
 *
 */
void Db::updLocVariable(const ELoc& loctype, int iech, int item, const EOperator& oper, double value)
{
  if (loctype == ELoc::UNKNOWN) return;
  if (!isSampleIndexValid(iech)) return;
  int icol = getColIdxByLocator(loctype, item);
  int internalAddress = _getAddress(iech, icol);

  double oldval = _array[internalAddress];
  double newval = modifyOperator(oper, oldval, value);
  _array[internalAddress] = newval;
}
void Db::updZVariable(int iech, int item, const EOperator& oper, double value)
{
  if (!isSampleIndexValid(iech)) return;
  int icol            = getColIdxByLocator(ELoc::Z, item);
  int internalAddress = _getAddress(iech, icol);

  double oldval           = _array[internalAddress];
  double newval           = modifyOperator(oper, oldval, value);
  _array[internalAddress] = newval;
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
    if (getLocNumber(ELoc::Z) != nvar)
    {
      messerr("This function requires %d variables (locator 'Z'). The 'Db' contains %d variables",
              nvar,getLocNumber(ELoc::Z));
      return false;
    }
  }
  else if (compare < 0)
  {
    if (! (getLocNumber(ELoc::Z) <= nvar))
    {
      messerr("This function requires nvar <= %d variables (locator 'Z'). The 'Db' contains %d variables",
              nvar,getLocNumber(ELoc::Z));
      return false;
    }
  }
  else
  {
    if (! (getLocNumber(ELoc::Z) > nvar))
    {
      messerr("This function requires nvar >= %d variables (locator 'Z'). The 'Db' contains %d variables",
              nvar,getLocNumber(ELoc::Z));
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
 *
 * @remark
 * The returned answer is false is there is no variable defined
 * or if the sample rank is not valid.
 * If 'nvar-max' is defined, the test is performed on the 'nvar_max'
 * first variables. Otherwise, it is performed on all ELOC.Z variables
 */
bool Db::isIsotopic(int iech, int nvar_max) const
{
  int nvar = getLocNumber(ELoc::Z);
  if (nvar_max > 0) nvar = MIN(nvar, nvar_max);
  if (nvar <= 0) return false;
  if (!isSampleIndexValid(iech)) return false;

  for (int ivar = 0; ivar < nvar; ivar++)
    if (FFFF(getZVariable(iech, ivar))) return false;
  return true;
}

/**
 * Check that all the active samples are isotopic
 */
bool Db::isAllIsotopic() const
{
  for (int iech = 0, nech = getSampleNumber(); iech < nech; iech++)
  {
    if (! isIsotopic(iech)) return false;
  }
  return true;
}

bool Db::isAllUndefined(int iech) const
{
  if (!isSampleIndexValid(iech)) return false;
  int nvar = getLocNumber(ELoc::Z);
  if (nvar <= 0) return false;

  for (int ivar = 0; ivar < nvar; ivar++)
    if (! FFFF(getZVariable(iech, ivar))) return true;
  return false;
}

bool Db::isAllUndefinedByType(const ELoc& loctype, int iech) const
{
  if (!isSampleIndexValid(iech)) return false;
  int natt = getLocNumber(loctype);
  if (natt <= 0) return false;

  for (int iatt = 0; iatt < natt; iatt++)
    if (! FFFF(getLocVariable(loctype, iech, iatt))) return true;
  return false;
}

int Db::getIntervalNumber() const
{
  return MAX(getLocNumber(ELoc::RKLOW), getLocNumber(ELoc::RKUP));
}

void Db::setInterval(int iech, int item, double rklow, double rkup)
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

void Db::setBound(int iech, int item, double lower, double upper)
{
  if (lower > upper)
  {
    messerr("Setting bounds: Lower (%lf) cannot be larger than upper (%lf)",
            lower,upper);
    return;
  }
  setLocVariable(ELoc::L,iech,item,lower);
  setLocVariable(ELoc::U,iech,item,upper);
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

VectorDouble Db::getGradient(int item, bool useSel) const
{
  if (!hasLocVariable(ELoc::G)) return VectorDouble();
  VectorDouble tab;

  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    if (useSel && ! isActive(iech)) continue;
    tab.push_back(getLocVariable(ELoc::G,iech,item));
  }
  return tab;
}

VectorDouble Db::getTangent(int item, bool useSel) const
{
  if (!hasLocVariable(ELoc::TGTE)) return VectorDouble();
  VectorDouble tab;

  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    if (useSel && ! isActive(iech)) continue;
    tab.push_back(getLocVariable(ELoc::TGTE,iech,item));
  }
  return tab;
}

/**
 * Return the Selection value at Sample 'iech'
 * @param iech Sample number
 * @return
 * @remark If the selection value if TEST, the sample is considered as masked off.
 */
int Db::getSelection(int iech) const
{
  if (!hasLocVariable(ELoc::SEL)) return 1;
  double value = getFromLocator(ELoc::SEL, iech, 0);
  if (FFFF(value)) return 0;
  int sel = (! isZero(value)) ? 1 :  0;
  return (sel);
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
  if (!hasLocVariable(ELoc::SEL)) return (getSampleNumber());

  /* Case when a selection is present */

  int count = 0;
  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    if (! isZero(getFromLocator(ELoc::SEL, iech, 0))) count++;
  }
  return count;
}

/**
 * Return the absolute rank of a sample from its relative rank
 * @param irel Relative rank
 * @return
 */
int Db::getRankRelativeToAbsolute(int irel) const
{
  if (! hasLocVariable(ELoc::SEL)) return irel;
  int nech = getSampleNumber(false);
  int jech = 0;
  for (int iabs = 0; iabs < nech; iabs++)
  {
    if (! isActive(iabs)) continue;
    if (irel == jech) return iabs;
    jech++;
  }
  return -1;
}

int Db::getRankAbsoluteToRelative(int iabs) const
{
  if (! hasLocVariable(ELoc::SEL)) return iabs;
  int nech = getSampleNumber(false);
  int irel = 0;
  for (int jabs = 0; jabs < nech; jabs++)
  {
    if (! isActive(jabs)) continue;
    if (jabs == iabs) return irel;
    irel++;
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
  if (!hasLocVariable(ELoc::SEL)) return _nech;

  if (!useSel) return _nech;
  int count = 0;
  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    if (!isZero(getFromLocator(ELoc::SEL, iech, 0))) count++;
  }
  return count;
}

/**
 * Returns the number of samples active and whose Z-value(item) is defined
 * @param item Rank of the Z-locator
 * @return
 */
int Db::getNumberActiveAndDefined(int item) const
{
  int count = 0;
  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    if (isActiveAndDefined(iech, item)) count++;
  }
  return count;
}

double Db::getWeight(int iech) const
{
  if (!hasLocVariable(ELoc::W)) return 1.;
  double w = getFromLocator(ELoc::W, iech, 0);
  if (FFFF(w)) w = 1.;
  if (w < 0) w = 0.;
  return (w);
}

VectorDouble Db::getWeights(bool useSel) const
{
  int icol = -1;
  int nech = getSampleNumber();
  VectorDouble sel;
  VectorDouble tab(nech);

  if (useSel) sel = getSelections();
  if (hasLocVariable(ELoc::W)) icol = getColIdxByLocator(ELoc::W, 0);

  int ecr = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    if (useSel && !sel.empty() && isZero(sel[iech])) continue;
    if (icol >= 0)
      tab[ecr] = getValueByColIdx(iech, icol);
    else
      tab[ecr] = 1.;
    ecr++;
  }
  tab.resize(ecr);
  return tab;
}

/****************************************************************************/
/*!
 **  Returns the list of Unique codes
 **
 ** \return  Pointer to the array containing a single occurence of each code
 **
 *****************************************************************************/
VectorDouble Db::getCodeList(void) const
{
  VectorDouble work(_nech);

  /* Load all the codes */

  int number = 0;
  for (int iech = 0; iech < _nech; iech++)
  {
    if (isActive(iech))
      work[number++] = getLocVariable(ELoc::C,iech,0);
  }

  /* Get the Unique occurrence */

  work.resize(number);
  VectorDouble tab = VH::unique(work);
  return (tab);
}

bool Db::isActiveDomain(int iech) const
{
  if (! hasLocVariable(ELoc::DOM)) return true;
  if (! GlobalEnvironment::getEnv()->isDomainReference()) return true;
  double value = getFromLocator(ELoc::DOM, iech, 0);
  if (FFFF(value)) return false;
  if (! GlobalEnvironment::getEnv()->matchDomainReference(value)) return true;
  return false;
}

/**
 * Returns the value of a simulation / variable for a given sample
 */
double Db::getSimvar(const ELoc& locatorType,
                     int iech,
                     int isimu,
                     int ivar,
                     int icase,
                     int nbsimu,
                     int nvar) const
{
  int item = getSimRank(isimu, ivar, icase, nbsimu, nvar);
  return getFromLocator(locatorType, iech, item);
}

/**
 * Set the value of a simulation / variable for a given sample
 */
void Db::setSimvar(const ELoc& locatorType,
                   int iech,
                   int isimu,
                   int ivar,
                   int icase,
                   int nbsimu,
                   int nvar,
                   double value)
{
  int item = getSimRank(isimu, ivar, icase, nbsimu, nvar);
  setFromLocator(locatorType, iech, item, value);
}

/**
 * Update the value of a simulation / variable for a given sample
 */
void Db::updSimvar(const ELoc& locatorType,
                   int iech,
                   int isimu,
                   int ivar,
                   int icase,
                   int nbsimu,
                   int nvar,
                   const EOperator& oper,
                   double value)
{
  int item = getSimRank(isimu, ivar, icase, nbsimu, nvar);

  // This direct addressing is meant to save time
  int icol = getColIdxByLocator(locatorType, item);
  if (icol < 0) return;
  int internalAddress = _getAddress(iech, icol);

  double oldval = _array[internalAddress];
  double newval = modifyOperator(oper, oldval, value);
  _array[internalAddress] = newval;
}

bool Db::isActive(int iech) const
{
  return (getSelection(iech) != 0 && isActiveDomain(iech));
}

bool Db::isActiveAndDefined(int iech, int item) const
{
  if (!isActive(iech)) return false;;
  return (! FFFF(getZVariable(iech, item)));
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
    if (FFFF(getZVariable(iech, item))) continue;
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
  if (number > size) return -1;
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
  int count = getFromLocatorNumber(locatorType);
  if (count <= 0) return namelist;
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
    if (icol < 0 || icol >= (int) _colNames.size()) continue;
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

VectorString Db::getName(const String& name) const
{
  return expandNameList(name);
}

VectorString Db::getNames(const VectorString& names) const
{
  return expandNameList(names);
}

VectorString Db::getAllNames(bool excludeRankAndCoordinates, bool verbose) const
{
  if (!excludeRankAndCoordinates) return _colNames;
  
  // From the list of all variables, exclude the following variables:
  // - the one named 'rank' (if any)
  // - the coordinates (if any)
  VectorString allnames = _colNames;
  VectorString names;
  for (int ivar = 0, nvar = (int) allnames.size(); ivar < nvar; ivar++)
  {
    // Exclude variable named 'rank'
    if (matchRegexp(allnames[ivar], "rank", false))
    {
      if (verbose) message("Excluding variable %s\n", allnames[ivar].c_str());
      continue;
    }

    // Exclude coordinates
    if (matchRegexp(allnames[ivar], "x*", false))
    {
      if (verbose) message("Excluding variable %s\n", allnames[ivar].c_str());
      continue;
    }

    // Add the names to the output list
    names.push_back(allnames[ivar]);
  }
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

void Db::setName(const VectorString& list, const String& name)
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
  int count = getFromLocatorNumber(locatorType);
  if (count <= 0) return;
  for (int i = 0; i < count; i++)
  {
    int icol = getColIdxByLocator(locatorType, i);
    if (icol < 0) continue;
    _colNames[icol] = incrementStringVersion(name, i+ 1);
  }
  correctNamesForDuplicates(_colNames);
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
  sstr << "Total number of samples      = " << getSampleNumber() << std::endl;
  if (hasLocVariable(ELoc::SEL))
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
    double vmin = VH::minimum(coor);
    double vmax = VH::maximum(coor);

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
 * @param mode 1 for basic statistics; 2 for class statistics
 * @param maxNClass Maximum number of printed classes
 * @return
 */
String Db::_summaryStats(VectorInt cols, int mode, int maxNClass) const
{
  std::stringstream sstr;

  int ncol = (cols.empty()) ? getColumnNumber() : static_cast<int> (cols.size());
  if (ncol <= 0) return sstr.str();

  sstr << toTitle(1, "Data Base Statistics");

  int nmask, ntest, nout;
  int nech = getSampleNumber(false);
  VectorDouble tab, wgt;

  // Loop on the columns

  for (int jcol = 0; jcol < ncol; jcol++)
  {
    int icol = (cols.empty()) ? jcol : cols[jcol];
    if (!isColIdxValid(icol)) continue;

    tab = getColumnByColIdx(icol, true);
    wgt = getWeights(true);
    StatResults stats = ut_statistics((int) tab.size(), tab.data(), NULL, wgt.data());

    sstr << icol + 1 << " - Name " << getNameByColIdx(icol) << " - Locator "
         << _getLocatorNameByColIdx(icol) << std::endl;
    sstr << " Nb of data          = " << toInt(nech) << std::endl;
    sstr << " Nb of active values = " << toInt(stats.nvalid) << std::endl;
    if (stats.nvalid <= 0) continue;

    /* Dispatch */

    if (mode == 1)
    {
      sstr << " Minimum value       = " << toDouble(stats.mini) << std::endl;
      sstr << " Maximum value       = " << toDouble(stats.maxi) << std::endl;
      sstr << " Mean value          = " << toDouble(stats.mean) << std::endl;
      sstr << " Standard Deviation  = " << toDouble(stats.stdv) << std::endl;
      sstr << " Variance            = " << toDouble(stats.stdv * stats.stdv) << std::endl;
    }
    else
    {
      double vmin = floor(stats.mini - 0.5);
      double vmax = ceil(stats.maxi + 0.5);
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
        sstr << " (" << toDouble(100. * classe[iclass] / stats.nvalid) << "%)";
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

  sstr << toMatrix(String(), colnames, VectorString(), true, number, ncol, tab);

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

VectorDouble Db::getSelections(void) const
{
  int nech = getSampleNumber();
  VectorDouble tab;

  if (!hasLocVariable(ELoc::SEL)) return tab;
  int icol = getColIdxByLocator(ELoc::SEL,0);
  if (!isColIdxValid(icol)) return tab;

  tab.resize(nech);
  for (int iech = 0; iech < nech; iech++)
    tab[iech] = getValueByColIdx(iech, icol);
  return tab;
}

/**
 * Returns a one_dimensional vector of values for valid samples for the set of
 * variables 'ivars'
 *
 * @param ivars   Vector giving the indices of the variables of interest
 * @param nbgh    Vector giving the ranks of the elligible samples (optional)
 * @param means   Vector of Means per variable (optional)
 * @param useSel  Discard the masked samples (if True)
 * @param useVerr Discard the samples where Verr (if it exists) is not correctly
 * defined
 *
 * @note: if the current 'db' has some Z-variable defined, only samples where
 * @note a variable is defined is considered (search for heterotopy).
 * @note: If argumennt 'Mean' is provided, the mean is subtracted from the output vector
 */

VectorDouble Db::getMultipleValuesActive(const VectorInt& ivars,
                                         const VectorInt& nbgh,
                                         const VectorDouble& means,
                                         bool useSel,
                                         bool useVerr) const
{
  VectorInt jvars = ivars;
  if (jvars.empty()) jvars = VH::sequence(getLocatorNumber(ELoc::Z));
  VectorDouble vec;
  const VectorVectorInt index = getMultipleRanksActive(jvars, nbgh, useSel, useVerr);
    
  int nvar = (int)jvars.size();
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int jvar = jvars[ivar];
    const VectorInt& local = index[ivar];
    for (int iech = 0, nech = (int)local.size(); iech < nech; iech++)
    {
      double value = getZVariable( iech, jvar);
      if (! means.empty()) value -= means[jvar];
      vec.push_back(value);
    }
  }
  return vec;
}

/**
 * Returns the list of indices 'index' for valid samples for the set of variables 'ivars'
 * as well as the count of samples (per variable)
 *
 * @param ivars   Vector giving the indices of the variables of interest
 * @param nbgh    Vector giving the ranks of the elligible samples (optional)
 * @param useSel  Discard the masked samples (if True)
 * @param useVerr Discard the samples where Verr (if it exists) is not correctly defined
 *
 * @note: if the current 'db' has some Z-variable defined, only samples where
 * @note a variable is defined is considered (search for heterotopy).
 */
VectorVectorInt Db::getMultipleRanksActive(const VectorInt &ivars,
                                           const VectorInt &nbgh,
                                           bool useSel,
                                           bool useVerr) const
{
  int nvar = (int) ivars.size();

  VectorVectorInt index(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int jvar = ivars[ivar];
    index[ivar] = getRanksActive(nbgh, jvar, useSel, useVerr);
  }
  return index;
}

VectorInt Db::getRanksActive(const VectorInt& nbgh, int item, bool useSel, bool useVerr) const
{
  double value;
  int nech_tot = getSampleNumber();

  // Create a vector of ranks of samples to be searched (using input 'nbgh' or not)
  VectorInt nbgh_init;
  if (nbgh.empty())
    nbgh_init = VH::sequence(nech_tot);
  else
    nbgh_init = nbgh;
  int nech_init = (int) nbgh_init.size();

  // Create the column index for the selection (only if 'useSel')
  int icol = (useSel) ? getColIdxByLocator(ELoc::SEL,0) : -1;

  // Update the search for variable, if no variable is defined
  if (getLocNumber(ELoc::Z) <= 0) item = -1;

  // Check the presence of variance of measurement error variable (only if 'useVerr')
  bool useV = false;
  if (useVerr && item >= 0)
  {
    if (getColIdxByLocator(ELoc::V, item) >= 0) useV = true;
  }

  // Constitute the resulting vector osf selected sample ranks
  VectorInt ranks;
  for (int jech = 0; jech < nech_init; jech++)
  {
    int iech = nbgh_init[jech];

    // Check against a possible selection
    if (icol >= 0)
    {
      value = getValueByColIdx(iech, icol);
      if (value <= 0) continue;
    }

    // Check against the existence of a target variable
    if (item >= 0)
    {
      value = getZVariable( iech, item);
      if (FFFF(value)) continue;
    }

    // Check against the validity of the Variance of Measurement Error variable
    if (useV)
    {
      value = getLocVariable(ELoc::V, iech, item);
      if (FFFF(value) || value < 0) continue;
    }

    // The sample is finally accepted
    ranks.push_back(iech);
  }
  return ranks;
}


/**
 *  Returns the column referred by its rank (0-based)
 *
 */
VectorDouble Db::getColumnByColIdx(int icol,
                                   bool useSel,
                                   bool flagCompress) const
{
  int nech = getSampleNumber(false);
  if (!isColIdxValid(icol)) return VectorDouble();

  VectorDouble tab(nech, TEST);
  VectorDouble sel;
  if (useSel) sel = getSelections();

  int ecr = 0;
  double value = TEST;
  for (int iech = 0; iech < nech; iech++)
  {
    bool defined = true;
    if (useSel && !sel.empty()) defined = (isOne(sel[iech]));
    if (! defined)
    {
      if (flagCompress) continue;
      value = TEST;
    }
    else
    {
      value = getValueByColIdx(iech, icol);
    }
    tab[ecr] = value;
    ecr++;
  }
  tab.resize(ecr);
  return tab;
}

/**
 * Returns a Column referred by its user-identification rank
 *
 */
VectorDouble Db::getColumnByUID(int iuid, bool useSel, bool flagCompress) const
{
  int icol = getColIdxByUID(iuid);
  if (icol < 0) return VectorDouble();
  return getColumnByColIdx(icol, useSel, flagCompress);
}

/**
 * Returns the contents of one Column identified by its locator type and item rank
 *
 */
VectorDouble Db::getColumnByLocator(const ELoc& locatorType,
                                   int locatorIndex,
                                   bool useSel,
                                   bool flagCompress) const
{
  int icol = getColIdxByLocator(locatorType, locatorIndex);
  if (icol < 0) return VectorDouble();
  return getColumnByColIdx(icol, useSel, flagCompress);
}

/**
 * Returns the contents of one Column identified by its name
 *
 */
VectorDouble Db::getColumn(const String &name,
                           bool useSel,
                           bool flagCompress) const
{
  VectorInt iuids = _ids(name, true);
  if (iuids.empty()) return VectorDouble();
  int icol = getColIdxByUID(iuids[0]);
  if (icol < 0) return VectorDouble();
  return getColumnByColIdx(icol, useSel, flagCompress);
}

/**
 * Returns the contents of a set of Columns identified by the locator type
 *
 */
VectorDouble Db::getColumnsByLocator(const ELoc &locatorType,
                                     bool useSel,
                                     bool flagCompress) const
{
  VectorString names = getNamesByLocator(locatorType);
  return getColumns(names, useSel, flagCompress);
}

/**
 * Returns the contents of a set of Columns identified by their user-identified ranks
 */
VectorDouble Db::getColumnsByUID(const VectorInt &iuids,
                                 bool useSel,
                                 bool flagCompress) const
{
  if (iuids.empty()) return VectorDouble();
  int nech = getSampleNumber(useSel);
  int nvar = static_cast<int> (iuids.size());
  VectorDouble retval(nvar * nech);

  /* Loop on the variables to be retrieved */

  int ecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    VectorDouble local = getColumnByUID(iuids[ivar], useSel, flagCompress);
    for (int iech = 0; iech < nech; iech++)
      retval[ecr++] = local[iech];
  }
  return retval;
}

/**
 * Returns the contents of a set of Columns specified by their ranks (0 based)
 *
 */
VectorDouble Db::getColumnsByColIdx(const VectorInt &icols,
                                    bool useSel,
                                    bool flagCompress) const
{
  int nech = getSampleNumber();
  int nvar = static_cast<int> (icols.size());
  VectorDouble retval(nvar * nech);

  /* Loop on the variables to be retrieved */

  int ecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    VectorDouble local = getColumnByColIdx(icols[ivar], useSel, flagCompress);
    if (local.empty()) continue;
    for (int iech = 0; iech < nech; iech++)
      retval[ecr++] = local[iech];
  }
  return retval;
}

/**
 * Returns the contents of a set of columns referred to by their rank interval (0 based)
 *
 */
VectorDouble Db::getColumnsByColIdxInterval(int icol_beg,
                                            int icol_end,
                                            bool useSel,
                                            bool flagCompress) const
{
  VectorInt icols;
  for (int icol = icol_beg; icol < icol_end; icol++)
    icols.push_back(icol);
  return getColumnsByColIdx(icols, useSel, flagCompress);
}

/**
 * Returns the contents of a set of columns specified by the interval of their user-identification ranks
 *
 */
VectorDouble Db::getColumnsByUIDRange(int iuid_beg,
                                      int iuid_end,
                                      bool useSel,
                                      bool flagCompress) const
{
  VectorInt iuids;
  for (int iuid = iuid_beg; iuid < iuid_end; iuid++)
    iuids.push_back(iuid);
  return getColumnsByUID(iuids, useSel, flagCompress);
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

VectorString Db::getItemNames(const VectorString& colnames) const
{
  return expandNameList(colnames);
}
VectorString Db::getItemNames(const String& colname) const
{
  return expandNameList(colname);
}
VectorString Db::getItemNames(const ELoc& locatorType) const
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

/**
 * Returns all the Columns contained in a Db
 *
 */
VectorDouble Db::getAllColumns(bool useSel, bool flagCompress) const
{
  VectorInt iuids = getAllUIDs();
  return getColumnsByUID(iuids, useSel, flagCompress);
}

/**
 * Setting the contents of all the Columns of a Db
 *
 * @param tabs Vector of vectors containing the values to be assigned
 */
void Db::setAllColumns(const VectorVectorDouble& tabs)
{
  VectorInt iuids = getAllUIDs();
  for (int iuid = 0; iuid < (int) iuids.size(); iuid++)
    setColumnByUID(tabs[iuid], iuids[iuid], false);
}

/**
 * Returns the contents of the Colmuns specified by their names
 *
 */
VectorDouble Db::getColumns(const VectorString &names,
                            bool useSel,
                            bool flagCompress) const
{
  if (names.empty()) return VectorDouble();
  VectorInt iuids =  _ids(names, false);
  return getColumnsByUID(iuids, useSel, flagCompress);
}

/**
 * Returns the contents of the Columns specified by their names (one variable per column)
 *
 */
VectorVectorDouble Db::getColumnsAsVVD(const VectorString &names,
                                       bool useSel,
                                       bool flagCompress) const
{
  VectorVectorDouble vec;
  if (names.empty()) return vec;
  VectorInt iuids =  _ids(names, false);

  for (int i = 0; i < (int) iuids.size(); i++)
    vec.push_back(getColumnByUID(iuids[i], useSel, flagCompress));
  return vec;
}

/**
 * Returns the contents of the columns specified by their names
 */
MatrixRectangular Db::getColumnsAsMatrix(const VectorString &names,
                                         bool useSel,
                                         bool flagCompress) const
{
  if (names.empty()) return MatrixRectangular();
  VectorInt iuids =  _ids(names, false);
  int nvar = (int) iuids.size();
  int nech = getSampleNumber(useSel && flagCompress);

  MatrixRectangular mat(nech,nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    VectorDouble X = getColumnByUID(iuids[ivar], useSel, flagCompress);
    mat.setColumn(ivar, X);
  }
  return mat;
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
  if (useSel) sel = getSelections();

  int ecr = 0;
  for (int iech = 0; iech < getSampleNumber(); iech++)
  {
    bool defined = true;
    if (useSel && !sel.empty()) defined = (isOne(sel[iech]));
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

  VectorInt iuids;
  for (unsigned int i = 0; i < names.size(); i++)
  {
    int icol = getRankInList(_colNames, names[i]);
    if (icol < 0) return VectorInt();
    int iuid = getUIDByColIdx(icol);
    if (iuid < 0) return VectorInt();
    iuids.push_back(iuid);
  }
  return iuids;
}

VectorInt Db::getUIDs(const VectorString& names) const
{
  if (names.empty()) return VectorInt();
  VectorInt iuids = _ids(names, false);
  return iuids;
}

VectorInt Db::getUIDsByLocator(const ELoc& locatorType) const
{
  VectorInt iuids;
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
  int ntab = static_cast<int>(tab.size()) / _nech;
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
 * Paint the column 'icol' with sample rank (1-based)
 * @param icol Index of the column to be painted (0-based)
 */
void Db::_createRank(int icol)
{
  int nech = getSampleNumber();
  for (int iech = 0; iech < nech; iech++)
    setArray(iech, icol, iech + 1);

  // Set the name

  _setNameByColIdx(icol, "rank");
}


/**
 * Create the sample rank variable (1-based) assuming that the Db is empty
 * @param nech Number of samples requested
 */
void Db::_addRank(int nech)
{
  if (getColumnNumber() > 0 || getSampleNumber() > 0)
  {
    messerr("Error: the Db should be empty in order to call _addRank. Nothing is done");
    return;
  }
  VectorDouble ranks = VH::sequence(1., (double) nech);
  addColumns(ranks, "rank");
}

void Db::_defineDefaultNames(int shift, const VectorString& names)
{
  int ncol = getColumnNumber() - shift;
  if (!names.empty())
  {
    if ((int) names.size() != ncol)
    {
      messerr("Argument 'names'(%d) must match the variables in 'tab'(%d)",
              (int) names.size(),ncol);
      messerr("Variables are not renamed");
    }
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
    my_throw("Error in the dimension of 'locatorNames'");

  ELoc locatorType;
  int locatorIndex, mult;
  for (int icol = 0; icol < ncol; icol++)
  {
    if (locatorIdentify(locatorNames[icol], &locatorType, &locatorIndex, &mult) == 0)
      setLocatorByUID(icol + shift, locatorType, locatorIndex);
  }
}

void Db::_defineDefaultLocatorsByNames(int shift, const VectorString& names)
{
  if (names.empty()) return;

  int ncol = getColumnNumber() - shift;
  if ((int) names.size() != ncol)
    my_throw("Error in the dimension of 'names'");

  ELoc locatorType;
  int locatorIndex, mult;
  for (int icol = 0; icol < ncol; icol++)
  {
    if (locatorIdentify(names[icol], &locatorType, &locatorIndex, &mult) == 0)
      setLocatorByUID(icol + shift, locatorType, locatorIndex);
  }
}

void Db::statisticsBySample(const VectorString& names,
                            const std::vector<EStatOption>& opers,
                            bool flagIso,
                            double proba,
                            double vmin,
                            double vmax,
                            const NamingConvention& namconv)
{
  DECLARE_UNUSED(flagIso);
  if (names.empty()) return;
  if (opers.empty()) return;

  VectorInt iuids = getUIDs(names);
  int noper = (int) opers.size();

  // Add the variables for PointWise statistics
  int iuidn = addColumnsByConstant(noper);
  if (iuidn < 0) return;

  VectorString nameloc = getNamesByUID(iuids);
  dbStatisticsVariables(this, nameloc, opers, iuidn, proba, vmin, vmax);

  namconv.setNamesAndLocators(this, iuidn);
  for (int i = 0; i < noper; i++)
  {
    const EStatOption& oper = opers[i];
    namconv.setNamesAndLocators(this, iuidn + i, oper.getKey());
  }
}

/**
 * The target variables are referred to by their user-designation ranks
 */
VectorDouble Db::statisticsMulti(const VectorString &names,
                                 bool flagIso,
                                 bool verbose,
                                 const String &title)
{
  if (names.empty()) return VectorDouble();
  Table table = dbStatisticsCorrel(this, names, flagIso);
  if (verbose)
  {
    table.setTitle(title);
    table.display();
  }
  return table.getValues();
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
int Db::getSimRank(int isimu, int ivar, int icase, int nbsimu, int nvar)
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
  VectorString locators = getLocators(true);
  VectorString names    = getName("*");
  
  bool ret = true;
  ret = ret && _recordWrite<int>(os, "Number of variables", ncol);
  ret = ret && _recordWrite<int>(os, "Number of samples", getSampleNumber());
  ret = ret && _recordWriteVec<String>(os, "Locators", locators);
  ret = ret && _recordWriteVec<String>(os, "Names", names);
  ret = ret && _commentWrite(os, "Array of values");
  for (int iech = 0, nech = getSampleNumber(); ret && iech < nech; iech++)
  {
    VectorDouble vals = getArrayBySample(iech);
    ret = ret && _recordWriteVec<double>(os, "", vals);
  }
  return ret;
}

bool Db::_deserialize(std::istream& is, bool /*verbose*/)
{
  int ncol = 0;
  int nech = 0;
  VectorString locators;
  VectorString names;

  // Read the file
  bool ret = true;
  ret = ret && _recordRead<int>(is, "Number of variables", ncol);
  ret = ret && _recordRead<int>(is, "Number of samples", nech);
  if (! ret) return ret;
  if (ncol > 0)
  {
    ret = ret && _recordReadVec<String>(is, "Locators", locators, ncol);
    ret = ret && _recordReadVec<String>(is, "Names", names, ncol);
  }

  VectorDouble allvalues(nech * ncol);
  VectorDouble::iterator it(allvalues.begin());
  for (int iech = 0; iech < nech && ret; iech++)
  {
    ret = ret && _recordReadVecInPlace<double>(is, "Array of values", it, ncol);
  }

  if (ret)
  {
    // Decode the locators
    std::vector<ELoc> tabloc;
    VectorInt tabnum;
    int  inum = 0, mult = 0;
    ELoc iloc;
    for (const auto& loc : locators)
    {
      if (locatorIdentify(loc, &iloc, &inum, &mult) != 0) return true;
      tabloc.push_back(iloc);
      tabnum.push_back(inum);
    }

    // Initialize the Db
    resetDims(ncol, nech);

    // Load the values
    _loadData(ELoadBy::SAMPLE, false, allvalues);

    // Update the column names and locators
    for (int i = 0; i < ncol; i++)
    {
      setNameByUID(i, names[i]);
      setLocatorByUID(i, tabloc[i], tabnum[i]);
    }
  }
  return ret;
}

void Db::_loadData(const ELoadBy& order, bool flagAddSampleRank, const VectorDouble& tab)
{
  // Preliminary check

  if (getColumnNumber() <= 0) return;
  int jcol = 0;

  // Add the rank (optional)

  if (flagAddSampleRank)
  {
    for (int iech = 0; iech < getSampleNumber(); iech++)
      setValueByColIdx(iech, jcol, iech + 1);
    setNameByUID(jcol, "rank");
    jcol++;
  }

  // Add the input array 'tab' (if provided)

  if (tab.empty()) return;
  int ntab = (flagAddSampleRank) ? getColumnNumber() - 1 : getColumnNumber();
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
}

bool Db::_isCountValid(const VectorInt& iuids, bool flagOne, bool verbose) const
{
  if (iuids.empty() && flagOne)
  {
    if (verbose) messerr("No variable name corresponding to your criterion");
    return false;
  }
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
    int ifac = (int) getZVariable(iech,0);
    if (ifac <= 0) continue;
    if (ifac > nfac) nfac = ifac;
  }
  return nfac;
}

VectorBool Db::getActiveArray() const
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

  VH::arrangeInPlace(0, rindex, xval, true, nech);

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
 * @param flagAddSampleRank true if the sample rank must be generated
 *
 * @remark A possible selection in 'dbin' will not be taken into account
 * @remark You can use either 'proportion' or 'number'
 */
int Db::resetSamplingDb(const Db* dbin,
                        double proportion,
                        int number,
                        const VectorString& names,
                        int seed,
                        bool verbose,
                        bool flagAddSampleRank)
{
  if (proportion <= 0. && number <= 0)
  {
    messerr("You must specify either 'proportion' or 'number'");
    return 1;
  }
  _clear();

  // Creating the vector of selected samples

  int nfrom = dbin->getSampleNumber();
  VectorInt ranks = VH::sampleRanks(nfrom, proportion, number, seed);
  _nech = static_cast<int> (ranks.size());
  if (verbose)
    message("From %d samples, the extraction concerns %d samples\n", nfrom,_nech);

  // Creating the vector of variables

  VectorString namloc = names;
  if (namloc.empty())
    namloc = dbin->getAllNames();
  _ncol = static_cast<int> (namloc.size());

  // Create the (empty) architecture

  int ncol = (flagAddSampleRank) ? _ncol + 1: _ncol;
  resetDims(ncol, _nech);

  if (flagAddSampleRank) _createRank(0);

  // Define the variables and the Locators

  _defineVariableAndLocators(dbin, namloc, (int) flagAddSampleRank);

  // Load samples

  _loadValues(dbin, namloc, ranks, (int) flagAddSampleRank);

  return 0;
}

/**
 * Define the Name and Locator of the variables, coming from the list of variables
 * contained in 'dbin' whose names are contained in 'names'.
 * @param dbin Input Db
 * @param names List of variables to be copied
 * @param shift Shift when storing the first attribute
 */
void Db::_defineVariableAndLocators(const Db* dbin, const VectorString& names, int shift)
{
  ELoc locatorType;
  int locatorIndex;
  for (int icol = 0, ncol=(int) names.size(); icol < ncol; icol++)
  {
    setNameByUID(icol + shift, names[icol]);
    if (dbin->getLocator(names[icol],&locatorType,&locatorIndex))
      setLocator(names[icol],locatorType,locatorIndex);
  }
}

/**
 * Load values of the variables 'names' of 'dbin' into the current Db
 * This copy is restricted to only active samples of 'dbin' (given by 'ranks')
 * @param db    Input Db
 * @param names List of variables to be copied
 * @param ranks Vector of ranks of the samples to be copied
 * @param shift Shift when storing the first vector of attributes
 */
void Db::_loadValues(const Db* db, const VectorString& names, const VectorInt& ranks, int shift)
{
  VectorDouble values(_nech);
  for (int icol = 0, ncol= (int) names.size(); icol < ncol; icol++)
  {
    int jcol = db->getColIdx(names[icol]);
    for (int iech = 0; iech < _nech; iech++)
      values[iech] = db->getValueByColIdx(ranks[iech],jcol);
    setColumnByColIdx(values, icol + shift);
  }
}

int Db::resetReduce(const Db *dbin,
                    const VectorString &names,
                    const VectorInt &ranks,
                    bool verbose)
{
  // Creating the vector of selected samples

  VectorInt ranksel = ranks;
  if (ranksel.empty())
  {
    if (dbin->hasLocVariable(ELoc::SEL))
      ranksel = dbin->getRanksActive();
    else
      ranksel = VH::sequence(dbin->getSampleNumber());
  }
  _nech = static_cast<int> (ranksel.size());
  bool flagMask = _nech != dbin->getSampleNumber();
  if (verbose)
    message("From %d samples, the extraction concerns %d samples\n", dbin->getSampleNumber(),_nech);

  // Creating the vector of variables

  VectorString namloc = names;
  if (namloc.empty())
    namloc = dbin->getAllNames();

  // Create the (empty) architecture

  _ncol = static_cast<int> (namloc.size());
  resetDims(_ncol, _nech);

  // Define the variables and the Locators

  _defineVariableAndLocators(dbin, namloc);

  // Load samples

  _loadValues(dbin, namloc, ranksel);

  // When the number of coordinates is 0 and if the input Db is a grid,
  // Create the coordinates before running the reduction.
  // Otherwise, the resulting Db (which is a 'point' Db) will have no coordinate
  // the coordinates are added before reduction

  if (getLocatorNumber(ELoc::X) <= 0)
  {
    // Extract vector of coordinates from input 'Db' (converted into a 'DbGrid')
    const DbGrid* dbgrid = dynamic_cast<const DbGrid*>(dbin);
    if (dbgrid != nullptr)
    {
      int ndim = dbin->getNDim();
      VectorVectorDouble coors = dbgrid->getAllCoordinates();
      VectorString namloc = generateMultipleNames("Coor", ndim);

      // Save the coordinates in the output file (after possible sample selection)
      for (int idim = 0; idim < ndim; idim++)
      {
        if (flagMask)
        {
          VectorDouble coor = VH::compress(coors[idim], ranksel);
          addColumns(coor, namloc[idim], ELoc::X, idim);
        }
        else
        {
          addColumns(coors[idim], namloc[idim], ELoc::X, idim);
        }
      }
    }
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
 ** \param[in]  flagAddSampleRank true if the Rank must be generated in the output Db
 **
 *****************************************************************************/
Db* Db::createFromDbGrid(int nech,
                         DbGrid* dbgrid,
                         int seed,
                         bool flag_exact,
                         bool flag_repulsion,
                         double range,
                         double beta,
                         bool flagAddSampleRank)
{
  Db* db = db_point_init(nech, VectorDouble(), VectorDouble(), dbgrid,
                         flag_exact, flag_repulsion, range, beta,
                         0., seed, flagAddSampleRank);
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
  {
    return;
  }

  if (combine == "not")
  {
    for (int iech = 0; iech < nech; iech++)
      sel[iech] = 1. - sel[iech];
    return;
  }

  // Read an already existing selection
  VectorDouble oldsel = getColumnByLocator(ELoc::SEL, 0);
  if (oldsel.empty()) return;

  if (combine == "or")
  {
    for (int iech = 0; iech < nech; iech++)
      sel[iech] = sel[iech] || oldsel[iech];
    return;
  }
  if (combine == "and")
  {
    for (int iech = 0; iech < nech; iech++)
      sel[iech] = sel[iech] && oldsel[iech];
    return;
  }
  if (combine == "xor")
  {
    for (int iech = 0; iech < nech; iech++)
      sel[iech] = !areEqual(sel[iech], oldsel[iech]);
    return;
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
                          bool flagAddSampleRank)
{
  Db* db = new Db;
  if (db->resetFromSamples(nech, order, tab, names, locatorNames, flagAddSampleRank) != 0)
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
                      bool flagAddSampleRank)
{
  Db* db = new Db;
  if (db->resetFromCSV(filename, verbose, csv, ncol_max, nrow_max, flagAddSampleRank) != 0)
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
                      double extend,
                      bool flagAddSampleRank)
{
  Db* db = db_point_init(nech, coormin, coormax, nullptr,
                         flag_exact, flag_repulsion, range, beta, extend, seed,
                         flagAddSampleRank);
  return db;
}

Db* Db::createFromOnePoint(const VectorDouble& tab, bool flagAddSampleRank)
{
  Db* db = new Db;
  if (db->resetFromOnePoint(tab, flagAddSampleRank) != 0)
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
                         bool verbose,
                         bool flagAddSampleRank)
{
  Db* db = new Db;
  if (db->resetSamplingDb(dbin, proportion, number, names, seed, verbose, flagAddSampleRank) != 0)
  {
    messerr("Error when clearing Db by Sampling another Db");
    delete db;
    return nullptr;
  }
  return db;
}

Db* Db::createReduce(const Db *dbin,
                     const VectorString &names,
                     const VectorInt &ranks,
                     bool verbose)
{
  Db* db = new Db;
  if (db->resetReduce(dbin, names, ranks, verbose) != 0)
  {
    db = dbin->clone();
  }
  return db;
}

bool Db::areSame(const String& name1,
                 const String& name2,
                 double eps,
                 bool useSel,
                 bool verbose) const
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

VectorInt Db::shrinkToValidRows(const VectorInt& rows) const
{
  if (rows.empty()) return rows;

  VectorInt new_rows;
  for (int i = 0; i < (int) rows.size(); i++)
  {
    if (rows[i] >= 0 && rows[i] < _nech) new_rows.push_back(rows[i]);
  }
  return new_rows;
}

VectorInt Db::shrinkToValidCols(const VectorInt& cols) const
{
  if (cols.empty()) return cols;

  VectorInt new_cols;
  for (int i = 0; i < (int) cols.size(); i++)
  {
    if (cols[i] >= 0 && cols[i] < _ncol) new_cols.push_back(cols[i]);
  }
  return new_cols;
}

/** Returns the ranks, within the exhaustive loop on variables then samples
 * to be used in the kriging system, knowing that we must discard the masked samples
 * and the samples whose value is not defined
 * @return
 */
VectorInt Db::getSampleRanks() const
{
  VectorInt vec;
  int nvar = getLocNumber(ELoc::Z);
  int nech = getSampleNumber();

  int lec = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int iech = 0; iech < nech; iech++, lec++)
    {
      if (! isActive(iech)) continue;
      if (FFFF(getZVariable( iech, ivar))) continue;
      vec.push_back(lec);
    }
  return vec;
}

/**
 * Creating a new Db loaded with random values
 * @param ndat Number of samples
 * @param ndim Dimension of the space
 * @param nvar Number of variables
 * @param nfex Number of external drift functions
 * @param ncode Number of codes (no code when 0)
 * @param varmax Maximum value for the measurement error
 * @param selRatio Percentage of samples that must be masked off
 * @param heteroRatio Vector of proportions of NA to be generated per variable
 * @param coormin Vector of minima of the rectangle containing data (0s if not defined)
 * @param coormax Vector of maxima of the rectangle containing data (1s if not defined)
 * @param seed Value for the Random Generator seed
 * @param flagAddSampleRank true if the sample rank must be generated
 * @return A pointer to the newly created Db
 *
 * @remarks
 * The coordinates are generated uniformly within [0,1]
 * The variance of measurement error is created only if 'varmax' is positive. Then a field is created
 * for each variable. this field is filled with random values uniformly generated in [0, varmax]
 * The external drift values are generated according to Normal distribution.
 */
Db* Db::createFillRandom(int ndat,
                         int ndim,
                         int nvar,
                         int nfex,
                         int ncode,
                         double varmax,
                         double selRatio,
                         const VectorDouble& heteroRatio,
                         const VectorDouble& coormin,
                         const VectorDouble& coormax,
                         int seed,
                         bool flagAddSampleRank)
{
  // Set the seed
  law_set_random_seed(seed);

  // Create the Db
  Db* db = Db::create();

  // Add the sample rank attribute
  if (flagAddSampleRank) db->_addRank(ndat);

  // Generate the vector of coordinates
  VectorVectorDouble coor(ndim);
  for (int idim = 0; idim < ndim; idim++)
  {
    double mini = ((int) coormin.size() == ndim) ? coormin[idim] : 0.;
    double maxi = ((int) coormax.size() == ndim) ? coormax[idim] : 1.;
    coor[idim] = VH::simulateUniform(ndat, mini, maxi);
  }
  db->addColumnsByVVD(coor, "x", ELoc::X);

  // Generate the Vectors of Variance of measurement error (optional)
  if (varmax > 0.)
  {
    VectorVectorDouble varm(nvar);
    for (int ivar = 0; ivar < nvar; ivar++)
      varm[ivar] = VH::simulateUniform(ndat, 0., varmax);
    db->addColumnsByVVD(varm, "v", ELoc::V);
  }

  // Generate the External Drift functions (optional)
  if (nfex > 0)
  {
    VectorVectorDouble fex(nfex);
    for (int ifex = 0; ifex < nfex; ifex++)
      fex[ifex] = VH::simulateGaussian(ndat);
    db->addColumnsByVVD(fex, "f", ELoc::F);
  }

  // Generate the selection (optional)
  if (selRatio > 0)
  {
    VectorDouble sel(ndat);
    VectorDouble rnd = VH::simulateUniform(ndat);
    for (int idat = 0; idat < ndat; idat++)
      sel[idat] = (rnd[idat] > selRatio) ? 1. : 0.;
    db->addColumns(sel, "sel", ELoc::SEL);
  }

  // Generate the variables
  bool flag_hetero = ((int) heteroRatio.size() == nvar);
  VectorVectorDouble vars(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    vars[ivar] = VH::simulateGaussian(ndat);
    if (flag_hetero)
    {
      VectorDouble rnd = VH::simulateUniform(ndat);
      for (int idat = 0; idat < ndat; idat++)
        if (rnd[idat] <= heteroRatio[ivar])
          vars[ivar][idat] = TEST;
    }
  }
  db->addColumnsByVVD(vars, "z", ELoc::Z);

  // Generate the code (optional)
  if (ncode > 0)
  {
    VectorDouble codes = VH::simulateUniform(ndat);
    for (int idat = 0; idat < ndat; idat++)
      codes[idat] = int(codes[idat] * (1.+ncode));
    db->addColumns(codes, "code", ELoc::C);
  }

  return db;
}

Table Db::printOneSample(int iech,
                         const VectorString& names,
                         bool excludeCoordinates,
                         bool skipTitle) const
{
  Table table;
  VectorString allNames = names;
  if (allNames.empty()) allNames = getAllNames(excludeCoordinates);
  VectorString localNames = expandNameList(allNames);

  const int nvar = (int)localNames.size();
  if (nvar <= 0) return table;
  if (! isSampleIndexValid(iech)) return table;

  table.reset(nvar, 1);
  table.setSkipDescription(true);
  if (! skipTitle)
    table.setTitle("Sample " + std::to_string(iech+1) + " / " + std::to_string(getSampleNumber()));

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    table.setRowName(ivar, localNames[ivar]);
    table.setValue(ivar, 0, getValue(localNames[ivar], iech));
  }
  return table;
}