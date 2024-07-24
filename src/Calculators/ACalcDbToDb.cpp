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
#include "Calculators/ACalcDbToDb.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/VectorHelper.hpp"

ACalcDbToDb::ACalcDbToDb(bool mustShareSameSpaceDimension)
    : ACalculator(),
      _mustShareSpaceDimension(mustShareSameSpaceDimension),
      _dbin(nullptr),
      _dbout(nullptr),
      _namconv(),
      _listVariablePermDbIn(),
      _listVariablePermDbOut(),
      _listVariableTempDbIn(),
      _listVariableTempDbOut()
{
}

ACalcDbToDb::~ACalcDbToDb()
{
}

int ACalcDbToDb::_getNDim() const
{
  if (_dbin != nullptr)
  {
    return _dbin->getNDim();
  }

  if (_dbout != nullptr)
  {
    return _dbout->getNDim();
  }
  return -1;
}

int ACalcDbToDb::_getNVar() const
{
  int nvar = 0;
  if (_dbin != nullptr)
  {
    if (nvar > 0)
    {
      if (nvar != _dbin->getLocNumber(ELoc::Z)) return -1;
    }
    else
    {
      nvar = _dbin->getLocNumber(ELoc::Z);
    }
  }
  return nvar;
}

bool ACalcDbToDb::_checkSpaceDimension()
{
  if (! _mustShareSpaceDimension) return true;

  /**************************************************/
   /* Cross-checking the Space Dimension consistency */
   /**************************************************/

   int ndim = 0;
   if (_dbin != nullptr)
   {
     if (ndim > 0)
     {
       if (ndim != _dbin->getNDim())
       {
         messerr("Inconsistent Space dimension:");
         messerr("- Current dimension = %d", ndim);
         messerr("- Space Dimension of 'dbin' = %d", _dbin->getNDim());
         return false;
       }
     }
     else
     {
       ndim = _dbin->getNDim();
     }
   }

   if (_dbout != nullptr)
   {
     if (ndim > 0)
     {
       if (ndim != _dbout->getNDim())
       {
         messerr("Inconsistent Space dimension:");
         messerr("- Current dimension = %d", ndim);
         messerr("- Space Dimension of 'dbout' = %d", _dbout->getNDim());
         return false;
       }
     }
     else
     {
//       ndim = _dbout->getNDim(); // Never reached
     }
   }
   return true;
}

bool ACalcDbToDb::_checkVariableNumber()
{
  int nvar = 0;
  if (_dbin != nullptr)
  {
    if (nvar > 0)
    {
      if (nvar != _dbin->getLocNumber(ELoc::Z))
      {
        messerr("Inconsistent the Variable Number:");
        messerr("- Current number = %d", nvar);
        messerr("- Number of variables in 'dbin' = %d",
                _dbin->getLocNumber(ELoc::Z));
        return false;
      }
    }
    else
    {
//      nvar = _dbin->getLocNumber(ELoc::Z); // Never reached
    }
  }
  return true;
}

bool ACalcDbToDb::_check()
{
  /**************************************************/
  /* Cross-checking the Space Dimension consistency */
  /**************************************************/

  if (! _checkSpaceDimension()) return false;

  /**************************************************/
  /* Cross-Checking the Variable Number consistency */
  /**************************************************/

  if (! _checkVariableNumber()) return false;

  return true;
}

/**
 * Returns a pointer to the relevant Db and issue a message if not defined
 * @param whichDb 1 for 'dbin'  and 2 for 'dbout'
 * @return A pointer to the Db or nullptr
 */
Db* ACalcDbToDb::_whichDb(int whichDb)
{
  Db *db;
  if (whichDb == 1)
    db = _dbin;
  else
    db = _dbout;
  if (db == nullptr)
  {
    messerr("Impossible to add variables in non-defined Db");
  }
  return db;
}

String ACalcDbToDb::_identifyVariable(int iuid) const
{
  return _dbin->getNameByUID(iuid);
}

/**
 * Store the IUID of the new variable in the relevant internal list
 * @param whichDb 1 for variable in 'dbin'; 2 for variable in 'dbout'
 * @param status  1 for variables to be stored; 2 for Temporary variable
 * @param iuids   Vector of UIDs of the new variable
 */
void ACalcDbToDb::_storeInVariableList(int whichDb,
                                       int status,
                                       const VectorInt &iuids)
{
  int number = (int) iuids.size();
  if (number <= 0) return;

  if (whichDb == 1)
  {
    if (status == 1)
    {
      for (int i = 0; i < number; i++)
        _listVariablePermDbIn.push_back(iuids[i]);
    }
    else
    {
      for (int i = 0; i < number; i++)
        _listVariableTempDbIn.push_back(iuids[i]);
    }
  }
  else
  {
    if (status == 1)
    {
      for (int i = 0; i < number; i++)
        _listVariablePermDbOut.push_back(iuids[i]);
    }
    else
    {
      for (int i = 0; i < number; i++)
        _listVariableTempDbOut.push_back(iuids[i]);
    }
  }
}
int ACalcDbToDb::_addVariableDb(int whichDb,
                                int status,
                                const ELoc& locatorType,
                                int locatorIndex,
                                int number,
                                double valinit)
{
  Db *db = _whichDb(whichDb);
  if (db == nullptr) return -1;
  int iuid = db->addColumnsByConstant(number, valinit, String(), locatorType,
                                      locatorIndex);
  if (iuid < 0) return -1;
  VectorInt iuids = VH::sequence(number, iuid);
  _storeInVariableList(whichDb, status, iuids);
  return iuid;
}

/**
 * Define the characteristics of the variables created by a Db2Db Calculator
 * @param whichDb     1 if the variable belongs to 'dbin'; 2 if it belongs to 'dbout'
 * @param names       Names of the variables in 'dbin' (or empty)
 * @param locatorType Locator for the names of input variables (or ELoc::UNKNOWN)
 * @param nvar        Number of variables (when constructed from locator)
 * @param iptr        IUID of the (first) variable to be renamed
 * @param qualifier   Name which will serve as 'qualifier' (when provided)
 * @param count       Number of variable named from the same basic name (using version number)
 * @param flagSetLocator True if the locator must be defined
 * @param locatorShift Shift to calculate the rank of the locator currently defined
 */
void ACalcDbToDb::_renameVariable(int whichDb,
                                  const VectorString& names,
                                  const ELoc& locatorType,
                                  int nvar,
                                  int iptr,
                                  const String& qualifier,
                                  int count,
                                  bool flagSetLocator,
                                  int locatorShift)
{
  if (whichDb == 1)
    _namconv.setNamesAndLocators(_dbin, names, locatorType, nvar,
                                 _dbin, iptr, qualifier, count, flagSetLocator,
                                 locatorShift);
  else
    _namconv.setNamesAndLocators(_dbin, names, locatorType, nvar,
                                 _dbout, iptr, qualifier, count, flagSetLocator,
                                 locatorShift);
}

void ACalcDbToDb::_cleanVariableDb(int status)
{
  // Dispatch

  if (status == 1)
  {
    // In 'dbin'
    if (!_listVariablePermDbIn.empty())
    {
      for (int i = 0; i < (int) _listVariablePermDbIn.size(); i++)
        _dbin->deleteColumnByUID(_listVariablePermDbIn[i]);
    }
    _listVariablePermDbIn.clear();

    // In 'dbout'
    if (!_listVariablePermDbOut.empty())
    {
      for (int i = 0; i < (int) _listVariablePermDbOut.size(); i++)
      {
        _dbout->deleteColumnByUID(_listVariablePermDbOut[i]);
      }
    }
    _listVariablePermDbOut.clear();
  }
  else
  {
    // In 'dbin'
    if (!_listVariableTempDbIn.empty())
    {
      for (int i = 0; i < (int) _listVariableTempDbIn.size(); i++)
        _dbin->deleteColumnByUID(_listVariableTempDbIn[i]);
    }
    _listVariableTempDbIn.clear();

    // In 'dbout'
    if (!_listVariableTempDbOut.empty())
    {
      for (int i = 0; i < (int) _listVariableTempDbOut.size(); i++)
        _dbout->deleteColumnByUID(_listVariableTempDbOut[i]);
    }
    _listVariableTempDbOut.clear();
  }
}

bool ACalcDbToDb::hasDbin(bool verbose) const
{
  if (_dbin == nullptr)
  {
    if (verbose) messerr("The argument 'dbin' must be defined");
    return false;
  }
  return true;
}
bool ACalcDbToDb::hasDbout(bool verbose) const
{
  if (_dbout == nullptr)
  {
    if (verbose) messerr("The argument 'dbout' must be defined");
    return false;
  }
  return true;
}

bool ACalcDbToDb::isGridIn(bool verbose) const
{
  if (!hasDbin(false)) return false;
  if (!_dbin->isGrid())
  {
    if (verbose) messerr("The argument 'dbin' should be a Grid File");
    return false;
  }
  return true;
}

bool ACalcDbToDb::isGridOut(bool verbose) const
{
  if (!hasDbout(false)) return false;
  if (!_dbout->isGrid())
  {
    if (verbose) messerr("The argument 'dbout' should be a Grid File");
    return false;
  }
  return true;
}

DbGrid* ACalcDbToDb::getGridin() const
{
  if (!hasDbin(false)) return nullptr;
  DbGrid *dbgrid = dynamic_cast<DbGrid*>(_dbin);
  return dbgrid;
}

DbGrid* ACalcDbToDb::getGridout() const
{
  if (!hasDbout(false)) return nullptr;
  DbGrid *dbgrid = dynamic_cast<DbGrid*>(_dbout);
  return dbgrid;
}

/*****************************************************************************/
/*!
 **  Derive some information from the Output db (if Grid) to the Input Db
 **
 ** \return  Error return code
 **
 ** \param[in]  mode        1 for allocation; -1 for deallocation
 ** \param[in]  locatorType Type of the pointer (ELoc)
 **
 ** \remark This function is only valid when the Output Db is a grid
 ** \remark However, in case of a Point output Db, this function should not
 ** \remark be used: the external drift functions should already be present
 ** \remark in the output Db.
 ** \remark If this is not the case, an error is issued.
 **
 ** \remark When called with mode=1, the new variables are added
 ** \remark When called with mode=-1, the variables are deleted (by type)
 **
 *****************************************************************************/
int ACalcDbToDb::_expandInformation(int mode, const ELoc& locatorType)
{
  if (getDbin() == nullptr || getDbout() == nullptr) return 0;

  // Check the number of fields to be expanded

  int ninfo;
  if (getDbout()->isGrid() && locatorType == ELoc::X)
    ninfo = getDbout()->getNDim();
  else
    ninfo = getDbout()->getFromLocatorNumber(locatorType);
  if (ninfo <= 0) return 0;

  // Check the corresponding number of variables in the Input File

  int ninfoIn = getDbin()->getFromLocatorNumber(locatorType);
  if (ninfo == ninfoIn) return 0;

  /* Case when the Output Db is not a grid */

  if (!getDbout()->isGrid())
  {
    messerr("The Output Db is not a Grid file");
    messerr("The Input Db does not contain the correct number of External Drifts");
    return 1;
  }
  DbGrid *dbgrid = dynamic_cast<DbGrid*>(getDbout());

  /* Dispatch */

  if (mode > 0)
  {
    // Here, the naming Convention is modified in order to anticipate the
    // locator of the newly created variables
    NamingConvention* namconv = NamingConvention::create("Migrate");
    namconv->setLocatorOutType(locatorType);
    int error = migrateByLocator(dbgrid, getDbin(), locatorType, 1,
                                 VectorDouble(), false, false, false, *namconv);
    delete namconv;
    if (error != 0) return 1;
  }
  else
  {
    getDbin()->deleteColumnsByLocator(locatorType);
  }
  return 0;
}

