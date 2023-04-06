/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f_private.h"

#include "Enum/ELoc.hpp"

#include "Basic/NamingConvention.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Calculators/ACalcDbToDb.hpp"
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"

#include <math.h>

CalcMigrate::CalcMigrate()
    : ACalcDbToDb(false),
      _iattOut(-1),
      _iuids(),
      _ldmax(1),
      _dmax(),
      _flagFill(false),
      _flagInter(false),
      _flagLocate(false),
      _locatorType(ELoc::Z)
{
}

CalcMigrate::~CalcMigrate()
{
}

bool CalcMigrate::_check()
{
  if (! ACalcDbToDb::_check()) return false;

  if (! hasDbin()) return false;
  if (! hasDbout()) return false;

  if (_iuids.empty())
  {
    messerr("At least one variable should be defined");
    return false;
  }
  if (_ldmax != 1 && _ldmax != 2)
  {
    messerr("Argument 'ldmax'(%d)  should be 1 (for L1 distance) or 2 (for L2 distance",_ldmax);
    return false;
  }
  return true;
}

bool CalcMigrate::_preprocess()
{
  int nvar = _getNVar();
  _iattOut = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, nvar, 0.);
  if (_iattOut < 0) return false;
  return true;
}

bool CalcMigrate::_postprocess()
{
  /* Free the temporary variables */
  _cleanVariableDb(2);

  int nvar = _getNVar();
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    _renameVariable(2, 1, _iattOut+ivar, _identifyVariable(_iuids[ivar]), 1, ! _flagLocate, ivar);
  }

  if (_flagLocate)
    getDbout()->setLocatorsByUID(nvar, _iattOut, _locatorType);

  return true;
}

void CalcMigrate::_rollback()
{
  _cleanVariableDb(1);
}

int CalcMigrate::_getNVar() const
{
  return (int) _iuids.size();
}

/****************************************************************************/
/*!
 **  Standard Kriging
 **
 ** \return  Error return code
 **
 *****************************************************************************/
bool CalcMigrate::_run()
{
  int nvar = _getNVar();

  // Perform the migrations

  for (int i = 0; i < nvar; i++)
  {
    int iatt1 = _iuids[i];
    int iatt2 = _iattOut + i;
    if (_migrate(getDbin(), getDbout(), iatt1, iatt2, _ldmax, _dmax, _flagFill,
                 _flagInter)) return false;
  }

  return true;
}

/*****************************************************************************/
/*!
 **  Migrates a variable from one Db to another one
 **
 ** \return  Error return code
 **
 ** \param[in]  db1        Descriptor of the input Db
 ** \param[in]  db2        Descriptor of the output Db
 ** \param[in]  name       Name of the attribute to be migrated
 ** \param[in]  ldmax      Type of distance for calculating maximum distance
 **                        1 for L1 and 2 for L2 distance
 ** \param[in]  dmax       Array of maximum distances (optional)
 ** \param[in]  flag_fill  Filling option
 ** \param[in]  flag_inter Interpolation
 ** \param[in]  namconv    Naming convention
 **
 *****************************************************************************/
int migrate(Db *db1,
            Db *db2,
            const String &name,
            int ldmax,
            const VectorDouble &dmax,
            int flag_fill,
            int flag_inter,
            const NamingConvention &namconv)
{
  CalcMigrate migrate;
  migrate.setDbin(db1);
  migrate.setDbout(db2);
  migrate.setNamingConvention(namconv);

  VectorInt iuids(1);
  iuids[0] = db1->getUID(name);
  migrate.setIuids(iuids);
  migrate.setLdmax(ldmax);
  migrate.setDmax(dmax);
  migrate.setFlagFill(flag_fill);
  migrate.setFlagInter(flag_inter);

  // Run the calculator
  int error = (migrate.run()) ? 0 : 1;

  return error;
}

/*****************************************************************************/
/*!
 **  Migrates a set of variables from one Db to another one
 **
 ** \return  Error return code
 **
 ** \param[in]  db1        Descriptor of the input Db
 ** \param[in]  db2        Descriptor of the output Db
 ** \param[in]  names      Name of the attribute to be migrated
 ** \param[in]  ldmax      Type of distance for calculating maximum distance
 **                        1 for L1 and 2 for L2 distance
 ** \param[in]  dmax       Array of maximum distances (optional)
 ** \param[in]  flag_fill  Filling option
 ** \param[in]  flag_inter Interpolation
 ** \param[in]  namconv    Naming convention
 **
 *****************************************************************************/
int migrateVariables(Db *db1,
                     Db *db2,
                     const VectorString &names,
                     int ldmax,
                     const VectorDouble &dmax,
                     int flag_fill,
                     int flag_inter,
                     const NamingConvention &namconv)
{
  CalcMigrate migrate;
  migrate.setDbin(db1);
  migrate.setDbout(db2);
  migrate.setNamingConvention(namconv);

  VectorInt iuids = db1->getUIDs(names);
  migrate.setIuids(iuids);
  migrate.setLdmax(ldmax);
  migrate.setDmax(dmax);
  migrate.setFlagFill(flag_fill);
  migrate.setFlagInter(flag_inter);

  // Run the calculator
  int error = (migrate.run()) ? 0 : 1;
  return error;
}

/*****************************************************************************/
/*!
 **  Migrates a variable from one Db to another one
 **
 ** \return  Error return code
 **
 ** \param[in]  db1        Descriptor of the input Db
 ** \param[in]  db2        Descriptor of the output Db
 ** \param[in]  atts       Array of attributes to be migrated
 ** \param[in]  ldmax      Type of distance for calculating maximum distance
 **                        1 for L1 and 2 for L2 distance
 ** \param[in]  dmax       Array of maximum distances (optional)
 ** \param[in]  flag_fill  Filling option
 ** \param[in]  flag_inter Interpolation
 ** \param[in]  namconv    Naming Convention
 **
 *****************************************************************************/
int migrateByAttribute(Db *db1,
                       Db *db2,
                       const VectorInt& atts,
                       int ldmax,
                       const VectorDouble &dmax,
                       int flag_fill,
                       int flag_inter,
                       const NamingConvention &namconv)
{
  CalcMigrate migrate;
  migrate.setDbin(db1);
  migrate.setDbout(db2);
  migrate.setNamingConvention(namconv);

  VectorInt iuids = atts;
  if (iuids.empty()) iuids = db1->getAllUIDs();

  migrate.setIuids(iuids);
  migrate.setLdmax(ldmax);
  migrate.setDmax(dmax);
  migrate.setFlagFill(flag_fill);
  migrate.setFlagInter(flag_inter);

  // Run the calculator
  int error = (migrate.run()) ? 0 : 1;
  return error;
}

/*****************************************************************************/
/*!
 **  Migrates all z-locator variables from one Db to another one
 **
 ** \return  Error return code
 **
 ** \param[in]  db1         Descriptor of the input Db
 ** \param[in]  db2         Descriptor of the output Db
 ** \param[in]  locatorType Locator Type
 ** \param[in]  ldmax       Type of distance for calculating maximum distance
 **                         1 for L1 and 2 for L2 distance
 ** \param[in]  dmax        Array of maximum distances (optional)
 ** \param[in]  flag_fill   Filling option
 ** \param[in]  flag_inter  Interpolation
 ** \param[in]  namconv     Naming convention
 **
 ** \remark The output variable receive the same locator as the input variables
 **
 *****************************************************************************/
int migrateByLocator(Db *db1,
                     Db *db2,
                     const ELoc& locatorType,
                     int ldmax,
                     const VectorDouble &dmax,
                     int flag_fill,
                     int flag_inter,
                     const NamingConvention &namconv)
{
  CalcMigrate migrate;
  migrate.setDbin(db1);
  migrate.setDbout(db2);
  migrate.setNamingConvention(namconv);

  VectorString names = db1->getNamesByLocator(locatorType);
  VectorInt iuids = db1->getUIDs(names);
  migrate.setIuids(iuids);
  migrate.setLdmax(ldmax);
  migrate.setDmax(dmax);
  migrate.setFlagFill(flag_fill);
  migrate.setFlagInter(flag_inter);
  migrate.setFlagLocate(true);
  migrate.setLocatorType(locatorType);

  // Run the calculator
  int error = (migrate.run()) ? 0 : 1;
  return error;
}

