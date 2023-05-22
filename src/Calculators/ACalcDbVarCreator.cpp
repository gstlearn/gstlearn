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
#include "Calculators/ACalcDbVarCreator.hpp"
#include "Calculators/ACalculator.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Basic/VectorHelper.hpp"

ACalcDbVarCreator::ACalcDbVarCreator()
    : ACalculator(),
      _db(nullptr),
      _namconv(),
      _listVariablePermDb(),
      _listVariableTempDb()
{
}

ACalcDbVarCreator::~ACalcDbVarCreator()
{
}

int ACalcDbVarCreator::_getNDim() const
{
  if (_db == nullptr) return -1;
  return  _db->getNDim();
}

int ACalcDbVarCreator::_getNVar() const
{
  if (_db == nullptr) return -1;
  return  _db->getLocatorNumber(ELoc::Z);
}

/**
 * Store the IUID of the new variable in the relevant internal list
 * @param status  1 for variables to be stored; 2 for Temporary variable
 * @param iuids   Vector of UIDs of the new variable
 */
void ACalcDbVarCreator::_storeInVariableList(int status, const VectorInt& iuids)
{
  int number = (int) iuids.size();
  if (number <= 0) return;

  if (status == 1)
  {
    for (int i = 0; i < number; i++)
      _listVariablePermDb.push_back(iuids[i]);
  }
  else
  {
    for (int i = 0; i < number; i++)
      _listVariableTempDb.push_back(iuids[i]);
  }
}

int ACalcDbVarCreator::_addVariableDb(int status,
                                      const ELoc& locatorType,
                                      int locatorIndex,
                                      int number,
                                      double valinit)
{
  if (_db == nullptr) return -1;
  int iuid = _db->addColumnsByConstant(number, valinit, String(), locatorType, locatorIndex);
  if (iuid < 0) return -1;
  VectorInt iuids = VH::sequence(number, iuid);
  _storeInVariableList(status, iuids);
  return iuid;
}

void ACalcDbVarCreator::_renameVariable(int nvar,
                                        int iptr,
                                        const String &name,
                                        int count)
{
  _namconv.setNamesAndLocators(_db, ELoc::Z, nvar, _db, iptr, name, count);
}

void ACalcDbVarCreator::_cleanVariableDb(int status)
{
  // Dispatch

  if (status == 1)
  {
    if (!_listVariablePermDb.empty())
    {
      for (int i = 0; i < (int) _listVariablePermDb.size(); i++)
        _db->deleteColumnByUID(_listVariablePermDb[i]);
    }
    _listVariablePermDb.clear();
  }
  else
  {
    if (!_listVariableTempDb.empty())
    {
      for (int i = 0; i < (int) _listVariableTempDb.size(); i++)
        _db->deleteColumnByUID(_listVariableTempDb[i]);
    }
    _listVariableTempDb.clear();
  }
}

bool ACalcDbVarCreator::hasDb(bool verbose) const
{
  if (_db == nullptr)
  {
    if (verbose)
      messerr("The argument 'db' must be defined");
    return false;
  }
  return true;
}
