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
#include "Calculators/ACalcDbVarCreator.hpp"
#include "Basic/VectorHelper.hpp"
#include "Calculators/ACalculator.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"

namespace gstlrn
{
ACalcDbVarCreator::ACalcDbVarCreator()
  : ACalculator()
  , _db(nullptr)
  , _namconv()
  , _listVariablePermDb()
  , _listVariableTempDb()
{
}

ACalcDbVarCreator::~ACalcDbVarCreator()
{
}

Id ACalcDbVarCreator::_getNDim() const
{
  if (_db == nullptr) return -1;
  return _db->getNDim();
}

Id ACalcDbVarCreator::_getNVar() const
{
  if (_db == nullptr) return -1;
  return _db->getNLoc(ELoc::Z);
}

/**
 * Store the IUID of the new variable in the relevant internal list
 * @param status  1 for variables to be stored; 2 for Temporary variable
 * @param iuids   Vector of UIDs of the new variable
 */
void ACalcDbVarCreator::_storeInVariableList(Id status, const VectorInt& iuids)
{
  Id number = static_cast<Id>(iuids.size());
  if (number <= 0) return;

  if (status == 1)
  {
    for (Id i = 0; i < number; i++)
      _listVariablePermDb.push_back(iuids[i]);
  }
  else
  {
    for (Id i = 0; i < number; i++)
      _listVariableTempDb.push_back(iuids[i]);
  }
}

Id ACalcDbVarCreator::_addVariableDb(Id status,
                                     const ELoc& locatorType,
                                     Id locatorIndex,
                                     Id number,
                                     double valinit)
{
  if (_db == nullptr) return -1;
  Id iuid = _db->addColumnsByConstant(number, valinit, String(), locatorType, locatorIndex);
  if (iuid < 0) return -1;
  VectorInt iuids = VH::sequence(number, iuid);
  _storeInVariableList(status, iuids);
  return iuid;
}

void ACalcDbVarCreator::_renameVariable(Id nvar,
                                        Id iptr,
                                        const ELoc& locatorInType,
                                        const String& qualifier,
                                        Id count)
{
  _namconv.setNamesAndLocators(_db, VectorString(), locatorInType, nvar, _db, iptr, qualifier, count);
}

void ACalcDbVarCreator::_cleanVariableDb(Id status)
{
  // Dispatch

  if (status == 1)
  {
    if (!_listVariablePermDb.empty())
    {
      for (Id i = 0; i < static_cast<Id>(_listVariablePermDb.size()); i++)
        _db->deleteColumnByUID(_listVariablePermDb[i]);
    }
    _listVariablePermDb.clear();
  }
  else
  {
    if (!_listVariableTempDb.empty())
    {
      for (Id i = 0; i < static_cast<Id>(_listVariableTempDb.size()); i++)
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
} // namespace gstlrn