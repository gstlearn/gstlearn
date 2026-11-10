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

  ACalcDbVarCreator::~ACalcDbVarCreator() {}

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
   * @param colids  Vector of Column descriptors of the new variable
   */
  void ACalcDbVarCreator::_storeInVariableList(Id status, const ColIDs& colids)
  {
    if (colids.empty()) return;
    auto& target = (status == 1) ? _listVariablePermDb : _listVariableTempDb;
    target.insert(target.end(), colids.begin(), colids.end());
  }

  Id ACalcDbVarCreator::_addVariableDb(
    Id status,
    const ELoc& locatorType,
    Id locatorIndex,
    Id number,
    Id nversion,
    double valinit)
  {
    DECLARE_UNUSED(locatorIndex);
    if (_db == nullptr) return -1;
    auto colids = _db->addMultiColumns(
      VectorString(), valinit, VectorDouble(), toERole(locatorType), nversion,
      false, false, 0, 0, number);
    if (colids.empty()) return -1;
    _storeInVariableList(status, colids);
    return colids[0].getUID();
  }

  void ACalcDbVarCreator::_renameVariable(
    Id nvar,
    Id iptr,
    const ELoc& locatorInType,
    const String& qualifier,
    Id count)
  {
    auto names = _db->getNames(toERole(locatorInType));
    _namconv.setOutput(names, nvar, _db, iptr, qualifier, count);
  }

  void ACalcDbVarCreator::_cleanVariableDb(Id status)
  {
    if (_db == nullptr) return;

    auto& list = (status == 1) ? _listVariablePermDb : _listVariableTempDb;
    // Deletion of the variables by reverse order to avoid any problem with the shifting of the columns
    for (auto colid = list.rbegin(); colid != list.rend(); ++colid)
      _db->deleteColumn(*colid);
    list.clear();
  }

  bool ACalcDbVarCreator::hasDb(bool verbose) const
  {
    if (_db == nullptr)
    {
      if (verbose) messerr("The argument 'db' must be defined");
      return false;
    }
    return true;
  }
} // namespace gstlrn
