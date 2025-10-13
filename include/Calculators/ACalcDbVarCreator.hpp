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
#pragma once

#include "gstlearn_export.hpp"

#include "Calculators/ACalculator.hpp"
#include "Db/Db.hpp"
#include "Basic/NamingConvention.hpp"

namespace gstlrn
{
class ELoc;

class GSTLEARN_EXPORT ACalcDbVarCreator: public ACalculator
{
public:
  ACalcDbVarCreator();
  ACalcDbVarCreator(const ACalcDbVarCreator &r) = delete;
  ACalcDbVarCreator& operator=(const ACalcDbVarCreator &r) = delete;
  virtual ~ACalcDbVarCreator();

  void setDb(Db* db) { _db = db; }
  void setNamingConvention(const NamingConvention& namconv) { _namconv = namconv; }

  Db*  getDb() const { return _db; }
  bool hasDb(bool verbose = false) const;

  const NamingConvention& getNamconv() const
  {
    return _namconv;
  }

protected:
  Id _getNDim() const;
  Id _getNVar() const;

  Id _addVariableDb(Id status,
                     const ELoc& locatorType,
                     Id locatorIndex = 0,
                     Id number = 1,
                     double valinit = 0.);
  void _renameVariable(Id nvar,
                       Id iptr,
                       const ELoc& locatorInType,
                       const String &qualifier,
                       Id count);
  void _storeInVariableList(Id status, const VectorInt& iuids);
  void _cleanVariableDb(Id status);

private:
  Db* _db;
  NamingConvention _namconv;
  VectorInt _listVariablePermDb;
  VectorInt _listVariableTempDb;
};
}