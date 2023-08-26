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
  int _getNDim() const;
  int _getNVar() const;

  int _addVariableDb(int status,
                     const ELoc& locatorType,
                     int locatorIndex = 0,
                     int number = 1,
                     double valinit = 0.);
  void _renameVariable(int nvar,
                       int iptr,
                       const String &name,
                       int count);
  void _storeInVariableList(int status, const VectorInt& iuids);
  void _cleanVariableDb(int status);

private:
  Db* _db;
  NamingConvention _namconv;
  VectorInt _listVariablePermDb;
  VectorInt _listVariableTempDb;
};
