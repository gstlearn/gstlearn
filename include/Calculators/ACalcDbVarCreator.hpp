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

#include "Basic/NamingConvention.hpp"
#include "Calculators/ACalculator.hpp"
#include "Db/Db.hpp"

namespace gstlrn
{
  class ELoc;

  class GSTLEARN_EXPORT ACalcDbVarCreator: public ACalculator
  {
  public:
    ACalcDbVarCreator();
    ACalcDbVarCreator(const ACalcDbVarCreator& r) = delete;
    ACalcDbVarCreator& operator=(const ACalcDbVarCreator& r) = delete;
    virtual ~ACalcDbVarCreator();

    void setDb(Db* db) { _db = db; }

    void setNamingConvention(const NamingConvention& namconv)
    {
      _namconv = namconv;
    }

    Db* getDb() const { return _db; }

    bool hasDb(bool verbose = false) const;

    const NamingConvention& getNamconv() const { return _namconv; }

  protected:
    Id _getNDim() const;
    Id _getNVar() const;

    void _renameVariable(
      Id nvar,
      Id iptr,
      const ELoc& locatorInType,
      const String& qualifier,
      Id count);
    void _storeInVariableList(Id status, const ColIDs& colids);
    void _cleanVariableDb(Id status);

  private:
    Db* _db;
    NamingConvention _namconv;
    ColIDs _listVariablePermDb;
    ColIDs _listVariableTempDb;
  };
} // namespace gstlrn
