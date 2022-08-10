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
#pragma once

#include "gstlearn_export.hpp"

#include "Calculators/ACalculator.hpp"
#include "Db/Db.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/Vector.hpp"

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
                     const ELoc &locatorType,
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
