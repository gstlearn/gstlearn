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
#include "Db/DbGrid.hpp"

namespace gstlrn
{
class ELoc;

class GSTLEARN_EXPORT ACalcDbToDb: public ACalculator
{
public:
  ACalcDbToDb(bool mustShareSameSpaceDimension = true);
  ACalcDbToDb(const ACalcDbToDb& r)            = delete;
  ACalcDbToDb& operator=(const ACalcDbToDb& r) = delete;
  virtual ~ACalcDbToDb();

  void setDbin(Db* dbin) { _dbin = dbin; }
  void setDbout(Db* dbout) { _dbout = dbout; }
  void setNamingConvention(const NamingConvention& namconv) { _namconv = namconv; }
  void setMustShareSpaceDimension(bool mustShareSpaceDimension)
  {
    _mustShareSpaceDimension = mustShareSpaceDimension;
  }

  Db* getDbin() const { return _dbin; }
  Db* getDbout() const { return _dbout; }
  DbGrid* getGridin() const;
  DbGrid* getGridout() const;
  const NamingConvention& getNamingConvention() const { return _namconv; }

  bool hasDbin(bool verbose = true) const;
  bool hasDbout(bool verbose = true) const;
  bool isGridIn(bool verbose = true) const;
  bool isGridOut(bool verbose = true) const;

protected:
  bool _check() override;
  bool _preprocess() override;
  Id _getNDim() const { return _ndim; }
  Id _getNVar() const { return _nvar; }
  bool _setNdim(Id ndim, bool flagForce = false);
  bool _setNvar(Id nvar, bool flagForce = false);

  Id _addVariableDb(Id whichDb,
                    Id status,
                    const ELoc& locatorType,
                    Id locatorIndex = 0,
                    Id number       = 1,
                    double valinit  = 0.);
  void _renameVariable(Id whichDb,
                       const VectorString& names,
                       const ELoc& locatorType,
                       Id nvar,
                       Id iptr,
                       const String& qualifier,
                       Id count,
                       bool flagSetLocator = true,
                       Id locatorShift     = 0);
  void _storeInVariableList(Id whichDb, Id status, const VectorInt& iuids);
  Id _expandInformation(Id mode, const ELoc& locatorType) const;
  void _cleanVariableDb(Id status);
  Db* _whichDb(Id whichDb);
  String _identifyVariable(Id iuid) const;

private:
  bool _checkSpaceDimension();
  bool _checkVariableNumber();

private:
  bool _mustShareSpaceDimension;
  Db* _dbin;
  Db* _dbout;
  NamingConvention _namconv;
  VectorInt _listVariablePermDbIn;
  VectorInt _listVariablePermDbOut;
  VectorInt _listVariableTempDbIn;
  VectorInt _listVariableTempDbOut;
  Id _ndim;
  Id _nvar;
};
} // namespace gstlrn