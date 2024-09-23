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
#include "Db/DbGrid.hpp"
#include "Basic/NamingConvention.hpp"

class ELoc;

class GSTLEARN_EXPORT ACalcDbToDb: public ACalculator
{
public:
  ACalcDbToDb(bool mustShareSameSpaceDimension = true);
  ACalcDbToDb(const ACalcDbToDb &r) = delete;
  ACalcDbToDb& operator=(const ACalcDbToDb &r) = delete;
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

  bool hasDbin(bool verbose = true) const;
  bool hasDbout(bool verbose = true) const;
  bool isGridIn(bool verbose = true) const;
  bool isGridOut(bool verbose = true) const;

protected:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  int _getNDim() const { return _ndim; }
  int _getNVar() const { return _nvar; }
  bool _setNdim(int ndim, bool flagForce = false);
  bool _setNvar(int nvar, bool flagForce = false);

  int _addVariableDb(int whichDb,
                     int status,
                     const ELoc& locatorType,
                     int locatorIndex = 0,
                     int number = 1,
                     double valinit = 0.);
  void _renameVariable(int whichDb,
                       const VectorString& names,
                       const ELoc& locatorType,
                       int nvar,
                       int iptr,
                       const String& qualifier,
                       int count,
                       bool flagSetLocator = true,
                       int locatorShift    = 0);
  void _storeInVariableList(int whichDb, int status, const VectorInt& iuids);
  int  _expandInformation(int mode, const ELoc& locatorType) const;
  void _cleanVariableDb(int status);
  Db*  _whichDb(int whichDb);
  String _identifyVariable(int iuid) const;

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
  int _ndim;
  int _nvar;
};
