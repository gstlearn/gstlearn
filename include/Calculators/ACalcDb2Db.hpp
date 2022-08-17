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

class GSTLEARN_EXPORT ACalcDb2Db: public ACalculator
{
public:
  ACalcDb2Db();
  ACalcDb2Db(const ACalcDb2Db &r) = delete;
  ACalcDb2Db& operator=(const ACalcDb2Db &r) = delete;
  virtual ~ACalcDb2Db();

  void setDbin(Db* dbin) { _dbin = dbin; }
  void setDbout(Db* dbout) { _dbout = dbout; }
  void setNamingConvention(const NamingConvention& namconv) { _namconv = namconv; }

  Db* getDbin() const { return _dbin; }
  Db* getDbout() const { return _dbout; }

  bool hasDbin(bool verbose = true) const;
  bool hasDbout(bool verbose = true) const;

protected:
  virtual bool _check() override;
  virtual int _getNDim() const;
  virtual int _getNVar() const;

  int _addVariableDb(int whichDb,
                     int status,
                     const ELoc &locatorType,
                     int number = 1,
                     double valinit = 0.);
  void _renameVariable(int nvar,
                       int iptr,
                       const String &name,
                       int count,
                       bool flagSetLocator = true);
  void _storeInVariableList(int whichDb, int status, const VectorInt& iuids);
  void _cleanVariableDb(int status);
  int  _expandInformation(int mode, const ELoc &locatorType);
  Db*  _whichDb(int whichDb);

private:
  Db* _dbin;
  Db* _dbout;
  NamingConvention _namconv;
  VectorInt _listVariablePermDbIn;
  VectorInt _listVariablePermDbOut;
  VectorInt _listVariableTempDbIn;
  VectorInt _listVariableTempDbOut;
};
