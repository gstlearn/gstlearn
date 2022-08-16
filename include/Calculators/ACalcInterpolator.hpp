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
#include "Model/Model.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/Vector.hpp"

class ELoc;

class GSTLEARN_EXPORT ACalcInterpolator: public ACalculator
{
public:
  ACalcInterpolator();
  ACalcInterpolator(const ACalcInterpolator &r) = delete;
  ACalcInterpolator& operator=(const ACalcInterpolator &r) = delete;
  virtual ~ACalcInterpolator();

  void setDbin(Db* dbin) { _dbin = dbin; }
  void setDbout(Db* dbout) { _dbout = dbout; }
  void setModel(Model *model) { _model = model; }
  void setNeighparam(ANeighParam *neighparam) { _neighparam = neighparam; }
  void setNamingConvention(const NamingConvention& namconv) { _namconv = namconv; }

  Db* getDbin() const { return _dbin; }
  Db* getDbout() const { return _dbout; }
  Model* getModel() const { return _model; }
  ANeighParam* getNeighparam() const { return _neighparam; }

  bool hasDbin(bool verbose = true) const;
  bool hasDbout(bool verbose = true) const;
  bool hasModel(bool verbose = true) const;
  bool hasNeighParam(bool verbose = true) const;

protected:
  virtual bool _check() override;
  int _getNDim() const;
  int _getNCova() const;
  virtual int _getNVar() const;

  int _expandInformation(int mode, const ELoc &locatorType);
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

  Db* _whichDb(int whichDb);

private:
  Db* _dbin;
  Db* _dbout;
  Model* _model;
  ANeighParam* _neighparam;
  NamingConvention _namconv;
  VectorInt _listVariablePermDbIn;
  VectorInt _listVariablePermDbOut;
  VectorInt _listVariableTempDbIn;
  VectorInt _listVariableTempDbOut;
};
