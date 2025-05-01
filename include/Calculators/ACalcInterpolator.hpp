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

#include "Model/ModelGeneric.hpp"
#include "gstlearn_export.hpp"

#include "Calculators/ACalcDbToDb.hpp"

class ELoc;
// class Model;
class ANeigh;

// TODO : Create InterpolatorParam ASpaceParam which inherits from ASPaceObject and AParam, which inherits from ASerializable, AStringable, IClonable
class GSTLEARN_EXPORT ACalcInterpolator: public ACalcDbToDb
{
public:
  ACalcInterpolator();
  ACalcInterpolator(const ACalcInterpolator &r) = delete;
  ACalcInterpolator& operator=(const ACalcInterpolator &r) = delete;
  virtual ~ACalcInterpolator();

  void setModel(ModelGeneric* model) { _model = model; }
  void setNeigh(ANeigh* neigh) { _neigh = neigh; }
  void setKrigopt(const KrigOpt& krigopt) { _krigopt = krigopt; }

  ModelGeneric*  getModel() const { return _model; }
  ANeigh* getNeigh() const { return _neigh; }
  const KrigOpt& getKrigopt() const { return _krigopt; }

  bool hasModel(bool verbose = true) const;
  bool hasNeigh(bool verbose = true) const;

protected:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  int  _getNCov() const { return _ncova; }
  bool _setNCov(int ncova);

  int _centerDataToGrid(DbGrid* dbgrid);

private:
  ModelGeneric*  _model;
  ANeigh* _neigh;
  KrigOpt _krigopt;
  int _ncova;
};
