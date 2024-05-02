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

#include "Calculators/ACalcDbToDb.hpp"

class ELoc;
class Model;
class ANeigh;

// TODO : Create InterpolatorParam ASpaceParam which inherits from ASPaceObject and AParam, which inherits from ASerializable, AStringable, IClonable
class GSTLEARN_EXPORT ACalcInterpolator: public ACalcDbToDb
{
public:
  ACalcInterpolator();
  ACalcInterpolator(const ACalcInterpolator &r) = delete;
  ACalcInterpolator& operator=(const ACalcInterpolator &r) = delete;
  virtual ~ACalcInterpolator();

  void setModel(Model *model)  { _model = model; }
  void setNeigh(ANeigh *neigh) { _neigh = neigh; }

  Model*  getModel() const { return _model; }
  ANeigh* getNeigh() const { return _neigh; }

  bool hasModel(bool verbose = true) const;
  bool hasNeigh(bool verbose = true) const;

protected:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual int  _getNDim() const override;
  virtual int  _getNVar() const override;
  virtual int  _getNCova() const;

  int _centerDataToGrid(DbGrid* dbgrid);

private:
  Model* _model;
  ANeigh* _neigh;
};
