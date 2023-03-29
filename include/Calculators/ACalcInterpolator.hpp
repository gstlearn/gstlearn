/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Calculators/ACalcDbToDb.hpp"
#include "Calculators/ACalculator.hpp"
#include "Model/Model.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Basic/NamingConvention.hpp"

class ELoc;

// TODO : Create InterpolatorParam ASpaceParam which inherits from ASPaceObject and AParam, which inherits from ASerializable, AStringable, IClonable
class GSTLEARN_EXPORT ACalcInterpolator: public ACalcDbToDb
{
public:
  ACalcInterpolator();
  ACalcInterpolator(const ACalcInterpolator &r) = delete;
  ACalcInterpolator& operator=(const ACalcInterpolator &r) = delete;
  virtual ~ACalcInterpolator();

  void setModel(Model *model) { _model = model; }
  void setNeighparam(ANeighParam *neighparam) { _neighparam = neighparam; }

  Model* getModel() const { return _model; }
  ANeighParam* getNeighparam() const { return _neighparam; }

  bool hasModel(bool verbose = true) const;
  bool hasNeighParam(bool verbose = true) const;

protected:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual int _getNDim() const override;
  virtual int _getNVar() const override;
  virtual int _getNCova() const;

  int _centerDataToGrid(DbGrid* dbgrid);

private:
  Model* _model;
  ANeighParam* _neighparam;
};
