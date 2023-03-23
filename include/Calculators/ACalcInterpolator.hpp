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
