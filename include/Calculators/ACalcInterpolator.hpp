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

protected:
  virtual bool _check() const override;
  int _expandInformation(int mode, const ELoc &locatorType);

private:
  Db* _dbin;
  Db* _dbout;
  Model* _model;
  ANeighParam* _neighparam;
};
