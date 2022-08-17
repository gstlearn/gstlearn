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

#include "geoslib_define.h"

#include "Calculators/ACalcDb2Db.hpp"

class Db;
class DbGrid;
class EStatOption;

class GSTLEARN_EXPORT CalcStatistics: public ACalcDb2Db
{
public:
  CalcStatistics();
  CalcStatistics(const CalcStatistics &r) = delete;
  CalcStatistics& operator=(const CalcStatistics &r) = delete;
  virtual ~CalcStatistics();

  void setRadius(int radius) { _radius = radius; }
  void setOper(const EStatOption &oper) { _oper = oper; }

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;

private:
  int    _iattOut;
  EStatOption _oper;
  int _radius;
};

GSTLEARN_EXPORT int dbStatisticsOnGrid(Db *db,
                                       DbGrid *dbgrid,
                                       const EStatOption &oper,
                                       int radius = 0,
                                       const NamingConvention &namconv = NamingConvention(
                                           "Stats"));
