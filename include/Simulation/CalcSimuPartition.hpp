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

#include "Basic/Plane.hpp"
#include "Simulation/ACalcSimulation.hpp"
#include "Simulation/SimuPartitionParam.hpp"

class Db;
class DbGrid;

typedef struct
{
  double valref;
  double valsim;
} Stack;

class GSTLEARN_EXPORT CalcSimuPartition: public ACalcSimulation
{
public:
  CalcSimuPartition(int mode, int nbsimu = 0, int seed = 4324324, bool verbose = false);
  CalcSimuPartition(const CalcSimuPartition &r) = delete;
  CalcSimuPartition& operator=(const CalcSimuPartition &r) = delete;
  virtual ~CalcSimuPartition();

  const SimuPartitionParam& getParparam() const { return _parparam; }
  void setParparam(const SimuPartitionParam& parparam) { _parparam = parparam; }

  void setMode(int mode) { _mode = mode; }
  void setVerbose(bool verbose) { _verbose = verbose; }

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;

  bool _voronoi();
  bool _poisson();
  double _stackSearch(const std::vector<Stack>& stacks, double valref);

private:
  int _mode;
  bool _verbose;
  int _iattOut;
  SimuPartitionParam _parparam;
};

GSTLEARN_EXPORT int tessellation_voronoi(DbGrid *dbgrid,
                                         Model *model,
                                         const SimuPartitionParam& parparam,
                                         int seed = 43243,
                                         int verbose = false,
                                         const NamingConvention& namconv = NamingConvention(
                                             "Voronoi"));
GSTLEARN_EXPORT int tessellation_poisson(DbGrid *dbgrid,
                                         Model *model,
                                         const SimuPartitionParam& parparam,
                                         int seed = 432432,
                                         int verbose = false,
                                         const NamingConvention& namconv = NamingConvention(
                                             "Poisson"));
