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

#include "Simulation/ACalcSimulation.hpp"
#include "Simulation/SimuSubstitutionParam.hpp"
#include "Basic/Plane.hpp"

class Db;
class DbGrid;

class GSTLEARN_EXPORT CalcSimuSubstitution: public ACalcSimulation
{
public:
  CalcSimuSubstitution(int nbsimu = 0, int seed = 4324324, bool verbose = false);
  CalcSimuSubstitution(const CalcSimuSubstitution &r) = delete;
  CalcSimuSubstitution& operator=(const CalcSimuSubstitution &r) = delete;
  virtual ~CalcSimuSubstitution();

  const SimuSubstitutionParam& getSubparam() const { return _subparam; }
  void setSubparam(const SimuSubstitutionParam &subparam) { _subparam = subparam; }

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;

  bool _simulate();
  void _calculValue(int ip, double factor, const VectorDouble& vector);
  VectorDouble _transToProp(const SimuSubstitutionParam& subparam,
                            bool verbose = false,
                            double eps = EPSILON5);

private:
  bool _verbose;
  int  _iattOut;
  SimuSubstitutionParam _subparam;
  std::vector<Plane> _planes;
};

GSTLEARN_EXPORT int substitution(DbGrid *dbgrid,
                                 SimuSubstitutionParam& subparam,
                                 int seed = 43242,
                                 int verbose = false,
                                 const NamingConvention& namconv = NamingConvention("SimSub"));
