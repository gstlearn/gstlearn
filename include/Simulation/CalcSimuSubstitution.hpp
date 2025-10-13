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

#include "geoslib_define.h"

#include "Basic/Plane.hpp"
#include "Simulation/ACalcSimulation.hpp"
#include "Simulation/SimuSubstitutionParam.hpp"

namespace gstlrn
{

class Db;
class DbGrid;

class GSTLEARN_EXPORT CalcSimuSubstitution: public ACalcSimulation
{
public:
  CalcSimuSubstitution(Id nbsimu = 0, Id seed = 4324324, bool verbose = false);
  CalcSimuSubstitution(const CalcSimuSubstitution& r)            = delete;
  CalcSimuSubstitution& operator=(const CalcSimuSubstitution& r) = delete;
  virtual ~CalcSimuSubstitution();

  const SimuSubstitutionParam& getSubparam() const { return _subparam; }
  void setSubparam(const SimuSubstitutionParam& subparam) { _subparam = subparam; }

private:
  bool _check() override;
  bool _preprocess() override;
  bool _run() override;
  bool _postprocess() override;
  void _rollback() override;

  bool _simulate();
  void _calculValue(Id ip, double factor, const VectorDouble& vector);
  static VectorDouble _transToProp(const SimuSubstitutionParam& subparam,
                                   bool verbose = false,
                                   double eps   = EPSILON5);

private:
  bool _verbose;
  Id _iattOut;
  SimuSubstitutionParam _subparam;
  std::vector<Plane> _planes;
};

GSTLEARN_EXPORT Id substitution(DbGrid* dbgrid,
                                 SimuSubstitutionParam& subparam,
                                 Id seed                        = 43242,
                                 Id verbose                     = false,
                                 const NamingConvention& namconv = NamingConvention("SimSub"));
} // namespace gstlrn