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

#include "Simulation/ACalcSimulation.hpp"
#include "Simulation/SimuPartitionParam.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
  class Db;
  class DbGrid;
  class Model;

  typedef struct
  {
    double valref;
    double valsim;
  } Stack;

  class GSTLEARN_EXPORT CalcSimuPartition: public ACalcSimulation
  {
  public:
    CalcSimuPartition(Id mode,
                      Id nbsimu = 0,
                      Id seed = 4324324,
                      bool verbose = false);
    CalcSimuPartition(const CalcSimuPartition& r) = delete;
    CalcSimuPartition& operator=(const CalcSimuPartition& r) = delete;
    virtual ~CalcSimuPartition();

    const SimuPartitionParam& getParparam() const { return _parparam; }

    void setParparam(const SimuPartitionParam& parparam)
    {
      _parparam = parparam;
    }

    void setMode(Id mode) { _mode = mode; }

    void setVerbose(bool verbose) { _verbose = verbose; }

  private:
    bool _check() override;
    bool _preprocess() override;
    bool _run() override;
    bool _postprocess() override;
    void _rollback() override;

    bool _voronoi();
    bool _poisson();
    static double _stackSearch(const std::vector<Stack>& stacks, double valref);

  private:
    Id _mode;
    bool _verbose;
    Id _iattOut;
    SimuPartitionParam _parparam;
    Model* _modelLocal;
  };

} // namespace gstlrn
