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

#include "Db/DbGrid.hpp"
#include "LinearOp/MultiGridSolver.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"
#include <memory>

namespace gstlrn
{

  class CovAniso;
  class MultiGridSolver;
  class DbGrid;
  class AMesh;
  class PrecisionOp;
  class ProjMatrix;

  class GSTLEARN_EXPORT MultiGridSPDE
  {

  public:
    MultiGridSPDE(
      const CovAniso* cov = nullptr,
      const DbGrid* gridfine = nullptr,
      MultiGridSolver* solver = nullptr);
    void buildGridHierarchy(Id nlevels_max = 2, Id n_rings = 1);

#ifndef SWIG
    static std::pair<DbGrid, bool> buildNextGrid(const DbGrid* dbfine);
#endif

    ProjMatrix buildProlongator(
      const DbGrid* dbfine,
      const DbGrid* dbcoarse,
      Id n_rings = 1);

  private:
    MultiGridSolver* _solver;
    std::unique_ptr<CovAniso> _cova;
    const DbGrid* _gridfine;
  };

  GSTLEARN_EXPORT MultiGridSolver* createMultiGridSolverSPDE(
    const CovAniso* cov,
    DbGrid* grid,
    Id nlevels_max = 2,
    Id n_rings = 1,
    Id niterPower = 100,
    double ratiochebmin = 0.1,
    double ratiochebmax = 1.1,
    Id smoothIter = 4);

} // namespace gstlrn
