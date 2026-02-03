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
#include "gstlearn_export.hpp"


namespace gstlrn
{

class CovAniso;
class MultiGridSolver;
class DbGrid;
class AMesh;

class GSTLEARN_EXPORT MultiGridSPDE
{

public:
  MultiGridSPDE(const CovAniso* cov = nullptr);

  void prepare(MultiGridSolver* solver,const DbGrid* grid);
  Id buildGridHierarchy(const DbGrid* dbfine, Id nlevels = 2);
  static std::pair<DbGrid, bool> buildNextGrid(const DbGrid* dbfine);
  MatrixSparse buildProlongator(const DbGrid* dbfine, const DbGrid* dbcoarse, Id n_rings = 1);

private:
  const CovAniso* _cova;

};

GSTLEARN_EXPORT MultiGridSolver* createMultiGridSolverSPDE(const CovAniso* cov, DbGrid* grid);

} // namespace gstlrn