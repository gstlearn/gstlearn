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

#include "LinearOp/MultiGridSPDE.hpp"
#include "Basic/OptCustom.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/RankHandler.hpp"
#include "Enum/ECalcMember.hpp"
#include "Enum/EKrigOpt.hpp"
#include "Estimation/KrigOpt.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "LinearOp/MultiGridSolver.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "LinearOp/PrecisionOpMatrix.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "LinearOp/IProj.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Matrix/NF_Triplet.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Space/SpacePoint.hpp"
#include "Space/SpaceRN.hpp"
#include "geoslib_define.h"
#include "Polynomials/APolynomial.hpp"
#include <memory>

namespace gstlrn
{

MultiGridSPDE::MultiGridSPDE(const CovAniso* cov, const DbGrid* gridfine, MultiGridSolver* solver)
  : _solver(solver)
  , _cova(std::make_unique<CovAniso>(*cov))
  , _gridfine(gridfine)
{
}

std::pair<DbGrid, bool> MultiGridSPDE::buildNextGrid(const DbGrid* dbfine)
{
  auto x0      = dbfine->getX0s();
  auto nxs     = dbfine->getNXs();
  auto dxs     = dbfine->getDXs();
  auto dim     = dbfine->getNDim();
  auto angles  = dbfine->getAngles();
  bool reduced = false;
  for (Id d = 0; d < dim; d++)
  {
    Id n_intervals = nxs[d] - 1;            // number of intervals in the fine grid
    Id nxcandidate = (n_intervals / 2) + 1; // number of nodes in the coarse grid if we keep every 2nd node
    if (nxcandidate > 5)
    {
      reduced  = true;
      nxs[d]   = n_intervals / 2 + 1; // number of intervals divided by 2 + 1 for the nodes
      double L = n_intervals * dxs[d];
      dxs[d]   = L / (nxs[d] - 1); // new grid spacing
    }
  }

  DbGrid next_db = DbGrid();
  next_db.reset(nxs, dxs, x0, angles, ELoadBy::fromKey("SAMPLE"), {}, {}, {}, false, false);
  return std::make_pair(next_db, reduced);
}
void MultiGridSPDE::buildGridHierarchy(Id nlevels_max, Id n_rings)
{
  n_rings = 1;
  const DbGrid* dbfine = _gridfine;
  std::pair<DbGrid, bool> current_coarse_pair;

  for (Id i = 0; i < nlevels_max - 1; i++)
  {
    auto next_pair = buildNextGrid(dbfine);
    
    if (!next_pair.second) break;

    auto mesh = MeshETurbo(dbfine);
    auto pop = std::make_unique<PrecisionOp>(&mesh, _cova.get(), false);
    auto prolongator = buildProlongator(dbfine, &next_pair.first, n_rings);

    _solver->addLevel(std::move(pop), std::make_unique<const ProjMatrix>(std::move(prolongator)));
    current_coarse_pair = std::move(next_pair);
    
    dbfine = &current_coarse_pair.first;
    
  }

  // 4. Last level(Coarse Solver)
  MeshETurbo meshlast(dbfine);
  auto pop = PrecisionOpMatrix(&meshlast, _cova.get());
  _solver->setCoarseSolver(*pop.getQ());
}
ProjMatrix MultiGridSPDE::buildProlongator(const DbGrid* dbfine, const DbGrid* dbcoarse, Id n_rings)
{
  MeshETurbo mesh_c(dbcoarse);
  MeshETurbo mesh_f(dbfine);

  // Part to make firstprivate
  MatrixSymmetric C;
  MatrixDense c0;
  VectorDouble weights;
  VectorDouble ones;
  VectorDouble s1;
  VectorDouble tabwork;
  RankHandler rkhref(dbcoarse);

  auto space = _cova->getSpace();
  SpacePoint pin(space), pout(space);
  //**********************************//
  auto nbthread = static_cast<I32>(OptCustom::query("ompthreads", 1));

  NF_Triplet triplet;
  std::vector<NF_Triplet> privateTriplets(nbthread);

  mesh_c.buildAdjacencyMatrix();
  _cova->optimizationPreProcessForData(dbcoarse);
  _cova->manage(dbcoarse, dbfine);
  CovAniso cova(*_cova);
  static RankHandler* rkh = nullptr;
  VectorInt indices_p1;
  VectorDouble lambdas;
  VectorInt all_parents;
#pragma omp threadprivate(rkh)
#pragma omp parallel for firstprivate(indices_p1, lambdas, all_parents, ones, pin, pout, tabwork, cova, C, c0, weights, s1) schedule(guided) num_threads(nbthread)
  for (Id i = 0; i < dbfine->getNSample(); i++)
  {
    indices_p1.clear();
    lambdas.clear();
    all_parents.clear();
    int tid = omp_get_thread_num();

    if (rkh == nullptr)
    {
      rkh = new RankHandler(dbcoarse);
      dbfine->initThread();
      dbcoarse->initThread();
      mesh_c.initThread();
      mesh_f.initThread();
    }

    for (Id idim = 0; idim < dbfine->getNDim(); idim++)
    {
      double value = dbfine->getCoordinate(i, idim);
      pin.setCoord(idim, value);
    }

    // Find the neighbors of the fine point in the coarse mesh

    mesh_c.getMeshFromCoordinates(pin.getCoords(), indices_p1, lambdas);

    for (const auto& idx: indices_p1)
    {
      auto voisins = mesh_c.getNRingsAdjacentApices(idx, n_rings);
      all_parents.insert(all_parents.end(), voisins.begin(), voisins.end());
    }

    std::sort(all_parents.begin(), all_parents.end());
    all_parents.erase(std::unique(all_parents.begin(), all_parents.end()), all_parents.end());

    // Compute the covariances

    rkh->defineSampleRanks(all_parents);

    // Most time consuming in the next line
    cova.evalCovMatOptimInPlace(C, dbcoarse, *rkh, EKrigOpt::POINT, ECalcMember::LHS, tabwork);
    c0.resize(all_parents.size(), 1);

    cova.evalCovVecRHSInPlace(c0.getViewOnColumnModify(0),
                              *rkh,
                              i,
                              EKrigOpt::POINT,
                              pin,
                              pout,
                              tabwork);

    // Algebra to compute ordinary kriging weights
    weights.resize(all_parents.size());

    CholeskyDense chol(C);
    chol.solve(c0.getViewOnColumn(0), weights);

    double sw = weights.sum();
    ones.resize(all_parents.size(), 1.);
    s1.resize(all_parents.size());
    chol.solve(ones, s1);
    double sc    = 1.0 / s1.sum();
    double ratio = sc * (1 - sw);
    s1 *= ratio;
    weights += s1;

    // Add the triplets for the prolongator matrix
    for (size_t j = 0; j < all_parents.size(); j++)
      if (std::abs(weights[j]) > 1e-5)
        privateTriplets[tid].add(i, all_parents[j], weights[j]);
  }

#pragma omp parallel num_threads(nbthread)
  {

    delete rkh;
    rkh = nullptr;
  }
  _cova->optimizationPostProcess();

  for (const auto& privateTriplet: privateTriplets)
    triplet.appendInPlace(privateTriplet);
  ProjMatrix result;
  result.resetFromTriplet(triplet);

  return result;
}

MultiGridSolver* createMultiGridSolverSPDE(const CovAniso* cov, DbGrid* grid, 
                                           Id nlevels_max, 
                                           Id n_rings, 
                                           Id niterPower, 
                                           double ratiochebmin, 
                                           double ratiochebmax, 
                                           Id smoothIter)
{
  auto* solver = new MultiGridSolver(ratiochebmin, ratiochebmax, smoothIter);
  MultiGridSPDE spde(cov, grid, solver);
  spde.buildGridHierarchy(nlevels_max, n_rings);
  solver->prepare(niterPower);
  return solver;
}

} // namespace gstlrn
