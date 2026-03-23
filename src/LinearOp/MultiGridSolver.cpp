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

#include "LinearOp/MultiGridSolver.hpp"

namespace gstlrn
{

  MultiGridSolver::MultiGridSolver() {}

  void MultiGridSolver::_vCycle(int lvl, constvect f, vect u) const
  {
    // 1. Base case : Solveur direct au niveau le plus grossier
    if (lvl == _nLevels - 1)
    {
      _coarseSolver.solve(f, u);
      return;
    }

    // 2. Pre-smoothing
    // auto cheby = Chebychev();
    // cheby.smoother(lvl, f, u);

    // 3. Calcul du résidu : r = f - Op * u
    // On utilise addToDest pour calculer le produit et soustra+++ire
    // Eigen::VectorXd r = f;
    //_operators[lvl]->addToDest(u, r); // Attention : gérer le signe selon ton implémentation de addToDest

    // 4. Restriction : r_c = P.point2mesh(r)
    // Eigen::VectorXd r_c = Eigen::VectorXd::Zero(_operators[lvl + 1]->getSize());
    //_transferOps[lvl]->point2mesh(r, r_c);

    // 5. Correction grossière récursive
    // Eigen::VectorXd delta_c = Eigen::VectorXd::Zero(r_c.size());
    //_vCycle(lvl + 1, r_c, delta_c);

    // 6. Prolongation : u = u + P.mesh2point(delta_c)
    // On utilise addMesh2point pour faire u += P * delta_c directement
    //_transferOps[lvl]->addMesh2point(delta_c, u);

    // 7. Post-smoothing
    //_applyChebyshev(lvl, f, u);
  }

} // namespace gstlrn
