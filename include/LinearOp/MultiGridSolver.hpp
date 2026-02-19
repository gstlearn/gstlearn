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

#include "LinearOp/ALinearOp.hpp"
#include "LinearOp/CholeskySparse.hpp"
#include "LinearOp/APreconditioner.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"
#include <Eigen/Core>
#include <functional>
#include <vector>

#pragma once

namespace gstlrn
{

class IProj;

class GSTLEARN_EXPORT MultiGridSolver : public APreconditioner
{

public:
  MultiGridSolver();




  // On prépare la hiérarchie
  MultiGridSolver& compute(const ALinearOp& /*default*/) { return *this; }

  // Méthode appelée par le CG d'Eigen
  template<typename Rhs>
  Rhs solve(const Rhs& b) const
  {
   // Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());
   // _vCycle(0, b, x);
    return b;
  }

  // Indispensable pour l'interface Eigen
  MultiGridSolver& analyzePattern(const ALinearOp& /*default*/) { return *this; }
  MultiGridSolver& factorize(const ALinearOp& /*default*/) { return *this; }
  

private:
  void _vCycle(int lvl, constvect f, vect u) const;

private:
  int _nLevels;
  //std::vector<const TLinOP*> _operators;  // La hiérarchie des Q
  std::vector<std::reference_wrapper<const IProj>> _transferOps; // La hiérarchie des transferts (spécialisations de IProj)
  std::vector<ALinearOp*> _operators;     // La hiérarchie des opérateurs
  VectorDouble _lMax;              // Rayons spectraux pour Chebyshev
  // Solveur direct pour le niv&eau L-1
  mutable CholeskySparse _coarseSolver;
};

} // namespace gstlrn