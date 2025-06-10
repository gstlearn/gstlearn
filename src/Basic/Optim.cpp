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

#include "geoslib_define.h"
#include "Basic/Optim.hpp"
#include <nlopt.h>
#include <vector>
#include <functional>
#include <memory>
#include <stdexcept>




Optim::Optim(opt_algorithm algo, int dim)
  : _opt(nlopt_create((nlopt_algorithm)algo, dim))
{
  if (!_opt) throw std::runtime_error("Échec de création de l'optimiseur NLopt");
}

Optim::~Optim()
{
  nlopt_destroy((nlopt_opt)_opt);
}

double Optim::callback(unsigned n, const double* x, double* grad, void* f_data)
{
  DECLARE_UNUSED(grad);
  auto* f = static_cast<std::function<double(const std::vector<double>&)>*>(f_data);
  return (*f)(std::vector<double>(x, x + n));
}
void Optim::setObjective(std::function<double(const std::vector<double>&)> objective)
{
  // Stocker la fonction dans un pointeur partagé pour que le callback y ait accès
  _objective = std::make_shared<std::function<double(const std::vector<double>&)>>(
    std::move(objective));

  // Définir le callback NLopt avec la fonction C statique
  nlopt_set_min_objective((nlopt_opt)_opt, &Optim::callback, _objective.get());
}

void Optim::setXtolRel(double tol)
{
  nlopt_set_xtol_rel((nlopt_opt)_opt, tol);
}

double Optim::optimize(std::vector<double>& x)
{
  double minf;
  nlopt_result res = nlopt_optimize((nlopt_opt)_opt, x.data(), &minf);
  if (res < 0) throw std::runtime_error("Échec de l'optimisation");
  return minf;
}

void Optim::setLowerBounds(const std::vector<double>& lb)
{
  nlopt_set_lower_bounds((nlopt_opt)_opt, lb.data());
}

void Optim::setUpperBounds(const std::vector<double>& ub)
{
  nlopt_set_upper_bounds((nlopt_opt)_opt, ub.data());
}



