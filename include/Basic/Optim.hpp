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


#include "geoslib_define.h"
#include <nlopt.h>
#include <vector>
#include <functional>
#include <memory>
#include <stdexcept>

class Optim
{
public:
  Optim(nlopt_algorithm algo, int dim)
    : _opt(nlopt_create(algo, dim))
  {
    if (!_opt) throw std::runtime_error("Échec de création de l'optimiseur NLopt");
  }

  ~Optim()
  {
    nlopt_destroy(_opt);
  }

  void setObjective(std::function<double(const std::vector<double>&)> objective)
  {
    // Stocker la fonction dans un pointeur partagé pour que le callback y ait accès
    _objective = std::make_shared<std::function<double(const std::vector<double>&)>>(
      std::move(objective));

    // Définir le callback NLopt avec la fonction C statique
    nlopt_set_min_objective(_opt, &Optim::callback, _objective.get());
  }

  void setXtolRel(double tol)
  {
    nlopt_set_xtol_rel(_opt, tol);
  }

  double optimize(std::vector<double>& x)
  {
    double minf;
    nlopt_result res = nlopt_optimize(_opt, x.data(), &minf);
    if (res < 0) throw std::runtime_error("Échec de l'optimisation");
    return minf;
  }

  void setLowerBounds(const std::vector<double>& lb)
  {
    nlopt_set_lower_bounds(_opt, lb.data());
  }

  void setUpperBounds(const std::vector<double>& ub)
  {
    nlopt_set_upper_bounds(_opt, ub.data());
  }

private:
  nlopt_opt _opt;
  std::shared_ptr<std::function<double(const std::vector<double>&)>> _objective;

  // Callback C compatible NLopt
  static double callback(unsigned n, const double* x, double* grad, void* f_data)
  {
    DECLARE_UNUSED(grad);
    auto* f = static_cast<std::function<double(const std::vector<double>&)>*>(f_data);
    return (*f)(std::vector<double>(x, x + n));
  }
};
