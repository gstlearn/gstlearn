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

#include "Basic/Optim.hpp"
#include "geoslib_define.h"
#include <nlopt.h>

Optim::Optim(opt_algorithm algo, int dim)
  : _opt(nlopt_create((nlopt_algorithm)algo, dim))
{
  if (!_opt) throw std::runtime_error("Échec de création de l'optimiseur NLopt");
}

Optim::~Optim()
{
  nlopt_destroy(_opt);
}

void Optim::setObjective(std::function<double(const std::vector<double>&)> objective)
{
  _objective = std::make_shared<std::function<double(const std::vector<double>&)>>(
    std::move(objective));
  nlopt_set_min_objective(_opt, &Optim::callback, this);
}

void Optim::setGradient(std::function<void(const std::vector<double>&, vect)> gradient)
{
  _gradient = std::make_shared<std::function<void(const std::vector<double>&, vect)>>(
    std::move(gradient));
}

void Optim::setGradientComponents(const std::vector<std::function<double(const std::vector<double>&)>>& partials)
{
  _gradientPartials = partials;
}

void Optim::setXtolRel(double tol)
{
  nlopt_set_xtol_rel(_opt, tol);
}

void Optim::setLowerBounds(const std::vector<double>& lb)
{
  nlopt_set_lower_bounds(_opt, lb.data());
}

void Optim::setUpperBounds(const std::vector<double>& ub)
{
  nlopt_set_upper_bounds(_opt, ub.data());
}

double Optim::optimize(std::vector<double>& x)
{
  double minf;
  nlopt_result res = nlopt_optimize(_opt, x.data(), &minf);
  if (res < 0) throw std::runtime_error("Échec de l'optimisation");
  return minf;
}

double Optim::callback(unsigned n, const double* x, double* grad, void* f_data)
{
  auto* that = static_cast<Optim*>(f_data);
  std::vector<double> xvec(x, x + n);

  if (grad != nullptr)
  {
    if (that->_gradient)
    {
      std::span<double> grad_span(grad, n);
      (*(that->_gradient))(xvec, grad_span);
    }
    else if (!that->_gradientPartials.empty())
    {
      if (that->_gradientPartials.size() != n)
        throw std::runtime_error("Incorrect number of gradient components");

      for (unsigned i = 0; i < n; ++i)
        grad[i] = that->_gradientPartials[i](xvec);
    }
    else 
    {
      const double eps = EPSILON8;
      std::vector<double> x_cur = xvec;
      for (unsigned i = 0; i < n; ++i)
      {
        
        x_cur[i] += eps;
        double f_plus = (*(that->_objective))(x_cur);

        x_cur[i] -= 2 * eps;
        double f_minus = (*(that->_objective))(x_cur);
        x_cur[i] += eps; // Restore original value
        grad[i] = (f_plus - f_minus) / (2 * eps);
      
      }
    }
  }

  return (*(that->_objective))(xvec);
}



