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
#include "Basic/AStringable.hpp"
#include <nlopt.h>

Optim::Optim(opt_algorithm algo, int dim)
  : _opt(nlopt_create((nlopt_algorithm)algo, dim))
  , _authorizedAnalyticalGradients(false)
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

void Optim::setGradient(std::function<void(const std::vector<double>&, vect)> gradient,
                        const std::vector<size_t>& dispatch,
                        const std::vector<size_t>& dispatchIndex)
{
  if (dispatch.empty())
  {
    _gradient = std::make_shared<std::function<void(const std::vector<double>&, vect)>>(
    std::move(gradient));
  }
  else 
  {

    _gradient = std::make_shared<std::function<void(const std::vector<double>&, vect)>>(
        [gradient, dispatch,dispatchIndex, this](const std::vector<double>& x, vect grad_reduced)
        {
            // Gradient complet, taille = dispatch.size() (i.e. taille initiale)
            this->_gradBuffer.resize(dispatch.size());
            gradient(x, this->_gradBuffer);

            // Initialiser le gradient réduit
            std::fill(grad_reduced.begin(), grad_reduced.end(), 0.0);

            // Agréger les contributions selon dispatch
            for (size_t i = 0; i < dispatch.size(); ++i)
               grad_reduced[dispatch[i]] += this->_gradBuffer[i];
        }
    );
  }
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

double Optim::minimize(std::vector<double>& x)
{
  double minf;
  nlopt_result res = nlopt_optimize(_opt, x.data(), &minf);
  if (res < 0) message("Warning, optimization return code is %d\n", res);
  return minf;
}

double Optim::callback(unsigned n, const double* x, double* grad, void* f_data) 
{
  auto* that = static_cast<Optim*>(f_data);
  bool gradAnalytic = that->_authorizedAnalyticalGradients;
  std::vector<double> xvec(x, x + n);
  double result = (*(that->_objective))(xvec);

  if (grad != nullptr)
  {
    if (that->_gradient && gradAnalytic)
    {
      vect grad_span(grad, n);
      (*(that->_gradient))(xvec, grad_span);
    }
    else if (!that->_gradientPartials.empty() && gradAnalytic)
    {
      if (that->_gradientPartials.size() != n)
        throw std::runtime_error("Incorrect number of gradient components");

      for (unsigned i = 0; i < n; ++i)
        grad[i] = that->_gradientPartials[i](xvec);
    }
    else
    {
      const double eps          = EPSILON8;
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
  return result;
}
