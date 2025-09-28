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
#define DEBUG_OPTIM
namespace gstlrn
{
Optim::Optim(opt_algorithm algo, Id dim)
  : _opt(nlopt_create((nlopt_algorithm)algo, static_cast<unsigned int>(dim)))
  , _authorizedAnalyticalGradients(true)
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

void Optim::setGradient(std::function<void(vect)> gradient,
                        const std::vector<size_t>& dispatch,
                        const std::vector<size_t>& dispatchIndex)
{
  // If no link between parameters, use directly the provided gradient
  if (dispatch.empty())
  {
    _gradient = std::make_shared<std::function<void(vect)>>(
      std::move(gradient));
  }
  else
  {
    _gradient = std::make_shared<std::function<void(vect)>>(
      [gradient, dispatch, dispatchIndex, this](vect grad_reduced)
      {
        // Full gradient, size = dispatch.size() (i.e. initial size of x)
        this->_gradBuffer.resize(dispatch.size());
        gradient(this->_gradBuffer);

        // Initialize reduced gradient
        std::fill(grad_reduced.begin(), grad_reduced.end(), 0.0);

        // Aggregate contributions according to dispatch
        for (size_t i = 0; i < dispatch.size(); ++i)
          grad_reduced[dispatch[i]] += this->_gradBuffer[i];
      });
  }
}

// Used for debugging purposes (not used in the effective code)
void Optim::evalGrad(vect res)
{
  if (_gradient)
    (*_gradient)(res);
  else
    throw std::runtime_error("Gradient function not set");
}

void Optim::setGradientComponents(const std::vector<std::function<double(const std::vector<double>&)>>& partials)
{
  _gradientPartials = partials;
}

void Optim::setXtolRel(double tol)
{
  nlopt_set_xtol_rel(_opt, tol);
}
void Optim::setLowerBounds(const std::vector<double>& lb,
                           const std::vector<size_t>& dispatch)
{
  if (dispatch.empty())
  {
    nlopt_set_lower_bounds(_opt, lb.data());
    return;
  }

  // Dispatch non vide → construire les bornes réduites
  Id n_reduced = 0;
  for (auto j: dispatch) n_reduced = std::max(n_reduced, (Id)j);
  n_reduced += 1; // indices 0-based

  std::vector<double> lb_reduced(n_reduced, -HUGE_VAL);

  for (size_t i = 0; i < lb.size(); ++i)
  {
    size_t j      = dispatch[i];
    lb_reduced[j] = std::max(lb_reduced[j], lb[i]);
  }

  nlopt_set_lower_bounds(_opt, lb_reduced.data());
}

void Optim::setUpperBounds(const std::vector<double>& ub,
                           const std::vector<size_t>& dispatch)
{
  if (dispatch.empty())
  {
    nlopt_set_upper_bounds(_opt, ub.data());
    return;
  }

  // Dispatch non vide → construire les bornes réduites
  Id n_reduced = 0;
  for (auto j: dispatch) n_reduced = std::max(n_reduced, (Id)j);
  n_reduced += 1; // indices 0-based

  std::vector<double> ub_reduced(n_reduced, HUGE_VAL);

  for (size_t i = 0; i < ub.size(); ++i)
  {
    size_t j      = dispatch[i];
    ub_reduced[j] = std::min(ub_reduced[j], ub[i]);
  }

  nlopt_set_upper_bounds(_opt, ub_reduced.data());
}

double Optim::minimize(std::vector<double>& x)
{

  // This part checks that the initial guess is within bounds
  std::vector<double> lb(nlopt_get_dimension(_opt));
  std::vector<double> ub(nlopt_get_dimension(_opt));

  nlopt_get_lower_bounds(_opt, lb.data());
  nlopt_get_upper_bounds(_opt, ub.data());

  size_t dim = nlopt_get_dimension(_opt);

  if (x.size() != dim)
  {
    messerr("Dimension mismatch: initial guess has size %d but optimizer expects %d\n",
            x.size(), dim);
  }

  for (size_t i = 0; i < dim; i++)
  {
    if (x[i] < lb[i] || x[i] > ub[i])
    {
      messerr("Initial guess x[%d]=%lf is out of bounds [%lf, %lf]\n",
              i, x[i], lb[i], ub[i]);
    }
  }

  double minf;
  nlopt_result res = nlopt_optimize(_opt, x.data(), &minf);
  if (res < 0) message("Warning, optimization return code is %d\n", res);
  return minf;
}

double Optim::callback(unsigned n, const double* x, double* grad, void* f_data)
{
  auto* that        = static_cast<Optim*>(f_data);
  bool gradAnalytic = that->_authorizedAnalyticalGradients;
  std::vector<double> xvec(x, x + n);

  // ---- Objectif ----
  double result = (*(that->_objective))(xvec);
  if (std::isnan(result) || std::isinf(result))
  {
#ifdef DEBUG_OPTIM
    std::cerr << "[WARN] Objective returned NaN/Inf -> HUGE_VAL\n";
#endif
    return HUGE_VAL;
  }

  // ---- Gradient ----
  if (grad != nullptr)
  {
    if (that->_gradient && gradAnalytic)
    {
      vect grad_span(grad, n);
      (*(that->_gradient))(grad_span);
    }
    else if (!that->_gradientPartials.empty() && gradAnalytic)
    {
      if (that->_gradientPartials.size() != n)
        throw std::runtime_error("Incorrect number of gradient components");

      for (unsigned i = 0; i < n; ++i)
      {
        grad[i] = that->_gradientPartials[i](xvec);
      }
    }
    else
    {
      const double eps          = EPSILON8;
      std::vector<double> x_cur = xvec;
      for (unsigned i = 0; i < n; ++i)
      {
        x_cur[i] += eps;
        double f_plus = (*(that->_objective))(x_cur);
        if (std::isnan(f_plus) || std::isinf(f_plus)) f_plus = result;

        x_cur[i] -= 2 * eps;
        double f_minus = (*(that->_objective))(x_cur);
        if (std::isnan(f_minus) || std::isinf(f_minus)) f_minus = result;

        x_cur[i] += eps; // Restore original value
        grad[i] = (f_plus - f_minus) / (2 * eps);
      }
    }

    // ---- Protection gradient ----
    for (unsigned i = 0; i < n; ++i)
    {
      if (std::isnan(grad[i]) || std::isinf(grad[i]))
      {
#ifdef DEBUG_OPTIM
        std::cerr << "[WARN] Gradient[" << i << "] = NaN/Inf -> forced to 0\n";
#endif
        grad[i] = 0.0;
      }
    }
  }

  return result;
}

} // namespace gstlrn
