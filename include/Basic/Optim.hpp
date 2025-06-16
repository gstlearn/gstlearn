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

#include "gstlearn_export.hpp"
#include "geoslib_define.h"
#include <vector>
#include <functional>
#include <memory>

typedef enum
{
  /* Naming conventions:

     NLOPT_{G/L}{D/N}_*
     = global/local derivative/no-derivative optimization,
     respectively

     *_RAND algorithms involve some randomization.

     *_NOSCAL algorithms are *not* scaled to a unit hypercube
     (i.e. they are sensitive to the units of x)
   */

  // NLOPT_GN_DIRECT = 0,
  // NLOPT_GN_DIRECT_L,
  // NLOPT_GN_DIRECT_L_RAND,
  // NLOPT_GN_DIRECT_NOSCAL,
  // NLOPT_GN_DIRECT_L_NOSCAL,
  // NLOPT_GN_DIRECT_L_RAND_NOSCAL,

  // NLOPT_GN_ORIG_DIRECT,
  // NLOPT_GN_ORIG_DIRECT_L,

  // NLOPT_GD_STOGO,
  // NLOPT_GD_STOGO_RAND,

  // NLOPT_LD_LBFGS_NOCEDAL,

  LBFGS = 11,

  // NLOPT_LN_PRAXIS,

  // NLOPT_LD_VAR1,
  // NLOPT_LD_VAR2,

  // NLOPT_LD_TNEWTON,
  // NLOPT_LD_TNEWTON_RESTART,
  // NLOPT_LD_TNEWTON_PRECOND,
  // NLOPT_LD_TNEWTON_PRECOND_RESTART,

  // NLOPT_GN_CRS2_LM,

  // NLOPT_GN_MLSL,
  // NLOPT_GD_MLSL,
  // NLOPT_GN_MLSL_LDS,
  // NLOPT_GD_MLSL_LDS,

  // NLOPT_LD_MMA,

  // NLOPT_LN_COBYLA,

  // NLOPT_LN_NEWUOA,
  // NLOPT_LN_NEWUOA_BOUND,

  NELDERMEAD = 28,
  // NLOPT_LN_SBPLX,

  // NLOPT_LN_AUGLAG,
  // NLOPT_LD_AUGLAG,
  // NLOPT_LN_AUGLAG_EQ,
  // NLOPT_LD_AUGLAG_EQ,

  // NLOPT_LN_BOBYQA,

  // NLOPT_GN_ISRES,

  // /* new variants that require local_optimizer to be set,
  //    not with older constants for backwards compatibility */
  // NLOPT_AUGLAG,
  // NLOPT_AUGLAG_EQ,
  // NLOPT_G_MLSL,
  // NLOPT_G_MLSL_LDS,

  // NLOPT_LD_SLSQP,

  // NLOPT_LD_CCSAQ,

  // NLOPT_GN_ESCH,

  // NLOPT_GN_AGS,

} opt_algorithm;

struct nlopt_opt_s;

class GSTLEARN_EXPORT Optim
{
public:
  Optim(opt_algorithm algo, int dim);
  Optim(const Optim&)            = delete;
  Optim& operator=(const Optim&) = delete;
  ~Optim();

  void setObjective(std::function<double(const std::vector<double>&)> objective);
  void setGradient(std::function<void(const std::vector<double>&, vect)> gradient);
  void setGradientComponents(const std::vector<std::function<double(const std::vector<double>&)>>& partials);

  void setXtolRel(double tol);
  void setLowerBounds(const std::vector<double>& lb);
  void setUpperBounds(const std::vector<double>& ub);
  double minimize(std::vector<double>& x);

private:
  static double callback(unsigned n, const double* x, double* grad, void* f_data);

  nlopt_opt_s* _opt;
  std::shared_ptr<std::function<double(const std::vector<double>&)>> _objective;
  std::shared_ptr<std::function<void(const std::vector<double>&, vect)>> _gradient;
  std::vector<std::function<double(const std::vector<double>&)>> _gradientPartials;
};
