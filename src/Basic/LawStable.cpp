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

#include "Basic/LawStable.hpp"
#include "Basic/Law.hpp"
#include "geoslib_define.h"
#include <cmath>

namespace gstlrn
{

  Id LawStable::_nbrTrials = 0;

  /**
   * Constructors and destructors for the stable distribution generators
   */

  LawStable::LawStable() {}

  LawStable::LawStable(const LawStable& r){
    DECLARE_UNUSED(r)} LawStable& LawStable::operator=(const LawStable& r)
  {
    DECLARE_UNUSED(r)
    return *this;
  }

  LawStable::~LawStable() {}

  double LawStable::_sinc(const double& u)
  {
    if (std::abs(u) > 0) return std::sin(u) / u;
    return 1.0;
  }

  double LawStable::_A(const double& u, const double& alpha)
  {
    double beta = 1 - alpha;
    double val = 0.0;
    if (u == 0.0)
    {
      val = alpha * std::log(alpha) + beta * std::log(beta);
    }
    else
    {
      val = alpha * std::log(std::sin(alpha * u));
      val += beta * std::log(std::sin(beta * u));
      val -= std::log(std::sin(u));
    }
    return std::exp(val / beta);
  }

  double LawStable::_BB0(const double& u, const double& alpha)
  {
    double beta = 1 - alpha;
    double val = 0.0;
    if (u != 0.0)
    {
      val = std::log(std::sin(u));
      val -= alpha * (std::log(std::sin(alpha * u)) - std::log(alpha));
      val -= beta * (std::log(std::sin(beta * u)) - std::log(beta));
    }
    return std::exp(val);
  }

  // functions for simulation

  double LawStable::law_stable_unilateral(const double alpha)
  {
    double val = 0.0;
    if (alpha == 1.0)
    {
      val = 1.0;
    }
    else
    {
      double U = GV_PI * law_uniform();
      double E = -std::log(law_uniform());
      double beta = (1 - alpha) / alpha;
      val = std::pow(LawStable::_A(U, alpha) / E, beta);
    }
    return val;
  }

  double LawStable::law_stable_bilateral(const double alpha)
  {
    double val = 0.0;
    if (alpha == 1.0)
    {
      val = LawStable::law_cauchy();
    }
    else
    {
      double U =
        GV_PI * law_uniform(-0.5, 0.5); // uniform distribution on [-pi/2, pi/2]
      double E = -std::log(law_uniform()); // exponential distribution
      double beta = (1 - alpha);
      val = std::pow(std::cos(beta * U) / E, beta) / std::cos(U);
      val = std::pow(val, 1 / alpha);
      val *= std::sin(alpha * U);
    }
    return val;
  }

  double LawStable::law_cauchy()
  {
    double val = law_gaussian() / law_gaussian();
    return val;
  }

  double LawStable::law_stable_unilateral_exptilt(
    const double alpha,
    const double lambda,
    bool doubleRejection)
  {
    double val = 0.0;
    if (doubleRejection)
    {
      val = _unilateral_exptilt_doubleRejection(alpha, lambda);
    }
    else
    {
      val = _unilateral_exptilt_singleRejection(alpha, lambda);
    }
    return val;
  }

  double LawStable::_unilateral_exptilt_singleRejection(
    const double alpha,
    const double lambda)
  {
    double val = 0.0;
    Id nbr = 0.0;
    do
    {
      nbr++;
      val = LawStable::law_stable_unilateral(alpha);
    } while (law_uniform() >= exp(-lambda * val));
    LawStable::setNbrTrial(nbr);
    return val;
  }

  double LawStable::_unilateral_exptilt_doubleRejection(
    const double alpha,
    const double lambda)
  {
    /*
    ** ------------------------------------------------------------
    ** set-up
    ** ------------------------------------------------------------
    */
    double beta = 1 - alpha;
    double b = beta / alpha;
    double gamma = std::pow(lambda, alpha) * alpha * beta;
    double xi = ((2 + std::sqrt(GV_PI / 2)) * std::sqrt(2 * gamma) + 1) / GV_PI;
    double psi = std::exp(-0.125 * gamma * GV_PI * GV_PI)
               * (2.0 + sqrt(GV_PI / 2.0)) * std::sqrt(gamma * GV_PI) / GV_PI;
    double w1 = xi * std::sqrt(GV_PI / (2.0 * gamma));
    double w2 = 2.0 * psi * std::sqrt(GV_PI);
    double w3 = xi * GV_PI;
    if (gamma >= 1.0)
    {
      w3 = 0.0;
    }
    else
    {
      w1 = 0.0;
    }
    double w_tot = w1 + w2 + w3;

    /*
    ** -----------------------------------------------
    ** generate U with density proportional to gG*
    ** -----------------------------------------------
    */
    Id nbr = 0.0;
    double U = 0.0;
    double Z = 0.0;
    double z = 0.0;
    do
    {
      Z = 0.0;
      z = 0.0;
      nbr++;
      // generate a proposal according to g**
      double W = law_uniform();
      if (W < w1 / w_tot)
      {
        U = law_gaussian();
        U = std::abs(U) / std::sqrt(gamma);
      }
      else if (W < (w1 + w2) / w_tot)
      {
        U = law_uniform();
        U = GV_PI * (1 - U * U);
      }
      else
      {
        U = law_uniform(0.0, GV_PI);
      }

      if (U < GV_PI)
      {
        double zeta = std::sqrt(_BB0(U, alpha));
        double phi = std::pow(std::sqrt(gamma) + alpha * zeta, 1 / alpha);
        z = phi / (phi - std::pow(gamma, 1 / (2 * alpha)));
        double g1 =
          std::exp(std::pow(lambda, alpha) * (1 - 1 / std::pow(zeta, 2)));
        g1 *= ((1 + sqrt(GV_PI / 2)) * std::sqrt(gamma) / zeta + z) / GV_PI;
        double g2 = psi / sqrt(GV_PI - U);
        if (gamma >= 1.0)
        {
          g2 = g2 + xi * std::exp(-0.5 * gamma * U * U);
        }
        else
        {
          g2 = g2 + xi;
        }
        Z = law_uniform() * g2 / g1;
      }
    } while (Z > 1.0);
    //  U has density g G*and Z is uniform on[0, 1]

    /*
    ** -----------------------------------------------
    ** generate X with density proportional to g(x, U)
    ** -----------------------------------------------
    */
    double a = _A(U, alpha);
    double m = std::pow(b * lambda / a, alpha);
    double delta = std::sqrt(m * alpha / a);
    double a1 = delta * std::sqrt(GV_PI / 2.0);
    double a2 = delta;
    double a3 = z / a;
    double s = a1 + a2 + a3;
    double X = 0.0;
    double val_test = 0.0;
    do
    {
      nbr++;
      X = 0.0;
      val_test = 0.0;
      double V = law_uniform();
      if (V < a1 / s)
      {
        double N_prime = law_gaussian();
        X = m - delta * std::abs(N_prime);
        val_test = -0.5 * N_prime * N_prime;
      }
      else if (V < (a1 + a2) / s)
      {
        double U_prime = law_uniform();
        X = m + delta * U_prime;
        val_test = 0.0;
      }
      else
      {
        double E_prime = -std::log(law_uniform());
        X = m + delta + a3 * E_prime;
        val_test = -E_prime;
      }
      val_test += a * (X - m) + lambda * (std::pow(X, -b) - std::pow(m, -b));
    } while ((X <= 0) | (val_test > -log(Z)));

    LawStable::setNbrTrial(nbr);
    return std::pow(X, -b);
  }

} // namespace gstlrn
