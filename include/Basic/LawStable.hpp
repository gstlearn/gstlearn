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
#include "gstlearn_export.hpp"
#include <cmath>

namespace gstlrn
{
  class GSTLEARN_EXPORT LawStable
  {
  public:
    LawStable();
    LawStable(const LawStable& r);
    LawStable& operator=(const LawStable& r);
    virtual ~LawStable();

    // auxiliary functions
    static double _sinc(const double& u);
    static double _A(const double& u, const double& alpha);
    static double _BB0(const double& u, const double& alpha);

    // simulation functions
    static double law_stable_unilateral(const double alpha);
    static double law_stable_bilateral(const double alpha);
    static double law_cauchy();

    static Id getNbrTrial() { return _nbrTrials; };

    static void setNbrTrial(Id val) { _nbrTrials = val; }

    static double law_stable_unilateral_exptilt(const double alpha,
                                                const double lambda = 0.0,
                                                bool doubleRejection = true);

  private:
    static double
      _unilateral_exptilt_singleRejection(const double alpha,
                                          const double lambda = 0.0);
    static double
      _unilateral_exptilt_doubleRejection(const double alpha,
                                          const double lambda = 0.0);

  private:
    static Id _nbrTrials; // number of trials
  };

  /*
  ** Old interface
  class GSTLEARN_EXPORT LawStable
  {
  public:
    LawStable(const double& alpha, const double& lambda = 0.0, const bool& option = true);
    LawStable(const LawStable& r);
    LawStable& operator=(const LawStable& r);
    virtual ~LawStable();
    double getAlpha() const { return _alpha; };
    double getLambda() const { return _lambda; };
    bool getOption() const { return _option; };
    double getNbrTrial() const { return _nbrTrials; };

    bool setAlpha(const double& alpha = 1.0);
    bool setLambda(const double& lambda = 0.0);
    bool setOption(const bool& option = true);
    static void setSeed(const Id& seed = 0);

    double simulate(bool isUni = true);

    static double _sinc(const double& u)
    {
      if (std::abs(u) > 0)
        return std::sin(u) / u;
      return 1.0;
    };
    static double _A(const double& u, const double& alpha)
    {
      double beta = 1 - alpha;
      double val  = 0.0;
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
    };
    static double _BB0(const double& u, const double& alpha)
    {
      double beta = 1 - alpha;
      double val  = 0.0;
      if (u != 0.0)
      {
        val = std::log(std::sin(u));
        val -= alpha * (std::log(std::sin(alpha * u)) - std::log(alpha));
        val -= beta * (std::log(std::sin(beta * u)) - std::log(beta));
      }
      return std::exp(val);
    };

    static double _simulate_cauchy()
    {
      double val = law_gaussian() / law_gaussian();
      return val;
    };

  private:
    double _simulate_uni_stable() const;
    double _simulate_bil_stable() const;
    double _simulate_exptilt_stable();
    double _sim_exptilt_stable_ini();
    double _sim_exptilt_stable_dbl();

  private:
    double _alpha;     // index of the stable distribution
    double _lambda;    // tilting coefficient
    bool _option;      // using double rejection method
    double _nbrTrials; // number of trials
  };
  */

} // namespace gstlrn
