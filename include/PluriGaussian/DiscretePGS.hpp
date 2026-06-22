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

#include "Basic/VectorNumT.hpp"

namespace gstlrn
{
  class GSTLEARN_EXPORT DiscretePGS
  {
  public:
    DiscretePGS(bool flag_cumul, Id nconf, Id ndisc, double cmin, double cmax);
    DiscretePGS(const DiscretePGS& r) = default;
    DiscretePGS& operator=(const DiscretePGS& r) = default;
    virtual ~DiscretePGS() = default;

    void printout(Id flagPrint) const;
    double calculate(Id iconf0, const double lows[2], const double ups[2]);
    double
      calculateByRank(Id iconf0, const double rklows[2], const double rkups[2])
        const;
    Id getRankFromProba(double gaussian);
    Id getCovrank(double cova, double cround[2]) const;

  private:
    Id _getNelem() const;
    double _getCoval(Id iconf) const;
    void _getRanks(double low, double up, Id* indmin, Id* indmax);
    void _manage(Id mode, Id rank, Id* nb_used, Id* nb_max) const;
    double _integrate2D(Id iconf0, Id idisc0, Id jdisc0) const;
    double _integrate3D(Id iconf0, Id idisc0, Id jdisc0, Id kdisc0) const;

  private:
    Id _nconf; // Number of covariance configurations
    Id _ndisc; // Number of discretization steps
    Id _flagCumul; // 1 if storing integer from -infinity to value
    // 0 if storing the value per discretized class
    double _cmin; // Minimum correlation value
    double _cmax; // Maximum correlation value
    double _dc; // Covariance class interval
    double _dp; // Probability quantum for discretization
    VectorDouble _v; // Vector of thresholds (Dim: ndisc+1)
    mutable VectorVectorDouble _res; // Dimension: [nconf][size]
  };
} // namespace gstlrn
