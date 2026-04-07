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

#include "Estimation/CalcGlobal.hpp"
#include "Estimation/CalcKriging.hpp"

#include "Basic/NamingConvention.hpp"
#include "geoslib_define.h"

#include "gstlearn_export.hpp"

namespace gstlrn
{
  class Db;
  class DbGrid;
  class ModelGeneric;
  class ANeigh;
  class Model;
  class NamingConvention;
  class AAnam;

  GSTLEARN_EXPORT Global_Result global_arithmetic(Db* dbin,
                                                  DbGrid* dbgrid,
                                                  ModelGeneric* model,
                                                  Id ivar0 = 0,
                                                  bool verbose = false);
  GSTLEARN_EXPORT Global_Result global_kriging(Db* dbin,
                                               Db* dbout,
                                               ModelGeneric* model,
                                               Id ivar0 = 0,
                                               bool verbose = false);
  GSTLEARN_EXPORT Id
    krimage(DbGrid* dbgrid,
            Model* model,
            ANeigh* neigh,
            bool flagFFT = false,
            bool verbose = false,
            Id seed = 13431,
            const NamingConvention& namconv = NamingConvention("Filtering"));
  GSTLEARN_EXPORT Id
    dbMorpho(DbGrid* dbgrid,
             const EMorpho& oper,
             double vmin = 0.,
             double vmax = 1.5,
             Id option = 0,
             const VectorInt& radius = VectorInt(),
             bool flagDistErode = false,
             bool verbose = false,
             const NamingConvention& namconv = NamingConvention("Morpho"));
  GSTLEARN_EXPORT Id
    dbSmoother(DbGrid* dbgrid,
               ANeigh* neigh,
               Id type = 1,
               double range = 1.,
               const NamingConvention& namconv = NamingConvention("Smooth"));

  GSTLEARN_EXPORT Id
    kriging(Db* dbin,
            Db* dbout,
            ModelGeneric* model,
            ANeigh* neigh = nullptr,
            bool flag_est = true,
            bool flag_std = true,
            bool flag_varz = false,
            const KrigOpt& krigopt = KrigOpt(),
            const NamingConvention& namconv = NamingConvention("Kriging"));
  GSTLEARN_EXPORT Id
    krigcell(Db* dbin,
             Db* dbout,
             ModelGeneric* model,
             ANeigh* neigh = nullptr,
             bool flag_est = true,
             bool flag_std = true,
             const KrigOpt& krigopt = KrigOpt(),
             const NamingConvention& namconv = NamingConvention("KrigCell"));
  GSTLEARN_EXPORT Id
    kribayes(Db* dbin,
             Db* dbout,
             ModelGeneric* model,
             ANeigh* neigh = nullptr,
             bool flag_est = true,
             bool flag_std = true,
             const NamingConvention& namconv = NamingConvention("Bayes"));
  GSTLEARN_EXPORT Id
    kriggam(Db* dbin,
            Db* dbout,
            ModelGeneric* model,
            ANeigh* neigh,
            AAnam* anam,
            const NamingConvention& namconv = NamingConvention("KrigGam"));
  GSTLEARN_EXPORT Krigtest_Res krigtest(Db* dbin,
                                        Db* dbout,
                                        ModelGeneric* model,
                                        ANeigh* neigh = nullptr,
                                        Id iech0 = 0,
                                        const KrigOpt& krigopt = KrigOpt(),
                                        bool verbose = true);
  GSTLEARN_EXPORT Id
    xvalid(Db* db,
           ModelGeneric* model,
           ANeigh* neigh = nullptr,
           bool flag_kfold = false,
           Id flag_xvalid_est = 1,
           Id flag_xvalid_std = 1,
           Id flag_xvalid_varz = 0,
           const KrigOpt& krigopt = KrigOpt(),
           const NamingConvention& namconv = NamingConvention("Xvalid"));
  GSTLEARN_EXPORT Id
    test_neigh(Db* dbin,
               Db* dbout,
               ModelGeneric* model,
               ANeigh* neigh = nullptr,
               const NamingConvention& namconv = NamingConvention("Neigh"));

  GSTLEARN_EXPORT Id
    krigingFactors(Db* dbin,
                   Db* dbout,
                   Model* model,
                   ANeigh* neigh,
                   bool flag_est = true,
                   bool flag_std = true,
                   const KrigOpt& krigopt = KrigOpt(),
                   const NamingConvention& namconv = NamingConvention("KD"));
  GSTLEARN_EXPORT Id krigingGradient(
    Db* dbin,
    Db* dbout,
    ModelGeneric* model,
    ANeigh* neigh,
    bool flag_est = true,
    bool flag_std = true,
    double ball_radius = 0.01,
    bool flagForceNumeric = false,
    const NamingConvention& namconv = NamingConvention("KrigGradient"));

  GSTLEARN_EXPORT Id inverseDistance(
    Db* dbin,
    Db* dbout,
    double exponent = 2.,
    bool flag_expand = true,
    double dmax = TEST,
    bool flag_est = true,
    bool flag_std = false,
    Model* model = nullptr,
    const NamingConvention& namconv = NamingConvention("InvDist"));
  GSTLEARN_EXPORT Id nearestNeighbor(
    Db* dbin,
    Db* dbout,
    bool flag_est = true,
    bool flag_std = false,
    Model* model = nullptr,
    const NamingConvention& namconv = NamingConvention("Nearest"));
  GSTLEARN_EXPORT Id
    movingAverage(Db* dbin,
                  Db* dbout,
                  ANeigh* neigh,
                  bool flag_est = true,
                  bool flag_std = false,
                  Model* model = nullptr,
                  const NamingConvention& namconv = NamingConvention("MovAve"));
  GSTLEARN_EXPORT Id
    movingMedian(Db* dbin,
                 Db* dbout,
                 ANeigh* neigh,
                 bool flag_est = true,
                 bool flag_std = false,
                 Model* model = nullptr,
                 const NamingConvention& namconv = NamingConvention("MovMed"));
  GSTLEARN_EXPORT Id
    leastSquares(Db* dbin,
                 Db* dbout,
                 ANeigh* neigh,
                 Id order = 0,
                 const NamingConvention& namconv = NamingConvention("LstSqr"));

} // namespace gstlrn