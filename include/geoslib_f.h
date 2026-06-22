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

#include "geoslib_d.h"

#include "Basic/NamingConvention.hpp"
#include "Covariances/CovAniso.hpp"
#include "Model/Model.hpp"
#include "Stats/Selectivity.hpp"
#include "Variogram/DirParam.hpp"
#include "Variogram/Vario.hpp"

namespace gstlrn
{
  class CovAniso;
  class Db;
  class VarioParam;
  class Model;
  class AAnam;
  class ANeigh;
  class Polygons;
  class RuleProp;
  class PCA;
  class ModelBoolean;
  class SimuBooleanParam;
  class SimuSphericalParam;
  class MeshSpherical;
  class SimuSubstitutionParam;
  class SimuRefineParam;

  /***********************/
  /* Functions for Model */
  /***********************/

  GSTLEARN_EXPORT Id model_auto_fit(
    Vario* vario,
    Model* model,
    bool verbose = false,
    const Option_AutoFit& mauto_arg = Option_AutoFit(),
    const Constraints& cons_arg = Constraints(),
    const Option_VarioFit& optvar_arg = Option_VarioFit());
  GSTLEARN_EXPORT Id vmap_auto_fit(
    const DbGrid* dbmap,
    Model* model,
    bool verbose = false,
    const Option_AutoFit& mauto_arg = Option_AutoFit(),
    const Constraints& cons_arg = Constraints(),
    const Option_VarioFit& optvar_arg = Option_VarioFit());

  /**********************************/
  /* High-level Interface Functions */
  /**********************************/

  GSTLEARN_EXPORT Id krigsum(
    Db* dbin,
    Db* dbout,
    Model* model,
    ANeigh* neigh,
    bool flag_positive = false,
    const NamingConvention& namconv = NamingConvention("KrigSum"));
  GSTLEARN_EXPORT Id declustering(
    Db* db,
    Model* model,
    Id method,
    ANeigh* neigh = nullptr,
    DbGrid* dbgrid = nullptr,
    const VectorDouble& radius = VectorDouble(),
    const VectorInt& ndisc = VectorInt(),
    Id flag_sel = false,
    bool verbose = false);
  GSTLEARN_EXPORT VectorDouble simsph_mesh(
    MeshSpherical* mesh,
    Model* model,
    const SimuSphericalParam& sphepar,
    Id seed = 54523,
    Id verbose = false);
  GSTLEARN_EXPORT MatrixDense fluid_extract(
    DbGrid* dbgrid,
    const String& name_facies,
    const String& name_fluid,
    const String& name_poro,
    const String& name_date,
    Id nfacies,
    Id nfluids,
    Id facies0,
    Id fluid0,
    Id ntime,
    double time0,
    double dtime,
    bool verbose = false);
  GSTLEARN_EXPORT Id db_proportion_estimate(
    Db* dbin,
    DbGrid* dbout,
    Model* model,
    Id niter = 100,
    bool verbose = false,
    const NamingConvention& namconv =
      NamingConvention("Prop", true, true, true, ELoc::fromKey("P")));
  GSTLEARN_EXPORT Id gibbs_sampler(
    Db* dbin,
    Model* model,
    Id nbsimu = 1,
    Id seed = 14341,
    Id gibbs_nburn = 10,
    Id gibbs_niter = 30,
    bool flag_moving = false,
    bool flag_norm = true,
    bool flag_multi_mono = false,
    bool flag_propagation = false,
    bool flag_sym_neigh = false,
    Id gibbs_optstats = 0,
    double percent = 5.,
    bool flag_ce = true,
    bool flag_cstd = false,
    bool verbose = false,
    const NamingConvention& namconv = NamingConvention("Gibbs"));

  GSTLEARN_EXPORT double ut_trace_discretize(
    const MatrixDense& trace,
    double disc,
    VectorDouble& xp,
    VectorDouble& yp,
    VectorDouble& dd,
    VectorDouble& del);
  GSTLEARN_EXPORT void ut_trace_sample(
    Db* db,
    const VectorDouble& xp,
    const VectorDouble& yp,
    const VectorDouble& dd,
    const ELoc& ptype,
    double radius,
    VectorDouble& xs,
    VectorDouble& ys,
    VectorInt& rks,
    VectorInt& lys,
    VectorInt& typ);
} // namespace gstlrn
