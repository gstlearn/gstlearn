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

#include "Covariances/CovAniso.hpp"
#include "geoslib_d.h"
#include "gstlearn_export.hpp"

#include "Basic/NamingConvention.hpp"
#include "Db/DbGrid.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Model/Constraints.hpp"
#include "Model/Model.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"
#include "Stats/Selectivity.hpp"
#include "Variogram/DirParam.hpp"

namespace gstlrn
{
class CovAniso;
class Db;
class Vario;
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

/***************************/
/* Functions for Variogram */
/***************************/

GSTLEARN_EXPORT Vario* variogram_pgs(Db* db,
                                     const VarioParam* varioparam,
                                     const RuleProp* ruleprop,
                                     Id flag_rho   = false,
                                     Id opt_correl = 2);

/***********************/
/* Functions for Model */
/***********************/

GSTLEARN_EXPORT Id model_auto_fit(Vario* vario,
                                   Model* model,
                                   bool verbose                      = false,
                                   const Option_AutoFit& mauto_arg   = Option_AutoFit(),
                                   const Constraints& cons_arg       = Constraints(),
                                   const Option_VarioFit& optvar_arg = Option_VarioFit());
GSTLEARN_EXPORT Id vmap_auto_fit(const DbGrid* dbmap,
                                  Model* model,
                                  bool verbose                      = false,
                                  const Option_AutoFit& mauto_arg   = Option_AutoFit(),
                                  const Constraints& cons_arg       = Constraints(),
                                  const Option_VarioFit& optvar_arg = Option_VarioFit());
GSTLEARN_EXPORT void set_test_discrete(bool flag_discret);
GSTLEARN_EXPORT Vario* model_pgs(Db* db,
                                 const VarioParam* varioparam,
                                 const RuleProp* ruleprop,
                                 const Model* model1,
                                 const Model* model2 = nullptr);

/**********************************/
/* High-level Interface Functions */
/**********************************/

GSTLEARN_EXPORT Id krigsum(Db* dbin,
                            Db* dbout,
                            Model* model,
                            ANeigh* neigh,
                            bool flag_positive              = false,
                            const NamingConvention& namconv = NamingConvention("KrigSum"));
GSTLEARN_EXPORT Id declustering(Db* db,
                                 Model* model,
                                 Id method,
                                 ANeigh* neigh              = nullptr,
                                 DbGrid* dbgrid             = nullptr,
                                 const VectorDouble& radius = VectorDouble(),
                                 const VectorInt& ndisc     = VectorInt(),
                                 Id flag_sel               = false,
                                 bool verbose               = false);
GSTLEARN_EXPORT Id simpgs(Db* dbin,
                           Db* dbout,
                           RuleProp* ruleprop,
                           Model* model1,
                           Model* model2                   = nullptr,
                           ANeigh* neigh                   = nullptr,
                           Id nbsimu                      = 1,
                           Id seed                        = 1321421,
                           Id flag_gaus                   = false,
                           Id flag_prop                   = false,
                           Id flag_check                  = false,
                           Id flag_show                   = false,
                           Id nbtuba                      = 100,
                           Id gibbs_nburn                 = 10,
                           Id gibbs_niter                 = 100,
                           double percent                  = 5.,
                           const NamingConvention& namconv = NamingConvention("Facies", true, true, true, ELoc::fromKey("FACIES")));
GSTLEARN_EXPORT Id simbipgs(Db* dbin,
                             Db* dbout,
                             RuleProp* ruleprop,
                             Model* model11,
                             Model* model12                  = nullptr,
                             Model* model21                  = nullptr,
                             Model* model22                  = nullptr,
                             ANeigh* neigh                   = nullptr,
                             Id nbsimu                      = 1,
                             Id seed                        = 43243,
                             Id flag_gaus                   = false,
                             Id flag_prop                   = false,
                             Id flag_check                  = false,
                             Id flag_show                   = false,
                             Id nbtuba                      = 100,
                             Id gibbs_nburn                 = 10,
                             Id gibbs_niter                 = 100,
                             double percent                  = 5.,
                             const NamingConvention& namconv = NamingConvention("Facies", true, true, true, ELoc::fromKey("FACIES")));
GSTLEARN_EXPORT VectorDouble simsph_mesh(MeshSpherical* mesh,
                                         Model* model,
                                         const SimuSphericalParam& sphepar,
                                         Id seed    = 54523,
                                         Id verbose = false);
GSTLEARN_EXPORT MatrixDense fluid_extract(DbGrid* dbgrid,
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
GSTLEARN_EXPORT Id db_proportion_estimate(Db* dbin,
                                           DbGrid* dbout,
                                           Model* model,
                                           Id niter                       = 100,
                                           bool verbose                    = false,
                                           const NamingConvention& namconv = NamingConvention("Prop", true, true, true, ELoc::fromKey("P")));
GSTLEARN_EXPORT Id gibbs_sampler(Db* dbin,
                                  Model* model,
                                  Id nbsimu,
                                  Id seed,
                                  Id gibbs_nburn,
                                  Id gibbs_niter,
                                  bool flag_moving,
                                  bool flag_norm,
                                  bool flag_multi_mono,
                                  bool flag_propagation,
                                  bool flag_sym_neigh,
                                  Id gibbs_optstats,
                                  double percent,
                                  bool flag_ce,
                                  bool flag_cstd,
                                  bool verbose                    = false,
                                  const NamingConvention& namconv = NamingConvention("Gibbs"));

} // namespace gstlrn