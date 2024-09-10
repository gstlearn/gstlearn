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
#include "gstlearn_export.hpp"
#include "geoslib_d.h"

#include "Basic/NamingConvention.hpp"
#include "Db/DbGrid.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Model/Constraints.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"
#include "Stats/Selectivity.hpp"
#include "Variogram/DirParam.hpp"

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

GSTLEARN_EXPORT Vario* variogram_pgs(Db *db,
                                     const VarioParam *varioparam,
                                     const RuleProp *ruleprop,
                                     int flag_rho = false,
                                     int opt_correl = 2);

/***********************/
/* Functions for Model */
/***********************/

GSTLEARN_EXPORT int model_auto_fit(Vario *vario,
                                   Model *model,
                                   bool verbose = false,
                                   const Option_AutoFit &mauto_arg = Option_AutoFit(),
                                   const Constraints &cons_arg = Constraints(),
                                   const Option_VarioFit &optvar_arg = Option_VarioFit());
GSTLEARN_EXPORT int vmap_auto_fit(const DbGrid* dbmap,
                                  Model* model,
                                  bool verbose                      = false,
                                  const Option_AutoFit& mauto_arg   = Option_AutoFit(),
                                  const Constraints& cons_arg       = Constraints(),
                                  const Option_VarioFit& optvar_arg = Option_VarioFit());
GSTLEARN_EXPORT int db_model_nostat(Db *db,
                                    Model *model,
                                    int icov = 0,
                                    const NamingConvention& namconv = NamingConvention("Nostat"));
GSTLEARN_EXPORT int db_cova_nostat(Db *db, ACov *cova,
                                   const NamingConvention& namconv = NamingConvention("Nostat"));
GSTLEARN_EXPORT void set_test_discrete(bool flag_discret);
GSTLEARN_EXPORT Vario* model_pgs(Db *db,
                                 const VarioParam *varioparam,
                                 const RuleProp *ruleprop,
                                 const Model *model1,
                                 const Model *model2 = nullptr);

/**********************/
/* Functions for SPDE */
/**********************/
#ifndef SWIG
GSTLEARN_EXPORT Cheb_Elem* spde_cheb_manage(int mode,
                                            int verbose,
                                            double power,
                                            const VectorDouble& blin,
                                            MatrixSparse *S,
                                            Cheb_Elem *cheb_old);
GSTLEARN_EXPORT int spde_chebychev_operate(MatrixSparse *S,
                                           Cheb_Elem *cheb_elem,
                                           const VectorDouble &lambda,
                                           const VectorDouble& x,
                                           VectorDouble& y);
#endif

/**********************************/
/* High-level Interface Functions */
/**********************************/

GSTLEARN_EXPORT int krigsum(Db *dbin,
                            Db *dbout,
                            Model *model,
                            ANeigh* neigh,
                            bool flag_positive = false,
                            const NamingConvention& namconv = NamingConvention("KrigSum"));
GSTLEARN_EXPORT int declustering(Db *db,
                                 Model *model,
                                 int method,
                                 ANeigh *neigh = nullptr,
                                 DbGrid *dbgrid = nullptr,
                                 const VectorDouble& radius = VectorDouble(),
                                 const VectorInt& ndisc = VectorInt(),
                                 int flag_sel = false,
                                 bool verbose = false);
GSTLEARN_EXPORT int simpgs(Db *dbin,
                           Db *dbout,
                           RuleProp *ruleprop,
                           Model *model1,
                           Model *model2 = nullptr,
                           ANeigh *neigh = nullptr,
                           int nbsimu = 1,
                           int seed = 1321421,
                           int flag_gaus = false,
                           int flag_prop = false,
                           int flag_check = false,
                           int flag_show = false,
                           int nbtuba = 100,
                           int gibbs_nburn = 10,
                           int gibbs_niter = 100,
                           double percent = 5.,
                           const NamingConvention& namconv = NamingConvention("Facies", true, true, true,
                                                                              ELoc::fromKey("FACIES")));
GSTLEARN_EXPORT int simbipgs(Db *dbin,
                             Db *dbout,
                             RuleProp *ruleprop,
                             Model *model11,
                             Model *model12 = nullptr,
                             Model *model21 = nullptr,
                             Model *model22 = nullptr,
                             ANeigh *neigh = nullptr,
                             int nbsimu = 1,
                             int seed = 43243,
                             int flag_gaus = false,
                             int flag_prop = false,
                             int flag_check = false,
                             int flag_show = false,
                             int nbtuba = 100,
                             int gibbs_nburn = 10,
                             int gibbs_niter = 100,
                             double percent = 5.,
                             const NamingConvention& namconv = NamingConvention("Facies", true, true, true,
                                                                                ELoc::fromKey("FACIES")));
GSTLEARN_EXPORT VectorDouble simsph_mesh(MeshSpherical *mesh,
                                         Model *model,
                                         const SimuSphericalParam& sphepar,
                                         int seed = 54523,
                                         int verbose = false);
GSTLEARN_EXPORT MatrixRectangular fluid_extract(DbGrid *dbgrid,
                                                const String& name_facies,
                                                const String& name_fluid,
                                                const String& name_poro,
                                                const String& name_date,
                                                int nfacies,
                                                int nfluids,
                                                int facies0,
                                                int fluid0,
                                                int ntime,
                                                double time0,
                                                double dtime,
                                                bool verbose = false);
GSTLEARN_EXPORT int simpgs_spde(Db* dbin,
                                Db* dbout,
                                RuleProp* ruleprop,
                                Model* model1,
                                Model* model2,
                                const String& triswitch,
                                const VectorDouble& gext,
                                int flag_gaus,
                                int flag_prop,
                                int flag_check,
                                int flag_show,
                                int nfacies,
                                int seed,
                                int nbsimu,
                                int gibbs_nburn,
                                int gibbs_niter,
                                int ngibbs_int,
                                int verbose,
                                double percent);
GSTLEARN_EXPORT int db_proportion_estimate(Db *dbin,
                                           DbGrid *dbout,
                                           Model *model,
                                           int niter = 100,
                                           bool verbose = false,
                                           const NamingConvention& namconv = NamingConvention("Prop", true, true, true,
                                                                                              ELoc::fromKey("P")));
GSTLEARN_EXPORT int gibbs_sampler(Db *dbin,
                                  Model *model,
                                  int nbsimu,
                                  int seed,
                                  int gibbs_nburn,
                                  int gibbs_niter,
                                  bool flag_moving,
                                  bool flag_norm,
                                  bool flag_multi_mono,
                                  bool flag_propagation,
                                  bool flag_sym_neigh,
                                  int gibbs_optstats,
                                  double percent,
                                  bool flag_ce,
                                  bool flag_cstd,
                                  bool verbose = false,
                                  const NamingConvention& namconv = NamingConvention("Gibbs"));

/********************************************/
/* Prototyping the functions in potential.c */
/********************************************/

GSTLEARN_EXPORT int potential_kriging(Db *db,
                                      Db *dbgrd,
                                      Db *dbtgt,
                                      DbGrid *dbout,
                                      Model *model,
                                      ANeigh *neigh,
                                      double nugget_grd = 0.,
                                      double nugget_tgt = 0.,
                                      bool flag_pot = true,
                                      bool flag_grad = false,
                                      bool flag_trans = false,
                                      bool flag_save_data = false,
                                      int opt_part = 0,
                                      bool verbose = false);
GSTLEARN_EXPORT int potential_cov(Model *model,
                                  bool verbose,
                                  int type1,
                                  const VectorDouble &x10,
                                  const VectorDouble &x1p,
                                  const VectorDouble &tx1,
                                  int type2,
                                  const VectorDouble &x20,
                                  const VectorDouble &x2p,
                                  const VectorDouble &tx2,
                                  VectorDouble &covtab);
GSTLEARN_EXPORT int potential_simulate(Db *dbiso,
                                       Db *dbgrd,
                                       Db *dbtgt,
                                       DbGrid *dbout,
                                       Model *model,
                                       ANeigh *neigh,
                                       double nugget_grd = 0.,
                                       double nugget_tgt = 0.,
                                       double dist_tempere = TEST,
                                       bool flag_trans = false,
                                       int seed = 135674,
                                       int nbsimu = 1,
                                       int nbtuba = 100,
                                       bool verbose = false);
GSTLEARN_EXPORT int potential_xvalid(Db *dbiso,
                                     Db *dbgrd,
                                     Db *dbtgt,
                                     Model *model,
                                     ANeigh *neigh,
                                     double nugget_grd = 0.,
                                     double nugget_tgt = 0.,
                                     bool flag_dist_conv = false,
                                     bool verbose = false);

