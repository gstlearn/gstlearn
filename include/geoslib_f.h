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

// WARNING: Make this include list as small as possible!
#include "gstlearn_export.hpp"
#include "geoslib_d.h"

#include "Enum/EKrigOpt.hpp"
#include "Enum/ECalcVario.hpp"
#include "Enum/ELoadBy.hpp"
#include "Enum/EConsElem.hpp"

#include "Basic/CSVformat.hpp"
#include "Basic/NamingConvention.hpp"
#include "Db/DbGrid.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Model/Constraints.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"
#include "Simulation/SimuBooleanParam.hpp"
#include "Simulation/SimuPartitionParam.hpp"
#include "Simulation/SimuFFTParam.hpp"
#include "Stats/Selectivity.hpp"
#include "Variogram/DirParam.hpp"

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

/**********************************************/
/* Prototyping the functions in acknowledge.c */
/**********************************************/

GSTLEARN_EXPORT void acknowledge_gstlearn(void);

/***********************/
/* Functions for Basic */
/***********************/

GSTLEARN_EXPORT VectorDouble util_set_array_double(int ntab,
                                                   const double *rtab);
GSTLEARN_EXPORT VectorInt util_set_array_integer(int ntab, const int *itab);
GSTLEARN_EXPORT VectorString util_set_array_char(int ntab, char **names);
GSTLEARN_EXPORT std::vector<char*> util_vs_to_vs(VectorString vs);

/*********************/
/* Functions for CSV */
/*********************/

GSTLEARN_EXPORT int csv_manage(const char *filename,
                               const CSVformat& csv,
                               int mode,
                               int nitem,
                               bool flag_integer = 0,
                               bool verbose = false);
GSTLEARN_EXPORT void csv_print_double(double value);

/***************************/
/* Functions for Data Base */
/***************************/

GSTLEARN_EXPORT VectorDouble db_get_grid_axis(DbGrid *dbgrid, int idim);
GSTLEARN_EXPORT VectorDouble db_get_attribute(Db *db,
                                              int iatt,
                                              bool verbose = false);
GSTLEARN_EXPORT VectorInt db_identify_variables_by_name(Db *db,
                                                        const String &pattern);
GSTLEARN_EXPORT int db_center_point_to_grid(Db *db_point,
                                            DbGrid *db_grid,
                                            double eps_random = EPSILON6);

/***************************/
/* Functions for Variogram */
/***************************/

GSTLEARN_EXPORT int variogram_direction_add(VarioParam *varioparam,
                                            int npas,
                                            int opt_code,
                                            int idate,
                                            double dpas,
                                            double toldis,
                                            double tolang,
                                            double bench,
                                            double cylrad,
                                            double tolcode,
                                            const VectorDouble &breaks,
                                            const VectorDouble &codir);
GSTLEARN_EXPORT int variogram_cloud(const Db *db,
                                    const VarioParam *varioparam,
                                    DbGrid *dbgrid,
                                    const NamingConvention& namconv = NamingConvention("Cloud"));
GSTLEARN_EXPORT DbGrid* db_variogram_cloud(Db *db,
                                           const VarioParam *varioparam,
                                           double lagmax = TEST,
                                           double varmax = TEST,
                                           int lagnb = 100,
                                           int varnb = 100,
                                           const NamingConvention& namconv = NamingConvention("Cloud"));
GSTLEARN_EXPORT Vario* variogram_pgs(Db *db,
                                     const VarioParam *varioparam,
                                     const RuleProp *ruleprop,
                                     int flag_rho = false,
                                     int opt_correl = 2);
GSTLEARN_EXPORT int vmap_compute(Db *db,
                                 DbGrid *dbmap,
                                 const ECalcVario &calcul_type = ECalcVario::fromKey("VARIOGRAM"),
                                 int radius = 0,
                                 bool flag_FFT = true,
                                 const NamingConvention& namconv = NamingConvention("VMAP"));
GSTLEARN_EXPORT DbGrid* db_vmap_compute(Db *db,
                                        const ECalcVario &calcul_type = ECalcVario::fromKey("VARIOGRAM"),
                                        const VectorInt& nxx = VectorInt(),
                                        const VectorDouble& dxx = VectorDouble(),
                                        int radius = 0.,
                                        bool flag_FFT = true,
                                        const NamingConvention& namconv = NamingConvention("VMAP"));
GSTLEARN_EXPORT Db* db_variogram(Db *db, const VarioParam *varioparam);
GSTLEARN_EXPORT int dbgrid_model(DbGrid *dbgrid,
                                 Model *model,
                                 const NamingConvention &namconv = NamingConvention("VMAP"));

/***********************/
/* Functions for Model */
/***********************/

GSTLEARN_EXPORT int model_auto_fit(Vario *vario,
                                   Model *model,
                                   bool verbose = false,
                                   const Option_AutoFit &mauto_arg = Option_AutoFit(),
                                   const Constraints &cons_arg = Constraints(),
                                   const Option_VarioFit &optvar_arg = Option_VarioFit());
GSTLEARN_EXPORT int vmap_auto_fit(const DbGrid *dbvmap,
                                  Model *model,
                                  bool verbose = false,
                                  const Option_AutoFit &mauto_arg = Option_AutoFit(),
                                  const Constraints &cons_arg = Constraints(),
                                  const Option_VarioFit &optvar_arg = Option_VarioFit());
GSTLEARN_EXPORT int db_model_nostat(Db *db,
                                    Model *model,
                                    int icov = 0,
                                    const NamingConvention& namconv = NamingConvention("Nostat"));
GSTLEARN_EXPORT int is_model_nostat_param(Model *model, const EConsElem &type0);
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
                                            int nblin,
                                            double *blin,
                                            cs *S,
                                            Cheb_Elem *cheb_old);
GSTLEARN_EXPORT int spde_chebychev_operate(cs *S,
                                           Cheb_Elem *cheb_elem,
                                           const VectorDouble &lambda,
                                           const double *x,
                                           double *y);
#endif

/**********************************/
/* High-level Interface Functions */
/**********************************/

GSTLEARN_EXPORT int db_grid_fill(DbGrid *dbgrid,
                                 int mode = 0,
                                 int seed = 34243,
                                 int radius = 1,
                                 bool verbose = false,
                                 const NamingConvention& namconv = NamingConvention("Fill"));
GSTLEARN_EXPORT int db_grid1D_fill(DbGrid *dbgrid,
                                   int mode = 0,
                                   int seed = 34243,
                                   const NamingConvention& namconv = NamingConvention("Fill"));
GSTLEARN_EXPORT int db_duplicate(Db *db,
                                 bool verbose = false,
                                 double *dist = nullptr,
                                 int opt_code = 0,
                                 double tolcode = 0.,
                                 const NamingConvention& namconv = NamingConvention("Duplicate", true, true, true,
                                                                                    ELoc::fromKey("SEL")));
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
                           int nboot = 10,
                           int niter = 100,
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
                             int nboot = 10,
                             int niter = 100,
                             double percent = 5.,
                             const NamingConvention& namconv = NamingConvention("Facies", true, true, true,
                                                                                ELoc::fromKey("FACIES")));
GSTLEARN_EXPORT int simbool(Db *dbin,
                            DbGrid *dbout,
                            ModelBoolean *tokens,
                            const SimuBooleanParam& boolparam = SimuBooleanParam(),
                            int seed = 432431,
                            bool flag_simu = true,
                            bool flag_rank = true,
                            bool verbose = false,
                            const NamingConvention& namconv = NamingConvention("Boolean"));
GSTLEARN_EXPORT int simsph(DbGrid *db,
                           Model *model,
                           const SimuSphericalParam& sphepar,
                           int seed,
                           bool verbose,
                           const NamingConvention& namconv = NamingConvention("SimSphe"));
GSTLEARN_EXPORT VectorDouble simsph_mesh(MeshSpherical *mesh,
                                         Model *model,
                                         const SimuSphericalParam& sphepar,
                                         int seed = 54523,
                                         int verbose = false);
GSTLEARN_EXPORT DbGrid* simfine(DbGrid *dbin,
                                Model *model,
                                const SimuRefineParam& param,
                                int seed);
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
GSTLEARN_EXPORT int simpgs_spde(Db *dbin,
                                Db *dbout,
                                RuleProp *ruleprop,
                                Model *model1,
                                Model *model2,
                                const String &triswitch,
                                const VectorDouble &gext,
                                int flag_gaus,
                                int flag_modif,
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
GSTLEARN_EXPORT Db* db_read_csv(const char *filename,
                                const CSVformat& csvfmt,
                                int verbose = 0,
                                int ncol_max = -1,
                                int nrow_max = -1,
                                int flag_add_rank = 0);
GSTLEARN_EXPORT int db_write_csv(Db *db,
                                 const char *filename,
                                 const CSVformat& csv,
                                 int flag_allcol = 1,
                                 int flag_coor = 1,
                                 bool flag_integer = false);
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

/*****************/
/* Various Tools */
/*****************/

GSTLEARN_EXPORT int db_tool_duplicate(Db *db1,
                                      Db *db2,
                                      bool flag_same,
                                      bool verbose,
                                      int opt_code,
                                      double tolcode,
                                      double *dist,
                                      double *sel);

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
                                     int flag_dist_conv = false,
                                     bool verbose = false);

