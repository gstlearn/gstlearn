/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

// WARNING: Make this include list as small as possible!
#include "gstlearn_export.hpp"
#include "geoslib_d.h"

#include "Basic/CSVformat.hpp"
#include "Db/ELoadBy.hpp"
#include "Db/DbGrid.hpp"
#include "Model/EConsElem.hpp"
#include "Model/Constraints.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Neigh/NeighImage.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Stats/Selectivity.hpp"

#include "Basic/NamingConvention.hpp"

class Db;
class Vario;
class VarioParam;
class Model;
class AAnam;
class ANeighParam;
class Polygons;
class RuleProp;
class ECalcVario;
class PCA;

/*************************/
/* Functions for License */
/*************************/

GSTLEARN_EXPORT int setup_license(const char *target_name);

/***********************/
/* Functions for Basic */
/***********************/

GSTLEARN_EXPORT VectorInt util_string_search(const VectorString &list_strings,
                                             const String &string,
                                             int verbose);

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

GSTLEARN_EXPORT Db* db_create_point(int nech,
                                    int ncol = 0,
                                    const ELoadBy &order = ELoadBy::COLUMN,
                                    int flag_add_rank = 0,
                                    const VectorDouble &tab = VectorDouble());
GSTLEARN_EXPORT DbGrid* db_create_grid_generic(int ndim,
                                               int ncol,
                                               const ELoadBy &order,
                                               int flag_add_rank,
                                               const VectorInt &nx,
                                               const VectorDouble &tab = VectorDouble());
GSTLEARN_EXPORT DbGrid* db_create_grid(int flag_g_rot,
                                       int ndim,
                                       int nvar,
                                       const ELoadBy &order,
                                       int flag_add_rank,
                                       const VectorInt &nx,
                                       const VectorDouble &x0,
                                       const VectorDouble &dx,
                                       const VectorDouble &angles = VectorDouble(),
                                       const VectorDouble &tab = VectorDouble());
GSTLEARN_EXPORT DbGrid* db_create_grid_2D(int flag_rot,
                                          int ncol,
                                          const ELoadBy &order,
                                          int flag_add_rank,
                                          int nx,
                                          int ny,
                                          double x0 = 0.,
                                          double y0 = 0.,
                                          double dx = 1.,
                                          double dy = 1.,
                                          double angle = 0.,
                                          const VectorDouble &tab = VectorDouble());
GSTLEARN_EXPORT DbGrid* db_create_grid_3D(int flag_rot,
                                          int ncol,
                                          const ELoadBy &order,
                                          int flag_add_rank,
                                          int nx,
                                          int ny,
                                          int nz,
                                          double x0 = 0.,
                                          double y0 = 0.,
                                          double z0 = 0.,
                                          double dx = 1.,
                                          double dy = 1.,
                                          double dz = 1.,
                                          double angle_z = 0.,
                                          double angle_y = 0.,
                                          double angle_x = 0.,
                                          const VectorDouble &tab = VectorDouble());
GSTLEARN_EXPORT VectorDouble db_get_grid_axis(DbGrid *dbgrid, int idim);
GSTLEARN_EXPORT VectorDouble db_get_attribute(Db *db,
                                              int iatt,
                                              bool verbose = false);
GSTLEARN_EXPORT VectorInt db_identify_variables_by_name(Db *db,
                                                        const String &pattern);
GSTLEARN_EXPORT void db_print(Db *db,
                              int flag_resume = 1,
                              int flag_vars = 1,
                              int flag_extend = 0,
                              int flag_stats = 0,
                              int flag_array = 0,
                              int nrank = 0,
                              int *ranks = NULL);
GSTLEARN_EXPORT void db_stats_print(const Db *db,
                                    const VectorInt &iatts = VectorInt(),
                                    const VectorString &opers = VectorString(),
                                    int flag_iso = 0,
                                    int flag_correl = 0,
                                    const String &title = String(),
                                    const String &radix = String());
GSTLEARN_EXPORT void db_stats_print(const Db *db,
                                    const VectorString &names,
                                    const VectorString &opers = VectorString(),
                                    int flag_iso = 0,
                                    int flag_correl = 0,
                                    const String &title = String(),
                                    const String &radix = String());

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
                                            const VectorDouble &codir,
                                            const VectorInt &grincr);
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
GSTLEARN_EXPORT void variogram_print(const Vario *vario);
GSTLEARN_EXPORT Vario* variogram_pgs(Db *db,
                                     const VarioParam *varioparam,
                                     const RuleProp *ruleprop,
                                     int flag_rho = false,
                                     int opt_correl = 2);
GSTLEARN_EXPORT int vmap_compute(Db *db,
                                 DbGrid *dbmap,
                                 const ECalcVario &calcul_type, // = ECalcVario::UNDEFINED,
                                 int radius = 0,
                                 bool flag_FFT = true,
                                 const NamingConvention& namconv = NamingConvention("VMAP"));
GSTLEARN_EXPORT DbGrid* db_vmap_compute(Db *db, const ECalcVario &calcul_type, // = ECalcVario::UNDEFINED,
                                        int nxx = 20,
                                        int nyy = 20,
                                        double dx = TEST,
                                        double dy = TEST,
                                        int radius = 0.,
                                        bool flag_FFT = true,
                                        const NamingConvention& namconv = NamingConvention("VMAP"));
GSTLEARN_EXPORT int pca_z2f(Db *db, const PCA *pca, bool flag_norm, bool flag_verbose,
                            const NamingConvention& namconv = NamingConvention("Z2F"));
GSTLEARN_EXPORT int pca_f2z(Db *db, const PCA *pca, bool flag_norm, bool flag_verbose,
                            const NamingConvention& namconv = NamingConvention("F2Z"));

/***********************/
/* Functions for Model */
/***********************/

GSTLEARN_EXPORT int model_auto_fit(const Vario *vario,
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

/******************************/
/* Functions for Anamorphosis */
/******************************/

GSTLEARN_EXPORT Selectivity anam_selectivity(AAnam *anam,
                                             int nclass,
                                             VectorDouble zcut,
                                             int flag_correct = 0,
                                             int verbose = 0);
GSTLEARN_EXPORT void selectivity_interpolate(int verbose,
                                             double *zcut,
                                             Selectivity& calest,
                                             Selectivity& calcut);
GSTLEARN_EXPORT int anamZToFactor(AAnam *anam,
                                 Db *db,
                                 const VectorInt &ifacs,
                                 const NamingConvention& namconv = NamingConvention("Factor"));
GSTLEARN_EXPORT int anamFactor2QT(Db *db,
                                  AAnam *anam,
                                  const VectorDouble& cutmine,
                                  double z_max,
                                  int flag_correct,
                                  const VectorInt& codes,
                                  const VectorInt& cols_est,
                                  const VectorInt& cols_std,
                                  VectorInt& qt_vars,
                                  bool verbose = false);

/******************************/
/* Functions for Neighborhood */
/******************************/

GSTLEARN_EXPORT int test_neigh(Db *dbin,
                               Db *dbout,
                               Model *model,
                               ANeighParam *neighparam,
                               const NamingConvention& namconv = NamingConvention("Neigh"));

/**********************/
/* Functions for SPDE */
/**********************/
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

/**********************************/
/* High-level Interface Functions */
/**********************************/

GSTLEARN_EXPORT int migrateByAttribute(Db *db1,
                                       Db *db2,
                                       const VectorInt &iatts = VectorInt(),
                                       int ldmax = 0,
                                       const VectorDouble &dmax = VectorDouble(),
                                       int flag_fill = false,
                                       int flag_inter = false,
                                       const NamingConvention& namconv = NamingConvention("Migrate"));
GSTLEARN_EXPORT int migrate(Db *db1,
                            Db *db2,
                            const String &name,
                            int ldmax = 0,
                            const VectorDouble &dmax = VectorDouble(),
                            int flag_fill = 0,
                            int flag_inter = 0,
                            const NamingConvention& namconv = NamingConvention("Migrate"));
GSTLEARN_EXPORT int migrateByLocator(Db *db1,
                                     Db *db2,
                                     const ELoc &locatorType,
                                     int ldmax = 0,
                                     const VectorDouble &dmax = VectorDouble(),
                                     int flag_fill = false,
                                     int flag_inter = false,
                                     const NamingConvention& namconv = NamingConvention("Migrate"));
GSTLEARN_EXPORT int db_selhull(Db *db1,
                               Db *db2,
                               bool verbose = false,
                               const NamingConvention& namconv = NamingConvention("Hull", ELoc::SEL));
GSTLEARN_EXPORT void db_polygon(Db *db,
                                Polygons *polygon,
                                int flag_sel = 0,
                                int flag_period = 0,
                                int flag_nested = 0,
                                const NamingConvention& namconv = NamingConvention("Polygon", ELoc::SEL));
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
                                 const NamingConvention& namconv = NamingConvention("Duplicate", ELoc::SEL));
GSTLEARN_EXPORT int kriging(Db *dbin,
                            Db *dbout,
                            Model *model,
                            ANeighParam *neighparam,
                            const EKrigOpt &calcul = EKrigOpt::PONCTUAL,
                            int flag_est = 1,
                            int flag_std = 1,
                            int flag_varz = 0,
                            VectorInt ndisc = VectorInt(),
                            VectorInt rank_colcok = VectorInt(),
                            VectorVectorDouble matCL = VectorVectorDouble(),
                            const NamingConvention& namconv = NamingConvention("Kriging"));
GSTLEARN_EXPORT int krigprof(Db *dbin,
                             Db *dbout,
                             Model *model,
                             ANeighParam *neighparam,
                             int flag_est = 1,
                             int flag_std = 1,
                             const NamingConvention& namconv = NamingConvention("KrigProf"));
GSTLEARN_EXPORT int krigdgm(Db *dbin,
                            Db *dbout,
                            Model *model,
                            ANeighParam *neighparam,
                            int flag_est = 1,
                            int flag_std = 1,
                            int flag_varz = 0,
                            double rval = 1.,
                            const NamingConvention& namconv = NamingConvention("DGM"));
GSTLEARN_EXPORT int krigcell(Db *dbin,
                             Db *dbout,
                             Model *model,
                             ANeighParam *neighparam,
                             int flag_est,
                             int flag_std,
                             VectorInt ndisc,
                             VectorInt rank_colcok = VectorInt(),
                             const NamingConvention& namconv = NamingConvention("KrigCell"));
GSTLEARN_EXPORT int krimage(DbGrid *dbgrid,
                            Model *model,
                            NeighImage *neighparam,
                            const NamingConvention& namconv = NamingConvention("Filtering"));
GSTLEARN_EXPORT int image_smoother(DbGrid *dbgrid,
                                   NeighImage *neighI,
                                   int type,
                                   double range,
                                   const NamingConvention& namconv = NamingConvention("Smoothing"));
GSTLEARN_EXPORT int xvalid(Db *db,
                           Model *model,
                           ANeighParam *neighparam,
                           int flag_code = 0,
                           int flag_est = 1,
                           int flag_std = 1,
                           int flag_varz = 0,
                           VectorInt rank_colcok = VectorInt(),
                           const NamingConvention& namconv = NamingConvention("Xvalid"));
GSTLEARN_EXPORT Krigtest_Res krigtest(Db *dbin,
                                      Db *dbout,
                                      Model *model,
                                      ANeighParam *neighparam,
                                      int iech0,
                                      const EKrigOpt &calcul = EKrigOpt::PONCTUAL,
                                      VectorInt ndisc = VectorInt());
GSTLEARN_EXPORT int kribayes(Db *dbin,
                             Db *dbout,
                             Model *model,
                             ANeighParam *neighparam,
                             const VectorDouble& dmean = VectorDouble(),
                             const VectorDouble& dcov = VectorDouble(),
                             int flag_est = 1,
                             int flag_std = 1,
                             const NamingConvention& namconv = NamingConvention("Bayes"));
GSTLEARN_EXPORT Global_Res global_kriging(Db *dbin,
                                          Db *dbout,
                                          Model *model,
                                          int ivar0 = 0,
                                          bool flag_verbose = false);
GSTLEARN_EXPORT Global_Res global_arithmetic(Db *dbin,
                                             DbGrid *dbgrid,
                                             Model *model,
                                             int ivar0 = 0,
                                             bool flag_verbose = false);
GSTLEARN_EXPORT int krigsum(Db *dbin,
                            Db *dbout,
                            Model *model,
                            NeighUnique* neighU,
                            bool flag_positive = false,
                            const NamingConvention& namconv = NamingConvention("KrigSum"));
GSTLEARN_EXPORT int kriggam(Db *dbin,
                            Db *dbout,
                            Model *model,
                            ANeighParam *neighparam,
                            AAnam *anam,
                            const NamingConvention& namconv = NamingConvention("KrigGam"));
GSTLEARN_EXPORT int declustering(Db *db,
                                 Model *model,
                                 int method,
                                 ANeighParam *neighparam = nullptr,
                                 DbGrid *dbgrid = nullptr,
                                 const VectorDouble& radius = VectorDouble(),
                                 const VectorInt& ndisc = VectorInt(),
                                 int flag_sel = false,
                                 bool verbose = false);
GSTLEARN_EXPORT int dk(Db* dbin,
                       DbGrid* dbsmu,
                       Model* model,
                       ANeighParam* neighparam,
                       AAnam* anam,
                       int nfactor,
                       const EKrigOpt &calcul,
                       const VectorInt &ndisc,
                       int flag_est = 1,
                       int flag_std = 1,
                       const NamingConvention& namconv = NamingConvention("KD"));
GSTLEARN_EXPORT int uc(Db *db,
                       AAnam *anam,
                       int att_est,
                       int att_var,
                       VectorDouble& cutmine,
                       double proba,
                       double var_bloc,
                       const VectorInt& codes,
                       int verbose,
                       VectorInt& qt_vars);
GSTLEARN_EXPORT int ce(Db *db,
                       AAnam *anam,
                       int att_est,
                       int att_std,
                       int flag_est,
                       int flag_std,
                       int flag_OK,
                       const VectorDouble& cutmine,
                       double proba,
                       const VectorInt& codes,
                       int nbsimu,
                       int verbose,
                       VectorInt& qt_vars);
GSTLEARN_EXPORT int simtub(Db *dbin,
                           Db *dbout,
                           Model *model,
                           ANeighParam *neighparam = nullptr,
                           int nbsimu = 1,
                           int seed = 43431,
                           int nbtuba = 100,
                           int flag_check = 0,
                           const NamingConvention& namconv = NamingConvention("Simu"));
GSTLEARN_EXPORT int simbayes(Db *dbin,
                             Db *dbout,
                             Model *model,
                             ANeighParam *neighparam,
                             int nbsimu = 1,
                             int seed = 132141,
                             const VectorDouble& dmean = VectorDouble(),
                             const VectorDouble& dcov = VectorDouble(),
                             int nbtuba = 100,
                             int flag_check = false,
                             const NamingConvention& namconv = NamingConvention("SimBayes"));
GSTLEARN_EXPORT int simpgs(Db *dbin,
                           Db *dbout,
                           RuleProp *ruleprop,
                           Model *model1,
                           Model *model2,
                           ANeighParam *neighparam = nullptr,
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
                           const NamingConvention& namconv = NamingConvention("Facies", ELoc::FACIES));
GSTLEARN_EXPORT int simbipgs(Db *dbin,
                             Db *dbout,
                             RuleProp *ruleprop,
                             Model *model11,
                             Model *model12,
                             Model *model21,
                             Model *model22,
                             ANeighParam *neighparam = nullptr,
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
                             const NamingConvention& namconv = NamingConvention("Facies", ELoc::FACIES));
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
                                           const NamingConvention& namconv = NamingConvention("Prop", ELoc::P));
GSTLEARN_EXPORT int gibbs_sampler(Db *db,
                                  Model *model,
                                  ANeighParam *neighparam,
                                  int nbsimu,
                                  int seed,
                                  int nboot,
                                  int niter,
                                  bool flag_norm,
                                  bool flag_multi_mono,
                                  bool flag_propagation,
                                  bool flag_sym_neigh,
                                  int gibbs_optstats,
                                  double percent,
                                  bool flag_ce,
                                  bool flag_cstd,
                                  bool verbose,
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

/***************************************/
/* Prototyping the functions in poly.c */
/***************************************/

GSTLEARN_EXPORT int polygon_inside(double xx,
                                   double yy,
                                   double zz,
                                   int flag_nested,
                                   Polygons *polygon);

