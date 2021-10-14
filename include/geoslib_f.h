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
#ifndef GEOSLIB_F_H
#define GEOSLIB_F_H

// WARNING: Make this include list as small as possible!
#include "geoslib_d.h"

#include "Db/ELoadBy.hpp"
#include "Model/Constraints.hpp"
#include "Basic/NamingConvention.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"

class Db;
class Vario;
class VarioParam;
class Model;
class Anam;
class Neigh;
class Polygons;
class RuleProp;
class ECalcVario;

/*************************/
/* Functions for License */
/*************************/

GEOSLIB_API int setup_license(const char *target_name);

/***********************/
/* Functions for Basic */
/***********************/

GEOSLIB_API void print_ivector(const char *title,
                                int flag_limit,
                                int ntab,
                                const int *itab);
GEOSLIB_API void print_ivector(const char *title,
                               int flag_limit,
                               int ntab,
                               const VectorInt& itab);
GEOSLIB_API void print_vector(const char *title,
                              int flag_limit,
                              int ntab,
                              const VectorDouble& tab);
GEOSLIB_API VectorInt util_string_search(const VectorString& list_strings,
                                         const String& string,
                                         int verbose);

GEOSLIB_API VectorDouble util_set_array_double(int ntab, const double *rtab);
GEOSLIB_API VectorInt util_set_array_integer(int ntab, const int *itab);
GEOSLIB_API VectorString util_set_array_char(int ntab, char **names);
GEOSLIB_API std::vector<char *> util_vs_to_vs(VectorString vs);

/*********************/
/* Functions for CSV */
/*********************/

GEOSLIB_API int csv_manage(const char *filename,
                           int mode,
                           int nitem,
                           bool flag_integer = 0,
                           const char *char_sep = ",",
                           const char *na_string = "NA",
                           bool verbose = false);
GEOSLIB_API void csv_print_string(const char *string);
GEOSLIB_API void csv_print_double(double value);
GEOSLIB_API void csv_print_eol(void);

/***************************/
/* Functions for Data Base */
/***************************/

GEOSLIB_API Db *db_create_point(int nech,
                                int ncol = 0,
                                const ELoadBy& order = ELoadBy::COLUMN,
                                int flag_add_rank = 0,
                                const VectorDouble& tab = VectorDouble());
GEOSLIB_API Db *db_create_grid_generic(int flag_rot,
                                       int ndim,
                                       int ncol,
                                       const ELoadBy& order,
                                       int flag_add_rank,
                                       const VectorInt& nx,
                                       const VectorDouble& tab = VectorDouble());
GEOSLIB_API Db *db_create_grid(int flag_g_rot,
                               int ndim,
                               int nvar,
                               const ELoadBy& order,
                               int flag_add_rank,
                               const VectorInt& nx,
                               const VectorDouble& x0,
                               const VectorDouble& dx,
                               const VectorDouble& angles = VectorDouble(),
                               const VectorDouble& tab = VectorDouble());
GEOSLIB_API Db *db_create_grid_2D(int flag_rot,
                                  int ncol,
                                  const ELoadBy& order,
                                  int flag_add_rank,
                                  int nx,
                                  int ny,
                                  double x0 = 0.,
                                  double y0 = 0.,
                                  double dx = 1.,
                                  double dy = 1.,
                                  double angle = 0.,
                                  const VectorDouble& tab = VectorDouble());
GEOSLIB_API Db *db_create_grid_3D(int flag_rot,
                                  int ncol,
                                  const ELoadBy& order,
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
                                  const VectorDouble& tab = VectorDouble());
GEOSLIB_API VectorDouble db_get_grid_axis(Db *dbgrid, int idim);
GEOSLIB_API VectorDouble db_get_attribute(Db *db, int iatt, bool verbose= false);
GEOSLIB_API VectorInt db_identify_variables_by_name(Db *db,
                                                    const String& pattern,
                                                    int verbose = 0);
GEOSLIB_API void db_print(Db *db,
                          int flag_resume = 1,
                          int flag_vars = 1,
                          int flag_extend = 0,
                          int flag_stats = 0,
                          int flag_array = 0,
                          int nrank = 0,
                          int *ranks = NULL);
GEOSLIB_API void db_stats_print(const Db *db,
                                const VectorInt& iatts = VectorInt(),
                                const VectorString& opers = VectorString(),
                                int flag_iso = 0,
                                int flag_correl = 0,
                                const String& title = String(),
                                const String& radix = String());
GEOSLIB_API void db_stats_print(const Db *db,
                                const VectorString& names,
                                const VectorString& opers = VectorString(),
                                int flag_iso = 0,
                                int flag_correl = 0,
                                const String& title = String(),
                                const String& radix = String());

/***************************/
/* Functions for Variogram */
/***************************/

GEOSLIB_API int variogram_direction_add(VarioParam *varioparam,
                                        int npas,
                                        int opt_code,
                                        int idate,
                                        double dpas,
                                        double toldis,
                                        double tolang,
                                        double bench,
                                        double cylrad,
                                        double tolcode,
                                        const VectorDouble& breaks,
                                        const VectorDouble& codir,
                                        const VectorInt&    grincr);
GEOSLIB_API int variogram_cloud(Db *db,
                                const VarioParam *varioparam,
                                Db *dbgrid,
                                NamingConvention namconv = NamingConvention("Cloud"));
GEOSLIB_API Db* db_variogram_cloud(Db *db,
                                   const VarioParam* varioparam,
                                   double lagmax = TEST,
                                   double varmax = TEST,
                                   int lagnb = 100,
                                   int varnb = 100,
                                   NamingConvention namconv = NamingConvention("Cloud"));
GEOSLIB_API void variogram_print(Vario *vario, int verbose = false);
GEOSLIB_API Vario* variogram_pgs(Db*               db,
                                 const VarioParam* varioparam,
                                 const RuleProp*   ruleprop,
                                 int               flag_rho = false,
                                 int               opt_correl = 2);
GEOSLIB_API int vmap_compute(Db*               db,
                             Db*               dbmap,
                             const ECalcVario& calcul_type,// = ECalcVario::UNDEFINED,
                             int               radius = 0,
                             bool              flag_FFT = true,
                             NamingConvention  namconv = NamingConvention("VMAP"));
GEOSLIB_API Db* db_vmap_compute(Db*                db,
                                const ECalcVario&  calcul_type,// = ECalcVario::UNDEFINED,
                                int                nxx = 20,
                                int                nyy = 20,
                                double             dx = TEST,
                                double             dy = TEST,
                                int                radius = 0.,
                                bool               flag_FFT = true,
                                NamingConvention   namconv = NamingConvention("VMAP"));

/***********************/
/* Functions for Model */
/***********************/

GEOSLIB_API Model *model_init(int ndim = 2,
                              int nvar = 1,
                              double field = 0,
                              int flag_linked = 0,
                              double ball_radius = 0.,
                              bool flag_gradient = false,
                              const VectorDouble& mean = VectorDouble(),
                              const VectorDouble& covar0 = VectorDouble());
GEOSLIB_API int model_auto_fit(Vario *vario,
                               Model *model,
                               bool verbose = false,
                               Option_AutoFit mauto = Option_AutoFit(),
                               const Constraints& consarg = Constraints(),
                               Option_VarioFit optvar = Option_VarioFit());
GEOSLIB_API int vmap_auto_fit(Db *dbvmap,
                              Model *model,
                              bool verbose = false,
                              Option_AutoFit mauto = Option_AutoFit(),
                              const Constraints& constraints = Constraints(),
                              Option_VarioFit optvar = Option_VarioFit());
GEOSLIB_API int db_model_nostat(Db *db,
                                Model *model,
                                int icov = 0,
                                NamingConvention namconv = NamingConvention("Nostat"));
GEOSLIB_API void set_test_discrete(bool flag_discret);
GEOSLIB_API Vario* model_pgs(Db *db,
                             const VarioParam *varioparam,
                             const RuleProp* ruleprop,
                             const Model* model1,
                             const Model* model2 = nullptr);

/******************************/
/* Functions for Anamorphosis */
/******************************/

GEOSLIB_API VectorDouble anam_selectivity(Anam *anam,
                                          int nclass,
                                          VectorDouble zcut,
                                          int flag_correct = 0,
                                          int verbose = 0);

/******************************/
/* Functions for Neighborhood */
/******************************/

GEOSLIB_API int test_neigh(Db *dbin, Db *dbout, Model *model, Neigh *neigh,
                           NamingConvention namconv = NamingConvention("Neigh"));

/**********************************/
/* High-level Interface Functions */
/**********************************/

GEOSLIB_API int migrateByAttribute(Db* db1,
                                   Db* db2,
                                   const VectorInt& iatts = VectorInt(),
                                   int ldmax = 0,
                                   const VectorDouble& dmax = VectorDouble(),
                                   int flag_fill = false,
                                   int flag_inter = false,
                                   NamingConvention namconv = NamingConvention("Migrate"));
GEOSLIB_API int migrate(Db* db1,
                        Db* db2,
                        const String& name,
                        int ldmax = 0,
                        const VectorDouble& dmax = VectorDouble(),
                        int flag_fill = 0,
                        int flag_inter = 0,
                        NamingConvention namconv = NamingConvention("Migrate"));
GEOSLIB_API int migrateByLocator(Db* db1,
                                 Db* db2,
                                 const ELoc& locatorType,
                                 int ldmax = 0,
                                 const VectorDouble& dmax = VectorDouble(),
                                 int flag_fill = false,
                                 int flag_inter = false,
                                 NamingConvention namconv = NamingConvention("Migrate"));
GEOSLIB_API int db_selhull(Db *db1,
                           Db *db2,
                           bool verbose = false,
                           NamingConvention namconv = NamingConvention("Hull",ELoc::SEL));
GEOSLIB_API void db_polygon(Db *db,
                            Polygons *polygon,
                            int flag_sel = 0,
                            int flag_period = 0,
                            int flag_nested = 0,
                            NamingConvention namconv = NamingConvention("Polygon",ELoc::SEL));
GEOSLIB_API int db_grid_fill(Db *dbgrid,
                             int mode = 0,
                             int seed = 34243,
                             int radius = 1,
                             bool verbose = false,
                             NamingConvention namconv = NamingConvention("Fill"));
GEOSLIB_API int db_grid1D_fill(Db *dbgrid,
                               int mode = 0,
                               int seed = 34243,
                               NamingConvention namconv = NamingConvention("Fill"));
GEOSLIB_API int db_duplicate(Db *db,
                             bool verbose = false,
                             double *dist = nullptr,
                             int opt_code = 0,
                             double tolcode = 0.,
                             NamingConvention namconv = NamingConvention("Duplicate",ELoc::SEL));
GEOSLIB_API int kriging(Db *dbin,
                        Db *dbout,
                        Model *model,
                        Neigh *neigh,
                        const EKrigOpt& calcul = EKrigOpt::PONCTUAL,
                        int flag_est = 1,
                        int flag_std = 1,
                        int flag_var = 0,
                        VectorInt ndisc = VectorInt(),
                        VectorInt rank_colcok = VectorInt(),
                        VectorDouble matCL = VectorDouble(),
                        NamingConvention namconv = NamingConvention("Kriging"));
GEOSLIB_API int xvalid(Db *db,
                       Model *model,
                       Neigh *neigh,
                       int flag_xvalid = 1,
                       int flag_code = 0,
                       int flag_est = 1,
                       int flag_std = 1,
                       VectorInt rank_colcok = VectorInt(),
                       NamingConvention namconv = NamingConvention("Xvalid"));
GEOSLIB_API int simtub(Db *dbin,
                       Db *dbout,
                       Model *model,
                       Neigh *neigh = nullptr,
                       int nbsimu = 1,
                       int seed = 43431,
                       int nbtuba = 100,
                       int flag_check = 0,
                       NamingConvention namconv = NamingConvention("Simu"));
GEOSLIB_API int simpgs(Db *dbin,
                       Db *dbout,
                       RuleProp *ruleprop,
                       Model *model1,
                       Model *model2,
                       Neigh *neigh,
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
                       double gibbs_eps = 1.e-3,
                       double delta = 1.,
                       NamingConvention namconv = NamingConvention("Facies",ELoc::FACIES));
GEOSLIB_API int simbipgs(Db *dbin,
                         Db *dbout,
                         RuleProp *ruleprop,
                         Model *model11,
                         Model *model12,
                         Model *model21,
                         Model *model22,
                         Neigh *neigh,
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
                         double gibbs_eps = 1.e-3,
                         NamingConvention namconv = NamingConvention("Facies",ELoc::FACIES));
GEOSLIB_API int simpgs_spde(Db *dbin,
                            Db *dbout,
                            RuleProp *ruleprop,
                            Model *model1,
                            Model *model2,
                            const String& triswitch,
                            const VectorDouble& gext,
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
GEOSLIB_API Db *db_read_csv(const char *file_name,
                            int verbose = 0,
                            int flag_header = 1,
                            int nskip = 0,
                            const char *char_sep = ",",
                            const char *char_dec = ".",
                            const char *na_string = "NA",
                            int ncol_max = -1,
                            int nrow_max = -1,
                            int flag_add_rank = 0);
GEOSLIB_API int db_write_csv(Db *db,
                             const char *filename,
                             int flag_header = 1,
                             int flag_allcol = 1,
                             int flag_coor = 1,
                             bool flag_integer = false,
                             const char *char_sep = ",",
                             const char *na_string = "NA");
GEOSLIB_API int db_proportion_estimate(Db *dbin,
                                       Db *dbout,
                                       Model *model,
                                       int niter = 100,
                                       bool verbose = false,
                                       NamingConvention namconv = NamingConvention("Prop",ELoc::P));
GEOSLIB_API int defineGeneralNeigh(int mode, Db* db, Model* model,
                                   Neigh* neigh);
GEOSLIB_API VectorInt getGeneralNeigh(Db* db, Neigh* neigh, int iech);
GEOSLIB_API int gibbs_sampler(Db* db,
                              Model* model,
                              Neigh* neigh,
                              int nbsimu,
                              int seed,
                              int nboot,
                              int niter,
                              bool flag_norm,
                              bool flag_multi_mono,
                              bool flag_propagation,
                              bool flag_sym_neigh,
                              bool flag_sym_Q,
                              int  gibbs_optstats,
                              double percent,
                              double gibbs_eps,
                              bool flag_ce,
                              bool flag_cstd,
                              bool verbose,
                              NamingConvention namconv = NamingConvention("Gibbs"));

/*****************/
/* Various Tools */
/*****************/
GEOSLIB_API int db_tool_duplicate(Db *db1,
                                  Db *db2,
                                  bool flag_same,
                                  bool verbose,
                                  int opt_code,
                                  double tolcode,
                                  double *dist,
                                  double *sel);

#endif
