/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/******************************************************************************/
#pragma once

#include "Covariances/CovAniso.hpp"
#include "gstlearn_export.hpp"
#include "geoslib_d.h"

#include "Enum/ECov.hpp"
#include "Enum/ELoc.hpp"

#include "Covariances/CovCalcMode.hpp"
#include "Basic/NamingConvention.hpp"
#include "Model/Constraints.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Variogram/DirParam.hpp"
#include "Variogram/Vario.hpp"

class AAnam;
class AnamDiscreteDD;
class AnamDiscreteIR;
class AnamEmpirical;
class AnamHermite;
class Rule;
class RuleShadow;
class MeshEStandard;
class MeshSpherical;
class CovInternal;
class Db;
class DbGrid;
class Model;
class Vario;
class VarioParam;
class ANeigh;
class NeighImage;
class NeighUnique;
class Polygons;
class PCA;
class Grid;
class SimuRefineParam;
class EStatOption;
class Faults;
class AMesh;
class SpaceTarget;

class cs;
class QChol;

/***************************************/
/* Prototyping the functions in math.c */
/***************************************/

GSTLEARN_EXPORT int foxleg_f(int ndat,
                             int npar,
                             int ncont,
                             const MatrixRectangular& acont,
                             VectorDouble &param,
                             VectorDouble &lower,
                             VectorDouble &upper,
                             VectorDouble &scale,
                             const Option_AutoFit &mauto,
                             int flag_title,
                             void (*func_evaluate)(int ndat,
                                                   int npar,
                                                   VectorDouble &param,
                                                   VectorDouble &work),
                             VectorDouble &tabexp,
                             VectorDouble &tabwgt);

/***************************************/
/* Prototyping the functions in util.c */
/***************************************/

GSTLEARN_EXPORT int* ut_split_into_two(int ncolor,
                                       int flag_half,
                                       int verbose,
                                       int *nposs);
GSTLEARN_EXPORT void projec_query(int *actif);
GSTLEARN_EXPORT void projec_print(void);
GSTLEARN_EXPORT void projec_toggle(int mode);
GSTLEARN_EXPORT void set_last_message(int mode, const char *string);
GSTLEARN_EXPORT void print_last_message(void);

GSTLEARN_EXPORT void ut_trace_discretize(int nseg,
                                         const double *trace,
                                         double disc,
                                         int *np_arg,
                                         double **xp_arg,
                                         double **yp_arg,
                                         double **dd_arg,
                                         double **del_arg,
                                         double *dist_arg);
GSTLEARN_EXPORT void ut_trace_sample(Db *db,
                                     const ELoc& ptype,
                                     int np,
                                     const double *xp,
                                     const double *yp,
                                     const double *dd,
                                     double radius,
                                     int *ns_arg,
                                     double **xs_arg,
                                     double **ys_arg,
                                     int **rks_arg,
                                     int **lys_arg,
                                     int **typ_arg);
GSTLEARN_EXPORT double ut_distance(int ndim, const double *tab1, const double *tab2);
GSTLEARN_EXPORT void ut_distance_allocated(int ndim,
                                           double **tab1,
                                           double **tab2);

/*****************************************/
/* Prototyping the functions in matrix.c */
/*****************************************/

GSTLEARN_EXPORT int matrix_invert(double *a, int neq, int rank);
GSTLEARN_EXPORT double matrix_determinant(int neq, const VectorDouble& b);
GSTLEARN_EXPORT int matrix_eigen(const double *a,
                                 int neq,
                                 double *value,
                                 double *vector);
GSTLEARN_EXPORT void matrix_product_safe(int n1,
                                         int n2,
                                         int n3,
                                         const double *v1,
                                         const double *v2,
                                         double *v3);
GSTLEARN_EXPORT int matrix_prod_norme(int transpose,
                                      int n1,
                                      int n2,
                                      const double *v1,
                                      const double *a,
                                      double *w);
GSTLEARN_EXPORT void matrix_transpose(int n1, int n2, VectorDouble& v1, VectorDouble& w1);
GSTLEARN_EXPORT int matrix_cholesky_decompose(const double *a,
                                              double *tl,
                                              int neq);
GSTLEARN_EXPORT void matrix_cholesky_product(int mode,
                                             int neq,
                                             int nrhs,
                                             const double *tl,
                                             const double *a,
                                             double *x);
GSTLEARN_EXPORT void matrix_cholesky_invert(int neq, const double *tl, double *xl);
GSTLEARN_EXPORT void matrix_combine(int nval,
                                    double coeffa,
                                    const double *a,
                                    double coeffb,
                                    const double *b,
                                    double *c);
GSTLEARN_EXPORT void matrix_product_by_diag(int mode,
                                            int neq,
                                            double *a,
                                            double *c,
                                            double *b);
GSTLEARN_EXPORT void matrix_triangle_to_square(int mode,
                                               int neq,
                                               const double *tl,
                                               double *a);
GSTLEARN_EXPORT int matrix_eigen_tridiagonal(const double *vecdiag,
                                             const double *vecinf,
                                             const double *vecsup,
                                             int neq,
                                             double *eigvec,
                                             double *eigval);

/*****************************************/
/* Prototyping the functions in morpho.c */
/*****************************************/

GSTLEARN_EXPORT int spill_point(DbGrid* dbgrid,
                                int ind_depth,
                                int ind_data,
                                int option,
                                bool flag_up,
                                int verbose_step,
                                double hmax,
                                double* h,
                                const double* th,
                                int* ix0,
                                int* iy0);

/****************************************/
/* Prototyping the functions in model.c */
/****************************************/

GSTLEARN_EXPORT int model_fitting_sills(Vario *vario,
                                        Model *model,
                                        const Constraints& constraints,
                                        const Option_AutoFit& mauto);
GSTLEARN_EXPORT int model_covmat_inchol(int verbose,
                                        Db *db,
                                        Model *model,
                                        double eta,
                                        int npivot_max,
                                        int nsize1,
                                        const int *ranks1,
                                        const double *center,
                                        int flag_sort,
                                        int *npivot_arg,
                                        int **Pret,
                                        double **Gret,
                                        const CovCalcMode* mode = nullptr);
GSTLEARN_EXPORT Model* model_duplicate_for_gradient(const Model *model,
                                       double ball_radius);
GSTLEARN_EXPORT void model_covupdt(Model *model,
                                   const double *c0,
                                   int flag_verbose,
                                   int *flag_nugget,
                                   double *nugget);
GSTLEARN_EXPORT void model_cova_characteristics(const ECov &type,
                                                char cov_name[STRING_LENGTH],
                                                int *flag_range,
                                                int *flag_param,
                                                int *min_order,
                                                int *max_ndim,
                                                int *flag_int_1d,
                                                int *flag_int_2d,
                                                int *flag_aniso,
                                                int *flag_rotation,
                                                double *scale,
                                                double *parmax);
GSTLEARN_EXPORT Model* model_combine(const Model *model1,
                                     const Model *model2,
                                     double r);

/*************************************/
/* Prototyping the functions in db.c */
/*************************************/

GSTLEARN_EXPORT void grid_iterator_init(Grid *grid,
                                        const VectorInt &order = VectorInt());
GSTLEARN_EXPORT VectorInt grid_iterator_next(Grid *grid);

GSTLEARN_EXPORT int db_name_identify(Db *db, const String &string);
GSTLEARN_EXPORT int db_locator_attribute_add(Db *db,
                                             const ELoc& locatorType,
                                             int number,
                                             int r_tem,
                                             double valinit,
                                             int *iptr);
GSTLEARN_EXPORT void db_locators_correct(VectorString &strings,
                                         const VectorInt &current,
                                         int flag_locnew);
GSTLEARN_EXPORT void db_grid_print(Db *db);

GSTLEARN_EXPORT int db_grid_define_coordinates(DbGrid *db);
GSTLEARN_EXPORT void db_sample_print(Db *db,
                                     int iech,
                                     int flag_ndim,
                                     int flag_nvar,
                                     int flag_nerr);
GSTLEARN_EXPORT int db_center(Db *db, double *center);
GSTLEARN_EXPORT void db_extension_rotated(Db *db,
                                          double *rotmat,
                                          VectorDouble& mini,
                                          VectorDouble& maxi,
                                          VectorDouble& delta);
GSTLEARN_EXPORT int db_selref(int ndim,
                              const int *nx,
                              const int *ref,
                              const double *tabin,
                              double *tabout);
GSTLEARN_EXPORT Db* db_regularize(Db *db, DbGrid *dbgrid, int flag_center);
GSTLEARN_EXPORT int compat_NDIM(Db *db1, Db *db2);
GSTLEARN_EXPORT double get_grid_value(DbGrid *dbgrid,
                                      int iptr,
                                      VectorInt& indg,
                                      int ix,
                                      int iy,
                                      int iz);
GSTLEARN_EXPORT void set_grid_value(DbGrid *dbgrid,
                                    int iptr,
                                    VectorInt& indg,
                                    int ix,
                                    int iy,
                                    int iz,
                                    double value);
GSTLEARN_EXPORT int get_LOCATOR_NITEM(const Db *db, const ELoc& locatorType);
GSTLEARN_EXPORT int is_grid_multiple(DbGrid *db1, DbGrid *db2);
GSTLEARN_EXPORT DbGrid* db_grid_reduce(DbGrid *db_grid,
                                       int iptr,
                                       const int *margin,
                                       const int *limmin,
                                       int flag_sel,
                                       int flag_copy,
                                       int verbose,
                                       double vmin,
                                       double vmax);
GSTLEARN_EXPORT double distance_inter(const Db *db1,
                                      const Db *db2,
                                      int iech1,
                                      int iech2,
                                      double *dist_vect);
GSTLEARN_EXPORT double distance_intra(const Db *db,
                                      int iech1,
                                      int iech2,
                                      double *dist_vect);
GSTLEARN_EXPORT double distance_grid(DbGrid *db,
                                     int flag_moins1,
                                     int iech1,
                                     int iech2,
                                     double *dist_vect);
GSTLEARN_EXPORT double* db_distances_general(Db *db1,
                                             Db *db2,
                                             int niso,
                                             int mode,
                                             int flag_same,
                                             int *n1,
                                             int *n2,
                                             double *dmin,
                                             double *dmax);
GSTLEARN_EXPORT int point_to_grid(const DbGrid *db,
                                  const double *coor,
                                  int flag_outside,
                                  int *indg);
GSTLEARN_EXPORT int point_to_bench(const DbGrid *db,
                                   double *coor,
                                   int flag_outside,
                                   int *indb);
GSTLEARN_EXPORT int index_point_to_grid(const Db *db,
                                        int iech,
                                        int flag_outside,
                                        const DbGrid *dbout,
                                        double *coor);
GSTLEARN_EXPORT int point_to_point(Db *db, const double *coor);
GSTLEARN_EXPORT int point_inside_grid(Db *db, int iech, const DbGrid *dbgrid);
GSTLEARN_EXPORT int db_gradient_components(DbGrid *dbgrid);
GSTLEARN_EXPORT int db_streamline(DbGrid *dbgrid,
                                  Db *dbpoint,
                                  int niter,
                                  double step,
                                  int flag_norm,
                                  int use_grad,
                                  int save_grad,
                                  int *nbline_loc,
                                  int *npline_loc,
                                  double **line_loc);
GSTLEARN_EXPORT void db_monostat(Db *db,
                                 int iatt,
                                 double *wtot,
                                 double *mean,
                                 double *var,
                                 double *mini,
                                 double *maxi);
GSTLEARN_EXPORT int db_gradient_update(Db *db);
GSTLEARN_EXPORT int surface(Db *db_point,
                            DbGrid *db_grid,
                            int icol,
                            double dlim,
                            double *dtab,
                            double *gtab);
GSTLEARN_EXPORT int db_edit(Db *db, int *flag_valid);
GSTLEARN_EXPORT int db_grid_copy(DbGrid *db1,
                                 DbGrid *db2,
                                 const int *ind1,
                                 const int *ind2,
                                 int ncol,
                                 int *cols);
GSTLEARN_EXPORT int db_grid_copy_dilate(DbGrid *db1,
                                        int iatt1,
                                        DbGrid *db2,
                                        int iatt2,
                                        int mode,
                                        const int *nshift);
GSTLEARN_EXPORT int db_proportion(Db *db,
                                  DbGrid *dbgrid,
                                  int nfac1max,
                                  int nfac2max,
                                  int *nclout);
GSTLEARN_EXPORT int db_merge(Db *db, int ncol, int *cols);
GSTLEARN_EXPORT int db_count_defined(Db *db, int icol);

GSTLEARN_EXPORT int db_prop_read(DbGrid *db, int ix, int iy, double *props);
GSTLEARN_EXPORT int db_prop_write(DbGrid *db, int ix, int iy, double *props);
GSTLEARN_EXPORT int db_resind(Db *db, int ivar, const VectorDouble& zcut);
GSTLEARN_EXPORT int db_gradient_modang_to_component(Db *db,
                                                    int ang_conv,
                                                    int iad_mod,
                                                    int iad_ang,
                                                    int iad_gx,
                                                    int iad_gy);
GSTLEARN_EXPORT int db_gradient_component_to_modang(Db *db,
                                                    int verbose,
                                                    int iad_gx,
                                                    int iad_gy,
                                                    int iad_mod,
                                                    int iad_ang,
                                                    double scale,
                                                    double ve);
GSTLEARN_EXPORT Db* db_point_init(int nech,
                                  const VectorDouble &coormin = VectorDouble(),
                                  const VectorDouble &coormax = VectorDouble(),
                                  DbGrid *dbgrid = nullptr,
                                  bool flag_exact = true,
                                  bool flag_repulsion = false,
                                  double range = 0.,
                                  double beta = 0.,
                                  double extend = 0.,
                                  int seed = 43241,
                                  bool flagAddSampleRank = true);
GSTLEARN_EXPORT int db_smooth_vpc(DbGrid *db, int width, double range);
GSTLEARN_EXPORT int db_grid2point_sampling(DbGrid *dbgrid,
                                           int nvar,
                                           int *vars,
                                           const int *npacks,
                                           int npcell,
                                           int nmini,
                                           int *nech,
                                           double **coor,
                                           double **data);
GSTLEARN_EXPORT int db_grid_patch(DbGrid* ss_grid,
                                  DbGrid* db_grid,
                                  int iptr_ss,
                                  int iptr_db,
                                  int iptr_rank,
                                  int new_rank,
                                  int oper,
                                  int verbose);

/****************************************/
/* Prototyping the functions in stats.c */
/****************************************/

GSTLEARN_EXPORT int stats_residuals(int verbose,
                                    int nech,
                                    double *tab,
                                    int ncut,
                                    double *zcut,
                                    int *nsorted,
                                    double *mean,
                                    double *residuals,
                                    double *T,
                                    double *Q);
GSTLEARN_EXPORT int db_upscale(DbGrid *dbgrid1,
                               DbGrid *dbgrid2,
                               int orient,
                               int verbose);
GSTLEARN_EXPORT int db_diffusion(DbGrid *dbgrid1,
                                 DbGrid *dbgrid2,
                                 int orient,
                                 int niter,
                                 int nseed,
                                 int seed,
                                 int verbose);


/****************************************/
/* Prototyping the functions in krige.c */
/****************************************/

GSTLEARN_EXPORT int is_flag_data_disc_defined(void);
GSTLEARN_EXPORT void set_DBIN(Db* dbin);
GSTLEARN_EXPORT void set_DBOUT(Db* dbout);
GSTLEARN_EXPORT int krige_koption_manage(int mode,
                                         int flag_check,
                                         const EKrigOpt &calcul,
                                         int flag_rand,
                                         const VectorInt& ndiscs = VectorInt());
GSTLEARN_EXPORT void krige_lhs_print(int nech,
                                     int neq,
                                     int nred,
                                     const int *flag,
                                     const double *lhs);
GSTLEARN_EXPORT void krige_rhs_print(int nvar,
                                     int nech,
                                     int neq,
                                     int nred,
                                     const int *flag,
                                     double *rhs);
GSTLEARN_EXPORT void krige_dual_print(int nech,
                                      int neq,
                                      int nred,
                                      const int *flag,
                                      double *dual);
GSTLEARN_EXPORT int bayes_simulate(Model *model,
                                   int nbsimu,
                                   const VectorDouble& rmean,
                                   const VectorDouble& rcov,
                                   VectorDouble& smean);
GSTLEARN_EXPORT int krigsampling_f(Db *dbin,
                                   Db *dbout,
                                   Model *model,
                                   double beta,
                                   VectorInt& ranks1,
                                   VectorInt& ranks2,
                                   bool flag_std,
                                   int verbose);
GSTLEARN_EXPORT int global_transitive(DbGrid* dbgrid,
                                      Model* model,
                                      int flag_verbose,
                                      int flag_regular,
                                      int ndisc,
                                      double* abundance,
                                      double* sse,
                                      double* cvtrans);
GSTLEARN_EXPORT int anakexp_f(DbGrid *db,
                              double *covdd,
                              double *covd0,
                              double top,
                              double bot,
                              int ncov_radius,
                              int neigh_radius,
                              int flag_sym,
                              int nfeq);
GSTLEARN_EXPORT int anakexp_3D(DbGrid *db,
                               double *cov_ref,
                               int cov_radius,
                               int neigh_ver,
                               int neigh_hor,
                               int flag_sym,
                               Model *model,
                               double nugget,
                               int nfeq,
                               int dbg_ix,
                               int dbg_iy);
GSTLEARN_EXPORT int sampling_f(Db* db,
                               Model* model,
                               double beta,
                               int method1,
                               int nsize1_max,
                               VectorInt& ranks1,
                               int method2,
                               int nsize2_max,
                               VectorInt& ranks2,
                               int verbose);
GSTLEARN_EXPORT int inhomogeneous_kriging(Db *dbdat,
                                          Db *dbsrc,
                                          Db *dbout,
                                          double power,
                                          int flag_source,
                                          Model *model_dat,
                                          Model *model_src);

/*****************************************/
/* Prototyping the functions in simtub.c */
/*****************************************/

GSTLEARN_EXPORT void simu_define_func_transf(void (*st_simu_transf)(Db*,
                                                                    int,
                                                                    int,
                                                                    int));
GSTLEARN_EXPORT void simu_define_func_update(void (*st_simu_update)(Db*,
                                                                    int,
                                                                    int,
                                                                    int));
GSTLEARN_EXPORT void simu_define_func_scale(void (*st_simu_scale)(Db*,
                                                                  int,
                                                                  int));
GSTLEARN_EXPORT void simu_func_categorical_transf(Db *db,
                                                  int verbose,
                                                  int isimu,
                                                  int nbsimu);
GSTLEARN_EXPORT void simu_func_continuous_update(Db *db,
                                                 int verbose,
                                                 int isimu,
                                                 int nbsimu);
GSTLEARN_EXPORT void simu_func_categorical_update(Db *db,
                                                  int verbose,
                                                  int isimu0,
                                                  int nbsimu0);
GSTLEARN_EXPORT void simu_func_continuous_scale(Db *db,
                                                int verbose,
                                                int nbsimu);
GSTLEARN_EXPORT void simu_func_categorical_scale(Db *db,
                                                 int verbose,
                                                 int nbsimu);

GSTLEARN_EXPORT void check_mandatory_attribute(const char *method,
                                               Db *db,
                                               const ELoc& locatorType);
GSTLEARN_EXPORT int simcond(Db *dbin,
                            Db *dbout,
                            Model *model,
                            int seed,
                            int nbsimu,
                            int nbtuba,
                            int gibbs_nburn,
                            int gibbs_niter,
                            int flag_check,
                            int flag_ce,
                            int flag_cstd,
                            int verbose);
GSTLEARN_EXPORT int simmaxstable(Db *dbout,
                                 Model *model,
                                 double ratio,
                                 int seed,
                                 int nbtuba,
                                 int flag_simu,
                                 int flag_rank,
                                 int verbose);
GSTLEARN_EXPORT int simRI(Db *dbout,
                          Model *model,
                          int ncut,
                          double *zcut,
                          double *wcut,
                          int seed,
                          int nbtuba,
                          int verbose);
GSTLEARN_EXPORT int simtub_constraints(Db* dbin,
                                       Db* dbout,
                                       Model* model,
                                       ANeigh* neigh,
                                       int seed,
                                       int nbtuba,
                                       int nbsimu_min,
                                       int nbsimu_quant,
                                       int niter_max,
                                       VectorInt& cols,
                                       int (*func_valid)(int flag_grid,
                                                         int nDim,
                                                         int nech,
                                                         int* nx,
                                                         double* dx,
                                                         double* x0,
                                                         double nonval,
                                                         double percent,
                                                         VectorDouble& tab));
GSTLEARN_EXPORT int db_simulations_to_ce(Db *db,
                                         const ELoc& locatorType,
                                         int nbsimu,
                                         int nvar,
                                         int *iptr_ce_arg,
                                         int *iptr_cstd_arg);

/*****************************************/
/* Prototyping the functions in simreg.c */
/*****************************************/

GSTLEARN_EXPORT int simfine_dim(DbGrid *dbin,
                                int nmult,
                                int *ndim,
                                int *ntot,
                                int *nx,
                                double *x0,
                                double *dx);
GSTLEARN_EXPORT int simfine_f(DbGrid *dbin,
                              Model *model,
                              const SimuRefineParam& param,
                              int seed,
                              VectorDouble &tab);

/*****************************************/
/* Prototyping the functions in thresh.c */
/*****************************************/

GSTLEARN_EXPORT int db_bounds_shadow(Db *db,
                                     Db *dbprop,
                                     RuleShadow *rule,
                                     Model *model,
                                     const VectorDouble &props,
                                     int flag_stat,
                                     int nfacies);

/*****************************************/
/* Prototyping the functions in geophy.c */
/*****************************************/

GSTLEARN_EXPORT int time_3db(double *HS,
                             double *T,
                             int NX,
                             int NY,
                             int NZ,
                             int BX,
                             int BY,
                             int BZ,
                             double XS,
                             double YS,
                             double ZS,
                             double HS_EPS_INIT,
                             int MSG);


/******************************************/
/* Prototyping the functions in mlayers.c */
/******************************************/
GSTLEARN_EXPORT int multilayers_vario(Db *dbin,
                                      DbGrid *dbout,
                                      Vario *vario,
                                      int nlayers,
                                      int flag_vel,
                                      int flag_ext,
                                      int irf_rank,
                                      int match_time,
                                      int colrefd,
                                      int colreft,
                                      int verbose);
GSTLEARN_EXPORT int multilayers_kriging(Db* dbin,
                                        DbGrid* dbout,
                                        Model* model,
                                        ANeigh* neigh,
                                        int flag_same,
                                        int flag_z,
                                        int flag_vel,
                                        int flag_cumul,
                                        int flag_ext,
                                        int flag_std,
                                        int flag_bayes,
                                        int irf_rank,
                                        int match_time,
                                        int dim_prior,
                                        double* prior_mean,
                                        double* prior_vars,
                                        int colrefd,
                                        int colreft,
                                        int colrefb,
                                        int verbose);
GSTLEARN_EXPORT int multilayers_get_prior(Db* dbin,
                                          DbGrid* dbout,
                                          Model* model,
                                          int flag_same,
                                          int flag_vel,
                                          int flag_ext,
                                          int irf_rank,
                                          int match_time,
                                          int colrefd,
                                          int colreft,
                                          int colrefb,
                                          int verbose,
                                          int* npar_arg,
                                          double** mean,
                                          double** vars);

/***************************************/
/* Prototyping the functions in spde.c */
/***************************************/
GSTLEARN_EXPORT QChol* qchol_manage(int mode, QChol *qchol);
GSTLEARN_EXPORT double spde_compute_correc(int ndim, double param);
GSTLEARN_EXPORT int spde_check(const Db *dbin,
                               const Db *dbout,
                               Model *model1,
                               Model *model2,
                               bool verbose,
                               const VectorDouble &gext,
                               bool mesh_dbin,
                               bool mesh_dbout,
                               bool flag_advanced,
                               bool flag_est,
                               bool flag_std,
                               bool flag_gibbs,
                               bool flag_modif);
GSTLEARN_EXPORT int spde_attach_model(Model *model);
GSTLEARN_EXPORT int m2d_gibbs_spde(Db *dbin,
                                   Db *dbout,
                                   Model *model,
                                   int flag_ed,
                                   int nlayer,
                                   int niter,
                                   int seed,
                                   int nbsimu,
                                   int icol_pinch,
                                   int flag_drift,
                                   int flag_ce,
                                   int flag_cstd,
                                   int verbose);
GSTLEARN_EXPORT SPDE_Option spde_option_alloc(void);
GSTLEARN_EXPORT void spde_option_update(SPDE_Option &s_option,
                                        const String &triswitch);
GSTLEARN_EXPORT int spde_prepar(Db *dbin,
                                Db *dbout,
                                const VectorDouble &gext,
                                SPDE_Option &s_option);
GSTLEARN_EXPORT int spde_posterior();
GSTLEARN_EXPORT int spde_process(Db *dbin,
                                 Db *dbout,
                                 SPDE_Option &s_option,
                                 int nbsimu,
                                 int ngibbs_nburn,
                                 int ngibbs_niter,
                                 int ngibbs_int);
#ifndef SWIG
GSTLEARN_EXPORT SPDE_Matelem& spde_get_current_matelem(int icov);
#endif
GSTLEARN_EXPORT AMesh* spde_mesh_load(Db *dbin,
                                      Db *dbout,
                                      const VectorDouble &gext,
                                      SPDE_Option &s_option,
                                      bool verbose = false);
GSTLEARN_EXPORT void spde_mesh_assign(AMesh *amesh,
                                      int ndim,
                                      int ncorner,
                                      int nvertex,
                                      int nmesh,
                                      const VectorInt& arg_meshes,
                                      const VectorDouble& arg_points,
                                      int verbose);
GSTLEARN_EXPORT int spde_build_matrices(Model *model, int verbose);
GSTLEARN_EXPORT int spde_eval(const VectorDouble& blin,
                              MatrixSparse *S,
                              const VectorDouble &Lambda,
                              const VectorDouble &TildeC,
                              double power,
                              VectorDouble& x,
                              VectorDouble& y);
GSTLEARN_EXPORT void spde_external_mesh_define(int icov0, AMesh *mesh);
GSTLEARN_EXPORT void spde_external_mesh_undefine(int icov0);
#ifndef SWIG
GSTLEARN_EXPORT int spde_external_copy(SPDE_Matelem &matelem, int icov0);
GSTLEARN_EXPORT MatrixSparse* spde_external_A_define(int icov0, MatrixSparse *A);
GSTLEARN_EXPORT MatrixSparse* spde_external_Q_define(int icov0, MatrixSparse *Q);
GSTLEARN_EXPORT MatrixSparse* spde_external_A_undefine(int icov0);
GSTLEARN_EXPORT MatrixSparse* spde_external_Q_undefine(int icov0);
#endif
GSTLEARN_EXPORT int kriging2D_spde(Db *dbin,
                                   Model *model,
                                   SPDE_Option &s_option,
                                   int verbose,
                                   int *nmesh_arg,
                                   int *nvertex_arg,
                                   VectorInt& meshes_arg,
                                   VectorDouble& points_arg);
#ifndef SWIG
GSTLEARN_EXPORT MatrixSparse* db_mesh_neigh(const Db *db,
                                            AMesh *amesh,
                                            double radius,
                                            int flag_exact,
                                            bool verbose,
                                            int *nactive,
                                            int **ranks);
#endif
GSTLEARN_EXPORT void spde_free_all(void);

/***************************************/
/* Prototyping the functions in math.c */
/***************************************/

GSTLEARN_EXPORT int db_trisurf(Db *db,
                               Model *model,
                               const String &triswitch,
                               int icode0,
                               int verbose,
                               int *ncode_arg,
                               int *ntri_arg,
                               int *npoint_arg,
                               double *codesel,
                               VectorInt &ntcode,
                               VectorInt &triangles,
                               VectorDouble &points);


/******************************************/
/* Prototyping the functions in cluster.c */
/******************************************/
GSTLEARN_EXPORT double* kclusters(double *data,
                                  int nvar,
                                  int nech,
                                  int nclusters,
                                  int npass,
                                  int mode,
                                  int verbose);
GSTLEARN_EXPORT int* kmedoids(double *data,
                              int nvar,
                              int nech,
                              int nclusters,
                              int npass,
                              int verbose);

