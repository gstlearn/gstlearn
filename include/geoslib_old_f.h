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

#include "Basic/NamingConvention.hpp"
#include "Covariances/CovAniso.hpp"
#include "Model/Constraints.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Variogram/DirParam.hpp"
#include "Variogram/Vario.hpp"
#include "geoslib_d.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
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

/***************************************/
/* Prototyping the functions in math.c */
/***************************************/

GSTLEARN_EXPORT Id foxleg_f(Id ndat,
                            Id npar,
                            Id ncont,
                            const MatrixDense& acont,
                            VectorDouble& param,
                            VectorDouble& lower,
                            VectorDouble& upper,
                            VectorDouble& scale,
                            const Option_AutoFit& mauto,
                            Id flag_title,
                            void (*func_evaluate)(Id ndat,
                                                  Id npar,
                                                  VectorDouble& param,
                                                  VectorDouble& work),
                            VectorDouble& tabexp,
                            VectorDouble& tabwgt);

/***************************************/
/* Prototyping the functions in util.c */
/***************************************/

GSTLEARN_EXPORT VectorInt ut_split_into_two(Id ncolor,
                                            Id flag_half,
                                            Id verbose,
                                            Id* nposs);

GSTLEARN_EXPORT void ut_trace_discretize(Id nseg,
                                         const double* trace,
                                         double disc,
                                         Id* np_arg,
                                         VectorDouble& xp,
                                         VectorDouble& yp,
                                         VectorDouble& dd,
                                         VectorDouble& del,
                                         double* dist_arg);
GSTLEARN_EXPORT void ut_trace_sample(Db* db,
                                     const ELoc& ptype,
                                     Id np,
                                     const double* xp,
                                     const double* yp,
                                     const double* dd,
                                     double radius,
                                     Id* ns_arg,
                                     VectorDouble& xs,
                                     VectorDouble& ys,
                                     VectorInt& rks,
                                     VectorInt& lys,
                                     VectorInt& typ);
GSTLEARN_EXPORT double ut_distance(Id ndim, const double* tab1, const double* tab2);
GSTLEARN_EXPORT void ut_distance_allocated(Id ndim,
                                           double** tab1,
                                           double** tab2);

/*****************************************/
/* Prototyping the functions in matrix.c */
/*****************************************/

GSTLEARN_EXPORT Id matrix_invert(double* a, Id neq, Id rank);
GSTLEARN_EXPORT double matrix_determinant(Id neq, const VectorDouble& b);
GSTLEARN_EXPORT void matrix_product_safe(Id n1,
                                         Id n2,
                                         Id n3,
                                         const double* v1,
                                         const double* v2,
                                         double* v3);
GSTLEARN_EXPORT Id matrix_prod_norme(Id transpose,
                                     Id n1,
                                     Id n2,
                                     const double* v1,
                                     const double* a,
                                     double* w);
GSTLEARN_EXPORT void matrix_transpose(Id n1, Id n2, VectorDouble& v1, VectorDouble& w1);
GSTLEARN_EXPORT Id matrix_cholesky_decompose(const double* a,
                                             double* tl,
                                             Id neq);

/*****************************************/
/* Prototyping the functions in morpho.c */
/*****************************************/

GSTLEARN_EXPORT Id spill_point(DbGrid* dbgrid,
                               Id ind_depth,
                               Id ind_data,
                               Id option,
                               bool flag_up,
                               Id verbose_step,
                               double hmax,
                               double* h,
                               const double* th,
                               Id* ix0,
                               Id* iy0);

/****************************************/
/* Prototyping the functions in model.c */
/****************************************/

GSTLEARN_EXPORT Id model_fitting_sills(Vario* vario,
                                       Model* model,
                                       const Constraints& constraints = Constraints(),
                                       const Option_VarioFit& optvar  = Option_VarioFit(),
                                       const Option_AutoFit& mauto    = Option_AutoFit());
GSTLEARN_EXPORT Id model_covmat_inchol(Id verbose,
                                       Db* db,
                                       Model* model,
                                       double eta,
                                       Id npivot_max,
                                       Id nsize1,
                                       const Id* ranks1,
                                       const double* center,
                                       Id flag_sort,
                                       Id* npivot_arg,
                                       VectorInt& pvec,
                                       VectorDouble& Gmatrix,
                                       const CovCalcMode* mode = nullptr);
GSTLEARN_EXPORT void model_cova_characteristics(const ECov& type,
                                                char cov_name[STRING_LENGTH],
                                                Id* flag_range,
                                                Id* flag_param,
                                                Id* min_order,
                                                Id* max_ndim,
                                                Id* flag_int_1d,
                                                Id* flag_int_2d,
                                                Id* flag_aniso,
                                                Id* flag_rotation,
                                                double* scale,
                                                double* parmax);
GSTLEARN_EXPORT Model* model_combine(const Model* model1,
                                     const Model* model2,
                                     double r);

/*************************************/
/* Prototyping the functions in db.c */
/*************************************/

GSTLEARN_EXPORT void grid_iterator_init(Grid* grid,
                                        const VectorInt& order = VectorInt());
GSTLEARN_EXPORT VectorInt grid_iterator_next(Grid* grid);

GSTLEARN_EXPORT Id db_name_identify(Db* db, const String& string);
GSTLEARN_EXPORT Id db_locator_attribute_add(Db* db,
                                            const ELoc& locatorType,
                                            Id number,
                                            Id r_tem,
                                            double valinit,
                                            Id* iptr);
GSTLEARN_EXPORT void db_locators_correct(VectorString& strings,
                                         const VectorInt& current,
                                         Id flag_locnew);
GSTLEARN_EXPORT void db_grid_print(Db* db);

GSTLEARN_EXPORT Id db_grid_define_coordinates(DbGrid* db);
GSTLEARN_EXPORT void db_sample_print(Db* db,
                                     Id iech,
                                     Id flag_ndim = 1,
                                     Id flag_nvar = 0,
                                     Id flag_nerr = 0,
                                     Id flag_blk  = 0);
GSTLEARN_EXPORT Id db_center(Db* db, double* center);
GSTLEARN_EXPORT Id db_selref(Id ndim,
                             const Id* nx,
                             const Id* ref,
                             const double* tabin,
                             double* tabout);
GSTLEARN_EXPORT Db* db_regularize(Db* db, DbGrid* dbgrid, Id flag_center);
GSTLEARN_EXPORT Id compat_NDIM(Db* db1, Db* db2);
GSTLEARN_EXPORT double get_grid_value(DbGrid* dbgrid,
                                      Id iptr,
                                      VectorInt& indg,
                                      Id ix,
                                      Id iy,
                                      Id iz);
GSTLEARN_EXPORT void set_grid_value(DbGrid* dbgrid,
                                    Id iptr,
                                    VectorInt& indg,
                                    Id ix,
                                    Id iy,
                                    Id iz,
                                    double value);
GSTLEARN_EXPORT Id get_LOCATOR_NITEM(const Db* db, const ELoc& locatorType);
GSTLEARN_EXPORT Id is_grid_multiple(DbGrid* db1, DbGrid* db2);
GSTLEARN_EXPORT DbGrid* db_grid_reduce(DbGrid* db_grid,
                                       Id iptr,
                                       const Id* margin,
                                       const Id* limmin,
                                       Id flag_sel,
                                       Id flag_copy,
                                       Id verbose,
                                       double vmin,
                                       double vmax);
GSTLEARN_EXPORT double distance_inter(const Db* db1,
                                      const Db* db2,
                                      Id iech1,
                                      Id iech2,
                                      double* dist_vect);
GSTLEARN_EXPORT double distance_intra(const Db* db,
                                      Id iech1,
                                      Id iech2,
                                      double* dist_vect);
GSTLEARN_EXPORT double distance_grid(DbGrid* db,
                                     Id flag_moins1,
                                     Id iech1,
                                     Id iech2,
                                     double* dist_vect);
GSTLEARN_EXPORT VectorDouble db_distances_general(Db* db1,
                                                  Db* db2,
                                                  Id niso,
                                                  Id mode,
                                                  Id flag_same,
                                                  Id* n1,
                                                  Id* n2,
                                                  double* dmin,
                                                  double* dmax);
GSTLEARN_EXPORT Id point_to_grid(const DbGrid* db,
                                 const double* coor,
                                 Id flag_outside,
                                 Id* indg);
GSTLEARN_EXPORT Id point_to_bench(const DbGrid* db,
                                  double* coor,
                                  Id flag_outside,
                                  Id* indb);
GSTLEARN_EXPORT Id index_point_to_grid(const Db* db,
                                       Id iech,
                                       Id flag_outside,
                                       const DbGrid* dbout,
                                       double* coor);
GSTLEARN_EXPORT Id point_to_point(Db* db, const double* coor);
GSTLEARN_EXPORT Id point_inside_grid(Db* db, Id iech, const DbGrid* dbgrid);
GSTLEARN_EXPORT Id db_gradient_components(DbGrid* dbgrid);
GSTLEARN_EXPORT Id db_streamline(DbGrid* dbgrid,
                                 Db* dbpoint,
                                 Id niter,
                                 double step,
                                 Id flag_norm,
                                 Id use_grad,
                                 Id save_grad,
                                 Id* nbline_loc,
                                 Id* npline_loc,
                                 VectorDouble& line);
GSTLEARN_EXPORT void db_monostat(Db* db,
                                 Id iatt,
                                 double* wtot,
                                 double* mean,
                                 double* var,
                                 double* mini,
                                 double* maxi);
GSTLEARN_EXPORT Id surface(Db* db_point,
                           DbGrid* db_grid,
                           Id icol,
                           double dlim,
                           double* dtab,
                           double* gtab);
GSTLEARN_EXPORT Id db_grid_copy(DbGrid* db1,
                                DbGrid* db2,
                                const Id* ind1,
                                const Id* ind2,
                                Id ncol,
                                Id* cols);
GSTLEARN_EXPORT Id db_grid_copy_dilate(DbGrid* db1,
                                       Id iatt1,
                                       DbGrid* db2,
                                       Id iatt2,
                                       Id mode,
                                       const Id* nshift);
GSTLEARN_EXPORT Id db_proportion(Db* db,
                                 DbGrid* dbgrid,
                                 Id nfac1max,
                                 Id nfac2max,
                                 Id* nclout);
GSTLEARN_EXPORT Id db_merge(Db* db, Id ncol, Id* cols);
GSTLEARN_EXPORT Id db_count_defined(Db* db, Id icol);

GSTLEARN_EXPORT Id db_prop_read(DbGrid* db, Id ix, Id iy, double* props);
GSTLEARN_EXPORT Id db_prop_write(DbGrid* db, Id ix, Id iy, double* props);
GSTLEARN_EXPORT Id db_resind(Db* db, Id ivar, const VectorDouble& zcut);
GSTLEARN_EXPORT Id db_gradient_modang_to_component(Db* db,
                                                   Id ang_conv,
                                                   Id iad_mod,
                                                   Id iad_ang,
                                                   Id iad_gx,
                                                   Id iad_gy);
GSTLEARN_EXPORT Id db_gradient_component_to_modang(Db* db,
                                                   Id verbose,
                                                   Id iad_gx,
                                                   Id iad_gy,
                                                   Id iad_mod,
                                                   Id iad_ang,
                                                   double scale,
                                                   double ve);
GSTLEARN_EXPORT Db* db_point_init(Id nech,
                                  const VectorDouble& coormin = VectorDouble(),
                                  const VectorDouble& coormax = VectorDouble(),
                                  DbGrid* dbgrid              = nullptr,
                                  bool flag_exact             = true,
                                  bool flag_repulsion         = false,
                                  double range                = 0.,
                                  double beta                 = 0.,
                                  double extend               = 0.,
                                  Id seed                     = 43241,
                                  bool flagAddSampleRank      = true);
GSTLEARN_EXPORT Id db_smooth_vpc(DbGrid* db, Id width, double range);
GSTLEARN_EXPORT Id db_grid2point_sampling(DbGrid* dbgrid,
                                          Id nvar,
                                          Id* vars,
                                          const Id* npacks,
                                          Id npcell,
                                          Id nmini,
                                          Id* nech,
                                          VectorDouble& coor,
                                          VectorDouble& data);
GSTLEARN_EXPORT Id db_grid_patch(DbGrid* ss_grid,
                                 DbGrid* db_grid,
                                 Id iptr_ss,
                                 Id iptr_db,
                                 Id iptr_rank,
                                 Id new_rank,
                                 Id oper,
                                 Id verbose);

/****************************************/
/* Prototyping the functions in stats.c */
/****************************************/

GSTLEARN_EXPORT Id stats_residuals(Id verbose,
                                   Id nech,
                                   double* tab,
                                   Id ncut,
                                   double* zcut,
                                   Id* nsorted,
                                   double* mean,
                                   double* residuals,
                                   double* T,
                                   double* Q);
GSTLEARN_EXPORT Id db_upscale(DbGrid* dbgrid1,
                              DbGrid* dbgrid2,
                              Id orient,
                              Id verbose);
GSTLEARN_EXPORT Id db_diffusion(DbGrid* dbgrid1,
                                DbGrid* dbgrid2,
                                Id orient,
                                Id niter,
                                Id nseed,
                                Id seed,
                                Id verbose);

/****************************************/
/* Prototyping the functions in krige.c */
/****************************************/

GSTLEARN_EXPORT void set_DBIN(Db* dbin);
GSTLEARN_EXPORT void set_DBOUT(Db* dbout);
GSTLEARN_EXPORT Id krige_koption_manage(Id mode,
                                        Id flag_check,
                                        const EKrigOpt& calcul,
                                        Id flag_rand,
                                        const VectorInt& ndiscs = VectorInt());
GSTLEARN_EXPORT void krige_lhs_print(Id nech,
                                     Id neq,
                                     Id nred,
                                     const Id* flag,
                                     const double* lhs);
GSTLEARN_EXPORT void krige_rhs_print(Id nvar,
                                     Id nech,
                                     Id neq,
                                     Id nred,
                                     const Id* flag,
                                     double* rhs);
GSTLEARN_EXPORT Id krigsampling_f(Db* dbin,
                                  Db* dbout,
                                  Model* model,
                                  double beta,
                                  VectorInt& ranks1,
                                  VectorInt& ranks2,
                                  bool flag_std,
                                  Id verbose);
GSTLEARN_EXPORT Id global_transitive(DbGrid* dbgrid,
                                     Model* model,
                                     Id flag_verbose,
                                     Id flag_regular,
                                     Id ndisc,
                                     double* abundance,
                                     double* sse,
                                     double* cvtrans);
GSTLEARN_EXPORT Id anakexp_f(DbGrid* db,
                             double* covdd,
                             double* covd0,
                             double top,
                             double bot,
                             Id ncov_radius,
                             Id neigh_radius,
                             Id flag_sym,
                             Id nfeq);
GSTLEARN_EXPORT Id anakexp_3D(DbGrid* db,
                              double* cov_ref,
                              Id cov_radius,
                              Id neigh_ver,
                              Id neigh_hor,
                              Id flag_sym,
                              Model* model,
                              double nugget,
                              Id nfeq,
                              Id dbg_ix,
                              Id dbg_iy);
GSTLEARN_EXPORT Id sampling_f(Db* db,
                              Model* model,
                              double beta,
                              Id method1,
                              Id nsize1_max,
                              VectorInt& ranks1,
                              Id method2,
                              Id nsize2_max,
                              VectorInt& ranks2,
                              Id verbose);
GSTLEARN_EXPORT Id inhomogeneous_kriging(Db* dbdat,
                                         Db* dbsrc,
                                         Db* dbout,
                                         double power,
                                         Id flag_source,
                                         Model* model_dat,
                                         Model* model_src);

/*****************************************/
/* Prototyping the functions in simtub.c */
/*****************************************/

GSTLEARN_EXPORT void simu_define_func_update(void (*st_simu_update)(Db*,
                                                                    Id,
                                                                    Id,
                                                                    Id));
GSTLEARN_EXPORT void simu_define_func_scale(void (*st_simu_scale)(Db*,
                                                                  Id,
                                                                  Id));
GSTLEARN_EXPORT void simu_func_categorical_transf(Db* db,
                                                  Id verbose,
                                                  Id isimu,
                                                  Id nbsimu);
GSTLEARN_EXPORT void simu_func_continuous_update(Db* db,
                                                 Id verbose,
                                                 Id isimu,
                                                 Id nbsimu);
GSTLEARN_EXPORT void simu_func_categorical_update(Db* db,
                                                  Id verbose,
                                                  Id isimu0,
                                                  Id nbsimu0);
GSTLEARN_EXPORT void simu_func_continuous_scale(Db* db,
                                                Id verbose,
                                                Id nbsimu);
GSTLEARN_EXPORT void simu_func_categorical_scale(Db* db,
                                                 Id verbose,
                                                 Id nbsimu);

GSTLEARN_EXPORT void check_mandatory_attribute(const char* method,
                                               Db* db,
                                               const ELoc& locatorType);
GSTLEARN_EXPORT Id simcond(Db* dbin,
                           Db* dbout,
                           Model* model,
                           Id seed,
                           Id nbsimu,
                           Id nbtuba,
                           Id gibbs_nburn,
                           Id gibbs_niter,
                           Id flag_check,
                           Id flag_ce,
                           Id flag_cstd,
                           Id verbose);
GSTLEARN_EXPORT Id simmaxstable(Db* dbout,
                                Model* model,
                                double ratio,
                                Id seed,
                                Id nbtuba,
                                Id flag_simu,
                                Id flag_rank,
                                Id verbose);
GSTLEARN_EXPORT Id simRI(Db* dbout,
                         Model* model,
                         Id ncut,
                         double* zcut,
                         double* wcut,
                         Id seed,
                         Id nbtuba,
                         Id verbose);
GSTLEARN_EXPORT Id simtub_constraints(Db* dbin,
                                      Db* dbout,
                                      Model* model,
                                      ANeigh* neigh,
                                      Id seed,
                                      Id nbtuba,
                                      Id nbsimu_min,
                                      Id nbsimu_quant,
                                      Id niter_max,
                                      VectorInt& cols,
                                      Id (*func_valid)(Id flag_grid,
                                                       Id nDim,
                                                       Id nech,
                                                       Id* nx,
                                                       double* dx,
                                                       double* x0,
                                                       double nonval,
                                                       double percent,
                                                       VectorDouble& tab));
GSTLEARN_EXPORT Id db_simulations_to_ce(Db* db,
                                        const ELoc& locatorType,
                                        Id nbsimu,
                                        Id nvar,
                                        Id* iptr_ce_arg,
                                        Id* iptr_cstd_arg);

/*****************************************/
/* Prototyping the functions in simreg.c */
/*****************************************/

GSTLEARN_EXPORT Id simfine_dim(DbGrid* dbin,
                               Id nmult,
                               Id* ndim,
                               Id* ntot,
                               Id* nx,
                               double* x0,
                               double* dx);
GSTLEARN_EXPORT Id simfine_f(DbGrid* dbin,
                             Model* model,
                             const SimuRefineParam& param,
                             Id seed,
                             VectorDouble& tab);

/*****************************************/
/* Prototyping the functions in thresh.c */
/*****************************************/

GSTLEARN_EXPORT Id db_bounds_shadow(Db* db,
                                    Db* dbprop,
                                    RuleShadow* rule,
                                    Model* model,
                                    const VectorDouble& props,
                                    Id flag_stat,
                                    Id nfacies);

/*****************************************/
/* Prototyping the functions in geophy.c */
/*****************************************/

GSTLEARN_EXPORT Id time_3db(double* HS,
                            double* T,
                            Id NX,
                            Id NY,
                            Id NZ,
                            Id BX,
                            Id BY,
                            Id BZ,
                            double XS,
                            double YS,
                            double ZS,
                            double HS_EPS_INIT,
                            Id MSG);

/******************************************/
/* Prototyping the functions in mlayers.c */
/******************************************/
GSTLEARN_EXPORT Id multilayers_vario(Db* dbin,
                                     DbGrid* dbout,
                                     Vario* vario,
                                     Id nlayers,
                                     Id flag_vel,
                                     Id flag_ext,
                                     Id irf_rank,
                                     Id match_time,
                                     Id colrefd,
                                     Id colreft,
                                     Id verbose);
GSTLEARN_EXPORT Id multilayers_kriging(Db* dbin,
                                       DbGrid* dbout,
                                       Model* model,
                                       ANeigh* neigh,
                                       Id flag_same,
                                       Id flag_z,
                                       Id flag_vel,
                                       Id flag_cumul,
                                       Id flag_ext,
                                       Id flag_std,
                                       Id flag_bayes,
                                       Id irf_rank,
                                       Id match_time,
                                       Id dim_prior,
                                       double* prior_mean,
                                       double* prior_vars,
                                       Id colrefd,
                                       Id colreft,
                                       Id colrefb,
                                       Id verbose);
GSTLEARN_EXPORT Id multilayers_get_prior(Db* dbin,
                                         DbGrid* dbout,
                                         Model* model,
                                         Id flag_same,
                                         Id flag_vel,
                                         Id flag_ext,
                                         Id irf_rank,
                                         Id match_time,
                                         Id colrefd,
                                         Id colreft,
                                         Id colrefb,
                                         Id verbose,
                                         Id* npar_arg,
                                         VectorDouble& mean,
                                         VectorDouble& vars);

/***************************************/
/* Prototyping the functions in spde.c */
/***************************************/
GSTLEARN_EXPORT Id m2d_gibbs_spde(Db* dbin,
                                  Db* dbout,
                                  Model* model,
                                  Id flag_ed,
                                  Id nlayer,
                                  Id niter,
                                  Id seed,
                                  Id nbsimu,
                                  Id icol_pinch,
                                  Id flag_drift,
                                  Id flag_ce,
                                  Id flag_cstd,
                                  Id verbose);

#ifndef SWIG
GSTLEARN_EXPORT SPDE_Matelem& spde_get_current_matelem(Id icov);
#endif
GSTLEARN_EXPORT AMesh* spde_mesh_load(Db* dbin,
                                      Db* dbout,
                                      const VectorDouble& gext,
                                      SPDE_Option& s_option,
                                      bool verbose);

/******************************************/
/* Prototyping the functions in cluster.c */
/******************************************/
GSTLEARN_EXPORT VectorDouble kclusters(const VectorDouble& data,
                                       Id nvar,
                                       Id nech,
                                       Id nclusters,
                                       Id npass,
                                       Id mode,
                                       Id verbose);
GSTLEARN_EXPORT VectorInt kmedoids(const VectorDouble& data,
                                   Id nvar,
                                   Id nech,
                                   Id nclusters,
                                   Id npass,
                                   Id verbose);
} // namespace gstlrn
