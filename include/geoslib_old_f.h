/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_d.h"

#include "Enum/EAnam.hpp"
#include "Enum/EJustify.hpp"
#include "Enum/ECalcMember.hpp"
#include "Enum/ECalcVario.hpp"
#include "Enum/ECov.hpp"
#include "Enum/ELoc.hpp"
#include "Enum/EDrift.hpp"
#include "Enum/EProcessOper.hpp"
#include "Enum/EConsElem.hpp"
#include "Enum/EConsType.hpp"
#include "Enum/ENeigh.hpp"

#include "Covariances/CovCalcMode.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/CSVformat.hpp"
#include "Model/Constraints.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Variogram/DirParam.hpp"

class AAnam;
class AnamDiscreteDD;
class AnamDiscreteIR;
class AnamEmpirical;
class AnamHermite;
class PropDef;
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

class cs;
class QChol;

/**********************************************/
/* Prototyping the functions in acknowledge.c */
/**********************************************/
GSTLEARN_EXPORT void inquire_gstlearn(char **release, char **date);

/***************************************/
/* Prototyping the functions in pile.c */
/***************************************/
GSTLEARN_EXPORT void pile_reset(int type);
GSTLEARN_EXPORT void piles_reset(void);
GSTLEARN_EXPORT int pile_next(int type);
GSTLEARN_EXPORT void pile_manage(int type, int rank, int mode, char *ptr);
GSTLEARN_EXPORT int pile_correct(int type, int rank, int mode);
GSTLEARN_EXPORT char* pile_get(int type, int rank);
GSTLEARN_EXPORT void piles_dump(void);

/**************************************/
/* Prototyping the functions in fft.c */
/**************************************/

GSTLEARN_EXPORT int fftn(int ndim,
                         const int dims[],
                         double Re[],
                         double Im[],
                         int iSign,
                         double scaling);

/***************************************/
/* Prototyping the functions in math.c */
/***************************************/

GSTLEARN_EXPORT int add_sill_constraints(Constraints& constraints,
                                              double constantSill);
GSTLEARN_EXPORT int add_unit_sill_constraints(Constraints& constraints);
GSTLEARN_EXPORT int foxleg_f(int ndat,
                             int npar,
                             int ncont,
                             const VectorDouble &acont,
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
GSTLEARN_EXPORT void set_keypair(const char *keyword,
                                 int origin,
                                 int nrow,
                                 int ncol,
                                 const double *values);
GSTLEARN_EXPORT void app_keypair(const char *keyword,
                                 int origin,
                                 int nrow,
                                 int ncol,
                                 double *values);
GSTLEARN_EXPORT void set_keypair_int(const char *keyword,
                                     int origin,
                                     int nrow,
                                     int ncol,
                                     int *values);
GSTLEARN_EXPORT void app_keypair_int(const char *keyword,
                                     int origin,
                                     int nrow,
                                     int ncol,
                                     int *values);
GSTLEARN_EXPORT double get_keypone(const char *keyword, double valdef);
GSTLEARN_EXPORT int get_keypair(const char *keyword,
                                int *nrow,
                                int *ncol,
                                double **values);
GSTLEARN_EXPORT int get_keypair_int(const char *keyword,
                                    int *nrow,
                                    int *ncol,
                                    int **values);
GSTLEARN_EXPORT void del_keypair(const char *keyword, int flag_exact);
GSTLEARN_EXPORT void print_keypair(int flag_short);
GSTLEARN_EXPORT void print_range(const char *title,
                                 int ntab,
                                 double *tab,
                                 double *sel);
GSTLEARN_EXPORT void ut_trace_discretize(int nseg,
                                         double *trace,
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
                                     double *xp,
                                     double *yp,
                                     double *dd,
                                     double radius,
                                     int *ns_arg,
                                     double **xs_arg,
                                     double **ys_arg,
                                     int **rks_arg,
                                     int **lys_arg,
                                     int **typ_arg);
GSTLEARN_EXPORT int solve_P2(double a, double b, double c, double *x);
GSTLEARN_EXPORT int solve_P3(double a, double b, double c, double d, double *x);
GSTLEARN_EXPORT double ut_distance(int ndim, const double *tab1, const double *tab2);
GSTLEARN_EXPORT void ut_distance_allocated(int ndim,
                                           double **tab1,
                                           double **tab2);

/*************************************/
/* Prototyping the functions in io.c */
/*************************************/

GSTLEARN_EXPORT void record_close(void);
#ifndef SWIG
GSTLEARN_EXPORT void redefine_message(void (*write_func)(const char*));
GSTLEARN_EXPORT void redefine_error(void (*warn_func)(const char*));
GSTLEARN_EXPORT void redefine_read(void (*read_func)(const char*, char*));
GSTLEARN_EXPORT void redefine_exit(void (*exit_func)(void));
#endif
GSTLEARN_EXPORT void mem_error(int nbyte);

GSTLEARN_EXPORT void message_extern(const char *string);
GSTLEARN_EXPORT void exit_extern();

GSTLEARN_EXPORT void string_strip_blanks(char *string, int flag_lead);
GSTLEARN_EXPORT void string_strip_quotes(char *string);

#if defined(_WIN32) || defined(_WIN64)
GSTLEARN_EXPORT char * strsep(char **stringp, const char* delim);
#endif
GSTLEARN_EXPORT void print_current_line(void);
GSTLEARN_EXPORT void file_dump(int ntab, double *tab);

/*****************************************/
/* Prototyping the functions in memory.c */
/*****************************************/

/* Overwriting memory management functions */

#define mem_free(tab)          mem_free_(__FILE__,__LINE__,tab)
#define mem_alloc(a,b)         mem_alloc_(__FILE__,__LINE__,a,b)
#define mem_calloc(a,b,c)      mem_calloc_(__FILE__,__LINE__,a,b,c)
#define mem_realloc(tab,a,b)   mem_realloc_(__FILE__,__LINE__,tab,a,b)
#define mem_copy(tab,a,b)      mem_copy_(__FILE__,__LINE__,tab,a,b)

GSTLEARN_EXPORT void memory_leak_set(int flag);
GSTLEARN_EXPORT void memory_leak_reset(void);
GSTLEARN_EXPORT void memory_leak_report(void);
GSTLEARN_EXPORT char* mem_alloc_(const char *call_file,
                                 unsigned int call_line,
                                 int size,
                                 int flag_fatal);
GSTLEARN_EXPORT char* mem_calloc_(const char *call_file,
                                  unsigned int call_line,
                                  int size_t,
                                  int size,
                                  int flag_fatal);
GSTLEARN_EXPORT char* mem_realloc_(const char *call_file,
                                   unsigned int call_line,
                                   char *tab,
                                   int size,
                                   int flag_fatal);
GSTLEARN_EXPORT char* mem_copy_(const char *call_file,
                                unsigned int call_line,
                                char *tabin,
                                int size,
                                int flag_fatal);
GSTLEARN_EXPORT char* mem_free_(const char *call_file,
                                unsigned int call_line,
                                char *tab);
GSTLEARN_EXPORT void mem_debug_set(int flag);
GSTLEARN_EXPORT void memory_status(const char *title);
GSTLEARN_EXPORT double** mem_tab_free(double **tab, int nvar);
GSTLEARN_EXPORT double** mem_tab_alloc(int nvar, int size, int flag_fatal);
GSTLEARN_EXPORT void time_start(void);
GSTLEARN_EXPORT void time_reset(void);
GSTLEARN_EXPORT void time_chunk_add(const char *call_name);
GSTLEARN_EXPORT void time_report(void);

/*****************************************/
/* Prototyping the functions in matrix.c */
/*****************************************/

GSTLEARN_EXPORT int matrix_get_extreme(int mode, int ntab, double *tab);
GSTLEARN_EXPORT void matrix_invsign(int ndim, double *a);
GSTLEARN_EXPORT int matrix_invert(double *a, int neq, int rank);
GSTLEARN_EXPORT int matrix_invert_triangle(int neq, double *tl, int rank);
GSTLEARN_EXPORT int matrix_invert_copy(const double *a, int neq, double *b);
GSTLEARN_EXPORT int matrix_invsym(double *a, int neq);
GSTLEARN_EXPORT int matrix_invgen(double *a,
                                  int neq,
                                  double *tabout,
                                  double *cond);
GSTLEARN_EXPORT int matrix_invsvdsym(double *a, int neq, int rank);
GSTLEARN_EXPORT int matrix_invreal(double *mat, int neq);
GSTLEARN_EXPORT int matrix_invreal_copy(const double *a, int neq, double *b);
GSTLEARN_EXPORT void matrix_svd_inverse(int neq,
                                        double *s,
                                        double *u,
                                        double *v,
                                        double *tabout);
GSTLEARN_EXPORT double matrix_determinant(int neq, const double *b);
GSTLEARN_EXPORT int matrix_cofactor(int neq, double *a, double *b);
GSTLEARN_EXPORT double matrix_cholesky_determinant(int neq, const double *tl);
GSTLEARN_EXPORT int matrix_eigen(const double *a,
                                 int neq,
                                 double *value,
                                 double *vector);
GSTLEARN_EXPORT int matrix_geigen(const double *a,
                                  const double *b,
                                  int neq,
                                  double *value,
                                  double *vector);
GSTLEARN_EXPORT void matrix_product(int n1,
                                    int n2,
                                    int n3,
                                    const double *v1,
                                    const double *v2,
                                    double *v3);
GSTLEARN_EXPORT void matrix_product_safe(int n1,
                                         int n2,
                                         int n3,
                                         const double *v1,
                                         const double *v2,
                                         double *v3);
GSTLEARN_EXPORT int matrix_prod_norme(int tranpose,
                                      int n1,
                                      int n2,
                                      const double *v1,
                                      const double *a,
                                      double *v2);
GSTLEARN_EXPORT void matrix_transpose(int n1, int n2, double *v1, double *w1);
GSTLEARN_EXPORT void matrix_int_transpose(int n1, int n2, int *v1, int *w1);
GSTLEARN_EXPORT void matrix_transpose_in_place(int n1, int n2, double *v1);
GSTLEARN_EXPORT void matrix_int_transpose_in_place(int n1, int n2, int *v1);
GSTLEARN_EXPORT int matrix_solve(int mode,
                                 const double *a,
                                 const double *b,
                                 double *x,
                                 int neq,
                                 int nrhs,
                                 int *pivot);
GSTLEARN_EXPORT double matrix_normA(double *b, double *a, int neq, int subneq);
GSTLEARN_EXPORT int matrix_cholesky_decompose(const double *a,
                                              double *tl,
                                              int neq);
GSTLEARN_EXPORT int matrix_LU_decompose(int neq,
                                        const double *a,
                                        double *tl,
                                        double *tu);
GSTLEARN_EXPORT int matrix_LU_solve(int neq,
                                    const double *tu,
                                    const double *tl,
                                    const double *b,
                                    double *x);
GSTLEARN_EXPORT int matrix_LU_invert(int neq, double* a);
GSTLEARN_EXPORT void matrix_cholesky_product(int mode,
                                             int neq,
                                             int nrhs,
                                             const double *tl,
                                             const double *a,
                                             double *x);
GSTLEARN_EXPORT int matrix_cholesky_solve(int neq,
                                          const double *tl,
                                          const double *b,
                                          double *x);
GSTLEARN_EXPORT int matrix_cholesky_to_invert(int neq, const double *tl, double *xl);
GSTLEARN_EXPORT void matrix_cholesky_invert(int neq, const double *tl, double *xl);
GSTLEARN_EXPORT void matrix_cholesky_norme(int mode,
                                           int neq,
                                           const double *tl,
                                           const double *a,
                                           double *b);
GSTLEARN_EXPORT void matrix_triangular_product(int neq,
                                               int mode,
                                               const double *al,
                                               const double *b,
                                               double *x);
GSTLEARN_EXPORT int is_matrix_definite_positive(int neq,
                                                const double *a,
                                                double *valpro,
                                                double *vecpro,
                                                int verbose);
GSTLEARN_EXPORT int is_matrix_non_negative(int nrow,
                                           int ncol,
                                           double *a,
                                           int verbose);
GSTLEARN_EXPORT int is_matrix_null(int nrow, int ncol, const double *a, int verbose);
GSTLEARN_EXPORT int is_matrix_symmetric(int neq, const double *a, int verbose);
GSTLEARN_EXPORT int is_matrix_correlation(int neq, double *a);
GSTLEARN_EXPORT int is_matrix_rotation(int neq, const double *a, int verbose);
GSTLEARN_EXPORT void matrix_produit_cholesky(int neq, const double *tl, double *a);
GSTLEARN_EXPORT VectorDouble matrix_produit_cholesky_VD(int neq, const double *tl);
GSTLEARN_EXPORT void matrix_set_identity(int neq, double *a);
GSTLEARN_EXPORT int is_matrix_product_identity(int neq,
                                               double *a,
                                               double *b,
                                               double *errmax);
GSTLEARN_EXPORT double* matrix_bind(int mode,
                                    int n11,
                                    int n12,
                                    double *a1,
                                    int n21,
                                    int n22,
                                    double *a2,
                                    int *n31,
                                    int *n32);
GSTLEARN_EXPORT void matrix_manage(int nrows,
                                   int ncols,
                                   int nr,
                                   int nc,
                                   int *rowsel,
                                   int *colsel,
                                   double *v1,
                                   double *v2);
GSTLEARN_EXPORT void matrix_combine(int nval,
                                    double coeffa,
                                    double *a,
                                    double coeffb,
                                    double *b,
                                    double *c);
GSTLEARN_EXPORT void matrix_fill_symmetry(int neq, double *a);
GSTLEARN_EXPORT double matrix_norminf(int neq, double *a);
GSTLEARN_EXPORT double matrix_norml1(int neq, double *a);
GSTLEARN_EXPORT void matrix_square(int neq, double *a, double *b);
GSTLEARN_EXPORT VectorDouble matrix_square_VD(int neq, const VectorDouble &a);
GSTLEARN_EXPORT void matrix_product_by_diag(int mode,
                                            int neq,
                                            double *a,
                                            double *c,
                                            double *b);
GSTLEARN_EXPORT void matrix_product_by_diag_VD(int mode,
                                               int neq,
                                               VectorDouble a,
                                               const VectorDouble &c);
GSTLEARN_EXPORT void matrix_triangle_to_square(int mode,
                                               int neq,
                                               const double *tl,
                                               double *a);
GSTLEARN_EXPORT void matrix_tri2sq(int neq, const double *tl, double *a);
GSTLEARN_EXPORT void matrix_square_to_triangle(int mode,
                                               int neq,
                                               const double *a,
                                               double *tl);
GSTLEARN_EXPORT void matrix_tl2tu(int neq, const double *tl, double *tu);
GSTLEARN_EXPORT void matrix_linear(int neq,
                                   double a1,
                                   double *a,
                                   double b1,
                                   double *b,
                                   double *x);
GSTLEARN_EXPORT int matrix_eigen_tridiagonal(const double *vecdiag,
                                             const double *vecinf,
                                             const double *vecsup,
                                             int neq,
                                             double *eigvec,
                                             double *eigval);
GSTLEARN_EXPORT int matrix_qo(int neq,
                              double *hmat,
                              double *gmat,
                              double *xmat);
GSTLEARN_EXPORT int matrix_qoc(int flag_invert,
                               int neq,
                               double *hmat,
                               double *gmat,
                               int na,
                               double *amat,
                               double *bmat,
                               double *xmat,
                               double *lambda);
GSTLEARN_EXPORT int matrix_qoci(int neq,
                                double *hmat,
                                double *gmat,
                                int nae,
                                double *aemat,
                                double *bemat,
                                int nai,
                                double *aimat,
                                double *bimat,
                                double *xmat);
GSTLEARN_EXPORT void matrix_range(int n1,
                                  int n2,
                                  double *v1,
                                  double *mini,
                                  double *maxi,
                                  double *norme1,
                                  double *norme2);

/****************************************/
/* Prototyping the functions in ascii.c */
/****************************************/

GSTLEARN_EXPORT void ascii_study_define(const char *study);
GSTLEARN_EXPORT void ascii_environ_read(char *file_name, int verbose);
GSTLEARN_EXPORT void ascii_filename(const char *type,
                                    int rank,
                                    int mode,
                                    char *filename);
GSTLEARN_EXPORT void ascii_simu_read(char *file_name,
                                     int verbose,
                                     int *nbsimu,
                                     int *nbtuba,
                                     int *seed);
GSTLEARN_EXPORT int ascii_option_defined(const char *file_name,
                                         int verbose,
                                         const char *option_name,
                                         int type,
                                         void *answer);

/*****************************************/
/* Prototyping the functions in morpho.c */
/*****************************************/

GSTLEARN_EXPORT int spill_point(DbGrid *dbgrid,
                                int ind_height,
                                int ind_data,
                                int option,
                                bool flag_up,
                                int verbose,
                                double hmax,
                                double *h,
                                double *th,
                                int *ix0,
                                int *iy0);

/****************************************/
/* Prototyping the functions in vario.c */
/****************************************/

GSTLEARN_EXPORT void vario_fix_codir(int ndim, VectorDouble &codir);

GSTLEARN_EXPORT int variogram_maximum_dist1D_reached(Db *db,
                                                     int iech,
                                                     int jech,
                                                     double maxdist);
GSTLEARN_EXPORT int geometry_compute(Db *db,
                                     Vario *vario,
                                     Vario_Order *vorder,
                                     int *npair);
GSTLEARN_EXPORT int variovect_compute(Db *db, Vario *vario, int ncomp);
GSTLEARN_EXPORT void variogram_extension(const Vario *vario,
                                         int ivar,
                                         int jvar,
                                         int idir0,
                                         int flag_norm,
                                         int flag_vars,
                                         double distmin,
                                         double distmax,
                                         double varmin,
                                         double varmax,
                                         int *flag_hneg,
                                         int *flag_gneg,
                                         double *c0,
                                         double *hmin,
                                         double *hmax,
                                         double *gmin,
                                         double *gmax);
GSTLEARN_EXPORT int code_comparable(const Db *db1,
                                    const Db *db2,
                                    int iech,
                                    int jech,
                                    int opt_code,
                                    int tolcode);
GSTLEARN_EXPORT int variogram_reject_pair(const Db *db,
                                          int iech,
                                          int jech,
                                          double dist,
                                          double psmin,
                                          double bench,
                                          double cylrad,
                                          const VectorDouble &codir,
                                          double *ps);
GSTLEARN_EXPORT bool variogram_reject_fault(const Db *db,
                                            int iech,
                                            int jech,
                                            const Faults *faults = nullptr);
GSTLEARN_EXPORT void variogram_scale(Vario *vario, int idir);
GSTLEARN_EXPORT int variogram_get_lag(const DirParam& dirparam,
                                      int idir,
                                      double ps,
                                      double psmin,
                                      double *dist,
                                      bool flag_asym);
GSTLEARN_EXPORT ECalcVario vario_identify_calcul_type(const String &cov_name);
GSTLEARN_EXPORT void vardir_print(Vario *vario, int idir, int verbose);
GSTLEARN_EXPORT void vardir_copy(VarioParam *vario_in,
                                 int idir_in,
                                 VarioParam *vario_out,
                                 int idir_out);
GSTLEARN_EXPORT void variogram_trans_cut(Vario *vario, int nh, double ycut);
GSTLEARN_EXPORT int correlation_f(Db *db1,
                                  Db *db2,
                                  DbGrid *dbgrid,
                                  int flag_same,
                                  int icol1,
                                  int icol2,
                                  int flag_verbose,
                                  double dmin,
                                  double dmax,
                                  double tolang,
                                  double slice_bench,
                                  double slice_radius,
                                  VectorDouble &codir,
                                  int flag_code,
                                  int tolcode,
                                  int *nindice,
                                  int **indices,
                                  double *correl);
GSTLEARN_EXPORT int correlation_ident(Db *db1,
                                      Db *db2,
                                      int icol1,
                                      int icol2,
                                      Polygons *polygon);
GSTLEARN_EXPORT int variogram_cloud_dim(Db *db,
                                        const VarioParam *varioparam,
                                        double *vmax);
GSTLEARN_EXPORT void variogram_cloud_ident(Db *db,
                                           DbGrid *dbgrid,
                                           Vario *vario,
                                           Polygons *polygon);
GSTLEARN_EXPORT void condexp(Db *db1,
                             Db *db2,
                             int icol1,
                             int icol2,
                             double mini,
                             double maxi,
                             int nclass,
                             int verbose,
                             int *ncond,
                             double *xcond,
                             double *ycond);
GSTLEARN_EXPORT int vario_extract(Vario *vario,
                                  ECalcVario *calcul_type,
                                  int *ndim,
                                  int *nvar,
                                  int *ndir,
                                  int *ndate,
                                  double *scale,
                                  double **dates);
GSTLEARN_EXPORT int vario_get_rank(Vario *vario, int idir, int idate);
GSTLEARN_EXPORT int variogram_y2z(Vario *vario, AAnam *anam, Model *model);

/****************************************/
/* Prototyping the functions in model.c */
/****************************************/

GSTLEARN_EXPORT void model_covtab_init(int flag_init,
                                       Model *model,
                                       double *covtab);
GSTLEARN_EXPORT Model* model_default(int ndim, int nvar);
GSTLEARN_EXPORT double model_calcul_basic(Model *model,
                                          int icov,
                                          const ECalcMember &member,
                                          const VectorDouble &d1);
GSTLEARN_EXPORT double model_calcul_stdev(Model *model,
                                          Db *db1,
                                          int iech1,
                                          Db *db2,
                                          int iech2,
                                          int verbose,
                                          double factor,
                                          const CovCalcMode* mode = nullptr);
GSTLEARN_EXPORT void model_calcul_drift(Model *model,
                                        const ECalcMember &member,
                                        const Db *db,
                                        int iech,
                                        double *drftab);
GSTLEARN_EXPORT void model_nostat_update(CovInternal *covint, Model *model);
GSTLEARN_EXPORT int model_add_cova(Model *model,
                                   const ECov &type,
                                   int flag_anisotropy,
                                   int flag_rotation,
                                   double range,
                                   double param,
                                   const VectorDouble &aniso_ranges,
                                   const VectorDouble &aniso_rotmat,
                                   const VectorDouble &coreg,
                                   double ball_radius);
GSTLEARN_EXPORT int model_sample(Vario *vario,
                                 Model *model,
                                 const CovCalcMode* mode);
GSTLEARN_EXPORT void model_calcul_cov(CovInternal *covint,
                                      Model *model,
                                      const CovCalcMode* mode,
                                      int flag_init,
                                      double weight,
                                      VectorDouble d1,
                                      double *covtab);
GSTLEARN_EXPORT int model_fitting_sills(Vario *vario,
                                        Model *model,
                                        const Constraints& constraints,
                                        const Option_AutoFit& mauto);
GSTLEARN_EXPORT int model_nfex(Model *model);
GSTLEARN_EXPORT int model_update_coreg(Model *model,
                                       double *aic,
                                       double *valpro,
                                       double *vecpro);
GSTLEARN_EXPORT int model_evaluate(Model *model,
                                   int ivar,
                                   int jvar,
                                   const CovCalcMode* mode,
                                   int nh,
                                   VectorDouble &codir,
                                   const double *h,
                                   double *g);
GSTLEARN_EXPORT int model_evaluate_nostat(Model *model,
                                          int ivar,
                                          int jvar,
                                          const CovCalcMode* mode,
                                          Db *db1,
                                          int iech1,
                                          Db *db2,
                                          int iech2,
                                          int nh,
                                          VectorDouble &codir,
                                          double *h,
                                          double *g);
GSTLEARN_EXPORT int model_grid(Model *model,
                               Db *db,
                               int ivar,
                               int jvar,
                               const CovCalcMode* mode,
                               double *g);
GSTLEARN_EXPORT double model_cxx(Model *model,
                                 Db *db1,
                                 Db *db2,
                                 int ivar,
                                 int jvar,
                                 int seed,
                                 double epsdist,
                                 const CovCalcMode* mode = nullptr);
GSTLEARN_EXPORT int model_covmat(Model *model,
                                  Db *db1,
                                  Db *db2,
                                  int ivar,
                                  int jvar,
                                  double *covmat,
                                  const CovCalcMode* mode = nullptr);
GSTLEARN_EXPORT MatrixSquareSymmetric model_covmatM(Model *model,
                                                    Db *db1,
                                                    Db *db2,
                                                    int ivar0,
                                                    int jvar0,
                                                    const CovCalcMode* mode = nullptr);
GSTLEARN_EXPORT double* model_covmat_by_ranks(Model *model,
                                              Db *db1,
                                              int nsize1,
                                              const int *ranks1,
                                              Db *db2,
                                              int nsize2,
                                              const int *ranks2,
                                              int ivar0 = -1,
                                              int jvar0 = -1,
                                              const CovCalcMode* mode = nullptr);
#ifndef SWIG
GSTLEARN_EXPORT cs* model_covmat_by_ranks_cs(Model *model,
                                             Db *db1,
                                             int nsize1,
                                             const int *ranks1,
                                             Db *db2,
                                             int nsize2,
                                             const int *ranks2,
                                             int ivar0 = -1,
                                             int jvar0 = -1,
                                             const CovCalcMode* mode = nullptr);
#endif
GSTLEARN_EXPORT int model_covmat_inchol(int verbose,
                                        Db *db,
                                        Model *model,
                                        double eta,
                                        int npivot_max,
                                        int nsize1,
                                        int *ranks1,
                                        double *center,
                                        int flag_sort,
                                        int *npivots,
                                        int **Pret,
                                        double **Gret,
                                        const CovCalcMode* mode = nullptr);
GSTLEARN_EXPORT int model_drift_mat(Model *model,
                                     const ECalcMember &member,
                                     Db *db,
                                     double *drfmat);
GSTLEARN_EXPORT int model_drift_vector(Model *model,
                                        const ECalcMember &member,
                                        Db *db,
                                        int iech,
                                        double *vector);
GSTLEARN_EXPORT void model_drift_filter(Model *model, int rank, int filter);
GSTLEARN_EXPORT Model* model_duplicate(const Model *model,
                                       double ball_radius,
                                       int mode);
GSTLEARN_EXPORT int model_stabilize(Model *model,
                                    int flag_verbose,
                                    double percent);
GSTLEARN_EXPORT int model_normalize(Model *model, int flag_verbose);
GSTLEARN_EXPORT void model_covupdt(Model *model,
                                   double *c0,
                                   int flag_verbose,
                                   int *flag_nugget,
                                   double *nugget);
GSTLEARN_EXPORT double model_drift_evaluate(int verbose,
                                            Model *model,
                                            const Db *db,
                                            int iech,
                                            int ivar,
                                            double *coef,
                                            double *drftab);
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
GSTLEARN_EXPORT double cova_get_scale_factor(const ECov &type, double param);
GSTLEARN_EXPORT Model* model_combine(const Model *model1,
                                     const Model *model2,
                                     double r);
GSTLEARN_EXPORT int model_get_nonugget_cova(Model *model);
GSTLEARN_EXPORT int model_regularize(Model *model,
                                     Vario *vario,
                                     DbGrid *dbgrid,
                                     const CovCalcMode* mode = nullptr);
GSTLEARN_EXPORT double constraints_get(const Constraints &constraints,
                                       const EConsType &icase,
                                       int igrf,
                                       int istr,
                                       const EConsElem &icons,
                                       int v1,
                                       int v2);
GSTLEARN_EXPORT void constraints_print(const Constraints &constraints);
GSTLEARN_EXPORT int modify_constraints_on_sill(Constraints &constraints);
GSTLEARN_EXPORT const CovInternal* get_external_covariance();

/***************************************/
/* Prototyping the functions in anam.c */
/***************************************/

GSTLEARN_EXPORT int anam_point_to_block(AAnam *anam,
                                        int verbose,
                                        double cvv,
                                        double coeff,
                                        double mu);

/*************************************/
/* Prototyping the functions in db.c */
/*************************************/

GSTLEARN_EXPORT void grid_iterator_init(Grid *grid,
                                        const VectorInt &order = VectorInt());
GSTLEARN_EXPORT VectorInt grid_iterator_next(Grid *grid);

GSTLEARN_EXPORT int* db_indg_alloc(const Db *db);
GSTLEARN_EXPORT int* db_indg_free(int *indice);
GSTLEARN_EXPORT double* db_sample_free(double *tab);
GSTLEARN_EXPORT double* db_sample_alloc(const Db *db, const ELoc& locatorType);
GSTLEARN_EXPORT int db_sample_load(Db *db,
                                   const ELoc& locatorType,
                                   int iech,
                                   double *tab);
GSTLEARN_EXPORT double* db_vector_free(double *tab);
GSTLEARN_EXPORT double* db_vector_alloc(const Db *db);
GSTLEARN_EXPORT int db_selection_get(const Db *db, int item, double *tab);
GSTLEARN_EXPORT int db_vector_get(Db *db,
                                  const ELoc& locatorType,
                                  int locatorIndex,
                                  double *tab);
GSTLEARN_EXPORT int db_vector_put(Db *db,
                                  const ELoc& locatorType,
                                  int locatorIndex,
                                  double *tab);
GSTLEARN_EXPORT int db_vector_get_att_sel_compress(Db *db,
                                                   int icol,
                                                   int *number,
                                                   double *tab);
GSTLEARN_EXPORT int db_vector_get_att(const Db *db, int iatt, double *tab);
GSTLEARN_EXPORT int db_vector_get_att_sel(Db *db, int iatt, double *tab);
GSTLEARN_EXPORT int db_name_set(Db *db, int iatt, const String &name);
GSTLEARN_EXPORT String db_name_get_by_att(const Db *db, int iatt);
GSTLEARN_EXPORT String db_name_get_by_col(Db *db, int icol);
GSTLEARN_EXPORT int db_name_identify(Db *db, const String &string);
GSTLEARN_EXPORT void db_attribute_del_mult(Db *db, int i_del, int n_del);
GSTLEARN_EXPORT void db_attribute_init(Db *db,
                                       int ncol,
                                       int iatt,
                                       double valinit);
GSTLEARN_EXPORT void db_attribute_copy(Db *db, int iatt_in, int iatt_out);
GSTLEARN_EXPORT int db_attribute_identify(const Db *db,
                                          const ELoc& locatorType,
                                          int locatorIndex);
GSTLEARN_EXPORT int db_sample_get_att(Db *db,
                                      int iech,
                                      int number,
                                      int iatt,
                                      double *tab);
GSTLEARN_EXPORT void db_sample_put_att(Db *db,
                                       int iech,
                                       int number,
                                       int iatt,
                                       double *tab);
GSTLEARN_EXPORT int db_locator_attribute_add(Db *db,
                                             const ELoc& locatorType,
                                             int number,
                                             int r_tem,
                                             double valinit,
                                             int *iptr);
GSTLEARN_EXPORT void db_locators_correct(VectorString &strings,
                                         const VectorInt &current,
                                         int flag_locnew);
GSTLEARN_EXPORT int db_coorvec_put(Db *db, int idim, double *tab);
GSTLEARN_EXPORT int db_coorvec_get(const Db *db, int idim, double *tab);
GSTLEARN_EXPORT int db_grid_match(DbGrid *db1, DbGrid *db2);
GSTLEARN_EXPORT int db_is_isotropic(const Db *db, int iech, double *data);
GSTLEARN_EXPORT void db_grid_print(Db *db);

GSTLEARN_EXPORT DbGrid* db_create_grid_multiple(DbGrid *dbin,
                                                const VectorInt &nmult,
                                                int flag_add_rank);
GSTLEARN_EXPORT DbGrid* db_create_grid_divider(DbGrid *dbin,
                                               const VectorInt &nmult,
                                               int flag_add_rank);
GSTLEARN_EXPORT DbGrid* db_create_grid_dilate(DbGrid *dbin,
                                              int mode,
                                              const VectorInt &nshift,
                                              int flag_add_rank);
GSTLEARN_EXPORT DbGrid* db_grid_sample(DbGrid *dbin, const VectorInt &nmult);
GSTLEARN_EXPORT int db_grid_define_coordinates(DbGrid *db);
GSTLEARN_EXPORT Db* db_create_from_target(const double *target,
                                          int ndim,
                                          int flag_add_rank);
GSTLEARN_EXPORT void db_sample_print(Db *db,
                                     int iech,
                                     int flag_ndim,
                                     int flag_nvar,
                                     int flag_nerr);
GSTLEARN_EXPORT int db_center(Db *db, double *center);
GSTLEARN_EXPORT void db_extension(const Db *db,
                                  VectorDouble& mini,
                                  VectorDouble& maxi,
                                  bool flag_preserve = false);
GSTLEARN_EXPORT void db_extension_rotated(Db *db,
                                          double *rotmat,
                                          VectorDouble& mini,
                                          VectorDouble& maxi,
                                          VectorDouble& delta);
GSTLEARN_EXPORT int db_attribute_range(const Db *db,
                                       int icol,
                                       double *mini,
                                       double *maxi,
                                       double *delta);
GSTLEARN_EXPORT int db_extension_diag(const Db *db, double *diag);
GSTLEARN_EXPORT double db_epsilon_distance(Db *db);
GSTLEARN_EXPORT int db_index_grid_to_sample(const DbGrid *db, const int *indg);
GSTLEARN_EXPORT void db_index_sample_to_grid(const DbGrid *db, int iech, int *indg);
GSTLEARN_EXPORT int db_index_sorted_in_grid(const DbGrid *db, int iech, int *indg);
GSTLEARN_EXPORT int db_selref(int ndim,
                              int *nx,
                              int *ref,
                              double *tabin,
                              double *tabout);
GSTLEARN_EXPORT Db* db_extract(Db *db, int *ranks);
GSTLEARN_EXPORT Db* db_regularize(Db *db, DbGrid *dbgrid, int flag_center);
GSTLEARN_EXPORT int compat_NDIM(Db *db1, Db *db2);
GSTLEARN_EXPORT double get_grid_value(DbGrid *dbgrid,
                                      int iptr,
                                      int *indg,
                                      int ix,
                                      int iy,
                                      int iz);
GSTLEARN_EXPORT void set_grid_value(DbGrid *dbgrid,
                                    int iptr,
                                    int *indg,
                                    int ix,
                                    int iy,
                                    int iz,
                                    double value);
GSTLEARN_EXPORT int get_LOCATOR_NITEM(const Db *db, const ELoc& locatorType);
GSTLEARN_EXPORT int exist_LOCATOR(Db *db, const ELoc& locatorType);
GSTLEARN_EXPORT double get_LOCATOR_ITEM(Db *db,
                                        const ELoc& locatorType,
                                        int locatorIndex,
                                        int iech);
GSTLEARN_EXPORT void set_LOCATOR_ITEM(Db *db,
                                      const ELoc& locatorType,
                                      int locatorIndex,
                                      int iech,
                                      double value);
GSTLEARN_EXPORT int db_get_rank_absolute_to_relative(Db *db, int iech0);
GSTLEARN_EXPORT int db_get_rank_relative_to_absolute(Db *db, int iech0);
GSTLEARN_EXPORT int is_grid_multiple(DbGrid *db1, DbGrid *db2);
GSTLEARN_EXPORT int db_grid_copy_params(DbGrid *dbin, int mode, DbGrid *dbout);
GSTLEARN_EXPORT DbGrid* db_grid_reduce(DbGrid *db_grid,
                                       int iptr,
                                       int *margin,
                                       int *limmin,
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
GSTLEARN_EXPORT double bench_distance(const Db *db, int iech1, int iech2);
GSTLEARN_EXPORT double cylinder_radius(const Db *db,
                                       int iech1,
                                       int iech2,
                                       const VectorDouble &codir);
GSTLEARN_EXPORT double db_grid_maille(Db *db);
GSTLEARN_EXPORT int point_to_grid(const DbGrid *db,
                                  double *coor,
                                  int flag_expand,
                                  int *indg);
GSTLEARN_EXPORT int point_to_bench(const DbGrid *db,
                                   double *coor,
                                   int flag_outside,
                                   int *indb);
GSTLEARN_EXPORT void grid_to_point(const DbGrid *db,
                                   int *indg,
                                   double *percent,
                                   double *coor);
GSTLEARN_EXPORT int index_point_to_grid(const Db *db,
                                        int iech,
                                        int flag_expand,
                                        const DbGrid *dbout,
                                        double *coor);
GSTLEARN_EXPORT int point_to_point(Db *db, double *coor);
GSTLEARN_EXPORT int point_inside_grid(Db *db, int iech, const DbGrid *dbgrid);
GSTLEARN_EXPORT int migrate_grid_to_coor(const DbGrid *db_grid,
                                         int iv_grid,
                                         const VectorVectorDouble& coords,
                                         VectorDouble& tab);
GSTLEARN_EXPORT int expand_point_to_coor(const Db *db1,
                                         int iatt,
                                         const VectorVectorDouble& coords,
                                         VectorDouble& tab);
GSTLEARN_EXPORT int expand_point_to_grid(Db *db_point,
                                         DbGrid *db_grid,
                                         int iatt,
                                         int iatt_time,
                                         int iatt_angle,
                                         int iatt_scaleu,
                                         int iatt_scalev,
                                         int iatt_scalew,
                                         int flag_index,
                                         int ldmax,
                                         const VectorDouble &dmax,
                                         VectorDouble &tab);

GSTLEARN_EXPORT int interpolate_variable_to_point(DbGrid *db_grid,
                                                  int iatt,
                                                  int np,
                                                  double *xp,
                                                  double *yp,
                                                  double *zp,
                                                  double *tab);
GSTLEARN_EXPORT int points_to_block(Db *dbpoint,
                                    DbGrid *dbgrid,
                                    int option,
                                    int flag_size,
                                    int iatt_time,
                                    int iatt_size,
                                    int iatt_angle,
                                    int iatt_scaleu,
                                    int iatt_scalev,
                                    int iatt_scalew);
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
GSTLEARN_EXPORT int manage_external_info(int mode,
                                         const ELoc& locatorType,
                                         Db *dbin,
                                         Db *dbout,
                                         int *istart);
GSTLEARN_EXPORT int manage_nostat_info(int mode,
                                       Model *model,
                                       Db *dbin,
                                       Db *dbout);
GSTLEARN_EXPORT int db_locate_in_grid(DbGrid *dbgrid, double *coor);
GSTLEARN_EXPORT void db_monostat(Db *db,
                                 int ivar,
                                 double *wtot,
                                 double *mean,
                                 double *var,
                                 double *mini,
                                 double *maxi);
GSTLEARN_EXPORT int db_normalize(Db *db,
                                 const char *oper,
                                 int ncol,
                                 int *cols,
                                 double center,
                                 double stdv);
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
                                 int *ind1,
                                 int *ind2,
                                 int ncol,
                                 int *cols);
GSTLEARN_EXPORT int db_grid_copy_dilate(DbGrid *db1,
                                        int iatt1,
                                        DbGrid *db2,
                                        int iatt2,
                                        int mode,
                                        int *nshift);
GSTLEARN_EXPORT int db_proportion(Db *db,
                                  DbGrid *dbgrid,
                                  int nfac1max,
                                  int nfac2max,
                                  int *nclout);
GSTLEARN_EXPORT int db_merge(Db *db, int ncol, int *cols);
GSTLEARN_EXPORT int db_count_defined(Db *db, int icol);

GSTLEARN_EXPORT int db_prop_read(DbGrid *db, int ix, int iy, double *props);
GSTLEARN_EXPORT int db_prop_write(DbGrid *db, int ix, int iy, double *props);
GSTLEARN_EXPORT int db_resind(Db *db, int ivar, int ncut, double *zcut);
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
GSTLEARN_EXPORT int db_compositional_transform(Db *db,
                                               int verbose,
                                               int mode,
                                               int type,
                                               int number,
                                               int *iatt_in,
                                               int *iatt_out,
                                               int *numout);
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
                                  int flag_add_rank = 1);
GSTLEARN_EXPORT int db_smooth_vpc(DbGrid *db, int width, double range);
GSTLEARN_EXPORT double* db_grid_sampling(DbGrid *dbgrid,
                                         double *x1,
                                         double *x2,
                                         int ndisc,
                                         int ncut,
                                         double *cuts,
                                         int *nval_ret);
GSTLEARN_EXPORT int db_grid2point_sampling(DbGrid *dbgrid,
                                           int nvar,
                                           int *vars,
                                           int *npacks,
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

/******************************************/
/* Prototyping the functions in spatial.c */
/******************************************/

GSTLEARN_EXPORT int cgi(Db *db,
                        int ivar,
                        double *center,
                        double *mvalue,
                        double *mvector,
                        double *inertia,
                        double *wztot);
GSTLEARN_EXPORT int spatial(Db *db,
                            double *totab,
                            double *parea,
                            double *eqarea);

/******************************************/
/* Prototyping the functions in convert.c */
/******************************************/

GSTLEARN_EXPORT DbGrid* db_grid_read_f2g(const char *filename, int verbose = 0);
GSTLEARN_EXPORT int db_grid_write_zycor(const char *filename, DbGrid *db, int icol);
GSTLEARN_EXPORT DbGrid* db_grid_read_zycor(const char* filename, int verbose=0);
GSTLEARN_EXPORT int db_grid_write_arcgis(const char *filename, DbGrid *db, int icol);
GSTLEARN_EXPORT int db_grid_write_XYZ(const char *filename, DbGrid *db, int icol);
GSTLEARN_EXPORT int db_write_vtk(const char *filename,
                                 DbGrid *db,
                                 const VectorInt &cols);
GSTLEARN_EXPORT int db_grid_write_bmp(const char *filename,
                                      DbGrid *db,
                                      int icol,
                                      int nsamplex = 1,
                                      int nsampley = 1,
                                      int nmult = 1,
                                      int ncolors = 1,
                                      int flag_low = 1,
                                      int flag_high = 1,
                                      double valmin = TEST,
                                      double valmax = TEST,
                                      int *red = nullptr,
                                      int *green = nullptr,
                                      int *blue = nullptr,
                                      int mask_red = 0,
                                      int mask_green = 0,
                                      int mask_blue = 0,
                                      int ffff_red = 232,
                                      int ffff_green = 232,
                                      int ffff_blue = 0,
                                      int low_red = 255,
                                      int low_green = 255,
                                      int low_blue = 255,
                                      int high_red = 255,
                                      int high_green = 0,
                                      int highblue = 0);
GSTLEARN_EXPORT DbGrid* db_grid_read_bmp(const char* filename, int verbose=0);
GSTLEARN_EXPORT int db_grid_write_irap(const char *filename,
                                       DbGrid *db,
                                       int icol,
                                       int nsamplex = 1,
                                       int nsampley = 1);
GSTLEARN_EXPORT int db_grid_write_ifpen(const char *filename,
                                       DbGrid *db,
                                       int ncol,
                                       int *icols);
GSTLEARN_EXPORT DbGrid* db_grid_read_ifpen(const char* filename, int verbose=0);
GSTLEARN_EXPORT int db_grid_write_eclipse(const char *filename,
                                          DbGrid *db,
                                          int icol);
GSTLEARN_EXPORT Db* db_well_read_las(const char *filename,
                                     double xwell,
                                     double ywell,
                                     double cwell,
                                     int verbose = 0);
GSTLEARN_EXPORT int csv_table_read(const String &filename,
                                   const CSVformat& csvfmt,
                                   int verbose,
                                   int ncol_max,
                                   int nrow_max,
                                   int *ncol_arg,
                                   int *nrow_arg,
                                   VectorString &names,
                                   VectorDouble &tab);

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
                                         VectorInt ndisc);
GSTLEARN_EXPORT void krige_lhs_print(int nech,
                                     int neq,
                                     int nred,
                                     int *flag,
                                     double *lhs);
GSTLEARN_EXPORT void krige_rhs_print(int nvar,
                                     int nech,
                                     int neq,
                                     int nred,
                                     int *flag,
                                     double *rhs);
GSTLEARN_EXPORT void krige_dual_print(int nech,
                                      int neq,
                                      int nred,
                                      int *flag,
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
                                   int nsize1,
                                   int *ranks1,
                                   int nsize2,
                                   int *ranks2,
                                   bool flag_std,
                                   int verbose);
GSTLEARN_EXPORT int global_transitive(DbGrid *dbgrid,
                                      Model *model,
                                      int flag_verbose,
                                      int flag_regular,
                                      int ndisc,
                                      double *zest,
                                      double *cve,
                                      double *cvtrans);
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
                               double *ref_var,
                               int ref_radius,
                               int neigh_ver,
                               int neigh_hor,
                               int flag_sym,
                               Model *model,
                               double nugget,
                               int nfeq,
                               int dbg_ix,
                               int dbg_iy);
GSTLEARN_EXPORT int sampling_f(Db *db,
                               Model *model,
                               double beta,
                               int method1,
                               int nsz1_max,
                               int nsize1,
                               int *ranks1,
                               int method2,
                               int nsz2_max,
                               int nsize2,
                               int *ranks2,
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

GSTLEARN_EXPORT int get_rank_from_propdef(PropDef *propdef, int ipgs, int igrf);
GSTLEARN_EXPORT void check_mandatory_attribute(const char *method,
                                               Db *db,
                                               const ELoc& locatorType);
GSTLEARN_EXPORT int simcond(Db *dbin,
                            Db *dbout,
                            Model *model,
                            int seed,
                            int nbsimu,
                            int nbtuba,
                            int nboot,
                            int niter,
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
GSTLEARN_EXPORT int simtub_constraints(Db *dbin,
                                       Db *dbout,
                                       Model *model,
                                       ANeigh *neigh,
                                       int seed,
                                       int nbtuba,
                                       int nbsimu,
                                       int nbquant,
                                       int niter_max,
                                       VectorInt &cols,
                                       int (*func_valid)(int flag_grid,
                                                         int ndim,
                                                         int nech,
                                                         int *nx,
                                                         double *dx,
                                                         double *x0,
                                                         double nonval,
                                                         double percent,
                                                         VectorDouble &tab));
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

/******************************************/
/* Prototyping the functions in simpart.c */
/******************************************/

GSTLEARN_EXPORT SubPlanes* poisson_manage_planes(int mode,
                                                 int np,
                                                 SubPlanes *splanes);
GSTLEARN_EXPORT int poisson_generate_planes(DbGrid *dbgrid, SubPlanes *splanes);

/*****************************************/
/* Prototyping the functions in thresh.c */
/*****************************************/

GSTLEARN_EXPORT Rule* rule_free(const Rule *rule);
GSTLEARN_EXPORT Model* model_rule_combine(const Model *model1,
                                          const Model *model2,
                                          const Rule *rule);
GSTLEARN_EXPORT int rule_thresh_define_shadow(PropDef *propdef,
                                              Db *dbin,
                                              const RuleShadow *rule,
                                              int facies,
                                              int iech,
                                              int isimu,
                                              int nbsimu,
                                              double *t1min,
                                              double *t1max,
                                              double *t2min,
                                              double *t2max,
                                              double *dsup,
                                              double *down);
GSTLEARN_EXPORT int rule_thresh_define(PropDef *propdef,
                                       Db *dbin,
                                       const Rule *rule,
                                       int facies,
                                       int iech,
                                       int isimu,
                                       int nbsimu,
                                       int flag_check,
                                       double *t1min,
                                       double *t1max,
                                       double *t2min,
                                       double *t2max);
GSTLEARN_EXPORT int db_rule_shadow(Db *db,
                                   Db *dbprop,
                                   RuleShadow *rule,
                                   Model *model1,
                                   const VectorDouble &propcst,
                                   int flag_stat,
                                   int nfacies);
GSTLEARN_EXPORT int db_bounds_shadow(Db *db,
                                     Db *dbprop,
                                     RuleShadow *rule,
                                     Model *model,
                                     const VectorDouble &propcst,
                                     int flag_stat,
                                     int nfacies);
GSTLEARN_EXPORT void proportion_rule_process(PropDef *propdef,
                                             const EProcessOper &mode);
GSTLEARN_EXPORT PropDef* proportion_manage(int mode,
                                           int flag_facies,
                                           int flag_stat,
                                           int ngrf1,
                                           int ngrf2,
                                           int nfac1,
                                           int nfac2,
                                           Db *db,
                                           const Db *dbprop,
                                           const VectorDouble &propcst,
                                           PropDef *proploc);
GSTLEARN_EXPORT void propdef_reset(PropDef *propdef);
GSTLEARN_EXPORT void proportion_print(PropDef *propdef);

/******************************************/
/* Prototyping the functions in seismic.c */
/******************************************/

GSTLEARN_EXPORT int seismic_estimate_XZ(DbGrid *db,
                                        Model *model,
                                        int nbench,
                                        int nv2max,
                                        int flag_ks,
                                        int flag_std,
                                        int flag_sort,
                                        int flag_stat);
GSTLEARN_EXPORT int seismic_simulate_XZ(DbGrid *db,
                                        Model *model,
                                        int nbench,
                                        int nv2max,
                                        int nbsimu,
                                        int seed,
                                        int flag_ks,
                                        int flag_sort,
                                        int flag_stat);
GSTLEARN_EXPORT int seismic_z2t_grid(int verbose,
                                     DbGrid *db_z,
                                     int iptr_v,
                                     int *nx,
                                     double *x0,
                                     double *dx);
GSTLEARN_EXPORT int seismic_t2z_grid(int verbose,
                                     DbGrid *db_t,
                                     int iptr_v,
                                     int *nx,
                                     double *x0,
                                     double *dx);
GSTLEARN_EXPORT int seismic_z2t_convert(DbGrid *db_z, int iptr_v, DbGrid *db_t);
GSTLEARN_EXPORT int seismic_t2z_convert(DbGrid *db_t, int iptr_v, DbGrid *db_z);
GSTLEARN_EXPORT int seismic_operate(DbGrid *db, int oper);
GSTLEARN_EXPORT int seismic_convolve(DbGrid *db,
                                     int flag_operate,
                                     int flag_contrast,
                                     int type,
                                     int ntw,
                                     int option,
                                     int tindex,
                                     double fpeak,
                                     double period,
                                     double amplitude,
                                     double distort,
                                     double val_before,
                                     double val_middle,
                                     double val_after,
                                     double *wavelet);

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

/*******************************************/
/* Prototyping the functions in variopgs.c */
/*******************************************/

GSTLEARN_EXPORT Vario_Order* vario_order_manage(int mode,
                                                int flag_dist,
                                                int size_aux,
                                                Vario_Order *vorder);
GSTLEARN_EXPORT int vario_order_add(Vario_Order *vorder,
                                    int iech,
                                    int jech,
                                    void *aux_iech,
                                    void *aux_jech,
                                    int ipas,
                                    int idir,
                                    double dist);
GSTLEARN_EXPORT Vario_Order* vario_order_final(Vario_Order *vorder, int *npair);
GSTLEARN_EXPORT void vario_order_print(Vario_Order *vorder,
                                       int idir_target,
                                       int ilag_target,
                                       int verbose);
GSTLEARN_EXPORT void vario_order_get_bounds(Vario_Order *vorder,
                                            int idir,
                                            int ipas,
                                            int *ifirst,
                                            int *ilast);
GSTLEARN_EXPORT void vario_order_get_indices(Vario_Order *vorder,
                                             int ipair,
                                             int *iech,
                                             int *jech,
                                             double *dist);
GSTLEARN_EXPORT void vario_order_get_auxiliary(Vario_Order *vorder,
                                               int ipair,
                                               char *aux_iech,
                                               char *aux_jech);

/******************************************/
/* Prototyping the functions in mlayers.c */
/******************************************/
GSTLEARN_EXPORT int variogram_mlayers(Db *db,
                                      int *seltab,
                                      Vario *vario,
                                      Vario_Order *vorder);
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
GSTLEARN_EXPORT int multilayers_kriging(Db *dbin,
                                        DbGrid *dbout,
                                        Model *model,
                                        ANeigh *neigh,
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
                                        double *prior_mean,
                                        double *prior_vars,
                                        int irefd,
                                        int ireft,
                                        int irefb,
                                        int verbose);
GSTLEARN_EXPORT int multilayers_get_prior(Db *dbin,
                                          DbGrid *dbout,
                                          Model *model,
                                          int flag_same,
                                          int flag_vel,
                                          int flag_ext,
                                          int irf_rank,
                                          int match_time,
                                          int irefd,
                                          int ireft,
                                          int irefb,
                                          int verbose,
                                          int *npar,
                                          double **mean,
                                          double **vars);

/*******************************************/
/* Prototyping the functions in delaunay.c */
/*******************************************/
GSTLEARN_EXPORT double* get_db_extension(Db *dbin, Db *dbout, int *nout);
GSTLEARN_EXPORT double* extend_grid(DbGrid *db, const double *gext, int *nout);
GSTLEARN_EXPORT double* extend_point(Db *db, const double *gext, int *nout);
GSTLEARN_EXPORT int MSS(int idim, int ipol, int icas, int icorn, int icoor);
GSTLEARN_EXPORT int meshes_2D_write(const char *file_name,
                                    const char *obj_name,
                                    int verbose,
                                    int ndim,
                                    int ncode,
                                    int ntri,
                                    int npoints,
                                    const VectorInt& ntcode,
                                    const VectorInt& triangles,
                                    const VectorDouble& points);
GSTLEARN_EXPORT AMesh* meshes_turbo_1D_grid_build(DbGrid *dbgrid);
GSTLEARN_EXPORT AMesh* meshes_turbo_2D_grid_build(DbGrid *dbgrid);
GSTLEARN_EXPORT AMesh* meshes_turbo_3D_grid_build(DbGrid *dbgrid);

GSTLEARN_EXPORT void mesh_stats(int ndim,
                                int ncorner,
                                int nmesh,
                                int *meshes,
                                double *points);

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
                                   int iptr_pinch,
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
                                 int gibbs_nburn,
                                 int gibbs_niter,
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
GSTLEARN_EXPORT int spde_build_stdev(double *vcur);
GSTLEARN_EXPORT int spde_f(Db *dbin,
                           Db *dbout,
                           Model *model,
                           const VectorDouble &gext,
                           SPDE_Option &s_option,
                           int mesh_dbin,
                           int mesh_dbout,
                           int seed,
                           int nbsimu,
                           int gibbs_nburn,
                           int gibbs_niter,
                           int ngibbs_int,
                           int flag_est,
                           int flag_std,
                           int flag_gibbs,
                           int flag_modif,
                           int verbose);
#ifndef SWIG
GSTLEARN_EXPORT int spde_eval(int nblin,
                              double *blin,
                              cs *S,
                              const VectorDouble &Lambda,
                              const VectorDouble &TildeC,
                              double power,
                              double *x,
                              double *y);
#endif
GSTLEARN_EXPORT void spde_external_mesh_define(int icov0, AMesh *mesh);
GSTLEARN_EXPORT void spde_external_mesh_undefine(int icov0);
#ifndef SWIG
GSTLEARN_EXPORT int spde_external_copy(SPDE_Matelem &matelem, int icov0);
GSTLEARN_EXPORT cs* spde_external_A_define(int icov0, cs *A);
GSTLEARN_EXPORT cs* spde_external_Q_define(int icov0, cs *Q);
GSTLEARN_EXPORT cs* spde_external_A_undefine(int icov0);
GSTLEARN_EXPORT cs* spde_external_Q_undefine(int icov0);
#endif
GSTLEARN_EXPORT int kriging2D_spde(Db *dbin,
                                   Model *model,
                                   SPDE_Option &s_option,
                                   int verbose,
                                   int *nmeshes,
                                   int *nvertex,
                                   VectorInt& meshes,
                                   VectorDouble& points);
#ifndef SWIG
GSTLEARN_EXPORT cs* db_mesh_neigh(const Db *db,
                                  AMesh *amesh,
                                  double radius,
                                  int flag_exact,
                                  int verbose,
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

GSTLEARN_EXPORT CTables* ct_tables_manage(int mode,
                                          int verbose,
                                          int flag_cumul,
                                          int nconf,
                                          int ndisc,
                                          double cmin,
                                          double cmax,
                                          CTables *ctables_old);
GSTLEARN_EXPORT void ct_tables_print(CTables *ctables, int flag_print);
GSTLEARN_EXPORT int ct_tableone_covrank(const CTables *ctables,
                                        double cova,
                                        double *cround);
GSTLEARN_EXPORT int ct_tableone_getrank_from_proba(CTables *ctables,
                                                   double gaussian);
GSTLEARN_EXPORT double ct_tableone_calculate(CTables *ctables,
                                             int iconf0,
                                             double *lows,
                                             double *ups);
GSTLEARN_EXPORT double ct_tableone_calculate_by_rank(CTables *ctables,
                                                     int iconf0,
                                                     double *rklows,
                                                     double *rkups);
GSTLEARN_EXPORT double ct_INTRES2(CTables *ctables,
                                  int iconf0,
                                  int idisc0,
                                  int jdisc0);
GSTLEARN_EXPORT double ct_INTRES3(CTables *ctables,
                                  int iconf0,
                                  int idisc0,
                                  int jdisc0,
                                  int kdisc0);

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


/***************************/
/* Sparse matrix inversion */
/***************************/
GSTLEARN_EXPORT int sparseinv(int n, int *Lp, int *Li, double *Lx, double *d, int *Up,
                              int *Uj, double *Ux, int *Zp, int *Zi, double *Zx, double *z,
                              int *Zdiagp, int *Lmunch);
