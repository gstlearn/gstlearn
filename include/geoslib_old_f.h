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

// Enums
#include "Covariances/ECov.hpp"
#include "Covariances/ECalcMember.hpp"
#include "Basic/EJustify.hpp"
#include "Db/ELoc.hpp"
#include "LithoRule/EProcessOper.hpp"
#include "Model/EConsElem.hpp"
#include "Model/EConsType.hpp"
#include "Anamorphosis/EAnam.hpp"
#include "Drifts/EDrift.hpp"
#include "Neigh/ENeigh.hpp"
#include "Variogram/ECalcVario.hpp"

// References
#include "Covariances/CovCalcMode.hpp"
#include "Model/Constraints.hpp"
#include "Model/Option_AutoFit.hpp"

class Anam;
class AnamDiscreteDD;
class AnamDiscreteIR;
class AnamEmpirical;
class AnamHermite;
class PropDef;
class Rule;
class RuleShadow;
class MeshEStandard;
class CovInternal;
class Db;
class Model;
class Vario;
class VarioParam;
class Neigh;
class Polygons;
class PCA;
class Grid;

class cs;
class QChol;
class tetgenio;

/**********************************************/
/* Prototyping the functions in acknowledge.c */
/**********************************************/
GSTLEARN_EXPORT void acknowledge_gstlearn(void);
GSTLEARN_EXPORT void inquire_gstlearn(char **release, char **date);

/******************************************/
/* Prototyping the functions in license.c */
/******************************************/
GSTLEARN_EXPORT int register_license_file(const char *file_name,
                                          const char *target_name);

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

GSTLEARN_EXPORT int opt_mauto_add_constraints(Option_AutoFit &mauto,
                                              double constantSill);
GSTLEARN_EXPORT int opt_mauto_add_unit_constraints(Option_AutoFit &mauto);
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

GSTLEARN_EXPORT double ut_deg2rad(double angle);
GSTLEARN_EXPORT double ut_rad2deg(double angle);
GSTLEARN_EXPORT int get_mirror_sample(int nx, int ix);
GSTLEARN_EXPORT void get_matrix(const char *title,
                                int flag_sym,
                                int flag_def,
                                int nx,
                                int ny,
                                double valmin,
                                double valmax,
                                double *tab);
GSTLEARN_EXPORT void get_rotation(const char *title,
                                  int flag_def,
                                  int ndim,
                                  double *rot);
GSTLEARN_EXPORT void ut_sort_double(int safe, int nech, int *ind, double *tab);
GSTLEARN_EXPORT void ut_sort_int(int safe, int nech, int *ind, int *tab);
GSTLEARN_EXPORT void ut_tab_unique(int ntab, double *tab, int *neff);
GSTLEARN_EXPORT void ut_statistics(int nech,
                                   double *tab,
                                   double *sel,
                                   double *wgt,
                                   int *nval,
                                   double *mini,
                                   double *maxi,
                                   double *delta,
                                   double *mean,
                                   double *stdv);
GSTLEARN_EXPORT double ut_cnp(int n, int k);
GSTLEARN_EXPORT int* ut_combinations(int n, int maxk, int *ncomb);
GSTLEARN_EXPORT int* ut_split_into_two(int ncolor,
                                       int flag_half,
                                       int verbose,
                                       int *nposs);
GSTLEARN_EXPORT double* ut_pascal(int ndim);
GSTLEARN_EXPORT double ut_median(double *tab, int ntab);
GSTLEARN_EXPORT void rgb2num(int r, int g, int b, int a, unsigned char *value);
GSTLEARN_EXPORT void num2rgb(unsigned char value,
                             int *r,
                             int *g,
                             int *b,
                             int *a);
GSTLEARN_EXPORT void ut_stats_mima(int nech,
                                   double *tab,
                                   double *sel,
                                   int *nvalid,
                                   double *mini,
                                   double *maxi);
GSTLEARN_EXPORT void ut_stats_mima_print(const char *title,
                                         int nech,
                                         double *tab,
                                         double *sel);
GSTLEARN_EXPORT void ut_facies_statistics(int nech,
                                          double *tab,
                                          double *sel,
                                          int *nval,
                                          int *mini,
                                          int *maxi);
GSTLEARN_EXPORT void ut_classify(int nech,
                                 double *tab,
                                 double *sel,
                                 int nclass,
                                 double start,
                                 double pas,
                                 int *nmask,
                                 int *ntest,
                                 int *nout,
                                 int *classe);
GSTLEARN_EXPORT void ut_normalize(int ntab, double *tab);
GSTLEARN_EXPORT void ut_rotation_sincos(double angle,
                                        double *cosa,
                                        double *sina);
GSTLEARN_EXPORT void ut_rotation_matrix_2D(double angle, double *rot);
GSTLEARN_EXPORT void ut_rotation_matrix_3D(double alpha,
                                           double beta,
                                           double gamma,
                                           double *rot);
GSTLEARN_EXPORT void ut_rotation_matrix(int ndim,
                                        const double *angles,
                                        double *rot);
GSTLEARN_EXPORT VectorDouble ut_rotation_matrix_VD(int ndim,
                                                   const VectorDouble &angles);
GSTLEARN_EXPORT void ut_rotation_init(int ndim, double *rot);
GSTLEARN_EXPORT int ut_rotation_check(double *rot, int ndim);
GSTLEARN_EXPORT void ut_rotation_copy(int ndim,
                                      const double *rotin,
                                      double *rotout);
GSTLEARN_EXPORT void ut_rotation_direction(double ct,
                                           double st,
                                           double *a,
                                           double *codir);
GSTLEARN_EXPORT int ut_angles_from_rotation_matrix(const double *rot,
                                                   int ndim,
                                                   double *angles);
GSTLEARN_EXPORT void ut_angles_from_codir(int ndim,
                                          int ndir,
                                          const VectorDouble &codir,
                                          VectorDouble &angles);
GSTLEARN_EXPORT void ut_angles_to_codir(int ndim,
                                        int ndir,
                                        const VectorDouble &angles,
                                        VectorDouble &codr);
GSTLEARN_EXPORT double ut_merge_extension(int ndim,
                                          double *mini_in,
                                          double *maxi_in,
                                          double *mini_out,
                                          double *maxi_out);
GSTLEARN_EXPORT void debug_reset(void);
GSTLEARN_EXPORT void debug_print(void);
GSTLEARN_EXPORT void debug_index(int rank);
GSTLEARN_EXPORT void debug_reference(int rank);
GSTLEARN_EXPORT int is_debug_reference_defined(void);
GSTLEARN_EXPORT void debug_define(const char *name, int status);
GSTLEARN_EXPORT int debug_query(const char *name);
GSTLEARN_EXPORT int debug_force(void);
GSTLEARN_EXPORT void string_to_uppercase(char *string);
GSTLEARN_EXPORT void string_to_lowercase(char *string);
GSTLEARN_EXPORT int string_compare(int flag_case,
                                   const char *string1,
                                   const char *string2);
GSTLEARN_EXPORT void projec_query(int *actif);
GSTLEARN_EXPORT void projec_print(void);
GSTLEARN_EXPORT void projec_toggle(int mode);
GSTLEARN_EXPORT void variety_define(int flag_sphere, double radius);
GSTLEARN_EXPORT void variety_query(int *flag_sphere);
GSTLEARN_EXPORT void variety_print(void);
GSTLEARN_EXPORT void variety_toggle(int mode);
GSTLEARN_EXPORT void variety_get_characteristics(double *radius);
GSTLEARN_EXPORT double ut_factorial(int k);
GSTLEARN_EXPORT void ut_log_factorial(int nbpoly, double *factor);
GSTLEARN_EXPORT double golden_search(double (*func_evaluate)(double test,
                                                             void *user_data),
                                     void *user_data,
                                     double tolstop,
                                     double a0,
                                     double c0,
                                     double *testval,
                                     double *niter);
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
                                     const ELoc &ptype,
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
GSTLEARN_EXPORT PL_Dist* pldist_manage(int mode,
                                       PL_Dist *pldist_loc,
                                       int ndim,
                                       int nvert);
GSTLEARN_EXPORT double distance_point_to_segment(double x0,
                                                 double y0,
                                                 double x1,
                                                 double y1,
                                                 double x2,
                                                 double y2,
                                                 double *xd,
                                                 double *yd,
                                                 int *nint);
GSTLEARN_EXPORT void distance_point_to_polyline(double x0,
                                                double y0,
                                                int nvert,
                                                const double *x,
                                                const double *y,
                                                PL_Dist *pldist);
GSTLEARN_EXPORT double distance_along_polyline(PL_Dist *pldist1,
                                               PL_Dist *pldist2,
                                               double *xk,
                                               double *yl);
GSTLEARN_EXPORT double distance_points_to_polyline(double ap,
                                                   double al,
                                                   double x1,
                                                   double y1,
                                                   double x2,
                                                   double y2,
                                                   int nvert,
                                                   double *x,
                                                   double *y);
GSTLEARN_EXPORT int db_unfold_polyline(Db *db,
                                       int nvert,
                                       double *xl,
                                       double *yl);
GSTLEARN_EXPORT int db_fold_polyline(Db *dbin,
                                     Db *dbout,
                                     int ncol,
                                     int *cols,
                                     int nvert,
                                     double *xl,
                                     double *yl);
GSTLEARN_EXPORT double ut_geodetic_angular_distance(double long1,
                                                    double lat1,
                                                    double long2,
                                                    double lat2);
GSTLEARN_EXPORT void ut_geodetic_angles(double long1,
                                        double lat1,
                                        double long2,
                                        double lat2,
                                        double long3,
                                        double lat3,
                                        double *a,
                                        double *b,
                                        double *c,
                                        double *ga,
                                        double *gb,
                                        double *gc);
GSTLEARN_EXPORT double ut_geodetic_triangle_perimeter(double long1,
                                                      double lat1,
                                                      double long2,
                                                      double lat2,
                                                      double long3,
                                                      double lat3);
GSTLEARN_EXPORT double ut_geodetic_triangle_surface(double long1,
                                                    double lat1,
                                                    double long2,
                                                    double lat2,
                                                    double long3,
                                                    double lat3);
GSTLEARN_EXPORT int is_in_spherical_triangle(double *coor,
                                             double surface,
                                             double *pts1,
                                             double *pts2,
                                             double *pts3,
                                             double *wgts);
GSTLEARN_EXPORT int is_in_spherical_triangle_optimized(double *coo0,
                                                       double *ptsa,
                                                       double *ptsb,
                                                       double *ptsc,
                                                       double *wgts);
GSTLEARN_EXPORT double ut_distance(int ndim, double *tab1, double *tab2);
GSTLEARN_EXPORT void ut_distance_allocated(int ndim,
                                           double **tab1,
                                           double **tab2);
GSTLEARN_EXPORT int segment_intersect(double xd1,
                                      double yd1,
                                      double xe1,
                                      double ye1,
                                      double xd2,
                                      double yd2,
                                      double xe2,
                                      double ye2,
                                      double *xint,
                                      double *yint);
GSTLEARN_EXPORT int ut_chebychev_coeffs(double (*func)(double,
                                                       double,
                                                       int,
                                                       double*),
                                        Cheb_Elem *cheb_elem,
                                        int nblin,
                                        double *blin);
GSTLEARN_EXPORT int ut_chebychev_count(double (*func)(double,
                                                      double,
                                                      int,
                                                      double*),
                                       Cheb_Elem *cheb_elem,
                                       double x,
                                       int nblin,
                                       double *blin);
GSTLEARN_EXPORT void ut_vandercorput(int n,
                                     int flag_sym,
                                     int flag_rot,
                                     int *ntri_arg,
                                     double **coor_arg);
GSTLEARN_EXPORT int ut_icosphere(int n,
                                 int flag_rot,
                                 int *ntri_arg,
                                 double **coor_arg);
GSTLEARN_EXPORT void ut_shuffle_array(int nrow, int ncol, double *tab);
GSTLEARN_EXPORT int ut_is_legendre_defined(void);
GSTLEARN_EXPORT void define_legendre(double (*legendre_sphPlm)(int,
                                                               int,
                                                               double),
                                     double (*legendre_Pl)(int, double));
GSTLEARN_EXPORT double ut_legendre(int flag_norm, int n, double v);
GSTLEARN_EXPORT double ut_flegendre(int flag_norm, int n, int k0, double theta);
GSTLEARN_EXPORT int* ut_name_decode(const char *name,
                                    int ndim,
                                    int *nx,
                                    int verbose);
GSTLEARN_EXPORT double* ut_rank_cells(int ndim,
                                      int *nx,
                                      int *order,
                                      int verbose);

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
GSTLEARN_EXPORT void constant_reset(void);
GSTLEARN_EXPORT void constant_define(const char *name, double value);
GSTLEARN_EXPORT void constant_print(void);
GSTLEARN_EXPORT double constant_query(const char *name);
GSTLEARN_EXPORT void mem_error(int nbyte);

GSTLEARN_EXPORT void message_extern(const char *string);
GSTLEARN_EXPORT void exit_extern();

GSTLEARN_EXPORT void mes_process(const char *string, int ntot, int rank);
GSTLEARN_EXPORT void string_strip_blanks(char *string, int flag_lead);
GSTLEARN_EXPORT void string_strip_quotes(char *string);

#if defined(_WIN32) || defined(_WIN64)
GSTLEARN_EXPORT char * strsep(char **stringp, const char* delim);
#endif
GSTLEARN_EXPORT void tab_prints(const char *title,
                                int ncol,
                                const EJustify &justify,
                                const char *string);
GSTLEARN_EXPORT void tab_printg(const char *title,
                                int ncol,
                                const EJustify &justify,
                                double value);
GSTLEARN_EXPORT void tab_printd(const char *title,
                                int ncol,
                                const EJustify &justify,
                                double value);
GSTLEARN_EXPORT void tab_printi(const char *title,
                                int ncol,
                                const EJustify &justify,
                                int value);
GSTLEARN_EXPORT void tab_print_rowname(const char *string, int taille);
GSTLEARN_EXPORT void tab_print_rc(const char *title,
                                  int ncol,
                                  const EJustify &justify,
                                  int mode,
                                  int value);
GSTLEARN_EXPORT void encode_printg(char *string,
                                   int ntcar,
                                   int ntdec,
                                   double value);
GSTLEARN_EXPORT void print_current_line(void);
GSTLEARN_EXPORT void print_matrix(const char *title,
                                  int flag_limit,
                                  int byrow,
                                  int nx,
                                  int ny,
                                  const double *sel,
                                  const double *tab);
GSTLEARN_EXPORT void print_trimat(const char *title,
                                  int mode,
                                  int neq,
                                  const double *tl);
GSTLEARN_EXPORT void print_imatrix(const char *title,
                                   int flag_limit,
                                   int bycol,
                                   int nx,
                                   int ny,
                                   const double *sel,
                                   const int *tab);
GSTLEARN_EXPORT void print_vector(const char *title,
                                  int flag_limit,
                                  int ntab,
                                  const double *tab);
GSTLEARN_EXPORT void print_names(int nx, int *ranks, VectorString names);
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

GSTLEARN_EXPORT void matrix_constant_define(int keywrd, double value);
GSTLEARN_EXPORT double matrix_constant_query(int keywrd);
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
GSTLEARN_EXPORT void matrix_svd_inverse(int neq,
                                        double *s,
                                        double *u,
                                        double *v,
                                        double *tabout);
GSTLEARN_EXPORT double matrix_determinant(int neq, const double *b);
GSTLEARN_EXPORT int matrix_cofactor(int neq, double *a, double *b);
GSTLEARN_EXPORT double matrix_cholesky_determinant(int neq, double *tl);
GSTLEARN_EXPORT int matrix_eigen(const double *a,
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
GSTLEARN_EXPORT void matrix_transpose_in_place(int n1, int n2, double *v1);
GSTLEARN_EXPORT void matrix_int_transpose_in_place(int n1, int n2, int *v1);
GSTLEARN_EXPORT int matrix_solve(int mode,
                                 const double *a,
                                 const double *b,
                                 double *x,
                                 int neq,
                                 int nrhs,
                                 int *pivot);
GSTLEARN_EXPORT double matrix_norm(double *a, int neq);
GSTLEARN_EXPORT double matrix_normA(double *b, double *a, int neq, int subneq);
GSTLEARN_EXPORT double inner_product(const double *a, const double *b, int neq);
GSTLEARN_EXPORT void vector_product(double *a, double *b, double *v);
GSTLEARN_EXPORT void vector_translate(int ndim,
                                      double *a,
                                      double *v,
                                      double *b);
GSTLEARN_EXPORT int matrix_cholesky_decompose(const double *a,
                                              double *tl,
                                              int neq);
GSTLEARN_EXPORT void matrix_cholesky_product(int mode,
                                             int neq,
                                             int nrhs,
                                             double *tl,
                                             double *a,
                                             double *x);
GSTLEARN_EXPORT int matrix_cholesky_solve(int neq,
                                          double *tl,
                                          double *b,
                                          double *x);
GSTLEARN_EXPORT int matrix_cholesky_to_invert(int neq, double *tl, double *xl);
GSTLEARN_EXPORT void matrix_cholesky_invert(int neq, double *tl, double *xl);
GSTLEARN_EXPORT void matrix_cholesky_norme(int mode,
                                           int neq,
                                           double *tl,
                                           double *a,
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
GSTLEARN_EXPORT int is_matrix_null(int nrow, int ncol, double *a, int verbose);
GSTLEARN_EXPORT int is_matrix_symmetric(int neq, const double *a, int verbose);
GSTLEARN_EXPORT int is_matrix_correlation(int neq, double *a);
GSTLEARN_EXPORT int is_matrix_rotation(int neq, const double *a, int verbose);
GSTLEARN_EXPORT void matrix_produit_lu(int neq, double *al, double *a);
GSTLEARN_EXPORT VectorDouble matrix_produit_lu_VD(int neq, double *tl);
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
                                               double *tl,
                                               double *a);
GSTLEARN_EXPORT void matrix_tri2sq(int neq, double *tl, double *a);
GSTLEARN_EXPORT void matrix_square_to_triangle(int mode,
                                               int neq,
                                               double *a,
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
GSTLEARN_EXPORT void ascii_external_filename(const char *filein,
                                             int mode,
                                             char *filename);
GSTLEARN_EXPORT void ascii_filename(const char *type,
                                    int rank,
                                    int mode,
                                    char *filename);
GSTLEARN_EXPORT int ascii_anam_write(const char *file_name,
                                     const Anam *anam,
                                     int verbose,
                                     int flag_calcul);
GSTLEARN_EXPORT int ascii_frac_write(const char *file_name,
                                     Frac_Environ *frac,
                                     int verbose);
GSTLEARN_EXPORT Db* ascii_db_read(const char *file_name,
                                  int must_grid,
                                  int verbose);
GSTLEARN_EXPORT Vario* ascii_vario_read(const char *file_name, bool verbose);
GSTLEARN_EXPORT Neigh* ascii_neigh_read(const char *file_name, int verbose);
GSTLEARN_EXPORT Model* ascii_model_read(const char *file_name, int verbose);
GSTLEARN_EXPORT void ascii_simu_read(char *file_name,
                                     int verbose,
                                     int *nbsimu,
                                     int *nbtuba,
                                     int *seed);
GSTLEARN_EXPORT Rule* ascii_rule_read(const char *file_name, int verbose);
GSTLEARN_EXPORT Anam* ascii_anam_read(const char *file_name, int verbose);
GSTLEARN_EXPORT Frac_Environ* ascii_frac_read(const char *file_name,
                                              int verbose);
GSTLEARN_EXPORT int ascii_option_defined(const char *file_name,
                                         int verbose,
                                         const char *option_name,
                                         int type,
                                         void *answer);

/*****************************************/
/* Prototyping the functions in morpho.c */
/*****************************************/

GSTLEARN_EXPORT int fluid_propagation(Db *dbgrid,
                                      int verbose,
                                      int seed,
                                      int niter,
                                      int ind_facies,
                                      int ind_fluid,
                                      int ind_perm,
                                      int ind_poro,
                                      int nfacies,
                                      int nfluids,
                                      int *speeds,
                                      int flag_show,
                                      double number_max,
                                      double volume_max);
GSTLEARN_EXPORT int fluid_extract(Db *dbgrid,
                                  int verbose,
                                  int ind_date,
                                  int ind_facies,
                                  int ind_fluid,
                                  int ind_poro,
                                  int nfacies,
                                  int nfluids,
                                  int facies0,
                                  int fluid0,
                                  int ntime,
                                  double time0,
                                  double dtime,
                                  double *tab);
GSTLEARN_EXPORT int spill_point(Db *dbgrid,
                                int ind_height,
                                int ind_data,
                                int flag_up,
                                int flag_cross,
                                int flag_unknown,
                                int flag_verbose,
                                double hmax,
                                double *h,
                                double *th,
                                int *ix0,
                                int *iy0);

/****************************************/
/* Prototyping the functions in vario.c */
/****************************************/

GSTLEARN_EXPORT void vario_fix_codir(int ndim, VectorDouble &codir);
GSTLEARN_EXPORT PCA* pca_free(PCA *pca);
GSTLEARN_EXPORT PCA* pca_alloc(int nvar);
GSTLEARN_EXPORT int pca_z2f(Db *db, PCA *pca, int flag_norm, int flag_verbose);
GSTLEARN_EXPORT int pca_f2z(Db *db, PCA *pca, int flag_norm, int flag_verbose);

GSTLEARN_EXPORT Vario* variogram_delete(Vario *vario);
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
GSTLEARN_EXPORT void variogram_scale(Vario *vario, int idir);
GSTLEARN_EXPORT int variogram_get_lag(Vario *vario,
                                      int idir,
                                      double ps,
                                      double psmin,
                                      double *dist);
GSTLEARN_EXPORT ECalcVario vario_identify_calcul_type(const String &cov_name);
GSTLEARN_EXPORT void vardir_print(Vario *vario, int idir, int verbose);
GSTLEARN_EXPORT void vardir_copy(VarioParam *vario_in,
                                 int idir_in,
                                 VarioParam *vario_out,
                                 int idir_out);
GSTLEARN_EXPORT void variogram_trans_cut(Vario *vario, int nh, double ycut);
GSTLEARN_EXPORT int correlation_f(Db *db1,
                                  Db *db2,
                                  Db *dbgrid,
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
                                           Db *dbgrid,
                                           Vario *vario,
                                           Polygons *polygon);
GSTLEARN_EXPORT int regression_f(Db *db1,
                                 Db *db2,
                                 int flag_mode,
                                 int icol,
                                 int ncol,
                                 int *icols,
                                 Model *model,
                                 int flag_one,
                                 int flag_verbose,
                                 int *count,
                                 double *coeff,
                                 double *variance,
                                 double *varres,
                                 double *correl);
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
GSTLEARN_EXPORT int maf_compute(Db *db,
                                int opt_code,
                                double tolcode,
                                VectorDouble &codir,
                                double tolang,
                                double bench,
                                double cylrad,
                                double h0,
                                double dh,
                                int verbose,
                                PCA *pca);
GSTLEARN_EXPORT int pca_compute(Db *db, int verbose, PCA *pca);
GSTLEARN_EXPORT int variogram_y2z(Vario *vario, Anam *anam, Model *model);

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
GSTLEARN_EXPORT double model_calcul_cov_ij(Model *model,
                                           const CovCalcMode &mode,
                                           int ivar,
                                           int jvar,
                                           const VectorDouble &d1);
GSTLEARN_EXPORT double model_calcul_stdev(Model *model,
                                          Db *db1,
                                          int iech1,
                                          Db *db2,
                                          int iech2,
                                          int verbose,
                                          double factor);
GSTLEARN_EXPORT void model_calcul_drift(Model *model,
                                        const ECalcMember &member,
                                        const Db *db,
                                        int iech,
                                        double *drftab);
GSTLEARN_EXPORT void model_variance0(Model *model,
                                     Koption *koption,
                                     double *covtab,
                                     double *var0);
GSTLEARN_EXPORT void model_variance0_nostat(Model *model,
                                            Koption *koption,
                                            CovInternal *covint,
                                            double *covtab,
                                            double *var0);
GSTLEARN_EXPORT Model* model_free(Model *model);
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
GSTLEARN_EXPORT int model_add_drift(Model *model,
                                    const EDrift &type,
                                    int rank_fex);
GSTLEARN_EXPORT int model_anamorphosis_set_factor(Model *model, int iclass);
GSTLEARN_EXPORT int model_sample(Vario *vario,
                                 Model *model,
                                 int flag_norm,
                                 int flag_cov);
GSTLEARN_EXPORT void model_calcul_cov(CovInternal *covint,
                                             Model *model,
                                             const CovCalcMode &mode,
                                             int flag_init,
                                             double weight,
                                             VectorDouble d1,
                                             double *covtab);
GSTLEARN_EXPORT int model_fitting_sills(const Vario *vario,
                                        Model *model,
                                        const Option_AutoFit& mauto);
GSTLEARN_EXPORT int model_nfex(Model *model);
GSTLEARN_EXPORT int model_update_coreg(Model *model,
                                       double *aic,
                                       double *valpro,
                                       double *vecpro);
GSTLEARN_EXPORT int model_evaluate(Model *model,
                                   int ivar,
                                   int jvar,
                                   int rank_sel,
                                   int flag_norm,
                                   int flag_cov,
                                   int nugget_opt,
                                   int nostd,
                                   int norder,
                                   const ECalcMember &member,
                                   int nh,
                                   VectorDouble &codir,
                                   double *h,
                                   double *g);
GSTLEARN_EXPORT int model_evaluate_nostat(Model *model,
                                          int ivar,
                                          int jvar,
                                          int rank_sel,
                                          int flag_norm,
                                          int flag_cov,
                                          int nugget_opt,
                                          int nostd,
                                          int norder,
                                          const ECalcMember &member,
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
                               int flag_norm,
                               int flag_cov,
                               double *g);
GSTLEARN_EXPORT double model_cxx(Model *model,
                                 Db *db1,
                                 Db *db2,
                                 int ivar,
                                 int jvar,
                                 int seed,
                                 double epsdist);
GSTLEARN_EXPORT void model_covmat(Model *model,
                                  Db *db1,
                                  Db *db2,
                                  int ivar,
                                  int jvar,
                                  int flag_norm,
                                  int flag_cov,
                                  double *covmat);
GSTLEARN_EXPORT double* model_covmat_by_ranks(Model *model,
                                              Db *db1,
                                              int nsize1,
                                              const int *ranks1,
                                              Db *db2,
                                              int nsize2,
                                              const int *ranks2,
                                              int ivar0 = -1,
                                              int jvar0 = -1,
                                              int flag_norm = 0,
                                              int flag_cov = 1);
GSTLEARN_EXPORT cs* model_covmat_by_ranks_cs(Model *model,
                                             Db *db1,
                                             int nsize1,
                                             const int *ranks1,
                                             Db *db2,
                                             int nsize2,
                                             const int *ranks2,
                                             int ivar0,
                                             int jvar0,
                                             int flag_norm,
                                             int flag_cov);
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
                                        double **Gret);
GSTLEARN_EXPORT void model_drift_mat(Model *model,
                                     const ECalcMember &member,
                                     Db *db,
                                     double *drfmat);
GSTLEARN_EXPORT void model_drift_vector(Model *model,
                                        const ECalcMember &member,
                                        Db *db,
                                        int iech,
                                        double *vector);
GSTLEARN_EXPORT void model_vector(Model *model,
                                  Db *db1,
                                  Db *db2,
                                  int ivar,
                                  int jvar,
                                  int iech,
                                  int flag_norm,
                                  int flag_cov,
                                  double *vector);
GSTLEARN_EXPORT void model_vector_nostat(Model *model,
                                         Db *db,
                                         int ivar,
                                         int jvar,
                                         int iech,
                                         double *vector);
GSTLEARN_EXPORT void model_vector_multivar(Model *model,
                                           Db *db,
                                           int ivar,
                                           int iech,
                                           int flag_norm,
                                           int flag_cov,
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
GSTLEARN_EXPORT int model_is_drift_defined(Model *model, const EDrift &type0);
GSTLEARN_EXPORT Model* input_model(int ndim,
                                   int nvar,
                                   int order,
                                   int flag_sill,
                                   int flag_norm,
                                   Model *model_in);
GSTLEARN_EXPORT int model_dimension(Model *model);
GSTLEARN_EXPORT double model_get_field(Model *model);
GSTLEARN_EXPORT int model_extract_cova(Model *model,
                                       int icov,
                                       ECov *cov_type,
                                       int *flag_aniso,
                                       double *param,
                                       VectorDouble &sill,
                                       VectorDouble &aniso_rotmat,
                                       VectorDouble &aniso_ranges);
GSTLEARN_EXPORT void model_extract_properties(Model *model, double *tape_range);
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
GSTLEARN_EXPORT double model_maximum_distance(Model *model);
GSTLEARN_EXPORT int model_maximum_order(Model *model);
GSTLEARN_EXPORT double model_scale2range(const ECov &type,
                                         double scale,
                                         double param);
GSTLEARN_EXPORT double model_range2scale(const ECov &type,
                                         double range,
                                         double param);
GSTLEARN_EXPORT double cova_get_scale_factor(const ECov &type, double param);
GSTLEARN_EXPORT Model* model_combine(const Model *model1,
                                     const Model *model2,
                                     double r);
GSTLEARN_EXPORT int model_get_nonugget_cova(Model *model);
GSTLEARN_EXPORT int model_regularize(Model *model,
                                     Vario *vario,
                                     Db *db,
                                     int opt_norm,
                                     double nug_ratio);
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

/****************************************/
/* Prototyping the functions in neigh.c */
/****************************************/

GSTLEARN_EXPORT int neigh_start(Db *dbin, Neigh *neigh);
GSTLEARN_EXPORT void neigh_stop(void);
GSTLEARN_EXPORT int neigh_select(Db *dbin,
                                 Db *dbout,
                                 int iech_out,
                                 Neigh *neigh,
                                 int flag_simu,
                                 int flag_no_var_check,
                                 int *nech,
                                 int *rank);
GSTLEARN_EXPORT Neigh* neigh_free(Neigh *neigh);
GSTLEARN_EXPORT Neigh* neigh_init_bench(int ndim,
                                        int flag_xvalid,
                                        double width);
GSTLEARN_EXPORT Neigh* neigh_init_unique(int ndim);
GSTLEARN_EXPORT Neigh* neigh_init_image(int ndim,
                                        int flag_xvalid,
                                        int skip,
                                        const VectorInt &nbgh_image = VectorInt());
GSTLEARN_EXPORT Neigh* neigh_init(int ndim,
                                  const ENeigh& type,
                                  int flag_xvalid,
                                  int flag_sector,
                                  int flag_aniso,
                                  int flag_rotation,
                                  int flag_continuous,
                                  int nmini,
                                  int nmaxi,
                                  int nsect,
                                  int nsmax,
                                  int skip,
                                  double width,
                                  double radius,
                                  double dist_cont,
                                  const VectorDouble &nbgh_radius = VectorDouble(),
                                  const VectorDouble &nbgh_rotmat = VectorDouble(),
                                  const VectorInt &nbgh_image = VectorInt());
GSTLEARN_EXPORT void neigh_print(const Neigh *neigh);
GSTLEARN_EXPORT void neigh_echo(Db *dbin,
                                Neigh *neigh,
                                int *rank,
                                int nsel,
                                double *tab);
GSTLEARN_EXPORT int neigh_extract(Neigh *neigh,
                                  ENeigh *type,
                                  int *nmini,
                                  int *nmaxi,
                                  int *nsect,
                                  int *nsmax,
                                  int *skip,
                                  int *flag_sector,
                                  int *flag_aniso,
                                  int *flag_rotation,
                                  int *flag_continuous,
                                  double *width,
                                  double *radius,
                                  double *dist_cont,
                                  VectorDouble &nbgh_rotmat,
                                  VectorDouble &nbgh_radius,
                                  VectorInt &nbgh_image);
GSTLEARN_EXPORT int* neigh_calc(Db *dbin,
                                Model *model,
                                Neigh *neigh,
                                double *target,
                                int *nech_out);
GSTLEARN_EXPORT double neigh_continuous_variance(Neigh *neigh,
                                                 Db *db1,
                                                 int rank1,
                                                 Db *db2,
                                                 int rankZ);

/***************************************/
/* Prototyping the functions in anam.c */
/***************************************/

GSTLEARN_EXPORT double anam_y2z(Anam *anam, double y, int flag_bound);
GSTLEARN_EXPORT void anam_update_hermitian(AnamHermite *anam_hermite,
                                           double pymin,
                                           double pzmin,
                                           double pymax,
                                           double pzmax,
                                           double aymin,
                                           double azmin,
                                           double aymax,
                                           double azmax,
                                           double r,
                                           const VectorDouble &psi_hn);
GSTLEARN_EXPORT void anam_update_discrete_DD(AnamDiscreteDD *anam_discrete_DD,
                                             int ncut,
                                             double scoef,
                                             double mu,
                                             const VectorDouble &zcut,
                                             const VectorDouble &pcaz2f,
                                             const VectorDouble &pcaf2z,
                                             const VectorDouble &stats);
GSTLEARN_EXPORT void anam_update_empirical(AnamEmpirical *anam_empirical,
                                           int ndisc,
                                           double pymin,
                                           double pzmin,
                                           double pymax,
                                           double pzmax,
                                           double aymin,
                                           double azmin,
                                           double aymax,
                                           double azmax,
                                           double sigma2e,
                                           const VectorDouble &tdisc);
GSTLEARN_EXPORT void anam_update_discrete_IR(AnamDiscreteIR *anam_discrste_IR,
                                             int ncut,
                                             double s,
                                             const VectorDouble &zcut,
                                             const VectorDouble &stats);
GSTLEARN_EXPORT int anam_discrete_DD_z2factor(Anam *anam,
                                              Db *db,
                                              int iptr,
                                              int nfact,
                                              VectorInt ifacs);
GSTLEARN_EXPORT int anam_discrete_IR_z2factor(Anam *anam,
                                              Db *db,
                                              int iptr,
                                              int nfact,
                                              VectorInt ifacs);
GSTLEARN_EXPORT int anam_discrete_z2factor(Anam *anam,
                                           Db *db,
                                           int nfact,
                                           const VectorInt &ifacs);
GSTLEARN_EXPORT int anam_point_to_block(Anam *anam,
                                        int verbose,
                                        double cvv,
                                        double coeff,
                                        double mu);
GSTLEARN_EXPORT double ce_compute_Z2(double krigest,
                                     double krigstd,
                                     const VectorDouble &phis);
GSTLEARN_EXPORT int anam_factor2qt(Db *db,
                                   Anam *anam,
                                   int ncutmine,
                                   double *cutmine,
                                   double z_max,
                                   int flag_correct,
                                   int nb_est,
                                   int *cols_est,
                                   int nb_std,
                                   int *cols_std,
                                   int ncode,
                                   int *codes,
                                   int *ncut,
                                   int *qt_vars);
GSTLEARN_EXPORT void selectivity_interpolate(int verbose,
                                             double *zcut,
                                             int nclass,
                                             double *calest,
                                             int ncut,
                                             double *calcut);
GSTLEARN_EXPORT int anam_get_r(Anam *anam, double cvv, double mu, double *r);
GSTLEARN_EXPORT int anam_vario_z2y(Anam *anam, double cvv, Vario *vario);

GSTLEARN_EXPORT int uc_f(Db *db,
                         Anam *anam,
                         int att_est,
                         int att_var,
                         int ncutmine,
                         double *cutmine,
                         double proba,
                         double var_bloc,
                         int ncode,
                         int *codes,
                         int verbose,
                         int *qt_vars);
GSTLEARN_EXPORT int ce_f(Db *db,
                         Anam *anam,
                         int att_est,
                         int att_std,
                         int flag_est,
                         int flag_std,
                         int flag_OK,
                         int ncutmine,
                         double *cutmine,
                         double proba,
                         int ncode,
                         int *codes,
                         int verbose,
                         int nbsimu,
                         int *qt_vars);

/*************************************/
/* Prototyping the functions in db.c */
/*************************************/

GSTLEARN_EXPORT void grid_iterator_init(Grid *grid,
                                        const VectorInt &order = VectorInt());
GSTLEARN_EXPORT VectorInt grid_iterator_next(Grid *grid);

GSTLEARN_EXPORT int* db_indg_alloc(const Db *db);
GSTLEARN_EXPORT int* db_indg_free(int *indice);
GSTLEARN_EXPORT double* db_sample_free(double *tab);
GSTLEARN_EXPORT double* db_sample_alloc(const Db *db, const ELoc &locatorType);
GSTLEARN_EXPORT int db_sample_load(Db *db,
                                   const ELoc &locatorType,
                                   int iech,
                                   double *tab);
GSTLEARN_EXPORT double* db_vector_free(double *tab);
GSTLEARN_EXPORT double* db_vector_alloc(const Db *db);
GSTLEARN_EXPORT int db_selection_get(const Db *db, int item, double *tab);
GSTLEARN_EXPORT int db_vector_get(Db *db,
                                  const ELoc &locatorType,
                                  int locatorIndex,
                                  double *tab);
GSTLEARN_EXPORT int db_vector_put(Db *db,
                                  const ELoc &locatorType,
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
                                          const ELoc &locatorType,
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
                                             const ELoc &locatorType,
                                             int number,
                                             int r_tem,
                                             double valinit,
                                             int *iptr);
GSTLEARN_EXPORT void db_locators_correct(VectorString &strings,
                                         const VectorInt &current,
                                         int flag_locnew);
GSTLEARN_EXPORT int db_coorvec_put(Db *db, int idim, double *tab);
GSTLEARN_EXPORT int db_coorvec_get(const Db *db, int idim, double *tab);
GSTLEARN_EXPORT Db* db_delete(Db *db);
GSTLEARN_EXPORT int db_grid_match(Db *db1, Db *db2);
GSTLEARN_EXPORT int db_is_isotropic(Db *db, int iech, double *data);
GSTLEARN_EXPORT void db_grid_print(Db *db);

GSTLEARN_EXPORT Db* db_create_grid_multiple(Db *dbin,
                                            const VectorInt &nmult,
                                            int flag_add_rank);
GSTLEARN_EXPORT Db* db_create_grid_divider(Db *dbin,
                                           const VectorInt &nmult,
                                           int flag_add_rank);
GSTLEARN_EXPORT Db* db_create_grid_dilate(Db *dbin,
                                          int mode,
                                          const VectorInt &nshift,
                                          int flag_add_rank);
GSTLEARN_EXPORT Db* db_grid_sample(Db *dbin, const VectorInt &nmult);
GSTLEARN_EXPORT int db_grid_define_coordinates(Db *db);
GSTLEARN_EXPORT Db* db_create_from_target(double *target,
                                          int ndim,
                                          int flag_add_rank);
GSTLEARN_EXPORT void db_sample_print(Db *db,
                                     int iech,
                                     int flag_ndim,
                                     int flag_nvar,
                                     int flag_nerr);
GSTLEARN_EXPORT int db_center(Db *db, double *center);
GSTLEARN_EXPORT int db_extension(Db *db,
                                 double *mini,
                                 double *maxi,
                                 double *delta);
GSTLEARN_EXPORT int db_extension_rotated(Db *db,
                                         double *rotmat,
                                         double *mini_arg,
                                         double *maxi_arg,
                                         double *delta_arg);
GSTLEARN_EXPORT int db_attribute_range(const Db *db,
                                       int icol,
                                       double *mini,
                                       double *maxi,
                                       double *delta);
GSTLEARN_EXPORT int db_extension_diag(const Db *db, double *diag);
GSTLEARN_EXPORT double db_epsilon_distance(Db *db);
GSTLEARN_EXPORT int db_index_grid_to_sample(const Db *db, const int *indg);
GSTLEARN_EXPORT void db_index_sample_to_grid(const Db *db, int iech, int *indg);
GSTLEARN_EXPORT int db_index_sorted_in_grid(const Db *db, int iech, int *indg);
GSTLEARN_EXPORT int db_selref(int ndim,
                              int *nx,
                              int *ref,
                              double *tabin,
                              double *tabout);
GSTLEARN_EXPORT Db* db_extract(Db *db, int *ranks);
GSTLEARN_EXPORT Db* db_regularize(Db *db, Db *dbgrid, int flag_center);
GSTLEARN_EXPORT int compat_NDIM(Db *db1, Db *db2);
GSTLEARN_EXPORT double get_grid_value(Db *dbgrid,
                                      int iptr,
                                      int *indg,
                                      int ix,
                                      int iy,
                                      int iz);
GSTLEARN_EXPORT void set_grid_value(Db *dbgrid,
                                    int iptr,
                                    int *indg,
                                    int ix,
                                    int iy,
                                    int iz,
                                    double value);
GSTLEARN_EXPORT int get_LOCATOR_NITEM(const Db *db, const ELoc &locatorType);
GSTLEARN_EXPORT int exist_LOCATOR(Db *db, const ELoc &locatorType);
GSTLEARN_EXPORT double get_LOCATOR_ITEM(Db *db,
                                        const ELoc &locatorType,
                                        int locatorIndex,
                                        int iech);
GSTLEARN_EXPORT void set_LOCATOR_ITEM(Db *db,
                                      const ELoc &locatorType,
                                      int locatorIndex,
                                      int iech,
                                      double value);
GSTLEARN_EXPORT int db_get_rank_absolute_to_relative(Db *db, int iech0);
GSTLEARN_EXPORT int db_get_rank_relative_to_absolute(Db *db, int iech0);
GSTLEARN_EXPORT int is_grid(const Db *db, bool verbose = false);
GSTLEARN_EXPORT int is_grid_multiple(Db *db1, Db *db2);
GSTLEARN_EXPORT int db_grid_copy_params(Db *dbin, int mode, Db *dbout);
GSTLEARN_EXPORT Db* db_grid_reduce(Db *db_grid,
                                   int iptr,
                                   int *margin,
                                   int *limmin,
                                   int flag_sel,
                                   int flag_copy,
                                   int verbose,
                                   double vmin,
                                   double vmax);
GSTLEARN_EXPORT double distance_inter(Db *db1,
                                      Db *db2,
                                      int iech1,
                                      int iech2,
                                      double *dist_vect);
GSTLEARN_EXPORT double distance_intra(const Db *db,
                                      int iech1,
                                      int iech2,
                                      double *dist_vect);
GSTLEARN_EXPORT double distance_grid(Db *db,
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
GSTLEARN_EXPORT int point_to_grid(const Db *db,
                                  double *coor,
                                  int flag_expand,
                                  int *indg);
GSTLEARN_EXPORT int point_to_bench(const Db *db,
                                   double *coor,
                                   int flag_outside,
                                   int *indb);
GSTLEARN_EXPORT void grid_to_point(const Db *db,
                                   int *indg,
                                   double *percent,
                                   double *coor);
GSTLEARN_EXPORT int index_point_to_grid(const Db *db,
                                        int iech,
                                        int flag_expand,
                                        const Db *dbout,
                                        double *coor);
GSTLEARN_EXPORT int point_to_point(Db *db, double *coor);
GSTLEARN_EXPORT int point_inside_grid(Db *db, int iech, Db *dbgrid);
GSTLEARN_EXPORT int migrate_grid_to_coor(const Db *db_grid,
                                         int iv_grid,
                                         int np,
                                         double *xp,
                                         double *yp,
                                         double *zp,
                                         double *tab);
GSTLEARN_EXPORT int expand_point_to_coor(const Db *db1,
                                         int iatt,
                                         int np,
                                         double *xp,
                                         double *yp,
                                         double *zp,
                                         double *tab);
GSTLEARN_EXPORT int expand_point_to_grid(Db *db_point,
                                         Db *db_grid,
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
GSTLEARN_EXPORT int db_center_point_to_grid(Db *db_point,
                                            Db *db_grid,
                                            double eps_random);
GSTLEARN_EXPORT int interpolate_variable_to_point(Db *db_grid,
                                                  int iatt,
                                                  int np,
                                                  double *xp,
                                                  double *yp,
                                                  double *zp,
                                                  double *tab);
GSTLEARN_EXPORT int points_to_block(Db *dbpoint,
                                    Db *dbgrid,
                                    int option,
                                    int flag_size,
                                    int iatt_time,
                                    int iatt_size,
                                    int iatt_angle,
                                    int iatt_scaleu,
                                    int iatt_scalev,
                                    int iatt_scalew);
GSTLEARN_EXPORT int db_gradient_components(Db *dbgrid);
GSTLEARN_EXPORT int db_streamline(Db *dbgrid,
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
                                         const ELoc &locatorType,
                                         Db *dbin,
                                         Db *dbout,
                                         int *istart);
GSTLEARN_EXPORT int manage_nostat_info(int mode,
                                       Model *model,
                                       Db *dbin,
                                       Db *dbout);
GSTLEARN_EXPORT int db_locate_in_grid(Db *dbgrid, double *coor);
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
                            Db *db_grid,
                            int icol,
                            double dlim,
                            double *dtab,
                            double *gtab);
GSTLEARN_EXPORT int db_edit(Db *db, int *flag_valid);
GSTLEARN_EXPORT int db_grid_copy(Db *db1,
                                 Db *db2,
                                 int *ind1,
                                 int *ind2,
                                 int ncol,
                                 int *cols);
GSTLEARN_EXPORT int db_grid_copy_dilate(Db *db1,
                                        int iatt1,
                                        Db *db2,
                                        int iatt2,
                                        int mode,
                                        int *nshift);
GSTLEARN_EXPORT int db_proportion(Db *db,
                                  Db *dbgrid,
                                  int nfac1max,
                                  int nfac2max,
                                  int *nclout);
GSTLEARN_EXPORT int db_merge(Db *db, int ncol, int *cols);
GSTLEARN_EXPORT int db_count_defined(Db *db, int icol);

GSTLEARN_EXPORT int db_prop_read(Db *db, int ix, int iy, double *props);
GSTLEARN_EXPORT int db_prop_write(Db *db, int ix, int iy, double *props);
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
GSTLEARN_EXPORT Db* db_point_init(int mode,
                                  int verbose,
                                  int ndim,
                                  int seed,
                                  double density,
                                  double range,
                                  double beta,
                                  Db *dbgrid,
                                  const VectorDouble &origin,
                                  const VectorDouble &extend);
GSTLEARN_EXPORT int db_smooth_vpc(Db *db, int width, double range);
GSTLEARN_EXPORT double* db_grid_sampling(Db *dbgrid,
                                         double *x1,
                                         double *x2,
                                         int ndisc,
                                         int ncut,
                                         double *cuts,
                                         int *nval_ret);
GSTLEARN_EXPORT int db_grid2point_sampling(Db *dbgrid,
                                           int nvar,
                                           int *vars,
                                           int *npacks,
                                           int npcell,
                                           int nmini,
                                           int *nech,
                                           double **coor,
                                           double **data);
GSTLEARN_EXPORT int db_grid_patch(Db *ss_grid,
                                  Db *db_grid,
                                  int iptr_ss,
                                  int iptr_db,
                                  int iptr_rank,
                                  int new_rank,
                                  int oper,
                                  int verbose);
GSTLEARN_EXPORT int db_polygon_distance(Db *db,
                                        Polygons *polygon,
                                        double dmax,
                                        int scale,
                                        int polin);

/****************************************/
/* Prototyping the functions in stats.c */
/****************************************/

GSTLEARN_EXPORT int stats_point_to_grid(Db *dbgrid,
                                        Db *db,
                                        const char *oper,
                                        int ivar,
                                        int jvar,
                                        int ncut,
                                        double *cuts,
                                        double *tab);
GSTLEARN_EXPORT int db_stats(Db *db,
                             const String &oper,
                             const VectorInt &cols,
                             int flag_mono,
                             int flag_verbose,
                             double *resta);
GSTLEARN_EXPORT int db_stats_grid(Db *db,
                                  Db *dbgrid,
                                  const char *oper,
                                  int ncol,
                                  int *cols,
                                  int radius);
GSTLEARN_EXPORT int stats_proportion(Db *dbin,
                                     Db *dbout,
                                     int pos,
                                     int nfacies,
                                     int radius);
GSTLEARN_EXPORT int stats_transition(Db *dbin,
                                     Db *dbout,
                                     int pos,
                                     int nfacies,
                                     int radius,
                                     int orient);
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
GSTLEARN_EXPORT int db_upscale(Db *dbgrid1,
                               Db *dbgrid2,
                               int orient,
                               int verbose);
GSTLEARN_EXPORT int db_diffusion(Db *dbgrid1,
                                 Db *dbgrid2,
                                 int orient,
                                 int niter,
                                 int nseed,
                                 int seed,
                                 int verbose);

/***************************************/
/* Prototyping the functions in skin.c */
/***************************************/

GSTLEARN_EXPORT Skin* skin_define(Db *db,
                                  int (*func_already_filled)(int ipos),
                                  int (*func_to_be_filled)(int ipos),
                                  double (*func_get_weight)(int ipos,
                                                            int idir));
GSTLEARN_EXPORT Skin* skin_undefine(Skin *skin);
GSTLEARN_EXPORT void skin_print(Skin *skin);
GSTLEARN_EXPORT int skin_init(Skin *skin, int verbose);
GSTLEARN_EXPORT int skin_remains(Skin *skin);
GSTLEARN_EXPORT void skin_next(Skin *skin, int *rank, int *ipos);
GSTLEARN_EXPORT int skin_unstack(Skin *skin, int rank, int ipos);
GSTLEARN_EXPORT int skin_grid_shift(Skin *skin, int lec, int dir, int *iad);

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

GSTLEARN_EXPORT int db_grid_read_zycor1(const char *filename,
                                        int verbose,
                                        int *nx,
                                        double *x0,
                                        double *dx);
GSTLEARN_EXPORT int db_grid_read_zycor2(const char *filename,
                                        int *nx,
                                        double *x0,
                                        double *dx,
                                        double *tab);
GSTLEARN_EXPORT int db_grid_read_bmp1(const char *filename,
                                      int verbose,
                                      int *nx,
                                      double *x0,
                                      double *dx);
GSTLEARN_EXPORT int db_grid_read_bmp2(const char *filename,
                                      int *nx,
                                      double *x0,
                                      double *dx,
                                      double *tab);
GSTLEARN_EXPORT int db_grid_read_prop1(const char *filename,
                                       int verbose,
                                       int *ncol,
                                       int *nx,
                                       double *x0,
                                       double *dx);
GSTLEARN_EXPORT int db_grid_read_prop2(const char *filename,
                                       int ncol_r,
                                       int *nx_r,
                                       double *x0_r,
                                       double *dx_r,
                                       double *tab);
GSTLEARN_EXPORT int db_grid_read_f2g(const char *filename,
                                     int verbose,
                                     int nx[2],
                                     double x0[3],
                                     double dx[3],
                                     double *angle,
                                     int *ncol,
                                     double **tab_arg);
GSTLEARN_EXPORT int db_grid_write_zycor(const char *filename, Db *db, int icol);
GSTLEARN_EXPORT int db_grid_write_XYZ(const char *filename, Db *db, int icol);
GSTLEARN_EXPORT int db_write_vtk(const char *filename,
                                 Db *db,
                                 const VectorInt &cols,
                                 const VectorString &names);
GSTLEARN_EXPORT int db_grid_write_bmp(const char *filename,
                                      Db *db,
                                      int icol,
                                      int nsamplex,
                                      int nsampley,
                                      int nmult,
                                      int ncolors,
                                      int flag_low,
                                      int flag_high,
                                      double valmin,
                                      double valmax,
                                      int *red,
                                      int *green,
                                      int *blue,
                                      int mask_red,
                                      int mask_green,
                                      int mask_blue,
                                      int ffff_red,
                                      int ffff_green,
                                      int ffff_blue,
                                      int low_red,
                                      int low_green,
                                      int low_blue,
                                      int high_red,
                                      int high_green,
                                      int highblue);
GSTLEARN_EXPORT int db_grid_write_irap(const char *filename,
                                       Db *db,
                                       int icol,
                                       int nsamplex,
                                       int nsampley);
GSTLEARN_EXPORT int db_grid_write_prop(const char *filename,
                                       Db *db,
                                       int ncol,
                                       int *icols);
GSTLEARN_EXPORT int db_grid_write_eclipse(const char *filename,
                                          Db *db,
                                          int icol);
GSTLEARN_EXPORT int db_well_read_las(const char *filename,
                                     int verbose,
                                     double xwell,
                                     double ywell,
                                     double cwell,
                                     int *nvarout,
                                     int *nechout,
                                     char ***var_names,
                                     double **tab);
GSTLEARN_EXPORT int csv_table_read(const String &filename,
                                   int verbose,
                                   int flag_header,
                                   int nskip,
                                   char char_sep,
                                   char char_dec,
                                   const String &na_string,
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
                                   double *rmean,
                                   double *rcov,
                                   double *smean);
GSTLEARN_EXPORT int image_smoother(Db *dbgrid,
                                   Neigh *neigh,
                                   int type,
                                   double range);
GSTLEARN_EXPORT int krigdgm_f(Db *dbin,
                              Db *dbout,
                              Model *model,
                              Neigh *neigh,
                              int flag_est,
                              int flag_std,
                              int flag_varz,
                              double rval);
GSTLEARN_EXPORT int krigcell_f(Db *dbin,
                               Db *dbout,
                               Model *model,
                               Neigh *neigh,
                               VectorInt ndisc,
                               int flag_est,
                               int flag_std,
                               VectorInt rank_colcok);
GSTLEARN_EXPORT int kriggam_f(Db *dbin,
                              Db *dbout,
                              Anam *anam,
                              Model *model,
                              Neigh *neigh);
GSTLEARN_EXPORT int krigprof_f(Db *dbin,
                               Db *dbout,
                               Model *model,
                               Neigh *neigh,
                               int ncode,
                               int flag_est,
                               int flag_std);
GSTLEARN_EXPORT int kribayes_f(Db *dbin,
                               Db *dbout,
                               Model *model,
                               Neigh *neigh,
                               double *dmean,
                               double *dcov,
                               int flag_est,
                               int flag_std);
GSTLEARN_EXPORT int krigsum_f(Db *dbin,
                              Db *dbout,
                              Model *model,
                              Neigh *neigh,
                              int flag_positive);
GSTLEARN_EXPORT int krigmvp_f(Db *dbin,
                              Db *db3grid,
                              Db *db2grid,
                              int fsum,
                              Model *model,
                              Neigh *neigh);

GSTLEARN_EXPORT int krigtest_dimension(Db *dbin,
                                       Db *dbout,
                                       Model *model,
                                       Neigh *neigh,
                                       int iech0,
                                       const EKrigOpt &calcul,
                                       VectorInt ndisc,
                                       int *ndim_ret,
                                       int *nech_ret,
                                       int *neq_ret,
                                       int *nrhs_ret);
GSTLEARN_EXPORT int krigtest_f(Db *dbin,
                               Db *dbout,
                               Model *model,
                               Neigh *neigh,
                               int iech0,
                               const EKrigOpt &calcul,
                               VectorInt ndisc,
                               int neq_out,
                               int nrhs_out,
                               double *xyz_out,
                               double *data_out,
                               double *lhs_out,
                               double *rhs_out,
                               double *wgt_out,
                               double *zam_out,
                               double *var_out);
GSTLEARN_EXPORT int krigsampling_f(Db *dbin,
                                   Db *dbout,
                                   Model *model,
                                   double beta,
                                   int nsize1,
                                   int *ranks1,
                                   int nsize2,
                                   int *ranks2,
                                   int flag_std,
                                   int verbose);
GSTLEARN_EXPORT int dk_f(Db *dbin,
                         Db *dbsmu,
                         Model *model,
                         Neigh *neigh,
                         int nfactor,
                         const VectorInt &nmult,
                         const VectorInt &ndisc,
                         int flag_est,
                         int flag_std);
GSTLEARN_EXPORT int global_arithmetic(Db *dbin,
                                      Db *dbgrid,
                                      Model *model,
                                      int ivar,
                                      int flag_verbose,
                                      int seed,
                                      double surface,
                                      double *zest,
                                      double *sse,
                                      double *cvgeo);
GSTLEARN_EXPORT int global_kriging(Db *dbin,
                                   Db *dbout,
                                   Model *model,
                                   int ivar,
                                   int flag_verbose,
                                   const EKrigOpt &calcul,
                                   int seed,
                                   double surface,
                                   double *zest,
                                   double *sse,
                                   double *cvgeo,
                                   double *weights);
GSTLEARN_EXPORT int global_transitive(Db *dbgrid,
                                      Model *model,
                                      int flag_verbose,
                                      int flag_regular,
                                      int ndisc,
                                      double *zest,
                                      double *cve,
                                      double *cvtrans);
GSTLEARN_EXPORT int invdist_f(Db *dbin,
                              Db *dbout,
                              int exponent,
                              int flag_expand,
                              double dmax);
GSTLEARN_EXPORT int anakexp_f(Db *db,
                              double *covdd,
                              double *covd0,
                              double top,
                              double bot,
                              int ncov_radius,
                              int neigh_radius,
                              int flag_sym,
                              int nfeq);
GSTLEARN_EXPORT int anakexp_3D(Db *db,
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
GSTLEARN_EXPORT int declustering_f(Db *db,
                                   Model *model,
                                   Neigh *neigh,
                                   Db *dbgrid,
                                   int method,
                                   double *radius,
                                   VectorInt ndisc,
                                   int flag_sel,
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
                                               const ELoc &locatorType);
GSTLEARN_EXPORT int simtub_workable(Model *model);
GSTLEARN_EXPORT int simdgm(Db *dbin,
                           Db *dbout,
                           Model *model,
                           Neigh *neigh,
                           double rval,
                           int seed,
                           int nbsimu,
                           int nbtuba,
                           int flag_check);
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
GSTLEARN_EXPORT int simbayes(Db *dbin,
                             Db *dbout,
                             Model *model,
                             Neigh *neigh,
                             double *dmean,
                             double *dcov,
                             int seed,
                             int nbsimu,
                             int nbtuba,
                             int flag_check);
GSTLEARN_EXPORT int simmaxstable(Db *dbout,
                                 Model *model,
                                 double ratio,
                                 int seed,
                                 int nbtuba,
                                 int flag_simu,
                                 int flag_rank,
                                 int verbose);
GSTLEARN_EXPORT int simtub_potential(Db *dbiso,
                                     Db *dbgrd,
                                     Db *dbtgt,
                                     Db *dbout,
                                     Model *model,
                                     int nbsimu,
                                     int nbtuba,
                                     double delta);
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
                                       Neigh *neigh,
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
                                         const ELoc &locatorType,
                                         int nbsimu,
                                         int nvar,
                                         int *iptr_ce_arg,
                                         int *iptr_cstd_arg);

/*****************************************/
/* Prototyping the functions in simfft.c */
/*****************************************/

GSTLEARN_EXPORT int simfft_f(Db *db,
                             Model *model,
                             int seed,
                             int nbsimu,
                             double percent,
                             int flag_aliasing);
GSTLEARN_EXPORT int simfft_support(Db *db,
                                   Model *model,
                                   double percent,
                                   int flag_aliasing,
                                   int nval,
                                   double *sigma,
                                   double *r2val,
                                   double *coeffs);

/*****************************************/
/* Prototyping the functions in simreg.c */
/*****************************************/

GSTLEARN_EXPORT int simfine_dim(Db *dbin,
                                int nmult,
                                int *ndim,
                                int *ntot,
                                int *nx,
                                double *x0,
                                double *dx);
GSTLEARN_EXPORT int simfine_f(Db *dbin,
                              Model *model,
                              int flag_ks,
                              int mult,
                              int seed,
                              VectorDouble &tab);

/*****************************************/
/* Prototyping the functions in simsub.c */
/*****************************************/

GSTLEARN_EXPORT int substitution(Db *dbgrid,
                                 int seed,
                                 int nfacies,
                                 int nstates,
                                 int flag_direct,
                                 int flag_coding,
                                 int flag_orient,
                                 int flag_auto,
                                 double intensity,
                                 double factor,
                                 double vector[3],
                                 double *trans,
                                 int colfac,
                                 int colang[3],
                                 int verbose);

/******************************************/
/* Prototyping the functions in simpart.c */
/******************************************/

GSTLEARN_EXPORT SubPlanes* poisson_manage_planes(int mode,
                                                 int np,
                                                 SubPlanes *splanes);
GSTLEARN_EXPORT int poisson_generate_planes(Db *dbgrid, SubPlanes *splanes);
GSTLEARN_EXPORT int tessellation_poisson(Db *dbgrid,
                                         Model *model,
                                         int seed,
                                         double intensity,
                                         int nbtuba,
                                         int verbose);
GSTLEARN_EXPORT int tessellation_voronoi(Db *dbgrid,
                                         Model *model,
                                         double *dilate,
                                         int seed,
                                         double intensity,
                                         int nbtuba,
                                         int verbose);

/*****************************************/
/* Prototyping the functions in simsph.c */
/*****************************************/
GSTLEARN_EXPORT int simsph_f(Db *db,
                             Model *model,
                             int seed,
                             int special,
                             int nbf,
                             int nfmax,
                             int verbose,
                             int flag_test,
                             int test_degree,
                             int test_order,
                             double test_phase);

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

GSTLEARN_EXPORT int seismic_estimate_XZ(Db *db,
                                        Model *model,
                                        int nbench,
                                        int nv2max,
                                        int flag_ks,
                                        int flag_std,
                                        int flag_sort,
                                        int flag_stat);
GSTLEARN_EXPORT int seismic_simulate_XZ(Db *db,
                                        Model *model,
                                        int nbench,
                                        int nv2max,
                                        int nbsimu,
                                        int seed,
                                        int flag_ks,
                                        int flag_sort,
                                        int flag_stat);
GSTLEARN_EXPORT int seismic_z2t_grid(int verbose,
                                     Db *db_z,
                                     int iptr_v,
                                     int *nx,
                                     double *x0,
                                     double *dx);
GSTLEARN_EXPORT int seismic_t2z_grid(int verbose,
                                     Db *db_t,
                                     int iptr_v,
                                     int *nx,
                                     double *x0,
                                     double *dx);
GSTLEARN_EXPORT int seismic_z2t_convert(Db *db_z, int iptr_v, Db *db_t);
GSTLEARN_EXPORT int seismic_t2z_convert(Db *db_t, int iptr_v, Db *db_z);
GSTLEARN_EXPORT int seismic_operate(Db *db, int oper);
GSTLEARN_EXPORT int seismic_convolve(Db *db,
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

/******************************************/
/* Prototyping the functions in simbool.c */
/******************************************/

GSTLEARN_EXPORT Tokens* tokens_free(Tokens *tokens);
GSTLEARN_EXPORT Tokens* tokens_create(int nb_tokens);
GSTLEARN_EXPORT int tokone_create(Tokens *token,
                                  int rank,
                                  int type,
                                  int npar,
                                  double prop,
                                  double factor_x2y,
                                  double factor_x2z,
                                  double factor_y2z,
                                  int *law,
                                  double *valarg);
GSTLEARN_EXPORT void tokone_print(Tokens *tokens, int rank);
GSTLEARN_EXPORT Tokens* tokens_input(void);
GSTLEARN_EXPORT void tokone_get_nbparams(Tokens *tokens,
                                         int rank,
                                         int *types,
                                         int *npar,
                                         double *prop);
GSTLEARN_EXPORT void tokone_get_params(Tokens *tokens,
                                       int rank,
                                       double *factor_x2y,
                                       double *factor_x2z,
                                       double *factor_y2z,
                                       int *law,
                                       double *valarg);
GSTLEARN_EXPORT void tokens_print(Tokens *tokens);
GSTLEARN_EXPORT int toktype_get_nbparams(int type);
GSTLEARN_EXPORT int simbool_f(Db *dbin,
                              Db *dbout,
                              Tokens *tokens,
                              int seed,
                              int nb_average,
                              int flag_stat,
                              int flag_simu,
                              int flag_rank,
                              double background,
                              double facies,
                              double *dilate,
                              double theta_cste,
                              double tmax,
                              int verbose);

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

/***************************************/
/* Prototyping the functions in poly.c */
/***************************************/

GSTLEARN_EXPORT Polygons* polygon_create(void);
GSTLEARN_EXPORT Polygons* polygon_free(Polygons *polygon);
GSTLEARN_EXPORT Polygons* polygon_add(Polygons *polygon,
                                      const VectorDouble &x,
                                      const VectorDouble &y,
                                      double zmin,
                                      double zmax);
GSTLEARN_EXPORT void polygon_print(Polygons *polygon, int flag_print);
GSTLEARN_EXPORT int polygon_inside(double xx,
                                   double yy,
                                   double zz,
                                   int flag_nested,
                                   Polygons *polygon);
GSTLEARN_EXPORT void polygon_extension(Polygons *polygon,
                                       double *xmin,
                                       double *xmax,
                                       double *ymin,
                                       double *ymax);
GSTLEARN_EXPORT double polygon_surface(Polygons *polygon);
GSTLEARN_EXPORT Polygons* input_polygon(void);
GSTLEARN_EXPORT Polygons* polygon_hull(const Db *db);
GSTLEARN_EXPORT int polygon_hull(const Db *db,
                                 VectorDouble &x,
                                 VectorDouble &y);

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

/**************************/
/* Prototyping fracture.c */
/**************************/
GSTLEARN_EXPORT Frac_Environ* fracture_alloc_environ(int nfamilies,
                                                     double xmax,
                                                     double ymax,
                                                     double deltax,
                                                     double deltay,
                                                     double mean,
                                                     double stdev);
GSTLEARN_EXPORT Frac_Environ* fracture_dealloc_environ(Frac_Environ *frac_environ);

GSTLEARN_EXPORT Frac_List* fracture_manage_list(int mode, Frac_List *frac_list);
GSTLEARN_EXPORT void fracture_update_family(Frac_Environ *frac_environ,
                                            int family,
                                            double orient,
                                            double dorient,
                                            double theta0,
                                            double alpha,
                                            double ratcst,
                                            double prop1,
                                            double prop2,
                                            double aterm,
                                            double bterm,
                                            double range);
GSTLEARN_EXPORT int fracture_add_fault(Frac_Environ *frac_environ,
                                       double fault_coord,
                                       double fault_orient);
GSTLEARN_EXPORT void fracture_update_fault(Frac_Environ *frac_environ,
                                           int ifault,
                                           int family,
                                           double thetal,
                                           double thetar,
                                           double rangel,
                                           double ranger);
GSTLEARN_EXPORT int fracture_simulate(Frac_Environ *frac_environ,
                                      int flag_sim_layer,
                                      int flag_sim_fract,
                                      int seed,
                                      int verbose,
                                      int nlayers_in,
                                      double *elevations,
                                      int *nlayers,
                                      int *ninfos,
                                      double **layinfo,
                                      Frac_List *frac_list);
GSTLEARN_EXPORT Frac_Environ* fracture_input(Frac_Environ *frac_def);
GSTLEARN_EXPORT void fracture_print(Frac_Environ *frac_environ);
GSTLEARN_EXPORT void fracture_list_print(const char *title,
                                         Frac_List *frac_list,
                                         int level);
GSTLEARN_EXPORT void fracture_export(Frac_List *frac_list,
                                     int *nfracs_arg,
                                     int *nbyfrac_arg,
                                     double **frac_segs_arg);
GSTLEARN_EXPORT Frac_List* fracture_import(int nval, double *frac_segs);
GSTLEARN_EXPORT double* fracture_extract_length(Frac_List *frac_list,
                                                int family,
                                                double cote,
                                                double dcote,
                                                int *ntab);
GSTLEARN_EXPORT double* fracture_extract_dist(Frac_List *frac_list,
                                              int family,
                                              double cote,
                                              double dcote,
                                              int *ntab);
GSTLEARN_EXPORT int fracture_to_block(Db *dbgrid,
                                      Frac_List *frac_list,
                                      double *locinfo,
                                      int n_layers,
                                      int nfamilies,
                                      double xmax,
                                      double *permtab,
                                      double perm_mat,
                                      double perm_bench,
                                      int ndisc);
GSTLEARN_EXPORT double* fracture_to_well(int nw_xy,
                                         double *well,
                                         Frac_List *frac_list,
                                         double xmax,
                                         double *permtab,
                                         int *nint,
                                         int *ncol);
GSTLEARN_EXPORT int fracture_well_to_block(Db *dbgrid,
                                           Frac_List *frac_list,
                                           int col_perm,
                                           int col_fluid,
                                           int flag_fluid,
                                           double val_fluid,
                                           double *wellout,
                                           int nw_xy,
                                           int ndisc,
                                           int verbose);

/******************************************/
/* Prototyping the functions in mlayers.c */
/******************************************/
GSTLEARN_EXPORT int variogram_mlayers(Db *db,
                                      int *seltab,
                                      Vario *vario,
                                      Vario_Order *vorder);
GSTLEARN_EXPORT int multilayers_vario(Db *dbin,
                                      Db *dbout,
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
                                        Db *dbout,
                                        Model *model,
                                        Neigh *neigh,
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
                                          Db *dbout,
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

/********************************************/
/* Prototyping the functions in potential.c */
/********************************************/
GSTLEARN_EXPORT int potential_cov(Model *model,
                                  int verbose,
                                  int type1,
                                  double *x10,
                                  double *x1p,
                                  double *tx1,
                                  int type2,
                                  double *x20,
                                  double *x2p,
                                  double *tx2,
                                  int *n1,
                                  int *n2,
                                  double *covtab);
GSTLEARN_EXPORT int potential_kriging(Db *db,
                                      Db *dbgrd,
                                      Db *dbtgt,
                                      Db *dbout,
                                      Model *model,
                                      Neigh *neigh,
                                      double nugget_grd,
                                      double nugget_tgt,
                                      int flag_grad,
                                      int flag_trans,
                                      int flag_drift,
                                      int verbose);
GSTLEARN_EXPORT int potential_simulate(Db *dbiso,
                                       Db *dbgrd,
                                       Db *dbtgt,
                                       Db *dbout,
                                       Model *model,
                                       Neigh *neigh,
                                       double nugget_grd,
                                       double nugget_tgt,
                                       double dist_tempere,
                                       int flag_trans,
                                       int seed,
                                       int nbsimu,
                                       int nbtuba,
                                       int verbose);
GSTLEARN_EXPORT int potential_xvalid(Db *dbiso,
                                     Db *dbgrd,
                                     Db *dbtgt,
                                     Model *model,
                                     Neigh *neigh,
                                     double nugget_grd,
                                     double nugget_tgt,
                                     int flag_dist_conv,
                                     int verbose);

/*******************************************/
/* Prototyping the functions in delaunay.c */
/*******************************************/
GSTLEARN_EXPORT int MSS(int idim, int ipol, int icas, int icorn, int icoor);
GSTLEARN_EXPORT Vercoloc* vercoloc_manage(int verbose,
                                          int mode,
                                          Db *dbin,
                                          Db *dbgrid,
                                          int mesh_dbin,
                                          Vercoloc *vercoloc);
GSTLEARN_EXPORT Vercoloc* vercoloc_from_external(int ndupl,
                                                 int *dupl_in,
                                                 int *dupl_out);
GSTLEARN_EXPORT Vertype* vertype_manage(int mode,
                                        Vertype *vertype,
                                        Vercoloc *vercoloc,
                                        int nvertex);
GSTLEARN_EXPORT int* vercoloc_get_dbin_indices(Vertype *vertype,
                                               Vercoloc *vercoloc,
                                               int *nbnodup);
GSTLEARN_EXPORT void triangulate(const char *triswitches,
                                 struct triangulateio *in,
                                 struct triangulateio *out,
                                 struct triangulateio *vorout);

GSTLEARN_EXPORT void meshes_1D_free(segmentio *t, int mode);
GSTLEARN_EXPORT void meshes_1D_init(segmentio *t);
GSTLEARN_EXPORT int meshes_1D_from_db(Db *db,
                                      int nmask,
                                      int *mask,
                                      segmentio *t);
GSTLEARN_EXPORT int meshes_1D_from_points(int nech, double *x, segmentio *t);
GSTLEARN_EXPORT void meshes_1D_default(Db *dbin, Db *dbout, segmentio *t);
GSTLEARN_EXPORT void meshes_1D_print(segmentio *t, int brief);
GSTLEARN_EXPORT int meshes_turbo_1D_grid_build(int verbose,
                                               Db *dbgrid,
                                               SPDE_Mesh *s_mesh);
GSTLEARN_EXPORT void meshes_1D_create(int verbose,
                                      struct segmentio *in,
                                      struct segmentio *out);
GSTLEARN_EXPORT void meshes_1D_load_vertices(segmentio *t,
                                             const char *name,
                                             int *ntab_arg,
                                             int *natt_arg,
                                             void **tab_arg);
GSTLEARN_EXPORT void meshes_1D_extended_domain(Db *dbout,
                                               const double *gext,
                                               segmentio *t);

GSTLEARN_EXPORT void meshes_2D_free(triangulateio *t, int mode);
GSTLEARN_EXPORT void meshes_2D_init(triangulateio *t);
GSTLEARN_EXPORT int meshes_2D_from_db(Db *db,
                                      int use_code,
                                      int nmask,
                                      int *mask,
                                      triangulateio *t);
GSTLEARN_EXPORT int meshes_2D_from_points(int nech,
                                          double *x,
                                          double *y,
                                          triangulateio *t);
GSTLEARN_EXPORT void meshes_2D_default(Db *dbin, Db *dbout, triangulateio *t);
GSTLEARN_EXPORT int meshes_2D_from_mem(int nseg,
                                       int ncol,
                                       int *segments,
                                       triangulateio *t);
GSTLEARN_EXPORT void meshes_2D_print(triangulateio *t, int brief);
GSTLEARN_EXPORT int meshes_2D_write(const char *file_name,
                                    const char *obj_name,
                                    int verbose,
                                    int ndim,
                                    int ncode,
                                    int ntri,
                                    int npoints,
                                    int *ntcode,
                                    int *triangles,
                                    double *points);
GSTLEARN_EXPORT int meshes_turbo_2D_grid_build(int verbose,
                                               Db *dbgrid,
                                               SPDE_Mesh *s_mesh);
GSTLEARN_EXPORT void meshes_2D_create(int verbose,
                                      const String &triswitches,
                                      struct triangulateio *in,
                                      struct triangulateio *out,
                                      struct triangulateio *vorout);
GSTLEARN_EXPORT void meshes_2D_load_vertices(triangulateio *t,
                                             const char *name,
                                             int *ntab_arg,
                                             int *natt_arg,
                                             void **tab_arg);
GSTLEARN_EXPORT void meshes_2D_extended_domain(Db *dbout,
                                               const double *gext,
                                               triangulateio *t);

GSTLEARN_EXPORT void meshes_3D_create(int verbose,
                                      const String &triswitch,
                                      tetgenio *in,
                                      tetgenio *out);
GSTLEARN_EXPORT int meshes_3D_from_db(Db *db,
                                      int nmask,
                                      int *mask,
                                      tetgenio *t);
GSTLEARN_EXPORT int meshes_3D_from_points(int nech,
                                          double *x,
                                          double *y,
                                          double *z,
                                          tetgenio *t);
GSTLEARN_EXPORT void meshes_3D_default(Db *dbin, Db *dbout, tetgenio *t);
GSTLEARN_EXPORT void meshes_3D_free(tetgenio *t);
GSTLEARN_EXPORT void meshes_3D_extended_domain(Db *dbout,
                                               const double *gext,
                                               tetgenio *t);
GSTLEARN_EXPORT int meshes_turbo_3D_grid_build(int verbose,
                                               Db *dbgrid,
                                               SPDE_Mesh *s_mesh);
GSTLEARN_EXPORT void meshes_3D_load_vertices(tetgenio *t,
                                             const char *name,
                                             int *ntab_arg,
                                             int *natt_arg,
                                             void **tab_arg);
GSTLEARN_EXPORT void meshes_3D_print(tetgenio *t, int brief);
GSTLEARN_EXPORT void mesh_stats(int ndim,
                                int ncorner,
                                int nmesh,
                                int *meshes,
                                double *points);
GSTLEARN_EXPORT void meshes_2D_sph_init(SphTriangle *t);
GSTLEARN_EXPORT void meshes_2D_sph_free(SphTriangle *t, int mode);
GSTLEARN_EXPORT int meshes_2D_sph_from_db(Db *db,
                                          int nmask,
                                          int *mask,
                                          SphTriangle *t);
GSTLEARN_EXPORT int meshes_2D_sph_from_points(int nech,
                                              double *x,
                                              double *y,
                                              SphTriangle *t);
GSTLEARN_EXPORT int meshes_2D_sph_from_auxiliary(const String &triswitch,
                                                 SphTriangle *t);
GSTLEARN_EXPORT void meshes_2D_sph_print(SphTriangle *t, int brief);
GSTLEARN_EXPORT int meshes_2D_sph_create(int verbose, SphTriangle *t);
GSTLEARN_EXPORT void meshes_2D_sph_load_vertices(SphTriangle *t,
                                                 const char *name,
                                                 int *ntab_arg,
                                                 int *natt_arg,
                                                 void **tab_arg);
GSTLEARN_EXPORT int trmesh_(int *n,
                            double *x,
                            double *y,
                            double *z__,
                            int *list,
                            int *lptr,
                            int *lend,
                            int *lnew,
                            int *near__,
                            int *next,
                            double *dist,
                            int *ier);
GSTLEARN_EXPORT int trlist_(int *n,
                            int *list,
                            int *lptr,
                            int *lend,
                            int *nrow,
                            int *nt,
                            int *ltri,
                            int *ier);
GSTLEARN_EXPORT void util_convert_sph2cart(double rlong,
                                           double rlat,
                                           double *x,
                                           double *y,
                                           double *z);
GSTLEARN_EXPORT void util_convert_cart2sph(double x,
                                           double y,
                                           double z,
                                           double *rlong,
                                           double *rlat);

/***************************************/
/* Prototyping the functions in spde.c */
/***************************************/
GSTLEARN_EXPORT QChol* qchol_manage(int mode, QChol *qchol);
GSTLEARN_EXPORT double spde_compute_correc(int ndim, double param);
GSTLEARN_EXPORT int spde_check(const Db *dbin,
                               const Db *dbout,
                               Model *model1,
                               Model *model2,
                               int verbose,
                               const VectorDouble &gext,
                               int mesh_dbin,
                               int mesh_dbout,
                               int flag_advanced,
                               int flag_est,
                               int flag_std,
                               int flag_gibbs,
                               int flag_modif);
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
GSTLEARN_EXPORT int spde_posterior(Db *dbin,
                                   Db *dbout,
                                   const VectorDouble &gext,
                                   SPDE_Option &s_option);
GSTLEARN_EXPORT int spde_process(Db *dbin,
                                 Db *dbout,
                                 SPDE_Option &s_option,
                                 int nbsimu,
                                 int gibbs_nburn,
                                 int gibbs_niter,
                                 int ngibbs_int);
GSTLEARN_EXPORT SPDE_Mesh* spde_mesh_manage(int mode, SPDE_Mesh *s_mesh_old);
GSTLEARN_EXPORT SPDE_Matelem& spde_get_current_matelem(int icov);
GSTLEARN_EXPORT int spde_mesh_load(SPDE_Mesh *s_mesh,
                                   int verbose,
                                   Db *dbin,
                                   Db *dbout,
                                   const VectorDouble &gext,
                                   SPDE_Option &s_option);
GSTLEARN_EXPORT void spde_mesh_assign(SPDE_Mesh *s_mesh,
                                      int ndim,
                                      int ncorner,
                                      int nvertex,
                                      int nmesh,
                                      int *meshes,
                                      double *points,
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
GSTLEARN_EXPORT int spde_eval(int nblin,
                              double *blin,
                              cs *S,
                              const VectorDouble &Lambda,
                              const VectorDouble &TildeC,
                              double power,
                              double *x,
                              double *y);
GSTLEARN_EXPORT int spde_external_mesh_define(int mode,
                                              int icov0,
                                              int ndim,
                                              int ncorner,
                                              int nvertex,
                                              int nmesh,
                                              int nbin,
                                              int nbout,
                                              int ndupl,
                                              int order,
                                              int *dupl_in,
                                              int *dupl_out,
                                              int *meshes,
                                              double *points);
GSTLEARN_EXPORT int spde_external_AQ_copy(SPDE_Matelem &matelem, int icov0);
GSTLEARN_EXPORT int spde_external_AQ_define(int mode,
                                            int icov0,
                                            int ndim,
                                            int nvertex,
                                            int nmesh,
                                            int nbin,
                                            int nbout,
                                            int ndupl,
                                            int order,
                                            int *dupl_in,
                                            int *dupl_out,
                                            cs *A,
                                            cs *Q);
GSTLEARN_EXPORT int spde_external_mesh_copy(SPDE_Mesh *s_mesh, int icov0);
GSTLEARN_EXPORT int kriging2D_spde(Db *dbin,
                                   Model *model,
                                   SPDE_Option &s_option,
                                   int verbose,
                                   int *ntri,
                                   int *npoint,
                                   int **triangles,
                                   double **points);
GSTLEARN_EXPORT cs* db_mesh_sparse(Db *db, MeshEStandard *amesh, int verbose);
GSTLEARN_EXPORT cs* db_mesh_neigh(const Db *db,
                                  SPDE_Mesh *s_mesh,
                                  double radius,
                                  int flag_exact,
                                  int verbose,
                                  int *nactive,
                                  int **ranks);
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
                               int **ntcode_arg,
                               int **triangle_arg,
                               double **points_arg);
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

