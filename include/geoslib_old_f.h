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
#ifndef GEOSLIB_OLDF_H
#define GEOSLIB_OLDF_H

#include "geoslib_d.h"
#include "csparse_d.h"
#include "csparse_f.h"
#include "Mesh/tetgen.h"
#include "segy.h"
#include "Neigh/Neigh.hpp"
#include "Variogram/Vario.hpp"
#include "Variogram/Dir.hpp"
#include "Model/Model.hpp"
#include "Model/Cova.hpp"
#include "Db/Db.hpp"
#include "Anamorphosis/Anam.hpp"
#include "Anamorphosis/AnamDiscreteDD.hpp"
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Anamorphosis/AnamEmpirical.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Anamorphosis/AnamUser.hpp"
#include "Model/CovNostatInternal.hpp"
#include "Model/Constraints.hpp"
#include "Stats/PCA.hpp"
#include "Mesh/MeshEStandard.hpp"
#include "Polygon/Polygons.hpp"
#include "Basic/NamingConvention.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"
#include "LithoRule/Rule.hpp"
#include "Model/ANoStat.hpp"

//#ifdef __cplusplus
//extern "C"
//{
//#endif

  /**********************************************/
  /* Prototyping the functions in acknowledge.c */
  /**********************************************/
  GEOSLIB_API void acknowledge_Geoslib(void);
  GEOSLIB_API void inquire_Geoslib(char **release, char **date);

  /******************************************/
  /* Prototyping the functions in license.c */
  /******************************************/
  GEOSLIB_API int register_license_file(const char *file_name,
                                        const char *target_name);

  /***************************************/
  /* Prototyping the functions in pile.c */
  /***************************************/
  GEOSLIB_API void pile_reset(int type);
  GEOSLIB_API void piles_reset(void);
  GEOSLIB_API int pile_next(int type);
  GEOSLIB_API void pile_manage(int type, int rank, int mode, char *ptr);
  GEOSLIB_API int pile_correct(int type, int rank, int mode);
  GEOSLIB_API char *pile_get(int type, int rank);
  GEOSLIB_API void piles_dump(void);

  /**************************************/
  /* Prototyping the functions in fft.c */
  /**************************************/

  GEOSLIB_API int fftn(int ndim,
                       const int dims[],
                       double Re[],
                       double Im[],
                       int iSign,
                       double scaling);

  /***************************************/
  /* Prototyping the functions in math.c */
  /***************************************/

  GEOSLIB_API int opt_mauto_add_constraints(Option_AutoFit& mauto,
                                            double constantSill);
  GEOSLIB_API int opt_mauto_add_unit_constraints(Option_AutoFit& mauto);
  GEOSLIB_API int foxleg_f(int ndat,
                           int npar,
                           int ncont,
                           const VectorDouble& acont,
                           VectorDouble& param,
                           VectorDouble& lower,
                           VectorDouble& upper,
                           VectorDouble& scale,
                           const Option_AutoFit& mauto,
                           int flag_title,
                           void (*func_evaluate)(int ndat,
                                                 int npar,
                                                 VectorDouble& param,
                                                 VectorDouble& work),
                           VectorDouble& tabexp,
                           VectorDouble& tabwgt);
  GEOSLIB_API int mvndst_infin(double low, double sup);
  GEOSLIB_API void mvndst(int n,
                          double *lower,
                          double *upper,
                          int *infin,
                          double *correl,
                          int maxpts,
                          double abseps,
                          double releps,
                          double *error,
                          double *value,
                          int *inform);
  GEOSLIB_API void mvndst2n(double *lower,
                            double *upper,
                            double *means,
                            double *correl,
                            int maxpts,
                            double abseps,
                            double releps,
                            double *error,
                            double *value,
                            int *inform);
  GEOSLIB_API void mvndst4(double *lower,
                           double *upper,
                           double *correl,
                           int maxpts,
                           double abseps,
                           double releps,
                           double *error,
                           double *value,
                           int *inform);

  /*****************************************/
  /* Prototyping the functions in bessel.c */
  /*****************************************/

  GEOSLIB_API int bessel_j(double x, double alpha, int nb, double *b);
  GEOSLIB_API int bessel_k(double x, double alpha, int nb, double *bk);

  /**************************************/
  /* Prototyping the functions in law.c */
  /**************************************/

  GEOSLIB_API int law_get_random_seed(void);
  GEOSLIB_API void law_set_random_seed(int seed);
  GEOSLIB_API double law_uniform(double mini, double maxi);
  GEOSLIB_API int law_int_uniform(int mini, int maxi);
  GEOSLIB_API double law_gaussian(void);
  GEOSLIB_API double law_exponential(void);
  GEOSLIB_API double law_gamma(double parameter);
  GEOSLIB_API int law_poisson(double parameter);
  GEOSLIB_API double law_stable_standard_agd(double alpha, double beta);
  GEOSLIB_API double law_stable_standard_a1gd(double beta);
  GEOSLIB_API double law_stable_standard_abgd(double alpha);
  GEOSLIB_API double law_stable_a(double alpha,
                                  double beta,
                                  double gamma,
                                  double delta);
  GEOSLIB_API double law_stable_a1(double beta, double gamma, double delta);
  GEOSLIB_API double law_stable(double alpha,
                                double beta,
                                double gamma,
                                double delta);
  GEOSLIB_API int law_binomial(int n, double p);
  GEOSLIB_API double law_beta1(double parameter1, double parameter2);
  GEOSLIB_API double law_beta2(double parameter1, double parameter2);
  GEOSLIB_API double law_df_gaussian(double value);
  GEOSLIB_API double law_dnorm(double value, double mean, double std);
  GEOSLIB_API double law_cdf_gaussian(double value);
  GEOSLIB_API double law_invcdf_gaussian(double value);
  GEOSLIB_API double law_gaussian_between_bounds(double binf, double bsup);
  GEOSLIB_API double law_df_bigaussian(double *vect,
                                       double *mean,
                                       double *corr);
  GEOSLIB_API double law_df_quadgaussian(double *vect, double *corr);
  GEOSLIB_API double law_df_multigaussian(int nvar, double *vect, double *corr);
  GEOSLIB_API void law_random_path(int nech, int *path);
  GEOSLIB_API double *law_exp_sample(double *tabin,
                                     int mode,
                                     int nvar,
                                     int nechin,
                                     int nechout,
                                     int niter,
                                     int nconst,
                                     double *consts,
                                     int seed,
                                     double percent);

  /***************************************/
  /* Prototyping the functions in util.c */
  /***************************************/

  GEOSLIB_API double ut_deg2rad(double angle);
  GEOSLIB_API double ut_rad2deg(double angle);
  GEOSLIB_API int get_mirror_sample(int nx, int ix);
  GEOSLIB_API void get_matrix(const char *title,
                              int flag_sym,
                              int flag_def,
                              int nx,
                              int ny,
                              double valmin,
                              double valmax,
                              double *tab);
  GEOSLIB_API void get_rotation(const char *title,
                                int flag_def,
                                int ndim,
                                double *rot);
  GEOSLIB_API void ut_sort_double(int safe, int nech, int *ind, double *tab);
  GEOSLIB_API void ut_sort_int(int safe, int nech, int *ind, int *tab);
  GEOSLIB_API void ut_tab_unique(int ntab, double *tab, int *neff);
  GEOSLIB_API void ut_statistics(int nech,
                                 double *tab,
                                 double *sel,
                                 double *wgt,
                                 int *nval,
                                 double *mini,
                                 double *maxi,
                                 double *delta,
                                 double *mean,
                                 double *stdv);
  GEOSLIB_API double ut_cnp(int n, int k);
  GEOSLIB_API int *ut_combinations(int n, int maxk, int *ncomb);
  GEOSLIB_API int *ut_split_into_two(int ncolor,
                                     int flag_half,
                                     int verbose,
                                     int *nposs);
  GEOSLIB_API double *ut_pascal(int ndim);
  GEOSLIB_API double ut_median(double *tab, int ntab);
  GEOSLIB_API void rgb2num(int r, int g, int b, int a, unsigned char *value);
  GEOSLIB_API void num2rgb(unsigned char value, int *r, int *g, int *b, int *a);
  GEOSLIB_API void ut_stats_mima(int nech,
                                 double *tab,
                                 double *sel,
                                 int *nvalid,
                                 double *mini,
                                 double *maxi);
  GEOSLIB_API void ut_stats_mima_print(const char *title,
                                       int nech,
                                       double *tab,
                                       double *sel);
  GEOSLIB_API void ut_facies_statistics(int nech,
                                        double *tab,
                                        double *sel,
                                        int *nval,
                                        int *mini,
                                        int *maxi);
  GEOSLIB_API void ut_classify(int nech,
                               double *tab,
                               double *sel,
                               int nclass,
                               double start,
                               double pas,
                               int *nmask,
                               int *ntest,
                               int *nout,
                               int *classe);
  GEOSLIB_API void ut_normalize(int ntab, double *tab);
  GEOSLIB_API void ut_rotation_sincos(double angle, double *cosa, double *sina);
  GEOSLIB_API void ut_rotation_matrix_2D(double angle, double *rot);
  GEOSLIB_API void ut_rotation_matrix_3D(double alpha,
                                         double beta,
                                         double gamma,
                                         double *rot);
  GEOSLIB_API void ut_rotation_matrix(int ndim,
                                      const double *angles,
                                      double *rot);
  GEOSLIB_API VectorDouble ut_rotation_matrix_VD(int ndim,
                                                 const VectorDouble& angles);
  GEOSLIB_API void ut_rotation_init(int ndim, double *rot);
  GEOSLIB_API int ut_rotation_check(double *rot, int ndim);
  GEOSLIB_API void ut_rotation_copy(int ndim,
                                    const double *rotin,
                                    double *rotout);
  GEOSLIB_API void ut_rotation_direction(double ct,
                                         double st,
                                         double *a,
                                         double *codir);
  GEOSLIB_API int ut_angles_from_rotation_matrix(const double *rot,
                                                 int ndim,
                                                 double *angles);
  GEOSLIB_API void ut_angles_from_codir(int ndim,
                                        int ndir,
                                        const VectorDouble& codir,
                                        VectorDouble& angles);
  GEOSLIB_API void ut_angles_to_codir(int ndim,
                                      int ndir,
                                      const VectorDouble& angles,
                                      VectorDouble& codr);
  GEOSLIB_API double ut_merge_extension(int ndim,
                                        double *mini_in,
                                        double *maxi_in,
                                        double *mini_out,
                                        double *maxi_out);
  GEOSLIB_API void debug_reset(void);
  GEOSLIB_API void debug_print(void);
  GEOSLIB_API void debug_index(int rank);
  GEOSLIB_API void debug_reference(int rank);
  GEOSLIB_API int is_debug_reference_defined(void);
  GEOSLIB_API void debug_define(const char *name, int status);
  GEOSLIB_API int debug_query(const char *name);
  GEOSLIB_API int debug_force(void);
  GEOSLIB_API void string_to_uppercase(char *string);
  GEOSLIB_API void string_to_lowercase(char *string);
  GEOSLIB_API int string_compare(int flag_case,
                                 const char *string1,
                                 const char *string2);
  GEOSLIB_API void projec_query(int *actif);
  GEOSLIB_API void projec_print(void);
  GEOSLIB_API void projec_toggle(int mode);
  GEOSLIB_API void variety_define(int flag_sphere, double radius);
  GEOSLIB_API void variety_query(int *flag_sphere);
  GEOSLIB_API void variety_print(void);
  GEOSLIB_API void variety_toggle(int mode);
  GEOSLIB_API void variety_get_characteristics(double *radius);
  GEOSLIB_API double ut_factorial(int k);
  GEOSLIB_API void ut_log_factorial(int nbpoly, double *factor);
  GEOSLIB_API double golden_search(double (*func_evaluate)(double test,
                                                           void *user_data),
                                   void *user_data,
                                   double tolstop,
                                   double a0,
                                   double c0,
                                   double *testval,
                                   double *niter);
  GEOSLIB_API double loggamma(double parameter);
  GEOSLIB_API void set_last_message(int mode, const char *string);
  GEOSLIB_API void print_last_message(void);
  GEOSLIB_API void set_keypair(const char *keyword,
                               int origin,
                               int nrow,
                               int ncol,
                               const double *values);
  GEOSLIB_API void app_keypair(const char *keyword,
                               int origin,
                               int nrow,
                               int ncol,
                               double *values);
  GEOSLIB_API void set_keypair_int(const char *keyword,
                                   int origin,
                                   int nrow,
                                   int ncol,
                                   int *values);
  GEOSLIB_API double get_keypone(const char *keyword, double valdef);
  GEOSLIB_API int get_keypair(const char *keyword,
                              int *nrow,
                              int *ncol,
                              double **values);
  GEOSLIB_API int get_keypair_int(const char *keyword,
                                  int *nrow,
                                  int *ncol,
                                  int **values);
  GEOSLIB_API void del_keypair(const char *keyword, int flag_exact);
  GEOSLIB_API void print_keypair(int flag_short);
  GEOSLIB_API void print_range(const char *title,
                               int ntab,
                               double *tab,
                               double *sel);
  GEOSLIB_API void ut_trace_discretize(int nseg,
                                       double *trace,
                                       double disc,
                                       int *np_arg,
                                       double **xp_arg,
                                       double **yp_arg,
                                       double **dd_arg,
                                       double **del_arg,
                                       double *dist_arg);
  GEOSLIB_API void ut_trace_sample(Db *db,
                                   ENUM_LOCS ptype,
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
  GEOSLIB_API int solve_P2(double a, double b, double c, double *x);
  GEOSLIB_API int solve_P3(double a, double b, double c, double d, double *x);
  GEOSLIB_API PL_Dist *pldist_manage(int mode,
                                     PL_Dist *pldist_loc,
                                     int ndim,
                                     int nvert);
  GEOSLIB_API double distance_point_to_segment(double x0,
                                               double y0,
                                               double x1,
                                               double y1,
                                               double x2,
                                               double y2,
                                               double *xd,
                                               double *yd,
                                               int *nint);
  GEOSLIB_API void distance_point_to_polyline(double x0,
                                              double y0,
                                              int nvert,
                                              const double *x,
                                              const double *y,
                                              PL_Dist *pldist);
  GEOSLIB_API double distance_along_polyline(PL_Dist *pldist1,
                                             PL_Dist *pldist2,
                                             double *xk,
                                             double *yl);
  GEOSLIB_API double distance_points_to_polyline(double ap,
                                                 double al,
                                                 double x1,
                                                 double y1,
                                                 double x2,
                                                 double y2,
                                                 int nvert,
                                                 double *x,
                                                 double *y);
  GEOSLIB_API int db_unfold_polyline(Db *db, int nvert, double *xl, double *yl);
  GEOSLIB_API int db_fold_polyline(Db *dbin,
                                   Db *dbout,
                                   int ncol,
                                   int *cols,
                                   int nvert,
                                   double *xl,
                                   double *yl);
  GEOSLIB_API double ut_geodetic_angular_distance(double long1,
                                                  double lat1,
                                                  double long2,
                                                  double lat2);
  GEOSLIB_API void ut_geodetic_angles(double long1,
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
  GEOSLIB_API double ut_geodetic_triangle_perimeter(double long1,
                                                    double lat1,
                                                    double long2,
                                                    double lat2,
                                                    double long3,
                                                    double lat3);
  GEOSLIB_API double ut_geodetic_triangle_surface(double long1,
                                                  double lat1,
                                                  double long2,
                                                  double lat2,
                                                  double long3,
                                                  double lat3);
  GEOSLIB_API int is_in_spherical_triangle(double *coor,
                                           double surface,
                                           double *pts1,
                                           double *pts2,
                                           double *pts3,
                                           double *wgts);
  GEOSLIB_API int is_in_spherical_triangle_optimized(double *coo0,
                                                     double *ptsa,
                                                     double *ptsb,
                                                     double *ptsc,
                                                     double *wgts);
  GEOSLIB_API double ut_distance(int ndim, double *tab1, double *tab2);
  GEOSLIB_API void ut_distance_allocated(int ndim,
                                         double **tab1,
                                         double **tab2);
  GEOSLIB_API int segment_intersect(double xd1,
                                    double yd1,
                                    double xe1,
                                    double ye1,
                                    double xd2,
                                    double yd2,
                                    double xe2,
                                    double ye2,
                                    double *xint,
                                    double *yint);
  GEOSLIB_API int ut_chebychev_coeffs(double (*func)(double,
                                                     double,
                                                     int,
                                                     double*),
                                      Cheb_Elem *cheb_elem,
                                      int nblin,
                                      double *blin);
  GEOSLIB_API int ut_chebychev_count(double (*func)(double,
                                                    double,
                                                    int,
                                                    double *),
                                     Cheb_Elem *cheb_elem,
                                     double x,
                                     int nblin,
                                     double *blin);
  GEOSLIB_API void ut_vandercorput(int n,
                                   int flag_sym,
                                   int flag_rot,
                                   int *ntri_arg,
                                   double **coor_arg);
  GEOSLIB_API int ut_icosphere(int n,
                               int flag_rot,
                               int *ntri_arg,
                               double **coor_arg);
  GEOSLIB_API void ut_shuffle_array(int nrow, int ncol, double *tab);
  GEOSLIB_API int ut_is_legendre_defined(void);
  GEOSLIB_API void define_legendre(double (*legendre_sphPlm)(int, int, double),
                                   double (*legendre_Pl)(int, double));
  GEOSLIB_API double ut_legendre(int flag_norm, int n, double v);
  GEOSLIB_API double ut_flegendre(int flag_norm, int n, int k0, double v);
  GEOSLIB_API int *ut_name_decode(const char *name,
                                  int ndim,
                                  int *nx,
                                  int verbose);
  GEOSLIB_API double *ut_rank_cells(int ndim, int *nx, int *order, int verbose);

  /*************************************/
  /* Prototyping the functions in io.c */
  /*************************************/

  //GEOSLIB_API FILE *file_open(const char *filename, int mode);
  GEOSLIB_API void record_close(void);
#ifndef SWIG
  GEOSLIB_API void redefine_message(void (*write_func)(const char *));
  GEOSLIB_API void redefine_error(void (*warn_func)(const char *));
  GEOSLIB_API void redefine_read(void (*read_func)(const char *, char *));
  GEOSLIB_API void redefine_exit(void (*exit_func)(void));
#endif
  GEOSLIB_API void constant_reset(void);
  GEOSLIB_API void constant_define(const char *name, double value);
  GEOSLIB_API void constant_print(void);
  GEOSLIB_API double constant_query(const char *name);
  GEOSLIB_API void mem_error(int nbyte);

  GEOSLIB_API void message_extern(const char *string);
  GEOSLIB_API void exit_extern();

  GEOSLIB_API void mes_process(const char *string, int ntot, int rank);
  GEOSLIB_API void string_strip_blanks(char *string, int flag_lead);
  GEOSLIB_API void string_strip_quotes(char *string);

#if defined(_WIN32) || defined(_WIN64)
  GEOSLIB_API char * strsep(char **stringp, const char* delim);
#endif
  GEOSLIB_API void tab_prints(const char *title,
                              int ncol,
                              int justify,
                              const char *string);
  GEOSLIB_API void tab_printg(const char *title,
                              int ncol,
                              int justify,
                              double value);
  GEOSLIB_API void tab_printd(const char *title,
                              int ncol,
                              int justify,
                              double value);
  GEOSLIB_API void tab_printi(const char *title,
                              int ncol,
                              int justify,
                              int value);
  GEOSLIB_API void tab_print_rowname(const char *string, int taille);
  GEOSLIB_API void tab_print_rc(const char *title,
                                int ncol,
                                int justify,
                                int mode,
                                int value);
  GEOSLIB_API void encode_printg(char *string,
                                 int ntcar,
                                 int ntdec,
                                 double value);
  GEOSLIB_API void print_current_line(void);
  GEOSLIB_API void print_matrix(const char *title,
                                int flag_limit,
                                int byrow,
                                int nx,
                                int ny,
                                const double *sel,
                                const double *tab);
  GEOSLIB_API void print_trimat(const char *title,
                                int mode,
                                int neq,
                                const double *tl);
  GEOSLIB_API void print_imatrix(const char *title,
                                 int flag_limit,
                                 int bycol,
                                 int nx,
                                 int ny,
                                 const double *sel,
                                 const int *tab);
  GEOSLIB_API void print_vector(const char *title,
                                int flag_limit,
                                int ntab,
                                const double *tab);
  GEOSLIB_API void print_names(int nx, int *ranks, VectorString names);
  GEOSLIB_API void file_dump(int ntab, double *tab);

  /*****************************************/
  /* Prototyping the functions in memory.c */
  /*****************************************/

  /* Overwriting memory management functions */

#define mem_free(tab)          mem_free_(__FILE__,__LINE__,tab)
#define mem_alloc(a,b)         mem_alloc_(__FILE__,__LINE__,a,b)
#define mem_calloc(a,b,c)      mem_calloc_(__FILE__,__LINE__,a,b,c)
#define mem_realloc(tab,a,b)   mem_realloc_(__FILE__,__LINE__,tab,a,b)
#define mem_duplicate(tab,a,b) mem_duplicate_(__FILE__,__LINE__,tab,a,b)
#define mem_copy(tab,a,b)      mem_copy_(__FILE__,__LINE__,tab,a,b)

  GEOSLIB_API void memory_leak_set(int flag);
  GEOSLIB_API void memory_leak_reset(void);
  GEOSLIB_API void memory_leak_report(void);
  GEOSLIB_API char *mem_alloc_(const char *call_file,
                               unsigned int call_line,
                               int size,
                               int flag_fatal);
  GEOSLIB_API char *mem_calloc_(const char *call_file,
                                unsigned int call_line,
                                int size_t,
                                int size,
                                int flag_fatal);
  GEOSLIB_API char *mem_realloc_(const char *call_file,
                                 unsigned int call_line,
                                 char *tab,
                                 int size,
                                 int flag_fatal);
  GEOSLIB_API char *mem_copy_(const char *call_file,
                              unsigned int call_line,
                              char *tabin,
                              int size,
                              int flag_fatal);
  GEOSLIB_API char *mem_free_(const char *call_file,
                              unsigned int call_line,
                              char *tab);
  GEOSLIB_API void mem_debug_set(int flag);
  GEOSLIB_API void memory_status(const char *title);
  GEOSLIB_API double **mem_tab_free(double **tab, int nvar);
  GEOSLIB_API double **mem_tab_alloc(int nvar, int size, int flag_fatal);
  GEOSLIB_API void time_start(void);
  GEOSLIB_API void time_reset(void);
  GEOSLIB_API void time_chunk_add(const char *call_name);
  GEOSLIB_API void time_report(void);

  /*****************************************/
  /* Prototyping the functions in matrix.c */
  /*****************************************/

  GEOSLIB_API void matrix_constant_define(int keywrd, double value);
  GEOSLIB_API double matrix_constant_query(int keywrd);
  GEOSLIB_API int matrix_get_extreme(int mode, int ntab, double *tab);
  GEOSLIB_API void matrix_invsign(int ndim, double *a);
  GEOSLIB_API int matrix_invert(double *a, int neq, int rank);
  GEOSLIB_API int matrix_invert_triangle(int neq,double *tl, int rank);
  GEOSLIB_API int matrix_invert_copy(const double *a, int neq, double *b);
  GEOSLIB_API int matrix_invsym(double *a, int neq);
  GEOSLIB_API int matrix_invgen(double *a,
                                int neq,
                                double *tabout,
                                double *cond);
  GEOSLIB_API int matrix_invsvdsym(double *a, int neq, int rank);
  GEOSLIB_API int matrix_invreal(double *mat, int neq);
  GEOSLIB_API void matrix_svd_inverse(int neq,
                                      double *s,
                                      double *u,
                                      double *v,
                                      double *tabout);
  GEOSLIB_API double matrix_determinant(int neq, const double *b);
  GEOSLIB_API int matrix_cofactor(int neq, double *a, double *b);
  GEOSLIB_API double matrix_cholesky_determinant(int neq, double *tl);
  GEOSLIB_API int matrix_eigen(const double *a,
                               int neq,
                               double *value,
                               double *vector);
  GEOSLIB_API void matrix_product(int n1,
                                  int n2,
                                  int n3,
                                  const double *v1,
                                  const double *v2,
                                  double *v3);
  GEOSLIB_API void matrix_product_safe(int n1,
                                       int n2,
                                       int n3,
                                       const double *v1,
                                       const double *v2,
                                       double *v3);
  GEOSLIB_API int matrix_prod_norme(int tranpose,
                                    int n1,
                                    int n2,
                                    const double *v1,
                                    const double *a,
                                    double *v2);
  GEOSLIB_API void matrix_transpose(int n1, int n2, double *v1, double *w1);
  GEOSLIB_API void matrix_transpose_in_place(int n1, int n2, double *v1);
  GEOSLIB_API void matrix_int_transpose_in_place(int n1, int n2, int *v1);
  GEOSLIB_API int matrix_solve(int mode,
                               const double *a,
                               const double *b,
                               double *x,
                               int neq,
                               int nrhs,
                               int *pivot);
  GEOSLIB_API double matrix_norm(double *a, int neq);
  GEOSLIB_API double matrix_normA(double *b, double *a, int neq, int subneq);
  GEOSLIB_API double inner_product(const double *a, const double *b, int neq);
  GEOSLIB_API void vector_product(double *a, double *b, double *v);
  GEOSLIB_API void vector_translate(int ndim, double *a, double *v, double *b);
  GEOSLIB_API int matrix_cholesky_decompose(const double *a, double *tl, int neq);
  GEOSLIB_API void matrix_cholesky_product(int mode,
                                           int neq,
                                           int nrhs,
                                           double *tl,
                                           double *a,
                                           double *x);
  GEOSLIB_API int matrix_cholesky_solve(int neq,
                                        double *tl,
                                        double *b,
                                        double *x);
  GEOSLIB_API int matrix_cholesky_to_invert(int neq, double *tl, double *xl);
  GEOSLIB_API void matrix_cholesky_invert(int neq, double *tl, double *xl);
  GEOSLIB_API void matrix_cholesky_norme(int mode,
                                         int neq,
                                         double *tl,
                                         double *a,
                                         double *b);
  GEOSLIB_API void matrix_triangular_product(int neq,
                                             int mode,
                                             const double *al,
                                             const double *b,
                                             double *x);
  GEOSLIB_API int is_matrix_definite_positive(int neq,
                                              const double *a,
                                              double *valpro,
                                              double *vecpro,
                                              int verbose);
  GEOSLIB_API int is_matrix_non_negative(int nrow,
                                         int ncol,
                                         double *a,
                                         int verbose);
  GEOSLIB_API int is_matrix_null(int nrow, int ncol, double *a, int verbose);
  GEOSLIB_API int is_matrix_symmetric(int neq, const double *a, int verbose);
  GEOSLIB_API int is_matrix_correlation(int neq, double *a);
  GEOSLIB_API int is_matrix_rotation(int neq, const double *a, int verbose);
  GEOSLIB_API void matrix_produit_lu(int neq, double *al, double *a);
  GEOSLIB_API VectorDouble matrix_produit_lu_VD(int neq, double *tl);
  GEOSLIB_API void matrix_set_identity(int neq, double *a);
  GEOSLIB_API int is_matrix_product_identity(int neq,
                                             double *a,
                                             double *b,
                                             double *errmax);
  GEOSLIB_API double *matrix_bind(int mode,
                                  int n11,
                                  int n12,
                                  double *a1,
                                  int n21,
                                  int n22,
                                  double *a2,
                                  int *n31,
                                  int *n32);
  GEOSLIB_API void matrix_manage(int nrows,
                                 int ncols,
                                 int nr,
                                 int nc,
                                 int *rowsel,
                                 int *colsel,
                                 double *v1,
                                 double *v2);
  GEOSLIB_API void matrix_combine(int nval,
                                  double coeffa,
                                  double *a,
                                  double coeffb,
                                  double *b,
                                  double *c);
  GEOSLIB_API void matrix_fill_symmetry(int neq, double *a);
  GEOSLIB_API double matrix_norminf(int neq, double *a);
  GEOSLIB_API double matrix_norml1(int neq, double *a);
  GEOSLIB_API void matrix_square(int neq, double *a, double *b);
  GEOSLIB_API VectorDouble matrix_square_VD(int neq,const VectorDouble& a);
  GEOSLIB_API void matrix_product_by_diag(int mode,
                                          int neq,
                                          double *a,
                                          double *c,
                                          double *b);
  GEOSLIB_API void matrix_product_by_diag_VD(int mode,
                                             int neq,
                                             VectorDouble a,
                                             const VectorDouble &c);
  GEOSLIB_API void matrix_triangle_to_square(int mode,
                                             int neq,
                                             double *tl,
                                             double *a);
  GEOSLIB_API void matrix_tri2sq(int neq, double *tl, double *a);
  GEOSLIB_API void matrix_tl2tu(int neq, const double *tl, double *tu);
  GEOSLIB_API void matrix_linear(int neq,
                                 double a1,
                                 double *a,
                                 double b1,
                                 double *b,
                                 double *x);
  GEOSLIB_API int matrix_eigen_tridiagonal(const double *vecdiag,
                                           const double *vecinf,
                                           const double *vecsup,
                                           int neq,
                                           double *eigvec,
                                           double *eigval);
  GEOSLIB_API int matrix_qo(int neq, double *hmat, double *gmat, double *xmat);
  GEOSLIB_API int matrix_qoc(int flag_invert,
                             int neq,
                             double *hmat,
                             double *gmat,
                             int na,
                             double *amat,
                             double *bmat,
                             double *xmat,
                             double *lambda);
  GEOSLIB_API int matrix_qoci(int neq,
                              double *hmat,
                              double *gmat,
                              int nae,
                              double *aemat,
                              double *bemat,
                              int nai,
                              double *aimat,
                              double *bimat,
                              double *xmat);
  GEOSLIB_API void matrix_range(int n1,
                                int n2,
                                double *v1,
                                double *mini,
                                double *maxi,
                                double *norme1,
                                double *norme2);

  /****************************************/
  /* Prototyping the functions in ascii.c */
  /****************************************/

  GEOSLIB_API void ascii_study_define(const char *study);
  GEOSLIB_API void ascii_environ_read(char *file_name, int verbose);
  GEOSLIB_API void ascii_external_filename(const char *filein,
                                           int mode,
                                           char *filename);
  GEOSLIB_API void ascii_filename(const char *type,
                                  int rank,
                                  int mode,
                                  char *filename);
  GEOSLIB_API int ascii_anam_write(const char *file_name,
                                   const Anam *anam,
                                   int verbose,
                                   int flag_calcul);
  GEOSLIB_API int ascii_frac_write(const char *file_name,
                                   Frac_Environ *frac,
                                   int verbose);
  GEOSLIB_API Db *ascii_db_read(const char *file_name,
                                int must_grid,
                                int verbose);
  GEOSLIB_API Vario *ascii_vario_read(const char *file_name, bool verbose);
  GEOSLIB_API Neigh *ascii_neigh_read(const char *file_name, int verbose);
  GEOSLIB_API Model *ascii_model_read(const char *file_name, int verbose);
  GEOSLIB_API void ascii_simu_read(char *file_name,
                                   int verbose,
                                   int *nbsimu,
                                   int *nbtuba,
                                   int *seed);
  GEOSLIB_API Rule *ascii_rule_read(const char *file_name, int verbose);
  GEOSLIB_API Anam *ascii_anam_read(const char *file_name, int verbose);
  GEOSLIB_API Frac_Environ *ascii_frac_read(const char *file_name, int verbose);
  GEOSLIB_API int ascii_option_defined(const char *file_name,
                                       int verbose,
                                       const char *option_name,
                                       int type,
                                       void *answer);
  GEOSLIB_API int ascii_table_write(const char *file_name,
                                    int verbose,
                                    int ntab,
                                    double *tab);
  GEOSLIB_API int ascii_tablei_write(const char *file_name,
                                     int verbose,
                                     int ntab,
                                     int *itab);
  GEOSLIB_API int ascii_table_read(const char *filename,
                                   int nskip,
                                   int ncol,
                                   int *nrow,
                                   double **tab_arg);

  /*****************************************/
  /* Prototyping the functions in morpho.c */
  /*****************************************/

  GEOSLIB_API int fluid_propagation(Db *dbgrid,
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
  GEOSLIB_API int fluid_extract(Db *dbgrid,
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
  GEOSLIB_API int spill_point(Db *dbgrid,
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

  GEOSLIB_API void vario_fix_codir(int ndim, VectorDouble& codir);
  GEOSLIB_API PCA *pca_free(PCA *pca);
  GEOSLIB_API PCA *pca_alloc(int nvar);
  GEOSLIB_API int pca_z2f(Db *db, PCA *pca, int flag_norm, int flag_verbose);
  GEOSLIB_API int pca_f2z(Db *db, PCA *pca, int flag_norm, int flag_verbose);

  GEOSLIB_API Vario *variogram_delete(Vario *vario);
  GEOSLIB_API void variogram_direction_del(Vario *vario, int idir);
  GEOSLIB_API int *variogram_sort(Db *db);
  GEOSLIB_API int variogram_maximum_dist1D_reached(Db *db,
                                                   int iech,
                                                   int jech,
                                                   double maxdist);
  GEOSLIB_API double variogram_maximum_distance(const Dir& dir);
  GEOSLIB_API int geometry_compute(Db *db,
                                  Vario *vario,
                                  Vario_Order *vorder,
                                  int *npair);
  GEOSLIB_API int variovect_compute(Db *db, Vario *vario, int ncomp);
  GEOSLIB_API void variogram_extension(Vario *vario,
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
  GEOSLIB_API int code_comparable(Db *db1,
                                  Db *db2,
                                  int iech,
                                  int jech,
                                  int opt_code,
                                  int tolcode);
  GEOSLIB_API int variogram_reject_pair(Db *db,
                                        int iech,
                                        int jech,
                                        double dist,
                                        double psmin,
                                        double bench,
                                        double cylrad,
                                        const VectorDouble& codir,
                                        double *ps);
  GEOSLIB_API void variogram_scale(Vario *vario, int idir);
  GEOSLIB_API int variogram_get_lag(Vario *vario,
                                    const Dir& dir,
                                    double ps,
                                    double psmin,
                                    double *dist);
  GEOSLIB_API int vario_identify_calcul_type(const String& cov_name);
  GEOSLIB_API void vardir_print(Vario *vario, int idir, int verbose);
  GEOSLIB_API void vardir_copy(Vario *vario_in,
                               int idir_in,
                               Vario *vario_out,
                               int idir_out);
  GEOSLIB_API void variogram_trans_cut(Vario *vario, int nh, double ycut);
  GEOSLIB_API int correlation_f(Db *db1,
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
                                VectorDouble& codir,
                                int flag_code,
                                int tolcode,
                                int *nindice,
                                int **indices,
                                double *correl);
  GEOSLIB_API int correlation_ident(Db *db1,
                                    Db *db2,
                                    int icol1,
                                    int icol2,
                                    Polygons *polygon);
  GEOSLIB_API int variogram_cloud_dim(Db *db, Vario *vario, double *vmax);
  GEOSLIB_API void variogram_cloud_ident(Db *db,
                                         Db *dbgrid,
                                         Vario *vario,
                                         Polygons *polygon);
  GEOSLIB_API int regression_f(Db *db1,
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
  GEOSLIB_API void condexp(Db *db1,
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
  GEOSLIB_API int vario_extract(Vario *vario,
                                int *cov_type,
                                int *ndim,
                                int *nvar,
                                int *ndir,
                                int *ndate,
                                double *scale,
                                double **dates);
  GEOSLIB_API int vardir_extract(Dir& dir,
                                 int ndim,
                                 int *flag_regular,
                                 int *npas,
                                 int *npatot,
                                 int *opt_code,
                                 int *size,
                                 double *dlag,
                                 double *toldis,
                                 double *tolang,
                                 double *slice_bench,
                                 double *cylrad,
                                 double *tolcode,
                                 double *codir,
                                 double *grincr);
  GEOSLIB_API int vario_get_rank(Vario *vario, int idir, int idate);
  GEOSLIB_API int vardir_dimension(Dir& dir);
  GEOSLIB_API void vardir_tab_extract(Vario *vario,
                                      int idir,
                                      int ivar,
                                      int jvar,
                                      int *count,
                                      int *center,
                                      int *rank,
                                      double *sw,
                                      double *gg,
                                      double *hh);
  GEOSLIB_API int maf_compute(Db *db,
                              int opt_code,
                              double tolcode,
                              VectorDouble& codir,
                              double tolang,
                              double bench,
                              double cylrad,
                              double h0,
                              double dh,
                              int verbose,
                              PCA *pca);
  GEOSLIB_API int pca_compute(Db *db, int verbose, PCA *pca);
  GEOSLIB_API int variogram_y2z(Vario *vario, Anam *anam, Model *model);

  /****************************************/
  /* Prototyping the functions in model.c */
  /****************************************/

  GEOSLIB_API void model_covtab_init(int flag_init,
                                     Model *model,
                                     double *covtab);
  GEOSLIB_API Model *model_default(int ndim, int nvar);
  GEOSLIB_API double model_calcul_basic(Model *model,
                                        int icov,
                                        int member,
                                        const VectorDouble& d1);
  GEOSLIB_API void model_calcul_cov(Model *model,
                                    CovCalcMode& mode,
                                    int flag_init,
                                    double weight,
                                    const VectorDouble& d1,
                                    double *covtab);
  GEOSLIB_API double model_calcul_cov_ij(Model *model,
                                         const CovCalcMode& mode,
                                         int ivar,
                                         int jvar,
                                         const VectorDouble& d1);
  GEOSLIB_API double model_calcul_stdev(Model *model,
                                        Db *db1,
                                        int iech1,
                                        Db *db2,
                                        int iech2,
                                        int verbose,
                                        double factor);
  GEOSLIB_API void model_calcul_cov_nostat(Model *model,
                                           CovCalcMode& mode,
                                           int flag_init,
                                           double weight,
                                           Db *db1,
                                           int iech1,
                                           Db *db2,
                                           int iech2,
                                           VectorDouble& d1,
                                           double *covtab);
  GEOSLIB_API void model_calcul_drift(Model *model,
                                      int member,
                                      Db *db,
                                      int iech,
                                      double *drftab);
  GEOSLIB_API void model_variance0(Model *model,
                                   Koption *koption,
                                   double *covtab,
                                   double *var0);
  GEOSLIB_API void model_variance0_nostat(Model *model,
                                          Koption *koption,
                                          Db *db0,
                                          int iech0,
                                          double *covtab,
                                          double *var0);
  GEOSLIB_API Model *model_free(Model *model);
  GEOSLIB_API void model_nostat_update(CovNostatInternal *cov_nostat,
                                       Model* model);
  GEOSLIB_API int model_add_cova(Model *model,
                                 int type,
                                 int flag_anisotropy,
                                 int flag_rotation,
                                 double range,
                                 double param,
                                 const VectorDouble& aniso_ranges,
                                 const VectorDouble& aniso_rotmat,
                                 const VectorDouble& coreg);
  GEOSLIB_API int model_add_drift(Model *model, int type, int rank_fex);
  GEOSLIB_API int model_add_no_property(Model *model);
  GEOSLIB_API int model_add_convolution(Model *model,
                                        int type,
                                        int idir,
                                        int count,
                                        double parameter);
  GEOSLIB_API int model_add_anamorphosis(Model *model,
                                         int anam_type,
                                         int anam_nclass,
                                         int anam_iclass,
                                         int anam_var,
                                         double anam_coefr,
                                         double anam_coefs,
                                         VectorDouble& anam_strcnt,
                                         VectorDouble& anam_stats);
  GEOSLIB_API int model_anamorphosis_set_factor(Model *model, int iclass);
  GEOSLIB_API int model_add_tapering(Model *model,
                                     int tape_type,
                                     double tape_range);
  GEOSLIB_API int model_sample(Vario *vario,
                               Model *model,
                               int flag_norm,
                               int flag_cov);
  GEOSLIB_API int model_setup(Model *model);
  GEOSLIB_API int model_fitting_sills(Vario *vario,
                                      Model *model,
                                      Option_AutoFit mauto);
  GEOSLIB_API int model_nfex(Model *model);
  GEOSLIB_API int model_update_coreg(Model *model,
                                     double *aic,
                                     double *valpro,
                                     double *vecpro);
  GEOSLIB_API int model_evaluate(Model *model,
                                 int ivar,
                                 int jvar,
                                 int rank_sel,
                                 int flag_norm,
                                 int flag_cov,
                                 int nugget_opt,
                                 int nostd,
                                 int norder,
                                 int member,
                                 int nh,
                                 VectorDouble& codir,
                                 double *h,
                                 double *g);
  GEOSLIB_API int model_evaluate_nostat(Model *model,
                                        int ivar,
                                        int jvar,
                                        int rank_sel,
                                        int flag_norm,
                                        int flag_cov,
                                        int nugget_opt,
                                        int nostd,
                                        int norder,
                                        int member,
                                        Db *db1,
                                        int iech1,
                                        Db *db2,
                                        int iech2,
                                        int nh,
                                        VectorDouble& codir,
                                        double *h,
                                        double *g);
  GEOSLIB_API int model_grid(Model *model,
                             Db *db,
                             int ivar,
                             int jvar,
                             int flag_norm,
                             int flag_cov,
                             double *g);
  GEOSLIB_API double model_cxx(Model *model,
                               Db *db1,
                               Db *db2,
                               int ivar,
                               int jvar,
                               int seed,
                               double epsdist);
  GEOSLIB_API void model_covmat(Model *model,
                                Db *db1,
                                Db *db2,
                                int ivar,
                                int jvar,
                                int flag_norm,
                                int flag_cov,
                                double *covmat);
  GEOSLIB_API double *model_covmat_by_ranks(Model *model,
                                            Db *db1,
                                            int nsize1,
                                            int *ranks1,
                                            Db *db2,
                                            int nsize2,
                                            int *ranks2,
                                            int ivar,
                                            int jvar,
                                            int flag_norm,
                                            int flag_cov);
  GEOSLIB_API int model_covmat_inchol(int verbose,
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
  GEOSLIB_API void model_covmat_nostat(Model *model,
                                       Db *db1,
                                       Db *db2,
                                       int ivar,
                                       int jvar,
                                       int flag_norm,
                                       int flag_cov,
                                       double *covmat);
  GEOSLIB_API void model_covmat_multivar(Model *model,
                                         Db *db,
                                         int flag_norm,
                                         int flag_cov,
                                         double *covmat);
  GEOSLIB_API void model_drift_mat(Model *model,
                                   int member,
                                   Db *db,
                                   double *drfmat);
  GEOSLIB_API void model_drift_vector(Model *model,
                                      int member,
                                      Db *db,
                                      int iech,
                                      double *vector);
  GEOSLIB_API void model_vector(Model *model,
                                Db *db1,
                                Db *db2,
                                int ivar,
                                int jvar,
                                int iech,
                                int flag_norm,
                                int flag_cov,
                                double *vector);
  GEOSLIB_API void model_vector_nostat(Model *model,
                                       Db *db,
                                       int ivar,
                                       int jvar,
                                       int iech,
                                       double *vector);
  GEOSLIB_API void model_vector_multivar(Model *model,
                                         Db *db,
                                         int ivar,
                                         int iech,
                                         int flag_norm,
                                         int flag_cov,
                                         double *vector);
  GEOSLIB_API void model_drift_filter(Model *model, int rank, int filter);
  GEOSLIB_API Model *model_duplicate(Model *model, double ball_radius,
                                     int mode);
//  GEOSLIB_API Model *model_modify(Model *model,
//                                  int new_nvar,
//                                  double *mean,
//                                  double *vars,
//                                  double *corr);
  GEOSLIB_API int model_stabilize(Model *model,
                                  int flag_verbose,
                                  double percent);
  GEOSLIB_API int model_normalize(Model *model, int flag_verbose);
  GEOSLIB_API void model_covupdt(Model *model,
                                 double *c0,
                                 int flag_verbose,
                                 int *flag_nugget,
                                 double *nugget);
  GEOSLIB_API double model_drift_evaluate(int verbose,
                                          Model *model,
                                          Db *db,
                                          int iech,
                                          int ivar,
                                          double *coef,
                                          double *drftab);
  GEOSLIB_API int model_is_drift_defined(Model *model, int type0);
  GEOSLIB_API Model *input_model(int ndim,
                                 int nvar,
                                 int order,
                                 int flag_sill,
                                 int flag_norm,
                                 Model *model_in);
  GEOSLIB_API int model_dimension(Model *model);
  GEOSLIB_API double model_get_field(Model *model);
  GEOSLIB_API int model_extract_cova(Model *model,
                                     int icov,
                                     int *cov_type,
                                     int *flag_aniso,
                                     double *param,
                                     VectorDouble& sill,
                                     VectorDouble& aniso_rotmat,
                                     VectorDouble& aniso_ranges);
  GEOSLIB_API void model_extract_properties(Model *model, double *tape_range);
  GEOSLIB_API void model_cova_characteristics(int rank,
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
  GEOSLIB_API double model_maximum_distance(Model *model);
  GEOSLIB_API int model_maximum_order(Model *model);
  GEOSLIB_API double model_scale2range(int type, double scale, double param);
  GEOSLIB_API double model_range2scale(int type, double range, double param);
  GEOSLIB_API double cova_get_scale_factor(int type, double param);
  GEOSLIB_API Model *model_combine(Model *model1, Model *model2, double r);
  GEOSLIB_API int model_get_nonugget_cova(Model *model);
  GEOSLIB_API int model_regularize(Model *model,
                                   Vario *vario,
                                   Db *db,
                                   int opt_norm,
                                   double nug_ratio);
  GEOSLIB_API double constraints_get(const Constraints& constraints,
                                     int type,
                                     int igrf,
                                     int istr,
                                     int elem,
                                     int v1,
                                     int v2);
  GEOSLIB_API void constraints_print(const Constraints& constraints);
  GEOSLIB_API int modify_constraints_on_sill(Constraints& constraints);
  GEOSLIB_API void fill_external_cov_model(External_Cov& E_Cov);

  /****************************************/
  /* Prototyping the functions in neigh.c */
  /****************************************/

  GEOSLIB_API int neigh_start(Db *dbin, Neigh *neigh);
  GEOSLIB_API void neigh_stop(void);
  GEOSLIB_API int neigh_select(Db *dbin,
                               Db *dbout,
                               int iech_out,
                               Neigh *neigh,
                               int flag_simu,
                               int *nech,
                               int *rank);
  GEOSLIB_API Neigh *neigh_free(Neigh *neigh);
  GEOSLIB_API Neigh *neigh_init_bench(int ndim, int flag_xvalid, double width);
  GEOSLIB_API Neigh *neigh_init_unique(int ndim);
  GEOSLIB_API Neigh *neigh_init_image(int ndim,
                                      int flag_xvalid,
                                      int skip,
                                      const VectorDouble& nbgh_image = VectorDouble());
  GEOSLIB_API Neigh *neigh_init(int ndim,
                                int type,
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
                                const VectorDouble& nbgh_radius = VectorDouble(),
                                const VectorDouble& nbgh_rotmat = VectorDouble(),
                                const VectorDouble& nbgh_image  = VectorDouble());
  GEOSLIB_API void neigh_print(const Neigh *neigh);
  GEOSLIB_API void neigh_echo(Db *dbin,
                              Neigh *neigh,
                              int *rank,
                              int nsel,
                              double *tab);
  GEOSLIB_API int neigh_extract(Neigh *neigh,
                                int *type,
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
                                VectorDouble& nbgh_rotmat,
                                VectorDouble& nbgh_radius,
                                VectorDouble& nbgh_image);
  GEOSLIB_API int *neigh_calc(Db *dbin,
                              Model *model,
                              Neigh *neigh,
                              double *target,
                              int *nech_out);
  GEOSLIB_API double neigh_continuous_variance(Neigh *neigh,
                                               Db *db1,
                                               int rank1,
                                               Db *db2,
                                               int rankZ);

  /***************************************/
  /* Prototyping the functions in anam.c */
  /***************************************/

  GEOSLIB_API double anam_y2z(Anam *anam, double y, int flag_bound);
  GEOSLIB_API void anam_update_hermitian(AnamHermite *anam_hermite,
                                         int nh,
                                         double pymin,
                                         double pzmin,
                                         double pymax,
                                         double pzmax,
                                         double aymin,
                                         double azmin,
                                         double aymax,
                                         double azmax,
                                         double r,
                                         double variance,
                                         const VectorDouble& psi_hn);
 GEOSLIB_API void anam_update_discrete_DD(AnamDiscreteDD   *anam_discrete_DD,
                                          int     ncut,
                                          double  scoef,
                                          double  mu,
                                          const VectorDouble& zcut,
                                          const VectorDouble& pcaz2f,
                                          const VectorDouble& pcaf2z,
                                          const VectorDouble& stats);
  GEOSLIB_API void anam_update_empirical(AnamEmpirical *anam_empirical,
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
                                        const VectorDouble& tdisc);
  GEOSLIB_API void anam_update_discrete_IR(AnamDiscreteIR *anam_discrste_IR,
                                           int ncut,
                                           double s,
                                           const VectorDouble& zcut,
                                           const VectorDouble& stats);
  GEOSLIB_API int anam_discrete_z2factor(Anam *anam,
                                         Db *db,
                                         int nfact,
                                         const VectorInt& ifacs);
  GEOSLIB_API int anam_point_to_block(Anam *anam,
                                      int verbose,
                                      double cvv,
                                      double coeff,
                                      double mu);
  GEOSLIB_API double ce_compute_Z2(double krigest,
                                   double krigstd,
                                   const VectorDouble& phis);
  GEOSLIB_API int anam_factor2qt(Db *db,
                                 Anam *anam,
                                 int ncutmine,
                                 double *cutmine,
                                 double z_max,
                                 int flag_correct,
                                 int nb_est,
                                 int* cols_est,
                                 int nb_std,
                                 int* cols_std,
                                 int ncode,
                                 int *codes,
                                 int *ncut,
                                 int *qt_vars);
  GEOSLIB_API void selectivity_interpolate(int verbose,
                                           double *zcut,
                                           int nclass,
                                           double *calest,
                                           int ncut,
                                           double *calcut);
  GEOSLIB_API int anam_get_r(Anam *anam, double cvv, double mu, double *r);
  GEOSLIB_API int anam_vario_z2y(Anam *anam, double cvv, Vario *vario);

  GEOSLIB_API int uc_f(Db *db,
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
  GEOSLIB_API int ce_f(Db *db,
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

  GEOSLIB_API void grid_iterator_init(GridC *grid,
                                      const VectorInt& order = VectorInt());
  GEOSLIB_API VectorInt grid_iterator_next(GridC *grid);

  GEOSLIB_API int *db_indg_alloc(Db *db);
  GEOSLIB_API int *db_indg_free(int *indice);
  GEOSLIB_API double *db_sample_free(double *tab);
  GEOSLIB_API double *db_sample_alloc(Db *db, ENUM_LOCS locatorType);
  GEOSLIB_API int db_sample_load(Db *db, ENUM_LOCS locatorType, int iech, double *tab);
  GEOSLIB_API double *db_vector_free(double *tab);
  GEOSLIB_API double *db_vector_alloc(Db *db);
  GEOSLIB_API int db_selection_get(Db *db, int item, double *tab);
  GEOSLIB_API int db_vector_get(Db *db, ENUM_LOCS locatorType, int locatorIndex, double *tab);
  GEOSLIB_API int db_vector_put(Db *db, ENUM_LOCS locatorType, int locatorIndex, double *tab);
  GEOSLIB_API int db_vector_get_att_sel_compress(Db *db,
                                                 int icol,
                                                 int *number,
                                                 double *tab);
  GEOSLIB_API int db_vector_get_att(Db *db, int iatt, double *tab);
  GEOSLIB_API int db_vector_get_att_sel(Db *db, int iatt, double *tab);
  GEOSLIB_API int db_name_set(Db *db, int iatt, const String& name);
  GEOSLIB_API String db_name_get_by_att(const Db *db, int iatt);
  GEOSLIB_API String db_name_get_by_col(Db *db, int icol);
  GEOSLIB_API int db_name_identify(Db *db, const String& string);
  GEOSLIB_API void db_attribute_del_mult(Db *db, int i_del, int n_del);
  GEOSLIB_API void db_attribute_init(Db *db,
                                     int ncol,
                                     int iatt,
                                     double valinit);
  GEOSLIB_API void db_attribute_copy(Db *db, int iatt_in, int iatt_out);
  GEOSLIB_API int db_attribute_identify(Db *db, ENUM_LOCS locatorType, int locatorIndex);
  GEOSLIB_API int db_sample_get_att(Db *db,
                                    int iech,
                                    int number,
                                    int iatt,
                                    double *tab);
  GEOSLIB_API void db_sample_put_att(Db *db,
                                     int iech,
                                     int number,
                                     int iatt,
                                     double *tab);
  GEOSLIB_API int db_locator_attribute_add(Db *db,
                                           ENUM_LOCS locatorType,
                                           int number,
                                           int r_tem,
                                           double valinit,
                                           int *iptr);
  GEOSLIB_API void db_locators_correct(VectorString& strings,
                                       const VectorInt& current,
                                       int flag_locnew);
  GEOSLIB_API int db_coorvec_put(Db *db, int idim, double *tab);
  GEOSLIB_API int db_coorvec_get(Db *db, int idim, double *tab);
  GEOSLIB_API Db *db_delete(Db *db);
  GEOSLIB_API int db_grid_match(Db *db1, Db *db2);
  GEOSLIB_API int db_is_isotropic(Db *db, int iech, double *data);
  GEOSLIB_API void db_grid_print(Db *db);

  GEOSLIB_API Db *db_create_grid_multiple(Db *dbin, int *nmult, int flag_add_rank);
  GEOSLIB_API Db *db_create_grid_divider(Db *dbin, int *nmult, int flag_add_rank);
  GEOSLIB_API Db *db_create_grid_dilate(Db *dbin, int mode, int *nshift,
                                        int flag_add_rank);
  GEOSLIB_API Db *db_grid_sample(Db *dbin, int *nmult);
  GEOSLIB_API int db_grid_define_coordinates(Db *db);
  GEOSLIB_API Db *db_create_from_target(double *target, int ndim, int flag_add_rank);
  GEOSLIB_API void db_sample_print(Db *db,
                                   int iech,
                                   int flag_ndim,
                                   int flag_nvar,
                                   int flag_nerr);
  GEOSLIB_API int db_center(Db *db, double *center);
  GEOSLIB_API int db_extension(Db *db,
                               double *mini,
                               double *maxi,
                               double *delta);
  GEOSLIB_API int db_extension_rotated(Db *db,
                                       double *rotmat,
                                       double *mini_arg,
                                       double *maxi_arg,
                                       double *delta_arg);
  GEOSLIB_API int db_attribute_range(Db *db,
                                     int icol,
                                     double *mini,
                                     double *maxi,
                                     double *delta);
  GEOSLIB_API int db_extension_diag(Db *db, double *diag);
  GEOSLIB_API double db_epsilon_distance(Db *db);
  GEOSLIB_API int db_index_grid_to_sample(Db *db, const int *indg);
  GEOSLIB_API void db_index_sample_to_grid(Db *db, int iech, int *indg);
  GEOSLIB_API int db_index_sorted_in_grid(Db *db, int iech, int *indg);
  GEOSLIB_API int db_selref(int ndim,
                            int *nx,
                            int *ref,
                            double *tabin,
                            double *tabout);
  GEOSLIB_API Db *db_extract(Db *db, int *ranks);
  GEOSLIB_API Db *db_regularize(Db *db, Db *dbgrid, int flag_center);
  GEOSLIB_API int compat_NDIM(Db *db1, Db *db2);
  GEOSLIB_API int same_mesh(Db *db1, Db *db2);
  GEOSLIB_API int same_rotation(Db *db1, Db *db2);
  GEOSLIB_API int get_NECH(const Db *db);
  GEOSLIB_API double get_ARRAY(const Db *db, int iech, int iatt);
  GEOSLIB_API double get_IDIM(Db *db, int iech, int idim);
  GEOSLIB_API double get_grid_IDIM(Db *db, int iech, int idim);
  GEOSLIB_API int    match_domain_ref(double value);
  GEOSLIB_API int    get_DOMAIN(Db *db, int iech);
  GEOSLIB_API void   domain_ref_define(int value, int verbose);
  GEOSLIB_API int    domain_ref_query(void);
  GEOSLIB_API void   domain_ref_print(void);
  GEOSLIB_API double get_grid_value(Db *dbgrid,
                                    int iptr,
                                    int *indg,
                                    int ix,
                                    int iy,
                                    int iz);
  GEOSLIB_API void set_grid_value(Db *dbgrid,
                                  int iptr,
                                  int *indg,
                                  int ix,
                                  int iy,
                                  int iz,
                                  double value);
  GEOSLIB_API int get_LOCATOR_NITEM(Db *db, ENUM_LOCS locatorType);
  GEOSLIB_API int exist_LOCATOR(Db *db, ENUM_LOCS locatorType);
  GEOSLIB_API double get_LOCATOR_ITEM(Db *db, ENUM_LOCS locatorType, int locatorIndex, int iech);
  GEOSLIB_API void set_LOCATOR_ITEM(Db *db,
                                    ENUM_LOCS locatorType,
                                    int locatorIndex,
                                    int iech,
                                    double value);
  GEOSLIB_API int db_get_rank_absolute_to_relative(Db *db, int iech0);
  GEOSLIB_API int db_get_rank_relative_to_absolute(Db *db, int iech0);
  GEOSLIB_API int is_grid(const Db *db, bool verbose=false);
  GEOSLIB_API int is_grid_multiple(Db *db1, Db *db2);
  GEOSLIB_API void get_grid_multiple(Db *db1,
                                     int *nmult,
                                     int flag_cell,
                                     VectorInt& nx,
                                     VectorDouble& dx,
                                     VectorDouble& x0);
  GEOSLIB_API void get_grid_divider(Db *db1,
                                    int *nmult,
                                    int flag_cell,
                                    VectorInt& nx,
                                    VectorDouble& dx,
                                    VectorDouble& x0);
  GEOSLIB_API int get_grid_dilate(Db *db,
                                  int mode,
                                  int *nshift,
                                  VectorInt& nx,
                                  VectorDouble& dx,
                                  VectorDouble& x0);
  GEOSLIB_API int db_grid_copy_params(Db *dbin, int mode, Db *dbout);
  GEOSLIB_API Db *db_grid_reduce(Db *db_grid,
                                 int iptr,
                                 int *margin,
                                 int *limmin,
                                 int flag_sel,
                                 int flag_copy,
                                 int verbose,
                                 double vmin,
                                 double vmax);
  GEOSLIB_API double distance_inter(Db *db1,
                                    Db *db2,
                                    int iech1,
                                    int iech2,
                                    double *dist_vect);
  GEOSLIB_API double distance_intra(Db *db,
                                    int iech1,
                                    int iech2,
                                    double *dist_vect);
  GEOSLIB_API double distance_grid(Db *db,
                                   int flag_moins1,
                                   int iech1,
                                   int iech2,
                                   double *dist_vect);
  GEOSLIB_API double *db_distances_general(Db *db1,
                                           Db *db2,
                                           int niso,
                                           int mode,
                                           int flag_same,
                                           int *n1,
                                           int *n2,
                                           double *dmin,
                                           double *dmax);
  GEOSLIB_API double cosdir(Db *db, int iech1, int iech2,
                            const VectorDouble& codir);
  GEOSLIB_API double bench_distance(Db *db, int iech1, int iech2);
  GEOSLIB_API double cylinder_radius(Db *db,
                                     int iech1,
                                     int iech2,
                                     const VectorDouble& codir);
  GEOSLIB_API double db_grid_maille(Db *db);
  GEOSLIB_API int point_to_grid(Db *db,
                                double *coor,
                                int flag_expand,
                                int *indg);
  GEOSLIB_API int point_to_bench(Db *db,
                                 double *coor,
                                 int flag_outside,
                                 int *indb);
  GEOSLIB_API void grid_to_point(Db *db,
                                 int *indg,
                                 double *percent,
                                 double *coor);
  GEOSLIB_API int index_point_to_grid(Db *db,
                                      int iech,
                                      int flag_expand,
                                      Db *dbout,
                                      double *coor);
  GEOSLIB_API int point_to_point(Db *db, double *coor);
  GEOSLIB_API int point_inside_grid(Db *db, int iech, Db *dbgrid);
  GEOSLIB_API int migrate_grid_to_coor(Db *db_grid,
                                       int iv_grid,
                                       int np,
                                       double *xp,
                                       double *yp,
                                       double *zp,
                                       double *tab);
  GEOSLIB_API int expand_point_to_coor(Db *db1,
                                       int iatt,
                                       int np,
                                       double *xp,
                                       double *yp,
                                       double *zp,
                                       double *tab);
  GEOSLIB_API int expand_point_to_grid(Db *db_point,
                                       Db *db_grid,
                                       int iatt,
                                       int iatt_time,
                                       int iatt_angle,
                                       int iatt_scaleu,
                                       int iatt_scalev,
                                       int iatt_scalew,
                                       int flag_index,
                                       int ldmax,
                                       const VectorDouble& dmax,
                                       VectorDouble& tab);
  GEOSLIB_API int db_center_point_to_grid(Db *db_point,
                                          Db *db_grid,
                                          double eps_random);
  GEOSLIB_API int interpolate_variable_to_point(Db *db_grid,
                                                int iatt,
                                                int np,
                                                double *xp,
                                                double *yp,
                                                double *zp,
                                                double *tab);
  GEOSLIB_API int points_to_block(Db *dbpoint,
                                  Db *dbgrid,
                                  int option,
                                  int flag_size,
                                  int iatt_time,
                                  int iatt_size,
                                  int iatt_angle,
                                  int iatt_scaleu,
                                  int iatt_scalev,
                                  int iatt_scalew);
  GEOSLIB_API int db_gradient_components(Db *dbgrid);
  GEOSLIB_API int db_streamline(Db *dbgrid,
                                Db *dbpoint,
                                int niter,
                                double step,
                                int flag_norm,
                                int use_grad,
                                int save_grad,
                                int *nbline_loc,
                                int *npline_loc,
                                double **line_loc);
  GEOSLIB_API int manage_external_info(int mode,
                                       ENUM_LOCS locatorType,
                                       Db *dbin,
                                       Db *dbout,
                                       int *istart);
  GEOSLIB_API int db_locate_in_grid(Db *dbgrid, double *coor);
  GEOSLIB_API void db_monostat(Db *db,
                               int ivar,
                               double *wtot,
                               double *mean,
                               double *var,
                               double *mini,
                               double *maxi);
  GEOSLIB_API int db_normalize(Db *db,
                               const char *oper,
                               int ncol,
                               int *cols,
                               double center,
                               double stdv);
  GEOSLIB_API int db_duplicate(Db *db1,
                               Db *db2,
                               int flag_same,
                               int flag_print,
                               int opt_code,
                               double tolcode,
                               double *dist,
                               double *sel);
  GEOSLIB_API int db_gradient_update(Db *db);
  GEOSLIB_API int surface(Db *db_point,
                          Db *db_grid,
                          int icol,
                          double dlim,
                          double *dtab,
                          double *gtab);
  GEOSLIB_API int db_edit(Db *db, int *flag_valid);
  GEOSLIB_API int db_grid_copy(Db *db1,
                               Db *db2,
                               int *ind1,
                               int *ind2,
                               int ncol,
                               int *cols);
  GEOSLIB_API int db_grid_copy_dilate(Db *db1,
                                      int iatt1,
                                      Db *db2,
                                      int iatt2,
                                      int mode,
                                      int *nshift);
  GEOSLIB_API int db_proportion(Db *db,
                                Db *dbgrid,
                                int nfac1max,
                                int nfac2max,
                                int *nclout);
  GEOSLIB_API int db_merge(Db *db, int ncol, int *cols);
  GEOSLIB_API int db_count_defined(Db *db, int icol);

  GEOSLIB_API int db_prop_read(Db *db, int ix, int iy, double *props);
  GEOSLIB_API int db_prop_write(Db *db, int ix, int iy, double *props);
  GEOSLIB_API int db_resind(Db *db, int ivar, int ncut, double *zcut);
  GEOSLIB_API int db_prop_thresh(Db *db,
                                 Db *dbprop,
                                 Rule *rule,
                                 Model *model,
                                 const VectorDouble& propcst,
                                 int flag_stat);
  GEOSLIB_API int db_gradient_modang_to_component(Db *db,
                                                  int ang_conv,
                                                  int iad_mod,
                                                  int iad_ang,
                                                  int iad_gx,
                                                  int iad_gy);
  GEOSLIB_API int db_gradient_component_to_modang(Db *db,
                                                  int verbose,
                                                  int iad_gx,
                                                  int iad_gy,
                                                  int iad_mod,
                                                  int iad_ang,
                                                  double scale,
                                                  double ve);
  GEOSLIB_API int db_compositional_transform(Db *db,
                                             int verbose,
                                             int mode,
                                             int type,
                                             int number,
                                             int *iatt_in,
                                             int *iatt_out,
                                             int *numout);
  GEOSLIB_API Db *db_point_init(int mode,
                                int verbose,
                                int ndim,
                                int seed,
                                double density,
                                double range,
                                double beta,
                                Db *dbgrid,
                                const VectorDouble& origin,
                                const VectorDouble& extend);
  GEOSLIB_API int db_smooth_vpc(Db *db, int width, double range);
  GEOSLIB_API double *db_grid_sampling(Db *dbgrid,
                                       double *x1,
                                       double *x2,
                                       int ndisc,
                                       int ncut,
                                       double *cuts,
                                       int *nval_ret);
  GEOSLIB_API int db_grid2point_sampling(Db *dbgrid,
                                         int nvar,
                                         int *vars,
                                         int *npacks,
                                         int npcell,
                                         int nmini,
                                         int *nech,
                                         double **coor,
                                         double **data);
  GEOSLIB_API int db_grid_patch(Db *ss_grid,
                                Db *db_grid,
                                int iptr_ss,
                                int iptr_db,
                                int iptr_rank,
                                int new_rank,
                                int oper,
                                int verbose);
  GEOSLIB_API int db_polygon_distance(Db *db,
                                      Polygons *polygon,
                                      double dmax,
                                      int scale,
                                      int polin);

  /****************************************/
  /* Prototyping the functions in stats.c */
  /****************************************/

  GEOSLIB_API int stats_point_to_grid(Db *dbgrid,
                                      Db *db,
                                      const char *oper,
                                      int ivar,
                                      int jvar,
                                      int ncut,
                                      double *cuts,
                                      double *tab);
  GEOSLIB_API int db_stats(Db *db,
                           const String& oper,
                           const VectorInt& cols,
                           int flag_mono,
                           int flag_verbose,
                           double *resta);
  GEOSLIB_API int db_stats_grid(Db *db,
                                Db *dbgrid,
                                const char *oper,
                                int ncol,
                                int *cols,
                                int radius);
  GEOSLIB_API int stats_proportion(Db *dbin,
                                   Db *dbout,
                                   int pos,
                                   int nfacies,
                                   int radius);
  GEOSLIB_API int stats_transition(Db *dbin,
                                   Db *dbout,
                                   int pos,
                                   int nfacies,
                                   int radius,
                                   int orient);
  GEOSLIB_API int stats_residuals(int verbose,
                                  int nech,
                                  double *tab,
                                  int ncut,
                                  double *zcut,
                                  int *nsorted,
                                  double *mean,
                                  double *residuals,
                                  double *T,
                                  double *Q);
  GEOSLIB_API int db_upscale(Db *dbgrid1, Db *dbgrid2, int orient, int verbose);
  GEOSLIB_API int db_diffusion(Db *dbgrid1,
                               Db *dbgrid2,
                               int orient,
                               int niter,
                               int nseed,
                               int seed,
                               int verbose);

  /***************************************/
  /* Prototyping the functions in skin.c */
  /***************************************/

  GEOSLIB_API Skin *skin_define(Db *db,
                                int (*func_already_filled)(int ipos),
                                int (*func_to_be_filled)(int ipos),
                                double (*func_get_weight)(int ipos, int idir));
  GEOSLIB_API Skin *skin_undefine(Skin *skin);
  GEOSLIB_API void skin_print(Skin *skin);
  GEOSLIB_API int  skin_init(Skin *skin, int verbose);
  GEOSLIB_API int  skin_remains(Skin *skin);
  GEOSLIB_API void skin_next(Skin *skin, int *rank, int *ipos);
  GEOSLIB_API int skin_unstack(Skin *skin, int rank, int ipos);
  GEOSLIB_API int skin_grid_shift(Skin *skin, int lec, int dir, int *iad);

  /******************************************/
  /* Prototyping the functions in spatial.c */
  /******************************************/

  GEOSLIB_API int cgi(Db *db,
                      int ivar,
                      double *center,
                      double *mvalue,
                      double *mvector,
                      double *inertia,
                      double *wztot);
  GEOSLIB_API int spatial(Db *db, double *totab, double *parea, double *eqarea);

  /******************************************/
  /* Prototyping the functions in convert.c */
  /******************************************/

  GEOSLIB_API int db_grid_read_zycor1(const char *filename,
                                      int verbose,
                                      int *nx,
                                      double *x0,
                                      double *dx);
  GEOSLIB_API int db_grid_read_zycor2(const char *filename,
                                      int *nx,
                                      double *x0,
                                      double *dx,
                                      double *tab);
  GEOSLIB_API int db_grid_read_bmp1(const char *filename,
                                    int verbose,
                                    int *nx,
                                    double *x0,
                                    double *dx);
  GEOSLIB_API int db_grid_read_bmp2(const char *filename,
                                    int *nx,
                                    double *x0,
                                    double *dx,
                                    double *tab);
  GEOSLIB_API int db_grid_read_prop1(const char *filename,
                                     int verbose,
                                     int *ncol,
                                     int *nx,
                                     double *x0,
                                     double *dx);
  GEOSLIB_API int db_grid_read_prop2(const char *filename,
                                     int ncol_r,
                                     int *nx_r,
                                     double *x0_r,
                                     double *dx_r,
                                     double *tab);
  GEOSLIB_API int db_grid_read_f2g(const char *filename,
                                   int verbose,
                                   int nx[2],
                                   double x0[3],
                                   double dx[3],
                                   double *angle,
                                   int *ncol,
                                   double **tab_arg);
  GEOSLIB_API int db_grid_write_zycor(const char *filename, Db *db, int icol);
  GEOSLIB_API int db_grid_write_XYZ(const char *filename, Db *db, int icol);
  GEOSLIB_API int db_write_vtk(const char *filename,
                               Db *db,
                               const VectorInt& cols,
                               const VectorString& names);
  GEOSLIB_API int db_grid_write_bmp(const char *filename,
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
  GEOSLIB_API int db_grid_write_irap(const char *filename,
                                     Db *db,
                                     int icol,
                                     int nsamplex,
                                     int nsampley);
  GEOSLIB_API int db_grid_write_prop(const char *filename,
                                     Db *db,
                                     int ncol,
                                     int *icols);
  GEOSLIB_API int db_grid_write_eclipse(const char *filename, Db *db, int icol);
  GEOSLIB_API int db_well_read_las(const char *filename,
                                   int verbose,
                                   double xwell,
                                   double ywell,
                                   double cwell,
                                   int *nvarout,
                                   int *nechout,
                                   char ***var_names,
                                   double **tab);
  GEOSLIB_API int csv_table_read(const char *filename,
                                 int verbose,
                                 int flag_header,
                                 int nskip,
                                 const char *char_sep,
                                 const char *char_dec,
                                 const char *na_string,
                                 int ncol_max,
                                 int nrow_max,
                                 int *ncol_arg,
                                 int *nrow_arg,
                                 VectorString& names,
                                 VectorDouble& tab);

  /****************************************/
  /* Prototyping the functions in krige.c */
  /****************************************/

  GEOSLIB_API int is_flag_data_disc_defined(void);
  GEOSLIB_API int krige_koption_manage(int mode,
                                       int flag_check,
                                       int calcul,
                                       int flag_rand,
                                       VectorInt ndisc);
  GEOSLIB_API void krige_lhs_print(int nech,
                                   int neq,
                                   int nred,
                                   int *flag,
                                   double *lhs);
  GEOSLIB_API void krige_rhs_print(int nvar,
                                   int nech,
                                   int neq,
                                   int nred,
                                   int *flag,
                                   double *rhs);
  GEOSLIB_API void krige_dual_print(int nech,
                                    int neq,
                                    int nred,
                                    int *flag,
                                    double *dual);
  GEOSLIB_API int bayes_simulate(Model *model,
                                 int nbsimu,
                                 double *rmean,
                                 double *rcov,
                                 double *smean);
  GEOSLIB_API int krimage_func(Db *dbgrid, Model *model, Neigh *neigh);
  GEOSLIB_API int image_smoother(Db *dbgrid,
                                 Neigh *neigh,
                                 int type,
                                 double range);
  GEOSLIB_API int krigdgm_f(Db *dbin,
                            Db *dbout,
                            Model *model,
                            Neigh *neigh,
                            int flag_est,
                            int flag_std,
                            int flag_varz,
                            double rval);
  GEOSLIB_API int krigcell_f(Db *dbin,
                             Db *dbout,
                             Model *model,
                             Neigh *neigh,
                             VectorInt ndisc,
                             int flag_est,
                             int flag_std,
                             VectorInt rank_colcok);
  GEOSLIB_API int kriggam_f(Db *dbin,
                            Db *dbout,
                            Anam *anam,
                            Model *model,
                            Neigh *neigh);
  GEOSLIB_API int krigprof_f(Db *dbin,
                             Db *dbout,
                             Model *model,
                             Neigh *neigh,
                             int ncode,
                             int flag_est,
                             int flag_std);
  GEOSLIB_API int kribayes_f(Db *dbin,
                             Db *dbout,
                             Model *model,
                             Neigh *neigh,
                             double *dmean,
                             double *dcov,
                             int flag_est,
                             int flag_std);
  GEOSLIB_API int krigsum_f(Db *dbin,
                            Db *dbout,
                            Model *model,
                            Neigh *neigh,
                            int flag_positive);
  GEOSLIB_API int krigmvp_f(Db *dbin,
                            Db *db3grid,
                            Db *db2grid,
                            int fsum,
                            Model *model,
                            Neigh *neigh);

  GEOSLIB_API int krigtest_dimension(Db *dbin,
                                     Db *dbout,
                                     Model *model,
                                     Neigh *neigh,
                                     int iech0,
                                     int calcul,
                                     VectorInt ndisc,
                                     int *ndim_ret,
                                     int *nech_ret,
                                     int *neq_ret,
                                     int *nrhs_ret);
  GEOSLIB_API int krigtest_f(Db *dbin,
                             Db *dbout,
                             Model *model,
                             Neigh *neigh,
                             int iech0,
                             int calcul,
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
  GEOSLIB_API int krigsampling_f(Db *dbin,
                                 Db *dbout,
                                 Model *model,
                                 double beta,
                                 int nsize1,
                                 int *ranks1,
                                 int nsize2,
                                 int *ranks2,
                                 int flag_std,
                                 int verbose);
  GEOSLIB_API int dk_f(Db *dbin,
                       Db *dbsmu,
                       Model *model,
                       Neigh *neigh,
                       int nfactor,
                       VectorInt nmult,
                       VectorInt ndisc,
                       int flag_est,
                       int flag_std);
  GEOSLIB_API int global_arithmetic(Db *dbin,
                                    Db *dbgrid,
                                    Model *model,
                                    int ivar,
                                    int flag_verbose,
                                    int seed,
                                    double surface,
                                    double *zest,
                                    double *sse,
                                    double *cvgeo);
  GEOSLIB_API int global_kriging(Db *dbin,
                                 Db *dbout,
                                 Model *model,
                                 int ivar,
                                 int flag_verbose,
                                 int calcul,
                                 int seed,
                                 double surface,
                                 double *zest,
                                 double *sse,
                                 double *cvgeo,
                                 double *weights);
  GEOSLIB_API int global_transitive(Db *dbgrid,
                                    Model *model,
                                    int flag_verbose,
                                    int flag_regular,
                                    int ndisc,
                                    double *zest,
                                    double *cve,
                                    double *cvtrans);
  GEOSLIB_API int invdist_f(Db *dbin,
                            Db *dbout,
                            int exponent,
                            int flag_expand,
                            double dmax);
  GEOSLIB_API int anakexp_f(Db *db,
                            double *covdd,
                            double *covd0,
                            double top,
                            double bot,
                            int ncov_radius,
                            int neigh_radius,
                            int flag_sym,
                            int nfeq);
  GEOSLIB_API int anakexp_3D(Db *db,
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
  GEOSLIB_API int sampling_f(Db *db,
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
  GEOSLIB_API int declustering_f(Db *db,
                                 Model *model,
                                 Neigh *neigh,
                                 Db *dbgrid,
                                 int method,
                                 double *radius,
                                 VectorInt ndisc,
                                 int flag_sel,
                                 int verbose);
  GEOSLIB_API int inhomogeneous_kriging(Db *dbdat,
                                        Db *dbsrc,
                                        Db *dbout,
                                        double power,
                                        int flag_source,
                                        Model *model_dat,
                                        Model *model_src);
  GEOSLIB_API void fill_external_cov_kriging(External_Cov& E_Cov);

  /*****************************************/
  /* Prototyping the functions in simtub.c */
  /*****************************************/

  GEOSLIB_API void simu_define_func_transf(void (*st_simu_transf)(Db *,
                                                                  int,
                                                                  int,
                                                                  int));
  GEOSLIB_API void simu_define_func_update(void (*st_simu_update)(Db *,
                                                                  int,
                                                                  int,
                                                                  int));
  GEOSLIB_API void simu_define_func_scale(void (*st_simu_scale)(Db *,
                                                                int,
                                                                int));
  GEOSLIB_API void simu_func_categorical_transf(Db *db,
                                                int verbose,
                                                int isimu,
                                                int nbsimu);
  GEOSLIB_API void simu_func_continuous_update(Db *db,
                                               int verbose,
                                               int isimu,
                                               int nbsimu);
  GEOSLIB_API void simu_func_categorical_update(Db *db,
                                                int verbose,
                                                int isimu0,
                                                int nbsimu0);
  GEOSLIB_API void simu_func_continuous_scale(Db *db, int verbose, int nbsimu);
  GEOSLIB_API void simu_func_categorical_scale(Db *db, int verbose, int nbsimu);

  GEOSLIB_API int get_rank_from_propdef(Props *propdef, int ipgs, int igrf);
  GEOSLIB_API void check_mandatory_attribute(const char *method,
                                             Db *db,
                                             ENUM_LOCS locatorType);
  GEOSLIB_API int simtub_workable(Model *model);
  GEOSLIB_API int simdgm(Db *dbin,
                         Db *dbout,
                         Model *model,
                         Neigh *neigh,
                         double rval,
                         int seed,
                         int nbsimu,
                         int nbtuba,
                         int flag_check);
  GEOSLIB_API int simcond(Db *dbin,
                          Db *dbout,
                          Model *model,
                          int seed,
                          int nbsimu,
                          int nbtuba,
                          int nboot,
                          int niter,
                          double gibbs_eps,
                          int flag_check,
                          int flag_ce,
                          int flag_cstd,
                          int verbose);
  GEOSLIB_API int simbayes(Db *dbin,
                           Db *dbout,
                           Model *model,
                           Neigh *neigh,
                           double *dmean,
                           double *dcov,
                           int seed,
                           int nbsimu,
                           int nbtuba,
                           int flag_check);
  GEOSLIB_API int simbipgs(Db *dbin,
                           Db *dbout,
                           Db *dbprop,
                           Rule *rule1,
                           Rule *rule2,
                           Model *model11,
                           Model *model12,
                           Model *model21,
                           Model *model22,
                           Neigh *neigh,
                           const VectorDouble& propcst,
                           int flag_stat,
                           int flag_gaus,
                           int flag_prop,
                           int flag_check,
                           int flag_show,
                           int nfac1,
                           int nfac2,
                           int seed,
                           int nbsimu,
                           int nbtuba,
                           int nboot,
                           int niter,
                           double percent,
                           double gibbs_eps);
  GEOSLIB_API int simmaxstable(Db *dbout,
                               Model *model,
                               double ratio,
                               int seed,
                               int nbtuba,
                               int flag_simu,
                               int flag_rank,
                               int verbose);
  GEOSLIB_API int simtub_potential(Db *dbiso,
                                   Db *dbgrd,
                                   Db *dbtgt,
                                   Db *dbout,
                                   Model *model,
                                   int nbsimu,
                                   int nbtuba,
                                   double delta);
  GEOSLIB_API int simRI(Db *dbout,
                        Model *model,
                        int ncut,
                        double *zcut,
                        double *wcut,
                        int seed,
                        int nbtuba,
                        int verbose);
  GEOSLIB_API int gibbs_sampler(Db *db,
                                Model *model,
                                int nbsimu,
                                int seed,
                                int nboot,
                                int niter,
                                int flag_norm,
                                int flag_propagation,
                                double percent,
                                double gibbs_eps,
                                int flag_ce,
                                int flag_cstd,
                                int verbose);
  GEOSLIB_API int gibbs_iter_monovariate(Props *propdef,
                                         Db *db,
                                         Model *model,
                                         double *covmat,
                                         int gibbs_nburn,
                                         int gibbs_niter,
                                         int isimu,
                                         int ipgs,
                                         int igrf,
                                         int nbsimu,
                                         double gibbs_eps,
                                         int verbose,
                                         double *mean);
  GEOSLIB_API int gibbs_iter_multivar(Props *propdef,
                                      Db *db,
                                      Model *model,
                                      int nboot,
                                      int maxiter,
                                      int isimu,
                                      int ipgs,
                                      double gibbs_eps,
                                      int verbose,
                                      double *mean);
  GEOSLIB_API int gibbs_iter_propagation(Props *propdef,
                                         Db *db,
                                         Model *model,
                                         int nboot,
                                         int maxiter,
                                         int isimu,
                                         int iptr,
                                         double gibbs_eps,
                                         int verbose,
                                         double *mean);
  GEOSLIB_API int simtub_constraints(Db *dbin,
                                     Db *dbout,
                                     Model *model,
                                     Neigh *neigh,
                                     int seed,
                                     int nbtuba,
                                     int nbsimu,
                                     int nbquant,
                                     int niter_max,
                                     VectorInt& cols,
                                     int (*func_valid)(int flag_grid,
                                                       int ndim,
                                                       int nech,
                                                       int *nx,
                                                       double *dx,
                                                       double *x0,
                                                       double nonval,
                                                       double percent,
                                                       VectorDouble& tab));
  GEOSLIB_API int db_simulations_to_ce(Db *db,
                                       ENUM_LOCS locatorType,
                                       int nbsimu,
                                       int nvar,
                                       int *iptr_ce_arg,
                                       int *iptr_cstd_arg);

  /*****************************************/
  /* Prototyping the functions in simfft.c */
  /*****************************************/

  GEOSLIB_API int simfft_f(Db *db,
                           Model *model,
                           int seed,
                           int nbsimu,
                           double percent,
                           int flag_aliasing);
  GEOSLIB_API int simfft_support(Db *db,
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

  GEOSLIB_API int simfine_dim(Db *dbin,
                              int nmult,
                              int *ndim,
                              int *ntot,
                              int *nx,
                              double *x0,
                              double *dx);
  GEOSLIB_API int simfine_f(Db *dbin,
                            Model *model,
                            int flag_ks,
                            int mult,
                            int seed,
                            VectorDouble& tab);

  /*****************************************/
  /* Prototyping the functions in simsub.c */
  /*****************************************/

  GEOSLIB_API int substitution(Db *dbgrid,
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

  GEOSLIB_API SubPlanes *poisson_manage_planes(int mode,
                                               int np,
                                               SubPlanes *splanes);
  GEOSLIB_API int poisson_generate_planes(Db *dbgrid, SubPlanes *splanes);
  GEOSLIB_API int tessellation_poisson(Db *dbgrid,
                                       Model *model,
                                       int seed,
                                       double intensity,
                                       int nbtuba,
                                       int verbose);
  GEOSLIB_API int tessellation_voronoi(Db *dbgrid,
                                       Model *model,
                                       double *dilate,
                                       int seed,
                                       double intensity,
                                       int nbtuba,
                                       int verbose);

  /*****************************************/
  /* Prototyping the functions in simsph.c */
  /*****************************************/
  GEOSLIB_API int simsph_f(Db *db,
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

  GEOSLIB_API Rule *rule_free(Rule *rule);
  GEOSLIB_API Rule *rule_init(int mode_rule,
                              double rho,
                              double slope,
                              double dinf,
                              double dsup,
                              const VectorDouble& shift,
                              const VectorInt& nodes,
                              int *nfacies,
                              int *ngrf,
                              int *ny1,
                              int *ny2);
  GEOSLIB_API Model *model_rule_combine(Model * model1,
                                        Model * model2,
                                        Rule * rule);
  GEOSLIB_API int rule_gaus2fac_result(Props *propdef,
                                       Db *dbout,
                                       Rule *rule,
                                       int *flag_used,
                                       int ipgs,
                                       int isimu,
                                       int nbsimu);
  GEOSLIB_API int rule_gaus2fac_result_shadow(Props *propdef,
                                              Db *dbout,
                                              Rule *rule,
                                              int *flag_used,
                                              int ipgs,
                                              int isimu,
                                              int nbsimu);
  GEOSLIB_API int rule_gaus2fac_data_shadow(Props *propdef,
                                            Db *dbin,
                                            Db *dbout,
                                            Rule *rule,
                                            int *flag_used,
                                            int ipgs,
                                            int isimu,
                                            int nbsimu);
  GEOSLIB_API int rule_gaus2fac_data(Props *propdef,
                                     Db *dbin,
                                     Db *dbout,
                                     Rule *rule,
                                     int *flag_used,
                                     int ipgs,
                                     int isimu,
                                     int nbsimu);
  GEOSLIB_API int rule_thresh_define_shadow(Props *propdef,
                                            Db *dbin,
                                            Rule *rule,
                                            int facies,
                                            int iech,
                                            int isimu,
                                            int nbsimu,
                                            int flag_check,
                                            double *t1min,
                                            double *t1max,
                                            double *t2min,
                                            double *t2max,
                                            double *dsup,
                                            double *down);
  GEOSLIB_API int rule_thresh_define(Props *propdef,
                                     Db *dbin,
                                     Rule *rule,
                                     int facies,
                                     int iech,
                                     int isimu,
                                     int nbsimu,
                                     int flag_check,
                                     double *t1min,
                                     double *t1max,
                                     double *t2min,
                                     double *t2max);
  GEOSLIB_API int rule_evaluate_bounds_shadow(Props *propdef,
                                              Db *dbin,
                                              Db *dbout,
                                              Rule *rule,
                                              int isimu,
                                              int igrf,
                                              int ipgs,
                                              int nbsimu,
                                              double delta);
  GEOSLIB_API int rule_evaluate_bounds(Props *propdef,
                                       Db *dbin,
                                       Db *dbout,
                                       Rule *rule,
                                       int isimu,
                                       int igrf,
                                       int ipgs,
                                       int nbsimu);
  GEOSLIB_API int db_rule_shadow(Db *db,
                                 Db *dbprop,
                                 Rule *rule,
                                 Model *model1,
                                 const VectorDouble& propcst,
                                 int flag_stat,
                                 int nfacies);
  GEOSLIB_API int db_bounds_shadow(Db *db,
                                   Db *dbprop,
                                   Rule *rule,
                                   Model *model,
                                   const VectorDouble& propcst,
                                   int flag_stat,
                                   int nfacies);
  GEOSLIB_API int db_bounds(Db *db,
                            Db *dbprop,
                            Rule *rule,
                            Model *model,
                            const VectorDouble& propcst,
                            int flag_stat,
                            int nfacies);
  GEOSLIB_API void proportion_rule_process(Props *propdef, int mode);
  GEOSLIB_API Props *proportion_manage(int mode,
                                       int flag_facies,
                                       int flag_stat,
                                       int ngrf1,
                                       int ngrf2,
                                       int nfac1,
                                       int nfac2,
                                       Db *db,
                                       Db *dbprop,
                                       const VectorDouble& propcst,
                                       Props *proploc);
  GEOSLIB_API void proportion_print(Props *propdef);

  /******************************************/
  /* Prototyping the functions in seismic.c */
  /******************************************/

  GEOSLIB_API int seismic_estimate_XZ(Db *db,
                                      Model *model,
                                      int nbench,
                                      int nv2max,
                                      int flag_ks,
                                      int flag_std,
                                      int flag_sort,
                                      int flag_stat);
  GEOSLIB_API int seismic_simulate_XZ(Db *db,
                                      Model *model,
                                      int nbench,
                                      int nv2max,
                                      int nbsimu,
                                      int seed,
                                      int flag_ks,
                                      int flag_sort,
                                      int flag_stat);
  GEOSLIB_API int seismic_z2t_grid(int verbose,
                                   Db *db_z,
                                   int iptr_v,
                                   int *nx,
                                   double *x0,
                                   double *dx);
  GEOSLIB_API int seismic_t2z_grid(int verbose,
                                   Db *db_t,
                                   int iptr_v,
                                   int *nx,
                                   double *x0,
                                   double *dx);
  GEOSLIB_API int seismic_z2t_convert(Db *db_z, int iptr_v, Db *db_t);
  GEOSLIB_API int seismic_t2z_convert(Db *db_t, int iptr_v, Db *db_z);
  GEOSLIB_API int seismic_operate(Db *db, int oper);
  GEOSLIB_API int seismic_convolve(Db *db,
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

  GEOSLIB_API Tokens *tokens_free(Tokens *tokens);
  GEOSLIB_API Tokens *tokens_create(int nb_tokens);
  GEOSLIB_API int tokone_create(Tokens *token,
                                int rank,
                                int type,
                                int npar,
                                double prop,
                                double factor_x2y,
                                double factor_x2z,
                                double factor_y2z,
                                int *law,
                                double *valarg);
  GEOSLIB_API void tokone_print(Tokens *tokens, int rank);
  GEOSLIB_API Tokens *tokens_input(void);
  GEOSLIB_API void tokone_get_nbparams(Tokens *tokens,
                                       int rank,
                                       int *types,
                                       int *npar,
                                       double *prop);
  GEOSLIB_API void tokone_get_params(Tokens *tokens,
                                     int rank,
                                     double *factor_x2y,
                                     double *factor_x2z,
                                     double *factor_y2z,
                                     int *law,
                                     double *valarg);
  GEOSLIB_API void tokens_print(Tokens *tokens);
  GEOSLIB_API int toktype_get_nbparams(int type);
  GEOSLIB_API int simbool_f(Db *dbin,
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

  GEOSLIB_API int time_3db(double *HS,
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

  GEOSLIB_API Polygons *polygon_create(void);
  GEOSLIB_API Polygons *polygon_free(Polygons *polygon);
  GEOSLIB_API Polygons *polygon_add(Polygons *polygon,
                                    const VectorDouble& x,
                                    const VectorDouble& y,
                                    double zmin,
                                    double zmax);
  GEOSLIB_API void polygon_print(Polygons *polygon, int flag_print);
  GEOSLIB_API int polygon_inside(double xx,
                                 double yy,
                                 double zz,
                                 int flag_nested,
                                 Polygons *polygon);
  GEOSLIB_API void polygon_extension(Polygons *polygon,
                                     double *xmin,
                                     double *xmax,
                                     double *ymin,
                                     double *ymax);
  GEOSLIB_API double polygon_surface(Polygons *polygon);
  GEOSLIB_API Polygons *input_polygon(void);
  GEOSLIB_API Polygons *polygon_hull(const Db *db);
  GEOSLIB_API int polygon_hull(const Db *db, VectorDouble& x, VectorDouble& y);

  /*******************************************/
  /* Prototyping the functions in variopgs.c */
  /*******************************************/

  GEOSLIB_API Vario_Order *vario_order_manage(int mode,
                                              int flag_dist,
                                              int size_aux,
                                              Vario_Order *vorder);
  GEOSLIB_API int vario_order_add(Vario_Order *vorder,
                                  int iech,
                                  int jech,
                                  void *aux_iech,
                                  void *aux_jech,
                                  int ipas,
                                  int idir,
                                  double dist);
  GEOSLIB_API Vario_Order *vario_order_final(Vario_Order *vorder, int *npair);
  GEOSLIB_API void vario_order_print(Vario_Order *vorder,
                                     int idir_target,
                                     int ilag_target,
                                     int verbose);
  GEOSLIB_API void vario_order_get_bounds(Vario_Order *vorder,
                                          int idir,
                                          int ipas,
                                          int *ifirst,
                                          int *ilast);
  GEOSLIB_API void vario_order_get_indices(Vario_Order *vorder,
                                           int ipair,
                                           int *iech,
                                           int *jech,
                                           double *dist);
  GEOSLIB_API void vario_order_get_auxiliary(Vario_Order *vorder,
                                             int ipair,
                                             char *aux_iech,
                                             char *aux_jech);
  GEOSLIB_API Rule *rule_auto(Db *db,
                              Db *dbprop,
                              Vario *vario,
                              Vario *varioind,
                              const VectorDouble& propcst,
                              int ncolor,
                              int ngrf,
                              int flag_stat,
                              int verbose);

  /**************************/
  /* Prototyping fracture.c */
  /**************************/
  GEOSLIB_API Frac_Environ *fracture_alloc_environ(int nfamilies,
                                                   double xmax,
                                                   double ymax,
                                                   double deltax,
                                                   double deltay,
                                                   double mean,
                                                   double stdev);
  GEOSLIB_API Frac_Environ *fracture_dealloc_environ(Frac_Environ *frac_environ);

  GEOSLIB_API Frac_List *fracture_manage_list(int mode, Frac_List *frac_list);
  GEOSLIB_API void fracture_update_family(Frac_Environ *frac_environ,
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
  GEOSLIB_API int fracture_add_fault(Frac_Environ *frac_environ,
                                     double fault_coord,
                                     double fault_orient);
  GEOSLIB_API void fracture_update_fault(Frac_Environ *frac_environ,
                                         int ifault,
                                         int family,
                                         double thetal,
                                         double thetar,
                                         double rangel,
                                         double ranger);
  GEOSLIB_API int fracture_simulate(Frac_Environ *frac_environ,
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
  GEOSLIB_API Frac_Environ *fracture_input(Frac_Environ *frac_def);
  GEOSLIB_API void fracture_print(Frac_Environ *frac_environ);
  GEOSLIB_API void fracture_list_print(const char *title,
                                       Frac_List *frac_list,
                                       int level);
  GEOSLIB_API void fracture_export(Frac_Environ *frac_environ,
                                   Frac_List *frac_list,
                                   int *nfracs_arg,
                                   int *nbyfrac_arg,
                                   double **frac_segs_arg);
  GEOSLIB_API Frac_List *fracture_import(int nval, double *frac_segs);
  GEOSLIB_API double *fracture_extract_length(Frac_List *frac_list,
                                              int family,
                                              double cote,
                                              double dcote,
                                              int *ntab);
  GEOSLIB_API double *fracture_extract_dist(Frac_List *frac_list,
                                            int family,
                                            double cote,
                                            double dcote,
                                            int *ntab);
  GEOSLIB_API int fracture_to_block(Db *dbgrid,
                                    Frac_List *frac_list,
                                    double *locinfo,
                                    int n_layers,
                                    int nfamilies,
                                    double xmax,
                                    double *permtab,
                                    int n_perm,
                                    double perm_mat,
                                    double perm_bench,
                                    int ndisc);
  GEOSLIB_API double *fracture_to_well(int nw_xy,
                                       double *well,
                                       Frac_List *frac_list,
                                       double xmax,
                                       double *permtab,
                                       int *nint,
                                       int *ncol);
  GEOSLIB_API int fracture_well_to_block(Db *dbgrid,
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
  GEOSLIB_API int variogram_mlayers(Db *db,
                                    int *seltab,
                                    Vario *vario,
                                    Vario_Order *vorder);
  GEOSLIB_API int multilayers_vario(Db *dbin,
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
  GEOSLIB_API int multilayers_kriging(Db *dbin,
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
  GEOSLIB_API int multilayers_get_prior(Db *dbin,
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
  GEOSLIB_API int potential_cov(Model *model,
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
  GEOSLIB_API int potential_kriging(Db *db,
                                    Db *dbgrd,
                                    Db *dbtgt,
                                    Db *dbout,
                                    Model *model,
                                    Neigh *neigh,
                                    double nugget_grd,
                                    double nugget_tgt,
                                    int z_number,
                                    double z_min,
                                    double _max,
                                    int flag_up,
                                    int flag_grad,
                                    int flag_trans,
                                    int flag_drift,
                                    int verbose);
  GEOSLIB_API int potential_simulate(Db *dbiso,
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
  GEOSLIB_API int potential_xvalid(Db *dbiso,
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
  GEOSLIB_API int MSS(int idim, int ipol, int icas, int icorn, int icoor);
  GEOSLIB_API Vercoloc *vercoloc_manage(int verbose,
                                        int mode,
                                        Db *dbin,
                                        Db *dbgrid,
                                        int mesh_dbin,
                                        Vercoloc *vercoloc);
  GEOSLIB_API Vercoloc *vercoloc_from_external(int ndupl,
                                               int *dupl_in,
                                               int *dupl_out);
  GEOSLIB_API Vertype *vertype_manage(int mode,
                                      Vertype *vertype,
                                      Vercoloc *vercoloc,
                                      int nvertex);
  GEOSLIB_API int *vercoloc_get_dbin_indices(Vertype *vertype,
                                             Vercoloc *vercoloc,
                                             int *nbnodup);
  GEOSLIB_API void triangulate(const char* triswitches,
                               struct triangulateio *in,
                               struct triangulateio *out,
                               struct triangulateio *vorout);

  GEOSLIB_API void meshes_1D_free(segmentio *t, int mode);
  GEOSLIB_API void meshes_1D_init(segmentio *t);
  GEOSLIB_API int meshes_1D_from_db(Db *db, int nmask, int *mask, segmentio *t);
  GEOSLIB_API int meshes_1D_from_points(int nech, double *x, segmentio *t);
  GEOSLIB_API void meshes_1D_default(Db *dbin, Db *dbout, segmentio *t);
  GEOSLIB_API void meshes_1D_print(segmentio *t, int brief);
  GEOSLIB_API int meshes_turbo_1D_grid_build(int verbose,
                                             Db *dbgrid,
                                             SPDE_Mesh *s_mesh);
  GEOSLIB_API void meshes_1D_create(int verbose,
                                    struct segmentio *in,
                                    struct segmentio *out);
  GEOSLIB_API void meshes_1D_load_vertices(segmentio *t,
                                           const char *name,
                                           int *ntab_arg,
                                           int *natt_arg,
                                           void **tab_arg);
  GEOSLIB_API void meshes_1D_extended_domain(Db *dbout,
                                             const double *gext,
                                             segmentio *t);

  GEOSLIB_API void meshes_2D_free(triangulateio *t, int mode);
  GEOSLIB_API void meshes_2D_init(triangulateio *t);
  GEOSLIB_API int meshes_2D_from_db(Db *db,
                                    int use_code,
                                    int nmask,
                                    int *mask,
                                    triangulateio *t);
  GEOSLIB_API int meshes_2D_from_points(int nech,
                                        double *x,
                                        double *y,
                                        triangulateio *t);
  GEOSLIB_API void meshes_2D_default(Db *dbin, Db *dbout, triangulateio *t);
  GEOSLIB_API int meshes_2D_from_mem(int nseg,
                                     int ncol,
                                     int *segments,
                                     triangulateio *t);
  GEOSLIB_API void meshes_2D_print(triangulateio *t, int brief);
  GEOSLIB_API int meshes_2D_write(const char *file_name,
                                  const char *obj_name,
                                  int verbose,
                                  int ndim,
                                  int ncode,
                                  int ntri,
                                  int npoints,
                                  int *ntcode,
                                  int *triangles,
                                  double *points);
  GEOSLIB_API int meshes_turbo_2D_grid_build(int verbose,
                                             Db *dbgrid,
                                             SPDE_Mesh *s_mesh);
  GEOSLIB_API void meshes_2D_create(int verbose,
                                    const String& triswitches,
                                    struct triangulateio *in,
                                    struct triangulateio *out,
                                    struct triangulateio *vorout);
  GEOSLIB_API void meshes_2D_load_vertices(triangulateio *t,
                                           const char *name,
                                           int *ntab_arg,
                                           int *natt_arg,
                                           void **tab_arg);
  GEOSLIB_API void meshes_2D_extended_domain(Db *dbout,
                                             const double *gext,
                                             triangulateio *t);

  GEOSLIB_API void meshes_3D_create(int verbose,
                                    const String& triswitch,
                                    tetgenio *in,
                                    tetgenio *out);
  GEOSLIB_API int meshes_3D_from_db(Db *db, int nmask, int *mask, tetgenio *t);
  GEOSLIB_API int meshes_3D_from_points(int nech,
                                        double *x,
                                        double *y,
                                        double *z,
                                        tetgenio *t);
  GEOSLIB_API void meshes_3D_default(Db *dbin, Db *dbout, tetgenio *t);
  GEOSLIB_API void meshes_3D_free(tetgenio *t);
  GEOSLIB_API void meshes_3D_extended_domain(Db *dbout,
                                             const double *gext,
                                             tetgenio *t);
  GEOSLIB_API int meshes_turbo_3D_grid_build(int verbose,
                                             Db *dbgrid,
                                             SPDE_Mesh *s_mesh);
  GEOSLIB_API void meshes_3D_load_vertices(tetgenio *t,
                                           const char *name,
                                           int *ntab_arg,
                                           int *natt_arg,
                                           void **tab_arg);
  GEOSLIB_API void meshes_3D_print(tetgenio *t, int brief);
  GEOSLIB_API void mesh_stats(int ndim,
                              int ncorner,
                              int nmesh,
                              int *meshes,
                              double *points);
  GEOSLIB_API void meshes_2D_sph_init(SphTriangle *t);
  GEOSLIB_API void meshes_2D_sph_free(SphTriangle *t, int mode);
  GEOSLIB_API int meshes_2D_sph_from_db(Db *db,
                                        int nmask,
                                        int *mask,
                                        SphTriangle *t);
  GEOSLIB_API int meshes_2D_sph_from_points(int nech,
                                            double *x,
                                            double *y,
                                            SphTriangle *t);
  GEOSLIB_API int meshes_2D_sph_from_auxiliary(const String& triswitch,
                                               SphTriangle *t);
  GEOSLIB_API void meshes_2D_sph_print(SphTriangle *t, int brief);
  GEOSLIB_API int meshes_2D_sph_create(int verbose, SphTriangle *t);
  GEOSLIB_API void meshes_2D_sph_load_vertices(SphTriangle *t,
                                               const char *name,
                                               int *ntab_arg,
                                               int *natt_arg,
                                               void **tab_arg);
  GEOSLIB_API int trmesh_(int *n,
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
  GEOSLIB_API int trlist_(int *n,
                          int *list,
                          int *lptr,
                          int *lend,
                          int *nrow,
                          int *nt,
                          int *ltri,
                          int *ier);
  GEOSLIB_API void util_convert_sph2cart(double rlong,
                                         double rlat,
                                         double *x,
                                         double *y,
                                         double *z);
  GEOSLIB_API void util_convert_cart2sph(double x,
                                         double y,
                                         double z,
                                         double *rlong,
                                         double *rlat);

  /***************************************/
  /* Prototyping the functions in spde.c */
  /***************************************/
  GEOSLIB_API QChol *qchol_manage(int mode, QChol *qchol);
  GEOSLIB_API double spde_compute_correc(int ndim, double param);
  GEOSLIB_API int spde_check(Db *dbin,
                             Db *dbout,
                             Model *model1,
                             Model *model2,
                             int verbose,
                             const VectorDouble& gext,
                             int mesh_dbin,
                             int mesh_dbout,
                             int flag_advanced,
                             int flag_est,
                             int flag_std,
                             int flag_gibbs,
                             int flag_modif);
  GEOSLIB_API int spde_attach_model(Model *model);
  GEOSLIB_API int m2d_gibbs_spde(Db *dbin,
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
  GEOSLIB_API SPDE_Option spde_option_alloc(void);
  GEOSLIB_API void spde_option_update(SPDE_Option& s_option, const String& triswitch);
  GEOSLIB_API int spde_prepar(Db *dbin,
                              Db *dbout,
                              const VectorDouble& gext,
                              SPDE_Option& s_option);
  GEOSLIB_API int spde_process(Db *dbin,
                               Db *dbout,
                               SPDE_Option& s_option,
                               int nbsimu,
                               int gibbs_nburn,
                               int gibbs_niter,
                               int ngibbs_int);
  GEOSLIB_API SPDE_Mesh *spde_mesh_manage(int mode, SPDE_Mesh *s_mesh_old);
  GEOSLIB_API SPDE_Matelem& spde_get_current_matelem(int icov);
  GEOSLIB_API int spde_mesh_load(SPDE_Mesh *s_mesh,
                                 int verbose,
                                 Db *dbin,
                                 Db *dbout,
                                 const VectorDouble& gext,
                                 SPDE_Option& s_option);
  GEOSLIB_API void spde_mesh_assign(SPDE_Mesh *s_mesh,
                                    int ndim,
                                    int ncorner,
                                    int nvertex,
                                    int nmesh,
                                    int *meshes,
                                    double *points,
                                    int verbose);
  GEOSLIB_API int spde_build_matrices(Model* model, int verbose);
  GEOSLIB_API int spde_build_stdev(double *vcur);
  GEOSLIB_API int spde_f(Db *dbin,
                         Db *dbout,
                         Model *model,
                         const VectorDouble& gext,
                         SPDE_Option& s_option,
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
  GEOSLIB_API int spde_eval(int nblin,
                            double *blin,
                            cs *S,
                            const VectorDouble& Lambda,
                            const VectorDouble& TildeC,
                            double power,
                            double *x,
                            double *y);
  GEOSLIB_API int spde_external_mesh_define(int mode,
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
  GEOSLIB_API int spde_external_AQ_define(int mode,
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
  GEOSLIB_API int spde_external_mesh_copy(SPDE_Mesh *s_mesh, int icov0);
  GEOSLIB_API int kriging2D_spde(Db *dbin,
                                 Model *model,
                                 SPDE_Option& s_option,
                                 int verbose,
                                 int *ntri,
                                 int *npoint,
                                 int **triangles,
                                 double **points);
  GEOSLIB_API cs *db_mesh_sparse(Db *db, MeshEStandard *amesh, int verbose);
  GEOSLIB_API cs *db_mesh_neigh(Db *db,
                                SPDE_Mesh *s_mesh,
                                double radius,
                                int flag_exact,
                                int verbose,
                                int *nactive,
                                int **ranks);
  GEOSLIB_API void spde_free_all(void);

  /***************************************/
  /* Prototyping the functions in math.c */
  /***************************************/

  GEOSLIB_API int db_trisurf(Db *db,
                             Model *model,
                             const String& triswitch,
                             int icode0,
                             int verbose,
                             int *ncode_arg,
                             int *ntri_arg,
                             int *npoint_arg,
                             double *codesel,
                             int **ntcode_arg,
                             int **triangle_arg,
                             double **points_arg);
  GEOSLIB_API CTables *ct_tables_manage(int mode,
                                        int verbose,
                                        int flag_cumul,
                                        int ndim,
                                        int nconf,
                                        int ndisc,
                                        double cmin,
                                        double cmax,
                                        CTables *ctables_old);
  GEOSLIB_API void ct_tables_print(CTables *ctables, int flag_print);
  GEOSLIB_API int ct_tableone_covrank(CTables *ctables,
                                      double cova,
                                      double *cround);
  GEOSLIB_API int ct_tableone_getrank_from_proba(CTables *ctables,
                                                 double gaussian);
  GEOSLIB_API double ct_tableone_calculate(CTables *ctables,
                                           int iconf0,
                                           double *lows,
                                           double *ups);
  GEOSLIB_API double ct_tableone_calculate_by_rank(CTables *ctables,
                                                   int iconf0,
                                                   double *rklows,
                                                   double *rkups);
  GEOSLIB_API double ct_INTRES2(CTables *ctables,
                                int iconf0,
                                int idisc0,
                                int jdisc0);
  GEOSLIB_API double ct_INTRES3(CTables *ctables,
                                int iconf0,
                                int idisc0,
                                int jdisc0,
                                int kdisc0);

  /******************************************/
  /* Prototyping the functions in cluster.c */
  /******************************************/
  GEOSLIB_API double *kclusters(double *data,
                                int nvar,
                                int nech,
                                int nclusters,
                                int npass,
                                int mode,
                                int verbose);
  GEOSLIB_API int *kmedoids(double *data,
                            int nvar,
                            int nech,
                            int nclusters,
                            int npass,
                            int verbose);

  /******************************************/
  /* Prototyping the functions in pthread.c */
  /******************************************/
  GEOSLIB_API Th_Rank th_delete();
  GEOSLIB_API Th_Rank th_create(void *(*func)(void *), void *arg);
  GEOSLIB_API int th_wait_for_all(int number, Th_Rank *th_ranks);
  GEOSLIB_API void th_mutex_init(Th_Mutex *th_mutex);
  GEOSLIB_API void th_mutex_lock(Th_Mutex *th_mutex);
  GEOSLIB_API void th_mutex_unlock(Th_Mutex *th_mutex);
  GEOSLIB_API void th_cond_init(Th_Cond *th_cond);
  GEOSLIB_API void th_cond_signal(Th_Cond *th_cond);
  GEOSLIB_API void th_cond_wait(Th_Cond *th_cond, Th_Mutex *th_mutex);
  GEOSLIB_API void th_cond_release_all(Th_Cond *th_cond);
  GEOSLIB_API int  ThC_get_ncores(void);
  GEOSLIB_API void ThC_create(ctpl::thread_pool &tp,
                              std::vector<std::future<void>> &futures,
                              int nTasks,
                              void *(*func)(int, int, void *),
                              void *arg);
  GEOSLIB_API void ThC_wait_for_all(std::vector<std::future<void>> &futures);

//#ifdef __cplusplus
//}
//#endif

#endif
