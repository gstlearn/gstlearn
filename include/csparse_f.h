#ifndef _CS_F_H
#define _CS_F_H

#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include "csparse_d.h"

GSTLEARN_EXPORT cs     *cs_add (const cs *A, const cs *B, double alpha, double beta) ;
GSTLEARN_EXPORT void    cs_add_cste(cs *A, double value);
GSTLEARN_EXPORT int     cs_cholsol (const cs *A, double *b, int order) ;
GSTLEARN_EXPORT int     cs_dupl (cs *A) ;
GSTLEARN_EXPORT int     cs_entry (cs *T, int i, int j, double x) ;
GSTLEARN_EXPORT int     cs_lusol (const cs *A, double *b, int order, double tol) ;
GSTLEARN_EXPORT int     cs_gaxpy (const cs *A, const double *x, double *y) ;
GSTLEARN_EXPORT cs     *cs_multiply (const cs *A, const cs *B) ;
GSTLEARN_EXPORT int     cs_qrsol (const cs *A, double *b, int order) ;
GSTLEARN_EXPORT cs     *cs_transpose (const cs *A, int values) ;
GSTLEARN_EXPORT void    cs_tmulvec(const cs *A, int nout, const double *x, double *y);
GSTLEARN_EXPORT void    cs_mulvec(const cs *A, int nout, const double *x, double *y);
GSTLEARN_EXPORT void    cs_vecmult(const cs *A, int nout, const double *x, double *y);
GSTLEARN_EXPORT void    cs_mulvec_uptri(const cs *A, int nout,
                                        const double *x, double *y, int flag_diag);
GSTLEARN_EXPORT void    cs_mulvec_lowtri(const cs *A, int nout,
                                         const double *x, double *y, int flag_diag);
GSTLEARN_EXPORT cs     *cs_matvecR(const cs *A, double *x, int oper);
GSTLEARN_EXPORT cs     *cs_matvecL(const cs *A, double *x, int oper);
GSTLEARN_EXPORT cs     *cs_matvecnorm(const cs *A, const double *x, int oper);
GSTLEARN_EXPORT void    cs_matvecnorm_inplace(cs *A, const double *x, int oper);
GSTLEARN_EXPORT cs     *cs_triplet (const cs *T) ;
GSTLEARN_EXPORT cs*     cs_diag(VectorDouble diag);
GSTLEARN_EXPORT double  cs_norm (const cs *A) ;
GSTLEARN_EXPORT int     cs_print (const cs *A, int brief) ;
GSTLEARN_EXPORT void    cs_print_only(const char *title, const cs *A,int nlimit);
GSTLEARN_EXPORT void    cs_print_nice (const char *title,const cs *A, int maxrow=-1, int maxcol=-1);
GSTLEARN_EXPORT cs     *cs_load (FILE *f) ;
GSTLEARN_EXPORT double *cs_col_sumrow(const cs *A,int *ncol,int *nrow);
GSTLEARN_EXPORT double  cs_maxsumabscol(const cs *A);
GSTLEARN_EXPORT void    cs_print_dim(const char *title,const cs *A);
GSTLEARN_EXPORT void    cs_print_short(const char *title, const cs *L, int nmax);
GSTLEARN_EXPORT void    cs_print_file(const char *radix, int rank, cs *A);
GSTLEARN_EXPORT cs     *cs_compress(cs *A);
GSTLEARN_EXPORT int    *cs_color_coding(cs *Q,int start,int *ncolor);
GSTLEARN_EXPORT cs     *cs_invert(const cs *A, int order, double epsilon = EPSILON6);

/* utilities */
GSTLEARN_EXPORT void   *cs_calloc  (int n, size_t size) ;
GSTLEARN_EXPORT void   *cs_free    (void *p) ;
GSTLEARN_EXPORT void   *cs_realloc (void *p, int n, size_t size, int *ok) ;
GSTLEARN_EXPORT cs     *cs_spalloc (int m, int n, int nzmax, int values, int triplet) ;
GSTLEARN_EXPORT cs     *cs_spfree  (cs *A) ;
GSTLEARN_EXPORT int     cs_sprealloc (cs *A, int nzmax) ;
GSTLEARN_EXPORT void   *cs_malloc  (int n, size_t size) ;

/* --- secondary CSparse routines and data structures ----------------------- */

GSTLEARN_EXPORT int    *cs_amd (const cs *A, int order) ;
GSTLEARN_EXPORT csn    *cs_chol (const cs *A, const css *S) ;
GSTLEARN_EXPORT csd    *cs_dmperm (const cs *A) ;
GSTLEARN_EXPORT int     cs_droptol (cs *A, double tol) ;
GSTLEARN_EXPORT int     cs_dropzeros (cs *A) ;
GSTLEARN_EXPORT int     cs_happly (const cs *V, int i, double beta, double *x) ;
GSTLEARN_EXPORT int     cs_ipvec (int n, const int *P, const double *b, double *x) ;
GSTLEARN_EXPORT int     cs_lsolve_lowtri( const cs *L, const double *x, double *y);
GSTLEARN_EXPORT int     cs_lsolve_uptri (const cs *L, const double *x, double *y);
GSTLEARN_EXPORT int     cs_lsolve (const cs *L, double *x) ;
GSTLEARN_EXPORT int     cs_ltsolve (const cs *L, double *x) ;
GSTLEARN_EXPORT csn    *cs_lu (const cs *A, const css *S, double tol) ;
GSTLEARN_EXPORT cs     *cs_permute (const cs *A, const int *P, const int *Q, int values) ;
GSTLEARN_EXPORT int    *cs_pinv (const int *P, int n) ;
GSTLEARN_EXPORT int     cs_pvec (int n, const int *P, const double *b, double *x) ;
GSTLEARN_EXPORT csn    *cs_qr (const cs *A, const css *S) ;
GSTLEARN_EXPORT css    *cs_schol (const cs *A, int order) ;
GSTLEARN_EXPORT css    *cs_sqr (const cs *A, int order, int qr) ;
GSTLEARN_EXPORT cs     *cs_symperm (const cs *A, const int *Pinv, int values) ;
GSTLEARN_EXPORT int     cs_usolve (const cs *U, double *x) ;
GSTLEARN_EXPORT int     cs_utsolve (const cs *U, double *x) ;
GSTLEARN_EXPORT int     cs_updown (cs *L, int sigma, const cs *C, const int *parent) ;
GSTLEARN_EXPORT css    *cs_sfree (css *S) ;
GSTLEARN_EXPORT csn    *cs_nfree (csn *N) ;
GSTLEARN_EXPORT csd    *cs_dfree (csd *D) ;

/* --- tertiary CSparse routines -------------------------------------------- */
GSTLEARN_EXPORT int    *cs_counts (const cs *A, const int *parent, const int *post, int ata) ;
GSTLEARN_EXPORT int     cs_cumsum (int *p, int *c, int n) ;
GSTLEARN_EXPORT int     cs_dfs (int j, cs *L, int top, int *xi, int *pstack, const int *Pinv) ;
GSTLEARN_EXPORT int    *cs_etree (const cs *A, int ata) ;
GSTLEARN_EXPORT int     cs_fkeep (cs *A, int (*fkeep) (int, int, double, void *), void *other) ;
GSTLEARN_EXPORT double  cs_house (double *x, double *beta, int n) ;
GSTLEARN_EXPORT int    *cs_maxtrans (const cs *A) ;
GSTLEARN_EXPORT int    *cs_post (int n, const int *parent) ;
GSTLEARN_EXPORT int     cs_reach (cs *L, const cs *B, int k, int *xi, const int *Pinv) ;
GSTLEARN_EXPORT csd    *cs_scc (cs *A) ;
GSTLEARN_EXPORT int     cs_scatter (const cs *A, int j, double beta, int *w, double *x,
                                    int mark, cs *C, int nz) ;
GSTLEARN_EXPORT int     cs_splsolve (cs *L, const cs *B, int k, int *xi, double *x,
                                     const int *Pinv) ;
GSTLEARN_EXPORT int     cs_tdfs (int j, int k, int *head, const int *next, int *post,
                                 int *stack) ;

/* utilities */
GSTLEARN_EXPORT csd    *cs_dalloc (int m, int n) ;
GSTLEARN_EXPORT cs     *cs_done (cs *C, void *w, void *x, int ok) ;
GSTLEARN_EXPORT int    *cs_idone (int *p, cs *C, void *w, int ok) ;
GSTLEARN_EXPORT csn    *cs_ndone (csn *N, cs *C, void *w, void *x, int ok) ;
GSTLEARN_EXPORT csd    *cs_ddone (csd *D, cs *C, void *w, int ok) ;
GSTLEARN_EXPORT void    cs_sparse_to_triplet (const cs *A, int flag_from_1,
                                              int *number,int **cols,int **rows,double ** vals);
GSTLEARN_EXPORT cs     *cs_arrays_to_sparse(int n, int nrow, int ncol,
                                            double *rows, double *cols, double *vals);
GSTLEARN_EXPORT cs     *cs_extract_submatrix(cs *C,
                                             int row_from, int row_length,
                                             int col_from, int col_length);
GSTLEARN_EXPORT cs     *cs_extract_submatrix_by_ranks(cs *C, int *row_array, int *col_array);
GSTLEARN_EXPORT cs     *cs_extract_submatrix_by_color(cs *C,int *colors,
                                                      int ref_color, int row_ok, int col_ok);
GSTLEARN_EXPORT void    cs_print_range(const char *title,const cs *C);
GSTLEARN_EXPORT cs     *cs_eye(int number,double value);
GSTLEARN_EXPORT cs     *cs_eye_tab(int number, double *values);
GSTLEARN_EXPORT cs     *cs_extract_diag(cs *C,int mode);
GSTLEARN_EXPORT void    cs_diag_suppress(cs *C);
GSTLEARN_EXPORT double *csd_extract_diag(cs *C,int mode);
GSTLEARN_EXPORT int     cs_sort_i(cs *C);
GSTLEARN_EXPORT int     sparseinv(int n, int *Lp, int *Li, double *Lx, double *d, int *Up,
                                  int *Uj, double *Ux, int *Zp, int *Zi, double *Zx, double *z,
                                  int *Zdiagp, int *Lmunch);
GSTLEARN_EXPORT void    cs_rowcol(const cs *A,int *nrows,int *ncols,int *count,double *percent);
GSTLEARN_EXPORT cs     *cs_duplicate(const cs *b1);
GSTLEARN_EXPORT cs     *cs_multiply_and_release(cs *b1,cs *b2,int flag_rel);
GSTLEARN_EXPORT cs     *cs_add_and_release(cs *b1, cs *b2, double alpha, double beta,
                                           int flag_rel);
GSTLEARN_EXPORT cs     *cs_normalize_by_diag_and_release(cs *Q, int flag_rel);
GSTLEARN_EXPORT cs     *cs_prod_norm(int mode, cs *A, cs *IhH);
GSTLEARN_EXPORT cs     *cs_prod_norm_single(int mode, cs *B);
GSTLEARN_EXPORT cs     *cs_prod_norm_diagonal(int mode, cs *B, VectorDouble diag);
GSTLEARN_EXPORT cs     *cs_prod_norm_and_release(cs *b1, cs *lambda, int flag_rel);
GSTLEARN_EXPORT int     cs_coarsening(cs *Q,int type,int **indCo,cs **L);
GSTLEARN_EXPORT cs     *cs_interpolate(cs *AA,cs *LL,int *indCo);
GSTLEARN_EXPORT cs     *cs_triangle(cs *A, int flag_upper, int flag_diag);
GSTLEARN_EXPORT void    cs_keypair(const char *key, cs *A, int flag_from_1);
GSTLEARN_EXPORT int     cs_scale(cs *C);
GSTLEARN_EXPORT int     cs_get_nrow(const cs *A);
GSTLEARN_EXPORT int     cs_get_ncol(const cs *A);
GSTLEARN_EXPORT int     cs_get_ncell(const cs *A);
GSTLEARN_EXPORT double  cs_get_value(const cs *A,int row, int col);
GSTLEARN_EXPORT void    cs_set_value(const cs *A,int row, int col, double value);
GSTLEARN_EXPORT double* cs_toArray(const cs *A);
GSTLEARN_EXPORT cs*     cs_strip(cs *A, double eps, int hypothesis = 3, bool verbose = false);
GSTLEARN_EXPORT int     cs_nnz(const cs* A);

GSTLEARN_EXPORT void cs_force_dimension(cs *T, int nrow, int ncol);

// Qchol operations
GSTLEARN_EXPORT int qchol_cholesky(int verbose,QChol *QC);
GSTLEARN_EXPORT void cs_chol_invert(QChol *qctt,double *xcr,double *rhs, double *work);
GSTLEARN_EXPORT void cs_chol_simulate(QChol	*qctt,double *simu,double *work);

// Multigrid operations
GSTLEARN_EXPORT cs_MGS *cs_multigrid_manage(cs_MGS *mgs,int mode,
                                            int nlevels, int path_type);
GSTLEARN_EXPORT void cs_multigrid_params(cs_MGS *mgs, int flag_cg,
                                         int type_coarse, int ngc, int nmg, int ngs,
                                         double tolgc, double tolnmg);
GSTLEARN_EXPORT void cs_multigrid_print(cs_MGS *mgs);
GSTLEARN_EXPORT int  cs_multigrid_get_nlevels(cs_MGS *mgs);
GSTLEARN_EXPORT int  cs_multigrid_setup(cs_MGS *mgs, QChol *Qctt,
                                        int flag_sel, int verbose, double **sel);
GSTLEARN_EXPORT int  cs_multigrid_process(cs_MGS *mgs, QChol *qctt, int verbose,
                                          double *x, double *b, double *work);
GSTLEARN_EXPORT void cs_multigrid_coarse2fine(cs_MGS *mgs,double *z,double *work);

GSTLEARN_EXPORT Triplet csToTriplet(const cs *A, bool flag_from_1 = false);
GSTLEARN_EXPORT String toStringDim(const String& title, const cs *A);
GSTLEARN_EXPORT String toStringRange(const String& title, const cs *C);
GSTLEARN_EXPORT bool cs_isSymmetric(const cs* A, bool verbose = false, bool detail = false);
GSTLEARN_EXPORT bool cs_isDiagonalDominant(cs *A, bool verbose = false, bool detail = false);
GSTLEARN_EXPORT bool cs_isDefinitePositive(cs* A, bool verbose = false);

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(Ap,j) (Ap [j] < 0)
#define CS_MARK(Ap,j) { Ap [j] = CS_FLIP (Ap [j]) ; }
#define CS_OVERFLOW(n,size) (n > INT_MAX / (int) size)
#endif
