#ifndef _CS_F_H
#define _CS_F_H

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif
#define CS_VER 1		    /* CSparse Version 1.2.0 */
#define CS_SUBVER 2
#define CS_SUBSUB 0
#define CS_DATE "Mar 6, 2006"	    /* CSparse release date */
#define CS_COPYRIGHT "Copyright (c) Timothy A. Davis, 2006"

#include "csparse_d.h"

cs     *cs_add (const cs *A, const cs *B, double alpha, double beta) ;
void    cs_add_cste(cs *A, double value);
int     cs_cholsol (const cs *A, double *b, int order) ;
int     cs_dupl (cs *A) ;
int     cs_entry (cs *T, int i, int j, double x) ;
int     cs_lusol (const cs *A, double *b, int order, double tol) ;
int     cs_gaxpy (const cs *A, const double *x, double *y) ;
cs     *cs_multiply (const cs *A, const cs *B) ;
int     cs_qrsol (const cs *A, double *b, int order) ;
cs     *cs_transpose (const cs *A, int values) ;
void    cs_tmulvec(const cs *A, const double *x, double *y);
void    cs_mulvec(const cs *A, int nout, const double *x, double *y);
void    cs_vecmult (const cs *A, const double *x, double *y);
void    cs_mulvec_uptri(const cs *A, int nout, 
                        const double *x, double *y, int flag_diag);
void    cs_mulvec_lowtri(const cs *A, int nout, 
                         const double *x, double *y, int flag_diag);
cs     *cs_matvecR(const cs *A, double *x, int oper);
cs     *cs_matvecL(const cs *A, double *x, int oper);
cs     *cs_matvecnorm(const cs *A, const double *x, int oper);
void    cs_matvecnorm_inplace(cs *A, const double *x, int oper);
cs     *cs_triplet (const cs *T) ;
double  cs_norm (const cs *A) ;
int     cs_print (const cs *A, int brief) ;
void    cs_print_only(const char *title, const cs *A,int nlimit);
void    cs_print_nice (const char *title,const cs *A, int maxrow, int maxcol);
cs     *cs_load (FILE *f) ;
double *cs_col_sumrow(const cs *A,int *ncol,int *nrow);
void    cs_print_dim(const char *title,const cs *A);
void    cs_print_short(const char *title, const cs *L, int nmax);
void    cs_print_file(const char *radix, int rank, cs *A);
cs     *cs_compress(cs *A);
int    *cs_color_coding(cs *Q,int start,int *ncolor);
cs     *cs_invert(const cs *A, int order, double epsilon = EPSILON6);

/* utilities */
void   *cs_calloc  (int n, size_t size) ;
void   *cs_free    (void *p) ;
void   *cs_realloc (void *p, int n, size_t size, int *ok) ;
cs     *cs_spalloc (int m, int n, int nzmax, int values, int triplet) ;
cs     *cs_spfree  (cs *A) ;
int     cs_sprealloc (cs *A, int nzmax) ;
void   *cs_malloc  (int n, size_t size) ;

/* --- secondary CSparse routines and data structures ----------------------- */

int    *cs_amd (const cs *A, int order) ;
csn    *cs_chol (const cs *A, const css *S) ;
csd    *cs_dmperm (const cs *A) ;
int     cs_droptol (cs *A, double tol) ;
int     cs_dropzeros (cs *A) ;
int     cs_happly (const cs *V, int i, double beta, double *x) ;
int     cs_ipvec (int n, const int *P, const double *b, double *x) ;
int     cs_lsolve_lowtri( const cs *L, const double *x, double *y);
int     cs_lsolve_uptri (const cs *L, const double *x, double *y);
int     cs_lsolve (const cs *L, double *x) ;
int     cs_ltsolve (const cs *L, double *x) ;
csn    *cs_lu (const cs *A, const css *S, double tol) ;
cs     *cs_permute (const cs *A, const int *P, const int *Q, int values) ;
int    *cs_pinv (const int *P, int n) ;
int     cs_pvec (int n, const int *P, const double *b, double *x) ;
csn    *cs_qr (const cs *A, const css *S) ;
css    *cs_schol (const cs *A, int order) ;
css    *cs_sqr (const cs *A, int order, int qr) ;
cs     *cs_symperm (const cs *A, const int *Pinv, int values) ;
int     cs_usolve (const cs *U, double *x) ;
int     cs_utsolve (const cs *U, double *x) ;
int     cs_updown (cs *L, int sigma, const cs *C, const int *parent) ;
css    *cs_sfree (css *S) ;
csn    *cs_nfree (csn *N) ;
csd    *cs_dfree (csd *D) ;

/* --- tertiary CSparse routines -------------------------------------------- */
int    *cs_counts (const cs *A, const int *parent, const int *post, int ata) ;
int     cs_cumsum (int *p, int *c, int n) ;
int     cs_dfs (int j, cs *L, int top, int *xi, int *pstack, const int *Pinv) ;
int    *cs_etree (const cs *A, int ata) ;
int     cs_fkeep (cs *A, int (*fkeep) (int, int, double, void *), void *other) ;
double  cs_house (double *x, double *beta, int n) ;
int    *cs_maxtrans (const cs *A) ;
int    *cs_post (int n, const int *parent) ;
int     cs_reach (cs *L, const cs *B, int k, int *xi, const int *Pinv) ;
csd    *cs_scc (cs *A) ;
int     cs_scatter (const cs *A, int j, double beta, int *w, double *x,
                    int mark, cs *C, int nz) ;
int     cs_splsolve (cs *L, const cs *B, int k, int *xi, double *x,
                     const int *Pinv) ;
int     cs_tdfs (int j, int k, int *head, const int *next, int *post,
                 int *stack) ;

/* utilities */
csd    *cs_dalloc (int m, int n) ;
cs     *cs_done (cs *C, void *w, void *x, int ok) ;
int    *cs_idone (int *p, cs *C, void *w, int ok) ;
csn    *cs_ndone (csn *N, cs *C, void *w, void *x, int ok) ;
csd    *cs_ddone (csd *D, cs *C, void *w, int ok) ;
void    cs_sparse_to_triplet (const cs *A, int flag_from_1,
                              int *number,int **cols,int **rows,double ** vals);
cs     *cs_arrays_to_sparse(int n, int nrow, int ncol,
                            double *rows, double *cols, double *vals);
cs     *cs_extract_submatrix(cs *C, 
                             int row_from, int row_length,
                             int col_from, int col_length);
cs     *cs_extract_submatrix_by_ranks(cs *C, int *row_array, int *col_array);
cs     *cs_extract_submatrix_by_color(cs *C,int *colors,
                                      int ref_color, int row_ok, int col_ok);
void    cs_print_range(const char *title,const cs *C);
cs     *cs_eye(int number,double value);
cs     *cs_eye_tab(int number, double *values);
cs     *cs_extract_diag(cs *C,int mode);
void    cs_diag_suppress(cs *C);
double *csd_extract_diag(cs *C,int mode);
int     cs_sort_i(cs *C);
int     sparseinv(int n, int *Lp, int *Li, double *Lx, double *d, int *Up, 
                  int *Uj, double *Ux, int *Zp, int *Zi, double *Zx, double *z, 
                  int *Zdiagp, int *Lmunch);
void    cs_rowcol(const cs *A,int *nrows,int *ncols,int *count,double *percent);
cs     *cs_duplicate(const cs *b1);
cs     *cs_multiply_and_release(cs *b1,cs *b2,int flag_rel);
cs     *cs_add_and_release(cs *b1, cs *b2, double alpha, double beta,
                           int flag_rel);
cs     *cs_normalize_by_diag_and_release(cs *Q, int flag_rel);
cs     *cs_prod_norm(int mode, cs *A, cs *IhH);
cs     *cs_prod_norm_and_release(cs *b1, cs *lambda, int flag_rel);
int     cs_coarsening(cs *Q,int type,int **indCo,cs **L);
cs     *cs_interpolate(cs *AA,cs *LL,int *indCo);
cs     *cs_triangle(cs *A, int flag_upper, int flag_diag);
void    cs_keypair(const char *key, cs *A, int flag_from_1);
int     cs_scale(cs *C);
int     cs_get_nrow(const cs *A);
int     cs_get_ncol(const cs *A);
double  cs_get_value(const cs *A,int row, int col);
void    cs_set_value(const cs *A,int row, int col, double value);


// Qchol operations
int  qchol_cholesky(int verbose,QChol *QC);
void cs_chol_invert(QChol *qctt,double *xcr,double *rhs, double *work);
void cs_chol_simulate(QChol	*qctt,double *simu,double *work);

// Multigrid operations
cs_MGS *cs_multigrid_manage(cs_MGS *mgs,int mode,
                                        int nlevels, int path_type);
void cs_multigrid_params(cs_MGS *mgs, int flag_cg,
                                     int type_coarse, int ngc, int nmg, int ngs,
                                     double tolgc, double tolnmg);
void cs_multigrid_print(cs_MGS *mgs);
int  cs_multigrid_get_nlevels(cs_MGS *mgs);
int  cs_multigrid_setup(cs_MGS *mgs, QChol *Qctt,
                                    int flag_sel, int verbose, double **sel);
int  cs_multigrid_process(cs_MGS *mgs, QChol *qctt, int verbose,
                                      double *x, double *b, double *work);
void cs_multigrid_coarse2fine(cs_MGS *mgs,double *z,double *work);

Triplet csToTriplet(const cs *A, bool flag_from_1 = false);
String toStringDim(const String& title, const cs *A);
String toStringRange(const String& title, const cs *C);
bool cs_isSymmetric(const cs* A, bool verbose = false, bool detail = false);
bool cs_isDiagonalDominant(cs *A, bool verbose = false, bool detail = false);
bool cs_isDefinitePositive(cs* A, bool verbose = false);

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(Ap,j) (Ap [j] < 0)
#define CS_MARK(Ap,j) { Ap [j] = CS_FLIP (Ap [j]) ; }
#define CS_OVERFLOW(n,size) (n > INT_MAX / (int) size)
#endif
