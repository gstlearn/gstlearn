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

#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"

class GSTLEARN_EXPORT Triplet
{
public:
    bool flagFromOne;
    int nrows;
    int ncols;
    VectorInt rows;
    VectorInt cols;
    VectorDouble values;

public:
    VectorDouble getValues() const { return values; }
    VectorInt getRows() const { return rows; }
    VectorInt getCols() const { return cols; }

    /// Has a specific implementation in the Target language
    DECLARE_TOTL;
};

#ifndef SWIG

class cs;
class css;
class csn;
class csd;

class GSTLEARN_EXPORT QChol
{
public:
    cs  *Q;
    css *S;
    csn *N;
};

class GSTLEARN_EXPORT cs_MG
{
public:
    int      nh;
    int      nH;
    double  *sumrow;
    cs      *IhH;
    QChol   *A;
};

class GSTLEARN_EXPORT cs_MGS
{
public:
    int       flag_cg;            /* Apply Conjugate-gradient */
    int       nlevels;            /* Number of multigrid levels */
    int       npath;              /* Number of paths */
    int       type_coarse;        /* Type of coarsening algorithm */
    int       ngc;                /* Maximum count of Conjugate-Gradient iters */
    int       nmg;                /* Maximum count of Multigrid operations*/
    int       ngs;                /* Number of Gauss_Siedel relaxation cycles */
    int       ncur;               /* Number of mesh vertices */
    int      *path;               /* Path description */
    double    tolnmg;             /* Tolerance for Multigrid */
    double    tolcg;              /* Tolerance for Conjugate-Gradient */
    double   *diag;               /* Normation diagonal */
    cs_MG   **mg;                 /* Array of cs_MG structures */
};

/// TODO : cs_*2 functions to be removed (encapsulation)
GSTLEARN_EXPORT cs     *cs_spfree2(cs *A);
GSTLEARN_EXPORT cs     *cs_spalloc2(int m, int n, int nzmax, int values, int triplet);
GSTLEARN_EXPORT int     cs_entry2(cs *T, int i, int j, double x);
GSTLEARN_EXPORT cs     *cs_triplet2(const cs *T);
GSTLEARN_EXPORT cs     *cs_transpose2(const cs *A, int values);
GSTLEARN_EXPORT cs     *cs_multiply2(const cs *A, const cs *B);
GSTLEARN_EXPORT int     cs_print2(const cs *A, int brief);

GSTLEARN_EXPORT void    cs_add_cste(cs *A, double value);
GSTLEARN_EXPORT void    cs_set_cste(cs *A, double value);

GSTLEARN_EXPORT void    cs_tmulvec(const cs *A, int nout, const double *x, double *y);
GSTLEARN_EXPORT void    cs_mulvec(const cs *A, int nout, const double *x, double *y);
GSTLEARN_EXPORT void    cs_vecmult(const cs *A, int nout, const double *x, double *y);
GSTLEARN_EXPORT void    cs_mulvec_uptri(const cs *A, int nout,
                                        const double *x, double *y, int flag_diag);
GSTLEARN_EXPORT void    cs_mulvec_lowtri(const cs *A, int nout,
                                         const double *x, double *y, int flag_diag);
GSTLEARN_EXPORT cs     *cs_matvecR(const cs *A, const double *x, int oper);
GSTLEARN_EXPORT cs     *cs_matvecL(const cs *A, const double *x, int oper);
GSTLEARN_EXPORT cs     *cs_matvecnorm(const cs *A, const double *x, int oper);
GSTLEARN_EXPORT void    cs_matvecnorm_inplace(cs *A, const double *x, int oper);
GSTLEARN_EXPORT double *cs_col_sumrow(const cs *A,int *ncol,int *nrow);
GSTLEARN_EXPORT double  cs_maxsumabscol(const cs *A);

GSTLEARN_EXPORT void    cs_print_only(const char *title, const cs *A,int nlimit);
GSTLEARN_EXPORT void    cs_print_nice (const char *title,const cs *A, int maxrow=-1, int maxcol=-1);
GSTLEARN_EXPORT void    cs_print_dim(const char *title,const cs *A);
GSTLEARN_EXPORT void    cs_print_short(const char *title, const cs *L, int nmax);
GSTLEARN_EXPORT void    cs_print_file(const char *radix, int rank, cs *A);
GSTLEARN_EXPORT cs     *cs_compress(cs *A);

GSTLEARN_EXPORT void    cs_force_dimension(cs *T, int nrow, int ncol);
GSTLEARN_EXPORT cs*     cs_diag(VectorDouble diag, double tol = EPSILON10);

GSTLEARN_EXPORT VectorInt cs_color_coding(cs *Q,int start,int *ncolor);
GSTLEARN_EXPORT cs     *cs_extract_submatrix_by_color(cs *C, const VectorInt& colors,
                                                      int ref_color, int row_ok, int col_ok);
GSTLEARN_EXPORT VectorDouble csd_extract_diag_VD(cs *C, int mode);
GSTLEARN_EXPORT cs     *cs_prod_norm_diagonal(int mode, const cs *B, VectorDouble diag);
GSTLEARN_EXPORT cs     *cs_arrays_to_sparse(int n, int nrow, int ncol,
                                            const double *rows,
                                            const double *cols,
                                            const double *vals);
GSTLEARN_EXPORT cs     *cs_vectors_to_sparse(int nrow,
                                             int ncol,
                                             const VectorDouble &rows,
                                             const VectorDouble &cols,
                                             const VectorDouble &values);
GSTLEARN_EXPORT int     cs_lsolve_lowtri( const cs *L, const double *x, double *y);
GSTLEARN_EXPORT int     cs_lsolve_uptri (const cs *L, const double *x, double *y);
GSTLEARN_EXPORT cs     *cs_invert(const cs *A, int order, double epsilon = EPSILON6);

// Multigrid operations
GSTLEARN_EXPORT cs_MGS *cs_multigrid_manage(cs_MGS *mgs,int mode,
                                            int nlevels, int path_type);
GSTLEARN_EXPORT void    cs_multigrid_params(cs_MGS *mgs, int flag_cg,
                                            int type_coarse, int ngc, int nmg, int ngs,
                                            double tolgc, double tolnmg);
GSTLEARN_EXPORT void    cs_multigrid_print(cs_MGS *mgs);
GSTLEARN_EXPORT int     cs_multigrid_get_nlevels(cs_MGS *mgs);
GSTLEARN_EXPORT int     cs_multigrid_setup(cs_MGS *mgs, QChol *Qctt,
                                           int flag_sel, int verbose, double **sel);
GSTLEARN_EXPORT int     cs_multigrid_process(cs_MGS *mgs, QChol *qctt, int verbose,
                                             double *x, double *b, double *work);
GSTLEARN_EXPORT void    cs_multigrid_coarse2fine(cs_MGS *mgs,double *z,double *work);

// Qchol operations
GSTLEARN_EXPORT int     qchol_cholesky(int verbose,QChol *QC);
GSTLEARN_EXPORT void    cs_chol_invert(QChol *qctt,double *xcr,double *rhs, double *work);
GSTLEARN_EXPORT void    cs_chol_simulate(QChol *qctt,double *simu,double *work);

// Multigrid operations
GSTLEARN_EXPORT cs_MGS *cs_multigrid_manage(cs_MGS *mgs,int mode,
                                            int nlevels, int path_type);
GSTLEARN_EXPORT void    cs_multigrid_params(cs_MGS *mgs, int flag_cg,
                                            int type_coarse, int ngc, int nmg, int ngs,
                                            double tolgc, double tolnmg);
GSTLEARN_EXPORT void    cs_multigrid_print(cs_MGS *mgs);
GSTLEARN_EXPORT int     cs_multigrid_get_nlevels(cs_MGS *mgs);
GSTLEARN_EXPORT int     cs_multigrid_setup(cs_MGS *mgs, QChol *Qctt,
                                           int flag_sel, int verbose, double **sel);
GSTLEARN_EXPORT int     cs_multigrid_process(cs_MGS *mgs, QChol *qctt, int verbose,
                                           double *x, double *b, double *work);
GSTLEARN_EXPORT void    cs_multigrid_coarse2fine(cs_MGS *mgs,double *z,double *work);

//
GSTLEARN_EXPORT Triplet csToTriplet(const cs *A, bool flag_from_1 = false);
GSTLEARN_EXPORT String  toStringDim(const String& title, const cs *A);
GSTLEARN_EXPORT String  toStringRange(const String& title, const cs *C);
GSTLEARN_EXPORT bool    cs_isSymmetric(const cs* A, bool verbose = false, bool detail = false);
GSTLEARN_EXPORT bool    cs_isDiagonalDominant(cs *A, bool verbose = false, bool detail = false);
GSTLEARN_EXPORT bool    cs_isDefinitePositive(cs* A, bool verbose = false);

GSTLEARN_EXPORT cs     *cs_extract_submatrix_by_ranks(cs *C, int *row_array, int *col_array);

GSTLEARN_EXPORT void    cs_sparse_to_triplet (const cs *A, int flag_from_1,
                                              int *number,int **cols,int **rows,double ** vals);
GSTLEARN_EXPORT cs     *cs_extract_submatrix(cs *C,
                                             int row_from, int row_length,
                                             int col_from, int col_length);
GSTLEARN_EXPORT void    cs_print_range(const char *title,const cs *C);
GSTLEARN_EXPORT cs     *cs_eye(int number,double value);
GSTLEARN_EXPORT cs     *cs_eye_tab(int number, double *values);
GSTLEARN_EXPORT cs     *cs_extract_diag(cs *C,int mode);
GSTLEARN_EXPORT void    cs_diag_suppress(cs *C);
GSTLEARN_EXPORT double *csd_extract_diag(cs *C,int mode);
GSTLEARN_EXPORT int     cs_sort_i(cs *C);

GSTLEARN_EXPORT void    cs_rowcol(const cs *A,int *nrows,int *ncols,int *count,double *percent);
GSTLEARN_EXPORT cs     *cs_duplicate(const cs *b1);
GSTLEARN_EXPORT cs     *cs_multiply_and_release(cs *b1,cs *b2,int flag_rel);
GSTLEARN_EXPORT cs     *cs_add_and_release(cs *b1, cs *b2, double alpha, double beta,
                                           int flag_rel);
GSTLEARN_EXPORT cs     *cs_normalize_by_diag_and_release(cs *Q, int flag_rel);
GSTLEARN_EXPORT cs     *cs_prod_norm(int mode, const cs *A, const cs *IhH);
GSTLEARN_EXPORT cs     *cs_prod_norm_single(int mode, cs *B);
GSTLEARN_EXPORT cs     *cs_prod_norm_and_release(cs *b1, cs *lambda, int flag_rel);
GSTLEARN_EXPORT int     cs_coarsening(cs *Q,int type,int **indCo,cs **L);
GSTLEARN_EXPORT cs     *cs_interpolate(cs *AA,cs *LL,int *indCo);
GSTLEARN_EXPORT cs     *cs_triangle(cs *A, int flag_upper, int flag_diag);
GSTLEARN_EXPORT void    cs_keypair(const char *key, cs *A, int flag_from_1);
GSTLEARN_EXPORT int     cs_scale(cs *C);
GSTLEARN_EXPORT int     cs_get_nrow(const cs *A);
GSTLEARN_EXPORT int     cs_get_ncol(const cs *A);
GSTLEARN_EXPORT int     cs_get_ncell(const cs *A);
GSTLEARN_EXPORT bool    cs_exist(const cs* A, int row, int col);
GSTLEARN_EXPORT double  cs_get_value(const cs *A,int row, int col);
GSTLEARN_EXPORT void    cs_set_value(const cs *A,int row, int col, double value);
GSTLEARN_EXPORT void    cs_add_value(const cs *A, int row, int col, double value);
GSTLEARN_EXPORT double* cs_toArray(const cs *A);
GSTLEARN_EXPORT cs*     cs_strip(cs *A, double eps, int hypothesis = 3, bool verbose = false);
GSTLEARN_EXPORT int     cs_nnz(const cs* A);
GSTLEARN_EXPORT cs*     cs_glue(const cs*A1, const cs* A2, bool shiftRow, bool shiftCol);

GSTLEARN_EXPORT void    cs_set_status_update_nonzero_value(int status = 2);
GSTLEARN_EXPORT int     cs_get_status_update_nonzero_value();

#endif // SWIG
