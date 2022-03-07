#ifndef _CS_D_H
#define _CS_D_H

#include "gstlearn_export.hpp"

/* --- primary CSparse routines and data structures ------------------------- */
class GSTLEARN_EXPORT cs // cs_sparse    /* matrix in compressed-column or triplet form */
{
public:
    int nzmax ;             /* maximum number of entries */
    int m ;                 /* number of rows */
    int n ;                 /* number of columns */
    int *p ;                /* column pointers (size n+1) or col indices (size nzmax) */
    int *i ;                /* row indices, size nzmax */
    double *x ;             /* numerical values, size nzmax */
    int nz ;                /* # of entries in triplet matrix, -1 for compressed-col */
};

/* --- secondary CSparse routines and data structures ----------------------- */
class GSTLEARN_EXPORT css // cs_symbolic  /* symbolic Cholesky, LU, or QR analysis */
{
public:
    int *Pinv ;            /* inverse row perm. for QR, fill red. perm for Chol */
    int *Q ;               /* fill-reducing column permutation for LU and QR */
    int *parent ;          /* elimination tree for Cholesky and QR */
    int *cp ;              /* column pointers for Cholesky, row counts for QR */
    int m2 ;               /* # of rows for QR, after adding fictitious rows */
    int lnz ;              /* # entries in L for LU or Cholesky; in V for QR */
    int unz ;              /* # entries in U for LU; in R for QR */
};

class GSTLEARN_EXPORT csn // cs_numeric   /* numeric Cholesky, LU, or QR factorization */
{
public:
    cs *L ;                /* L for LU and Cholesky, V for QR */
    cs *U ;                /* U for LU, R for QR, not used for Cholesky */
    int *Pinv ;            /* partial pivoting for LU */
    double *B ;            /* beta [0..n-1] for QR */
};

class GSTLEARN_EXPORT csd // cs_dmperm_results    /* cs_dmperm or cs_scc output */
{
public:
    int *P ;            /* size m, row permutation */
    int *Q ;            /* size n, column permutation */
    int *R ;            /* size nb+1, block k is rows R[k] to R[k+1]-1 in A(P,Q) */
    int *S ;            /* size nb+1, block k is cols S[k] to S[k+1]-1 in A(P,Q) */
    int nb ;            /* # of blocks in fine dmperm decomposition */
    int rr [5] ;        /* coarse row decomposition */
    int cc [5] ;        /* coarse column decomposition */
};

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

class GSTLEARN_EXPORT Triplet
{
public:
    bool flagFromOne;
    int nrows;
    int ncols;
    VectorInt rows;
    VectorInt cols;
    VectorDouble values;
};

#endif
