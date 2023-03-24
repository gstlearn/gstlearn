/*
                                      csparse

Original Author: Timothy Davis
Website: https://people.math.sc.edu/Burkardt/c_src/csparse/csparse.html
License: see doc/csparse_license.txt
*/

/*
Author: Timothy Davis

License:

CSPARSE: a Concise Sparse matrix package.
Copyright (c) 2006, Timothy A. Davis.
http://www.cise.ufl.edu/research/sparse/CSparse

CSPARSE is free software; you can redistribute it and/or modify it under the terms of
the GNU Lesser General Public License as published by the Free Software Foundation;
either version 2.1 of the License, or (at your option) any later version.

CSPARSE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with
this Module; if not, write to the Free Software Foundation,
Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
*/

/*
 Modified by MINES PARIS / PSL (2022)
*/
#ifndef _CS_F_H
#define _CS_F_H

#include "csparse_d.h"

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

MY_EXPORT cs     *cs_add (const cs *A, const cs *B, double alpha, double beta) ;
MY_EXPORT int     cs_cholsol (const cs *A, double *b, int order) ;
MY_EXPORT int     cs_dupl (cs *A) ;
MY_EXPORT int     cs_entry (cs *T, int i, int j, double x) ;
MY_EXPORT int     cs_lusol (const cs *A, double *b, int order, double tol) ;
MY_EXPORT int     cs_gaxpy (const cs *A, const double *x, double *y) ;
MY_EXPORT cs     *cs_multiply (const cs *A, const cs *B) ;
MY_EXPORT int     cs_qrsol (const cs *A, double *b, int order) ;
MY_EXPORT cs     *cs_transpose (const cs *A, int values) ;
MY_EXPORT cs     *cs_triplet (const cs *T) ;
MY_EXPORT double  cs_norm (const cs *A) ;
MY_EXPORT int     cs_print (const cs *A, int brief) ;
MY_EXPORT cs     *cs_load (FILE *f) ;

/* utilities */
MY_EXPORT void   *cs_calloc  (int n, size_t size) ;
MY_EXPORT void   *cs_free    (void *p) ;
MY_EXPORT void   *cs_realloc (void *p, int n, size_t size, int *ok) ;
MY_EXPORT cs     *cs_spalloc (int m, int n, int nzmax, int values, int triplet) ;
MY_EXPORT cs     *cs_spfree  (cs *A) ;
MY_EXPORT int     cs_sprealloc (cs *A, int nzmax) ;
MY_EXPORT void   *cs_malloc  (int n, size_t size) ;

/* --- secondary CSparse routines and data structures ----------------------- */

MY_EXPORT int    *cs_amd (const cs *A, int order) ;
MY_EXPORT csn    *cs_chol (const cs *A, const css *S) ;
MY_EXPORT csd    *cs_dmperm (const cs *A) ;
MY_EXPORT int     cs_droptol (cs *A, double tol) ;
MY_EXPORT int     cs_dropzeros (cs *A) ;
MY_EXPORT int     cs_happly (const cs *V, int i, double beta, double *x) ;
MY_EXPORT int     cs_ipvec (int n, const int *P, const double *b, double *x) ;
MY_EXPORT int     cs_lsolve (const cs *L, double *x) ;
MY_EXPORT int     cs_ltsolve (const cs *L, double *x) ;
MY_EXPORT csn    *cs_lu (const cs *A, const css *S, double tol) ;
MY_EXPORT cs     *cs_permute (const cs *A, const int *P, const int *Q, int values) ;
MY_EXPORT int    *cs_pinv (const int *P, int n) ;
MY_EXPORT int     cs_pvec (int n, const int *P, const double *b, double *x) ;
MY_EXPORT csn    *cs_qr (const cs *A, const css *S) ;
MY_EXPORT css    *cs_schol (const cs *A, int order) ;
MY_EXPORT css    *cs_sqr (const cs *A, int order, int qr) ;
MY_EXPORT cs     *cs_symperm (const cs *A, const int *Pinv, int values) ;
MY_EXPORT int     cs_usolve (const cs *U, double *x) ;
MY_EXPORT int     cs_utsolve (const cs *U, double *x) ;
MY_EXPORT int     cs_updown (cs *L, int sigma, const cs *C, const int *parent) ;
MY_EXPORT css    *cs_sfree (css *S) ;
MY_EXPORT csn    *cs_nfree (csn *N) ;
MY_EXPORT csd    *cs_dfree (csd *D) ;

/* --- tertiary CSparse routines -------------------------------------------- */
MY_EXPORT int    *cs_counts (const cs *A, const int *parent, const int *post, int ata) ;
MY_EXPORT int     cs_cumsum (int *p, int *c, int n) ;
MY_EXPORT int     cs_dfs (int j, cs *L, int top, int *xi, int *pstack, const int *Pinv) ;
MY_EXPORT int    *cs_etree (const cs *A, int ata) ;
MY_EXPORT int     cs_fkeep (cs *A, int (*fkeep) (int, int, double, void *), void *other) ;
MY_EXPORT double  cs_house (double *x, double *beta, int n) ;
MY_EXPORT int    *cs_maxtrans (const cs *A) ;
MY_EXPORT int    *cs_post (int n, const int *parent) ;
MY_EXPORT int     cs_reach (cs *L, const cs *B, int k, int *xi, const int *Pinv) ;
MY_EXPORT csd    *cs_scc (cs *A) ;
MY_EXPORT int     cs_scatter (const cs *A, int j, double beta, int *w, double *x,
                                    int mark, cs *C, int nz) ;
MY_EXPORT int     cs_splsolve (cs *L, const cs *B, int k, int *xi, double *x,
                                     const int *Pinv) ;
MY_EXPORT int     cs_tdfs (int j, int k, int *head, const int *next, int *post,
                                 int *stack) ;

/* utilities */
MY_EXPORT csd    *cs_dalloc (int m, int n) ;
MY_EXPORT cs     *cs_done (cs *C, void *w, void *x, int ok) ;
MY_EXPORT int    *cs_idone (int *p, cs *C, void *w, int ok) ;
MY_EXPORT csn    *cs_ndone (csn *N, cs *C, void *w, void *x, int ok) ;
MY_EXPORT csd    *cs_ddone (csd *D, cs *C, void *w, int ok) ;

/* inversion */
MY_EXPORT int     sparseinv(int n, int *Lp, int *Li, double *Lx, double *d, int *Up,
                            int *Uj, double *Ux, int *Zp, int *Zi, double *Zx, double *z,
                            int *Zdiagp, int *Lmunch);

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(Ap,j) (Ap [j] < 0)
#define CS_MARK(Ap,j) { Ap [j] = CS_FLIP (Ap [j]) ; }
#define CS_OVERFLOW(n,size) (n > INT_MAX / (int) size)

#endif
