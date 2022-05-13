/*
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
 Modified by ARMINES-MINES Paris (2022)
*/

#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Basic/File.hpp"
#include "Basic/String.hpp"
#include "Basic/OptDbg.hpp"
#include "csparse_d.h"
#include "csparse_f.h"

#include <math.h>
#include <limits.h>
#include <set>

/* Global symbols for Csparse */

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif
#define CS_VER 1                   /* CSparse Version 1.2.0 */
#define CS_SUBVER 2
#define CS_SUBSUB 0
#define CS_DATE "Mar 6, 2006"     /* CSparse release date */
#define CS_COPYRIGHT "Copyright (c) Timothy A. Davis, 2006"
#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(Ap,j) (Ap [j] < 0)
#define CS_MARK(Ap,j) { Ap [j] = CS_FLIP (Ap [j]) ; }
#define CS_OVERFLOW(n,size) (n > INT_MAX / (int) size)

#define MAX_NEIGH 100
#define XCR(ilevel,i)       (xcr[(ilevel) * ncur + (i)])
#define RHS(ilevel,i)       (rhs[(ilevel) * ncur + (i)])
#define MAT(i,j)            (mat[(i) * n + (j)])
#define DEBUG 0


/*
 Purpose:

 CS_ADD computes C = alpha*A + beta*B for sparse A and B.

 Reference:

 Timothy Davis,
 Direct Methods for Sparse Linear Systems,
 SIAM, Philadelphia, 2006.
 */
cs* cs_add(const cs *A, const cs *B, double alpha, double beta)
{
  int p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values;
  double *x, *Bx, *Cx;
  cs *C;
  if (!A || !B) return (NULL); /* check inputs */
  m = A->m;
  anz = A->p[A->n];
  n = B->n;
  Bp = B->p;
  Bx = B->x;
  bnz = Bp[n];
  w = (int*) cs_calloc(m, sizeof(int));
  values = (A->x != NULL) && (Bx != NULL);
  x = values ? (double*) cs_malloc(m, sizeof(double)) : NULL;
  C = cs_spalloc(m, n, anz + bnz, values, 0);
  if (!C || !w || (values && !x))
  {
    messerr("Core allocation problem in CSparse Library (%d x %d)", m, n);
    return (cs_done(C, w, x, 0));
  }
  Cp = C->p;
  Ci = C->i;
  Cx = C->x;
  for (j = 0; j < n; j++)
  {
    Cp[j] = nz; /* column j of C starts here */
    nz = cs_scatter(A, j, alpha, w, x, j + 1, C, nz); /* alpha*A(:,j)*/
    nz = cs_scatter(B, j, beta, w, x, j + 1, C, nz); /* beta*B(:,j) */
    if (values) for (p = Cp[j]; p < nz; p++)
      Cx[p] = x[Ci[p]];
  }
  Cp[n] = nz; /* finalize the last column of C */
  cs_sprealloc(C, 0); /* remove extra space from C */
  return (cs_done(C, w, x, 1)); /* success; free workspace, return C */
}

/*
 Purpose:

 CS_WCLEAR clears W.

 Reference:

 Timothy Davis,
 Direct Methods for Sparse Linear Systems,
 SIAM, Philadelphia, 2006.
 */
static int cs_wclear(int mark, int lemax, int *w, int n)
{
  int k;
  if (mark < 2 || (mark + lemax < 0))
  {
    for (k = 0; k < n; k++)
      if (w[k] != 0) w[k] = 1;
    mark = 2;
  }
  return (mark); /* at this point, w [0..n-1] < mark holds */
}

/* keep off-diagonal entries; drop diagonal entries */
static int cs_diag(int i, int j, double /*aij*/, void* /*other*/)
{
  return (i != j);
}

/* p = amd(A+A') if symmetric is true, or amd(A'A) otherwise */
/*
 Purpose:

 CS_AMD carries out the approximate minimum degree algorithm.

 Reference:

 Timothy Davis,
 Direct Methods for Sparse Linear Systems,
 SIAM, Philadelphia, 2006.

 Parameters:

 Input, int ORDER:
 -1:natural, <other
 0:Cholesky,
 1:LU,
 2:QR
 */
int* cs_amd(const cs *A, int order)
{
  cs *C, *A2, *AT;
  int *Cp, *Ci, *last, *ww, *len, *nv, *next, *P, *head, *elen, *degree, *w,
      *hhead, *ATp, *ATi, d, dk, dext, lemax = 0, e, elenk, eln, i, j, k, k1,
      k2, k3, jlast, ln, dense, nzmax, mindeg = 0, nvi, nvj, nvk, mark, wnvi,
      ok, cnz, nel = 0, p, p1, p2, p3, p4, pj, pk, pk1, pk2, pn, q, n, m;
  unsigned int h;
  /* --- Construct matrix C ----------------------------------------------- */
  if (!A || order < 0) return (NULL); /* check inputs; quick return */
  AT = cs_transpose(A, 0); /* compute A' */
  if (!AT) return (NULL);
  m = A->m;
  n = A->n;
  dense = CS_MAX(16, (int) (10 * sqrt((double) n))); /* find dense threshold */
  dense = CS_MIN(n - 2, dense);
  if (order == 0 && n == m)
  {
    C = cs_add(A, AT, 0, 0); /* C = A+A' */
  }
  else if (order == 1)
  {
    ATp = AT->p; /* drop dense columns from AT */
    ATi = AT->i;
    for (p2 = 0, j = 0; j < m; j++)
    {
      p = ATp[j]; /* column j of AT starts here */
      ATp[j] = p2; /* new column j starts here */
      if (ATp[j + 1] - p > dense) continue; /* skip dense col j */
      for (; p < ATp[j + 1]; p++)
        ATi[p2++] = ATi[p];
    }
    ATp[m] = p2; /* finalize AT */
    A2 = cs_transpose(AT, 0); /* A2 = AT' */
    C = A2 ? cs_multiply(AT, A2) : NULL; /* C=A'*A with no dense rows */
    cs_spfree(A2);
  }
  else
  {
    C = cs_multiply(AT, A); /* C=A'*A */
  }
  cs_spfree(AT);
  if (!C) return (NULL);
  P = (int*) cs_malloc(n + 1, sizeof(int)); /* allocate result */
  ww = (int*) cs_malloc(8 * (n + 1), sizeof(int));/* get workspace */
  len = ww;
  nv = ww + (n + 1);
  next = ww + 2 * (n + 1);
  head = ww + 3 * (n + 1);
  elen = ww + 4 * (n + 1);
  degree = ww + 5 * (n + 1);
  w = ww + 6 * (n + 1);
  hhead = ww + 7 * (n + 1);
  last = P; /* use P as workspace for last */
  cs_fkeep(C, &cs_diag, NULL); /* drop diagonal entries */
  Cp = C->p;
  cnz = Cp[n];
  if (!cs_sprealloc(C, cnz + cnz / 5 + 2 * n)) return (cs_idone(P, C, ww, 0));
  /* --- Initialize quotient graph ---------------------------------------- */
  for (k = 0; k < n; k++)
    len[k] = Cp[k + 1] - Cp[k];
  len[n] = 0;
  nzmax = C->nzmax;
  Ci = C->i;
  for (i = 0; i <= n; i++)
  {
    head[i] = -1; /* degree list i is empty */
    last[i] = -1;
    next[i] = -1;
    hhead[i] = -1; /* hash list i is empty */
    nv[i] = 1; /* node i is just one node */
    w[i] = 1; /* node i is alive */
    elen[i] = 0; /* Ek of node i is empty */
    degree[i] = len[i]; /* degree of node i */
  }
  mark = cs_wclear(0, 0, w, n); /* clear w */
  elen[n] = -2; /* n is a dead element */
  Cp[n] = -1; /* n is a root of assembly tree */
  w[n] = 0; /* n is a dead element */
  /* --- Initialize degree lists ------------------------------------------ */
  for (i = 0; i < n; i++)
  {
    d = degree[i];
    if (d == 0) /* node i is empty */
    {
      elen[i] = -2; /* element i is dead */
      nel++;
      Cp[i] = -1; /* i is a root of assemby tree */
      w[i] = 0;
    }
    else if (d > dense) /* node i is dense */
    {
      nv[i] = 0; /* absorb i into element n */
      elen[i] = -1; /* node i is dead */
      nel++;
      Cp[i] = CS_FLIP(n);
      nv[n]++;
    }
    else
    {
      if (head[d] != -1) last[head[d]] = i;
      next[i] = head[d]; /* put node i in degree list d */
      head[d] = i;
    }
  }
  while (nel < n) /* while (selecting pivots) do */
  {
    /* --- Select node of minimum approximate degree -------------------- */
    for (k = -1; mindeg < n && (k = head[mindeg]) == -1; mindeg++)
      ;
    if (next[k] != -1) last[next[k]] = -1;
    head[mindeg] = next[k]; /* remove k from degree list */
    elenk = elen[k]; /* elenk = |Ek| */
    nvk = nv[k]; /* # of nodes k represents */
    nel += nvk; /* nv[k] nodes of A eliminated */
    /* --- Garbage collection ------------------------------------------- */
    if (elenk > 0 && cnz + mindeg >= nzmax)
    {
      for (j = 0; j < n; j++)
      {
        if ((p = Cp[j]) >= 0) /* j is a live node or element */
        {
          Cp[j] = Ci[p]; /* save first entry of object */
          Ci[p] = CS_FLIP(j); /* first entry is now CS_FLIP(j) */
        }
      }
      for (q = 0, p = 0; p < cnz;) /* scan all of memory */
      {
        if ((j = CS_FLIP(Ci[p++])) >= 0) /* found object j */
        {
          Ci[q] = Cp[j]; /* restore first entry of object */
          Cp[j] = q++; /* new pointer to object j */
          for (k3 = 0; k3 < len[j] - 1; k3++)
            Ci[q++] = Ci[p++];
        }
      }
      cnz = q; /* Ci [cnz...nzmax-1] now free */
    }
    /* --- Construct new element ---------------------------------------- */
    dk = 0;
    nv[k] = -nvk; /* flag k as in Lk */
    p = Cp[k];
    pk1 = (elenk == 0) ? p : cnz; /* do in place if elen[k] == 0 */
    pk2 = pk1;
    for (k1 = 1; k1 <= elenk + 1; k1++)
    {
      if (k1 > elenk)
      {
        e = k; /* search the nodes in k */
        pj = p; /* list of nodes starts at Ci[pj]*/
        ln = len[k] - elenk; /* length of list of nodes in k */
      }
      else
      {
        e = Ci[p++]; /* search the nodes in e */
        pj = Cp[e];
        ln = len[e]; /* length of list of nodes in e */
      }
      for (k2 = 1; k2 <= ln; k2++)
      {
        i = Ci[pj++];
        if ((nvi = nv[i]) <= 0) continue; /* node i dead, or seen */
        dk += nvi; /* degree[Lk] += size of node i */
        nv[i] = -nvi; /* negate nv[i] to denote i in Lk*/
        Ci[pk2++] = i; /* place i in Lk */
        if (next[i] != -1) last[next[i]] = last[i];
        if (last[i] != -1) /* remove i from degree list */
        {
          next[last[i]] = next[i];
        }
        else
        {
          head[degree[i]] = next[i];
        }
      }
      if (e != k)
      {
        Cp[e] = CS_FLIP(k); /* absorb e into k */
        w[e] = 0; /* e is now a dead element */
      }
    }
    if (elenk != 0) cnz = pk2; /* Ci [cnz...nzmax] is free */
    degree[k] = dk; /* external degree of k - |Lk\i| */
    Cp[k] = pk1; /* element k is in Ci[pk1..pk2-1] */
    len[k] = pk2 - pk1;
    elen[k] = -2; /* k is now an element */
    /* --- Find set differences ----------------------------------------- */
    mark = cs_wclear(mark, lemax, w, n); /* clear w if necessary */
    for (pk = pk1; pk < pk2; pk++) /* scan 1: find |Le\Lk| */
    {
      i = Ci[pk];
      if ((eln = elen[i]) <= 0) continue;/* skip if elen[i] empty */
      nvi = -nv[i]; /* nv [i] was negated */
      wnvi = mark - nvi;
      for (p = Cp[i]; p <= Cp[i] + eln - 1; p++) /* scan Ei */
      {
        e = Ci[p];
        if (w[e] >= mark)
        {
          w[e] -= nvi; /* decrement |Le\Lk| */
        }
        else if (w[e] != 0) /* ensure e is a live element */
        {
          w[e] = degree[e] + wnvi; /* 1st time e seen in scan 1 */
        }
      }
    }
    /* --- Degree update ------------------------------------------------ */
    for (pk = pk1; pk < pk2; pk++) /* scan2: degree update */
    {
      i = Ci[pk]; /* consider node i in Lk */
      p1 = Cp[i];
      p2 = p1 + elen[i] - 1;
      pn = p1;
      for (h = 0, d = 0, p = p1; p <= p2; p++) /* scan Ei */
      {
        e = Ci[p];
        if (w[e] != 0) /* e is an unabsorbed element */
        {
          dext = w[e] - mark; /* dext = |Le\Lk| */
          if (dext > 0)
          {
            d += dext; /* sum up the set differences */
            Ci[pn++] = e; /* keep e in Ei */
            h += e; /* compute the hash of node i */
          }
          else
          {
            Cp[e] = CS_FLIP(k); /* aggressive absorb. e->k */
            w[e] = 0; /* e is a dead element */
          }
        }
      }
      elen[i] = pn - p1 + 1; /* elen[i] = |Ei| */
      p3 = pn;
      p4 = p1 + len[i];
      for (p = p2 + 1; p < p4; p++) /* prune edges in Ai */
      {
        j = Ci[p];
        if ((nvj = nv[j]) <= 0) continue; /* node j dead or in Lk */
        d += nvj; /* degree(i) += |j| */
        Ci[pn++] = j; /* place j in node list of i */
        h += j; /* compute hash for node i */
      }
      if (d == 0) /* check for mass elimination */
      {
        Cp[i] = CS_FLIP(k); /* absorb i into k */
        nvi = -nv[i];
        dk -= nvi; /* |Lk| -= |i| */
        nvk += nvi; /* |k| += nv[i] */
        nel += nvi;
        nv[i] = 0;
        elen[i] = -1; /* node i is dead */
      }
      else
      {
        degree[i] = CS_MIN(degree[i], d); /* update degree(i) */
        Ci[pn] = Ci[p3]; /* move first node to end */
        Ci[p3] = Ci[p1]; /* move 1st el. to end of Ei */
        Ci[p1] = k; /* add k as 1st element in of Ei */
        len[i] = pn - p1 + 1; /* new len of adj. list of node i */
        h %= n; /* finalize hash of i */
        next[i] = hhead[h]; /* place i in hash bucket */
        hhead[h] = i;
        last[i] = h; /* save hash of i in last[i] */
      }
    } /* scan2 is done */
    degree[k] = dk; /* finalize |Lk| */
    lemax = CS_MAX(lemax, dk);
    mark = cs_wclear(mark + lemax, lemax, w, n); /* clear w */
    /* --- Supernode detection ------------------------------------------ */
    for (pk = pk1; pk < pk2; pk++)
    {
      i = Ci[pk];
      if (nv[i] >= 0) continue; /* skip if i is dead */
      h = last[i]; /* scan hash bucket of node i */
      i = hhead[h];
      hhead[h] = -1; /* hash bucket will be empty */
      for (; i != -1 && next[i] != -1; i = next[i], mark++)
      {
        ln = len[i];
        eln = elen[i];
        for (p = Cp[i] + 1; p <= Cp[i] + ln - 1; p++)
          w[Ci[p]] = mark;
        jlast = i;
        for (j = next[i]; j != -1;) /* compare i with all j */
        {
          ok = (len[j] == ln) && (elen[j] == eln);
          for (p = Cp[j] + 1; ok && p <= Cp[j] + ln - 1; p++)
          {
            if (w[Ci[p]] != mark) ok = 0; /* compare i and j*/
          }
          if (ok) /* i and j are identical */
          {
            Cp[j] = CS_FLIP(i); /* absorb j into i */
            nv[i] += nv[j];
            nv[j] = 0;
            elen[j] = -1; /* node j is dead */
            j = next[j]; /* delete j from hash bucket */
            next[jlast] = j;
          }
          else
          {
            jlast = j; /* j and i are different */
            j = next[j];
          }
        }
      }
    }
    /* --- Finalize new element------------------------------------------ */
    for (p = pk1, pk = pk1; pk < pk2; pk++) /* finalize Lk */
    {
      i = Ci[pk];
      if ((nvi = -nv[i]) <= 0) continue;/* skip if i is dead */
      nv[i] = nvi; /* restore nv[i] */
      d = degree[i] + dk - nvi; /* compute external degree(i) */
      d = CS_MIN(d, n - nel - nvi);
      if (head[d] != -1) last[head[d]] = i;
      next[i] = head[d]; /* put i back in degree list */
      last[i] = -1;
      head[d] = i;
      mindeg = CS_MIN(mindeg, d); /* find new minimum degree */
      degree[i] = d;
      Ci[p++] = i; /* place i in Lk */
    }
    nv[k] = nvk; /* # nodes absorbed into k */
    if ((len[k] = p - pk1) == 0) /* length of adj list of element k*/
    {
      Cp[k] = -1; /* k is a root of the tree */
      w[k] = 0; /* k is now a dead element */
    }
    if (elenk != 0) cnz = p; /* free unused space in Lk */
  }
  /* --- Postordering ----------------------------------------------------- */
  for (i = 0; i < n; i++)
    Cp[i] = CS_FLIP(Cp[i]);/* fix assembly tree */
  for (j = 0; j <= n; j++)
    head[j] = -1;
  for (j = n; j >= 0; j--) /* place unordered nodes in lists */
  {
    if (nv[j] > 0) continue; /* skip if j is an element */
    next[j] = head[Cp[j]]; /* place j in list of its parent */
    head[Cp[j]] = j;
  }
  for (e = n; e >= 0; e--) /* place elements in lists */
  {
    if (nv[e] <= 0) continue; /* skip unless e is an element */
    if (Cp[e] != -1)
    {
      next[e] = head[Cp[e]]; /* place e in list of its parent */
      head[Cp[e]] = e;
    }
  }
  for (k = 0, i = 0; i <= n; i++) /* postorder the assembly tree */
  {
    if (Cp[i] == -1) k = cs_tdfs(i, k, head, next, P, w);
  }
  return (cs_idone(P, C, ww, 1));
}

/* compute nonzero pattern of L(k,:) */
static int cs_ereach(const cs *A,
                     int k,
                     const int *parent,
                     int *s,
                     int *w,
                     double *x,
                     int top)
{
  int i, p, len, *Ap = A->p, *Ai = A->i;
  double *Ax = A->x;
  for (p = Ap[k]; p < Ap[k + 1]; p++) /* get pattern of L(k,:) */
  {
    i = Ai[p]; /* A(i,k) is nonzero */
    if (i > k) continue; /* only use upper triangular part of A */
    x[i] = Ax[p]; /* x(i) = A(i,k) */
    for (len = 0; w[i] != k; i = parent[i]) /* traverse up etree */
    {
      s[len++] = i; /* L(k,i) is nonzero */
      w[i] = k; /* mark i as visited */
    }
    while (len > 0)
      s[--top] = s[--len]; /* push path onto stack */
  }
  return (top); /* s [top..n-1] contains pattern of L(k,:)*/
}

/* L = chol (A, [Pinv parent cp]), Pinv is optional */
csn* cs_chol(const cs *A, const css *S)
{
  double d, lki, *Lx, *x;
  int top, i, p, k, n, *Li, *Lp, *cp, *Pinv, *w, *s, *c, *parent;
  cs *L, *C, *E;
  csn *N;
  if (!A || !S || !S->cp || !S->parent) return (NULL); /* check inputs */
  n = A->n;
  N = (csn*) cs_calloc(1, sizeof(csn));
  w = (int*) cs_malloc(3 * n, sizeof(int));
  s = w + n, c = w + 2 * n;
  x = (double*) cs_malloc(n, sizeof(double));
  cp = S->cp;
  Pinv = S->Pinv;
  parent = S->parent;
  C = Pinv ? cs_symperm(A, Pinv, 1) : ((cs*) A);
  E = Pinv ? C : NULL;
  if (!N || !w || !x || !C) return (cs_ndone(N, E, w, x, 0));
  N->L = L = cs_spalloc(n, n, cp[n], 1, 0);
  if (!L)
  {
    messerr("Core allocation problem in CSparse Library (%d x %d)", n, n);
    return (cs_ndone(N, E, w, x, 0));
  }

  Lp = L->p;
  Li = L->i;
  Lx = L->x;
  for (k = 0; k < n; k++)
  {
    /* --- Nonzero pattern of L(k,:) ------------------------------------ */
    Lp[k] = c[k] = cp[k]; /* column k of L starts here */
    x[k] = 0; /* x (0:k) is now zero */
    w[k] = k; /* mark node k as visited */
    top = cs_ereach(C, k, parent, s, w, x, n); /* find row k of L*/
    d = x[k]; /* d = C(k,k) */
    x[k] = 0; /* clear workspace for k+1st iteration */
    /* --- Triangular solve --------------------------------------------- */
    for (; top < n; top++) /* solve L(0:k-1,0:k-1) * x = C(:,k) */
    {
      i = s[top]; /* s [top..n-1] is pattern of L(k,:) */
      lki = x[i] / Lx[Lp[i]]; /* L(k,i) = x (i) / L(i,i) */
      x[i] = 0; /* clear workspace for k+1st iteration */
      for (p = Lp[i] + 1; p < c[i]; p++)
      {
        x[Li[p]] -= Lx[p] * lki;
      }
      d -= lki * lki; /* d = d - L(k,i)*L(k,i) */
      p = c[i]++;
      Li[p] = k; /* store L(k,i) in column i */
      Lx[p] = lki;
    }
    /* --- Compute L(k,k) ----------------------------------------------- */
    if (d <= 0)
    {
      messerr("Pivot line %d/%d should be positive (d= %lf)", k + 1, n, d);
      return (cs_ndone(N, E, w, x, 0)); /* not positive definite */
    }
    p = c[k]++;
    Li[p] = k; /* store L(k,k) = sqrt(d) in column k */
    Lx[p] = sqrt(d);
  }
  Lp[n] = cp[n]; /* finalize L */
  return (cs_ndone(N, E, w, x, 1)); /* success: free E,w,x; return N */
}

/* x=A\b where A is symmetric positive definite; b overwritten with solution */
int cs_cholsol(const cs *A, double *b, int order)
{
  double *x;
  css *S;
  csn *N;
  int n, ok;
  if (!A || !b) return (0); /* check inputs */
  n = A->n;
  S = cs_schol(A, order); /* ordering and symbolic analysis */
  N = cs_chol(A, S); /* numeric Cholesky factorization */
  x = (double*) cs_malloc(n, sizeof(double));
  ok = (S && N && x);
  if (ok)
  {
    cs_ipvec(n, S->Pinv, b, x); /* x = P*b */
    cs_lsolve(N->L, x); /* x = L\x */
    cs_ltsolve(N->L, x); /* x = L'\x */
    cs_pvec(n, S->Pinv, x, b); /* b = P'*x */
  }
  cs_free(x);
  cs_sfree(S);
  cs_nfree(N);
  return (ok);
}

/* process edge (j,i) of the matrix */
static void cs_cedge(int j,
                     int i,
                     const int *first,
                     int *maxfirst,
                     int *delta,
                     int *prevleaf,
                     int *ancestor)
{
  int q, s, sparent, jprev;
  if (i <= j || first[j] <= maxfirst[i]) return;
  maxfirst[i] = first[j]; /* update max first[j] seen so far */
  jprev = prevleaf[i]; /* j is a leaf of the ith subtree */
  delta[j]++; /* A(i,j) is in the skeleton matrix */
  if (jprev != -1)
  {
    /* q = least common ancestor of jprev and j */
    for (q = jprev; q != ancestor[q]; q = ancestor[q])
      ;
    for (s = jprev; s != q; s = sparent)
    {
      sparent = ancestor[s]; /* path compression */
      ancestor[s] = q;
    }
    delta[q]--; /* decrement to account for overlap in q */
  }
  prevleaf[i] = j; /* j is now previous leaf of ith subtree */
}

/* colcount = column counts of LL'=A or LL'=A'A, given parent & post ordering */
int* cs_counts(const cs *A, const int *parent, const int *post, int ata)
{
  int i, j, k, p, n, m, ii, s, *ATp, *ATi, *maxfirst, *prevleaf, *ancestor,
      *head = NULL, *next = NULL, *colcount, *w, *first, *delta;
  cs *AT;
  if (!A || !parent || !post) return (NULL); /* check inputs */
  m = A->m;
  n = A->n;
  s = 4 * n + (ata ? (n + m + 1) : 0);
  w = (int*) cs_malloc(s, sizeof(int));
  first = w + 3 * n; /* get workspace */
  ancestor = w;
  maxfirst = w + n;
  prevleaf = w + 2 * n;
  delta = colcount = (int*) cs_malloc(n, sizeof(int)); /* allocate result */
  AT = cs_transpose(A, 0);
  if (!AT || !colcount || !w) return (cs_idone(colcount, AT, w, 1));
  for (k = 0; k < s; k++)
    w[k] = -1; /* clear workspace w [0..s-1] */
  for (k = 0; k < n; k++) /* find first [j] */
  {
    j = post[k];
    delta[j] = (first[j] == -1) ? 1 : 0; /* delta[j]=1 if j is a leaf */
    for (; j != -1 && first[j] == -1; j = parent[j])
      first[j] = k;
  }
  ATp = AT->p;
  ATi = AT->i;
  if (ata)
  {
    head = w + 4 * n;
    next = w + 5 * n + 1;
    for (k = 0; k < n; k++)
      w[post[k]] = k; /* invert post */
    for (i = 0; i < m; i++)
    {
      k = n; /* k = least postordered column in row i */
      for (p = ATp[i]; p < ATp[i + 1]; p++)
        k = CS_MIN(k, w[ATi[p]]);
      next[i] = head[k]; /* place row i in link list k */
      head[k] = i;
    }
  }
  for (i = 0; i < n; i++)
    ancestor[i] = i; /* each node in its own set */
  for (k = 0; k < n; k++)
  {
    j = post[k]; /* j is the kth node in postordered etree */
    if (parent[j] != -1) delta[parent[j]]--; /* j is not a root */
    if (ata)
    {
      for (ii = head[k]; ii != -1; ii = next[ii])
      {
        for (p = ATp[ii]; p < ATp[ii + 1]; p++)
          cs_cedge(j, ATi[p], first, maxfirst, delta, prevleaf, ancestor);
      }
    }
    else
    {
      for (p = ATp[j]; p < ATp[j + 1]; p++)
        cs_cedge(j, ATi[p], first, maxfirst, delta, prevleaf, ancestor);
    }
    if (parent[j] != -1) ancestor[j] = parent[j];
  }
  for (j = 0; j < n; j++) /* sum up delta's of each child */
  {
    if (parent[j] != -1) colcount[parent[j]] += colcount[j];
  }
  return (cs_idone(colcount, AT, w, 1)); /* success: free workspace */
}

/* p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c */
int cs_cumsum(int *p, int *c, int n)
{
  int i, nz = 0;
  if (!p || !c) return (-1); /* check inputs */
  for (i = 0; i < n; i++)
  {
    p[i] = nz;
    nz += c[i];
    c[i] = p[i];
  }
  p[n] = nz;
  return (nz); /* return sum (c [0..n-1]) */
}

/* depth-first-search of the graph of a matrix, starting at node j */
int cs_dfs(int j, cs *L, int top, int *xi, int *pstack, const int *Pinv)
{
  int i, p, p2, done, jnew, head = 0, *Lp, *Li;
  if (!L || !xi || !pstack) return (-1);
  Lp = L->p;
  Li = L->i;
  xi[0] = j; /* initialize the recursion stack */
  while (head >= 0)
  {
    j = xi[head]; /* get j from the top of the recursion stack */
    jnew = Pinv ? (Pinv[j]) : j;
    if (!CS_MARKED(Lp, j))
    {
      CS_MARK(Lp, j); /* mark node j as visited */
      pstack[head] = (jnew < 0) ? 0 : CS_UNFLIP(Lp[jnew]);
    }
    done = 1; /* node j done if no unvisited neighbors */
    p2 = (jnew < 0) ? 0 : CS_UNFLIP(Lp[jnew + 1]);
    for (p = pstack[head]; p < p2; p++) /* examine all neighbors of j */
    {
      i = Li[p]; /* consider neighbor node i */
      if (CS_MARKED(Lp, i)) continue; /* skip visited node i */
      pstack[head] = p; /* pause depth-first search of node j */
      xi[++head] = i; /* start dfs at node i */
      done = 0; /* node j is not done */
      break; /* break, to start dfs (i) */
    }
    if (done) /* depth-first search at node j is done */
    {
      head--; /* remove j from the recursion stack */
      xi[--top] = j; /* and place in the output stack */
    }
  }
  return (top);
}

/* breadth-first search for coarse decomposition (C0,C1,R1 or R0,R3,C3) */
static int cs_bfs(const cs *A,
                  int n,
                  int *wi,
                  int *wj,
                  int *queue,
                  const int *imatch,
                  const int *jmatch,
                  int mark)
{
  int *Ap, *Ai, head = 0, tail = 0, j, i, p, j2;
  cs *C;
  for (j = 0; j < n; j++) /* place all unmatched nodes in queue */
  {
    if (imatch[j] >= 0) continue; /* skip j if matched */
    wj[j] = 0; /* j in set C0 (R0 if transpose) */
    queue[tail++] = j; /* place unmatched col j in queue */
  }
  if (tail == 0) return (1); /* quick return if no unmatched nodes */
  C = (mark == 1) ? ((cs*) A) : cs_transpose(A, 0);
  if (!C) return (0); /* bfs of C=A' to find R0,R3,C3 */
  Ap = C->p;
  Ai = C->i;
  while (head < tail) /* while queue is not empty */
  {
    j = queue[head++]; /* get the head of the queue */
    for (p = Ap[j]; p < Ap[j + 1]; p++)
    {
      i = Ai[p];
      if (wi[i] >= 0) continue; /* skip if i is marked */
      wi[i] = mark; /* i in set R1 (C3 if transpose) */
      j2 = jmatch[i]; /* traverse alternating path to j2 */
      if (wj[j2] >= 0) continue;/* skip j2 if it is marked */
      wj[j2] = mark; /* j2 in set C1 (R3 if transpose) */
      queue[tail++] = j2; /* add j2 to queue */
    }
  }
  if (mark != 1) cs_spfree(C); /* free A' if it was created */
  return (1);
}

/* collect matched rows and columns into P and Q */
static void cs_matched(int m,
                       const int *wi,
                       const int *jmatch,
                       int *P,
                       int *Q,
                       int *cc,
                       int *rr,
                       int set,
                       int mark)
{
  int kc = cc[set], i;
  int kr = rr[set - 1];
  for (i = 0; i < m; i++)
  {
    if (wi[i] != mark) continue; /* skip if i is not in R set */
    P[kr++] = i;
    Q[kc++] = jmatch[i];
  }
  cc[set + 1] = kc;
  rr[set] = kr;
}

/*
 Purpose:

 CS_UNMATCHED collects unmatched rows into the permutation vector P.
 */
static void cs_unmatched(int m, const int *wi, int *P, int *rr, int set)
{
  int i, kr = rr[set];
  for (i = 0; i < m; i++)
    if (wi[i] == 0) P[kr++] = i;
  rr[set + 1] = kr;
}

/* return 1 if row i is in R2 */
static int cs_rprune(int i, int /*j*/, double /*aij*/, void *other)
{
  int *rr = (int*) other;
  return (i >= rr[1] && i < rr[2]);
}

/* Given A, find coarse dmperm */
csd* cs_dmperm(const cs *A)
{
  int m, n, i, j, k, p, cnz, nc, *jmatch, *imatch, *wi, *wj, *Pinv, *Cp, *Ci,
      *Ps, *Rs, nb1, nb2, *P, *Q, *cc, *rr, *R, *S, ok;
  cs *C;
  csd *D, *scc;
  /* --- Maximum matching ------------------------------------------------- */
  if (!A) return (NULL); /* check inputs */
  m = A->m;
  n = A->n;
  D = cs_dalloc(m, n); /* allocate result */
  if (!D) return (NULL);
  P = D->P;
  Q = D->Q;
  R = D->R;
  S = D->S;
  cc = D->cc;
  rr = D->rr;
  jmatch = cs_maxtrans(A); /* max transversal */
  imatch = jmatch + m; /* imatch = inverse of jmatch */
  if (!jmatch) return (cs_ddone(D, NULL, jmatch, 0));
  /* --- Coarse decomposition --------------------------------------------- */
  wi = R;
  wj = S; /* use R and S as workspace */
  for (j = 0; j < n; j++)
    wj[j] = -1; /* unmark all cols for bfs */
  for (i = 0; i < m; i++)
    wi[i] = -1; /* unmark all rows for bfs */
  cs_bfs(A, n, wi, wj, Q, imatch, jmatch, 1); /* find C0, C1, R1 */
  ok = cs_bfs(A, m, wj, wi, P, jmatch, imatch, 3); /* find R0, R3, C3 */
  if (!ok) return (cs_ddone(D, NULL, jmatch, 0));
  cs_unmatched(n, wj, Q, cc, 0); /* unmatched set C0 */
  cs_matched(m, wi, jmatch, P, Q, cc, rr, 1, 1); /* set R1 and C1 */
  cs_matched(m, wi, jmatch, P, Q, cc, rr, 2, -1); /* set R2 and C2 */
  cs_matched(m, wi, jmatch, P, Q, cc, rr, 3, 3); /* set R3 and C3 */
  cs_unmatched(m, wi, P, rr, 3); /* unmatched set R0 */
  cs_free(jmatch);
  /* --- Fine decomposition ----------------------------------------------- */
  Pinv = cs_pinv(P, m); /* Pinv=P' */
  if (!Pinv) return (cs_ddone(D, NULL, NULL, 0));
  C = cs_permute(A, Pinv, Q, 0);/* C=A(P,Q) (it will hold A(R2,C2)) */
  cs_free(Pinv);
  if (!C) return (cs_ddone(D, NULL, NULL, 0));
  Cp = C->p;
  Ci = C->i;
  nc = cc[3] - cc[2]; /* delete cols C0, C1, and C3 from C */
  if (cc[2] > 0) for (j = cc[2]; j <= cc[3]; j++)
    Cp[j - cc[2]] = Cp[j];
  C->n = nc;
  if (rr[2] - rr[1] < m) /* delete rows R0, R1, and R3 from C */
  {
    cs_fkeep(C, cs_rprune, rr);
    cnz = Cp[nc];
    if (rr[1] > 0) for (p = 0; p < cnz; p++)
      Ci[p] -= rr[1];
  }
  C->m = nc;
  scc = cs_scc(C); /* find strongly-connected components of C*/
  if (!scc) return (cs_ddone(D, C, NULL, 0));
  /* --- Combine coarse and fine decompositions --------------------------- */
  Ps = scc->P; /* C(Ps,Ps) is the permuted matrix */
  Rs = scc->R; /* kth block is Rs[k]..Rs[k+1]-1 */
  nb1 = scc->nb; /* # of blocks of A(*/
  for (k = 0; k < nc; k++)
    wj[k] = Q[Ps[k] + cc[2]]; /* combine */
  for (k = 0; k < nc; k++)
    Q[k + cc[2]] = wj[k];
  for (k = 0; k < nc; k++)
    wi[k] = P[Ps[k] + rr[1]];
  for (k = 0; k < nc; k++)
    P[k + rr[1]] = wi[k];
  nb2 = 0; /* create the fine block partitions */
  R[0] = 0;
  S[0] = 0;
  if (cc[2] > 0) nb2++; /* leading coarse block A (R1, [C0 C1]) */
  for (k = 0; k < nb1; k++) /* coarse block A (R2,C2) */
  {
    R[nb2] = Rs[k] + rr[1]; /* A (R2,C2) splits into nb1 fine blocks */
    S[nb2] = Rs[k] + cc[2];
    nb2++;
  }
  if (rr[2] < m)
  {
    R[nb2] = rr[2]; /* trailing coarse block A ([R3 R0], C3) */
    S[nb2] = cc[3];
    nb2++;
  }
  R[nb2] = m;
  S[nb2] = n;
  D->nb = nb2;
  cs_dfree(scc);
  return (cs_ddone(D, C, NULL, 1));
}

static int cs_tol(int /*i*/, int /*j*/, double aij, void *tol)
{
  return (fabs(aij) > *((double*) tol));
}
int cs_droptol(cs *A, double tol)
{
  return (cs_fkeep(A, &cs_tol, &tol)); /* keep all large entries */
}

static int cs_nonzero(int /*i*/, int /*j*/, double aij, void* /*other*/)
{
  return (aij != 0);
}
int cs_dropzeros(cs *A)
{
  return (cs_fkeep(A, &cs_nonzero, NULL)); /* keep all nonzero entries */
}

/*
 Purpose:

 CS_DUPL removes duplicate entries from A.

 Reference:

 Timothy Davis,
 Direct Methods for Sparse Linear Systems,
 SIAM, Philadelphia, 2006.
 */
int cs_dupl(cs *A)
{
  int i, j, p, q, nz = 0, n, m, *Ap, *Ai, *w;
  double *Ax;
  if (!A) return (0); /* check inputs */
  m = A->m;
  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  w = (int*) cs_malloc(m, sizeof(int)); /* get workspace */
  if (!w) return (0); /* out of memory */
  for (i = 0; i < m; i++)
    w[i] = -1; /* row i not yet seen */
  for (j = 0; j < n; j++)
  {
    q = nz; /* column j will start at q */
    for (p = Ap[j]; p < Ap[j + 1]; p++)
    {
      i = Ai[p]; /* A(i,j) is nonzero */
      if (w[i] >= q)
      {
        Ax[w[i]] += Ax[p]; /* A(i,j) is a duplicate */
      }
      else
      {
        w[i] = nz; /* record where row i occurs */
        Ai[nz] = i; /* keep A(i,j) */
        Ax[nz++] = Ax[p];
      }
    }
    Ap[j] = q; /* record start of column j */
  }
  Ap[n] = nz; /* finalize A */
  cs_free(w); /* free workspace */
  return (cs_sprealloc(A, 0)); /* remove extra space from A */
}

/* add an entry to a triplet matrix; return 1 if ok, 0 otherwise */
int cs_entry(cs *T, int i, int j, double x)
{
  if (!T || (T->nz >= T->nzmax && !cs_sprealloc(T, 2 * (T->nzmax)))) return (0);
  if (T->x) T->x[T->nz] = x;
  T->i[T->nz] = i;
  T->p[T->nz++] = j;
  T->m = CS_MAX(T->m, i + 1);
  T->n = CS_MAX(T->n, j + 1);
  return (1);
}

/* compute the etree of A (using triu(A), or A'A without forming A'A */
int* cs_etree(const cs *A, int ata)
{
  int i, k, p, m, n, inext, *Ap, *Ai, *w, *parent, *ancestor, *prev;
  if (!A) return (NULL); /* check inputs */
  m = A->m;
  n = A->n;
  Ap = A->p;
  Ai = A->i;
  parent = (int*) cs_malloc(n, sizeof(int));
  w = (int*) cs_malloc(n + (ata ? m : 0), sizeof(int));
  ancestor = w;
  prev = w + n;
  if (!w || !parent) return (cs_idone(parent, NULL, w, 0));
  if (ata) for (i = 0; i < m; i++)
    prev[i] = -1;
  for (k = 0; k < n; k++)
  {
    parent[k] = -1; /* node k has no parent yet */
    ancestor[k] = -1; /* nor does k have an ancestor */
    for (p = Ap[k]; p < Ap[k + 1]; p++)
    {
      i = ata ? (prev[Ai[p]]) : (Ai[p]);
      for (; i != -1 && i < k; i = inext) /* traverse from i to k */
      {
        inext = ancestor[i]; /* inext = ancestor of i */
        ancestor[i] = k; /* path compression */
        if (inext == -1) parent[i] = k; /* no anc., parent is k */
      }
      if (ata) prev[Ai[p]] = k;
    }
  }
  return (cs_idone(parent, NULL, w, 1));
}

/* drop entries for which fkeep(A(i,j)) is false; return nz if OK, else -1 */
int cs_fkeep(cs *A, int (*fkeep)(int, int, double, void*), void *other)
{
  int j, p, nz = 0, n, *Ap, *Ai;
  double *Ax;
  if (!A || !fkeep) return (-1); /* check inputs */
  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  for (j = 0; j < n; j++)
  {
    p = Ap[j]; /* get current location of col j */
    Ap[j] = nz; /* record new location of col j */
    for (; p < Ap[j + 1]; p++)
    {
      if (fkeep(Ai[p], j, Ax ? Ax[p] : 1, other))
      {
        if (Ax) Ax[nz] = Ax[p]; /* keep A(i,j) */
        Ai[nz++] = Ai[p];
      }
    }
  }
  return (Ap[n] = nz); /* finalize A and return nnz(A) */
}

/* y = A*x+y */
int cs_gaxpy(const cs *A, const double *x, double *y)
{
  int p, j, n, *Ap, *Ai;
  double *Ax;
  if (!A || !x || !y) return (0); /* check inputs */
  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  for (j = 0; j < n; j++)
  {
    for (p = Ap[j]; p < Ap[j + 1]; p++)
    {
      y[Ai[p]] += Ax[p] * x[j];
    }
  }
  return (1);
}

/* apply the ith Householder vector to x */
int cs_happly(const cs *V, int i, double beta, double *x)
{
  int p, *Vp, *Vi;
  double *Vx, tau = 0;
  if (!V || !x) return (0); /* check inputs */
  Vp = V->p;
  Vi = V->i;
  Vx = V->x;
  for (p = Vp[i]; p < Vp[i + 1]; p++) /* tau = v'*x */
  {
    tau += Vx[p] * x[Vi[p]];
  }
  tau *= beta; /* tau = beta*(v'*x) */
  for (p = Vp[i]; p < Vp[i + 1]; p++) /* x = x - v*tau */
  {
    x[Vi[p]] -= Vx[p] * tau;
  }
  return (1);
}

/* create a Householder reflection [v,beta,s]=house(x), overwrite x with v,
 * where (I-beta*v*v')*x = s*x.  See Algo 5.1.1, Golub & Van Loan, 3rd ed. */
double cs_house(double *x, double *beta, int n)
{
  double s, sigma = 0;
  int i;
  if (!x || !beta) return (-1); /* check inputs */
  for (i = 1; i < n; i++)
    sigma += x[i] * x[i];
  if (sigma == 0)
  {
    s = fabs(x[0]); /* s = |x(0)| */
    (*beta) = (x[0] <= 0) ? 2 : 0;
    x[0] = 1;
  }
  else
  {
    s = sqrt(x[0] * x[0] + sigma); /* s = norm (x) */
    x[0] = (x[0] <= 0) ? (x[0] - s) : (-sigma / (x[0] + s));
    (*beta) = -1. / (s * x[0]);
  }
  return (s);
}

/* x(P) = b, for dense vectors x and b; P=NULL denotes identity */
int cs_ipvec(int n, const int *P, const double *b, double *x)
{
  int k;
  if (!x || !b) return (0); /* check inputs */
  for (k = 0; k < n; k++)
    x[P ? P[k] : k] = b[k];
  return (1);
}

/*
 Purpose:

 CS_LOAD loads a triplet matrix from a file.

 Reference:

 Timothy Davis,
 Direct Methods for Sparse Linear Systems,
 SIAM, Philadelphia, 2006.
 */
cs* cs_load(FILE *f)
{
  int i, j;
  double x;
  cs *T;
  if (!f) return (NULL);
  T = cs_spalloc(0, 0, 1, 1, 1);
  while (gslFScanf(f, "%d %d %lg\n", &i, &j, &x) == 3)
  {
    if (!cs_entry(T, i, j, x)) return (cs_spfree(T));
  }
  return (T);
}

/*
 Purpose:

 CS_LSOLVE solves L*x=b.

 Discussion:

 On input, X contains the right hand side, and on output, the solution.

 Reference:

 Timothy Davis,
 Direct Methods for Sparse Linear Systems,
 SIAM, Philadelphia, 2006.
 */
int cs_lsolve(const cs *L, double *x)
{
  int p, j, n, *Lp, *Li, Lpj, Lpjp1;
  double *Lx, xval;
  if (!L || !x) return (0); /* check inputs */
  n = L->n;
  Lp = L->p;
  Li = L->i;
  Lx = L->x;

  Lpj = Lpjp1 = Lp[0];
  for (j = 0; j < n; j++)
  {
    Lpj = Lpjp1;
    Lpjp1 = Lp[j + 1];
    x[j] /= Lx[Lpj];
    xval = x[j];
    for (p = Lpj + 1; p < Lpjp1; p++)
    {
      x[Li[p]] -= Lx[p] * xval;
    }
  }
  return (1);
}

/*
 Purpose:

 CS_LTSOLVE solves L'*x=b.

 Discussion:

 On input, X contains the right hand side, and on output, the solution.

 Reference:

 Timothy Davis,
 Direct Methods for Sparse Linear Systems,
 SIAM, Philadelphia, 2006.
 */
int cs_ltsolve(const cs *L, double *x)
{
  int p, j, n, *Lp, *Li, Lpj, Lpjp1;
  double *Lx;
  if (!L || !x) return (0); /* check inputs */
  n = L->n;
  Lp = L->p;
  Li = L->i;
  Lx = L->x;

  Lpj = Lpjp1 = Lp[n];
  for (j = n - 1; j >= 0; j--)
  {
    Lpjp1 = Lpj;
    Lpj = Lp[j];
    for (p = Lpj + 1; p < Lpjp1; p++)
    {
      x[j] -= Lx[p] * x[Li[p]];
    }
    x[j] /= Lx[Lpj];
  }
  return (1);
}

/* [L,U,Pinv]=lu(A, [Q lnz unz]). lnz and unz can be guess */
csn* cs_lu(const cs *A, const css *S, double tol)
{
  cs *L, *U;
  csn *N;
  double pivot, *Lx, *Ux, *x, a, t;
  int *Lp, *Li, *Up, *Ui, *Pinv, *xi, *Q, n, ipiv, k, top, p, i, col, lnz, unz;
  if (!A || !S) return (NULL); /* check inputs */
  n = A->n;
  Q = S->Q;
  lnz = S->lnz;
  unz = S->unz;
  x = (double*) cs_malloc(n, sizeof(double));
  xi = (int*) cs_malloc(2 * n, sizeof(int));
  N = (csn*) cs_calloc(1, sizeof(csn));
  if (!x || !xi || !N) return (cs_ndone(N, NULL, xi, x, 0));
  N->L = L = cs_spalloc(n, n, lnz, 1, 0); /* initial L and U */
  N->U = U = cs_spalloc(n, n, unz, 1, 0);
  N->Pinv = Pinv = (int*) cs_malloc(n, sizeof(int));
  if (!L || !U || !Pinv)
  {
    messerr("Core allocation problem in CSparse Library (%d x %d)", n, n);
    return (cs_ndone(N, NULL, xi, x, 0));
  }
  Lp = L->p;
  Up = U->p;
  for (i = 0; i < n; i++)
    x[i] = 0; /* clear workspace */
  for (i = 0; i < n; i++)
    Pinv[i] = -1; /* no rows pivotal yet */
  for (k = 0; k <= n; k++)
    Lp[k] = 0; /* no cols of L yet */
  lnz = unz = 0;
  for (k = 0; k < n; k++) /* compute L(:,k) and U(:,k) */
  {
    /* --- Triangular solve --------------------------------------------- */
    Lp[k] = lnz; /* L(:,k) starts here */
    Up[k] = unz; /* U(:,k) starts here */
    if ((lnz + n > L->nzmax && !cs_sprealloc(L, 2 * L->nzmax + n)) || (unz + n
        > U->nzmax
                                                                       && !cs_sprealloc(
                                                                           U,
                                                                           2 * U->nzmax + n)))
    {
      return (cs_ndone(N, NULL, xi, x, 0));
    }
    Li = L->i;
    Lx = L->x;
    Ui = U->i;
    Ux = U->x;
    col = Q ? (Q[k]) : k;
    top = cs_splsolve(L, A, col, xi, x, Pinv); /* x = L\A(:,col) */
    /* --- Find pivot --------------------------------------------------- */
    ipiv = -1;
    a = -1;
    for (p = top; p < n; p++)
    {
      i = xi[p]; /* x(i) is nonzero */
      if (Pinv[i] < 0) /* row i is not pivotal */
      {
        if ((t = fabs(x[i])) > a)
        {
          a = t; /* largest pivot candidate so far */
          ipiv = i;
        }
      }
      else /* x(i) is the entry U(Pinv[i],k) */
      {
        Ui[unz] = Pinv[i];
        Ux[unz++] = x[i];
      }
    }
    if (ipiv == -1 || a <= 0) return (cs_ndone(N, NULL, xi, x, 0));
    if (Pinv[col] < 0 && fabs(x[col]) >= a * tol) ipiv = col;
    /* --- Divide by pivot ---------------------------------------------- */
    pivot = x[ipiv]; /* the chosen pivot */
    Ui[unz] = k; /* last entry in U(:,k) is U(k,k) */
    Ux[unz++] = pivot;
    Pinv[ipiv] = k; /* ipiv is the kth pivot row */
    Li[lnz] = ipiv; /* first entry in L(:,k) is L(k,k) = 1 */
    Lx[lnz++] = 1;
    for (p = top; p < n; p++) /* L(k+1:n,k) = x / pivot */
    {
      i = xi[p];
      if (Pinv[i] < 0) /* x(i) is an entry in L(:,k) */
      {
        Li[lnz] = i; /* save unpermuted row in L */
        Lx[lnz++] = x[i] / pivot; /* scale pivot column */
      }
      x[i] = 0; /* x [0..n-1] = 0 for next k */
    }
  }
  /* --- Finalize L and U ------------------------------------------------- */
  Lp[n] = lnz;
  Up[n] = unz;
  Li = L->i; /* fix row indices of L for final Pinv */
  for (p = 0; p < lnz; p++)
    Li[p] = Pinv[Li[p]];
  cs_sprealloc(L, 0); /* remove extra space from L and U */
  cs_sprealloc(U, 0);
  return (cs_ndone(N, NULL, xi, x, 1)); /* success */
}

/* x=A\b where A is unsymmetric; b overwritten with solution */
int cs_lusol(const cs *A, double *b, int order, double tol)
{
  double *x;
  css *S;
  csn *N;
  int n, ok;
  if (!A || !b) return (0); /* check inputs */
  n = A->n;
  S = cs_sqr(A, order, 0); /* ordering and symbolic analysis */
  N = cs_lu(A, S, tol); /* numeric LU factorization */
  x = (double*) cs_malloc(n, sizeof(double));
  ok = (S && N && x);
  if (ok)
  {
    cs_ipvec(n, N->Pinv, b, x); /* x = P*b */
    cs_lsolve(N->L, x); /* x = L\x */
    cs_usolve(N->U, x); /* x = U\x */
    cs_ipvec(n, S->Q, x, b); /* b = Q*x */
  }
  cs_free(x);
  cs_sfree(S);
  cs_nfree(N);
  return (ok);
}

#ifdef MATLAB_MEX_FILE
  #define malloc mxMalloc
  #define free mxFree
  #define realloc mxRealloc
  #define calloc mxCalloc
  #endif

/* wrapper for malloc */
void* cs_malloc(int n, size_t size)
{
  return (CS_OVERFLOW(n, size) ?
  NULL : mem_alloc(CS_MAX(n, 1) * static_cast<int>(size), 0));
}

/* wrapper for calloc */
void* cs_calloc(int n, size_t size)
{
  return (CS_OVERFLOW(n, size) ?
  NULL : (void*) mem_calloc(CS_MAX(n, 1), static_cast<int>(size), 0));
}

/* wrapper for free */
void* cs_free(void *p)
{
  if (p)
  {
    mem_free ((char *) p); /* free p if it is not already NULL */
  }
  return (NULL); /* return NULL to simplify the use of cs_free */
}

/* wrapper for realloc */
void* cs_realloc(void *p, int n, size_t size, int *ok)
{
  void *p2;
  *ok = !CS_OVERFLOW(n, size); /* guard against int overflow */
  if (!(*ok)) return (p); /* p unchanged if n too large */
  /* realloc the block */
  p2 = (void*) mem_realloc((char* ) p, CS_MAX(n, 1) * static_cast<int>(size),
                           0);
  *ok = (p2 != NULL);
  return ((*ok) ? p2 : p); /* return original p if failure */
}

/* find an augmenting path starting at column k and extend the match if found */
static void cs_augment(int k,
                       const cs *A,
                       int *jmatch,
                       int *cheap,
                       int *w,
                       int *js,
                       int *is,
                       int *ps)
{
  int found = 0, p, i = -1, *Ap = A->p, *Ai = A->i, head = 0, j;
  js[0] = k; /* start with just node k in jstack */
  while (head >= 0)
  {
    /* --- Start (or continue) depth-first-search at node j ------------- */
    j = js[head]; /* get j from top of jstack */
    if (w[j] != k) /* 1st time j visited for kth path */
    {
      w[j] = k; /* mark j as visited for kth path */
      for (p = cheap[j]; p < Ap[j + 1] && !found; p++)
      {
        i = Ai[p]; /* try a cheap assignment (i,j) */
        found = (jmatch[i] == -1);
      }
      cheap[j] = p; /* start here next time j is traversed*/
      if (found)
      {
        is[head] = i; /* column j matched with row i */
        break; /* end of augmenting path */
      }
      ps[head] = Ap[j]; /* no cheap match: start dfs for j */
    }
    /* --- Depth-first-search of neighbors of j ------------------------- */
    for (p = ps[head]; p < Ap[j + 1]; p++)
    {
      i = Ai[p]; /* consider row i */
      if (w[jmatch[i]] == k) continue; /* skip jmatch [i] if marked */
      ps[head] = p + 1; /* pause dfs of node j */
      is[head] = i; /* i will be matched with j if found */
      js[++head] = jmatch[i]; /* start dfs at column jmatch [i] */
      break;
    }
    if (p == Ap[j + 1]) head--; /* node j is done; pop from stack */
  } /* augment the match if path found: */
  if (found) for (p = head; p >= 0; p--)
    jmatch[is[p]] = js[p];
}

/* find a maximum transveral */
int* cs_maxtrans(const cs *A) /* returns jmatch [0..m-1]; imatch [0..n-1] */
{
  int i, j, k, n, m, p, n2 = 0, m2 = 0, *Ap, *jimatch, *w, *cheap, *js, *is,
      *ps, *Ai, *Cp, *jmatch, *imatch;
  cs *C;
  if (!A) return (NULL); /* check inputs */
  n = A->n;
  m = A->m;
  Ap = A->p;
  Ai = A->i;
  w = jimatch = (int*) cs_calloc(m + n, sizeof(int)); /* allocate result */
  if (!jimatch) return (NULL);
  for (j = 0; j < n; j++) /* count non-empty rows and columns */
  {
    n2 += (Ap[j] < Ap[j + 1]);
    for (p = Ap[j]; p < Ap[j + 1]; p++)
      w[Ai[p]] = 1;
  }
  for (i = 0; i < m; i++)
    m2 += w[i];
  C = (m2 < n2) ? cs_transpose(A, 0) : ((cs*) A); /* transpose if needed */
  if (!C) return (cs_idone(jimatch, (m2 < n2) ? C : NULL, NULL, 0));
  n = C->n;
  m = C->m;
  Cp = C->p;
  jmatch = (m2 < n2) ? jimatch + n : jimatch;
  imatch = (m2 < n2) ? jimatch : jimatch + m;
  w = (int*) cs_malloc(5 * n, sizeof(int)); /* allocate workspace */
  if (!w) return (cs_idone(jimatch, (m2 < n2) ? C : NULL, w, 0));
  cheap = w + n;
  js = w + 2 * n;
  is = w + 3 * n;
  ps = w + 4 * n;
  for (j = 0; j < n; j++)
    cheap[j] = Cp[j]; /* for cheap assignment */
  for (j = 0; j < n; j++)
    w[j] = -1; /* all columns unflagged */
  for (i = 0; i < m; i++)
    jmatch[i] = -1; /* nothing matched yet */
  for (k = 0; k < n; k++)
    cs_augment(k, C, jmatch, cheap, w, js, is, ps);
  for (j = 0; j < n; j++)
    imatch[j] = -1; /* find row match */
  for (i = 0; i < m; i++)
    if (jmatch[i] >= 0) imatch[jmatch[i]] = i;
  return (cs_idone(jimatch, (m2 < n2) ? C : NULL, w, 1));
}

/* C = A * B */
cs* cs_multiply(const cs *A, const cs *B)
{
  int p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values, *Bi;
  double *x, *Bx, *Cx;
  cs *C;
  if (!A || !B)
  {
    messerr("Both matrices should be defined");
    return (NULL); /* check inputs */
  }
  m = A->m;
  anz = A->p[A->n];
  n = B->n;
  Bp = B->p;
  Bi = B->i;
  Bx = B->x;
  bnz = Bp[n];
  w = (int*) cs_calloc(m, sizeof(int));
  values = (A->x != NULL) && (Bx != NULL);
  x = values ? (double*) cs_malloc(m, sizeof(double)) : NULL;
  C = cs_spalloc(m, n, anz + bnz, values, 0);
  if (!C || !w || (values && !x))
  {
    messerr("Core allocation problem in CSparse Library (%d x %d)", m, n);
    return (cs_done(C, w, x, 0));
  }
  Cp = C->p;
  for (j = 0; j < n; j++)
  {
    if (nz + m > C->nzmax && !cs_sprealloc(C, 2 * (C->nzmax) + m))
    {
      messerr("Out of memory in CSparse Library %d", 2 * (C->nzmax) + m);
      return (cs_done(C, w, x, 0)); /* out of memory */
    }
    Ci = C->i;
    Cx = C->x; /* C may have been reallocated */
    Cp[j] = nz; /* column j of C starts here */
    for (p = Bp[j]; p < Bp[j + 1]; p++)
    {
      nz = cs_scatter(A, Bi[p], Bx ? Bx[p] : 1, w, x, j + 1, C, nz);
    }
    if (values) for (p = Cp[j]; p < nz; p++)
      Cx[p] = x[Ci[p]];
  }
  Cp[n] = nz; /* finalize the last column of C */
  cs_sprealloc(C, 0); /* remove extra space from C */
  return (cs_done(C, w, x, 1)); /* success; free workspace, return C */
}

/* 1-norm of a sparse matrix = max (sum (abs (A))), largest column sum */
double cs_norm(const cs *A)
{
  int p, j, n, *Ap;
  double *Ax, norm = 0, s;
  if (!A || !A->x) return (-1); /* check inputs */
  n = A->n;
  Ap = A->p;
  Ax = A->x;
  for (j = 0; j < n; j++)
  {
    for (s = 0, p = Ap[j]; p < Ap[j + 1]; p++)
      s += fabs(Ax[p]);
    norm = CS_MAX(norm, s);
  }
  return (norm);
}

/* C = A(P,Q) where P and Q are permutations of 0..m-1 and 0..n-1. */
cs* cs_permute(const cs *A, const int *Pinv, const int *Q, int values)
{
  int p, j, k, nz = 0, m, n, *Ap, *Ai, *Cp, *Ci;
  double *Cx, *Ax;
  cs *C;
  if (!A) return (NULL); /* check inputs */
  m = A->m;
  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  C = cs_spalloc(m, n, Ap[n], values && Ax != NULL, 0);
  if (!C)
  {
    messerr("Core allocation problem in CSparse Library (%d x %d)", m, n);
    return (cs_done(C, NULL, NULL, 0)); /* out of memory */
  }
  Cp = C->p;
  Ci = C->i;
  Cx = C->x;
  for (k = 0; k < n; k++)
  {
    Cp[k] = nz; /* column k of C is column Q[k] of A */
    j = Q ? (Q[k]) : k;
    for (p = Ap[j]; p < Ap[j + 1]; p++)
    {
      if (Cx) Cx[nz] = Ax[p]; /* row i of A is row Pinv[i] of C */
      Ci[nz++] = Pinv ? (Pinv[Ai[p]]) : Ai[p];
    }
  }
  Cp[n] = nz; /* finalize the last column of C */
  return (cs_done(C, NULL, NULL, 1));
}

/* Pinv = P', or P = Pinv' */
int* cs_pinv(int const *P, int n)
{
  int k, *Pinv;
  if (!P) return (NULL); /* P = NULL denotes identity */
  Pinv = (int*) cs_malloc(n, sizeof(int)); /* allocate resuult */
  if (!Pinv) return (NULL); /* out of memory */
  for (k = 0; k < n; k++)
    Pinv[P[k]] = k;/* invert the permutation */
  return (Pinv); /* return result */
}

/* post order a forest */
int* cs_post(int n, const int *parent)
{
  int j, k = 0, *post, *w, *head, *next, *stack;
  if (!parent) return (NULL); /* check inputs */
  post = (int*) cs_malloc(n, sizeof(int)); /* allocate result */
  w = (int*) cs_malloc(3 * n, sizeof(int)); /* 3*n workspace */
  head = w;
  next = w + n;
  stack = w + 2 * n;
  if (!w || !post) return (cs_idone(post, NULL, w, 0));
  for (j = 0; j < n; j++)
    head[j] = -1; /* empty link lists */
  for (j = n - 1; j >= 0; j--) /* traverse nodes in reverse order*/
  {
    if (parent[j] == -1) continue; /* j is a root */
    next[j] = head[parent[j]]; /* add j to list of its parent */
    head[parent[j]] = j;
  }
  for (j = 0; j < n; j++)
  {
    if (parent[j] != -1) continue; /* skip j if it is not a root */
    k = cs_tdfs(j, k, head, next, post, stack);
  }
  return (cs_idone(post, NULL, w, 1)); /* success; free w, return post */
}

/* print a sparse matrix */
int cs_print(const cs *A, int brief)
{
  int p, j, m, n, nzmax, nz, *Ap, *Ai;
  double *Ax;
  if (!A)
  {
    message("(null)\n");
    return (0);
  }
  m = A->m;
  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  nzmax = A->nzmax;
  nz = A->nz;
  message("CSparse Version %d.%d.%d, %s.  %s\n", CS_VER, CS_SUBVER, CS_SUBSUB,
          CS_DATE, CS_COPYRIGHT);
  if (nz < 0)
  {
    message("%d-by-%d, nzmax: %d nnz: %d, 1-norm: %g\n", m, n, nzmax, Ap[n],
            cs_norm(A));
    for (j = 0; j < n; j++)
    {
      message("    col %d : locations %d to %d\n", j, Ap[j], Ap[j + 1] - 1);
      for (p = Ap[j]; p < Ap[j + 1]; p++)
      {
        message("     %d : %g\n", Ai[p], Ax ? Ax[p] : 1);
        if (brief && p > 20)
        {
          message("  ...\n");
          return (1);
        }
      }
    }
  }
  else
  {
    message("triplet: %d-by-%d, nzmax: %d nnz: %d\n", m, n, nzmax, nz);
    for (p = 0; p < nz; p++)
    {
      message("     %d %d : %g\n", Ai[p], Ap[p], Ax ? Ax[p] : 1);
      if (brief && p > 20)
      {
        message("  ...\n");
        return (1);
      }
    }
  }
  return (1);
}

/* x = b(P), for dense vectors x and b; P=NULL denotes identity */
int cs_pvec(int n, const int *P, const double *b, double *x)
{
  int k;
  if (!x || !b) return (0); /* check inputs */
  for (k = 0; k < n; k++)
    x[k] = b[P ? P[k] : k];
  return (1);
}

/* sparse QR factorization [V,beta,p,R] = qr (A) */
csn* cs_qr(const cs *A, const css *S)
{
  double *Rx, *Vx, *Ax, *Beta, *x;
  int i, k, p, m, n, vnz, p1, top, m2, len, col, rnz, *s, *leftmost, *Ap, *Ai,
      *parent, *Rp, *Ri, *Vp, *Vi, *w, *Pinv, *Q;
  cs *R, *V;
  csn *N;
  if (!A || !S || !S->parent || !S->Pinv) return (NULL); /* check inputs */
  m = A->m;
  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  Q = S->Q;
  parent = S->parent;
  Pinv = S->Pinv;
  m2 = S->m2;
  vnz = S->lnz;
  rnz = S->unz;
  leftmost = Pinv + m + n;
  w = (int*) cs_malloc(m2 + n, sizeof(int));
  x = (double*) cs_malloc(m2, sizeof(double));
  N = (csn*) cs_calloc(1, sizeof(csn));
  if (!w || !x || !N) return (cs_ndone(N, NULL, w, x, 0));
  s = w + m2; /* size n */
  for (k = 0; k < m2; k++)
    x[k] = 0; /* clear workspace x */
  N->L = V = cs_spalloc(m2, n, vnz, 1, 0); /* allocate V */
  N->U = R = cs_spalloc(m2, n, rnz, 1, 0); /* allocate R, m2-by-n */
  N->B = Beta = (double*) cs_malloc(n, sizeof(double));
  if (!R || !V || !Beta)
  {
    messerr("Core allocation problem in CSparse Library (%d x %d)", m2, n);
    return (cs_ndone(N, NULL, w, x, 0));
  }
  Rp = R->p;
  Ri = R->i;
  Rx = R->x;
  Vp = V->p;
  Vi = V->i;
  Vx = V->x;
  for (i = 0; i < m2; i++)
    w[i] = -1; /* clear w, to mark nodes */
  rnz = 0;
  vnz = 0;
  for (k = 0; k < n; k++) /* compute V and R */
  {
    Rp[k] = rnz; /* R(:,k) starts here */
    Vp[k] = p1 = vnz; /* V(:,k) starts here */
    w[k] = k; /* add V(k,k) to pattern of V */
    Vi[vnz++] = k;
    top = n;
    col = Q ? Q[k] : k;
    for (p = Ap[col]; p < Ap[col + 1]; p++) /* find R(:,k) pattern */
    {
      i = leftmost[Ai[p]]; /* i = min(find(A(i,Q))) */
      for (len = 0; w[i] != k; i = parent[i]) /* traverse up to k */
      {
        s[len++] = i;
        w[i] = k;
      }
      while (len > 0)
        s[--top] = s[--len]; /* push path on stack */
      i = Pinv[Ai[p]]; /* i = permuted row of A(:,col) */
      x[i] = Ax[p]; /* x (i) = A(.,col) */
      if (i > k && w[i] < k) /* pattern of V(:,k) = x (k+1:m) */
      {
        Vi[vnz++] = i; /* add i to pattern of V(:,k) */
        w[i] = k;
      }
    }
    for (p = top; p < n; p++) /* for each i in pattern of R(:,k) */
    {
      i = s[p]; /* R(i,k) is nonzero */
      cs_happly(V, i, Beta[i], x); /* apply (V(i),Beta(i)) to x */
      Ri[rnz] = i; /* R(i,k) = x(i) */
      Rx[rnz++] = x[i];
      x[i] = 0;
      if (parent[i] == k) vnz = cs_scatter(V, i, 0, w, NULL, k, V, vnz);
    }
    for (p = p1; p < vnz; p++) /* gather V(:,k) = x */
    {
      Vx[p] = x[Vi[p]];
      x[Vi[p]] = 0;
    }
    Ri[rnz] = k; /* R(k,k) = norm (x) */
    Rx[rnz++] = cs_house(Vx + p1, Beta + k, vnz - p1); /* [v,beta]=house(x) */
  }
  Rp[n] = rnz; /* finalize R */
  Vp[n] = vnz; /* finalize V */
  return (cs_ndone(N, NULL, w, x, 1)); /* success */
}

/* x=A\b where A can be rectangular; b overwritten with solution */
int cs_qrsol(const cs *A, double *b, int order)
{
  double *x;
  css *S;
  csn *N;
  cs *AT = NULL;
  int k, m, n, ok;
  if (!A || !b) return (0); /* check inputs */
  n = A->n;
  m = A->m;
  if (m >= n)
  {
    S = cs_sqr(A, order, 1); /* ordering and symbolic analysis */
    N = cs_qr(A, S); /* numeric QR factorization */
    x = (double*) cs_calloc(S ? S->m2 : 1, sizeof(double));
    ok = (S && N && x);
    if (ok)
    {
      cs_ipvec(m, S->Pinv, b, x); /* x(0:m-1) = P*b(0:m-1) */
      for (k = 0; k < n; k++) /* apply Householder refl. to x */
      {
        cs_happly(N->L, k, N->B[k], x);
      }
      cs_usolve(N->U, x); /* x = R\x */
      cs_ipvec(n, S->Q, x, b); /* b(0:n-1) = Q*x (permutation) */
    }
  }
  else
  {
    AT = cs_transpose(A, 1); /* Ax=b is underdetermined */
    S = cs_sqr(AT, order, 1); /* ordering and symbolic analysis */
    N = cs_qr(AT, S); /* numeric QR factorization of A' */
    x = (double*) cs_calloc(S ? S->m2 : 1, sizeof(double));
    ok = (AT && S && N && x);
    if (ok)
    {
      cs_pvec(m, S->Q, b, x); /* x(0:m-1) = Q'*b (permutation) */
      cs_utsolve(N->U, x); /* x = R'\x */
      for (k = m - 1; k >= 0; k--) /* apply Householder refl. to x */
      {
        cs_happly(N->L, k, N->B[k], x);
      }
      cs_pvec(n, S->Pinv, x, b); /* b (0:n-1) = P'*x */
    }
  }
  cs_free(x);
  cs_sfree(S);
  cs_nfree(N);
  cs_spfree(AT);
  return (ok);
}

/* xi [top...n-1] = nodes reachable from graph of L*P' via nodes in B(:,k).
 * xi [n...2n-1] used as workspace */
int cs_reach(cs *L, const cs *B, int k, int *xi, const int *Pinv)
{
  int p, n, top, *Bp, *Bi, *Lp;
  if (!L || !B || !xi) return (-1);
  n = L->n;
  Bp = B->p;
  Bi = B->i;
  Lp = L->p;
  top = n;
  for (p = Bp[k]; p < Bp[k + 1]; p++)
  {
    if (!CS_MARKED(Lp, Bi[p])) /* start a dfs at unmarked node i */
    {
      top = cs_dfs(Bi[p], L, top, xi, xi + n, Pinv);
    }
  }
  for (p = top; p < n; p++)
    CS_MARK(Lp, xi[p]); /* restore L */
  return (top);
}

/* x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse */
int cs_scatter(const cs *A,
               int j,
               double beta,
               int *w,
               double *x,
               int mark,
               cs *C,
               int nz)
{
  int i, p, *Ap, *Ai, *Ci;
  double *Ax;
  if (!A || !w || !C) return (-1); /* ensure inputs are valid */
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  Ci = C->i;
  for (p = Ap[j]; p < Ap[j + 1]; p++)
  {
    i = Ai[p]; /* A(i,j) is nonzero */
    if (w[i] < mark)
    {
      w[i] = mark; /* i is new entry in column j */
      Ci[nz++] = i; /* add i to pattern of C(:,j) */
      if (x) x[i] = beta * Ax[p]; /* x(i) = beta*A(i,j) */
    }
    else if (x) x[i] += beta * Ax[p]; /* i exists in C(:,j) already */
  }
  return (nz);
}

/* find the strongly connected components of a square matrix */
csd* cs_scc(cs *A) /* matrix A temporarily modified, then restored */
{
  int n, i, k, b = 0, top, *xi, *pstack, *P, *R, *Ap, *ATp;
  cs *AT;
  csd *D;
  if (!A) return (NULL);
  n = A->n;
  Ap = A->p;
  D = cs_dalloc(n, 0);
  AT = cs_transpose(A, 0); /* AT = A' */
  xi = (int*) cs_malloc(2 * n, sizeof(int)); /* allocate workspace */
  pstack = xi + n;
  if (!D || !AT || !xi) return (cs_ddone(D, AT, xi, 0));
  P = D->P;
  R = D->R;
  ATp = AT->p;
  top = n;
  for (i = 0; i < n; i++) /* first dfs(A) to find finish times (xi) */
  {
    if (!CS_MARKED(Ap, i)) top = cs_dfs(i, A, top, xi, pstack, NULL);
  }
  for (i = 0; i < n; i++)
    CS_MARK(Ap, i); /* restore A; unmark all nodes*/
  top = n;
  b = n;
  for (k = 0; k < n; k++) /* dfs(A') to find strongly connnected comp. */
  {
    i = xi[k]; /* get i in reverse order of finish times */
    if (CS_MARKED(ATp, i)) continue; /* skip node i if already ordered */
    R[b--] = top; /* node i is the start of a component in P */
    top = cs_dfs(i, AT, top, P, pstack, NULL);
  }
  R[b] = 0; /* first block starts at zero; shift R up */
  for (k = b; k <= n; k++)
    R[k - b] = R[k];
  D->nb = R[n + 1] = b = n - b; /* b = # of strongly connected components */
  return (cs_ddone(D, AT, xi, 1));
}

/* ordering and symbolic analysis for a Cholesky factorization */
css* cs_schol(const cs *A, int order)
{
  int n, *c, *post, *P;
  cs *C;
  css *S;
  if (!A) return (NULL); /* check inputs */
  n = A->n;
  S = (css*) cs_calloc(1, sizeof(css)); /* allocate symbolic analysis */
  if (!S) return (NULL); /* out of memory */
  P = cs_amd(A, order); /* P = amd(A+A'), or natural */
  S->Pinv = cs_pinv(P, n); /* find inverse permutation */
  cs_free(P);
  if (order >= 0 && !S->Pinv) return (cs_sfree(S));
  C = cs_symperm(A, S->Pinv, 0); /* C = spones(triu(A(P,P))) */
  S->parent = cs_etree(C, 0); /* find etree of C */
  post = cs_post(n, S->parent); /* postorder the etree */
  c = cs_counts(C, S->parent, post, 0); /* find column counts of chol(C) */
  cs_free(post);
  cs_spfree(C);
  S->cp = (int*) cs_malloc(n + 1, sizeof(int)); /* find column pointers for L */
  S->unz = S->lnz = cs_cumsum(S->cp, c, n);
  cs_free(c);
  return ((S->lnz >= 0) ? S : cs_sfree(S));
}

/* solve Lx=b(:,k), leaving pattern in xi[top..n-1], values scattered in x. */
int cs_splsolve(cs *L, const cs *B, int k, int *xi, double *x, const int *Pinv)
{
  int j, jnew, p, px, top, n, *Lp, *Li, *Bp, *Bi;
  double *Lx, *Bx;
  if (!L || !B || !xi || !x) return (-1);
  Lp = L->p;
  Li = L->i;
  Lx = L->x;
  n = L->n;
  Bp = B->p;
  Bi = B->i;
  Bx = B->x;
  top = cs_reach(L, B, k, xi, Pinv); /* xi[top..n-1]=Reach(B(:,k)) */
  for (p = top; p < n; p++)
    x[xi[p]] = 0;/* clear x */
  for (p = Bp[k]; p < Bp[k + 1]; p++)
    x[Bi[p]] = Bx[p]; /* scatter B */
  for (px = top; px < n; px++)
  {
    j = xi[px]; /* x(j) is nonzero */
    jnew = Pinv ? (Pinv[j]) : j; /* j is column jnew of L */
    if (jnew < 0) continue; /* column jnew is empty */
    for (p = Lp[jnew] + 1; p < Lp[jnew + 1]; p++)
    {
      x[Li[p]] -= Lx[p] * x[j]; /* x(i) -= L(i,j) * x(j) */
    }
  }
  return (top); /* return top of stack */
}

/* compute vnz, Pinv, leftmost, m2 from A and parent */
static int* cs_vcount(const cs *A, const int *parent, int *m2, int *vnz)
{
  int i, k, p, pa, n = A->n, m = A->m, *Ap = A->p, *Ai = A->i;
  int *Pinv = (int*) cs_malloc(2 * m + n, sizeof(int)), *leftmost = Pinv + m
                                                                    + n;
  int *w = (int*) cs_malloc(m + 3 * n, sizeof(int));
  int *next = w, *head = w + m, *tail = w + m + n, *nque = w + m + 2 * n;
  if (!Pinv || !w) return (cs_idone(Pinv, NULL, w, 0));
  for (k = 0; k < n; k++)
    head[k] = -1; /* queue k is empty */
  for (k = 0; k < n; k++)
    tail[k] = -1;
  for (k = 0; k < n; k++)
    nque[k] = 0;
  for (i = 0; i < m; i++)
    leftmost[i] = -1;
  for (k = n - 1; k >= 0; k--)
  {
    for (p = Ap[k]; p < Ap[k + 1]; p++)
    {
      leftmost[Ai[p]] = k; /* leftmost[i] = min(find(A(i,:)))*/
    }
  }
  for (i = m - 1; i >= 0; i--) /* scan rows in reverse order */
  {
    Pinv[i] = -1; /* row i is not yet ordered */
    k = leftmost[i];
    if (k == -1) continue; /* row i is empty */
    if (nque[k]++ == 0) tail[k] = i; /* first row in queue k */
    next[i] = head[k]; /* put i at head of queue k */
    head[k] = i;
  }
  (*vnz) = 0;
  (*m2) = m;
  for (k = 0; k < n; k++) /* find row permutation and nnz(V)*/
  {
    i = head[k]; /* remove row i from queue k */
    (*vnz)++; /* count V(k,k) as nonzero */
    if (i < 0) i = (*m2)++; /* add a fictitious row */
    Pinv[i] = k; /* associate row i with V(:,k) */
    if (--nque[k] <= 0) continue; /* skip if V(k+1:m,k) is empty */
    (*vnz) += nque[k]; /* nque [k] = nnz (V(k+1:m,k)) */
    if ((pa = parent[k]) != -1) /* move all rows to parent of k */
    {
      if (nque[pa] == 0) tail[pa] = tail[k];
      next[tail[k]] = head[pa];
      head[pa] = next[i];
      nque[pa] += nque[k];
    }
  }
  for (i = 0; i < m; i++)
    if (Pinv[i] < 0) Pinv[i] = k++;
  return (cs_idone(Pinv, NULL, w, 1));
}

/* symbolic analysis for QR or LU */
css* cs_sqr(const cs *A, int order, int qr)
{
  int n, k, ok = 1, *post;
  css *S;
  if (!A) return (NULL); /* check inputs */
  n = A->n;
  S = (css*) cs_calloc(1, sizeof(css)); /* allocate symbolic analysis */
  if (!S) return (NULL); /* out of memory */
  S->Q = cs_amd(A, order); /* fill-reducing ordering */
  if (order >= 0 && !S->Q) return (cs_sfree(S));
  if (qr) /* QR symbolic analysis */
  {
    cs *C = (order >= 0) ? cs_permute(A, NULL, S->Q, 0) : ((cs*) A);
    S->parent = cs_etree(C, 1); /* etree of C'*C, where C=A(:,Q) */
    post = cs_post(n, S->parent);
    S->cp = cs_counts(C, S->parent, post, 1); /* col counts chol(C'*C) */
    cs_free(post);
    ok = C && S->parent && S->cp;
    if (ok) S->Pinv = cs_vcount(C, S->parent, &(S->m2), &(S->lnz));
    ok = ok && S->Pinv;
    if (ok) for (S->unz = 0, k = 0; k < n; k++)
      S->unz += S->cp[k];
    if (order >= 0) cs_spfree(C);
  }
  else
  {
    S->unz = 4 * (A->p[n]) + n; /* for LU factorization only, */
    S->lnz = S->unz; /* guess nnz(L) and nnz(U) */
  }
  return (ok ? S : cs_sfree(S));
}

/* C = A(p,p) where A and C are symmetric the upper part stored, Pinv not P */
cs* cs_symperm(const cs *A, const int *Pinv, int values)
{
  int i, j, p, q, i2, j2, n, *Ap, *Ai, *Cp, *Ci, *w;
  double *Cx, *Ax;
  cs *C;
  if (!A) return (NULL);
  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  C = cs_spalloc(n, n, Ap[n], values && (Ax != NULL), 0);
  w = (int*) cs_calloc(n, sizeof(int));
  if (!C || !w)
  {
    messerr("Core allocation problem in CSparse Library (%d x %d)", n, n);
    return (cs_done(C, w, NULL, 0)); /* out of memory */
  }
  Cp = C->p;
  Ci = C->i;
  Cx = C->x;
  for (j = 0; j < n; j++) /* count entries in each column of C */
  {
    j2 = Pinv ? Pinv[j] : j; /* column j of A is column j2 of C */
    for (p = Ap[j]; p < Ap[j + 1]; p++)
    {
      i = Ai[p];
      if (i > j) continue; /* skip lower triangular part of A */
      i2 = Pinv ? Pinv[i] : i; /* row i of A is row i2 of C */
      w[CS_MAX(i2, j2)]++; /* column count of C */
    }
  }
  cs_cumsum(Cp, w, n); /* compute column pointers of C */
  for (j = 0; j < n; j++)
  {
    j2 = Pinv ? Pinv[j] : j; /* column j of A is column j2 of C */
    for (p = Ap[j]; p < Ap[j + 1]; p++)
    {
      i = Ai[p];
      if (i > j) continue; /* skip lower triangular part of A*/
      i2 = Pinv ? Pinv[i] : i; /* row i of A is row i2 of C */
      Ci[q = w[CS_MAX(i2, j2)]++] = CS_MIN(i2, j2);
      if (Cx) Cx[q] = Ax[p];
    }
  }
  return (cs_done(C, w, NULL, 1)); /* success; free workspace, return C */
}

/* depth-first search and postorder of a tree rooted at node j */
int cs_tdfs(int j, int k, int *head, const int *next, int *post, int *stack)
{
  int i, p, top = 0;
  if (!head || !next || !post || !stack) return (-1); /* check inputs */
  stack[0] = j; /* place j on the stack */
  while (top >= 0) /* while (stack is not empty) */
  {
    p = stack[top]; /* p = top of stack */
    i = head[p]; /* i = youngest child of p */
    if (i == -1)
    {
      top--; /* p has no unordered children left */
      post[k++] = p; /* node p is the kth postordered node */
    }
    else
    {
      head[p] = next[i]; /* remove i from children of p */
      stack[++top] = i; /* start dfs on child node i */
    }
  }
  return (k);
}

/* C = A' */
cs* cs_transpose(const cs *A, int values)
{
  int p, q, j, *Cp, *Ci, n, m, *Ap, *Ai, *w;
  double *Cx, *Ax;
  cs *C;
  if (!A) return (NULL);
  m = A->m;
  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  C = cs_spalloc(n, m, Ap[n], values && Ax, 0); /* allocate result */
  w = (int*) cs_calloc(m, sizeof(int));
  if (!C || !w || m <= 0 || n <= 0)
  {
    messerr("Problem when transposing a matrix in CSparse Library (%d x %d)", m,
            n);
    return (cs_done(C, w, NULL, 0));
  }
  Cp = C->p;
  Ci = C->i;
  Cx = C->x;
  for (p = 0; p < Ap[n]; p++)
    w[Ai[p]]++; /* row counts */
  cs_cumsum(Cp, w, m); /* row pointers */
  for (j = 0; j < n; j++)
  {
    for (p = Ap[j]; p < Ap[j + 1]; p++)
    {
      Ci[q = w[Ai[p]]++] = j; /* place A(i,j) as entry C(j,i) */
      if (Cx) Cx[q] = Ax[p];
    }
  }
  return (cs_done(C, w, NULL, 1)); /* success; free w and return C */
}

/* C = compressed-column form of a triplet matrix T */
cs* cs_triplet(const cs *T)
{
  int m, n, nz, p, k, *Cp, *Ci, *w, *Ti, *Tj;
  double *Cx, *Tx;
  cs *C;
  if (!T) return (NULL); /* check inputs */
  m = T->m;
  n = T->n;
  Ti = T->i;
  Tj = T->p;
  Tx = T->x;
  nz = T->nz;
  C = cs_spalloc(m, n, nz, Tx != NULL, 0); /* allocate result */
  w = (int*) cs_calloc(n, sizeof(int)); /* get workspace */
  if (!C || !w)
  {
    messerr("Core allocation problem in CSparse Library (%d x %d)", m, n);
    return (cs_done(C, w, NULL, 0)); /* out of memory */
  }
  Cp = C->p;
  Ci = C->i;
  Cx = C->x;
  for (k = 0; k < nz; k++)
    w[Tj[k]]++; /* count rows in column Tj[k]  */
  cs_cumsum(Cp, w, n); /* column pointers */
  for (k = 0; k < nz; k++)
  {
    Ci[p = w[Tj[k]]++] = Ti[k]; /* A(i,j) is the pth entry in C */
    if (Cx) Cx[p] = Tx[k];
  }
  return (cs_done(C, w, NULL, 1)); /* success; free w and return C */
}

/* sparse Cholesky update/downdate, L*L' + sigma*w*w' (sigma = +1 or -1) */
int cs_updown(cs *L, int sigma, const cs *C, const int *parent)
{
  int p, f, j, *Lp, *Li, *Cp, *Ci;
  double *Lx, *Cx, alpha, beta = 1, delta, gamma, w1, w2, *w, n, beta2 = 1;
  if (!L || !C || !parent) return (0);
  Lp = L->p;
  Li = L->i;
  Lx = L->x;
  n = L->n;
  Cp = C->p;
  Ci = C->i;
  Cx = C->x;
  if ((p = Cp[0]) >= Cp[1]) return (1); /* return if C empty */
  w = (double*) cs_malloc((int) n, sizeof(double));
  if (!w) return (0);
  f = Ci[p];
  for (; p < Cp[1]; p++)
    f = CS_MIN(f, Ci[p]); /* f = min (find (C)) */
  for (j = f; j != -1; j = parent[j])
    w[j] = 0; /* clear workspace w */
  for (p = Cp[0]; p < Cp[1]; p++)
    w[Ci[p]] = Cx[p]; /* w = C */
  for (j = f; j != -1; j = parent[j]) /* walk path f up to root */
  {
    p = Lp[j];
    alpha = w[j] / Lx[p]; /* alpha = w(j) / L(j,j) */
    beta2 = beta * beta + sigma * alpha * alpha;
    if (beta2 <= 0) break; /* not positive definite */
    beta2 = sqrt(beta2);
    delta = (sigma > 0) ? (beta / beta2) : (beta2 / beta);
    gamma = sigma * alpha / (beta2 * beta);
    Lx[p] = delta * Lx[p] + ((sigma > 0) ? (gamma * w[j]) : 0);
    beta = beta2;
    for (p++; p < Lp[j + 1]; p++)
    {
      w1 = w[Li[p]];
      w[Li[p]] = w2 = w1 - alpha * Lx[p];
      Lx[p] = delta * Lx[p] + gamma * ((sigma > 0) ? w1 : w2);
    }
  }
  cs_free(w);
  return (beta2 > 0);
}

/* solve Ux=b where x and b are dense.   x=b on input, solution on output. */
int cs_usolve(const cs *U, double *x)
{
  int p, j, n, *Up, *Ui;
  double *Ux;
  if (!U || !x) return (0); /* check inputs */
  n = U->n;
  Up = U->p;
  Ui = U->i;
  Ux = U->x;
  for (j = n - 1; j >= 0; j--)
  {
    x[j] /= Ux[Up[j + 1] - 1];
    for (p = Up[j]; p < Up[j + 1] - 1; p++)
    {
      x[Ui[p]] -= Ux[p] * x[j];
    }
  }
  return (1);
}

/* allocate a sparse matrix (triplet form or compressed-column form) */
cs* cs_spalloc(int m, int n, int nzmax, int values, int triplet)
{
  cs *A = (cs*) cs_calloc(1, sizeof(cs)); /* allocate the cs struct */
  if (!A) return (NULL); /* out of memory */
  A->m = m; /* define dimensions and nzmax */
  A->n = n;
  A->nzmax = nzmax = CS_MAX(nzmax, 1);
  A->nz = triplet ? 0 : -1; /* allocate triplet or comp.col */
  A->p = (int*) cs_malloc(triplet ? nzmax : n + 1, sizeof(int));
  A->i = (int*) cs_malloc(nzmax, sizeof(int));
  A->x = values ? (double*) cs_malloc(nzmax, sizeof(double)) : NULL;
  return ((!A->p || !A->i || (values && !A->x)) ? cs_spfree(A) : A);
}

/* change the max # of entries sparse matrix */
int cs_sprealloc(cs *A, int nzmax)
{
  int ok, oki, okj = 1, okx = 1;
  if (!A) return (0);
  nzmax = (nzmax <= 0) ? (A->p[A->n]) : nzmax;
  A->i = (int*) cs_realloc(A->i, nzmax, sizeof(int), &oki);
  if (A->nz >= 0) A->p = (int*) cs_realloc(A->p, nzmax, sizeof(int), &okj);
  if (A->x) A->x = (double*) cs_realloc(A->x, nzmax, sizeof(double), &okx);
  ok = (oki && okj && okx);
  if (ok) A->nzmax = nzmax;
  return (ok);
}

/* free a sparse matrix */
cs* cs_spfree(cs *A)
{
  if (!A) return (NULL); /* do nothing if A already NULL */
  cs_free(A->p);
  cs_free(A->i);
  cs_free(A->x);
  return ((cs*) cs_free(A)); /* free the cs struct and return NULL */
}

/* free a numeric factorization */
csn* cs_nfree(csn *N)
{
  if (!N) return (NULL); /* do nothing if N already NULL */
  cs_spfree(N->L);
  cs_spfree(N->U);
  cs_free(N->Pinv);
  cs_free(N->B);
  return ((csn*) cs_free(N)); /* free the csn struct and return NULL */
}

/* free a symbolic factorization */
css* cs_sfree(css *S)
{
  if (!S) return (NULL); /* do nothing if S already NULL */
  cs_free(S->Pinv);
  cs_free(S->Q);
  cs_free(S->parent);
  cs_free(S->cp);
  return ((css*) cs_free(S)); /* free the css struct and return NULL */
}

/* allocate a cs_dmperm or cs_scc result */
csd* cs_dalloc(int m, int n)
{
  csd *D;
  D = (csd*) cs_calloc(1, sizeof(csd));
  if (!D) return (NULL);
  D->P = (int*) cs_malloc(m, sizeof(int));
  D->R = (int*) cs_malloc(m + 6, sizeof(int));
  D->Q = (int*) cs_malloc(n, sizeof(int));
  D->S = (int*) cs_malloc(n + 6, sizeof(int));
  return ((!D->P || !D->R || !D->Q || !D->S) ? cs_dfree(D) : D);
}

/* free a cs_dmperm or cs_scc result */
csd* cs_dfree(csd *D)
{
  if (!D) return (NULL); /* do nothing if D already NULL */
  cs_free(D->P);
  cs_free(D->Q);
  cs_free(D->R);
  cs_free(D->S);
  return ((csd*) cs_free(D));
}

/* free workspace and return a sparse matrix result */
cs* cs_done(cs *C, void *w, void *x, int ok)
{
  cs_free(w); /* free workspace */
  cs_free(x);
  return (ok ? C : cs_spfree(C)); /* return result if OK, else free it */
}

/* free workspace and return int array result */
int* cs_idone(int *p, cs *C, void *w, int ok)
{
  cs_spfree(C); /* free temporary matrix */
  cs_free(w); /* free workspace */
  return (ok ? p : (int*) cs_free(p)); /* return result if OK, else free it */
}

/* free workspace and return a numeric factorization (Cholesky, LU, or QR) */
csn* cs_ndone(csn *N, cs *C, void *w, void *x, int ok)
{
  cs_spfree(C); /* free temporary matrix */
  cs_free(w); /* free workspace */
  cs_free(x);
  return (ok ? N : cs_nfree(N)); /* return result if OK, else free it */
}

/* free workspace and return a csd result */
csd* cs_ddone(csd *D, cs *C, void *w, int ok)
{
  cs_spfree(C); /* free temporary matrix */
  cs_free(w); /* free workspace */
  return (ok ? D : cs_dfree(D)); /* return result if OK, else free it */
}

/* solve U'x=b where x and b are dense.  x=b on input, solution on output. */
int cs_utsolve(const cs *U, double *x)
{
  int p, j, n, *Up, *Ui;
  double *Ux;
  if (!U || !x) return (0); /* check inputs */
  n = U->n;
  Up = U->p;
  Ui = U->i;
  Ux = U->x;
  for (j = 0; j < n; j++)
  {
    for (p = Up[j]; p < Up[j + 1] - 1; p++)
    {
      x[j] -= Ux[p] * x[Ui[p]];
    }
    x[j] /= Ux[p];
  }
  return (1);
}

/* sparseinv: computes the sparse inverse subset, using Takahashi's equations.

 On input, the pattern of Z must be equal to the symbolic Cholesky
 factorization of A+A', where A=(L+I)*(U+I).   The pattern of L+U must be a
 subset of Z.  Z must have zero-free diagonal.  These conditions are
 difficult to check, so they are assumed to hold.  Results will be completely
 wrong if the conditions do not hold.

 This function performs the same amount of work as the initial LU
 factorization, assuming that the pattern of P*A*Q is symmetric.  For large
 matrices, this function can take a lot more time than LU in MATLAB, even if
 P*A*Q is symmetric.   This is because LU is a multifrontal method, whereas
 this sparseinv function is based on gather/scatter operations.

 The basic integer type is an Int, or ptrdiff_t, which is 32 bits on a 32
 bits and 64 bits on a 64 bit system.  The function returns the flop count as
 an Int.  This will not overflow on a 64 bit system but might on a 32 bit.
 The total work is flops + O(n + nnz(Z)).  Since flops > n and flops > nnz(Z),
 this is O(flops).

 Copyright 2011, Timothy A. Davis, http://www.suitesparse.com
 */
int sparseinv /* returns -1 on error, or flop count if OK */
(
/* inputs, not modified on output: */
int n, /* L, U, D, and Z are n-by-n */
 int *Lp, /* L is sparse, lower triangular, stored by column */
 int *Li, /* the row indices of L must be sorted */
 double *Lx, /* diagonal of L, if present, is ignored */
 double *d, /* diagonal of D, of size n */
 int *Up, /* U is sparse, upper triangular, stored by row */
 int *Uj, /* the column indices of U need not be sorted */
 double *Ux, /* diagonal of U, if present, is ignored */
 int *Zp, /* Z is sparse, stored by column */
 int *Zi, /* the row indices of Z must be sorted */

 /* output, not defined on input: */
 double *Zx,

 /* workspace: */
 double *z, /* size n, zero on input, restored as such on output */
 int *Zdiagp, /* size n */
 int *Lmunch /* size n */
 )
{
  double ljk, zkj;
  int j, i, k, p, znz, pdiag, up, zp, flops = n;

  /* ---------------------------------------------------------------------- */
  /* initializations */
  /* ---------------------------------------------------------------------- */

  /* clear the numerical values of Z */
  znz = Zp[n];
  for (p = 0; p < znz; p++)
  {
    Zx[p] = 0;
  }

  /* find the diagonal of Z and initialize it */
  for (j = 0; j < n; j++)
  {
    pdiag = -1;
    for (p = Zp[j]; p < Zp[j + 1] && pdiag == -1; p++)
    {
      if (Zi[p] == j)
      {
        pdiag = p;
        Zx[p] = 1 / d[j];
      }
    }
    Zdiagp[j] = pdiag;
    if (pdiag == -1) return (-1); /* Z must have a zero-free diagonal */
  }

  /* Lmunch [k] points to the last entry in column k of L */
  for (k = 0; k < n; k++)
  {
    Lmunch[k] = Lp[k + 1] - 1;
  }

  /* ---------------------------------------------------------------------- */
  /* compute the sparse inverse subset */
  /* ---------------------------------------------------------------------- */

  for (j = n - 1; j >= 0; j--)
  {

    /* ------------------------------------------------------------------ */
    /* scatter Z (:,j) into z workspace */
    /* ------------------------------------------------------------------ */

    /* only the lower triangular part is needed, since the upper triangular
     part is all zero */
    for (p = Zdiagp[j]; p < Zp[j + 1]; p++)
    {
      z[Zi[p]] = Zx[p];
    }

    /* ------------------------------------------------------------------ */
    /* compute the strictly upper triangular part of Z (:,j) */
    /* ------------------------------------------------------------------ */

    /* for k = (j-1):-1:1 but only for the entries Z(k,j) */
    for (p = Zdiagp[j] - 1; p >= Zp[j]; p--)
    {
      /* Z (k,j) = - U (k,k+1:n) * Z (k+1:n,j) */
      k = Zi[p];
      zkj = 0;
      flops += (Up[k + 1] - Up[k]);
      for (up = Up[k]; up < Up[k + 1]; up++)
      {
        /* skip the diagonal of U, if present */
        i = Uj[up];
        if (i > k)
        {
          zkj -= Ux[up] * z[i];
        }
      }
      z[k] = zkj;
    }

    /* ------------------------------------------------------------------ */
    /* left-looking update to lower triangular part of Z */
    /* ------------------------------------------------------------------ */

    /* for k = (j-1):-1:1 but only for the entries Z(k,j) */
    for (p = Zdiagp[j] - 1; p >= Zp[j]; p--)
    {
      k = Zi[p];

      /* ljk = L (j,k) */
      if (Lmunch[k] < Lp[k] || Li[Lmunch[k]] != j)
      {
        /* L (j,k) is zero, so there is no work to do */
        continue;
      }
      ljk = Lx[Lmunch[k]--];

      /* Z (k+1:n,k) = Z (k+1:n,k) - Z (k+1:n,j) * L (j,k) */
      flops += (Zp[k + 1] - Zdiagp[k]);
      for (zp = Zdiagp[k]; zp < Zp[k + 1]; zp++)
      {
        Zx[zp] -= z[Zi[zp]] * ljk;
      }
    }

    /* ------------------------------------------------------------------ */
    /* gather Z (:,j) back from z workspace */
    /* ------------------------------------------------------------------ */

    for (p = Zp[j]; p < Zp[j + 1]; p++)
    {
      i = Zi[p];
      Zx[p] = z[i];
      z[i] = 0;
    }
  }
  return (flops);
}

/* compressed-column form into arrays */
/* number: Number of non-zero terms in the sparse matrix A */
void cs_sparse_to_triplet(const cs *A,
                          int flag_from_1,
                          int *number,
                          int **cols,
                          int **rows,
                          double **vals)
{
  int p, j, n, nz, nnz, ecr, *Ap, *Ai, ok;
  double *Ax;

  *number = 0;
  *cols = *rows = nullptr;
  *vals = nullptr;
  if (!A) return;
  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  nz = A->nz;
  if (nz >= 0) return;

  nnz = Ap[n];
  (*cols) = (int*) cs_malloc(nnz, sizeof(int));
  (*rows) = (int*) cs_malloc(nnz, sizeof(int));
  (*vals) = (double*) cs_malloc(nnz, sizeof(double));

  ecr = 0;
  for (j = 0; j < n; j++)
    for (p = Ap[j]; p < Ap[j + 1]; p++)
    {
      (*cols)[ecr] = (flag_from_1) ? j + 1 : j;
      (*rows)[ecr] = (flag_from_1) ? Ai[p] + 1 : Ai[p];
      (*vals)[ecr] = Ax ? Ax[p] : 1;
      if (ABS((*vals)[ecr]) <= 0) continue;
      ecr++;
    }

  if (ecr < nnz)
  {
    (*cols) = (int*) cs_realloc(*cols, ecr, sizeof(int), &ok);
    (*rows) = (int*) cs_realloc(*rows, ecr, sizeof(int), &ok);
    (*vals) = (double*) cs_realloc(*vals, ecr, sizeof(double), &ok);
  }
  *number = ecr;
  return;
}

cs* cs_arrays_to_sparse(int n,
                        int nrow,
                        int ncol,
                        double *rows,
                        double *cols,
                        double *vals)
{
  cs *Q, *Qtriplet;
  int ip1, ip2, error, row_max, col_max;

  error = 1;
  Q = Qtriplet = nullptr;

  row_max = col_max = -1;
  Qtriplet = cs_spalloc(0, 0, 1, 1, 1);
  for (int i = 0; i < n; i++)
  {
    ip1 = (int) rows[i];
    ip2 = (int) cols[i];
    if (ip1 > row_max) row_max = ip1;
    if (ip2 > col_max) col_max = ip2;
    if (!cs_entry(Qtriplet, ip1, ip2, vals[i])) goto label_end;
  }

  if (nrow > 0 && row_max >= nrow)
  {
    messerr("Inconsistency between number of rows (%d)", nrow);
    messerr("and maximum row index in argument 'rows' (%d)", row_max);
    goto label_end;
  }
  if (ncol > 0 && col_max >= ncol)
  {
    messerr("Inconsistency between number of columns (%d)", ncol);
    messerr("and maximum column index in argument 'cols' (%d)", col_max);
    goto label_end;
  }
  if (row_max < nrow - 1 || col_max < ncol - 1)
  {
    /* Add a fictitious entry to ensure the dimension of the sparse matrix */

    if (!cs_entry(Qtriplet, nrow - 1, ncol - 1, 0.)) goto label_end;
  }
  Q = cs_triplet(Qtriplet);
  if (Q == nullptr) goto label_end;

  error = 0;

  label_end: Qtriplet = cs_spfree(Qtriplet);
  if (error) Q = cs_spfree(Q);
  return (Q);
}

/* Extract a sparse submatrix */

// row_from and col_from are given strating from 0
cs* cs_extract_submatrix(cs *C,
                         int row_from,
                         int row_length,
                         int col_from,
                         int col_length)
{
  cs *Atriplet, *A;
  int *cols, *rows, number, ir, ic;
  double *vals;

  /* Initializations */

  cols = rows = nullptr;
  vals = nullptr;
  A = Atriplet = nullptr;

  /* Convert the contents of the sparse matrix into columns */

  cs_sparse_to_triplet(C, 0, &number, &cols, &rows, &vals);

  /* Fill the new sparse triplet */

  Atriplet = cs_spalloc(0, 0, 1, 1, 1);
  if (Atriplet == nullptr) goto label_end;
  for (int i = 0; i < number; i++)
  {
    ic = cols[i] - col_from;
    if (ic < 0 || ic >= col_length) continue;
    ir = rows[i] - row_from;
    if (ir < 0 || ir >= row_length) continue;
    if (!cs_entry(Atriplet, ir, ic, vals[i])) goto label_end;
  }

  A = cs_triplet(Atriplet);

  label_end: Atriplet = cs_spfree(Atriplet);
  cols = (int*) cs_free((char*) cols);
  rows = (int*) cs_free((char*) rows);
  vals = (double*) cs_free((char*) vals);
  return (A);
}

/* Extract a sparse submatrix */
/* 'rank_rows' and 'rank_cols' must have same dimension as C */
/* The arrays 'rank_rows' and 'rank_cols' may be absent */
/* Their value gives the rank of the saved element or -1 */
cs* cs_extract_submatrix_by_ranks(cs *C, int *rank_rows, int *rank_cols)
{
  cs *Atriplet, *A;
  int *cols, *rows, number, old_row, old_col, new_row, new_col;
  double *vals;

  /* Initializations */

  cols = rows = nullptr;
  vals = nullptr;
  A = Atriplet = nullptr;

  /* Convert the contents of the sparse matrix into columns */

  cs_sparse_to_triplet(C, 0, &number, &cols, &rows, &vals);

  /* Initialize the output matrix */

  Atriplet = cs_spalloc(0, 0, 1, 1, 1);
  if (Atriplet == nullptr) goto label_end;

  /* Fill the new sparse triplet */

  for (int i = 0; i < number; i++)
  {
    old_row = rows[i];
    old_col = cols[i];
    new_row = (rank_rows != nullptr) ? rank_rows[old_row] : old_row;
    new_col = (rank_cols != nullptr) ? rank_cols[old_col] : old_col;
    if (new_row < 0 || new_col < 0) continue;
    if (!cs_entry(Atriplet, new_row, new_col, vals[i])) goto label_end;
  }

  A = cs_triplet(Atriplet);

  label_end: cols = (int*) cs_free((char*) cols);
  rows = (int*) cs_free((char*) rows);
  vals = (double*) cs_free((char*) vals);
  Atriplet = cs_spfree(Atriplet);
  return (A);
}

/* Extract a sparse submatrix */
/* The array 'colors' has the same dimension as C */
/* The element of 'C' must be kept if: */
/* - the color of its row number is equal to 'ref_color' if 'row_ok'==TRUE */
/*   or different if 'row_ok'==FALSE */
/* and if */
/* - the color of its column number is equal to 'ref-color' if 'col_ok'==TRUE*/
/*   or different if 'col_ok'== FALSE */
cs* cs_extract_submatrix_by_color(cs *C,
                                  int *colors,
                                  int ref_color,
                                  int row_ok,
                                  int col_ok)
{
  cs *Atriplet, *A;
  int *cols, *rows, *u_row, *u_col, number, ir, ic, n;
  double *vals;

  /* Initializations */

  cols = rows = u_row = u_col = nullptr;
  vals = nullptr;
  A = Atriplet = nullptr;

  /* Convert the contents of the sparse matrix into columns */

  cs_sparse_to_triplet(C, 0, &number, &cols, &rows, &vals);

  /* Initialize the output matrix */

  Atriplet = cs_spalloc(0, 0, 1, 1, 1);
  if (Atriplet == nullptr) goto label_end;

  /* Core allocation */

  n = C->n;
  u_row = (int*) mem_alloc(sizeof(int) * n, 1);
  u_col = (int*) mem_alloc(sizeof(int) * n, 1);

  ir = 0;
  for (int i = 0; i < n; i++)
  {
    u_row[i] = -1;
    if (row_ok && colors[i] != ref_color) continue;
    if (!row_ok && colors[i] == ref_color) continue;
    u_row[i] = ir++;
  }

  ic = 0;
  for (int i = 0; i < n; i++)
  {
    u_col[i] = -1;
    if (col_ok && colors[i] != ref_color) continue;
    if (!col_ok && colors[i] == ref_color) continue;
    u_col[i] = ic++;
  }

  /* Fill the new sparse triplet */

  for (int i = 0; i < number; i++)
  {
    ir = u_row[rows[i]];
    ic = u_col[cols[i]];
    if (ir < 0 || ic < 0) continue;
    if (!cs_entry(Atriplet, ir, ic, vals[i])) goto label_end;
  }

  A = cs_triplet(Atriplet);

  label_end: u_row = (int*) mem_free((char* ) u_row);
  u_col = (int*) mem_free((char* ) u_col);
  cols = (int*) cs_free((char*) cols);
  rows = (int*) cs_free((char*) rows);
  vals = (double*) cs_free((char*) vals);
  Atriplet = cs_spfree(Atriplet);
  return (A);
}

int cs_get_nrow(const cs *A)
{
  if (A == nullptr) return 0;
  cs *AT = cs_transpose(A, 1);
  if (AT == nullptr) return 0;
  int nrow = AT->n;
  AT = cs_spfree(AT);
  return (nrow);
}

int cs_get_ncol(const cs *A)
{
  if (A == nullptr) return 0;
  int ncol = A->n;
  return (ncol);
}

int cs_get_ncell(const cs *A)
{
  if (A == nullptr) return 0;
  int ncol = A->n;
  cs *AT = cs_transpose(A, 1);
  int nrow = AT->n;
  AT = cs_spfree(AT);
  return (nrow * ncol);
}

void cs_print_dim(const char *title, const cs *A)
{
  cs *AT;
  int n1, n2;

  n1 = n2 = 0;
  if (A == nullptr) return;
  n1 = A->n;
  AT = cs_transpose(A, 1);
  if (AT != nullptr) n2 = AT->n;
  message("%s: Nrow=%d Ncol=%d\n", title, n2, n1);
  if (AT != nullptr) AT = cs_spfree(AT);
}

String toStringDim(const String &title, const cs *A)
{
  std::stringstream sstr;
  cs *AT;
  int n1, n2;

  n1 = n2 = 0;
  if (A == nullptr) return sstr.str();
  n1 = A->n;
  AT = cs_transpose(A, 1);
  if (AT != nullptr) n2 = AT->n;
  if (!title.empty()) sstr << title << " : ";
  sstr << "Nrows=" << n2 << " - Ncols=" << n1 << std::endl;
  if (AT != nullptr) AT = cs_spfree(AT);
  return sstr.str();
}

void cs_print_range(const char *title, const cs *C)
{
  int *cols, *rows, number, nvalid;
  double *vals, mini, maxi;

  /* Initializations */

  if (C == nullptr) return;
  cols = rows = nullptr;
  vals = nullptr;

  /* Convert the contents of the sparse matrix into columns */

  cs_sparse_to_triplet(C, 0, &number, &cols, &rows, &vals);

  /* Calculate the extreme values */

  nvalid = 0;
  mini = maxi = TEST;
  ut_stats_mima(number, vals, NULL, &nvalid, &mini, &maxi);

  /* Printout */

  if (title != NULL)
    message("%s\n", title);
  else
    message("Sparse matrix\n");
  message(" Descr: m=%d n=%d nzmax=%d\n", C->m, C->n, C->nzmax);
  message(" Range: [%lf ; %lf] (%d/%d)\n", mini, maxi, nvalid, number);

  /* Core deallocation */

  cols = (int*) cs_free((char*) cols);
  rows = (int*) cs_free((char*) rows);
  vals = (double*) cs_free((char*) vals);
}

String toStringRange(const String &title, const cs *C)
{
  std::stringstream sstr;
  int *cols, *rows, number, nvalid;
  double *vals, mini, maxi;

  /* Initializations */

  if (C == nullptr) return sstr.str();
  cols = rows = nullptr;
  vals = nullptr;

  /* Convert the contents of the sparse matrix into columns */

  cs_sparse_to_triplet(C, 0, &number, &cols, &rows, &vals);

  /* Calculate the extreme values */

  nvalid = 0;
  mini = maxi = TEST;
  ut_stats_mima(number, vals, NULL, &nvalid, &mini, &maxi);

  /* Printout */

  if (!title.empty()) sstr << title << std::endl;

  sstr << " Descr: m=" << C->m << " - n=" << C->n << " - nzmax=" << C->nzmax
  << std::endl;
  sstr << " Range: [" << mini << " ; " << maxi << "] (" << nvalid << " / "
       << number << ")" << std::endl;

  /* Core deallocation */

  cols = (int*) cs_free((char*) cols);
  rows = (int*) cs_free((char*) rows);
  vals = (double*) cs_free((char*) vals);

  return sstr.str();
}

/* Construct the sparse diagonal matrix */
cs* cs_eye(int number, double value)
{
  cs *Atriplet, *A;

  /* Initializations */

  A = nullptr;

  /* Fill the new sparse triplet */

  Atriplet = cs_spalloc(0, 0, 1, 1, 1);
  if (Atriplet == nullptr) goto label_end;
  for (int i = 0; i < number; i++)
  {
    if (!cs_entry(Atriplet, i, i, value)) goto label_end;
  }

  A = cs_triplet(Atriplet);

  label_end: Atriplet = cs_spfree(Atriplet);
  return (A);
}

// Build a sparse matrix containing the diagonal of a sparse matrix
// mode: Operation on the diagonal term
//  1 for diagonal
// -1 for inverse diagonal
//  2 for squared diagonal
// -2 for inverse square root of diagonal

cs* cs_extract_diag(cs *C, int mode)
{
  cs *Atriplet, *A;
  double *Cx, value;
  int *Cp, *Ci;

  /* Initializations */

  A = nullptr;
  Atriplet = cs_spalloc(0, 0, 1, 1, 1);
  if (Atriplet == nullptr) goto label_end;
  Cp = C->p;
  Ci = C->i;
  Cx = C->x;

  /* Loop on the rows */

  for (int j = 0; j < C->n; j++)
  {
    for (int p = Cp[j]; p < Cp[j + 1]; p++)
    {
      if (Ci[p] != j) continue;
      value = Cx[p];

      switch (mode)
      {
        case 1:
          if (!cs_entry(Atriplet, j, j, value)) goto label_end;
          break;

        case -1:
          if (!cs_entry(Atriplet, j, j, 1. / value)) goto label_end;
          break;

        case 2:
          if (!cs_entry(Atriplet, j, j, value * value)) goto label_end;
          break;

        case -2:
          if (!cs_entry(Atriplet, j, j, 1. / sqrt(value))) goto label_end;
          break;
      }
    }
  }
  A = cs_triplet(Atriplet);

  label_end: Atriplet = cs_spfree(Atriplet);
  return (A);
}

// Extract the (transformed) diagonal of a sparse matrix
// mode: Operation on the diagonal term
//  1 for diagonal
// -1 for inverse diagonal
//  2 for squared diagonal
// -2 for inverse square root of diagonal

double* csd_extract_diag(cs *C, int mode)
{
  int *Cp, *Ci, size;
  double *Cx, *diag, value;

  /* Initializations */

  diag = nullptr;
  size = C->n;
  Cp = C->p;
  Ci = C->i;
  Cx = C->x;

  /* Core allocation */

  diag = (double*) mem_alloc(sizeof(double) * size, 0);
  if (diag == nullptr) goto label_end;
  for (int i = 0; i < size; i++)
    diag[i] = 0.;

  /* Loop on the rows */

  for (int j = 0; j < C->n; j++)
  {
    for (int p = Cp[j]; p < Cp[j + 1]; p++)
    {
      if (Ci[p] != j) continue;
      value = Cx[p];

      switch (mode)
      {
        case 1:
          diag[j] = value;
          break;

        case -1:
          diag[j] = 1. / value;
          break;

        case 2:
          diag[j] = value * value;
          break;

        case -2:
          diag[j] = 1. / sqrt(value);
          break;
      }
    }
  }

  label_end: return (diag);
}

void cs_diag_suppress(cs *C)
{
  int *Ci, *Cp;
  double *Cx;

  /* Initializations */

  Cp = C->p;
  Ci = C->i;
  Cx = C->x;

  /* Loop on the rows */

  for (int j = 0; j < C->n; j++)
  {
    for (int p = Cp[j]; p < Cp[j + 1]; p++)
    {
      if (Ci[p] == j) Cx[p] = 0.;
    }
  }
  return;
}

int cs_sort_i(cs *C)
{
  int *rank, size, n, j, i, ecr;

  /* Core allocation */

  n = C->n;
  size = MAX(C->m, n);
  rank = (int*) mem_alloc(sizeof(int) * size, 0);
  if (rank == nullptr) return (1);

  for (j = 0; j < n; j++)
  {
    for (i = C->p[j], ecr = 0; i < C->p[j + 1]; i++)
      rank[ecr++] = C->i[i];
    ut_sort_int(0, ecr, NULL, rank);
    for (i = C->p[j], ecr = 0; i < C->p[j + 1]; i++)
      C->i[i] = rank[ecr++];
  }

  rank = (int*) mem_free((char* ) rank);
  return (0);
}

/* Return the number of rows and columns */
/* as well as the percentage of filled terms */
void cs_rowcol(const cs *A, int *nrows, int *ncols, int *count, double *percent)
{
  cs *AT;
  int *Ap;
  double *Ax;

  (*nrows) = (*ncols) = (*count) = 0;
  (*percent) = 0.;
  if (!A) return;

  Ap = A->p;
  Ax = A->x;
  if (A->nz >= 0) return;

  for (int j = 0; j < A->n; j++)
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
    {
      if (ABS(Ax[p]) > 0) (*count)++;
    }

  *ncols = A->n;
  AT = cs_transpose(A, 1);
  *nrows = AT->n;
  AT = cs_spfree(AT);

  if ((*nrows) > 0 && (*ncols) > 0)
    (*percent) = ((100. * (double) (*count))
        / ((double) (*nrows) * (double) (*ncols)));
}

/* Print a nice sparse matrix */
/* The format is copied from Matrix package */
void cs_print_nice(const char *title, const cs *A, int maxrow, int maxcol)
{
  int p, j, m, n, *Ap, *Ai, npass, jdeb, jfin, found, num_line;
  double *Ax;
  int nbypass = 7;

  if (!A)
  {
    message("(null)\n");
    return;
  }
  m = A->m;
  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  if (A->nz >= 0) return;
  if (maxcol >= 0) n = maxcol;
  if (maxrow >= 0) m = maxrow;

  npass = (int) ceil((double) n / (double) nbypass);

  /* Print the title (optional) */

  if (title != NULL)
    message("%s", title);
  else
    message("Print Sparse Matrix");
  if (maxrow >= 0) message(" nrows<=%d", maxrow);
  if (maxcol >= 0) message(" ncols<=%d", maxcol);
  message("\n");

  /* Loop on the passes */

  for (int ipass = 0; ipass < npass; ipass++)
  {
    jdeb = ipass * nbypass;
    jfin = MIN(jdeb + nbypass, n);

    /* Title of the columns */

    message("      ");
    for (j = jdeb; j < jfin; j++)
      message("    [,%3d]", j + 1);
    message("\n");

    /* Loop on the lines */

    for (int i = 0; i < m; i++)
    {
      message("[%3d,] ", i + 1);

      /* Loop on the columns */

      for (j = jdeb; j < jfin; j++)
      {

        /* Search for the correct line number */

        found = -1;
        for (p = Ap[j]; p < Ap[j + 1] && found < 0; p++)
        {
          num_line = Ai[p];
          if (num_line == i) found = p;
        }

        if (found < 0)
          message(" .        ");
        else
          message("%9.4lf ", Ax[found]);
      }
      message("\n");
    }
    message("\n");
  }
  return;
}

/* Print the only non-zero terms of a sparse matrix */
void cs_print_only(const char *title, const cs *A, int nlimit)
{
  int p, j, n, *Ap, *Ai, ecr;
  double *Ax, value;

  if (!A)
  {
    message("(null)\n");
    return;
  }
  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;

  /* Print the title (optional) */

  if (title != NULL) message("Only non-zero terms in %s\n", title);

  /* Loop on the elements */

  ecr = 0;
  for (j = 0; j < n; j++)
  {
    for (p = Ap[j]; p < Ap[j + 1]; p++)
    {
      value = Ax[p];
      if (ABS(value) <= 1.e-6) continue;
      message("i=%5d j=%5d Value = %lf\n", Ai[p], j, value);
      ecr++;
      if (nlimit > 0 && ecr > nlimit) return;
    }
  }
  return;
}

void cs_print_short(const char *title, const cs *L, int nmax)
{
  int p, j, i, n, *Lp, *Li;
  double *Lx, value;

  n = L->n;
  Lp = L->p;
  Li = L->i;
  Lx = L->x;

  /* Print the title (optional) */

  if (title != NULL) message("\n%s\n", title);

  for (j = 0; j < MIN(n, nmax); j++)
  {
    message("[%d] - ", j + 1);
    for (p = Lp[j]; p < Lp[j + 1]; p++)
    {
      i = Li[p];
      value = Lx[p];
      if (ABS(value) > 1.e-10) message("[%d] %7.4lf ", i + 1, Lx[p]);
    }
    message("\n");
  }
}

/* Construct the sparse diagonal matrix from a vector of values */
cs* cs_eye_tab(int number, double *values)
{
  cs *Atriplet, *A;

  /* Initializations */

  A = nullptr;

  /* Fill the new sparse triplet */

  Atriplet = cs_spalloc(0, 0, 1, 1, 1);
  if (Atriplet == nullptr) goto label_end;
  for (int i = 0; i < number; i++)
  {
    if (!cs_entry(Atriplet, i, i, values[i])) goto label_end;
  }

  A = cs_triplet(Atriplet);

  label_end: Atriplet = cs_spfree(Atriplet);
  return (A);
}

cs* cs_multiply_and_release(cs *b1, cs *b2, int flag_release)
{
  cs *bres;

  bres = cs_multiply(b1, b2);
  if (bres == nullptr) goto label_end;
  if (flag_release) b1 = cs_spfree(b1);

  label_end: return (bres);
}

cs* cs_add_and_release(cs *b1,
                       cs *b2,
                       double alpha,
                       double beta,
                       int flag_release)
{
  cs *bres;

  bres = cs_add(b1, b2, alpha, beta);
  if (bres == nullptr) goto label_end;
  if (flag_release) b1 = cs_spfree(b1);

  label_end: return (bres);
}

cs* cs_duplicate(const cs *b1)
{
  cs *bres;
  bres = cs_add(b1, b1, 1., 0.);
  return (bres);
}

cs* cs_prod_norm_and_release(cs *b1, cs *lambda, int flag_release)
{
  cs *bres1, *bres2;

  bres1 = bres2 = nullptr;

  bres1 = cs_multiply(lambda, b1);
  if (bres1 == nullptr) goto label_end;
  bres2 = cs_multiply(bres1, lambda);
  if (bres2 == nullptr) goto label_end;
  if (flag_release) b1 = cs_spfree(b1);

  label_end: bres1 = cs_spfree(bres1);
  return (bres2);
}

/* Return per column, the sum along the row */
/* Note: Check existence of returned argument */
/* Note: The returned array must be freed by calling function */
double* cs_col_sumrow(const cs *A, int *ncol, int *nrow)
{
  int *Ap, error;
  cs *AT;      // Will store t(A)
  double *Ax, *vect, sumval;

  error = 1;
  *ncol = *nrow = 0;
  AT = nullptr;
  vect = nullptr;
  if (!A) goto label_end;
  if (A->nz >= 0) goto label_end;
  Ap = A->p;
  Ax = A->x;

  /* Transpose A into AT to follow ordering per column */

  AT = cs_transpose(A, 1);
  if (AT == nullptr) goto label_end;

  /* Core allocation */

  vect = (double*) mem_alloc(sizeof(double) * A->n, 0);
  if (vect == nullptr) goto label_end;

  /* Loop on the rows */

  for (int j = 0; j < A->n; j++)
  {
    sumval = 0.;
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
      sumval += Ax[p];
    vect[j] = sumval;
  }

  *ncol = A->n;
  *nrow = AT->n;
  error = 0;

  label_end: AT = cs_spfree(AT);
  if (error) vect = (double*) mem_free((char* ) vect);
  return (vect);
}

/* Operate the product of a vector by a sparse matrix */
/* y = A %*% x */
void cs_vecmult(const cs *A, int nout, const double *x, double *y)
{
  int *Ap, *Ai, n;
  double *Ax, value;

  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;

  for (int j = 0; j < nout; j++)
    y[j] = 0.;

  for (int j = 0; j < n; j++)
  {
    value = 0.;
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
      value += Ax[p] * x[Ai[p]];
    y[j] = value;
  }
}

/* Operate the product of a vector by a sparse matrix */
/* y = t(A) %*% x */
void cs_tmulvec(const cs *A, int nout, const double *x, double *y)
{
  int *Ap, *Ai, n;
  double *Ax, value;

  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;

  for (int j = 0; j < nout; j++)
    y[j] = 0.;

  for (int j = 0; j < n; j++)
  {
    value = 0.;
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
      value += Ax[p] * x[Ai[p]];
    y[j] = value;
  }
}

/* Operate the product of a vector by a sparse matrix */
/* y = x %*% A */
void cs_mulvec(const cs *A, int nout, const double *x, double *y)
{
  int *Ap, *Ai, n;
  double *Ax;

  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;

  for (int j = 0; j < nout; j++)
    y[j] = 0.;

  for (int j = 0; j < n; j++)
  {
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
    {
      y[Ai[p]] += Ax[p] * x[j];
    }
  }
}

cs* cs_normalize_by_diag_and_release(cs *Q, int flag_release)
{
  cs *diag, *Qp;

  Qp = diag = nullptr;

  diag = cs_extract_diag(Q, -2);
  if (diag == nullptr) goto label_end;
  Qp = cs_prod_norm_and_release(Q, diag, flag_release);
  if (Qp == nullptr) goto label_end;

  label_end: diag = cs_spfree(diag);
  return (Qp);
}

static double func0(double x)
{
  return (x);
}
static double func1(double x)
{
  return (1. / x);
}
static double func2(double x)
{
  return (1 / sqrt(x));
}
static double func3(double x)
{
  return (sqrt(x));
}

/* Operate right-product of sparse matrix A by diagonal matrix X */
/* (entered as a vector): B = A %*% oper(X) */
/* 'oper' is Identity(0), Inverse(1), Inverse Square Root(2), Square Root(3) */

cs* cs_matvecR(const cs *A, double *x, int oper)
{
  int *Ap, *Ai, n;
  double *Ax, *Bx;
  cs *B;
  double (*oper_func)(double);

  oper_func = NULL;
  if (oper == 0)
    oper_func = func0;
  else if (oper == 1)
    oper_func = func1;
  else if (oper == 2)
    oper_func = func2;
  else if (oper == 3)
    oper_func = func3;
  else
    messageAbort("Function not found");

  B = cs_duplicate(A);
  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  Bx = B->x;

  for (int j = 0; j < n; j++)
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
      Bx[p] = Ax[p] * oper_func(x[Ai[p]]);
  return (B);
}

/* Operate left-product of sparse matrix A by diagonal matrix X */
/* (entered as a vector): B = oper(X) %*% A */
/* 'oper' is Identity(0), Inverse(1), Inverse Square Root(2), Square Root(3) */

cs* cs_matvecL(const cs *A, double *x, int oper)
{
  int *Ap, n;
  double *Ax, *Bx;
  cs *B;
  double (*oper_func)(double);

  oper_func = NULL;
  if (oper == 0)
    oper_func = func0;
  else if (oper == 1)
    oper_func = func1;
  else if (oper == 2)
    oper_func = func2;
  else if (oper == 3)
    oper_func = func3;
  else
    messageAbort("Function not found");

  B = cs_duplicate(A);
  n = A->n;
  Ap = A->p;
  Ax = A->x;
  Bx = B->x;

  for (int j = 0; j < n; j++)
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
      Bx[p] = Ax[p] * oper_func(x[j]);
  return (B);
}

/* Operate norm-product of sparse matrix A by diagonal matrix X */
/* (entered as a vector): B = oper(X) %*% A %*% oper(X) */
/* 'oper': Identity(0), Inverse(1), Inverse Square Root(2), Square Root(3) */

cs* cs_matvecnorm(const cs *A, const double *x, int oper)
{
  int *Ap, *Ai, n;
  double *Ax, *Bx;
  cs *B;
  double (*oper_func)(double);

  oper_func = NULL;
  if (oper == 0)
    oper_func = func0;
  else if (oper == 1)
    oper_func = func1;
  else if (oper == 2)
    oper_func = func2;
  else if (oper == 3)
    oper_func = func3;
  else
    messageAbort("Function not found");

  B = cs_duplicate(A);
  if (B == nullptr) return (B);
  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  Bx = B->x;

  for (int j = 0; j < n; j++)
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
      Bx[p] = Ax[p] * oper_func(x[j]) * oper_func(x[Ai[p]]);
  return (B);
}

/* Same of cs_matvecnorm ... but in place                                   */
/* 'oper': Identity(0), Inverse(1), Inverse Square Root(2), Square Root(3) */

void cs_matvecnorm_inplace(cs *A, const double *x, int oper)
{
  int *Ap, *Ai, n;
  double *Ax;
  double (*oper_func)(double);

  oper_func = NULL;
  if (oper == 0)
    oper_func = func0;
  else if (oper == 1)
    oper_func = func1;
  else if (oper == 2)
    oper_func = func2;
  else if (oper == 3)
    oper_func = func3;
  else
    messageAbort("Function not found");

  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;

  for (int j = 0; j < n; j++)
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
      Ax[p] = Ax[p] * oper_func(x[j]) * oper_func(x[Ai[p]]);
}

static void st_get_FiCo(cs *L,
                        cs *Lt,
                        int *lambda,
                        int *indUd,
                        int *indFi,
                        int *indCo)
{
  int n, nU, icol, newCo, *Ltp, *Lti, *LLp, *LLi, nFi, newFi;
  double *Ltx;

  std::map<int, std::set<int> > tab;
  std::map<int, std::set<int> >::reverse_iterator rit;
  std::map<int, std::set<int> >::iterator it;
  std::set<int> *vec;

  /* Initializations */

  Ltp = L->p;
  Ltx = L->x;
  Lti = L->i;
  LLp = Lt->p;
  LLi = Lt->i;

  n = L->n;
  nU = 0;
  for (int i = 0; i < n; i++)
    nU += indUd[i];

  for (int i = 0; i < n; i++)
  {
    it = tab.find(lambda[i]);
    if (it == tab.end()) tab[lambda[i]] = std::set<int>();
    tab[lambda[i]].insert(i);
  }

  // Update indices

  while (nU > 0)
  {
    // Find the sample mostly connected

    rit = tab.rbegin();
    newCo = *rit->second.begin();
    indCo[newCo] = 1;
    indUd[newCo] = 0;
    nU--;
    vec = &rit->second;
    vec->erase(newCo);
    if (vec->empty()) tab.erase(lambda[newCo]);

    // Add the samples to the fine set 

    nFi = 0;
    for (int p = Ltp[newCo]; p < Ltp[newCo + 1]; p++)
    {
      icol = Lti[p];
      if (Ltx[p] == 1 && indUd[icol] == 1)
      {
        indFi[nFi++] = icol;
        indCo[icol] = 0;
        indUd[icol] = 0;
        nU--;
        it = tab.find(lambda[icol]);
        if (it == tab.end())
        {
          messageAbort("%d not found for %d (%d)", lambda[icol], icol, nU);
        }
        vec = &it->second;
        vec->erase(icol);
        if (vec->empty()) tab.erase(lambda[icol]);
      }
    }

    // Update Lambda (due to newCo)

    for (int p = LLp[newCo]; p < LLp[newCo + 1]; p++)
    {
      if (indUd[LLi[p]])
      {
        it = tab.find(lambda[LLi[p]]);
        if (it == tab.end())
        {
          messageAbort("Fi : %lf not found for %d", lambda[LLi[p]], LLi[p]);
        }
        vec = &it->second;
        vec->erase(LLi[p]);
        if (vec->empty()) tab.erase(lambda[LLi[p]]);
        lambda[LLi[p]]--;
        it = tab.find(lambda[LLi[p]]);
        if (it == tab.end()) // means lambda[LLi[p]] is not in tab as a key.
        tab[lambda[LLi[p]]] = std::set<int>();
        tab[lambda[LLi[p]]].insert(LLi[p]);
      }
    }

    // Update Lambda (due to newFi)

    if (nFi > 0)
    {
      for (int k = 0; k < nFi; k++)
      {
        newFi = indFi[k];
        for (int p = LLp[newFi]; p < LLp[newFi + 1]; p++)
        {
          if (indUd[LLi[p]])
          {
            it = tab.find(lambda[LLi[p]]);
            vec = &it->second;
            vec->erase(LLi[p]);
            if (vec->empty()) tab.erase(lambda[LLi[p]]);
            lambda[LLi[p]]++;
            it = tab.find(lambda[LLi[p]]);
            if (it == tab.end()) tab[lambda[LLi[p]]] = std::set<int>();
            tab[lambda[LLi[p]]].insert(LLi[p]);
          }
        }
      }
    }
  }
}

static int st_update_neigh(int *n_arg, int indloc, int *n_tab, int *r_tab)
{
  int n;

  n = *n_arg;

  // Look for already registered index value

  for (int i = 0; i < n; i++)
  {
    if (r_tab[i] != indloc) continue;
    n_tab[i]++;
    return (0);
  }

  // New index to be registered

  if (n > MAX_NEIGH) return (1);
  r_tab[n] = indloc;
  n_tab[n] = 1;
  n++;
  *n_arg = n;
  return (0);
}

static int st_coarse_type0(cs *Q,
                           int *indUd,
                           int *indFi,
                           int *indCo,
                           cs **Lret,
                           cs **Ltret)
{
  int *lambda, *Qp, *Qi, n, ip, error;
  cs *L, *Lt, *Ltriplet;
  double *Qx, u, minu;
  double eps = 0.25;

  /* Initializations */

  error = 1;
  Qp = Qi = lambda = nullptr;
  L = Lt = Ltriplet = nullptr;
  Qx = nullptr;
  n = Q->n;
  Qp = Q->p;
  Qx = Q->x;
  Qi = Q->i;

  /* Core allocation */

  Ltriplet = cs_spalloc(0, 0, 1, 1, 1);
  if (Ltriplet == nullptr) goto label_end;
  lambda = (int*) mem_alloc(sizeof(int) * n, 0);
  if (lambda == nullptr) goto label_end;

  for (int i = 0; i < n; i++)
    lambda[i] = 0;

  /* Building Ltriplet sparse matrix */

  for (int j = 0; j < n; j++)
  {
    minu = 0.;
    for (int p = Qp[j]; p < Qp[j + 1]; p++)
    {
      u = Qx[p];
      if (u < minu) minu = u;
    }
    for (int p = Qp[j]; p < Qp[j + 1]; p++)
    {
      u = Qx[p];
      if (u < eps * minu)
      {
        ip = Qi[p];
        if (!cs_entry(Ltriplet, j, ip, 1.)) goto label_end;
        lambda[ip]++;
      }
    }
  }

  /* Convert from triplet to sparse matrix */

  L = cs_triplet(Ltriplet);
  Ltriplet = cs_spfree(Ltriplet);

  // Transpose L

  Lt = cs_transpose(L, 1);
  if (Lt == nullptr) goto label_end;

  // Transform triplet into sparse matrix

  st_get_FiCo(L, Lt, lambda, indUd, indFi, indCo);

  // Set the error return code

  error = 0;
  *Lret = L;
  *Ltret = Lt;

  label_end: lambda = (int*) mem_free((char* ) lambda);
  Ltriplet = cs_spfree(Ltriplet);
  if (error)
  {
    L = cs_spfree(L);
    Lt = cs_spfree(Lt);
  }
  return (error);
}

static int st_coarse_typen(cs* /*L*/,
                           cs* Lt,
                           int type,
                           int* indUd,
                           int* indFi,
                           int* indCo)
{
  cs *Lout, *Loutt, *Ltriplet;
  double *Lx;
  int *lambda, *Lp, *Li, n, ip, iq, error, n_neigh, nb;
  int n_neigh_tab[MAX_NEIGH], r_neigh_tab[MAX_NEIGH];

  /* Initializations */

  error = 1;
  Lp = Li = lambda = nullptr;
  Lout = Loutt = Ltriplet = nullptr;
  Lx = nullptr;
  n = Lt->n;
  Lp = Lt->p;
  Lx = Lt->x;
  Li = Lt->i;

  /* Core allocation */

  Ltriplet = cs_spalloc(0, 0, 1, 1, 1);
  if (Ltriplet == nullptr) goto label_end;
  lambda = (int*) mem_alloc(sizeof(int) * n, 0);
  if (lambda == nullptr) goto label_end;

  for (int i = 0; i < n; i++)
    lambda[i] = 0;

  /* Building Ltriplet sparse matrix */

  for (int j = 0; j < n; j++)
  {
    // Skip if 'j' is fine
    if (indCo[j] != 1) continue;
    n_neigh = 0;

    // Loop on the nodes connected to row 'j'
    for (int p = Lp[j]; p < Lp[j + 1]; p++)
    {
      // Skip the non connected node
      if (Lx[p] == 0) continue;
      // Get the column index
      ip = Li[p];

      // Skip if 'ip' is coarse
      if (indCo[ip] == 1) continue;

      // Loop on the neighbors of 'iq'
      for (int q = Lp[ip]; q < Lp[ip + 1]; q++)
      {
        // Skip the non connected node
        if (Lx[q] == 0) continue;
        iq = Li[q];
        if (iq == j) continue;
        if (indCo[iq] == 0) continue;
        if (st_update_neigh(&n_neigh, iq, n_neigh_tab, r_neigh_tab))
          goto label_end;
      }
    }

    for (int k = 0; k < n_neigh; k++)
    {
      iq = r_neigh_tab[k];
      nb = n_neigh_tab[k];
      if (nb >= type)
      {
        if (!cs_entry(Ltriplet, j, iq, 1.)) goto label_end;
        lambda[iq]++;
      }
    }
  }

  /* Convert from triplet to sparse matrix */

  Lout = cs_triplet(Ltriplet);
  Ltriplet = cs_spfree(Ltriplet);

  // Transpose

  Loutt = cs_transpose(Lout, 1);

  // Transform triplet into sparse matrix

  st_get_FiCo(Lout, Loutt, lambda, indUd, indFi, indCo);

  /* Set the error return code */

  error = 0;

  label_end: lambda = (int*) mem_free((char* ) lambda);
  Lout = cs_spfree(Lout);
  Loutt = cs_spfree(Loutt);
  Ltriplet = cs_spfree(Ltriplet);
  return (error);
}

//
// L:    List of samples strongly connected:
//    for i Qij<0 and ABS(Qij)>eps * max(Qij<0)
// type: 0 for standard coarsening
//    1 for aggressive (type 1)
//    2 for aggressive (type 2)
//
int cs_coarsening(cs *Q, int type, int **indCo_ret, cs **L_ret)
{
  int *indUd, *indCo, *indFi;
  int n, error;
  cs *L, *Lt;

  // Initializations

  error = 1;
  L = Lt = nullptr;
  indUd = indCo = indFi = nullptr;
  n = Q->n;

  // Core allocation

  indUd = (int*) mem_alloc(sizeof(int) * n, 0);
  if (indUd == nullptr) goto label_end;
  indCo = (int*) mem_alloc(sizeof(int) * n, 0);
  if (indCo == nullptr) goto label_end;
  indFi = (int*) mem_alloc(sizeof(int) * n, 0);
  if (indFi == nullptr) goto label_end;

  // Construct L (for type = 0)

  for (int i = 0; i < n; i++)
  {
    indUd[i] = 1;
    indFi[i] = 0;
    indCo[i] = 0;
  }
  if (st_coarse_type0(Q, indUd, indFi, indCo, &L, &Lt)) goto label_end;

  if (type == 0)
  {
    error = 0;
    goto label_end;
  }

  /* Construct L for aggressive coarsening */

  for (int i = 0; i < n; i++)
    indUd[i] = indCo[i];

  if (st_coarse_typen(L, Lt, type, indUd, indFi, indCo)) goto label_end;

  /* Set the error return code */

  error = 0;

  label_end: indUd = (int*) mem_free((char* ) indUd);
  indFi = (int*) mem_free((char* ) indFi);
  L = cs_spfree(L);
  if (error)
  {
    indCo = (int*) mem_free((char* ) indCo);
    Lt = cs_spfree(Lt);
  }
  *indCo_ret = indCo;
  *L_ret = Lt;
  return (error);
}

cs* cs_interpolate(cs *AA, cs *Lt, int *Co)
{
  cs *IH, *IHtriplet;
  double *u, *AAx, *Ltx, sunip, sunim, supim, supip, alpha, beta, fact, val;
  int *Nip, *Pip, *AAp, *AAi, *Ltp, *Lti, n, error, rCo, ip, *indCo, s;

  // Initialization
  error = 1;
  IH = IHtriplet = nullptr;
  u = nullptr;
  Nip = Pip = nullptr;
  AAp = AA->p;
  AAx = AA->x;
  AAi = AA->i;
  n = AA->n;
  indCo = nullptr;

  // Core allocation */
  u = (double*) mem_alloc(sizeof(double) * n, 0);
  if (u == nullptr) goto label_end;
  indCo = (int*) mem_alloc(sizeof(int) * n, 0);
  if (indCo == nullptr) goto label_end;

  Nip = (int*) mem_alloc(sizeof(int) * n, 0);
  if (Nip == nullptr) goto label_end;
  Pip = (int*) mem_alloc(sizeof(int) * n, 0);
  if (Pip == nullptr) goto label_end;

  // Transpose
  IHtriplet = cs_spalloc(0, 0, 1, 1, 1);
  if (IHtriplet == nullptr) goto label_end;
  Ltp = Lt->p;
  Ltx = Lt->x;
  Lti = Lt->i;

  s = 0;
  for (int i = 0; i < n; i++)
  {
    indCo[i] = 0;
    if (!Co[i]) continue;
    indCo[i] = s;
    s++;
  }

  // Loop on the Fine elements (not Coarse)

  for (int i = 0; i < n; i++)
  {
    if (Co[i] == 1) continue;

    // Initializations
    sunim = supim = sunip = supip = 0.;

    // Connected neighbors of 'i' ('i' excluded)
    for (int p = AAp[i]; p < AAp[i + 1]; p++)
    {
      ip = AAi[p];
      Nip[ip] = 0;
      val = u[ip] = AAx[p];
      if (val == 0.) continue;
      if (ip != i)
      {
        if (val > 0)
        {
          Nip[ip] = 1;
          sunip += val;
        }
        else
        {
          sunim += val;
        }
      }
    }

    // Strongly connected neighbors on Coarse
    for (int p = AAp[i]; p < AAp[i + 1]; p++)
    {
      ip = AAi[p];
      Pip[ip] = 0;
    }

    for (int p = Ltp[i]; p < Ltp[i + 1]; p++)
    {
      ip = Lti[p];
      val = u[ip];

      if (val == 0.) continue;
      if (Co[ip] && Ltx[p] != 0.)
      {
        if (val > 0)
        {
          Pip[ip] = 1;
          supip += val;
        }
        else
        {
          Pip[ip] = -1;
          supim += val;
        }
      }
    }

    alpha = sunim / supim;
    fact = u[i];
    if (supip > 0.)
    {
      sunip = supip = 0;
      for (int p = AAp[i]; p < AAp[i + 1]; p++)
      {
        ip = AAi[p];
        if (Nip[ip]) sunip += AAx[p];
        if (Pip[ip] > 0) supip += AAx[p];
      }
      beta = sunip / supip;
    }
    else
    {
      fact += sunip;
      beta = 0.;
    }

    // Add the entries relative to 'Pip'

    for (int p = AAp[i]; p < AAp[i + 1]; p++)
    {
      ip = AAi[p];
      if (!Co[ip]) continue;
      if (Pip[ip] != 0.)
      {
        val = (Pip[ip] > 0) ? beta : alpha;
        val *= -u[ip] / fact;
        if (!cs_entry(IHtriplet, indCo[ip], i, val)) goto label_end;
      }
    }
  }

  // Add the entries corresponding to the diagonal

  rCo = -1;
  for (int p = 0; p < n; p++)
  {
    if (!Co[p]) continue;
    rCo++;
    if (!cs_entry(IHtriplet, rCo, p, 1.)) goto label_end;
  }

  // Transform from triplet to sparse
  IH = cs_triplet(IHtriplet);
  IHtriplet = cs_spfree(IHtriplet);

  // Set the error return code
  error = 0;

  label_end: u = (double*) mem_free((char* ) u);
  indCo = (int*) mem_free((char* ) indCo);
  Nip = (int*) mem_free((char* ) Nip);
  Pip = (int*) mem_free((char* ) Pip);
  IHtriplet = cs_spfree(IHtriplet);
  if (error) IH = cs_spfree(IH);
  return (IH);
}

// Calculate the norm as follows:
// mode=1: t(B) %*% A %*% B (initial form)
// mode=2: B %*% A %*% t(B)
//
cs* cs_prod_norm(int mode, cs *A, cs *B)

{
  cs *Bt, *BtA, *BA, *Res;

  // Initializations
  Bt = BtA = BA = Res = nullptr;

  // Transpose matrix B

  Bt = cs_transpose(B, 1);
  if (Bt == nullptr) goto label_end;

  // Dispatch 

  if (mode == 1)
  {
    BtA = cs_multiply(Bt, A);
    if (BtA == nullptr) goto label_end;
    Res = cs_multiply(BtA, B);
    if (Res == nullptr) goto label_end;
  }
  else
  {
    BA = cs_multiply(B, A);
    if (BA == nullptr) goto label_end;
    Res = cs_multiply(BA, Bt);
    if (Res == nullptr) goto label_end;
  }
  label_end: Bt = cs_spfree(Bt);
  BA = cs_spfree(BA);
  BtA = cs_spfree(BtA);
  return (Res);
}

// flag_upper is 1 -> Keep upper triangular part; otherwise lower triangle
// flag_diag is 0  -> Diagonal is cleared
cs* cs_triangle(cs *A, int flag_upper, int flag_diag)
{
  cs *At;
  int *Ati, *Atp, i;
  double *Atx;

  At = cs_duplicate(A);
  Atp = At->p;
  Atx = At->x;
  Ati = At->i;

  // Loop on the columns 'j'
  for (int j = 0; j < At->n; j++)
  {
    // Loop on the non-zero part of the column
    for (int p = Atp[j]; p < Atp[j + 1]; p++)
    {
      // Row index 'i'
      i = Ati[p];

      if (flag_upper)
      {
        if (i > j) Atx[p] = 0.;
      }
      else
      {
        if (j > i) Atx[p] = 0.;
      }
      if (!flag_diag && i == j) Atx[p] = 0.;
    }
  }
  return (At);
}

int cs_lsolve_lowtri(const cs *L, const double *x, double *y)
/*
 Purpose:

 CS_LSOLVE solves LowerTri(L)*y=x
 Where LowerTri(L) is the lower triangular part of L (symmetric)
 */
{
  int j, n, p, *Lp, *Li;
  double *Lx, pivot, value;
  if (!L || !x) return (0); /* check inputs */
  n = L->n;
  Lp = L->p;
  Li = L->i;
  Lx = L->x;

  for (int i = 0; i < n; i++)
  {
    value = x[i];
    pivot = 0;
    for (p = Lp[i]; p < Lp[i + 1]; p++)
    {
      j = Li[p];
      if (j < i)
        value -= Lx[p] * y[j];
      else if (i == j) pivot = Lx[p];
    }

    y[i] = value / pivot;
  }
  return (1);
}

int cs_lsolve_uptri(const cs *L, const double *x, double *y)
/*
 Purpose:

 CS_LSOLVE solves UpperTri(L)*y=x
 Where UpperTri(L) is the upper triangular part of L (symmetric)
 */
{
  int j, n, p, *Lp, *Li;
  double *Lx, pivot, value;
  if (!L || !x) return (0); /* check inputs */
  n = L->n;
  Lp = L->p;
  Li = L->i;
  Lx = L->x;

  for (int i = n - 1; i >= 0; i--)
  {
    value = x[i];
    pivot = 0;
    for (p = Lp[i]; p < Lp[i + 1]; p++)
    {
      j = Li[p];
      if (j > i)
        value -= Lx[p] * y[j];
      else if (i == j) pivot = Lx[p];
    }

    y[i] = value / pivot;
  }
  return (1);
}

/* Operate product of a vector by triangular upper part of sparse matrix */
/* y = A %*% x */
void cs_mulvec_uptri(const cs *A,
                     int nout,
                     const double *x,
                     double *y,
                     int flag_diag)
{
  int *Ap, *Ai, i, n;
  double *Ax, value;

  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;

  for (int j = 0; j < nout; j++)
    y[j] = 0.;

  for (int j = 0; j < n; j++)
  {
    value = x[j];
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
    {
      i = Ai[p];
      if (!flag_diag && i == j) continue;
      if (j < i) continue;
      y[i] += Ax[p] * value;
    }
  }
}

/* Operate product of a vector by triangular lower part of sparse matrix */
/* y = A %*% x */
void cs_mulvec_lowtri(const cs *A,
                      int nout,
                      const double *x,
                      double *y,
                      int flag_diag)
{
  int *Ap, *Ai, i, n;
  double *Ax, value;

  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;

  for (int j = 0; j < nout; j++)
    y[j] = 0.;

  for (int j = 0; j < n; j++)
  {
    value = x[j];
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
    {
      i = Ai[p];
      if (!flag_diag && i == j) continue;
      if (j > i) continue;
      y[i] += Ax[p] * value;
    }
  }
}

cs* cs_compress(cs *A)
{
  int *cols, *rows, number;
  double *vals, value;
  cs *Q, *Qtriplet;

  cs_sparse_to_triplet(A, 0, &number, &cols, &rows, &vals);

  Q = Qtriplet = nullptr;
  Qtriplet = cs_spalloc(0, 0, 1, 1, 1);
  for (int i = 0; i < number; i++)
  {
    value = vals[i];
    if (ABS(value) < 1.e-10) continue;
    (void) cs_entry(Qtriplet, rows[i], cols[i], value);
  }
  Q = cs_triplet(Qtriplet);
  Qtriplet = cs_spfree(Qtriplet);
  return (Q);
}

void cs_keypair(const char *key, cs *A, int flag_from_1)
{
  int *rows, *cols, number;
  double *vals;
  char name[100];

  number = 0;
  cols = rows = nullptr;
  vals = nullptr;
  if (A == nullptr) return;

  cs_sparse_to_triplet(A, flag_from_1, &number, &cols, &rows, &vals);

  (void) gslSPrintf(name, "%s.cols", key);
  set_keypair_int(name, 1, number, 1, cols);
  (void) gslSPrintf(name, "%s.rows", key);
  set_keypair_int(name, 1, number, 1, rows);
  (void) gslSPrintf(name, "%s.vals", key);
  set_keypair(name, 1, number, 1, vals);

  rows = (int*) mem_free((char* ) rows);
  cols = (int*) mem_free((char* ) cols);
  vals = (double*) mem_free((char* ) vals);
}

void cs_print_file(const char *radix, int rank, cs *A)
{
  int *rows, *cols, number;
  double *vals;
  FILE *file;
  char filename[100];

  number = 0;
  cols = rows = nullptr;
  vals = nullptr;
  if (A == nullptr) return;

  if (!IFFFF(rank))
    (void) gslSPrintf(filename, "%s-%d", radix, rank);
  else
    (void) gslStrcpy(filename, radix);

  file = gslFopen(filename, "w");
  if (file == nullptr) return;

  cs_sparse_to_triplet(A, 0, &number, &cols, &rows, &vals);

  for (int i = 0; i < number; i++)
  {
    fprintf(file, "%10d %10d %20.10lf\n", cols[i], rows[i], vals[i]);
  }

  (void) fclose(file);

  rows = (int*) mem_free((char* ) rows);
  cols = (int*) mem_free((char* ) cols);
  vals = (double*) mem_free((char* ) vals);
}

/****************************************************************************/
/*!
 **  Finalize the construction of the QChol structure.
 **  Perform the Cholesky decomposition
 **
 ** \return  Error return code
 **
 ** \param[in]  verbose   Verbose flag
 ** \param[in]  QC   QChol structure to be finalized
 **
 ** \remarks In case of problem the message is issued in this function
 ** \remarks If the decomposition is already performed, nothing is done
 **
 *****************************************************************************/
int qchol_cholesky(int verbose, QChol *QC)

{
  int nmax = 8;

  /* Check that the Q matrix has already been defined */

  if (QC->Q == nullptr) return (1);

  /* Cholesky decomposition is only valid for square matric */

  if (QC->Q->m != QC->Q->n)
  {
    messerr("You wish to perform a Cholesky Decomposition of a Matrix");
    messerr("which is not square: %d x %d", QC->Q->m, QC->Q->n);
    messerr("This must be an error");
    return (1);
  }

  /* Perform the Cholesky decomposition */

  if (verbose) message("  Cholesky Decomposition... ");

  if (QC->S == nullptr)
  {
    if (verbose) message("Ordering... ");
    QC->S = cs_schol(QC->Q, 0);
    if (QC->S == nullptr)
    {
      messerr("Error in cs_schol function");
      goto label_err;
    }
  }

  if (QC->N == nullptr)
  {
    if (verbose) message("Factorization... ");
    QC->N = cs_chol(QC->Q, QC->S);
    if (QC->N == nullptr)
    {
      messerr("Error in cs_chol function");
      goto label_err;
    }
  }
  if (verbose) message("Finished\n");

  /* Optional printout */

  if (OptDbg::query(EDbg::KRIGING) || OptDbg::force())
  {
    message("Q Sparse Matrix\n");
    cs_print(QC->Q, 1);
    cs_print_range("Q", QC->Q);
  }
  return (0);

  label_err: if (verbose)
    cs_print_nice("Cholesky Decomposition of QC", QC->Q, nmax, nmax);
  QC->N = cs_nfree(QC->N);
  QC->S = cs_sfree(QC->S);
  return (1);
}

/****************************************************************************/
/*!
 **  Define the path according to the number of levels
 **
 ** \param[in]  mgs    cs_MGS structure
 ** \param[in]  nlevels    Number of coarsening levels
 ** \param[in]  path_type Type of the path
 **
 *****************************************************************************/
static void st_path_define(cs_MGS *mgs, int nlevels, int path_type)
{
  int ecr, n;
  static int niw = 6;
  static int iw[] = { 1, 1, -1, 1, -1, -1 };

  /* Dispatch */

  switch (path_type)
  {
    case 1:      // V type
      mgs->npath = 2 * nlevels;
      mgs->path = (int*) mem_alloc(sizeof(int) * mgs->npath, 1);

      ecr = 0;
      for (int i = 0; i < nlevels; i++)
        mgs->path[ecr++] = 1;
      for (int i = 0; i < nlevels; i++)
        mgs->path[ecr++] = -1;
      break;

    case 2:      // W type
      if (nlevels == 1)
      {
        mgs->npath = 2 * nlevels;
        mgs->path = (int*) mem_alloc(sizeof(int) * mgs->npath, 1);
        mgs->path[0] = 1;
        mgs->path[1] = -1;
      }
      else
      {
        mgs->npath = niw;
        mgs->path = (int*) mem_alloc(sizeof(int) * mgs->npath, 1);
        for (int i = 0; i < niw; i++)
          mgs->path[i] = iw[i];

        for (int ilevel = 2; ilevel < nlevels; ilevel++)
        {
          n = mgs->npath;
          mgs->npath = 2 + 2 * n;
          mgs->path = (int*) mem_realloc((char* ) mgs->path,
                                         sizeof(int) * mgs->npath, 1);
          for (int i = n - 1; i >= 0; i--)
          {
            mgs->path[i + 1] = mgs->path[i];
            mgs->path[i + 1 + n] = mgs->path[i];
          }
          mgs->path[0] = 1;
          mgs->path[2 * n + 1] = -1;
        }
      }
      break;

    case 3:      // F type
      mgs->npath = nlevels;
      mgs->path = (int*) mem_alloc(sizeof(int) * mgs->npath, 1);
      for (int i = 0; i < nlevels; i++)
        mgs->path[i] = 1;

      for (int ilevel = 1; ilevel < nlevels; ilevel++)
      {
        n = mgs->npath;
        mgs->npath += 2 * ilevel;
        mgs->path = (int*) mem_realloc((char* ) mgs->path,
                                       sizeof(int) * mgs->npath, 1);
        for (int i = 0; i < ilevel; i++)
        {
          mgs->path[n + i] = -1;
          mgs->path[n + i + ilevel] = 1;
        }
      }

      n = mgs->npath;
      mgs->npath += nlevels;
      mgs->path = (int*) mem_realloc((char* ) mgs->path,
                                     sizeof(int) * mgs->npath, 1);
      for (int i = 0; i < nlevels; i++)
        mgs->path[i + n] = -1;
      break;
  }
}

/****************************************************************************/
/*!
 **  Define the default parameters for the Multigrid option
 **
 ** \param[in]  mgs    cs_MGS structure
 **
 *****************************************************************************/
static void st_multigrid_set_default_params(cs_MGS *mgs)
{
  mgs->flag_cg = 1;
  mgs->type_coarse = 1;
  mgs->ngc = 50;
  mgs->nmg = 4;
  mgs->ngs = 1;
  mgs->tolcg = 1.e-07;
  mgs->tolnmg = 1.e-07;
}

/****************************************************************************/
/*!
 **  Define the parameters for the multigrid manipulation
 **
 ** \param[in]  mgs     cs_MGS structure
 ** \param[in]  flag_cg     1 for Conjugate-Gradient use; 0 otherwise
 ** \param[in]  type_coarse Type of coarsening algorithm
 ** \param[in]  ngc     Maximum number of Conjugate-Gradient iterations
 ** \param[in]  nmg     Maximum number of mutligrid iterations
 ** \param[in]  ngs     Number of Gauss-Siedel relaxation cycles
 ** \param[in]  tolcg     Tolerance for the Conjugate-Gradient algorithm
 ** \param[in]  tolnmg     Tolerance for the Multigrid algorithm
 **
 *****************************************************************************/
void cs_multigrid_params(cs_MGS *mgs,
                         int flag_cg,
                         int type_coarse,
                         int ngc,
                         int nmg,
                         int ngs,
                         double tolcg,
                         double tolnmg)
{
  if (mgs == nullptr) return;
  mgs->flag_cg = flag_cg;
  if (mgs->nlevels == 0) mgs->flag_cg = 0;
  mgs->type_coarse = type_coarse;
  mgs->ngc = ngc;
  mgs->nmg = nmg;
  mgs->ngs = ngs;
  mgs->tolcg = tolcg;
  mgs->tolnmg = tolnmg;
}

/****************************************************************************/
/*!
 **  Allocate one sub-structure for the multigrid manipulation
 **
 ** \param[in]  mode     1 for allocation; -1 for deallocation
 ** \param[in]  mg     cs_MG structure (used for deallocation)
 **
 *****************************************************************************/
static cs_MG* st_monogrid_manage(int mode, cs_MG *mg)
{
  /* Dispatch */

  if (mode > 0)
  {
    mg = (cs_MG*) mem_alloc(sizeof(cs_MG), 1);
    mg->nh = 0;
    mg->nH = 0;
    mg->sumrow = nullptr;
    mg->IhH = nullptr;
    mg->A = qchol_manage(1, NULL);
  }
  else
  {
    if (mg != nullptr)
    {
      mg->sumrow = (double*) mem_free((char* ) mg->sumrow);
      mg->IhH = cs_spfree(mg->IhH);
      mg->A = qchol_manage(-1, mg->A);
      mg = (cs_MG*) mem_free((char* ) mg);
    }
  }
  return (mg);
}

/****************************************************************************/
/*!
 **  Allocate the structure for the multigrid manipulation
 **
 ** \param[in]  mgs     cs_MGS structure
 ** \param[in]  mode     1 for allocation; -1 for deallocation
 ** \param[in]  nlevels     Number of levels of the multigrid (only for mode > 0)
 ** \param[in]  path_type   Type of the Path (1:V; 2:W, 3:F) (only for mode > 0)
 **
 *****************************************************************************/
cs_MGS* cs_multigrid_manage(cs_MGS *mgs, int mode, int nlevels, int path_type)
{
  /* Dispatch */

  if (mode > 0)
  {
    mgs = (cs_MGS*) mem_alloc(sizeof(cs_MGS), 1);
    mgs->nlevels = nlevels;
    mgs->diag = nullptr;
    st_path_define(mgs, nlevels, path_type);
    st_multigrid_set_default_params(mgs);

    mgs->mg = (cs_MG**) mem_alloc(sizeof(cs_MG*) * (1 + nlevels), 1);
    for (int i = 0; i <= nlevels; i++)
      mgs->mg[i] = st_monogrid_manage(1, NULL);
  }
  else
  {
    if (mgs == nullptr) return (mgs);
    mgs->path = (int*) mem_free((char* ) mgs->path);
    mgs->diag = (double*) mem_free((char* ) mgs->diag);
    for (int i = 0; i <= nlevels; i++)
      mgs->mg[i] = st_monogrid_manage(-1, mgs->mg[i]);
    mgs->mg = (cs_MG**) mem_free((char* ) mgs->mg);
    mgs = (cs_MGS*) mem_free((char* ) mgs);
  }
  return (mgs);
}

/****************************************************************************/
/*!
 **  Print the path
 **
 *****************************************************************************/
static void st_path_print(int nlevels, int npath, int *path)
{
  if (npath == 0) return;
  message("MultiGrid Path =");
  for (int i = 0; i < npath; i++)
    message(" %d", path[i]);
  message(" -> Number of levels = %d\n", nlevels);
}

/****************************************************************************/
/*!
 **  Print one cs_MG structure
 **
 *****************************************************************************/
static void st_mg_print(cs_MGS *mgs, int rank)
{
  cs_MG *mg;

  mg = mgs->mg[rank];
  mestitle(2, "Contents of the MG structure for level %d", rank);
  if (mg->nh > 0 && mg->nH > 0)
    message("Transition between %d and %d vertices\n", mg->nh, mg->nH);
  if (mg->IhH != NULL) cs_print_range("Range of IhH", mg->IhH);
  if (mg->A != NULL) cs_print_range("Range of A", mg->A->Q);
  if (mg->sumrow != NULL)
    print_range("Range of sumrow", mg->nh, mg->sumrow, NULL);
}

/****************************************************************************/
/*!
 **  Print the cs_MGS structure
 **
 ** \param[in] mgs     cs_MGS structure
 **
 *****************************************************************************/
void cs_multigrid_print(cs_MGS *mgs)
{
  mestitle(1, "Multigrid Levels");
  st_path_print(mgs->nlevels, mgs->npath, mgs->path);
  print_range("Range of diag", mgs->ncur, mgs->diag, NULL);
  for (int i = 0; i <= mgs->nlevels; i++)
    st_mg_print(mgs, i);
}

/****************************************************************************/
/*!
 **  Returns the number of multigrid levels
 **
 ** \param[in] mgs     cs_MGS structure
 **
 *****************************************************************************/
int cs_multigrid_get_nlevels(cs_MGS *mgs)
{
  if (mgs == nullptr)
    return (0);
  else
    return (mgs->nlevels);
}

/****************************************************************************/
/*!
 **  Normalize / Denormalize the initial solution and RHS
 **  for the multigrid kriging
 **
 ** \param[in] mgs    cs_MGS structure
 ** \param[in] mode    1 : Normalize; -1 : Denormalize
 ** \param[in] z      Solution vector
 ** \param[in] b      Right-hand side
 **
 ** \remarks: When mode>0, z and b are divided by diagonal
 ** \remarks: When mode<0, z is multiplied by diagonal
 **
 *****************************************************************************/
static void st_multigrid_scale(cs_MGS *mgs, int mode, double *z, double *b)
{

  /* Dispatch */

  if (mgs->diag == nullptr) return;
  if (mode > 0)
  {
    for (int icur = 0; icur < mgs->ncur; icur++)
    {
      b[icur] *= mgs->diag[icur];
      z[icur] /= mgs->diag[icur];
    }
  }
  else
  {
    for (int icur = 0; icur < mgs->ncur; icur++)
      z[icur] *= mgs->diag[icur];
  }
}

/****************************************************************************/
/*!
 **  Ascent step
 **
 ** \param[in]  mgs     cs_MGS structure
 ** \param[in]  level     Current level
 ** \param[in]  flag_init   1 if the output vector must be initialized
 ** \param[in]  flag_scale  1 if the output vector must be scaled
 ** \param[in]  zin     Input vector
 **
 ** \param[out] zout     Output vector
 ** \param[out] work     Working array
 **
 ** \remark Arguments 'zin' and 'zout' may coincide (if correctly dimensionned)
 **
 *****************************************************************************/
static void st_multigrid_ascent(cs_MGS *mgs,
                                int level,
                                int flag_init,
                                int flag_scale,
                                double *zin,
                                double *zout,
                                double *work)
{
  cs_MG *mg;

  if (DEBUG)
    message("Ascending from %d to %d (init=%d scale=%d)\n", level + 1, level,
            flag_init, flag_scale);
  mg = mgs->mg[level];
  cs_tmulvec(mg->IhH, mg->nh, zin, work);
  if (flag_init)
    for (int icur = 0; icur < mg->nh; icur++)
      zout[icur] = work[icur];
  else
    for (int icur = 0; icur < mg->nh; icur++)
      zout[icur] += work[icur];

  if (flag_scale) for (int icur = 0; icur < mg->nh; icur++)
    if (mg->sumrow[icur] != 0.) zout[icur] /= mg->sumrow[icur];
}

/****************************************************************************/
/*!
 **  Convert results from coarsest to final scale
 **
 ** \param[in]    mgs      cs_MGS structure
 ** \param[in,out] z      Input vector
 **
 ** \param[out]    work      Working array
 **
 ** \remark Arguments 'z' and 'work' maust be dimensioned to the finest size
 **
 *****************************************************************************/
void cs_multigrid_coarse2fine(cs_MGS *mgs, double *z, double *work)
{
  for (int ilevel = mgs->nlevels - 1; ilevel >= 0; ilevel--)
    st_multigrid_ascent(mgs, ilevel, 1, 1, z, z, work);
}

/****************************************************************************/
/*!
 **  Descent step
 **
 ** \param[in]  mgs     cs_MGS structure
 ** \param[in]  level     Current level
 ** \param[in]  zin     Input vector
 ** \param[in]  rhsin     Input R.H.S. vector
 **
 ** \param[out] rhsout     Output R.H.S. vector
 ** \param[out] work     Working array
 **
 *****************************************************************************/
static void st_multigrid_descent(cs_MGS *mgs,
                                 int level,
                                 double *zin,
                                 double *rhsin,
                                 double *rhsout,
                                 double *work)
{
  cs_MG *mg;

  if (DEBUG) message("Descending from %d to %d\n", level - 1, level);
  mg = mgs->mg[level - 1];
  cs_mulvec(mg->A->Q, mg->nh, zin, work);
  for (int icur = 0; icur < mg->nh; icur++)
    work[icur] = rhsin[icur] - work[icur];
  cs_mulvec(mg->IhH, mg->nH, work, rhsout);
}

/****************************************************************************/
/*!
 **  Inversion using Cholesky
 **
 ** \param[in]  qctt     Qchol structure
 ** \param[in,out]  xcr Current vector
 ** \param[in]  rhs     Current R.H.S. vector
 **
 ** \param[out] work     Working array
 **
 *****************************************************************************/
void cs_chol_invert(QChol *qctt, double *xcr, double *rhs, double *work)
{
  int n;

  if (DEBUG) message("Cholesky Inversion\n");
  n = qctt->Q->n;
  cs_ipvec(n, qctt->S->Pinv, rhs, work);
  cs_lsolve(qctt->N->L, work);
  cs_ltsolve(qctt->N->L, work);
  cs_pvec(n, qctt->S->Pinv, work, xcr);
}

/****************************************************************************/
/*!
 **  Simulate using Cholesky
 **
 ** \param[in]  qctt     Qchol structure
 **
 ** \param[out] simu     Simulated array
 ** \param[out] work     Working array
 **
 *****************************************************************************/
void cs_chol_simulate(QChol *qctt, double *simu, double *work)
{
  int n;

  if (DEBUG) message("Cholesky Simulation\n");
  n = qctt->Q->n;

  cs_ltsolve(qctt->N->L, work);
  cs_pvec(n, qctt->S->Pinv, work, simu);
}

/****************************************************************************/
/*!
 **  Relaxation step
 **
 ** \param[in]  mgs     cs_MGS structure
 ** \param[in]  level     Current level
 ** \param[in]  mode     1 for Descending and -1 for Ascending
 ** \param[in,out]  xcr     Current vector
 ** \param[in]  rhs     Current R.H.S. vector
 **
 ** \param[out] work     Working array
 **
 *****************************************************************************/
static void st_relaxation(cs_MGS *mgs,
                          int level,
                          int mode,
                          double *xcr,
                          double *rhs,
                          double *work)
{
  cs_MG *mg;

  mg = mgs->mg[level];

  if (mode > 0)
  {

    /* Descending case */

    if (DEBUG) message("Relaxation Descending\n");

    if (level != 0) for (int icur = 0; icur < mg->nh; icur++)
      xcr[icur] = 0.;

    for (int i = 0; i < mgs->ngs; i++)
    {
      cs_mulvec_uptri(mg->A->Q, mg->nh, xcr, work, 0);
      for (int icur = 0; icur < mg->nh; icur++)
        xcr[icur] = rhs[icur] - work[icur];
      cs_lsolve_lowtri(mg->A->Q, xcr, work);
      for (int icur = 0; icur < mg->nh; icur++)
        xcr[icur] = work[icur];
    }
  }
  else
  {

    /* Ascending case */

    if (DEBUG) message("Relaxation Ascending\n");
    for (int i = 0; i < mgs->ngs; i++)
    {
      cs_mulvec_lowtri(mg->A->Q, mg->nh, xcr, work, 0);
      for (int icur = 0; icur < mg->nh; icur++)
        xcr[icur] = rhs[icur] - work[icur];
      cs_lsolve_uptri(mg->A->Q, xcr, work);
      for (int icur = 0; icur < mg->nh; icur++)
        xcr[icur] = work[icur];
    }
  }
}

/****************************************************************************/
/*!
 **  Iterative phases for the multigrid kriging
 **
 ** \return  Error return code
 **
 ** \param[in]    mgs     cs_MGS structure
 ** \param[in]    verbose  Verbose flag
 ** \param[in,out] x     Input vector
 ** \param[in]    b     R.H.S. vector
 **
 ** \param[out] work     Working array
 **
 ** \remark 'verbose' is passed as argument as it may be different from one
 ** \remark case to the other when calling this function
 **
 *****************************************************************************/
static int st_multigrid_kriging_prec(cs_MGS *mgs,
                                     int verbose,
                                     double *x,
                                     double *b,
                                     double *work)
{
  double *xcr, *rhs, *scores, norm, delta, score;
  int error, nlevels, level, ncur, niter, mode, flag_sym, nfois;

  /* Initializations */

  error = 1;
  xcr = rhs = scores = nullptr;
  ncur = mgs->ncur;
  nlevels = mgs->nlevels;
  norm = matrix_norm(b, ncur);
  flag_sym = (int) get_keypone("MG_Flag_Symmetric", 1.);
  if (verbose)
    message("Pre-conditioning phase (Niter=%d Tol=%15.10lf)\n", mgs->nmg,
            mgs->tolnmg);

  /* Core allocation */

  xcr = (double*) mem_alloc(sizeof(double) * ncur * (1 + nlevels), 0);
  if (xcr == nullptr) goto label_end;
  rhs = (double*) mem_alloc(sizeof(double) * ncur * (1 + nlevels), 0);
  if (rhs == nullptr) goto label_end;
  scores = (double*) mem_alloc(sizeof(double) * mgs->nmg, 0);
  if (scores == nullptr) goto label_end;

  /* Initialize the internal arrays */

  for (level = 0; level <= nlevels; level++)
    for (int icur = 0; icur < ncur; icur++)
    {
      XCR(level,icur) = (level == 0) ? x[icur] : 0.;
      RHS(level,icur) = (level == 0) ? b[icur] : 0.;
    }

  /* Loop on the multigrid iterations */

  niter = 0;
  nfois = (flag_sym) ? 2 : 1;

  for (int iter = 0; iter < mgs->nmg; iter++)
  {
    level = 0;

    if (mgs->nlevels > 0)
    {
      /* Loop for symmetrization */

      for (int ifois = 0; ifois < nfois; ifois++)
      {

        /* Loop on the path levels */

        for (int k = 0; k < mgs->npath; k++)
        {
          if (level != nlevels)
          {
            mode = (k > 0) ? mgs->path[k - 1] : 1;
            if (flag_sym) mode = 2 * ifois - 1;
            st_relaxation(mgs, level, mode, &XCR(level, 0), &RHS(level, 0),
                          work);
          }
          else
            cs_chol_invert(mgs->mg[level]->A, &XCR(level, 0), &RHS(level, 0),
                           work);

          level += mgs->path[k];

          if (mgs->path[k] > 0)
            st_multigrid_descent(mgs, level, &XCR(level - 1, 0),
                                 &RHS(level - 1, 0), &RHS(level, 0), work);
          else
            st_multigrid_ascent(mgs, level, 0, 0, &XCR(level + 1, 0),
                                &XCR(level, 0), work);
        }
        mode = (!flag_sym) ? -1 : 1;
        st_relaxation(mgs, 0, mode, &XCR(level, 0), &RHS(level, 0), work);
      }
    }
    else
    {
      cs_chol_invert(mgs->mg[level]->A, &XCR(level, 0), &RHS(level, 0), work);
    }

    // Calculate the score

    score = 0.;
    cs_mulvec(mgs->mg[0]->A->Q, mgs->mg[0]->A->Q->n, &XCR(0, 0), work);
    for (int icur = 0; icur < ncur; icur++)
    {
      delta = (b[icur] - work[icur]);
      score += delta * delta;
    }
    score = sqrt(score / norm);
    scores[niter++] = score;
    if (verbose)
      message("Iteration %3d -> Score = %15.10lf\n", iter + 1, score);
    if (score < mgs->tolnmg) break;
  }

  // Save the result 

  for (int icur = 0; icur < ncur; icur++)
    x[icur] = XCR(0, icur);

  // Store the scores 

  set_keypair("Multigrid_Scores", 1, 1, niter, scores);

  // Set the error return code

  error = 0;

  label_end: xcr = (double*) mem_free((char* ) xcr);
  rhs = (double*) mem_free((char* ) rhs);
  scores = (double*) mem_free((char* ) scores);
  return (error);
}

/****************************************************************************/
/*!
 **  Conjugate Gradient algorithm for the multigrid kriging
 **
 ** \return  Error return code
 **
 ** \param[in]    mgs    cs_MGS structure
 ** \param[in]    verbose Verbose flag
 ** \param[in,out] x    Input vector
 ** \param[in]    b    R.H.S. vector
 **
 ** \param[out] work    Working array
 **
 *****************************************************************************/
static int st_multigrid_kriging_cg(cs_MGS *mgs,
                                   int verbose,
                                   double *x,
                                   double *b,
                                   double *work)
{
  cs_MG *mg;
  double *resid, *p, *z, *scores, *temp, sn, s, alpha, beta, norm, score;
  int ncur, error, niter;

  // Initializations 

  error = 1;
  ncur = mgs->ncur;
  norm = matrix_norm(b, ncur);
  resid = p = z = scores = temp = nullptr;
  if (verbose)
    message("Conjugate-Gradient Phase (Nmax=%d Tol=%15.10lf)\n", mgs->ngc,
            mgs->tolcg);

  // Core allocation

  p = (double*) mem_alloc(sizeof(double) * ncur, 0);
  if (p == nullptr) goto label_end;
  z = (double*) mem_alloc(sizeof(double) * ncur, 0);
  if (z == nullptr) goto label_end;
  resid = (double*) mem_alloc(sizeof(double) * ncur, 0);
  if (resid == nullptr) goto label_end;
  temp = (double*) mem_alloc(sizeof(double) * ncur, 0);
  if (temp == nullptr) goto label_end;
  scores = (double*) mem_alloc(sizeof(double) * mgs->ngc, 0);
  if (scores == nullptr) goto label_end;

  // Calculate initial residual

  mg = mgs->mg[0];
  cs_mulvec(mg->A->Q, mg->nh, x, work);
  for (int icur = 0; icur < mg->nh; icur++)
    resid[icur] = b[icur] - work[icur];
  for (int icur = 0; icur < ncur; icur++)
    z[icur] = 0.;
  st_multigrid_kriging_prec(mgs, 0, z, resid, work);
  matrix_product(1, ncur, 1, resid, z, &sn);
  for (int icur = 0; icur < ncur; icur++)
    p[icur] = z[icur];

  // Loop 

  niter = 0;
  for (int iter = 0; iter < mgs->ngc; iter++)
  {
    s = sn;
    cs_mulvec(mg->A->Q, mg->nh, p, temp);
    matrix_product(1, ncur, 1, p, temp, &alpha);
    alpha = s / alpha;

    for (int icur = 0; icur < ncur; icur++)
    {
      x[icur] += alpha * p[icur];
      resid[icur] -= alpha * temp[icur];
    }

    score = sqrt(matrix_norm(resid, ncur) / norm);
    scores[niter++] = score;
    if (verbose)
      message("Iteration Gradient %3d -> Score = %15.10lf\n", iter + 1, score);
    if (score < mgs->tolcg) break;

    for (int icur = 0; icur < ncur; icur++)
      z[icur] = 0.;
    st_multigrid_kriging_prec(mgs, 0, z, resid, work);

    matrix_product(1, ncur, 1, resid, z, &sn);
    matrix_product(1, ncur, 1, temp, z, &beta);
    beta *= -alpha / s;
    for (int icur = 0; icur < ncur; icur++)
      p[icur] = z[icur] + beta * p[icur];
  }

  // Store the scores 

  set_keypair("Multigrid_Gradient_Scores", 1, 1, niter, scores);

  // Set the error return code

  error = 0;

  label_end: p = (double*) mem_free((char* ) p);
  z = (double*) mem_free((char* ) z);
  resid = (double*) mem_free((char* ) resid);
  temp = (double*) mem_free((char* ) temp);
  scores = (double*) mem_free((char* ) scores);
  return (error);
}

/****************************************************************************/
/*!
 **  Perform the multigrid kriging
 **
 ** \return  Error return code
 **
 ** \param[in]  mgs     cs_MGS structure
 ** \param[in]  qctt     QChol matrix of the upper level
 ** \param[in]  verbose     Verbose flag
 ** \param[in]  x0     Input vector
 ** \param[in]  b     R.H.S. vector
 **
 ** \param[out] work     Working array
 **
 *****************************************************************************/
int cs_multigrid_process(cs_MGS *mgs,
                         QChol *qctt,
                         int verbose,
                         double *x0,
                         double *b,
                         double *work)
{
  // Perform the setup (if not already done)

  if (mgs->diag == nullptr)
  {
    if (cs_multigrid_setup(mgs, qctt, 0, verbose, NULL)) return (1);
  }
  else
  {
    if (mgs->ncur != qctt->Q->n)
      messageAbort("Check that multigrid has been setup up correctly");
  }

  // Initial normalization

  st_multigrid_scale(mgs, 1, x0, b);

  // Dispatch

  if (mgs->flag_cg)
  {

    /* Multigrid using Conjugate-gradient */

    if (st_multigrid_kriging_cg(mgs, verbose, x0, b, work)) return (1);
  }
  else
  {

    /* Iterative multigrid processing */

    if (st_multigrid_kriging_prec(mgs, verbose, x0, b, work)) return (1);
  }

  // Denormalize the results

  st_multigrid_scale(mgs, -1, x0, NULL);

  return (0);
}

/****************************************************************************/
/*!
 **  Update the selection vector
 **
 ** \param[in] ncur   Dimension of array 'sel'
 ** \param[in] sel   Array containing the current selection
 ** \param[in] indCo   Array of selected samples
 **
 ** \remarks The array 'indCo' is relative to the input selection
 **
 *****************************************************************************/
static void st_selection_update(int ncur, double *sel, int *indCo)
{
  int ecr, nval;

  ecr = nval = 0;
  for (int iech = 0; iech < ncur; iech++)
  {
    if (sel[iech] == 0) continue;
    sel[iech] = indCo[ecr++];
    nval += static_cast<int>(sel[iech]);
  }
}

/****************************************************************************/
/*!
 **  Setup the Multigrid system
 **
 ** \return  Error return code
 **
 ** \param[in] mgs    cs_MGS structure
 ** \param[in] qctt    QChol matrix of the upper level
 ** \param[in] flag_sel    If 1, an output selection is created
 ** \param[in] verbose    Verbose flag
 **
 ** \param[out] sel_arg    Vector of selection (in double as used in db)
 **
 *****************************************************************************/
int cs_multigrid_setup(cs_MGS *mgs,
                       QChol *qctt,
                       int flag_sel,
                       int verbose,
                       double **sel_arg)
{
  int *indCo, error, flag_print;
  cs *L;
  double *sel;
  cs_MG *mg, *mgold;

  // Initializations

  error = 1;
  flag_print = (int) get_keypone("MG_Flag_Print", 0.);
  indCo = nullptr;
  L = nullptr;
  sel = nullptr;
  if (verbose) mestitle(1, "Coarsening %d levels", mgs->nlevels);
  if (flag_print) cs_print_file("QTT_avant", ITEST, qctt->Q);

  // Define the size of the system

  mgs->ncur = qctt->Q->n;

  // Initialize the Selection vector 

  if (flag_sel)
  {
    sel = (double*) mem_alloc(sizeof(double) * mgs->ncur, 0);
    if (sel == nullptr) goto label_end;
    for (int i = 0; i < mgs->ncur; i++)
      sel[i] = 1;
  }

  // Normalize the initial Qctt matrix

  if (mgs->nlevels > 0)
  {
    mgs->diag = csd_extract_diag(qctt->Q, -2);
    if (mgs->diag == nullptr) return (1);
    qctt->Q = cs_normalize_by_diag_and_release(qctt->Q, 1);
  }
  if (flag_print) cs_print_file("QTT_apres", ITEST, qctt->Q);

  // Loop on the levels 

  for (int ilevel = 0; ilevel <= mgs->nlevels; ilevel++)
  {
    mg = mgs->mg[ilevel];

    // Store A

    if (ilevel == 0)
    {
      mg->A->Q = cs_duplicate(qctt->Q);
      mg->nh = mg->A->Q->n;
    }
    else
    {
      mgold = mgs->mg[ilevel - 1];
      mg->A->Q = cs_prod_norm(2, mgold->A->Q, mgold->IhH);
      if (mg->A->Q == nullptr) goto label_end;
    }
    if (flag_print) cs_print_file("A", ilevel, mg->A->Q);

    if (ilevel < mgs->nlevels)
    {
      // Extract the coarse grid

      if (cs_coarsening(mg->A->Q, mgs->type_coarse, &indCo, &L)) goto label_end;
      if (flag_print) cs_print_file("L", ilevel, L);
      if (flag_print) for (int ik = 0; ik < mg->A->Q->n; ik++)
        message("indco[%d] = %d\n", ik, indCo[ik]);

      // Interpolation 

      mg->IhH = cs_interpolate(mg->A->Q, L, indCo);
      if (mg->IhH == nullptr) goto label_end;
      if (flag_print) cs_print_file("IhH", ilevel, mg->IhH);
      L = cs_spfree(L);

      // Calculate the sum of rows per column 

      mg->sumrow = cs_col_sumrow(mg->IhH, &mg->nh, &mg->nH);
      if (mg->sumrow == nullptr) goto label_end;

      // Update the selection

      if (flag_sel) st_selection_update(mgs->ncur, sel, indCo);
      indCo = (int*) mem_free((char* ) indCo);
    }

    // Cholesky of the last level 

    if (ilevel == mgs->nlevels)
    {
      if (qchol_cholesky(verbose, mg->A)) goto label_end;
    }
  }

  // Optional printout

  if (verbose) cs_multigrid_print(mgs);

  // Set the error return code

  if (flag_sel) *sel_arg = sel;
  error = 0;

  label_end: if (error) sel = (double*) mem_free((char* ) sel);
  indCo = (int*) mem_free((char* ) indCo);
  L = cs_spfree(L);
  return (error);
}

static int st_find_color(int *Qp,
                         int *Qi,
                         double *Qx,
                         int imesh,
                         int ncolor,
                         int *colors,
                         int *temp)
{
  int i;

  /* Initializations */

  for (int j = 0; j < ncolor; j++)
    temp[j] = 0;

  /* Checks the colors of the connected nodes */

  for (int p = Qp[imesh]; p < Qp[imesh + 1]; p++)
  {
    if (ABS(Qx[p]) <= 0.) continue;
    i = Qi[p];
    if (!IFFFF(colors[i])) temp[colors[i] - 1]++;
  }

  /* Look for a free color */

  for (int j = 0; j < ncolor; j++)
  {
    if (temp[j] == 0) return (j + 1);
  }
  return (-1);
}

/****************************************************************************/
/*!
 **  Creating the color coding for mesh nodes
 **
 ** \return  Error return code
 **
 ** \param[in] Q      cs structure
 ** \param[in] start    Starting value for colors ranking (usually 0 or 1)
 **
 ** \param[out] ncols    Number of colors
 **
 ** \remarks The returned array must be freed by the calling function
 **
 *****************************************************************************/
int* cs_color_coding(cs *Q, int start, int *ncols)
{
  int *colors, *temp, error, ncolor, nmesh, next_col;
  double *Qx;
  int *Qp, *Qi;

  /* Initializations */

  error = 1;
  ncolor = *ncols = 0;
  colors = temp = nullptr;
  Qp = Q->p;
  Qi = Q->i;
  Qx = Q->x;
  nmesh = Q->n;

  /* Core allocation */

  colors = (int*) mem_alloc(sizeof(int) * nmesh, 0);
  if (colors == nullptr) goto label_end;
  temp = (int*) mem_alloc(sizeof(int) * nmesh, 0);
  if (temp == nullptr) goto label_end;
  for (int i = 0; i < nmesh; i++)
    colors[i] = ITEST;

  /* Loop on the nodes of the mesh */

  for (int j = 0; j < nmesh; j++)
  {
    next_col = st_find_color(Qp, Qi, Qx, j, ncolor, colors, temp);
    if (next_col < 0)
    {
      ncolor++;
      colors[j] = ncolor;
    }
    else
    {
      colors[j] = next_col;
    }
  }

  /* Update the colors ranking fixing the starting value */

  for (int j = 0; j < nmesh; j++)
    colors[j] = colors[j] + start - 1;

  /* Set the error return code */

  *ncols = ncolor;
  error = 0;

  label_end: temp = (int*) mem_free((char* ) temp);
  if (error) colors = (int*) mem_free((char* ) colors);
  return (colors);
}

int cs_scale(cs *A)
{
  int *Ap, *Ai, n, j, p, error;
  double *Ax, *diag;

  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;

  error = 1;
  diag = csd_extract_diag(A, 1);
  if (diag == nullptr) goto label_end;

  for (j = 0; j < n; j++)
    for (p = Ap[j]; p < Ap[j + 1]; p++)
      Ax[p] = -Ax[p] / diag[Ai[p]];

  error = 0;

  label_end: diag = (double*) mem_free((char* ) diag);
  return (error);
}

/* Get a value from the Sparse matrix */
double cs_get_value(const cs *A, int row, int col)
{
  int *Ap, *Ai;
  double *Ax;

  if (!A)
  {
    return TEST;
  }
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;

  /* Loop on the elements */

  for (int p = Ap[row]; p < Ap[row + 1]; p++)
  {
    if (Ai[p] == col) return Ax[p];
  }
  return 0.;
}

/* Set a value from the Sparse matrix */
/* This function can only update a non-zero value; otherwise nothing is done */
void cs_set_value(const cs *A, int row, int col, double value)
{
  int *Ap, *Ai;
  double *Ax;

  if (!A)
  {
    return;
  }
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;

  /* Loop on the elements */

  for (int p = Ap[row]; p < Ap[row + 1]; p++)
  {
    if (Ai[p] == col)
    {
      Ax[p] = value;
      return;
    }
  }
  if (value != 0.)
  {
    my_throw("Sparse matrix: cannot modify a zero-value by a non-zero one");
  }
}

void cs_add_cste(cs *A, double value)
{
  int *Ap, n;
  double *Ax;

  if (!A)
  {
    return;
  }
  Ap = A->p;
  Ax = A->x;
  n = A->n;

  for (int j = 0; j < n; j++)
  {
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
    {
      if (Ax[p] != 0.) Ax[p] += value;
    }
  }
}

/* Calculate the inverse of sparse matrix A and store the result in sparse B */
cs* cs_invert(const cs *A, int order, double epsilon)
{
  double *x, *b;
  cs *Bt;
  css *S;
  csn *N;
  cs *B = NULL;
  int n, ok;
  if (!A) return (0); /* check inputs */
  n = A->n;
  Bt = cs_spalloc(0, 0, 1, 1, 1);
  S = cs_schol(A, order); /* ordering and symbolic analysis */
  N = cs_chol(A, S); /* numeric Cholesky factorization */
  x = (double*) cs_malloc(n, sizeof(double));
  b = (double*) cs_malloc(n, sizeof(double));
  ok = (S && N && x && b);
  if (ok)
  {
    // Loop on all columns of the Identity matrix
    for (int icol = 0; icol < n; icol++)
    {

      // Fill the column
      for (int j = 0; j < n; j++)
        b[j] = 0.;
      b[icol] = 1.;

      // Solve the linear system and returns the result in 'x'
      cs_ipvec(n, S->Pinv, b, x); /* x = P*b */
      cs_lsolve(N->L, x); /* x = L\x */
      cs_ltsolve(N->L, x); /* x = L'\x */
      cs_pvec(n, S->Pinv, x, b); /* b = P'*x */

      // Store the elements in the output sparse matrix (triplet format)
      for (int j = 0; j < n; j++)
      {
        if (ABS(b[j]) > epsilon) cs_entry(Bt, icol, j, b[j]);
      }
    }
    B = cs_triplet(Bt);
  }
  cs_free(Bt);
  cs_free(x);
  cs_free(b);
  cs_sfree(S);
  cs_nfree(N);
  return (B);
}

Triplet csToTriplet(const cs *A, bool flagFrom1)
{
  Triplet triplet;

  if (!A) return triplet;
  int ncols = A->n;

  cs *AT = cs_transpose(A, 1);
  int nrows = 0;
  if (AT != nullptr)
  {
    nrows = AT->n;
    AT = cs_spfree(AT);
  }

  triplet.nrows = nrows;
  triplet.ncols = ncols;

  int *Ap = A->p;
  int *Ai = A->i;
  double *Ax = A->x;
  int nz = A->nz;
  if (nz >= 0) return triplet;

  int nnz = Ap[ncols];
  triplet.rows.resize(nnz);
  triplet.cols.resize(nnz);
  triplet.values.resize(nnz);

  int ecr = 0;
  for (int j = 0; j < ncols; j++)
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
    {
      triplet.cols[ecr] = (flagFrom1) ? j + 1 : j;
      triplet.rows[ecr] = (flagFrom1) ? Ai[p] + 1 : Ai[p];
      triplet.values[ecr] = Ax ? Ax[p] : 1;
      if (ABS(triplet.values[ecr]) <= 0) continue;
      ecr++;
    }

  if (ecr < nnz)
  {
    triplet.cols.resize(ecr);
    triplet.rows.resize(ecr);
    triplet.values.resize(ecr);
  }
  return triplet;
}

bool cs_isSymmetric(const cs *A, bool verbose, bool detail)
{
  int nrows, ncols, count;
  double percent;
  cs_rowcol(A, &nrows, &ncols, &count, &percent);
  if (nrows != ncols)
  {
    messerr("The sparse matrix is not square (%d x %d)", nrows, ncols);
    return false;
  }

  if (verbose) message("Testing if Matrix is Symmetric:\n");
  int numError = 0;
  for (int irow = 0; irow < nrows; irow++)
    for (int icol = irow; icol < ncols; icol++)
    {
      double aij = cs_get_value(A, irow, icol);
      double aji = cs_get_value(A, icol, irow);
      if (ABS(ABS(aij) - ABS(aji)) > EPSILON5)
      {
        if (verbose && detail)
          messerr("Element (%d,%d)=%lf is different from (%d,%d)=%lf", irow,
                  icol, aij, icol, irow, aji);
        numError++;
      }
    }
  if (verbose)
  {
    if (numError > 0)
      messerr("-> Matrix is not Symmetric");
    else
      message("-> Test successful\n");
  }
  return numError <= 0;
}

bool cs_isDiagonalDominant(cs *A, bool verbose, bool detail)
{
  int nrows, ncols, count;
  double percent;
  cs_rowcol(A, &nrows, &ncols, &count, &percent);
  if (nrows != ncols)
  {
    messerr("The sparse matrix is not square (%d x %d)", nrows, ncols);
    return false;
  }

  if (verbose) message("Testing if Matrix is Diagonal Dominant\n");
  int numError = 0;
  for (int irow = 0; irow < nrows; irow++)
  {
    double row_total = 0.;
    double row_pivot = 0.;
    for (int icol = 0; icol < ncols; icol++)
    {
      double value = cs_get_value(A, irow, icol);
      if (icol == irow)
        row_pivot = ABS(value);
      else
        row_total += ABS(value);
    }
    if (row_total > row_pivot)
    {
      if (verbose && detail)
        messerr("Error in Row (%d): Sum of abs-values=%lf Pivot=%lf", irow,
                row_total, row_pivot);
      numError++;
    }
  }
  if (verbose)
  {
    if (numError > 0)
      messerr("-> Matrix is not Diagonal Dominant");
    else
      message("-> Test successful");
  }
  return numError <= 0;
}

bool cs_isDefinitePositive(cs *A, bool verbose)
{
  int error = 1;
  css *S = nullptr;
  csn *N = nullptr;

  if (verbose) message("Testing if Matrix is Definite Positive:\n");

  S = cs_schol(A, 0);
  if (S == nullptr) goto label_end;

  N = cs_chol(A, S);
  if (N == nullptr) goto label_end;

  error = 0;

  // Free memory

  label_end: N = cs_nfree(N);
  S = cs_sfree(S);

  // Optional message

  if (verbose)
  {
    if (error)
      messerr("-> Matrix is not Definite Positive");
    else
      message("-> Test successful\n");
  }
  return error <= 0;
}

/* Convert a sparse matrix to a complete matrix
 returned as vector of double */
double* cs_toArray(const cs *A)
{
  if (!A)
  {
    message("(null)\n");
    return nullptr;
  }
  int n = A->n;
  int *Ap = A->p;
  int *Ai = A->i;
  double *Ax = A->x;

  // Core allocation

  int size = cs_get_ncell(A);
  double *mat = (double*) mem_alloc(sizeof(double) * size, 1);
  for (int i = 0; i < size; i++)
    mat[i] = 0.;

  /* Loop on the elements */

  for (int j = 0; j < n; j++)
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
      MAT(j,Ai[p]) = Ax[p];
  return mat;
}

int cs_nnz(const cs *A)
{
  if (!A) return 0;
  int n = A->n;
  int *Ap = A->p;
  return (Ap[n]);
}

/**
 * Strip off elements of the input Sparse Matrix whose absolute value is smaller than eps.
 * Return a new compressed sparse matrix.
 * Consequently, strip off the vector P[in/out] which contains the order of the samples
 * @param A       Input sparse matrix
 * @param eps     Tolerance on absolute value of the input sparse matrix elements
 * @param hypothesis Stripping hypothesis
 * @param verbose Verbose flag
 * @return A pointer to the new sparse matrix (it must be freed by calling function)
 *
 * @note Note that in the current version, the input matrix A is also modified
 */
cs* cs_strip(cs *A, double eps, int hypothesis, bool verbose)
{
  cs *Q, *Qtriplet;
  if (!A) return nullptr;
  int ntotal = cs_nnz(A);

  if (eps <= 0)
  {
    Q = cs_duplicate(A);

    if (verbose)
    {
      message("No Stripping Sparse Matrix:\n");
      message("- Dimension of the sparse matrix   = %d\n", Q->n);
      message("- Number of non-zero terms         = %d\n", ntotal);
    }

    return Q;
  }
  int Apj, Apjp1;
  double epsloc;
  bool rescale = false;

  int error = 1;
  int n = A->n;
  int *Ap = A->p;
  int *Ai = A->i;
  double *Ax = A->x;

  Q = Qtriplet = nullptr;
  Qtriplet = cs_spalloc(0, 0, 1, 1, 1);

  /* Loop on the elements */

  Apj = Apjp1 = Ap[0];
  for (int j = 0; j < n; j++)
  {
    Apj = Apjp1;
    Apjp1 = Ap[j + 1];

    // Find the value of the diagonal term

    int p0 = -1;
    for (int p = Apj; p < Apjp1 && p0 < 0; p++)
    {
      if (Ai[p] == j) p0 = p;
    }

    // Adapt the value for 'epsloc' according to the choice 'hyp'

    if (hypothesis == 1)
    {
      epsloc = eps;
    }
    else if (hypothesis == 2)
    {
      epsloc = eps * j / n;
    }
    else if (hypothesis == 3)
    {
      epsloc = eps * Ax[p0];
      rescale = false;
    }
    else
    {
      epsloc = eps * Ax[p0];
      rescale = true;
    }

    // Establish the correcting value for non-zero terms

    double ratio = 1.;
    if (rescale)
    {
      double totpos = 0.;
      double totnul = 0.;
      for (int p = Apj; p < Apjp1; p++)
      {
        if (Ai[p] == p0) continue;
        double value = Ax[p];
        if (ABS(value) < epsloc)
        {
          totnul += value * value;
        }
        else
        {
          totpos += value * value;
        }
      }
      ratio = sqrt((totpos + totnul) / totpos);
    }

    // Review all elements of the row, keeping only the ones larger than 'epsloc'

    for (int p = Apj; p < Apjp1; p++)
    {
      double value = Ax[p];
      if (ABS(value) < epsloc)
      {
        Ax[p] = 0.;
      }
      else
      {
        double coeff = (p == p0) ? 1. : ratio;
        if (!cs_entry(Qtriplet, Ai[p], j, value * coeff)) goto label_end;
      }
    }
  }

  // Transform triplet into Sparse matrix

  Q = cs_triplet(Qtriplet);
  if (Q == nullptr) goto label_end;

  // Clean the P vector

  if (verbose)
  {
    int nremain = cs_nnz(Q);
    message("Stripping Sparse Matrix:\n");
    message("- Tolerance = %lf\n", eps);
    message("- Filtering Hypothesis = %d\n", hypothesis);
    message("- Dimension of the sparse matrix    = %d\n", n);
    message("- Initial Number of non-zero values = %d\n", ntotal);
    message("- Final number of non-zero values   = %d\n", nremain);
    double reduc = 100. * (double) (ntotal - nremain) / (double) ntotal;
    message("- Reduction percentage              = %6.2lf\n", reduc);
  }

  error = 0;

  label_end: Qtriplet = cs_spfree(Qtriplet);
  if (error) Q = cs_spfree(Q);
  return (Q);
}
