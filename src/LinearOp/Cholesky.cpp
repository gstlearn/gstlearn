/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "geoslib_old_f.h"
#include "geoslib_f.h"

#include "LinearOp/Cholesky.hpp"
#include "LinearOp/Identity.hpp"
#include "Basic/AException.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/VectorHelper.hpp"
#include "Matrix/LinkMatrixSparse.hpp"

#include <iostream>

// External library /// TODO : Dependency to csparse to be removed
#include "csparse_d.h"
#include "csparse_f.h"

Cholesky::Cholesky(const cs *mat, bool flagDecompose)
    : ALinearOp(),
      _mat(nullptr),
      _matS(nullptr),
      _matN(nullptr),
      _work()
{
  reset(mat, flagDecompose);
}

Cholesky::Cholesky(const Cholesky &m)
    : ALinearOp(m),
      _mat(nullptr),
      _matS(nullptr),
      _matN(nullptr),
      _work()
{
  reset(m._mat, m._isDecomposed());
}

Cholesky& Cholesky::operator=(const Cholesky &m)
{
  if (this != &m)
  {
    ALinearOp::operator =(m);
    reset(m._mat, m._isDecomposed());
  }
  return *this;
}

Cholesky::~Cholesky()
{
  _clean();
}

void Cholesky::_clean()
{
  _matS = cs_sfree(_matS);
  _matN = cs_nfree(_matN);
}

int Cholesky::reset(const cs* mat, bool flagDecompose)
{
  if (mat == nullptr) return 0;

  // Clear the already existing contents
  _clean();

  // Check that the input matrix is square and symmetric
  if (! cs_isSymmetric(mat)) return 1;

  // Duplicate the sparse matrix
  _mat = mat;

  // Perform the Cholesky decomposition
  if (flagDecompose) _decompose();

  return 0;
}

int Cholesky::getSize() const
{
  if (! _isDefined())
    return 0;
  else
    return _mat->n;
}

/*****************************************************************************/
/*!
**  Evaluate the product: 'outv' = MAT^{-1} * 'inv'
**
** \param[in]  inv     Array of input values
**
** \param[out] outv    Array of output values
**
*****************************************************************************/
void Cholesky::evalInverse(const VectorDouble &inv, VectorDouble &outv) const
{
  if (!_isDefined())
  {
    messerr("The sparse matrix should be defined beforehand");
    return;
  }

  // Perform Cholesky factorization (if needed)
  _decompose();

  int n = getSize();
  for (int i = 0; i < n; i++) _work[i] = 0.;
  cs_ipvec(n, _matS->Pinv, inv.data(), _work.data());
  cs_lsolve(_matN->L, _work.data());
  cs_ltsolve(_matN->L, _work.data());
  cs_pvec(n, _matS->Pinv, _work.data(), outv.data());
}

/*****************************************************************************/
/*!
**  Operate the operation: 'outv' = MAT * 'inv'
**
** \param[in]  inv       Array of input values
**
** \param[out] outv      Array of output values
**
*****************************************************************************/
void Cholesky::_evalDirect(const VectorDouble &inv, VectorDouble &outv) const
{
  if (! _isDefined())
  {
    messerr("The sparse matrix should be defined beforehand");
    return;
  }
  int n = getSize();
  cs_vecmult(_mat, n, inv.data(), outv.data());
}

void Cholesky::printout(const char *title, bool verbose) const
{
  if (! _isDefined()) return;

  if (title != NULL) message("%s\n", title);

  int nrows = 0;
  int ncols = 0;
  int count = 0;
  double percent = 0.;
  cs_rowcol(_mat, &nrows, &ncols, &count, &percent);
  message("- Nrows(%d) x Ncols(%d) - Non-zeros(%d) [%6.2lf (percent)]", nrows,
          ncols, count, percent);
  if (_matS != nullptr || _matN != nullptr) message(" (Cholesky)");
  message("\n");

  if (verbose)
    cs_print_nice("Symmetric Matrix", _mat);

  message("\n");
}

/****************************************************************************/
/*!
 **  Finalize the construction of the QChol structure.
 **  Perform the Cholesky decomposition
 **
 ** \param[in]  verbose   Verbose flag
 **
 ** \remarks In case of problem the message is issued in this function
 ** \remarks If the decomposition is already performed, nothing is done
 **
 *****************************************************************************/
void Cholesky::_decompose(bool verbose) const
{
  if (!_isDefined()) return;
  if (_isDecomposed()) return;

  /* Perform the Cholesky decomposition */

  if (verbose) message("  Cholesky Decomposition... ");

  if (verbose) message("Ordering... ");
  _matS = cs_schol(_mat, 0);
  if (_matS == nullptr)
  {
    messerr("Error in cs_schol function");
    _matN = cs_nfree(_matN);
    return;
  }

  if (verbose) message("Factorization... ");
  _matN = cs_chol(_mat, _matS);
  if (_matN == nullptr)
  {
    messerr("Error in cs_chol function");
    _matS = cs_sfree(_matS);
    return;
  }

  _work.resize(getSize());

  if (verbose) message("Finished\n");
  return;
}

/****************************************************************************/
/*!
 **  Simulate using Cholesky
 **
 ** \param[out] inv      input Vector
 ** \param[out] outv     Simulated output vector
 **
 *****************************************************************************/
void Cholesky::simulate(VectorDouble& inv, VectorDouble& outv)
{
  if (!_isDefined())
  {
    messerr("The sparse matrix should be defined beforehand");
    return;
  }

  // Perform Cholesky factorization (if needed)
  _decompose();

  int n = getSize();
  cs_ltsolve(_matN->L, inv.data());
  cs_pvec(n, _matS->Pinv, inv.data(), outv.data());
}


/****************************************************************************/
/*!
 **  Perform the calculation of the Standard Deviation of Estimation Error
 **
 ** \param[out] vcur     Output array
 ** \param[in]  flagStDev FALSE for a variance calculation, True for StDev.
 **
 *****************************************************************************/
void Cholesky::stdev(VectorDouble& vcur, bool flagStDev)
{
  if (!_isDefined())
  {
    messerr("The sparse matrix should be defined beforehand");
    return;
  }

  // Perform Cholesky factorization (if needed)
   _decompose();

  VectorDouble z;
  VectorDouble wz;
  VectorInt wZdiagp;
  VectorInt wLmunch;
  VectorDouble d2;
  VectorDouble diag;

  cs* Dinv = nullptr;
  cs* LDinv = nullptr;
  cs* TLDinv = nullptr;
  cs* Pattern = nullptr;

  int ntarget = getSize();
  int nzmax = 0;

  /* Pre-processing */

  d2 = csd_extract_diag_VD(_matN->L, 2);
  Dinv = cs_extract_diag(_matN->L, -1);
  if (Dinv == nullptr) goto label_end;
  LDinv = cs_multiply(_matN->L, Dinv);
  if (LDinv == nullptr) goto label_end;
  TLDinv = cs_transpose(LDinv, 1);
  if (TLDinv == nullptr) goto label_end;
  Pattern = cs_add(LDinv, TLDinv, 1, 1);
  if (Pattern == nullptr) goto label_end;
  if (cs_sort_i(Pattern)) goto label_end;
  if (cs_sort_i(LDinv)) goto label_end;

  /* Core allocation */

  nzmax = Pattern->nzmax;
  z.resize(ntarget, 0);
  wz.resize(nzmax, 0);
  wZdiagp.resize(nzmax, 0);
  wLmunch.resize(nzmax, 0);

  if (sparseinv(ntarget, LDinv->p, LDinv->i, LDinv->x, d2.data(), LDinv->p, LDinv->i,
                LDinv->x, Pattern->p, Pattern->i, Pattern->x, wz.data(),
                wZdiagp.data(), wLmunch.data()) == -1) goto label_end;

  /* Extracting the diagonal of wz */

  diag = csd_extract_diag_VD(Pattern, 1);
  cs_pvec(ntarget, _matS->Pinv, diag.data(), z.data());

  if (flagStDev)
    for (int iech = 0; iech < ntarget; iech++)
      vcur[iech] = sqrt(z[iech]);
  else
    for (int iech = 0; iech < ntarget; iech++)
      vcur[iech] = z[iech];

  /* Set the error return code */

  label_end:
  Dinv = cs_spfree(Dinv);
  LDinv = cs_spfree(LDinv);
  TLDinv = cs_spfree(TLDinv);
  Pattern = cs_spfree(Pattern);
}

double Cholesky::computeLogDet() const
{
  if (!_isDefined())
  {
    messerr("The sparse matrix should be defined beforehand");
    return TEST;
  }

  // Perform Cholesky factorization (if needed)
   _decompose();

  VectorDouble diag = csd_extract_diag_VD(_matN->L, 1);
  double det = 0.;
  for (int i = 0; i < (int) diag.size(); i++)
    det += log(diag[i]);

  return 2. * det;
}
