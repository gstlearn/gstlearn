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
#include "geoslib_old_f.h"

#include "LinearOp/Cholesky.hpp"
#include "Matrix/LinkMatrixSparse.hpp"

#include "csparse_f.h"
#include <Eigen/src/Core/Matrix.h>

Cholesky::Cholesky(const MatrixSparse* mat)
    : ALinearOp(),
      _S(nullptr),
      _N(nullptr),
      _matCS(mat)
{
  _compute();
}

Cholesky::~Cholesky()
{
  _clean();
}

void Cholesky::_clean()
{
  if (_matCS == nullptr) return;
  if (_matCS->isFlagEigen()) // TODO : We should not be aware of _matCS internal storage
  {
    // Nothing to be done
  }
  else
  {
    _S = cs_sfree2(_S);
    _S = nullptr;
    _N = cs_nfree2(_N);
    _N = nullptr;
  }
  _matCS = nullptr;
}

int Cholesky::getSize() const
{
  if (_matCS == nullptr)
    return 0;
  return _matCS->getNRows();
}

/*****************************************************************************/
/*!
**  Evaluate the product: 'outv' = MAT^{-1} * 'inv' (ALinearOp heritage)
**
** \param[in]  vecin   Array of input values
**
** \param[out] vecout  Array of output values
**
*****************************************************************************/
void Cholesky::evalInverse(const VectorDouble &vecin, VectorDouble &vecout) const
{
  if (! isValid()) return;

  solve(vecin, vecout);
}

/*****************************************************************************/
/*!
**  Operate the operation: 'outv' = MAT * 'inv' (ALinearOp heritage)
**
** \param[in]  inv       Array of input values
**
** \param[out] outv      Array of output values
**
*****************************************************************************/
int  Cholesky::_addToDest(const Eigen::VectorXd& inv,
                           Eigen::VectorXd& outv) const
{
  if (!isValid()) return 1;

  // Map Eigen Vector to VectorDouble arguments
  // TODO : VectorXd => VectorDouble = Memory copy !!
  VectorDouble einv(inv.data(), inv.data() + inv.size());
  VectorDouble eoutv(outv.size());
  // TODO : add add to dest
  _matCS->prodMatVecInPlace(einv, eoutv);
  // TODO : VectorDouble => Existing preallocated VectorXd = Memory copy !!
  outv = Eigen::Map<Eigen::VectorXd>(eoutv.data(), eoutv.size());
  return 0;
}

/****************************************************************************/
/*!
 **  Finalize the construction of the QChol structure.
 **  Perform the Cholesky decomposition
 **
 ** \remarks In case of problem the message is issued in this function
 ** \remarks If the decomposition is already performed, nothing is done
 **
 *****************************************************************************/
void Cholesky::_compute()
{
  if (_matCS == nullptr)
  {
    messerr("The argument '_matCS' must be defined");
    return;
  }
  if (_matCS->isFlagEigen())
  {
    _cholSolver.compute(_matCS->getEigenMatrix());
  }
  else
  {
    _S = cs_schol(_matCS->getCS(), 0);
    if (_S == nullptr)
    {
      messerr("Error in cs_schol function");
      _clean();
      return;
    }

    _N = cs_chol(_matCS->getCS(), _S);
    if (_N == nullptr)
    {
      messerr("Error in cs_chol function");
      _clean();
      return;
    }
  }
}

int Cholesky::solve(const Eigen::VectorXd& b, Eigen::VectorXd& x) const
{
  if (! isValid()) return 1;

  if (_matCS->isFlagEigen())
  {
    x = _cholSolver.solve(b);
  }
  else
  {
    int size = _matCS->getNRows();
    VectorDouble work(size, 0.);
    cs_ipvec(size, _S->Pinv, b.data(), work.data());
    cs_lsolve(_N->L, work.data());
    cs_ltsolve(_N->L, work.data());
    cs_pvec(size, _S->Pinv, work.data(), x.data());
  }
  return 0;
}

int Cholesky::solve(const VectorDouble& b, VectorDouble& x) const
{
  if (! isValid()) return 1;

  if (_matCS->isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> bm(b.data(), _matCS->getNCols());
    Eigen::Map<Eigen::VectorXd> xm(x.data(), _matCS->getNRows());
    xm = _cholSolver.solve(bm);
  }
  else
  {
    int size = _matCS->getNRows();
    VectorDouble work(size, 0.);
    cs_ipvec(size, _S->Pinv, b.data(), work.data());
    cs_lsolve(_N->L, work.data());
    cs_ltsolve(_N->L, work.data());
    cs_pvec(size, _S->Pinv, work.data(), x.data());
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Simulate using Cholesky
 **
 ** \param[out] b   Input Vector
 ** \param[out] x   Simulated output vector
 **
 *****************************************************************************/
int Cholesky::simulate(const VectorDouble& b, VectorDouble& x) const
{
  if (! isValid()) return 1;
  int size = _matCS->getNRows();

  if (_matCS->isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> bm(b.data(), b.size());
    Eigen::Map<Eigen::VectorXd> xm(x.data(), x.size());

    Eigen::ArrayXd Ddm = 1.0 / _cholSolver.vectorD().array().sqrt();
    Eigen::VectorXd DW = ((bm.array()) * Ddm).matrix();
    Eigen::VectorXd Y = _cholSolver.matrixU().solve(DW);
    xm = _cholSolver.permutationPinv() * Y;
  }
  else
  {
    VectorDouble work = b; // We must work on a copy of b in order to preserve constness
    cs_ltsolve(_N->L, work.data());
    cs_pvec(size, _S->Pinv, work.data(), x.data());
  }

  return 0;
}

int Cholesky::simulate(const Eigen::VectorXd& b, Eigen::VectorXd& x) const
{
  if (! isValid()) return 1;
  int size = _matCS->getNRows();

  if (_matCS->isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> bm(b.data(), b.size());
    Eigen::Map<Eigen::VectorXd> xm(x.data(), x.size());

    Eigen::ArrayXd Ddm = 1.0 / _cholSolver.vectorD().array().sqrt();
    Eigen::VectorXd DW = ((bm.array()) * Ddm).matrix();
    Eigen::VectorXd Y = _cholSolver.matrixU().solve(DW);
    xm = _cholSolver.permutationPinv() * Y;
  }
  else
  {
    Eigen::VectorXd work = b; // We must work on a copy of b in order to preserve constness
    cs_ltsolve(_N->L, work.data());
    cs_pvec(size, _S->Pinv, work.data(), x.data());
  }

  return 0;
}
/****************************************************************************/
/*!
 **  Perform the calculation of the Standard Deviation of Estimation Error
 **
 ** \param[out] vcur     Output array
 ** \param[in]  flagStDev FALSE for a variance calculation, True for StDev.
 **
 *****************************************************************************/
int Cholesky::stdev(VectorDouble& vcur, bool flagStDev) const
{
  if (! isValid()) return 1;

  if (_matCS->isFlagEigen())
  {
    /// TODO : calculate stdev when eigen
    messerr("The calculation of 'stdev' is not yet performed with Eigen Library");
    return 1;
  }
  VectorDouble z;
  VectorDouble wz;
  VectorInt wZdiagp;
  VectorInt wLmunch;
  VectorDouble d2;
  VectorDouble diag;

  cs *Dinv = nullptr;
  cs *LDinv = nullptr;
  cs *TLDinv = nullptr;
  cs *Pattern = nullptr;

  int ntarget = getSize();
  int nzmax = 0;

  /* Pre-processing */

  d2 = csd_extract_diag_VD(_N->L, 2);
  Dinv = cs_extract_diag(_N->L, -1);
  if (Dinv == nullptr) goto label_end;
  LDinv = cs_multiply(_N->L, Dinv);
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
  nzmax = Pattern->nzmax;
  z.resize(ntarget, 0);
  wz.resize(nzmax, 0);
  wZdiagp.resize(nzmax, 0);
  wLmunch.resize(nzmax, 0);

  if (sparseinv(ntarget, LDinv->p, LDinv->i, LDinv->x, d2.data(), LDinv->p,
                LDinv->i, LDinv->x, Pattern->p, Pattern->i, Pattern->x,
                wz.data(), wZdiagp.data(), wLmunch.data())
      == -1) goto label_end;

  /* Extracting the diagonal of wz */

  diag = csd_extract_diag_VD(Pattern, 1);
  cs_pvec(ntarget, _S->Pinv, diag.data(), z.data());

  if (flagStDev)
    for (int iech = 0; iech < ntarget; iech++)
      vcur[iech] = sqrt(z[iech]);
  else
    for (int iech = 0; iech < ntarget; iech++)
      vcur[iech] = z[iech];

  /* Set the error return code */

  label_end:
  cs_spfree2(Dinv);
  cs_spfree2(LDinv);
  cs_spfree2(TLDinv);
  cs_spfree2(Pattern);
  return 0;
}

double Cholesky::getLogDeterminant() const
{
  if (! isValid()) return TEST;
  if (_matCS->isFlagEigen())
  {
    double det = 0.;
    const auto& diag = _cholSolver.vectorD();
    for (int i = 0; i < _matCS->getNRows(); ++i)
      det += log(diag[i]);
    det *= 2;
    return det;
  }
  VectorDouble diag = csd_extract_diag_VD(_N->L, 1);
  double det = 0.;
  for (int i = 0; i < (int) diag.size(); i++)
    det += log(diag[i]);
  return 2. * det;
}
