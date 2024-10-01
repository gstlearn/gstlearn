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
#include "LinearOp/CholeskySparse.hpp"
#include "Matrix/LinkMatrixSparse.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Core/SparseInv.hpp"

#include "csparse_f.h"
#include <Eigen/src/Core/Matrix.h>
#include <vector>

CholeskySparse::CholeskySparse(const MatrixSparse* mat, bool inverse)
  : ACholesky(mat, inverse)
  , _flagEigen(false)
  , _S(nullptr)
  , _N(nullptr)
  , _factor(nullptr)
{
  const MatrixSparse* matCS = dynamic_cast<const MatrixSparse*>(mat);
  _flagEigen                = matCS->isFlagEigen();

  _prepare();
}

CholeskySparse::~CholeskySparse()
{
  _clean();
}

void CholeskySparse::_clean()
{
  if (_flagEigen)
  {
    delete _factor;
    _factor = nullptr;
  }
  else
  {
    _S = cs_sfree2(_S);
    _S = nullptr;
    _N = cs_nfree2(_N);
    _N = nullptr;
  }
}

int CholeskySparse::_addSolveX(const constvect b, vect x) const
{
  if (_flagEigen)
  { 
    Eigen::Map<const Eigen::VectorXd> bm(b.data(),b.size());
    Eigen::Map<Eigen::VectorXd> outm(x.data(),x.size());
    outm += _factor->solve(bm);
  }
  else
  {
    VectorDouble work(_size, 0.);
    cs_ipvec(_size, _S->Pinv, b.data(), work.data());
    cs_lsolve(_N->L, work.data());
    cs_ltsolve(_N->L, work.data());
    add_cs_pvec(_size, _S->Pinv, work.data(), x.data());
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
int CholeskySparse::stdev(VectorDouble& vcur, bool flagStDev) const
{
  if (_flagEigen)
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

double CholeskySparse::computeLogDeterminant() const
{
  if (_flagEigen)
  {
    double det = 0.;
    const auto& diag = _factor->vectorD(); //Diagonal of the LDL^t decomposition (don't multiply by 2.!)
    for (int i = 0; i < _size; ++i)
      det += log(diag[i]);
    return det;
  }
  VectorDouble diag = csd_extract_diag_VD(_N->L, 1);
  double det = 0.;
  for (int i = 0; i < (int) diag.size(); i++)
    det += log(diag[i]);
  return 2. * det;
}

void CholeskySparse::_prepare() const
{
  const MatrixSparse* matCS = dynamic_cast<const MatrixSparse*>(_mat);
  if (_factor == nullptr)
  {
    if (_flagEigen)
    {
      _factor = new Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>;
      _factor->compute(matCS->getEigenMatrix());
    }
    else
    {
      _S = cs_schol(matCS->getCS(), 0);
      if (_S == nullptr)
      {
        messerr("Error in cs_schol function");
        return;
      }

      _N = cs_chol(matCS->getCS(), _S);
      if (_N == nullptr)
      {
        messerr("Error in cs_chol function");
        return;
      }
    }
  }
}

int CholeskySparse::_addInvLtX(const constvect b, vect x) const
{
  if (_flagEigen)
  {
    Eigen::VectorXd temp(x.size());
    std::fill(temp.data(), temp.data() + temp.size(), 0.0);
    Eigen::Map<const Eigen::VectorXd> bm(b.data(), b.size());
    Eigen::Map<Eigen::VectorXd> xm(x.data(), x.size());

    Eigen::ArrayXd Ddm = 1.0 / _factor->vectorD().array().sqrt();
    Eigen::VectorXd DW = ((bm.array()) * Ddm).matrix();
    Eigen::VectorXd Y  = _factor->matrixU().solve(DW);
    temp               = _factor->permutationPinv() * Y;
    xm += temp;
  }
  else
  {
    std::vector<double> work(b.size());
    work.assign(
      b.begin(),
      b.end()); // We must work on a copy of b in order to preserve constness
    cs_ltsolve(_N->L, work.data());
    add_cs_pvec(_size, _S->Pinv, work.data(), x.data());
  }
  return 0;
}

int CholeskySparse::_addLX(const constvect inv, vect outv) const
{
  if (_flagEigen)
  {
    Eigen::Map<const Eigen::VectorXd> invm(inv.data(), inv.size());
    Eigen::Map<Eigen::VectorXd> outvm(outv.data(), outv.size());
    Eigen::VectorXd temp(invm.size());
    std::fill(temp.data(), temp.data() + temp.size(), 0.0);

    Eigen::ArrayXd Ddm = _factor->vectorD().array().sqrt();
    Eigen::VectorXd DW = invm.array() * Ddm;
    Eigen::VectorXd Y  = _factor->matrixL() * DW;
    temp               = _factor->permutationPinv() * Y;
    outvm += temp;
  }
  else
  {
    messerr("This option has not been programmed yet");
  }
  return 0;
}
