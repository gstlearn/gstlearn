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
#include "LinearOp/CholeskySparseInv.hpp"
#include "Matrix/LinkMatrixSparse.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Core/SparseInv.hpp"

#include "csparse_f.h"
#include <Eigen/src/Core/Matrix.h>
#include <vector>

CholeskySparse::CholeskySparse(const MatrixSparse* mat)
  : ACholesky(mat)
  , _flagEigen(false)
  , _S(nullptr)
  , _N(nullptr)
  , _factor(nullptr)
{
  const MatrixSparse* matCS = dynamic_cast<const MatrixSparse*>(mat);
  _flagEigen                = matCS->isFlagEigen();

  (void)_prepare();
}

CholeskySparse::CholeskySparse(const CholeskySparse& m)
  : ACholesky(m)
  , _flagEigen(m._flagEigen)
  , _S(nullptr)
  , _N(nullptr)
  , _factor(nullptr)
{
  if (_flagEigen)
  {
    _S = m._S;
    _N = m._N;
  }
  else
  {
    if (m._factor != nullptr)
    {
      _factor = new Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>;
      _factor = m._factor;
    }
  }
}

CholeskySparse& CholeskySparse::operator=(const CholeskySparse& m)
{
  if (this != &m)
  {
    ACholesky::operator=(m);
    _flagEigen = m._flagEigen;
    if (_flagEigen)
    {
      _S = m._S;
      _N = m._N;
    }
    else
    {
      if (m._factor != nullptr)
      {
        _factor = new Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>;
        _factor = m._factor;
      }
    }
  }
  return *this;
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
  if (_mat == nullptr) return 1;
  int ntarget = getSize();
  vcur.resize(ntarget, 0);
  if (_flagEigen)
  {
    if (_stdevEigen(vcur)) return 1;
  }
  else
  {
    if (_stdevOld(vcur)) return 1;
  }

  if (flagStDev)
    for (int iech = 0; iech < ntarget; iech++)
      vcur[iech] = sqrt(vcur[iech]);
  return 0;
}

/****************************************************************************/
/*!
 **  Perform the calculation of the Standard Deviation of Estimation Error
 **
 ** \param[out] vcur     Output array
 **
 *****************************************************************************/
int CholeskySparse::_stdevOld(VectorDouble& vcur) const
{
  VectorDouble wz;
  VectorInt wZdiagp;
  VectorInt wLmunch;
  VectorDouble d2;
  VectorDouble diag;

  cs* Dinv    = nullptr;
  cs* LDinv   = nullptr;
  cs* TLDinv  = nullptr;
  cs* Pattern = nullptr;

  int ntarget = getSize();
  int nzmax   = 0;

  /* Pre-processing */

  d2   = csd_extract_diag_VD(_N->L, 2);
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
  wz.resize(nzmax, 0);
  wZdiagp.resize(nzmax, 0);
  wLmunch.resize(nzmax, 0);

  if (sparseinv(ntarget, LDinv->p, LDinv->i, LDinv->x, d2.data(), LDinv->p,
                LDinv->i, LDinv->x, Pattern->p, Pattern->i, Pattern->x,
                wz.data(), wZdiagp.data(), wLmunch.data()) == -1) goto label_end;

  /* Extracting the diagonal of wz */

  diag = csd_extract_diag_VD(Pattern, 1);
  cs_pvec(ntarget, _S->Pinv, diag.data(), vcur.data());

label_end:
  cs_spfree2(Dinv);
  cs_spfree2(LDinv);
  cs_spfree2(TLDinv);
  cs_spfree2(Pattern);
  return 0;
}

double CholeskySparse::computeLogDeterminant() const
{
  if (!isReady()) return TEST;
  if (_flagEigen)
  {
    double det       = 0.;
    const auto& diag = _factor->vectorD(); // Diagonal of the LDL^t decomposition (don't multiply by 2.!)
    for (int i = 0; i < _size; ++i)
      det += log(diag[i]);
    return det;
  }
  VectorDouble diag = csd_extract_diag_VD(_N->L, 1);
  double det        = 0.;
  for (int i = 0; i < (int)diag.size(); i++)
    det += log(diag[i]);
  return 2. * det;
}

int CholeskySparse::setMatrix(const MatrixSparse* mat)
{
  _mat  = mat;
  _size = mat->getNRows();
  return _prepare();
}

int CholeskySparse::_prepare() const
{
  if (_mat == nullptr) return 1;
  const MatrixSparse* matCS = dynamic_cast<const MatrixSparse*>(_mat);

  if (_flagEigen)
  {
    if (_factor != nullptr) return 0;

    _factor = new Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>;
    _factor->compute(matCS->getEigenMatrix());
    if (_factor == nullptr)
    {
      messerr("Error when computing Cholesky Decomposition");
      return 1;
    }
  }
  else
  {
    if (_S != nullptr && _N != nullptr) return 0;
    _S = cs_schol(matCS->getCS(), 0);
    if (_S == nullptr)
    {
      messerr("Error in cs_schol function");
      return 1;
    }

    _N = cs_chol(matCS->getCS(), _S);
    if (_N == nullptr)
    {
      messerr("Error in cs_chol function");
      return 1;
    }
  }
  _setReady();
  return 0;
}

int CholeskySparse::addSolveX(const constvect vecin, vect vecout) const
{
  if (!isReady()) return 1;
  if (_flagEigen)
  {
    Eigen::Map<const Eigen::VectorXd> bm(vecin.data(), vecin.size());
    Eigen::Map<Eigen::VectorXd> outm(vecout.data(), vecout.size());
    outm += _factor->solve(bm);
  }
  else
  {
    VectorDouble work(_size, 0.);
    cs_ipvec(_size, _S->Pinv, vecin.data(), work.data());
    cs_lsolve(_N->L, work.data());
    cs_ltsolve(_N->L, work.data());
    add_cs_pvec(_size, _S->Pinv, work.data(), vecout.data());
  }
  return 0;
}

int CholeskySparse::addInvLtX(const constvect vecin, vect vecout) const
{
  if (!isReady()) return 1;
  if (_flagEigen)
  {
    Eigen::VectorXd temp(vecout.size());
    std::fill(temp.data(), temp.data() + temp.size(), 0.0);
    Eigen::Map<const Eigen::VectorXd> mvecin(vecin.data(), vecin.size());
    Eigen::Map<Eigen::VectorXd> mvecout(vecout.data(), vecout.size());

    Eigen::ArrayXd Ddm = 1.0 / _factor->vectorD().array().sqrt();
    Eigen::VectorXd DW = ((mvecin.array()) * Ddm).matrix();
    Eigen::VectorXd Y  = _factor->matrixU().solve(DW);
    temp               = _factor->permutationPinv() * Y;
    mvecout += temp;
  }
  else
  {
    std::vector<double> work(vecin.size());
    work.assign(
      vecin.begin(),
      vecin.end()); // We must work on a copy of b in order to preserve constness
    cs_ltsolve(_N->L, work.data());
    add_cs_pvec(_size, _S->Pinv, work.data(), vecout.data());
  }
  return 0;
}

int CholeskySparse::addLtX(const constvect vecin, vect vecout) const
{
  if (!isReady()) return 1;
  if (_flagEigen)
  {
    Eigen::VectorXd temp(vecout.size());
    std::fill(temp.data(), temp.data() + temp.size(), 0.0);
    Eigen::Map<const Eigen::VectorXd> mvecin(vecin.data(), vecin.size());
    Eigen::Map<Eigen::VectorXd> mvecout(vecout.data(), vecout.size());

    temp               = _factor->permutationP() * mvecin;
    Eigen::VectorXd Y  = _factor->matrixU() * temp;
    Eigen::ArrayXd Ddm = _factor->vectorD().array().sqrt();
    Eigen::VectorXd DW = Y.array() * Ddm;

    mvecout += DW;
  }
  else
  {
    messerr("This option has not been programmed yet");
  }
  return 0;
}

int CholeskySparse::addLX(const constvect vecin, vect vecout) const
{
  if (!isReady()) return 1;
  if (_flagEigen)
  {
    Eigen::Map<const Eigen::VectorXd> mvecin(vecin.data(), vecin.size());
    Eigen::Map<Eigen::VectorXd> mvecout(vecout.data(), vecout.size());
    Eigen::VectorXd temp(mvecin.size());
    std::fill(temp.data(), temp.data() + temp.size(), 0.0);

    Eigen::ArrayXd Ddm = _factor->vectorD().array().sqrt();
    Eigen::VectorXd DW = mvecin.array() * Ddm;
    Eigen::VectorXd Y  = _factor->matrixL() * DW;
    temp               = _factor->permutationPinv() * Y;
    mvecout += temp;
  }
  else
  {
    messerr("This option has not been programmed yet");
  }
  return 0;
}

int CholeskySparse::addInvLX(const constvect vecin, vect vecout) const
{
  if (!isReady()) return 1;
  if (_flagEigen)
  {
    Eigen::Map<const Eigen::VectorXd> mvecin(vecin.data(), vecin.size());
    Eigen::Map<Eigen::VectorXd> mvecout(vecout.data(), vecout.size());
    Eigen::VectorXd temp(mvecin.size());
    std::fill(temp.data(), temp.data() + temp.size(), 0.0);

    temp               = _factor->permutationP() * mvecin;
    Eigen::VectorXd Y  = _factor->matrixL().solve(temp);
    Eigen::ArrayXd Ddm = 1.0 / _factor->vectorD().array().sqrt();
    Eigen::VectorXd DW = ((Y.array()) * Ddm).matrix();

    mvecout += DW;
  }
  else
  {
    messerr("This option has not been programmed yet");
  }
  return 0;
}

/**
 * @brief Compute the inverse of the 'this' matrix
 * 
 * @param vcur Storing the diagonal of the inverse matrix
 * @return int 
 *
 * @note: The method 'partial_inverse' used assumes a LTT decomposition
 * (which is not the decomposition of _factor [LDLT]). Hence a local
 * decomposition is performed again here.
 */
int CholeskySparse::_stdevEigen(VectorDouble& vcur) const
{
  Eigen::Map<Eigen::VectorXd> vcurm(vcur.data(), vcur.size());

  // Find the pointor on the initial matrix
  const auto a = dynamic_cast<const MatrixSparse*>(_mat)->getEigenMatrix();

  // Construct a SimplicialLLT matrix (instead of the LDLT stored in '_factor')
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt;
  llt.compute(a);

  Eigen::SparseMatrix Qinv = partial_inverse(llt, a);
  vcurm = Qinv.diagonal();
  return 0;
}