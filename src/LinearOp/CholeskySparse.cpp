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
#include "Matrix/MatrixSparse.hpp"
#include <Eigen/src/Core/Matrix.h>

namespace gstlrn
{
CholeskySparse::CholeskySparse(const MatrixSparse& mat)
  : ACholesky(mat)
  , _factor(nullptr)
{
  (void)_prepare(mat);
}

CholeskySparse::CholeskySparse(const CholeskySparse& m)
  : ACholesky(m)
  , _factor(nullptr)
{
  if (m._factor != nullptr)
  {
    _factor = new Eigen::SimplicialLDLT<Sp>;
    _factor = m._factor;
  }
}

CholeskySparse& CholeskySparse::operator=(const CholeskySparse& m)
{
  if (this != &m)
  {
    ACholesky::operator=(m);
    if (m._factor != nullptr)
    {
      _factor = new Eigen::SimplicialLDLT<Sp>;
      _factor = m._factor;
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
  delete _factor;
  _factor = nullptr;
}

/****************************************************************************/
/*!
 **  Perform the calculation of the Standard Deviation of Estimation Error
 **
 ** \param[out] vcur     Output array
 ** \param[out] proj Projection to the final output dimension
 ** \param[in]  flagStDev FALSE for a variance calculation, True for StDev.
 **
 *****************************************************************************/
Id CholeskySparse::stdev(VectorDouble& vcur,
                         const MatrixSparse* proj,
                         bool flagStDev) const
{
  if (_stdev(vcur, proj)) return 1;

  if (flagStDev)
    for (Id iech = 0, ntarget = static_cast<Id>(vcur.size()); iech < ntarget; iech++)
      vcur[iech] = sqrt(vcur[iech]);
  return 0;
}

double CholeskySparse::computeLogDeterminant() const
{
  if (!isReady()) return TEST;
  double det       = 0.;
  const auto& diag = _factor->vectorD(); // Diagonal of the LDL^t decomposition (don't multiply by 2.!)
  for (Id i = 0; i < _size; ++i)
    det += log(diag[i]);
  return det;
}

Id CholeskySparse::setMatrix(const MatrixSparse& mat)
{
  _size = mat.getNRows();
  return _prepare(mat);
}

Id CholeskySparse::_prepare(const MatrixSparse& mat) const
{
  if (_factor != nullptr) return 0;

  _factor = new Eigen::SimplicialLDLT<Sp>;
  _factor->compute(mat.eigenMat());
  if (_factor == nullptr)
  {
    messerr("Error when computing Cholesky Decomposition");
    return 1;
  }
  _setReady();
  return 0;
}

Id CholeskySparse::addSolveX(const constvect vecin, vect vecout) const
{
  if (!isReady()) return 1;
  Eigen::Map<const Eigen::VectorXd> bm(vecin.data(), vecin.size());
  Eigen::Map<Eigen::VectorXd> outm(vecout.data(), vecout.size());
  outm += _factor->solve(bm);
  return 0;
}

Id CholeskySparse::addInvLtX(const constvect vecin, vect vecout) const
{
  if (!isReady()) return 1;
  Eigen::VectorXd temp(vecout.size());
  std::fill(temp.data(), temp.data() + temp.size(), 0.0);
  Eigen::Map<const Eigen::VectorXd> mvecin(vecin.data(), vecin.size());
  Eigen::Map<Eigen::VectorXd> mvecout(vecout.data(), vecout.size());

  Eigen::ArrayXd Ddm = 1.0 / _factor->vectorD().array().sqrt();
  Eigen::VectorXd DW = ((mvecin.array()) * Ddm).matrix();
  Eigen::VectorXd Y  = _factor->matrixU().solve(DW);
  temp               = _factor->permutationPinv() * Y;
  mvecout += temp;
  return 0;
}

Id CholeskySparse::addLtX(const constvect vecin, vect vecout) const
{
  if (!isReady()) return 1;
  Eigen::VectorXd temp(vecout.size());
  std::fill(temp.data(), temp.data() + temp.size(), 0.0);
  Eigen::Map<const Eigen::VectorXd> mvecin(vecin.data(), vecin.size());
  Eigen::Map<Eigen::VectorXd> mvecout(vecout.data(), vecout.size());

  temp               = _factor->permutationP() * mvecin;
  Eigen::VectorXd Y  = _factor->matrixU() * temp;
  Eigen::ArrayXd Ddm = _factor->vectorD().array().sqrt();
  Eigen::VectorXd DW = Y.array() * Ddm;

  mvecout += DW;
  return 0;
}

Id CholeskySparse::addLX(const constvect vecin, vect vecout) const
{
  if (!isReady()) return 1;
  Eigen::Map<const Eigen::VectorXd> mvecin(vecin.data(), vecin.size());
  Eigen::Map<Eigen::VectorXd> mvecout(vecout.data(), vecout.size());
  Eigen::VectorXd temp(mvecin.size());
  std::fill(temp.data(), temp.data() + temp.size(), 0.0);

  Eigen::ArrayXd Ddm = _factor->vectorD().array().sqrt();
  Eigen::VectorXd DW = mvecin.array() * Ddm;
  Eigen::VectorXd Y  = _factor->matrixL() * DW;
  temp               = _factor->permutationPinv() * Y;
  mvecout += temp;
  return 0;
}

Id CholeskySparse::addInvLX(const constvect vecin, vect vecout) const
{
  if (!isReady()) return 1;
  Eigen::Map<const Eigen::VectorXd> mvecin(vecin.data(), vecin.size());
  Eigen::Map<Eigen::VectorXd> mvecout(vecout.data(), vecout.size());
  Eigen::VectorXd temp(mvecin.size());
  std::fill(temp.data(), temp.data() + temp.size(), 0.0);

  temp               = _factor->permutationP() * mvecin;
  Eigen::VectorXd Y  = _factor->matrixL().solve(temp);
  Eigen::ArrayXd Ddm = 1.0 / _factor->vectorD().array().sqrt();
  Eigen::VectorXd DW = ((Y.array()) * Ddm).matrix();

  mvecout += DW;
  return 0;
}

/**
 * @brief Returns the diagonal of the inverse of 'this' matrix
 *
 * @param vcur Output vector
 * @param proj Projection sparse matrix
 *
 * @note: The method 'partial_inverse' used assumes a LTT decomposition
 * (which is not the decomposition of _factor [LDLT]). Hence a local
 * decomposition is performed again here.
 * This should be optimally replaced by a more clever version of
 * the original Takahashi algorithm (see sparseinv in old code)
 */
Id CholeskySparse::_stdev(VectorDouble& vcur,
                          const MatrixSparse* proj) const
{
  Eigen::Map<Eigen::VectorXd> vcurm(vcur.data(), vcur.size());
  vcurm.setZero();

  // Conversion de P en RowMajor (une seule fois)
  const auto& Pcol = proj->eigenMat(); // ℓ×k, col‑major
  using SpRowMat   = Eigen::SparseMatrix<double, Eigen::RowMajor>;
  SpRowMat P       = Pcol; // copie creuse -> row‑major ; coût négligeable si P très creuse

  // Alias ----------------------------------------------------------------------
  using SpRowMat = Eigen::SparseMatrix<double, Eigen::RowMajor>; // P (row‑major)
  using SpVec    = Eigen::SparseVector<double>;
  using DenseVec = Eigen::VectorXd;

  const Id k = P.cols();
  const Id l = P.rows();
  SpVec p_i(k);

  for (Id i = 0; i < l; ++i)
  {
    p_i.setZero();

    // 1) Extraction de la ligne i de P (row‑major) ------------------------
    for (SpRowMat::InnerIterator it(P, i); it; ++it)
      p_i.coeffRef(it.col()) = it.value();

    // 2) Résolution :  Q y = p_iᵀ   --------------------------------------
    DenseVec y = _factor->solve(p_i);

    // 3) Contribution diagonale : p_i · y  -------------------------------
    vcurm(i) = p_i.dot(y);
  }
  return 0;
}
} // namespace gstlrn