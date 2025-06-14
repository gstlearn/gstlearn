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
// Taken from:
// Simpson, Dan. 2024. "Random C++ Part 2: Sparse Partial Inverses in Eigen."
// September 5, 2024. https://dansblog.netlify.app/posts/2024-05-08-laplace/laplace.html.
//
#pragma once

#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

typedef Eigen::SparseMatrix<double>::StorageIndex StorageIndex;

template<typename SpMat>
class MatchPattern
{
  using T = typename SpMat::value_type;
  StorageIndex* m_outer;
  StorageIndex* m_inner;
  T* m_val;
  StorageIndex m_cols;
  StorageIndex m_nnz;

public:
  MatchPattern(const SpMat& A, const SpMat& pattern)
  {
    m_cols = pattern.cols();
    m_nnz  = pattern.nonZeros();

    m_outer = new StorageIndex[m_cols + 1];
    std::copy(pattern.outerIndexPtr(), pattern.outerIndexPtr() + m_cols + 1, m_outer);
    m_inner = new StorageIndex[m_nnz];
    std::copy(pattern.innerIndexPtr(), pattern.innerIndexPtr() + m_nnz, m_inner);
    m_val = new T[m_nnz];

    T* valptr = m_val;
    for (int j = 0; j < m_cols; ++j)
    {
      typename SpMat::InnerIterator Acol(A, j);
      for (typename SpMat::InnerIterator pattern_col(pattern, j);
           pattern_col; ++pattern_col)
      {
        while (Acol && (Acol.row() < pattern_col.row()))
        {
          ++Acol;
        }
        *valptr++ = Acol.value();
        ++Acol;
      }
    }
  }

  // Specialization for rank-1 matrices A = bc^T
  MatchPattern(
    const typename Eigen::Matrix<T, 1, Eigen::Dynamic>& b,
    const typename Eigen::Matrix<T, 1, Eigen::Dynamic>& c,
    const SpMat& pattern)
  {
    m_cols  = pattern.cols();
    m_nnz   = pattern.nonZeros();
    m_outer = new StorageIndex[m_cols + 1];
    std::copy(pattern.outerIndexPtr(), pattern.outerIndexPtr() + m_cols + 1, m_outer);
    m_inner = new StorageIndex[m_nnz];
    std::copy(pattern.innerIndexPtr(), pattern.innerIndexPtr() + m_nnz, m_inner);
    m_val = new T[m_nnz];

    T* valptr = m_val;
    for (int j = 0; j < m_cols; ++j)
    {
      for (typename SpMat::InnerIterator pattern_col(pattern, j);
           pattern_col; ++pattern_col)
      {
        *valptr++ = b.coeff(pattern_col.row()) * c.coeff(j);
      }
    }
  }

  ~MatchPattern()
  {
    delete[] m_inner;
    delete[] m_outer;
    delete[] m_val;
  }

  SpMat operator()()
  {
    return Eigen::Map<SpMat>(
      m_cols,
      m_cols,
      m_nnz,
      m_outer,
      m_inner,
      m_val);
  }
};

template<typename SpChol, typename SpMat>
typename SpChol::MatrixType partial_inverse(
  const SpChol& llt,
  const SpMat& pattern)
{
  typedef typename SpMat::ReverseInnerIterator reverse_it;
  StorageIndex ncols = llt.cols();
  const SpMat& L     = llt.matrixL();
  SpMat Qinv         = L.template selfadjointView<Eigen::Lower>();

  for (int i = ncols - 1; i >= 0; --i)
  {
    reverse_it QinvcolI(Qinv, i);
    for (reverse_it LcolI_slow(L, i); LcolI_slow; --LcolI_slow)
    {
      // inner sum iterators
      reverse_it LcolI(L, i);
      reverse_it QinvcolJ(Qinv, LcolI_slow.row());

      // Initialize Qinv[j,i]
      QinvcolI.valueRef() = 0.0;

      // Inner-most sum
      while (LcolI.row() > i)
      {
        // First up, sync the iterators
        while (QinvcolJ && (LcolI.row() < QinvcolJ.row()))
        {
          --QinvcolJ;
        }
        if (QinvcolJ && (QinvcolJ.row() == LcolI.row()))
        {
          QinvcolI.valueRef() -= LcolI.value() * QinvcolJ.value();
          --QinvcolJ;
        }
        --LcolI;
      }
      // At this point LcolI is the diagonal value
      if (i == LcolI_slow.row())
      {
        QinvcolI.valueRef() += 1 / LcolI.value();
        QinvcolI.valueRef() /= LcolI.value();
      }
      else
      {
        QinvcolI.valueRef() /= LcolI.value();
        // Set Qinv[i,j] = Qinv[j,i]
        while (QinvcolJ.row() > i)
        {
          --QinvcolJ;
        }
        QinvcolJ.valueRef() = QinvcolI.value();
      }
      --QinvcolI;
    }
  }

  // Undo the permutation
  Qinv = Qinv.twistedBy(llt.permutationP().inverse());

  // Return the non-zero elements of Qinv corresponding to the non-zero
  // elements of Q
  return MatchPattern(Qinv, pattern)();
}
