#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
 
#include "Matrix/MatrixSparse.hpp"

class MatrixReplacement;
using Eigen::SparseMatrix;

// MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
template<>
struct Eigen::internal::traits<MatrixReplacement> :  public Eigen::internal::traits<Eigen::SparseMatrix<double> >
{};
 
// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement> {
public:

  // Required typedefs, constants, and method:
  typedef double Scalar;
  typedef double RealScalar;
  typedef int StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor = false
  };

  Eigen::Index rows() const { return mp_mat->getNRows(); }
  Eigen::Index cols() const { return mp_mat->getNCols(); }
 
  template<typename Rhs>
  Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
    return Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
  }
 
  // Custom API:
  MatrixReplacement() : mp_mat(0) {}
 
  void attachMyMatrix(const MatrixSparse* mat) {
    mp_mat = mat;
  }
  const MatrixSparse* my_matrix() const { return mp_mat; }
 
private:
  const MatrixSparse* mp_mat;
};

// Implementation of MatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
template<typename Rhs>
struct Eigen::internal::generic_product_impl<MatrixReplacement, Rhs, Eigen::SparseShape, Eigen::DenseShape, Eigen::GemvProduct> // GEMV stands for matrix-vector
: Eigen::internal::generic_product_impl_base<MatrixReplacement,Rhs,Eigen::internal::generic_product_impl<MatrixReplacement,Rhs> >
{
  typedef typename Product<MatrixReplacement,Rhs>::Scalar Scalar;

  template<typename Dest>
  static void scaleAndAddTo(Dest& dst, const MatrixReplacement& lhs, const Rhs& rhs, const Scalar& alpha)
  {
    // This method should implement "dst += alpha * lhs * rhs" inplace,
    // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
    assert(alpha==Scalar(1) && "scaling is not implemented");
    EIGEN_ONLY_USED_FOR_DEBUG(alpha);

    // Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
    // but let's do something fancier (and less efficient):
    for(Index i=0; i<lhs.cols(); ++i)
      dst += rhs(i) * lhs.my_matrix()->getEigenMatrix().col(i);
  }
};
 
