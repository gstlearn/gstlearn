#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
 
//#include "LinearOp/ALinearOp.hpp"
#include "Matrix/MatrixSparse.hpp"

class ALinearOpEigenCG;
using Eigen::SparseMatrix;
 
namespace Eigen {
namespace internal {
  // ALinearOpEigenCG looks-like a SparseMatrix, so let's inherits its traits:
  template<>
  struct traits<ALinearOpEigenCG> :  public Eigen::internal::traits<Eigen::SparseMatrix<double> >
  {};
}
}
 
// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class ALinearOpEigenCG : public Eigen::EigenBase<ALinearOpEigenCG> {
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
/*
  Index rows() const { return mp_mat->rows(); }
  Index cols() const { return mp_mat->cols(); }
 
  template<typename Rhs>
  Eigen::Product<ALinearOpEigenCG,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
    return Eigen::Product<ALinearOpEigenCG,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
  }
 
  // Custom API:
  ALinearOpEigenCG() : mp_mat(0) {}
 
  void attachMyMatrix(const SparseMatrix<double> &mat) {
    mp_mat = &mat;
  }
  const SparseMatrix<double> my_matrix() const { return *mp_mat; }
 
private:
  const SparseMatrix<double> *mp_mat;
*/
///*
  Index rows() const { return mp_mat->getNRows(); }
  Index cols() const { return mp_mat->getNCols(); }
 
  template<typename Rhs>
  Eigen::Product<ALinearOpEigenCG,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
    return Eigen::Product<ALinearOpEigenCG,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
  }
 
  // Custom API:
  ALinearOpEigenCG() : mp_mat(0) {}
 
  void attachMyMatrix(const MatrixSparse* mat) {
    mp_mat = mat;
  }
  const MatrixSparse* my_matrix() const { return mp_mat; }
 
private:
  const MatrixSparse* mp_mat;
//*/
};
 
 
// Implementation of ALinearOpEigenCG * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
namespace internal {
 
  template<typename Rhs>
  struct generic_product_impl<ALinearOpEigenCG, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
  : generic_product_impl_base<ALinearOpEigenCG,Rhs,generic_product_impl<ALinearOpEigenCG,Rhs> >
  {
    typedef typename Product<ALinearOpEigenCG,Rhs>::Scalar Scalar;
 
    template<typename Dest>
    static void scaleAndAddTo(Dest& dst, const ALinearOpEigenCG& lhs, const Rhs& rhs, const Scalar& alpha)
    {
      // This method should implement "dst += alpha * lhs * rhs" inplace,
      // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
      assert(alpha==Scalar(1) && "scaling is not implemented");
      EIGEN_ONLY_USED_FOR_DEBUG(alpha);
 
      // Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
      // but let's do something fancier (and less efficient):
      for(Index i=0; i<lhs.cols(); ++i)
        dst += rhs(i) * lhs.my_matrix()->getEigenMatrix().col(i);
        //dst += rhs(i) * lhs.my_matrix().col(i);
    }
  };
 
}
}
