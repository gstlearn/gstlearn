/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                               */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "geoslib_old_f.h"

#include "../doc/dev/sandbox/CholeskyEigenCG.hpp"
#include "Basic/AException.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/Timer.hpp"
#include "Basic/VectorHelper.hpp"
#include "Core/SparseInv.hpp"
#include "Matrix/LinkMatrixSparse.hpp"

#include "csparse_f.h"

CholeskyEigenCG::CholeskyEigenCG(const MatrixSparse* mat)
  : ALinearOp()
  , _S(nullptr)
  , _N(nullptr)
  , _matCS(mat)
{
  _compute();
}

CholeskyEigenCG::~CholeskyEigenCG()
{
  _clean();
}

void CholeskyEigenCG::_clean()
{
  if (_matCS == nullptr) return;
  _matCS = nullptr;
}

int CholeskyEigenCG::getSize() const
{
  if (_matCS == nullptr)
    return 0;
  else
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
void CholeskyEigenCG::evalInverse(const VectorDouble& vecin, VectorDouble& vecout) const
{
  if (!isValid()) return;

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
void CholeskyEigenCG::_evalDirect(const VectorDouble& inv, VectorDouble& outv) const
{
  if (!isValid()) return;
  _matCS->prodMatVecInPlace(inv, outv);
}

/****************************************************************************/
/*!
 **  Finalize the construction of the QChol structure.
 **  Perform the CholeskyEigenCG decomposition
 **
 ** \remarks In case of problem the message is issued in this function
 ** \remarks If the decomposition is already performed, nothing is done
 **
 *****************************************************************************/
void CholeskyEigenCG::_compute()
{
  if (_matCS == nullptr)
  {
    messerr("The argument '_matCS' must be defined");
    return;
  }
  _cholSolver.compute(_matCS->eigenMat());
}

int CholeskyEigenCG::solve(const VectorDouble& b, VectorDouble& x) const
{
  if (!isValid()) return 1;

  Eigen::Map<const Eigen::VectorXd> bm(b.data(), _matCS->getNCols());
  Eigen::Map<Eigen::VectorXd> xm(x.data(), _matCS->getNRows());
  xm = _cholSolver.solve(bm);
  return 0;
}

/****************************************************************************/
/*!
 **  Simulate using CholeskyEigenCG
 **
 ** \param[out] b   Input Vector
 ** \param[out] x   Simulated output vector
 **
 *****************************************************************************/
int CholeskyEigenCG::simulate(const VectorDouble& b, VectorDouble& x) const
{
  if (!isValid()) return 1;
  int size = _matCS->getNRows();

  Eigen::Map<const Eigen::VectorXd> bm(b.data(), b.size());
  Eigen::Map<Eigen::VectorXd> xm(x.data(), x.size());

  Eigen::ArrayXd Ddm = 1.0 / _cholSolver.vectorD().array().sqrt();
  Eigen::VectorXd DW = ((bm.array()) * Ddm).matrix();
  Eigen::VectorXd Y  = _cholSolver.matrixU().solve(DW);
  xm                 = _cholSolver.permutationPinv() * Y;

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
int CholeskyEigenCG::stdev(VectorDouble& vcur, bool flagStDev) const
{
  if (!isValid()) return 1;

  messerr("The calculation of 'stdev' is not yet performed with Eigen Library");
  return 1;
}

double CholeskyEigenCG::getLogDeterminant() const
{
  if (!isValid()) return TEST;
  //    return log(_cholSolver.determinant()); // This should be avoided to prevent overflow for large matrix
  double det = 0.;
  auto& diag = _cholSolver.vectorD();
  for (int i = 0; i < _matCS->getNRows(); ++i)
    det += log(diag[i]);
  return 2. * det;
}
