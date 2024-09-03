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
#include "Polynomials/ClassicalPolynomial.hpp"
#include "Basic/AStringable.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include <Eigen/src/Core/Matrix.h>

ClassicalPolynomial::ClassicalPolynomial()
{
  // TODO Auto-generated constructor stub
}


ClassicalPolynomial::ClassicalPolynomial(const VectorDouble& coeffs)
{
  init(coeffs);
}

ClassicalPolynomial::~ClassicalPolynomial()
{
  // TODO Auto-generated destructor stub
}

double ClassicalPolynomial::eval(double x) const
{
  double result=_coeffs.back();
  for (int i = static_cast<int> (_coeffs.size())-2; i >= 0; i--)
  {
    result *= x;
    result +=_coeffs[i];
  }
  return result;
}

// Horner scheme starting from the lowest degree
// (since it adds the result to the input vector, the classical scheme can t be used)
#ifndef SWIG
void ClassicalPolynomial::evalOpCumul(MatrixSparse* Op, const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const
{
  int n = static_cast<int> (inv.size());
  Eigen::VectorXd work(n);
  Eigen::VectorXd work2(n);
  Eigen::VectorXd *swap1,*swap2,*swap3;

  swap1 = &work;
  swap2 = &work2;

  for (int i = 0; i<n ; i++)
  {
    outv[i] += _coeffs[0] * inv[i];
  }

  Op->prodMatVecInPlace(inv, *swap1);

  for (int j = 1; j < (int) _coeffs.size(); j++)
  {
    for (int i = 0; i < n; i++)
    {
      outv[i] += _coeffs[j] * (*swap1)[i];
    }

    if (j < (int) _coeffs.size() - 1)
    {
      Op->prodMatVecInPlace(*swap1, *swap2);
      swap3 = swap1;
      swap1 = swap2;
      swap2 = swap3;
    }
  }
}

void ClassicalPolynomial::addEvalOp(ALinearOp* Op,const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const
{
  int n = static_cast<int> (inv.size());

  Eigen::VectorXd *swap1,*swap2,*swap3;

  if (_work.size()!= n)
  {
    _work.resize(n);
  }
  if (_work2.size()!= n)
  {
    _work2.resize(n);
  }

  swap1 = &_work;
  swap2 = &_work2;

  for (int i = 0; i<n ; i++)
  {
    outv[i] += _coeffs[0] * inv[i];
  }

  Op->evalDirect(inv, *swap1);

  for (int j = 1; j < (int) _coeffs.size(); j++)
  {
    for (int i = 0; i < n; i++)
    {
      outv[i] += _coeffs[j] * (*swap1)[i];
    }

    if (j < (int) _coeffs.size() - 1)
    {
      Op->evalDirect(*swap1, *swap2);
      swap3 = swap1;
      swap1 = swap2;
      swap2 = swap3;
    }
  }
}

// Classical Hörner scheme starting from the highest degree
void ClassicalPolynomial::evalOp(MatrixSparse* Op,
                                 const Eigen::VectorXd& inv,
                                 Eigen::VectorXd& outv) const
{
  int n = static_cast<int>(inv.size());
  Eigen::VectorXd work(n);

  for (int i = 0; i < n; i++)
    outv[i] = _coeffs.back() * inv[i];

  for (int j = static_cast<int>(_coeffs.size()) - 2; j >= 0; j--)
  {
    Op->prodMatVecInPlace(outv, work);
    for (int i = 0; i < n; i++)
    {
      outv[i] = _coeffs[j] * inv[i] + work[i];
    }
  }
}

// Classical Hörner scheme starting from the highest degree
void ClassicalPolynomial::evalOpTraining(MatrixSparse *Op,
                                         const Eigen::VectorXd &inv,
                                         std::vector<Eigen::VectorXd> &store,
                                         Eigen::VectorXd &work) const
{
  int n = static_cast<int>(inv.size());

  if (work.size() == 0)
  {
    work.resize(n);
  }

  for (int i = 0; i < n; i++)
  {
    store[_coeffs.size() - 1][i] = _coeffs.back() * inv[i];
  }

  for (int j = (int) _coeffs.size() - 2; j >= 0; j--)
  {
    Op->prodMatVecInPlace(store[j + 1], work);
    for (int i = 0; i < n; i++)
    {
      store[j][i] = _coeffs[j] * inv[i] + work[i];
    }
  }
}
#endif
void ClassicalPolynomial::evalDerivOp(ShiftOpCs* shiftOp,
                                      const Eigen::VectorXd& inv,
                                      Eigen::VectorXd& outv,
                                      int iapex,
                                      int igparam)const
{
  int n = static_cast<int> (inv.size());
  int degree = static_cast<int> (getCoeffs().size());
  ClassicalPolynomial* polycur = this->clone();
  Eigen::VectorXd work(n);
  Eigen::VectorXd work2(n);
  Eigen::VectorXd *swap1,*swap2,*swap3;
  MatrixSparse* Op = shiftOp->getS();
  MatrixSparse* derivOp = shiftOp->getSGrad(iapex,igparam);

  swap1 = &work;
  swap2 = &work2;

  VectorEigen::fill(outv, 0.);

  for(int i = 0 ; i< n; i++)
  {
    work[i] = inv[i];
  }

  auto coeffsCur = polycur->getCoeffs();

  for(int i = 0; i < degree - 1 ;i++)
  {
    derivOp->prodMatVecInPlace(*swap1,*swap2);
    coeffsCur.erase(coeffsCur.begin());
    polycur->init(coeffsCur);
    polycur->evalOpCumul(Op,*swap2,outv);

    if(i<degree-2) // to avoid useless and time consuming computation since it prepares next iteration
    {
      Op->prodMatVecInPlace(*swap1,*swap2);
      swap3 = swap1;
      swap1 = swap2;
      swap2 = swap3;
    }
   }
   delete polycur;
}

void ClassicalPolynomial::evalDerivOp(ShiftOpCs* /*shiftOp*/,
                                      const VectorDouble& /*inv*/,
                                      VectorDouble& /*outv*/,
                                      int /*iapex*/,
                                      int /*igparam*/)
{
  //TODO Call the Eigen::VectorXd function
  messerr("evalDerivOp is not implemented for vectorsDouble");
}

//void ClassicalPolynomial::evalDerivOpOptim(ShiftOpCs* shiftOp,
//                                           const Eigen::VectorXd& in1,
//                                           Eigen::VectorXd& in2,
//                                           Eigen::VectorXd& outv,
//                                           const std::vector<Eigen::VectorXd> workpoly,
//                                           int iapex,
//                                           int igparam) const
//{
//  int n = static_cast<int> (in1.size());
//  Eigen::VectorXd work1(n);
//  Eigen::VectorXd work2(n);
//  Eigen::VectorXd work3(n);
//  Eigen::VectorXd deriv(n);
//  Op->prodMatVecInPlace(in2,work2);
//
//    for(int i = 0; i < n ;i++)
//    {
//       work1[i] = _coeffs.back() * in1[i];
//       deriv[i] = 0.;
//    }
//
//    for(int j = static_cast<int> (_coeffs.size())-2; j >= 0; j--)
//    {
//      Op->prodMatVecInPlace(work1,work1);
//      for (int i = 0; i<n ; i++)
//      {
//          work1[i] = _coeffs[j] * in1[i] + work1[i];
//      }
//    }
//}

void ClassicalPolynomial::evalDerivOpOptim(ShiftOpCs* /*shiftOp*/,
                                           VectorDouble& /*temp1*/,
                                           VectorDouble& /*temp2*/,
                                           VectorDouble& /*outv*/,
                                           const VectorVectorDouble& /*workpoly*/,
                                           int /*iapex*/,
                                           int /*igparam*/)
{
   //TODO Call the Eigen::VectorXd function (try to put it in the mother class Polynomial)
   messerr("evalDerivOpOptim is not implemented for vectorsDouble");
}


void ClassicalPolynomial::evalDerivOpOptim(ShiftOpCs* shiftOp,
                                           Eigen::VectorXd& temp1,
                                           Eigen::VectorXd& temp2,
                                           Eigen::VectorXd& outv,
                                           const std::vector<Eigen::VectorXd>& workpoly,
                                           int iapex,
                                           int igparam) const
{
  int degree = (int) _coeffs.size();

  MatrixSparse* S = shiftOp->getS();
  MatrixSparse* gradS = shiftOp->getSGrad(iapex,igparam);

  shiftOp->getSGrad(iapex, igparam)->prodMatVecInPlace(workpoly[degree - 1], outv);

  for (int i = degree - 3; i >= 0; i--)
  {
    S->prodMatVecInPlace(outv, temp1);
    gradS->prodMatVecInPlace(workpoly[i + 1], temp2);
    VectorEigen::addInPlace(temp1, temp2, outv);
  }
}
