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
#include "Basic/VectorHelper.hpp"
#include "LinearOp/ShiftOpCs.hpp"

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
void ClassicalPolynomial::evalOpCumul(MatrixSparse* Op, const VectorDouble& inv, VectorDouble& outv) const
{
  int n = static_cast<int> (inv.size());
  VectorDouble work(n);
  VectorDouble work2(n);
  VectorDouble *swap1,*swap2,*swap3;

  swap1 = &work;
  swap2 = &work2;

  for (int i = 0; i<n ; i++)
  {
    outv[i] += _coeffs[0] * inv[i];
  }

  Op->prodMatVec(inv, *swap1);

  for (int j = 1; j < (int) _coeffs.size(); j++)
  {
    for (int i = 0; i < n; i++)
    {
      outv[i] += _coeffs[j] * (*swap1)[i];
    }

    if (j < (int) _coeffs.size() - 1)
    {
      Op->prodMatVec(*swap1, *swap2);
      swap3 = swap1;
      swap1 = swap2;
      swap2 = swap3;
    }
  }
}

// Classical Hörner scheme starting from the highest degree
void ClassicalPolynomial::evalOp(MatrixSparse* Op,
                                 const VectorDouble& inv,
                                 VectorDouble& outv) const
{
  int n = static_cast<int>(inv.size());
  VectorDouble work(n);

  for (int i = 0; i < n; i++)
    outv[i] = _coeffs.back() * inv[i];

  for (int j = static_cast<int>(_coeffs.size()) - 2; j >= 0; j--)
  {
    Op->prodMatVec(outv, work);
    for (int i = 0; i < n; i++)
    {
      outv[i] = _coeffs[j] * inv[i] + work[i];
    }
  }
}

// Classical Hörner scheme starting from the highest degree
void ClassicalPolynomial::evalOpTraining(MatrixSparse *Op,
                                         const VectorDouble &inv,
                                         VectorVectorDouble &store,
                                         VectorDouble &work) const
{
  int n = static_cast<int>(inv.size());

  if (work.empty())
  {
    work.resize(n);
  }

  for (int i = 0; i < n; i++)
  {
    store[_coeffs.size() - 1][i] = _coeffs.back() * inv[i];
  }

  for (int j = (int) _coeffs.size() - 2; j >= 0; j--)
  {
    Op->prodMatVec(store[j + 1], work);
    for (int i = 0; i < n; i++)
    {
      store[j][i] = _coeffs[j] * inv[i] + work[i];
    }
  }
}
#endif
void ClassicalPolynomial::evalDerivOp(ShiftOpCs* shiftOp,
                                      const VectorDouble& inv,
                                      VectorDouble& outv,
                                      int iapex,
                                      int igparam)const
{
  int n = static_cast<int> (inv.size());
  int degree = static_cast<int> (getCoeffs().size());
  ClassicalPolynomial* polycur = this->clone();
  VectorDouble work(n);
  VectorDouble work2(n);
  VectorDouble *swap1,*swap2,*swap3;
  MatrixSparse* Op = shiftOp->getS();
  MatrixSparse* derivOp = shiftOp->getSGrad(iapex,igparam);

  swap1 = &work;
  swap2 = &work2;

  for(auto &e : outv)
  {
    e = 0;
  }

  for(int i = 0 ; i< n; i++)
  {
    work[i] = inv[i];
  }

  auto coeffsCur = polycur->getCoeffs();

  for(int i = 0; i < degree - 1 ;i++)
  {
    derivOp->prodMatVec(*swap1,*swap2);
    coeffsCur.erase(coeffsCur.begin());
    polycur->init(coeffsCur);
    polycur->evalOpCumul(Op,*swap2,outv);

    if(i<degree-2) // to avoid useless and time consuming computation since it prepares next iteration
    {
      Op->prodMatVec(*swap1,*swap2);
      swap3 = swap1;
      swap1 = swap2;
      swap2 = swap3;
    }
   }
   delete polycur;
}

//void ClassicalPolynomial::evalDerivOpOptim(ShiftOpCs* shiftOp,
//                                           const VectorDouble& in1,
//                                           VectorDouble& in2,
//                                           VectorDouble& outv,
//                                           const VectorVectorDouble workpoly,
//                                           int iapex,
//                                           int igparam) const
//{
//  int n = static_cast<int> (in1.size());
//  VectorDouble work1(n);
//  VectorDouble work2(n);
//  VectorDouble work3(n);
//  VectorDouble deriv(n);
//  Op->prodMatVec(in2,work2);
//
//    for(int i = 0; i < n ;i++)
//    {
//       work1[i] = _coeffs.back() * in1[i];
//       deriv[i] = 0.;
//    }
//
//    for(int j = static_cast<int> (_coeffs.size())-2; j >= 0; j--)
//    {
//      Op->prodMatVec(work1,work1);
//      for (int i = 0; i<n ; i++)
//      {
//          work1[i] = _coeffs[j] * in1[i] + work1[i];
//      }
//    }
//}

void ClassicalPolynomial::evalDerivOpOptim(ShiftOpCs* shiftOp,
                                           VectorDouble& temp1,
                                           VectorDouble& temp2,
                                           VectorDouble& outv,
                                           const VectorVectorDouble workpoly,
                                           int iapex,
                                           int igparam) const
{
  int degree = (int) _coeffs.size();

  MatrixSparse* S = shiftOp->getS();
  MatrixSparse* gradS = shiftOp->getSGrad(iapex,igparam);

  shiftOp->getSGrad(iapex, igparam)->prodMatVec(workpoly[degree - 1], outv);

  for (int i = degree - 3; i >= 0; i--)
  {
    S->prodMatVec(outv, temp1);
    gradS->prodMatVec(workpoly[i + 1], temp2);
    VH::addInPlace(temp1, temp2, outv);
  }
}
