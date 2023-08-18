/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "Polynomials/ClassicalPolynomial.hpp"
#include "Basic/VectorHelper.hpp"
#include "LinearOp/ShiftOpCs.hpp"

#include "Matrix/LinkMatrixSparse.hpp"
#include "Matrix/LinkMatrixSparse.hpp"

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
void ClassicalPolynomial::evalOpCumul(cs* Op, const VectorDouble& inv, VectorDouble& outv) const
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

  cs_vecmult(Op, (int) swap1->size(), inv.data(), swap1->data());

  for (int j = 1; j < (int) _coeffs.size(); j++)
  {
    for (int i = 0; i < n; i++)
    {
      outv[i] += _coeffs[j] * (*swap1)[i];
    }

    if (j < (int) _coeffs.size() - 1)
    {
      cs_vecmult(Op, (int) swap2->size(), swap1->data(), swap2->data());
      swap3 = swap1;
      swap1 = swap2;
      swap2 = swap3;
    }
  }
}

// Classical Hörner scheme starting from the highest degree
void ClassicalPolynomial::evalOp(cs* Op,
                                 const VectorDouble& inv,
                                 VectorDouble& outv) const
{
  int n = static_cast<int>(inv.size());
  VectorDouble work(n);

  for (int i = 0; i < n; i++)
  {
    outv[i] = _coeffs.back() * inv[i];
  }

  int nout = (int) work.size();
  for (int j = static_cast<int>(_coeffs.size()) - 2; j >= 0; j--)
  {
    cs_vecmult(Op, nout, outv.data(), work.data());
    for (int i = 0; i < n; i++)
    {
      outv[i] = _coeffs[j] * inv[i] + work[i];
    }
  }
}

// Classical Hörner scheme starting from the highest degree
void ClassicalPolynomial::evalOpTraining(cs* Op, const VectorDouble& inv,VectorVectorDouble& store,VectorDouble& work) const
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

  int nout = (int) work.size();
  for (int j = (int) _coeffs.size() - 2; j >= 0; j--)
  {
    cs_vecmult(Op, nout, store[j + 1].data(), work.data());
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
  cs* Op = shiftOp->getS();
  cs* derivOp = shiftOp->getSGrad(iapex,igparam);

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

  int nout = (int) swap2->size();
  for(int i = 0; i < degree - 1 ;i++)
  {
    cs_vecmult(derivOp, nout, swap1->data(),swap2->data());
    coeffsCur.erase(coeffsCur.begin());
    polycur->init(coeffsCur);
    polycur->evalOpCumul(Op,*swap2,outv);

    if(i<degree-2) // to avoid useless and time consuming computation since it prepares next iteration
    {
      cs_vecmult(Op, nout, swap1->data(),swap2->data());
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
//  cs_vecmult(Op,work2.size(),in2.data(),work2.data());
//
//    for(int i = 0; i < n ;i++)
//    {
//       work1[i] = _coeffs.back() * in1[i];
//       deriv[i] = 0.;
//    }
//
//    for(int j = static_cast<int> (_coeffs.size())-2; j >= 0; j--)
//    {
//      cs_vecmult(Op,work1.size(),work1.data(),work1.data());
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

  cs* S = shiftOp->getS();
  cs* gradS = shiftOp->getSGrad(iapex,igparam);

  cs_vecmult(shiftOp->getSGrad(iapex, igparam), (int) outv.size(),
             workpoly[degree - 1].data(), outv.data());

  for (int i = degree - 3; i >= 0; i--)
  {
    cs_vecmult(S, (int) temp1.size(), outv.data(), temp1.data());
    cs_vecmult(gradS, (int) temp2.size(),
               workpoly[i + 1].data(), temp2.data());
    VH::addInPlace(temp1, temp2, outv);
  }
}
