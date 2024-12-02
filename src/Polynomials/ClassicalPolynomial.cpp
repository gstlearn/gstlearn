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
#include "Basic/VectorNumT.hpp"
#include "Enum/EOperator.hpp"
#include "LinearOp/ALinearOp.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/NF_Triplet.hpp"
#include "geoslib_define.h"

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
void ClassicalPolynomial::evalOpCumul(MatrixSparse* Op,
                                      const constvect inv,
                                      vect outv) const
{
  int n = static_cast<int> (inv.size());
  VectorDouble work(n);
  VectorDouble work2(n);
  vect *swap1,*swap2,*swap3;

  vect ws(work);
  vect w2s(work2);
  swap1 = &ws;
  swap2 = &w2s;

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

void ClassicalPolynomial::addEvalOp(ALinearOp* Op,
                                    const constvect inv,
                                    vect outv) const
{
  int n = static_cast<int> (inv.size());

  vect swap1,swap2,swap3;

  if ((int)_work.size()!= n)
  {
    _work.resize(n);
  }
  if ((int)_work2.size()!= n)
  {
    _work2.resize(n);
  }

  swap1 = vect(_work);
  swap2 = vect(_work2);

  for (int i = 0; i<n ; i++)
  {
    outv[i] += _coeffs[0] * inv[i];
  }

  Op->evalDirect(inv, swap1);

  for (int j = 1; j < (int) _coeffs.size(); j++)
  {
    for (int i = 0; i < n; i++)
    {
      outv[i] += _coeffs[j] * (swap1)[i];
    }

    if (j < (int) _coeffs.size() - 1)
    {
      Op->evalDirect(swap1, swap2);
      swap3 = swap1;
      swap1 = swap2;
      swap2 = swap3;
    }
  }
}

// Classical Hörner scheme starting from the highest degree
void ClassicalPolynomial::evalOp(MatrixSparse* Op,
                                 const constvect inv,
                                 vect outv) const
{
  int n = static_cast<int>(inv.size());
  std::vector<double> work(n);
  vect ws(work);
  for (int i = 0; i < n; i++)
    outv[i] = _coeffs.back() * inv[i];

  int degree = (int) _coeffs.size();
  for (int j = degree - 2; j >= 0; j--)
  {
    Op->prodMatVecInPlace(outv, ws);
    for (int i = 0; i < n; i++)
    {
      outv[i] = _coeffs[j] * inv[i] + work[i];
    }
  }
}

/**
 * @brief Returns the rank-th term of the Diagonal of 'Op'
 * in its Polynomail expression through Horner mechanism
 * It is similar to the method 'evalOp' but targets the diagonal only
 * 
 * @param S Target Sparse matrix (possibly not even concretized)
 * @param rank Rank of the target
 * @return double 
 */
double ClassicalPolynomial::evalOpByRank(MatrixSparse* S, int rank) const
{
  int degree = (int)_coeffs.size();
  
  MatrixSparse* work = S->getColumnAsMatrixSparse(rank, _coeffs.back());
  MatrixSparse* outv = nullptr;
  for (int j = degree - 2; j >= 0; j--)
  {
    if (j != degree - 2) work->prodMatMatInPlace(S, outv);
    if (j == 0) break;
    delete outv;
    outv = work->clone();
    outv->setValue(rank, 0, outv->getValue(rank, 0) + _coeffs[j]);
  }

  double retval = work->getValue(rank,0) + _coeffs[0];
  delete outv;
  delete work;
  return retval;
}

// Classical Hörner scheme starting from the highest degree
void ClassicalPolynomial::evalOpTraining(
  MatrixSparse* Op,
  const constvect inv,
  std::vector<std::vector<double>>& store,
  std::vector<double>& work) const
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
    constvect stores(store[j+1]);
    vect ws(work);
    Op->prodMatVecInPlace(stores, ws);
    for (int i = 0; i < n; i++)
    {
      store[j][i] = _coeffs[j] * inv[i] + work[i];
    }
  }
}
#endif
// void ClassicalPolynomial::evalDerivOp(ShiftOpMatrix* shiftOp,
//                                       const constvect& inv,
//                                       vect& outv,
//                                       int iapex,
//                                       int igparam)const
// {
//   int n = static_cast<int> (inv.size());
//   int degree = static_cast<int> (getCoeffs().size());
//   ClassicalPolynomial* polycur = this->clone();
//   VectorDouble work(n);
//   VectorDouble work2(n);
//   vect ws(work);
//   vect w2s(work2);
//   vect *swap1,*swap2,*swap3;
//   MatrixSparse* Op = shiftOp->getS();
//   MatrixSparse* derivOp = shiftOp->getSGrad(iapex,igparam);

//   swap1 = &ws;
//   swap2 = &w2s;

//   std::fill(outv.begin(),outv.end(),0.);

//   for(int i = 0 ; i< n; i++)
//   {
//     work[i] = inv[i];
//   }

//   auto coeffsCur = polycur->getCoeffs();

//   for(int i = 0; i < degree - 1 ;i++)
//   {
//     derivOp->prodMatVecInPlace(*swap1,*swap2);
//     coeffsCur.erase(coeffsCur.begin());
//     polycur->init(coeffsCur);
//     polycur->evalOpCumul(Op,*swap2,outv);

//     if(i<degree-2) // to avoid useless and time consuming computation since it prepares next iteration
//     {
//       Op->prodMatVecInPlace(*swap1,*swap2);
//       swap3 = swap1;
//       swap1 = swap2;
//       swap2 = swap3;
//     }
//    }
//    delete polycur;
// }

// void ClassicalPolynomial::evalDerivOp(ShiftOpMatrix* shiftOp,
//                                       const VectorDouble& inv,
//                                       VectorDouble& outv,
//                                       int iapex,
//                                       int igparam)
// {
//   DECLARE_UNUSED(shiftOp);
//   DECLARE_UNUSED(inv);
//   DECLARE_UNUSED(outv);
//   DECLARE_UNUSED(iapex);
//   DECLARE_UNUSED(igparam);
//  //TODO Call the Eigen::VectorXd function
//  messerr("evalDerivOp is not implemented for vectorsDouble");
// }

//void ClassicalPolynomial::evalDerivOpOptim(ShiftOpMatrix* shiftOp,
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

// void ClassicalPolynomial::evalDerivOpOptim(ShiftOpMatrix* shiftOp,
//                                            VectorDouble& temp1,
//                                            VectorDouble& temp2,
//                                            VectorDouble& outv,
//                                            const VectorVectorDouble& workpoly,
//                                            int iapex,
//                                            int igparam)
// {
//   DECLARE_UNUSED(shiftOp);
//   DECLARE_UNUSED(temp1);
//   DECLARE_UNUSED(temp2);
//   DECLARE_UNUSED(outv);
//   DECLARE_UNUSED(workpoly);
//   DECLARE_UNUSED(iapex);
//   DECLARE_UNUSED(igparam);
//    //TODO Call the Eigen::VectorXd function (try to put it in the mother class Polynomial)
//    messerr("evalDerivOpOptim is not implemented for vectorsDouble");
// }


// void ClassicalPolynomial::evalDerivOpOptim(ShiftOpMatrix* shiftOp,
//                                            vect& temp1,
//                                            vect& temp2,
//                                            vect& outv,
//                                            const VectorVectorDouble& workpoly,
//                                            int iapex,
//                                            int igparam) 
// {
//   DECLARE_UNUSED(shiftOp,temp1,temp2,outv,workpoly,iapex,igparam)
//   // int degree = (int) _coeffs.size();

//   // MatrixSparse* S = shiftOp->getS();
//   // MatrixSparse* gradS = shiftOp->getSGrad(iapex,igparam);

//   // shiftOp->getSGrad(iapex, igparam)->prodMatVecInPlace(workpoly[degree - 1], outv);

//   // for (int i = degree - 3; i >= 0; i--)
//   // {
//   //   S->prodMatVecInPlace(outv, temp1);
//   //   gradS->prodMatVecInPlace(workpoly[i + 1], temp2);
//   //   VectorEigen::addInPlace(temp1, temp2, outv);
//   //}
// }
