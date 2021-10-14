/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/

#include "Polynomials/ClassicalPolynomial.hpp"


#include "Basic/Vector.hpp"
#include "LinearOp/ShiftOpCs.hpp"

#include "csparse_d.h"
#include "csparse_f.h"

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
//(since it adds the result to the input vector, the classical scheme can t be used)
void ClassicalPolynomial::evalOpCumul(cs* Op, const VectorDouble& in, VectorDouble& out) const
{
  int n = static_cast<int> (in.size());
  VectorDouble work(n);
  VectorDouble work2(n);
  VectorDouble *swap1,*swap2,*swap3;

  swap1 = &work;
  swap2 = &work2;


  for (int i = 0; i<n ; i++)
  {
    out[i] += _coeffs[0] * in[i];
  }

  cs_vecmult(Op,in.data(),swap1->data());

  for(int j=1;j<(int)_coeffs.size();j++)
  {

    for (int i = 0; i<n ; i++)
    {
        out[i] += _coeffs[j] * (*swap1)[i];
    }

    if(j <(int) _coeffs.size()-1)
    {
         cs_vecmult(Op,swap1->data(),swap2->data());
         swap3 = swap1;
         swap1 = swap2;
         swap2 = swap3;
    }
  }
}

// Classical Hörner scheme starting from the highest degree
void ClassicalPolynomial::evalOp(cs* Op, const VectorDouble& in, VectorDouble& out) const
{
  int n = static_cast<int> (in.size());
  VectorDouble work(n);


  for(int i = 0; i < n ;i++)
  {
     out[i] = _coeffs.back() * in[i];
  }

  for(int j = static_cast<int> (_coeffs.size())-2; j >= 0; j--)
  {
    cs_vecmult(Op,out.data(),work.data());
    for (int i = 0; i<n ; i++)
    {
        out[i] = _coeffs[j] * in[i] + work[i];
    }
  }
}

// Classical Hörner scheme starting from the highest degree
void ClassicalPolynomial::evalOpTraining(cs* Op, const VectorDouble& in,VectorVectorDouble& store) const
{
  int n = static_cast<int> (in.size());

  VectorDouble work(n);


  for(int i = 0; i < n ;i++)
  {
     store[0][i] = _coeffs.back() * in[i];
  }

  for(int j = static_cast<int> (_coeffs.size())-2; j >= 0; j--)
  {
    cs_vecmult(Op,store[j].data(),work.data());
    for (int i = 0; i<n ; i++)
    {
        store[j][i] = _coeffs[j] * in[i] + work[i];
    }
  }
}



void ClassicalPolynomial::evalDerivOp(ShiftOpCs* shiftOp,
                                      const VectorDouble& in,
                                      VectorDouble& out,
                                      int iapex,
                                      int igparam)const
{
  int n = static_cast<int> (in.size());
  int degree = static_cast<int> (getCoeffs().size());
  ClassicalPolynomial* polycur = (ClassicalPolynomial*)this->clone();
  VectorDouble work(n);
  VectorDouble work2(n);
  VectorDouble *swap1,*swap2,*swap3;
  cs* Op = shiftOp->getS();
  cs* derivOp = shiftOp->getSGrad(iapex,igparam);

  swap1 = &work;
  swap2 = &work2;

  for(auto &e : out)
  {
    e = 0;
  }

  for(int i = 0 ; i< n; i++)
  {
    work[i] = in[i];
  }

  auto coeffsCur = polycur->getCoeffs();

  for(int i = 0; i < degree - 1 ;i++)
  {
    cs_vecmult(derivOp,swap1->data(),swap2->data());
    coeffsCur.erase(coeffsCur.begin());
    polycur->init(coeffsCur);
    polycur->evalOpCumul(Op,*swap2,out);

    if(i<degree-2) // to avoid useless and time consuming computation since it prepares next iteration
    {
      cs_vecmult(Op,swap1->data(),swap2->data());
      swap3 = swap1;
      swap1 = swap2;
      swap2 = swap3;
    }
   }
   delete polycur;
}



void ClassicalPolynomial::evalDerivOpOptim(cs* Op,
                                      const VectorDouble& in1,
                                      const VectorDouble& in2,
                                      VectorDouble& out
                                      )const
{

  int n = static_cast<int> (in1.size());
  VectorDouble work1(n);
  VectorDouble work2(n);
  VectorDouble work3(n);
  VectorDouble deriv(n);
  cs_vecmult(Op,in2.data(),work2.data());

    for(int i = 0; i < n ;i++)
    {
       work1[i] = _coeffs.back() * in1[i];
       deriv[i] = 0.;
    }

    for(int j = static_cast<int> (_coeffs.size())-2; j >= 0; j--)
    {
      cs_vecmult(Op,work1.data(),work1.data());
      for (int i = 0; i<n ; i++)
      {
          work1[i] = _coeffs[j] * in1[i] + work1[i];
      }
    }

}
