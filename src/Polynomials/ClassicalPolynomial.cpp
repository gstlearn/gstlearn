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
  for (int i = _coeffs.size()-2; i >= 0; i--)
  {
    result *= x;
    result +=_coeffs[i];
  }
  return result;
}

// Horner scheme starting from the lowest degree (since it add the result to the input vector)
void ClassicalPolynomial::evalOpCumul(cs* Op, const VectorDouble& in, VectorDouble& out) const
{
  int n = in.size();
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
  int n = in.size();
  VectorDouble work(n);


  for(int i = 0; i < n ;i++)
  {
     out[i] = _coeffs.back() * in[i];
  }

  for(int j=_coeffs.size()-2;j>=0;j--)
  {
    cs_vecmult(Op,out.data(),work.data());
    for (int i = 0; i<n ; i++)
    {
        out[i] = _coeffs[j] * in[i] + work[i];
    }
  }
}



void ClassicalPolynomial::evalDerivOp(ShiftOpCs* shiftOp,
                                      const VectorDouble& in,
                                      VectorDouble& out,
                                      int iapex,
                                      int igparam)const
{
  int n = in.size();
  int degree = getCoeffs().size();
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
    coeffsCur.erase(coeffsCur.begin(),coeffsCur.begin()+1);
    polycur->display();
    polycur->init(coeffsCur);
    polycur->evalOpCumul(Op,*swap2,out);

    if(i<degree-2)
    {
      cs_vecmult(Op,swap1->data(),swap2->data());
    }

    swap3 = swap1;
    swap1 = swap2;
    swap2 = swap3;
   }

   delete polycur;
}
