/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* Created on: 18 dec. 2020 by N. Desassis                                    */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "LinearOp/PrecisionOpMultiConditional.hpp"
#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Polynomials/Chebychev.hpp"

#include <functional>

#include <math.h>

PrecisionOpMultiConditional::PrecisionOpMultiConditional()
  :_multiPrecisionOp(std::vector<PrecisionOp*>())
  ,_multiProjData(std::vector<IProjMatrix*>())
  ,_varianceData()
  ,_ndat(0)
  ,_ncova(0)
  ,_work1(VectorDouble())
  ,_work1bis(VectorDouble())
  ,_work1ter(VectorDouble())
  ,_workdata(VectorDouble())
  ,_work2(VectorVectorDouble())
  ,_work3(VectorVectorDouble())
{

}


VectorVectorDouble PrecisionOpMultiConditional::computeRhs(const VectorDouble& datVal) const
{
  VectorVectorDouble rhs(sizes());
  for(int i = 0; i< sizes(); i++)
  {
    rhs[i].resize(size(i));
  }
  computeRhsInPlace(datVal,rhs);
  return rhs;

}

void PrecisionOpMultiConditional::computeRhsInPlace(const VectorDouble& datVal, VectorVectorDouble& rhs) const
{
  VectorDouble temp = datVal;
  for(int i = 0; i < static_cast<int>(datVal.size()) ; i++)
  {
    temp[i] /= getVarianceData(i);
  }

  for(int i = 0; i < sizes(); i++)
  {
    _multiProjData[i]->point2mesh(temp,rhs[i]);
  }
}

void PrecisionOpMultiConditional::push_back(PrecisionOp* pmatElem,
                                            IProjMatrix* projDataElem)
{
  if (sizes() == 0 && projDataElem != nullptr)
  {
    _ndat = projDataElem->getPointNumber(); //TODO Vérifier la cohérence. _ndat doit coïncider pour tous les projDataElem.
    _work1.resize(_ndat);
    _workdata.resize(_ndat);
  }
  _multiPrecisionOp.push_back(pmatElem);
  _work2.push_back(VectorDouble(pmatElem->getSize()));
  _multiProjData.push_back(projDataElem);
  _updated();
  _ncova++;
}

std::pair<double,double> PrecisionOpMultiConditional::rangeEigenValQ() const
{
  std::pair<double,double> result = _multiPrecisionOp[0]->getRangeEigenVal();

  for (int i = 1; i < (int)_multiPrecisionOp.size(); i++)
  {
    std::pair<double,double> vals =_multiPrecisionOp[i]->getRangeEigenVal();
    result.first  = MIN(result.first ,vals.first);
    result.second = MAX(result.second,vals.second);
  }

  return result;
}

/* Evaluate the max along columns of the sum along lines of AtA. */
/*Since the terms of A are positive, we can compute AtA * 1_n  and take the max */

double PrecisionOpMultiConditional::getMaxEigenValProj() const
{

  _allocate(3);
  fillVal(_work3,1.);
  AtA(_work3,_work2);
  return ALinearOpMulti::max(_work2);
}

std::pair<double,double> PrecisionOpMultiConditional::computeRangeEigenVal() const
{
  std::pair<double,double> result = rangeEigenValQ();
  result.second += getMaxEigenValProj();
  return result;
}


 void PrecisionOpMultiConditional::preparePoly(Chebychev& logPoly) const
{
  std::pair<double,double> ranges = computeRangeEigenVal();
  double a = ranges.first;
  double b = ranges.second;
  logPoly.setA(a);
  logPoly.setB(b);
  logPoly.setNcMax(1500);
  std::function<double(double)> f;
  f = [] (double val){return log(val);};
  logPoly.fit(f,a,b,2*EPSILON4/(a+b));
}

double PrecisionOpMultiConditional::computeLogDetOp(int nsimus, int seed) const
{
  Chebychev logPoly;
  preparePoly(logPoly);
  
  VectorVectorDouble gauss(_ncova);
  VectorVectorDouble out(_ncova);
  law_set_random_seed(seed);

  for (int i = 0; i<_ncova; i++)
  {
    int size = _multiPrecisionOp[i]->getSize();
    gauss[i] =  VectorDouble(size);
    out[i] = VectorDouble(size);
  }

  double val = 0.;

  for (int i = 0; i < nsimus; i++)
  {
    for (auto &e : gauss)
    {
      VH::simulateGaussianInPlace(e);
    }

    logPoly.evalOp(this, gauss, _work3);
    val += innerProduct(gauss, _work3);
  }

  return val / nsimus;

}

double PrecisionOpMultiConditional::computeLogDetQ(int nsimus, int seed) const
{
  double result = 0.;
  for (auto &e : _multiPrecisionOp)
  {
    result += e->computeLogDet(nsimus,seed);
  }
  return result;
}


double PrecisionOpMultiConditional::sumLogVar() const
{
  double s = 0.;
  for (auto &e : _varianceData)
  {
    s += log(e);
  }
  return s;
}

// We use the fact that log|Sigma| = log |Q + A^t diag^(-1) (sigma) A|- log|Q| + Sum(log sigma_i^2)

double PrecisionOpMultiConditional::computeTotalLogDet(int nsimus , int seed ) const
{
  double result = 0;
  result += computeLogDetOp(nsimus,seed);
  result -= computeLogDetQ(nsimus,seed);
  result += sumLogVar();
  return result;
}

double PrecisionOpMultiConditional::computeQuadratic(const VectorDouble& x) const
{
  evalInvCov(x,_work1ter);
  return VH::innerProduct(x,_work1ter);
}


void PrecisionOpMultiConditional::AtA(const VectorVectorDouble& inv, VectorVectorDouble& outv) const
{

  for(auto &e : _workdata)
  {
    e = 0.;
  }

  for (int imod = 0; imod < sizes(); imod++)
  {
    _multiProjData[imod]->mesh2point(inv[imod], _work1);
    VH::addInPlace(_workdata,_work1);
  }

   VH::divideInPlace(_workdata, _varianceData);

  for (int imod = 0; imod < sizes(); imod++)
  {
    _multiProjData[imod]->point2mesh(_workdata, outv[imod]);
  }

}

/*****************************************************************************/
/*!
** Compute diag(Q1,...,Qncova) x + 1/nugget [A1,...,Ancova]^t [A1,...,Ancova] x
** in a block form where ncova is the number of basic structures excluding the
** nugget effect. Qi are the precision matrices associated to each structure and Ai
** are the projection matrices from the meshing vertices to the data locations.
**
** \param[in]  inv     Array of input values
**
** \param[out] outv    Array of output values
**
*******************************************************************************/
void PrecisionOpMultiConditional::_evalDirect(const VectorVectorDouble& inv,
                                              VectorVectorDouble& outv) const
{
  _init();
  AtA(inv,_work2);
  for (int imod = 0; imod < sizes(); imod++)
    _multiPrecisionOp[imod]->eval(inv[imod], outv[imod]);
  sum(_work2,outv, outv);

}

void PrecisionOpMultiConditional::simulateOnMeshing(VectorDouble& gauss,VectorVectorDouble& result,int icov) const
{
    _multiPrecisionOp[icov]->simulateOneInPlace(gauss,result[icov]);
}

void PrecisionOpMultiConditional::simulateOnDataPointFromMeshings(const VectorVectorDouble& simus,
                                                                  VectorDouble& result) const
{
  VH::fill(result,0.,_ndat);

  for(int icov = 0; icov <  sizes(); icov++)
  {
    _multiProjData[icov]->mesh2point(simus[icov],_work1);
    VH::addInPlace(result,_work1);
  }

  for(int idat = 0; idat < _ndat; idat++)
  {
    result[idat]+= sqrt(_varianceData[idat]) * law_gaussian();
  }

}

void PrecisionOpMultiConditional::_allocate(int i) const
{
  if (i == 1)
  {
    if(_work1bis.empty())
    {
      _work1bis.resize(_ndat);
    }
  }

  if (i == 0)
  {
    if(_work1.empty())
    {
      _work1.resize(_ndat);
    }
  }

  if (i == 3)
  {
    _work3.resize(sizes());

    for(int j = 0; j<sizes(); j++)
    {
      if(_work3[j].empty())
      {
        _work3[j].resize(size(j));
      }
    }
  }
  if (i == 2)
  {
    if(_workdata.empty())
    {
      _workdata.resize(_ndat);
    }
  }
  if (i == 4)
  {
    if(_work1ter.empty())
    {
      _work1ter.resize(_ndat);
    }
  }

}


void PrecisionOpMultiConditional::evalInvCov(const VectorDouble& inv, VectorDouble& result) const
{

  _allocate(0);
  _allocate(1);
  _allocate(2);
  _allocate(3);
  _allocate(4);

  for(int idat = 0; idat < _ndat; idat++)
  {
    result[idat] = inv[idat]/_varianceData[idat];
  }

  for(int icov = 0; icov < sizes(); icov++)
  {
    _multiProjData[icov]->point2mesh(result,_work2[icov]);
  }

  evalInverse(_work2,_work3);

  for(int icov = 0; icov < sizes(); icov ++)
  {
     _multiProjData[icov]->mesh2point(_work3[icov],_work1bis);

     for(int idat = 0; idat < _ndat; idat++)
     {
        result[idat] -=  1./_varianceData[idat] * _work1bis[idat];
     }
  }
}

VectorDouble PrecisionOpMultiConditional::computeCoeffs(const VectorDouble& Y, const VectorVectorDouble& X) const
{

  _allocate(4);
  int xsize = static_cast<int>(X.size());
  VectorDouble XtInvSigmaZ(static_cast<int>(xsize));
  MatrixSquareSymmetric XtInvSigmaX(xsize,false);
  VectorDouble result(xsize);

  for(int i = 0; i< xsize; i++)
  {
    evalInvCov(X[i],_work1ter);
    XtInvSigmaZ[i] = VH::innerProduct(Y,_work1ter);

    for(int j = i; j < xsize;j++)
    {
      XtInvSigmaX.setValue(i,j, VH::innerProduct(X[j],_work1ter));
    }
  }

  XtInvSigmaX.solve(XtInvSigmaZ,result);



  return result;
}

PrecisionOpMultiConditional::~PrecisionOpMultiConditional(){}
