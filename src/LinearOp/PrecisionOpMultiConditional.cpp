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
#include "LinearOp/PrecisionOpMultiConditional.hpp"
#include "Basic/Law.hpp"
#include "Basic/VectorNumT.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Polynomials/Chebychev.hpp"
#include "geoslib_define.h"
#include <functional>

#include <math.h>
#include <vector>
#include <algorithm>

PrecisionOpMultiConditional::PrecisionOpMultiConditional()
  :_multiPrecisionOp(std::vector<PrecisionOp*>())
  ,_multiProjData(std::vector<IProj*>())
  ,_varianceData()
  ,_ndat(0)
  ,_ncova(0)
  ,_work1()
  ,_work1bis()
  ,_work1ter()
  ,_workdata()
  ,_work2()
  ,_work3()
{

}

PrecisionOpMultiConditional::~PrecisionOpMultiConditional()
{
}


std::vector<std::vector<double>> PrecisionOpMultiConditional::computeRhs(const std::vector<double>& datVal) const
{
  std::vector<std::vector<double>>rhs(sizes());
  for(int i = 0, n = sizes(); i < n; i++)
  {
    rhs[i].resize(size(i));
  }
  computeRhsInPlace(datVal,rhs);
  return rhs;
}

void PrecisionOpMultiConditional::computeRhsInPlace(const std::vector<double>& datVal, std::vector<std::vector<double>>& rhs) const
{
  std::vector<double> temp(datVal.size());

  for(int i = 0; i < static_cast<int>(datVal.size()) ; i++)
  {
    temp[i] = datVal[i] / getVarianceData(i);
  }

  for(int i = 0, n = sizes(); i < n; i++)
  {
    constvect tempm(temp);
    vect rhsm(rhs[i]);
    _multiProjData[i]->point2mesh(tempm, rhsm);
  }
}


int PrecisionOpMultiConditional::push_back(PrecisionOp *pmatElem,
                                           IProj *projDataElem)
{
  if (sizes() == 0 && projDataElem != nullptr)
  {
    // From the first element of the list, get the number of data and set _ndat
    _ndat = projDataElem->getNPoint(); //TODO Vérifier la cohérence. _ndat doit coïncider pour tous les projDataElem.
    _work1.resize(_ndat);
    _workdata.resize(_ndat);
  }

  // Check that '_ndat' is the same for all IProj

  for (int i = 0, n = sizes(); i < n; i++)
  {
    int ndatloc = _multiProjData[i]->getNPoint();
    if (ndatloc != _ndat)
    {
      messerr("The Projection matrix for element %d refers to %d data",i,ndatloc);
      messerr("It should be %d as for the others", _ndat);
      return 1;
    }
  }
  _multiPrecisionOp.push_back(pmatElem);
  _work2.push_back(std::vector<double>(pmatElem->getSize()));
  _multiProjData.push_back(projDataElem);
  _updated();
  _ncova++;
  return 0;
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
/* Since the terms of A are positive, we can compute AtA * 1_n  and take the max */
double PrecisionOpMultiConditional::getMaxEigenValProj() const
{

  _allocate(3);
  for (auto &e: _work3)
  {
    std::fill(e.begin(),e.end(), 1.);
  }
  _AtA(_work3,_work2);
  return VectorHelper::maximum(_work2);
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

double PrecisionOpMultiConditional::computeLogDetOp(int nbsimu) const
{
  Chebychev logPoly;
  preparePoly(logPoly);  
  double val = 0.;
  for (int i = 0; i < nbsimu; i++)
  {
    for (int j = 0; j < _ncova; j++)
    { 
      VH::simulateGaussianInPlace(_work1);
      std::fill(_work1bis.begin(),_work1bis.end(),0.);
      vect w1s(_work1);
      vect w1bis(_work1bis);
      logPoly.addEvalOp(_multiPrecisionOp[j], w1s, w1bis);
      val += VH::innerProduct(_work1bis,_work1);
    }
  }
  return val / nbsimu;
}

double PrecisionOpMultiConditional::computeLogDetQ(int nbsimu) const
{
  double result = 0.;
  for (const auto &e : _multiPrecisionOp)
  {
    result += e->getLogDeterminant(nbsimu);
  }
  return result;
}

double PrecisionOpMultiConditional::sumLogVar() const
{
  double s = 0.;
  for (const auto &e : _varianceData)
  {
    s += log(e);
  }
  return s;
}

// We use the fact that log|Sigma| = log |Q + A^t diag^(-1) (sigma) A|- log|Q| + Sum(log sigma_i^2)

double PrecisionOpMultiConditional::computeTotalLogDet(int nbsimu ) const
{
  double a1 = computeLogDetOp(nbsimu);
  double a2 = computeLogDetQ(nbsimu);
  double a3 = sumLogVar();
  return a1 - a2 + a3;
}

double PrecisionOpMultiConditional::computeQuadratic(const std::vector<double>& x) const
{
  evalInvCov(x,_work1ter);
  return VH::innerProduct(_work1ter,x);
}

void PrecisionOpMultiConditional::_AtA(const std::vector<std::vector<double>>& inv, std::vector<std::vector<double>>& outv) const
{
  std::fill(_workdata.begin(),_workdata.end(),0.);

  for (int imod = 0; imod < sizes(); imod++)
  {
    constvect invs(inv[imod]);
    vect w1s(_work1);
    _multiProjData[imod]->mesh2point(invs, w1s);
    VH::addInPlace(_workdata,_work1);
  }

  VH::divideInPlace(_workdata,_varianceData);

  for (int imod = 0; imod < sizes(); imod++)
  {
    constvect wds(_workdata);
    vect outs(outv[imod]);
    _multiProjData[imod]->point2mesh(wds, outs);
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
void PrecisionOpMultiConditional::_evalDirect(const std::vector<std::vector<double>>& inv,
                                              std::vector<std::vector<double>>& outv) const
{
  prepare();
  _AtA(inv,_work2);
  for (int imod = 0; imod < sizes(); imod++)
  {
    constvect invm(inv[imod]);
    vect outm(outv[imod]);
    _multiPrecisionOp[imod]->evalDirect(invm, outm);
  }
  VH::addInPlace(_work2, outv, outv);
}

void PrecisionOpMultiConditional::simulateOnMeshings(std::vector<std::vector<double>> &result) const
{
  for (int icov = 0, ncov = (int) _multiPrecisionOp.size(); icov < ncov; icov++)
    simulateOnMeshing(result[icov], icov);
}

void PrecisionOpMultiConditional::simulateOnMeshing(std::vector<double> &result,
                                                    int icov) const
{
  VectorDouble gauss(_multiPrecisionOp[icov]->getSize());
  VH::simulateGaussianInPlace(gauss);
  constvect gaussm(gauss);
  vect resultm(result);

  _multiPrecisionOp[icov]->evalSimulate(gaussm, resultm);
}

void PrecisionOpMultiConditional::simulateOnDataPointFromMeshings(const std::vector<std::vector<double>> &simus,
                                                                  std::vector<double>& result) const
{
  result.resize(_ndat);
  std::fill(result.begin(),result.end(),0.);

  for(int icov = 0; icov <  sizes(); icov++)
  {
    constvect simuss(simus[icov]);
    vect w1s(_work1);
    _multiProjData[icov]->mesh2point(simuss,w1s);
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
    if(_work1bis.size() == 0)
    {
      _work1bis.resize(_ndat);
    }
  }

  if (i == 0)
  {
    if(_work1.size() == 0)
    {
      _work1.resize(_ndat);
    }
  }

  if (i == 3)
  {
    _work3.resize(sizes());

    for(int j = 0; j<sizes(); j++)
    {
      if(_work3[j].size() == 0)
      {
        _work3[j].resize(size(j));
      }
    }
  }
  if (i == 2)
  {
    if(_workdata.size() == 0)
    {
      _workdata.resize(_ndat);
    }
  }
  if (i == 4)
  {
    if(_work1ter.size() == 0)
    {
      _work1ter.resize(_ndat);
    }
  }
}

void PrecisionOpMultiConditional::evalInvCov(const constvect inv,
                                             std::vector<double>& result) const
{
  _allocate(0);
  _allocate(1);
  _allocate(2);
  _allocate(3);
  _allocate(4);

  for(int idat = 0; idat < _ndat; idat++)
  {
    result[idat] = inv[idat] / _varianceData[idat];
  }

  for(int icov = 0; icov < sizes(); icov++)
  {
    constvect results(result);
    vect w2s(_work2[icov]);
    _multiProjData[icov]->point2mesh(results,w2s);
  }
  evalInverse(_work2,_work3);

  for(int icov = 0; icov < sizes(); icov ++)
  {
     constvect w3s(_work3[icov]);
     vect w1bis(_work1bis);
     _multiProjData[icov]->mesh2point(w3s,w1bis);

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
  MatrixSquareSymmetric XtInvSigmaX(xsize);
  VectorDouble result(xsize);

  for(int i = 0; i< xsize; i++)
  {
    constvect xm(X[i].data(),X[i].size());
    evalInvCov(xm,_work1ter);

    constvect Ys(Y);
    constvect w1i(_work1ter);
    XtInvSigmaZ[i] = VH::innerProduct(Ys,w1i);

    for(int j = i; j < xsize;j++)
    {
      constvect xmj(X[j].data(),X[j].size());
      double prod = VH::innerProduct(xmj,w1i);
      XtInvSigmaX.setValue(i,j,prod);
    }
  }

  XtInvSigmaX.solve(XtInvSigmaZ,result);

  return result;
}

