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
#include "Matrix/VectorEigen.hpp"
#include "Polynomials/Chebychev.hpp"
#include <Eigen/src/Core/Map.h>
#include <Eigen/src/Core/Matrix.h>
#include <functional>

#include <math.h>

PrecisionOpMultiConditional::PrecisionOpMultiConditional()
  :_multiPrecisionOp(std::vector<PrecisionOp*>())
  ,_multiProjData(std::vector<IProjMatrix*>())
  ,_varianceData()
  ,_ndat(0)
  ,_ncova(0)
  ,_work1(Eigen::VectorXd())
  ,_work1bis(Eigen::VectorXd())
  ,_work1ter(Eigen::VectorXd())
  ,_workdata(Eigen::VectorXd())
  ,_work2(std::vector<Eigen::VectorXd>())
  ,_work3(std::vector<Eigen::VectorXd>())
{

}

PrecisionOpMultiConditional::~PrecisionOpMultiConditional()
{
}

VectorVectorDouble PrecisionOpMultiConditional::computeRhs(const VectorDouble& datVal) const
{
  Eigen::Map<const Eigen::VectorXd> datvalm(datVal.data(),datVal.size());
  auto res = computeRhs(datvalm);
  VectorVectorDouble result;
  for (int i = 0; i < (int)res.size(); i++)
  {
    auto v = VectorEigen::copyIntoVD(res[i]);
    result.push_back(v);
  }
  return result;
}

std::vector<Eigen::VectorXd> PrecisionOpMultiConditional::computeRhs(const Eigen::VectorXd& datVal) const
{
  std::vector<Eigen::VectorXd>rhs(sizes());
  for(int i = 0, n = sizes(); i < n; i++)
  {
    rhs[i].resize(size(i));
  }
  computeRhsInPlace(datVal,rhs);
  return rhs;
}

void PrecisionOpMultiConditional::computeRhsInPlace(const Eigen::VectorXd& datVal, std::vector<Eigen::VectorXd>& rhs) const
{
  Eigen::VectorXd temp(datVal.size());

  for(int i = 0; i < static_cast<int>(datVal.size()) ; i++)
  {
    temp[i] = datVal[i] / getVarianceData(i);
  }

  for(int i = 0; i < sizes(); i++)
  {
    _multiProjData[i]->point2mesh(temp, rhs[i]);
  }
}


int PrecisionOpMultiConditional::push_back(PrecisionOp *pmatElem,
                                           IProjMatrix *projDataElem)
{
  if (sizes() == 0 && projDataElem != nullptr)
  {
    // From the first element of the list, get the number of data and set _ndat
    _ndat = projDataElem->getPointNumber(); //TODO Vérifier la cohérence. _ndat doit coïncider pour tous les projDataElem.
    _work1.resize(_ndat);
    _workdata.resize(_ndat);
  }

  // Check that '_ndat' is the same for all IProjMatrix

  for (int i = 0, n = sizes(); i < n; i++)
  {
    int ndatloc = _multiProjData[i]->getPointNumber();
    if (ndatloc != _ndat)
    {
      messerr("The Projection matrix for element %d refers to %d data",i,ndatloc);
      messerr("It should be %d as for the others", _ndat);
      return 1;
    }
  }
  _multiPrecisionOp.push_back(pmatElem);
  _work2.push_back(Eigen::VectorXd(pmatElem->getSize()));
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
    VectorEigen::fill(e, 1.);
  }
  _AtA(_work3,_work2);
  return VectorEigen::maximum(_work2);
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

double PrecisionOpMultiConditional::computeLogDetOp(int nbsimu, int seed) const
{
  Chebychev logPoly;
  preparePoly(logPoly);

  law_set_random_seed(seed);
  
  double val = 0.;
  for (int i = 0; i < nbsimu; i++)
  {
    for (int j = 0; j < _ncova; j++)
    { 
      VectorEigen::simulateGaussianInPlace(_work1);
      VectorEigen::fill(_work1bis,0.);
      logPoly.addEvalOp(_multiPrecisionOp[j], _work1, _work1bis);
      val += _work1.adjoint() * _work1bis;
    }
  }
  return val / nbsimu;
}

double PrecisionOpMultiConditional::computeLogDetQ(int nbsimu, int seed) const
{
  double result = 0.;
  for (const auto &e : _multiPrecisionOp)
  {
    result += e->getLogDeterminant(nbsimu,seed);
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

double PrecisionOpMultiConditional::computeTotalLogDet(int nbsimu , int seed ) const
{
  double a1 = computeLogDetOp(nbsimu,seed);
  double a2 = computeLogDetQ(nbsimu,seed);
  double a3 = sumLogVar();
  return a1 - a2 + a3;
}

double PrecisionOpMultiConditional::computeQuadratic(const Eigen::VectorXd& x) const
{
  evalInvCov(x,_work1ter);
  return x.adjoint() *_work1ter;
}

void PrecisionOpMultiConditional::_AtA(const std::vector<Eigen::VectorXd>& inv, std::vector<Eigen::VectorXd>& outv) const
{
  VectorEigen::fill(_workdata,0.);

  for (int imod = 0; imod < sizes(); imod++)
  {
    _multiProjData[imod]->mesh2point(inv[imod], _work1);
    VectorEigen::addInPlace(_work1, _workdata);
  }

  Eigen::Map<const Eigen::VectorXd> vm(_varianceData.data(), _varianceData.size());
  VectorEigen::divideInPlace(vm,_workdata);

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
void PrecisionOpMultiConditional::_evalDirect(const std::vector<Eigen::VectorXd>& inv,
                                              std::vector<Eigen::VectorXd>& outv) const
{
  prepare();
  _AtA(inv,_work2);
  for (int imod = 0; imod < sizes(); imod++)
    _multiPrecisionOp[imod]->evalDirect(inv[imod], outv[imod]);
  VectorEigen::addInPlace(_work2, outv, outv);
}

void PrecisionOpMultiConditional::simulateOnMeshings(std::vector<Eigen::VectorXd> &result) const
{
  for (int icov = 0, ncov = (int) _multiPrecisionOp.size(); icov < ncov; icov++)
    simulateOnMeshing(result[icov], icov);
}

void PrecisionOpMultiConditional::simulateOnMeshing(Eigen::VectorXd &result,
                                                    int icov) const
{
  Eigen::VectorXd gauss(_multiPrecisionOp[icov]->getSize());
  VectorEigen::simulateGaussianInPlace(gauss);
  _multiPrecisionOp[icov]->evalSimulate(gauss, result);
}

void PrecisionOpMultiConditional::simulateOnDataPointFromMeshings(const std::vector<Eigen::VectorXd> &simus,
                                                                  Eigen::VectorXd& result) const
{
  result.resize(_ndat);
  VectorEigen::fill(result,0.);

  for(int icov = 0; icov <  sizes(); icov++)
  {
    _multiProjData[icov]->mesh2point(simus[icov],_work1);
    VectorEigen::addInPlace(_work1,result);
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

void PrecisionOpMultiConditional::evalInvCov(const Eigen::VectorXd& inv, Eigen::VectorXd& result) const
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
  MatrixSquareSymmetric XtInvSigmaX(xsize);
  VectorDouble result(xsize);

  for(int i = 0; i< xsize; i++)
  {
    Eigen::Map<const Eigen::VectorXd> xm(X[i].data(),X[i].size());
    evalInvCov(xm,_work1ter);

    Eigen::Map<const Eigen::VectorXd> ym(Y.data(),Y.size());
    XtInvSigmaZ[i] = ym.adjoint() * _work1ter;

    for(int j = i; j < xsize;j++)
    {
      Eigen::Map<const Eigen::VectorXd> xmj(X[j].data(),X[j].size());
      XtInvSigmaX.setValue(i,j,  xmj.adjoint() * _work1ter);
    }
  }

  XtInvSigmaX.solve(XtInvSigmaZ,result);

  return result;
}

