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

double PrecisionOpMultiConditional::getMaxEigenValProj() const
{
  double result = _multiProjData[0]->getMaxAtDinvA(_varianceData);

  for (int i = 1; i < (int)_multiProjData.size(); i++)
  {
    double vals =_multiProjData[i]->getMaxAtDinvA(_varianceData);
    result = MAX(result,vals);
  }

  return result;
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
  std::function<double(double)> f;
  f = [this](double val){return log(val);};
  logPoly.fit(f,a,b);
}

double PrecisionOpMultiConditional::computeLogDetOp(int nsimus, int seed) const
{
  Chebychev logPoly;
  preparePoly(logPoly);
  VectorVectorDouble gauss(_ncova);
  VectorVectorDouble out(_ncova);

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
      ut_vector_simulate_gaussian_inplace(e);
    }

//    logPoly.evalOp()
//    val += gauss[i] * result[i];
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


double PrecisionOpMultiConditional::computeTotalLogDet(int nsimus , int seed ) const
{
  double result = computeLogDetOp(nsimus,seed);
  result += computeLogDetQ(nsimus,seed);
  return result;
}

/*****************************************************************************/
/*!
** Compute diag(Q1,...,Qncova) x + 1/nugget [A1,...,Ancova]^t [A1,...,Ancova] x
** in a block form where ncova is the number of basic structures excluding the
** nugget effect. Qi are the precision matrices associated to each structure and Ai
** are the projection matrices from the meshing vertices to the data locations.
**
** \param[in]  in     Array of input values
**
** \param[out] out    Array of output values
**
*******************************************************************************/
void PrecisionOpMultiConditional::_evalDirect(const VectorVectorDouble& in,
                                              VectorVectorDouble& out) const
{
  _init();

  for(auto &e : _workdata)
  {
    e = 0.;
  }

  for (int imod = 0; imod < sizes(); imod++)
  {
    _multiProjData[imod]->mesh2point(in[imod], _work1);
    ut_vector_add_inplace(_workdata,_work1);
  }

  ut_vector_divide_vec(_workdata, _varianceData);

  for (int imod = 0; imod < sizes(); imod++)
  {
    _multiPrecisionOp[imod]->eval(in[imod], out[imod]);
     _multiProjData[imod]->point2mesh(_workdata, _work2[imod]);
   }

  _linearComb(1., _work2, 1., out, out);

}

void PrecisionOpMultiConditional::simulateOnMeshing(const VectorDouble& gauss,VectorVectorDouble& result) const
{
  for(int icov=0;icov< sizes();icov++)
  {
    _multiPrecisionOp[icov]->eval(gauss,result[icov]);
  }
}

void PrecisionOpMultiConditional::simulateOnDataPointFromMeshings(const VectorVectorDouble& simus,
                                                                  VectorDouble& result) const
{
  ut_vector_fill(result,0.,_ndat);

  for(int icov = 0; icov <  sizes(); icov++)
  {
    _multiProjData[icov]->mesh2point(simus[icov],_work1);
    ut_vector_add_inplace(result,_work1);
  }

  for(int idat = 0; idat < _ndat; idat++)
  {
    result[idat]+= sqrt(_varianceData[idat]) * law_gaussian();
  }

}

void PrecisionOpMultiConditional::evalInvCov(const VectorDouble& in, VectorDouble& result) const
{
  if(static_cast<int>(_work3.size()) <= 0)
  {
    _work3.resize(sizes());
  }

  if(static_cast<int>(_work1ter.size()) <= 0)
  {
      _work1ter.resize(_ndat);
  }

  for(int i = 0; i<sizes();i++)
  {
    if(_work3[i].empty())
    {
    _work3[i].resize(size(i));
    }
  }

  if(_work1bis.empty())
  {
    _work1bis.resize(_ndat);
  }

  for(int idat = 0; idat < _ndat; idat++)
  {
    result[idat] = in[idat]/_varianceData[idat];
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

  int xsize = static_cast<int>(X.size());
  VectorDouble XtInvSigmaZ(static_cast<int>(xsize));
  MatrixSquareSymmetric XtInvSigmaX(xsize,false);
  VectorDouble result(xsize);

  for(int i = 0; i< xsize; i++)
  {
    evalInvCov(X[i],_work1ter);
    XtInvSigmaZ[i] = ut_vector_inner_product(Y,_work1ter);

    for(int j = i; j < xsize;j++)
    {
      XtInvSigmaX.setValue(i,j,ut_vector_inner_product(X[j],_work1ter));
    }
  }

//  std::cout<< XtInvSigmaX.getValue(0,0)<< "  "<< XtInvSigmaX.getValue(0,1) <<std::endl;
//  std::cout<< XtInvSigmaX.getValue(1,0)<< "  "<< XtInvSigmaX.getValue(1,1) <<std::endl;
//  std::cout <<" ---------------------" <<std::endl;
//  std::cout<< XtInvSigmaZ[0] << "  "<< XtInvSigmaZ[1] <<std::endl;


  XtInvSigmaX.solve(XtInvSigmaZ,result);



  return result;
}

PrecisionOpMultiConditional::~PrecisionOpMultiConditional(){}
