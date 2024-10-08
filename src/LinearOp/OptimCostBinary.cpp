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
#include "LinearOp/LinearOpCGSolver.hpp"

#include "LinearOp/OptimCostBinary.hpp"
#include "LinearOp/HessianOp.hpp"
#include "LinearOp/IOptimCost.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "Basic/Law.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/VectorHelper.hpp"
#include <math.h>

OptimCostBinary::OptimCostBinary()
    : IOptimCost(),
      _isInitialized(false),
      _flagSeismic(false),
      _meanPropRaw(0.),
      _meanPropGaus(0.),
      _pMat(nullptr),
      _projData(nullptr),
      _projSeis(nullptr),
      _propSeis(),
      _varSeis(),
      _grad(),
      _workp(),
      _workx(),
      _workv(),
      _lambdav(),
      _works()
{
}

OptimCostBinary::OptimCostBinary(const OptimCostBinary &m)
    : IOptimCost(),
      _isInitialized(m._isInitialized),
      _flagSeismic(m._flagSeismic),
      _meanPropRaw(m._meanPropRaw),
      _meanPropGaus(m._meanPropGaus),
      _pMat(m._pMat),
      _projData(m._projData),
      _projSeis(m._projSeis),
      _propSeis(m._propSeis),
      _varSeis(m._varSeis),
      _grad(),
      _workp(),
      _workx(),
      _workv(),
      _lambdav()
{

}

OptimCostBinary& OptimCostBinary::operator = (const OptimCostBinary &m)
{
  if (this != &m)
  {
    _isInitialized = m._isInitialized;
    _flagSeismic = m._flagSeismic;
    _meanPropRaw = m._meanPropRaw;
    _meanPropGaus = m._meanPropGaus;
    _pMat = m._pMat;
    _projData = m._projData;
    _projSeis = m._projSeis;
    _propSeis = m._propSeis;
    _varSeis = m._varSeis;
  }
  return *this;
}

OptimCostBinary::~OptimCostBinary() 
{
}

/*****************************************************************************/
/*!
**  Initialize the Binary Cost Operator
**
** \param[in]  pmat     The precision matrix to be optimized
** \param[in]  projdata The Projection operator between Data and Meshing
** \param[in]  projseis The Projection operator between Seismic and Meshing
** \param[in]  propseis Array of facies proportions
** \param[in]  varseis  Array of variance attached to the seismic
**
*****************************************************************************/
void OptimCostBinary::reset(PrecisionOp* pmat,
                            const ProjMatrix* projdata,
                            const ProjMatrix* projseis,
                            const VectorDouble& propseis,
                            const VectorDouble& varseis)
{
  // Assignment of pointers
  _pMat     = pmat;
  _projData = projdata;
  _projSeis = projseis;
  _propSeis = propseis;
  _varSeis  = varseis;
  
  // Copy the memory chunks passed as arguments
  int nvertex = _projData->getApexNumber();
  int npoint  = _projData->getPointNumber();
  
  // Particular case of the Seismic

  _flagSeismic = (projseis != (ProjMatrix *) NULL && projseis->getPointNumber() > 0);
  if (_flagSeismic)
  {
    int nseis = _projSeis->getPointNumber();
    _works.resize(nseis);
  }

  // Auxiliary working arrays
  _grad.resize(nvertex);
  _workp.resize(npoint);
  _workx.resize(npoint);
  _workv.resize(nvertex);
  _lambdav.resize(nvertex);

  // Set the initialization flag
  _isInitialized = true;
}

/*****************************************************************************/
/*!
**  Perform the minimization
**
**  \return The array of facies proportions (Dimension: nvertex)
**
** \param[in]  indic      Array containing the Facies indicators (see remarks)
**                        (Dimension: npoint)
** \param[in]  verbose    Verbose flag
** \param[in]  maxiter    Maximum number of iterations for Optimization algo.
** \param[in]  eps        Tolerance for Optimization algorithm
**
** \remarks The argument 'indic' should contain 0 or 1 for active constraints
** \remarks The inactive constraints should be set to TEST.
**
*****************************************************************************/
VectorDouble OptimCostBinary::minimize(VectorDouble& indic,
                                       bool verbose,
                                       int maxiter,
                                       double eps)
{
  double normgrad, costv;
  bool flagContinue;
  int  iter;
  HessianOp *hess = nullptr;
  VectorDouble propfac;
  // Statistics on the input data (only for verbose option)

  if (verbose)
  {
    message("Mean proportion (provided as input) = %lf\n",_meanPropRaw);
    VH::displayStats("Proportions calculated on Data",indic);
  }

  try 
  {
    if (! _isInitialized) 
      my_throw("'OptimCostBinary' must be initialized beforehand");
    int nvertex = _pMat->getSize();

    // Instantiate the Hessian

    hess = new HessianOp();
    if (hess->init(_pMat, _projData, _projSeis, indic, _propSeis, _varSeis)) 
      my_throw("Problem in Hessian init() method");

    // Core allocation

    VectorDouble lambdat(nvertex, 0.);
    VectorDouble step(nvertex, 0.);
    propfac.resize(nvertex);
    propfac.fill(0.);
  
    // Evaluate the reference cost value

    costv = _evaluateCost(indic,propfac);

    // Iterations to minimize the Cost function

    iter = 0;
    flagContinue = true;
    while(flagContinue)
    {
      iter++;
      _evaluateGrad(indic,propfac,&normgrad);
      if (OptDbg::query(EDbg::CONVERGE))
        message("Iteration #%d (max=%d) - Cost=%lf - NormGrad=%lf (eps=%lg)\n",
                iter,maxiter,costv,normgrad,eps);

      flagContinue = (normgrad > eps && iter < maxiter);

      // Perform the Conjugate Gradient on the Hessian

      hess->setLambda(propfac);
      LinearOpCGSolver<HessianOp> solver(hess);
      solver.solve(_grad,step);
      
      bool flagSortie = false;
      while (! flagSortie)
      {
        for (int i=0; i<nvertex; i++) lambdat[i] = propfac[i] - step[i];
        double costt = _evaluateCost(indic,lambdat);
        if (costt < costv * 1.000001)
        {
          costv = costt;
          for (int i=0; i<nvertex; i++) propfac[i] = lambdat[i];
          flagSortie = true;
        }
        for (int i=0; i<nvertex; i++) step[i] /= 2.;
      }
    }

    // Final conversion from Gaussian to Proportion scale

    for (int i=0; i<nvertex; i++)
      propfac[i] = 1. - law_cdf_gaussian(propfac[i]);

    // Calculate the Proportion statistics (verbose option)

    if (verbose)
    {
      VH::displayStats("Calculated Proportions",propfac);
    }
  }
  catch(const std::string& str)
  {
    // TODO : Check if std::exception can be used
    messerr("%s", str.c_str());
  }

  delete(hess);
  return propfac;
}

/*****************************************************************************/
/*!
**  Calculate the Gradient
**
** \param[in]  indic      Array containing the Facies indicators (see remarks)
**                        (Dimension: npoint)
** \param[in]  lambda     Array of input values
**
** \param[out] out        Array of output gradients
**
*****************************************************************************/
void OptimCostBinary::calculateGradient(const VectorDouble& indic,
                                        const VectorDouble& lambda,
                                        double* out)
{
  double normgrad;

  _evaluateGrad(indic,lambda,&normgrad);

  for (int i=0; i<_projData->getApexNumber(); i++) out[i] = _grad[i];
}

/*****************************************************************************/
/*!
**  Returns the Number of Data Points
**
*****************************************************************************/
int OptimCostBinary::getNPoint() const
{
  if (! _isInitialized) 
    my_throw("'OptimCostBinary' must be initialized beforehand");
  return _projData->getPointNumber();
}

/*****************************************************************************/
/*!
**  Returns the Number of Meshing Vertices 
**
*****************************************************************************/
int OptimCostBinary::getNVertex() const
{
  if (! _isInitialized) 
    my_throw("'OptimCostBinary' must be initialized beforehand");
  return _projData->getApexNumber();
}

/*****************************************************************************/
/*!
**  Toggle the use of the Seismic constraint
**
** \param[in]  status     Status assigned to the Seismic Constraints
**
** \remarks When the Seismic is not defined (i.e. 'projSeis not defined)
** \remarks this function is useless
**
*****************************************************************************/
void OptimCostBinary::toggleSeismic(bool status)
{
  // When Seismic is not defined, this precedure is useless
  if (_projSeis == (ProjMatrix*)NULL || _projSeis->getPointNumber() <= 0)
    return;

  _flagSeismic = status;
}

/*****************************************************************************/
/*!
**  Set the Mean proportion for the indicator
**
** \param[in]  meanprop     Value of the mean proportion (raw scale)
**
*****************************************************************************/
int OptimCostBinary::setMeanProportion(double meanprop)
{
  if (meanprop < 0 || meanprop > 1)
  {
    messerr("The argument 'meanprop' should lie between 0 and 1 (%lf)");
    return 1;
  }
  _meanPropRaw  = meanprop;
  _meanPropGaus = law_invcdf_gaussian(1. - meanprop);
  return 0;
}

/*****************************************************************************/
/*!
**  Internal function to evaluate the Cost
**
** \param[in]  indic      Array containing the Facies indicators 
**                        (Dimension: npoint)
** \param[in]  lambda     Array of input values
**
*****************************************************************************/
double OptimCostBinary::_evaluateCost(const VectorDouble& indic,
                                      const VectorDouble& lambda)
{
  double result = 0.;

  // Contribution of the Data

  _projData->mesh2point(lambda,_workp);

  double sum_pos = 0.;
  double sum_neg = 0.;
  for (int i=0; i<_projData->getPointNumber(); i++)
  {
    if (FFFF(indic[i])) continue;
    if (indic[i] > 0.)
      sum_pos -= log(1. - law_cdf_gaussian(_workp[i]));
    else
      sum_neg -= log(     law_cdf_gaussian(_workp[i]));
  }
  result += sum_pos + sum_neg;

  // Contribution of the spatial structure 

  for (int i=0; i<_projData->getApexNumber(); i++) 
    _lambdav[i] = lambda[i] - _meanPropGaus;
  _pMat->evalDirect(_lambdav, _workv);

  double sum_str = 0.;
  for (int i=0; i<_projData->getApexNumber(); i++)
    sum_str += 0.5 * _lambdav[i] * _workv[i];
  result += sum_str;

  // Contribution of the seismic (optional)

  if (_flagSeismic)
  {
    _contributeSeismic(lambda);

    double sum_seis = 0.;
    for (int i=0; i<_projSeis->getPointNumber(); i++) 
      sum_seis += 0.5 * _works[i] * _varSeis[i] * _works[i];

    result += sum_seis;
  }
  return result;
}

/*****************************************************************************/
/*!
**  Internal function to evaluate the Gradient
**
** \param[in]  indic      Array containing the Facies indicators (see remarks)
**                        (Dimension: npoint)
** \param[in]  lambda     Array of input values
**
** \param[out] normgrad   Returnewd normalized gradient value
**
*****************************************************************************/
void OptimCostBinary::_evaluateGrad(const VectorDouble& indic,
                                    const VectorDouble& lambda,
                                    double* normgrad)
{

  // Contribution of the spatial structure
  // Array 'lambda' must be corrected from the contribution of the 'mean'

  for (int i=0; i<_projData->getApexNumber(); i++) 
    _lambdav[i] = lambda[i] - _meanPropGaus;
  _pMat->evalDirect(_lambdav,_grad);

  // Contribution of the Data

  _projData->mesh2point(lambda,_workp);
  for (int i=0; i<_projData->getPointNumber(); i++)
  {
    if (FFFF(indic[i]))
      _workp[i] = 0.;           // TODO: to be checked
    else
      _workp[i] = 
        law_df_gaussian(_workp[i]) / (indic[i] - law_cdf_gaussian(_workp[i]));
  }
  _projData->point2mesh(_workp, _workv);
  for (int i=0; i<_projData->getApexNumber(); i++) _grad[i] += _workv[i];
  
  // Contribution of the Seismic

  if (_flagSeismic)
  {
    _contributeSeismicDerivative(lambda);
    for (int i=0; i<_projData->getApexNumber(); i++) _grad[i] += _workv[i];
  }

  // Evaluate the norm of the gradient
  (*normgrad) = 0.;
  for (int i=0; i<_projData->getApexNumber(); i++)
    (*normgrad) += _grad[i] * _grad[i];
}

/*****************************************************************************/
/*!
**  Internal function to evaluate the seismic contribution to cost function
**
** \param[in]  lambda     Array of input values
**
*****************************************************************************/
void OptimCostBinary::_contributeSeismic(const VectorDouble& lambda)
{
  for (int i=0; i<_projSeis->getApexNumber(); i++) 
    _workv[i] = law_cdf_gaussian(lambda[i]);
  _projSeis->mesh2point(_workv, _works);
  for (int i=0; i<_projSeis->getPointNumber(); i++) 
    _works[i] -= _propSeis[i];
}

/*****************************************************************************/
/*!
**  Internal function to evaluate the seismic contribution to gradient function
**
** \param[in]  lambda     Array of input values
**
*****************************************************************************/
void OptimCostBinary::_contributeSeismicDerivative(const VectorDouble& lambda)
{
  _contributeSeismic(lambda);

  for (int i=0; i<_projSeis->getPointNumber(); i++)
    _works[i] *= _varSeis[i];
  _projSeis->point2mesh(_works, _workv);
  for (int i=0; i<_projSeis->getApexNumber(); i++) 
    _workv[i] *= law_df_gaussian(lambda[i]);
}
