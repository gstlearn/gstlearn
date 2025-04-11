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

#include "Covariances/CorMatern.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/CorAniso.hpp"
#include "Covariances/CovContext.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpaceComposite.hpp"
#include "Space/SpacePoint.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "geoslib_define.h"
#include "Basic/MathFunc.hpp"

#include <cmath>
#include <vector>

CorMatern::CorMatern(const VectorDouble &ranges,
                     const VectorDouble &angles,
                     const VectorDouble& coeffScales, 
                     const VectorDouble& params ,
                     bool flagRange)
  : ACov()
  , _nVar(params.size())
  , _corRef(CorAniso::createAnisotropic(CovContext(1, ranges.size()),ECov::MATERN, ranges, params[0], angles, flagRange))
  , _corMatern(*_corRef)
  , _coeffScales(coeffScales)
  , _params(params)
  , _corMax(_nVar)
  , _angles(angles)
{
  if ((int)_coeffScales.size() != _nVar - 1)
  {
    messerr("CorMatern: inconsistent size between coeffScales and params");
    messerr("CorMatern: coeffScales size = %d, params size = %d", _coeffScales.size(), _nVar);
    _nVar = 0;
    _corMax = MatrixSymmetric(0);
    _coeffScales = VectorDouble();
    _params = VectorDouble();
    return;
  }
  _coeffScales.push_front(1.);
  setContext(_corMatern.getContext());
  _ctxt.setNVar(_nVar);
  for (int ivar = 0; ivar < _nVar; ivar++)
  {
    _corMax.setValue(ivar, ivar, 1.);
    for (int jvar = ivar + 1; jvar < _nVar; jvar++)
    {
      double scalei = _coeffScales[ivar];
      double scalej = _coeffScales[jvar];
      double scaleij = computeScale(ivar,jvar);
      double nui = _params[ivar];
      double nuj = _params[jvar];
      double nuij = computeParam(ivar,jvar);
      double gni =  exp(loggamma(nui));
      double gnj =  exp(loggamma(nuj));
      double gnij = exp(loggamma(nuij));
      double ratioi = gni / pow(scalei, 2. * nui);
      double ratioj = gnj / pow(scalej, 2. * nuj);
      double ratioij = gnij / pow(scaleij, 2. * nuij);
      double val = ratioij / sqrt(ratioi * ratioj);
      _corMax.setValue(ivar, jvar, val);
    }
  }
}

CorMatern::CorMatern(const CorMatern& r)
  : ACov(r)
  , _corMatern(r._corMatern)
{   
    _nVar = r._nVar;
    _corRef = r._corRef;
    _coeffScales = r._coeffScales;
    _params = r._params;
    _corMax = r._corMax;
    _angles = r._angles;
}

CorMatern& CorMatern::operator=(const CorMatern& r)
{
  if (this != &r)
  {
    ACov::operator=(r);
    _nVar = r._nVar;
    _corMatern = r._corMatern;
    _coeffScales = r._coeffScales;
    _params = r._params;
    _corMax = r._corMax;
    _angles = r._angles;
  }
  return *this;
}

double CorMatern::computeScale(int ivar, int jvar) const
{ 
    if (ivar == jvar)
    {
        return  _coeffScales[ivar];
    }
    double ci2 =  _coeffScales[ivar] * _coeffScales[ivar];
    double cj2 =  _coeffScales[jvar] * _coeffScales[jvar];
   
    return sqrt( 0.5 * (ci2+ cj2));
} 

double CorMatern::computeParam(int ivar, int jvar) const
{
  if (ivar == jvar)
  {
    return _params[ivar];
  }

  return 0.5 * (_params[ivar] + _params[jvar]);
} 

CorMatern::~CorMatern()
{
  delete _corRef;
}

void CorMatern::_optimizationSetTarget(SpacePoint& pt) const
{
  DECLARE_UNUSED(pt)
}

void CorMatern::_optimizationPreProcess(int mode, const std::vector<SpacePoint>& ps) const
{
  DECLARE_UNUSED(mode)
  DECLARE_UNUSED(ps)
  // _covS->_optimizationPreProcess(p);
  // _covTemp->_optimizationPreProcess(p);
}

void CorMatern::_optimizationPostProcess() const
{
  //_covS->optimizationPostProcess();
  //_covTemp->optimizationPostProcess();
}

double CorMatern::_eval(const SpacePoint& p1,
                        const SpacePoint& p2,
                        int ivar,
                        int jvar,
                        const CovCalcMode* mode) const
{
    _corMatern.setParam(computeParam(ivar, jvar));

    VectorDouble scales = _corRef->getScales();
  
    double correcScale = computeScale(ivar, jvar);
    for (int idim = 0; idim < (int)_space->getNDim(); idim++)
    {
        scales[idim] *= correcScale;
    }
    _corMatern.setRotationAnglesAndRadius(_angles,VectorDouble(),scales);
    return _corMax.getValue(ivar, jvar) * _corMatern.evalCov(p1, p2, 0 , 0 , mode);
}

