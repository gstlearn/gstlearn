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
#include "Covariances/ACovFunc.hpp"

#include "Model/Cova.hpp"

#include "geoslib_f.h"
#include "geoslib_f_private.h"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"

ACovFunc::ACovFunc(const ENUM_COVS& type, const CovContext& ctxt)
: _type(type),
  _ctxt(ctxt),
  _param(TEST)
{
  if (!isConsistent())
    my_throw ("Cannot create such covariance function in that context");
}

ACovFunc::ACovFunc(const ACovFunc &r)
: _type(r._type),
  _ctxt(r._ctxt),
  _param(r._param)
{
}

ACovFunc& ACovFunc::operator=(const ACovFunc &r)
{
  if (this != &r)
  {
    _type = r._type;
    _ctxt = r._ctxt;
    _param = r._param;
  }
  return *this;
}

ACovFunc::~ACovFunc()
{
}

void ACovFunc::setParam(double param)
{
  /// TODO : Do not throw in setter. Check range and build the error message here.
  double max = getParMax();
  if (param < 0. || (!FFFF(max) && param > max))
    my_throw("Wrong third parameter value");
  _param = param;
}

void ACovFunc::setField(double field)
{
  if (field < EPSILON10)
    my_throw("Cannot scale with zero");
  _ctxt.setField(field);
}

double ACovFunc::evalCov(double h) const
{
  return _evaluateCov(h);
}
double ACovFunc::evalCovDerivative(int degree, double h) const
{
  return _evaluateCovDerivate(degree, h);
}

VectorDouble ACovFunc::evalCovVec(const VectorDouble& vech) const
{
  VectorDouble vec;
  for (const auto& h : vech)
    vec.push_back(evalCov(h));
  return vec;
}
VectorDouble ACovFunc::evalCovDerivativeVec(int degree,
                                            const VectorDouble& vech) const
{
  VectorDouble vec;
  for (const auto& i : vech)
    vec.push_back(evalCovDerivative(degree, i));
  return vec;
}
std::string ACovFunc::toString(int level) const
{
  std::stringstream sstr;
  sstr << getCovName();
  if (hasParam())
    sstr << " (Third Parameter = " << getParam() << ")";
  sstr << std::endl;
  return sstr.str();
}

/// Test consistency toward the current context
bool ACovFunc::isConsistent() const
{
  unsigned int maxndim = getMaxNDim();
  if ((maxndim > 0 && (maxndim < _ctxt.getNDim())))
    /// TODO : Test irfDegree vs getMinOrder in CovElem because zonal anisotropies
       return false;
  return true;
}

bool ACovFunc::hasInt1D() const
{
  return (getMaxNDim() >= 1 && getMinOrder() <= 0);
}
bool ACovFunc::hasInt2D() const
{
  return (getMaxNDim() >= 2 && getMinOrder() <= 0);
}
 /**
  * Calculate covariance derivatives, i.e.
  * Degree 1: C^1(r) / r
  * degree 2: C^2(r)
  * Degree 3: C^3(r)
  * Degree 4: C^4(r)
  * @param degree Level of derivation
  * @param h Normalized distance
  * @return
  */
double ACovFunc::_evaluateCovDerivate(int degree, double h) const
{
  my_throw("Undefined derivative for this covariance");
  return 0.;
}
