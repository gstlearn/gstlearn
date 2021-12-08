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
#include "Covariances/CovContext.hpp"

#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Space/ASpace.hpp"
#include "Basic/Vector.hpp"
#include "Variogram/Vario.hpp"
#include "Space/SpaceRN.hpp"
#include "Db/Db.hpp"
#include "geoslib_f.h"

CovContext::CovContext(int nvar,
                       const ASpace* space,
                       int irfMaxDegree,
                       double field)

: ASpaceObject(space)
,  _nVar(nvar)
,  _irfMaxDegree(irfMaxDegree)
,  _field(field)
,  _ballRadius(0.)
,  _mean()
,  _covar0()
{
  _update();
}

CovContext::CovContext(int nvar,
                       int ndim,
                       int irfMaxDegree,
                       double field)
: ASpaceObject(SpaceRN(ndim))
,  _nVar(nvar)
,  _irfMaxDegree(irfMaxDegree)
,  _field(field)
,  _ballRadius(0.)
,  _mean()
,  _covar0()
{
  _update();
}

CovContext::CovContext(const Db *db, int irfMaxDegree, const ASpace* space)
: ASpaceObject(space)
, _nVar(0)
, _irfMaxDegree(irfMaxDegree)
, _field(0.)
, _ballRadius(0.)
, _mean()
, _covar0()
{
  /// TODO : check Db dimension vs provided space
  _nVar = db->getVariableNumber();
  // As it does not make sense not to have any variable, this number is set to 1 at least
  if (_nVar <= 1) _nVar = 1;
  _field = db->getFieldSize();
  _update();
}

CovContext::CovContext(const Vario* vario, int irfMaxDegree, const ASpace* space)
: ASpaceObject(space)
, _nVar(0)
, _irfMaxDegree(irfMaxDegree)
, _field(0.)
, _ballRadius(0.)
, _mean()
, _covar0()
{
  /// TODO : check vario dimension vs provided space
  _nVar = vario->getVariableNumber();
  _field = vario->getHmax();
  _update();
}

CovContext::CovContext(const CovContext &r)
: ASpaceObject(r)
, _nVar(r._nVar)
, _irfMaxDegree(r._irfMaxDegree)
, _field(r._field)
, _ballRadius(r._ballRadius)
, _mean(r._mean)
, _covar0(r._covar0)
{
}

CovContext& CovContext::operator=(const CovContext &r)
{
  if (this != &r)
  {
    ASpaceObject::operator =(r);
    _nVar = r._nVar;
    _irfMaxDegree = r._irfMaxDegree;
    _field = r._field;
    _ballRadius = r._ballRadius;
    _mean = r._mean;
    _covar0 = r._covar0;
  }
  return *this;
}

CovContext::~CovContext()
{
}

String CovContext::toString(int level) const
{
  std::stringstream sstr;
  sstr << ASpaceObject::toString(level);
  sstr << "Nb Variables       = "       << _nVar << std::endl;
  sstr << "Maximum IRF Degree = "       << _irfMaxDegree << std::endl;
  sstr << "Field Size         = "       << _field << std::endl;
  sstr << "Ball Radius for Gradient = " << _ballRadius << std::endl;
  sstr << "Mean(s)            = "       << ut_vector_string(_mean);
  sstr << "Covariance (0)     = "       << ut_vector_string(_covar0);
  return sstr.str();
}

bool CovContext::isConsistent(const ASpace* /*space*/) const
{
  /// TODO: Consistency of CovContext toward a space: Possible duplicate:
  /// - CovFatory::_isValid
  /// - ACovFunc::isConsistent
  return true;
}

bool CovContext::isEqual(const CovContext &r) const
{
  return (_nVar == r.getNVar()                 &&
          _irfMaxDegree == r.getIrfMaxDegree() &&
          _field == r.getField()               &&
          _space->isEqual(r.getSpace())        &&
          _ballRadius == r._ballRadius         &&
          ut_vector_same(_mean, r._mean)       &&
          ut_vector_same(_covar0, r._covar0)
          );
}

double CovContext::getMean(int ivar) const
{
  if (ivar < 0 || ivar >= (int) _mean.size())
    throw("Invalid argument in _getMean");
  return _mean[ivar];
}

double CovContext::getCovar0(int ivar, int jvar) const
{
  int rank = _getIndex(ivar, jvar);
  if (rank < 0 || rank >= (int) _covar0.size())
    throw("Invalid argument in _setCovar0");
  return _covar0[rank];
}

void CovContext::setMean(const VectorDouble& mean)
{
  if (_mean.size() == mean.size())
    _mean = mean;
}

/**
 * Define the Mean for one variable
 * @param ivar Rank of the variable (starting from 0)
 * @param mean Value for the mean
 */
void CovContext::setMean(int ivar, const double mean)
{
  if (ivar < 0 || ivar >= (int) _mean.size())
    throw("Invalid argument in _setMean");
  _mean[ivar] = mean;
}

/**
 * Define the covariance at the origin
 * @param covar0 Values
 */
void CovContext::setCovar0(const VectorDouble& covar0)
{
  if (_covar0.size() == covar0.size())
    _covar0 = covar0;
}

void CovContext::setCovar0(int ivar, int jvar, double covar0)
{
  int rank = _getIndex(ivar, jvar);
  if (rank < 0 || rank >= (int) _covar0.size())
    throw("Invalid argument in _setCovar0");
  _covar0[rank] = covar0;
}

int CovContext::_getIndex(int ivar, int jvar) const
{
  return ivar * getNVar() + jvar;
}

void CovContext::_update()
{
  _mean.resize(_nVar, 0.);
  MatrixSquareSymmetric Id(_nVar);
  Id.setIdentity();
  _covar0 = Id.getValues();
}
