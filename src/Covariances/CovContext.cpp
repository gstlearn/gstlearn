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

#include "MatrixC/MatrixCSSym.hpp"
#include "Space/ASpace.hpp"
#include "Basic/Vector.hpp"
#include "geoslib_f.h"

CovContext::CovContext(int nvar,
                       int irfMaxDegree,
                       double field,
                       const ASpace* space)
: _nVar(nvar),
  _irfMaxDegree(irfMaxDegree),
  _field(field),
  _ballRadius(0.),
  _mean(),
  _covar0(),
  _space(space)
{
  _update();
  if (_space == nullptr)
    _space = ASpaceObject::getGlobalSpace();
}

CovContext::CovContext(const Db *db,
                       int irfMaxDegree,
                       const ASpace* space)
: _nVar(0),
  _irfMaxDegree(irfMaxDegree),
  _field(0.),
  _ballRadius(0.),
  _mean(),
  _covar0(),
  _space(space)
{
  _nVar = db->getVariableNumber();
  // As it does not make sense not to have any variable, this number is set to 1 at least
  if (_nVar <= 1) _nVar = 1;
  _field = db->getFieldSize();
  _update();
  if (_space == nullptr)
    _space = ASpaceObject::getGlobalSpace();
}


CovContext::CovContext(const CovContext &r)
: _nVar(r._nVar),
  _irfMaxDegree(r._irfMaxDegree),
  _field(r._field),
  _ballRadius(r._ballRadius),
  _mean(r._mean),
  _covar0(r._covar0),
  _space(r._space)
  {
}

CovContext& CovContext::operator=(const CovContext &r)
{
  if (this != &r)
  {
    _nVar = r._nVar;
    _irfMaxDegree = r._irfMaxDegree;
    _field = r._field;
    _ballRadius = r._ballRadius;
    _mean = r._mean;
    _covar0 = r._covar0;
    _space = r._space;
  }
  return *this;
}

CovContext::~CovContext()
{
}

std::string CovContext::toString(int level) const
{
  std::stringstream sstr;
  sstr << _space->toString() << std::endl;
  sstr << "Nb Variables       = "       << _nVar << std::endl;
  sstr << "Maximum IRF Degree = "       << _irfMaxDegree << std::endl;
  sstr << "Field Size         = "       << _field << std::endl;
  sstr << "Ball Radius for Gradient = " << _ballRadius << std::endl;
  sstr << "Mean(s)            = "       << ut_vector_string(_mean);
  sstr << "Covariance (0)     = "       << ut_vector_string(_covar0);
  return sstr.str();
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

void CovContext::_update()
{
  _mean.resize(_nVar, 0.);
  MatrixCSSym Id = MatrixCSSym(_nVar);
  Id.setIdentity();
  _covar0 = Id.getValues();
}
