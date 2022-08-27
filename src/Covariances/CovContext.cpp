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

/**
 * Create a covariances context giving the number dimensions of a predefined space RN
 *
 * @param nvar         Number of variables
 * @param space        Space definition
 * @param field        Maximum field distance (used for covariances having no sill)
 */
CovContext::CovContext(int nvar, const ASpace *space, double field)

    : ASpaceObject(space),
      _nVar(nvar),
      _field(field),
      _mean(),
      _covar0()
{
  _update();
}

/**
 * Create a covariances context giving the number dimensions of a predefined space RN
 *
 * @param nvar         Number of variables
 * @param ndim         Number of dimension of the euclidean space (RN)
 * @param field        Maximum field distance (used for covariances having no sill)
 * @param mean         Vector of Means
 * @param covar0       Vector of variance-covariance
 */
CovContext::CovContext(int nvar,
                       int ndim,
                       double field,
                       const VectorDouble &mean,
                       const VectorDouble &covar0)
    : ASpaceObject(SpaceRN(ndim)),
      _nVar(nvar),
      _field(field),
      _mean(mean),
      _covar0(covar0)
{
  _update();
}

CovContext::CovContext(const Db *db, const ASpace* space)
    : ASpaceObject(space),
      _nVar(0),
      _field(1.),
      _mean(),
      _covar0()
{
  /// TODO : check Db dimension vs provided space
  _nVar = db->getVariableNumber();
  // As it does not make sense not to have any variable, this number is set to 1 at least
  if (_nVar <= 1) _nVar = 1;
  _field = db->getColumnSize();
  _update();
}

CovContext::CovContext(const Vario *vario, const ASpace *space)
    : ASpaceObject(space),
      _nVar(0),
      _field(1.),
      _mean(),
      _covar0()
{
  /// TODO : check vario dimension vs provided space
  _nVar = vario->getVariableNumber();
  _field = vario->getHmax();
  _update();
}

CovContext::CovContext(const CovContext &r)
    : ASpaceObject(r),
      _nVar(r._nVar),
      _field(r._field),
      _mean(r._mean),
      _covar0(r._covar0)
{
}

CovContext& CovContext::operator=(const CovContext &r)
{
  if (this != &r)
  {
    ASpaceObject::operator =(r);
    _nVar = r._nVar;
    _field = r._field;
    _mean = r._mean;
    _covar0 = r._covar0;
  }
  return *this;
}

CovContext::~CovContext()
{
}

String CovContext::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  sstr << ASpaceObject::toString(strfmt);
  sstr << "Nb Variables       = "       << _nVar << std::endl;
  sstr << "Field Size         = "       << _field << std::endl;
  sstr << "Mean(s)            = "       << ut_vector_string(_mean);
  sstr << "Covariance (0)     = "       << ut_vector_string(_covar0);
  return sstr.str();
}

CovContext* CovContext::create(int nvar, int ndim, double field)
{
  CovContext* ctxt = new CovContext(nvar, ndim, field);
  return ctxt;
}

bool CovContext::isConsistent(const ASpace* space) const
{
  /// TODO: Consistency of CovContext toward a space: Possible duplicate:
  /// - CovFactory::_isValid
  /// - ACovFunc::isConsistent
  return (_space->isEqual(space));
}

/**
 * Checks that two CovContext are 'similar'
 * @param r Secondary CovContext to be compared with this
 * @return
 */
bool CovContext::isEqual(const CovContext &r) const
{
  return (_nVar == r.getNVar()                 &&
          _space->isEqual(r.getSpace())
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
  if (_mean.empty())
    _mean.resize(_nVar, 0.);
  if (_covar0.empty())
  {
    MatrixSquareSymmetric Id(_nVar);
    Id.setIdentity();
    _covar0 = Id.getValues();
  }
}

/**
 * This operation sets the contents of the current CovContext class
 * by copying the information from a source CovContext
 * @param ctxt Source CovContext
 *
 * @remark: This operation does not allow changing the number of variables
 */
void CovContext::copyCovContext(const CovContext& ctxt)
{
  if (ctxt._nVar != _nVar)
  {
    messerr("The update of a CovContext does not allow modifying");
    messerr("the number of variables");
    messerr("Operation is cancelled");
    return;
  }
  _field  = ctxt._field;
  _mean   = ctxt._mean;
  _covar0 = ctxt._covar0;
}
