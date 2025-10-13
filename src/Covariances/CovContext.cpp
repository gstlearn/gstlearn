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
#include "Covariances/CovContext.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/VectorNumT.hpp"
#include "Db/Db.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpaceRN.hpp"
#include "Variogram/Vario.hpp"

namespace gstlrn
{
/**
 * Create a covariances context giving the number dimensions of a predefined space RN
 *
 * @param nvar         Number of variables
 * @param space        Space definition
 */
CovContext::CovContext(Id nvar, const ASpaceSharedPtr& space)

  : ASpaceObject(space)
  , _nVar(nvar)
  , _field(TEST)
  , _covar0()
{
  _update();
}

/**
 * Create a covariances context giving the number dimensions of a predefined space RN
 *
 * @param nvar         Number of variables
 * @param ndim         Number of dimension of the euclidean space (RN)
 * @param covar0       Vector of variance-covariance
 */
CovContext::CovContext(Id nvar,
                       Id ndim,
                       const VectorDouble& covar0)
  : ASpaceObject(SpaceRN::create(ndim))
  , _nVar(nvar)
  , _field(TEST)
  , _covar0(covar0)
{
  _update();
}

CovContext::CovContext(const Db* db, const ASpaceSharedPtr& space)
  : ASpaceObject(space)
  , _nVar(0)
  , _field(TEST)
  , _covar0()
{
  /// TODO : check Db dimension vs provided space
  _nVar = db->getNLoc(ELoc::Z);
  // As it does not make sense not to have any variable, this number is set to 1 at least
  if (_nVar <= 1) _nVar = 1;
  _update();
}

CovContext::CovContext(const Vario* vario)
  : ASpaceObject(vario->getSpace())
  , _nVar(0)
  , _field(TEST)
  , _covar0()
{
  /// TODO : check vario dimension vs provided space
  _nVar  = vario->getNVar();
  _field = vario->getHmax();
  _update();
}

CovContext::CovContext(const CovContext& r)
  : ASpaceObject(r)
  , _nVar(r._nVar)
  , _field(r._field)
  , _covar0(r._covar0)
{
}

CovContext& CovContext::operator=(const CovContext& r)
{
  if (this != &r)
  {
    ASpaceObject::operator=(r);
    _nVar   = r._nVar;
    _field  = r._field;
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
  sstr << "Nb Variables       = " << _nVar << std::endl;
  if (!FFFF(_field))
    sstr << "Field Size         = " << _field << std::endl;
  sstr << "Covariance (0)     = " << VH::toStringAsVD(_covar0);
  return sstr.str();
}

CovContext* CovContext::create(Id nvar, Id ndim)
{
  auto* ctxt = new CovContext(nvar, ndim);
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
bool CovContext::isEqual(const CovContext& r) const
{
  return (_nVar == r.getNVar() && _space->isEqual(r.getSpace().get()));
}

double CovContext::getCovar0(Id ivar, Id jvar) const
{
  auto rank = _getIndex(ivar, jvar);
  if (rank < 0 || rank >= static_cast<Id>(_covar0.size()))
    my_throw("Invalid argument in _setCovar0");
  return _covar0[rank];
}

/**
 * Define the covariance at the origin
 * @param covar0 Values
 */
void CovContext::setCovar0s(const VectorDouble& covar0)
{
  if (_covar0.size() == covar0.size())
    _covar0 = covar0;
}

void CovContext::setCovar0(Id ivar, Id jvar, double covar0)
{
  auto rank = _getIndex(ivar, jvar);
  if (rank < 0 || rank >= static_cast<Id>(_covar0.size()))
    my_throw("Invalid argument in _setCovar0");
  _covar0[rank] = covar0;
}

Id CovContext::_getIndex(Id ivar, Id jvar) const
{
  return ivar * getNVar() + jvar;
}

void CovContext::_update()
{
  if (_nVar * _nVar != static_cast<Id>(_covar0.size()))
  {
    MatrixSymmetric Id(_nVar);
    Id.setIdentity();
    _covar0 = Id.getValues();
  }
}

/**
 * This operation sets the contents of the current CovContext class
 * by copying the information from a source CovContext
 * @param ctxt Source CovContext
 * @param severe When severe, does not allow modifying the number of variables
 */
void CovContext::copyCovContext(const CovContext& ctxt, bool severe)
{
  if (severe && ctxt._nVar != _nVar)
  {
    messerr("The update of a CovContext does not allow modifying the number of variables (old=%d -> New=%d)",
            _nVar, ctxt._nVar);
    messerr("Operation is cancelled");
    return;
  }
  _nVar   = ctxt._nVar;
  _field  = ctxt._field;
  _covar0 = ctxt._covar0;
  _update();
}

const CovContext* CovContext::createReduce(const VectorInt& validVars) const
{
  Id ecr, lec;
  Id nvar   = static_cast<Id>(validVars.size());
  auto ndim = getNDim();
  VectorBool valids(_nVar, false);
  for (Id ivar = 0; ivar < nvar; ivar++) valids[validVars[ivar]] = true;

  VectorDouble covar0(nvar * nvar, 0);
  ecr = 0;
  lec = 0;
  for (Id ivar = 0; ivar < _nVar; ivar++)
    for (Id jvar = 0; jvar < _nVar; jvar++)
    {
      if (valids[ivar] && valids[jvar])
      {
        covar0[ecr++] = _covar0[lec];
      }
      lec++;
    }

  auto* newctxt = new CovContext(nvar, static_cast<Id>(ndim), covar0);
  return newctxt;
}
} // namespace gstlrn