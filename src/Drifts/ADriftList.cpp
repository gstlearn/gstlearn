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
#include "Space/ASpace.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Drifts/ADriftList.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"

ADriftList::ADriftList(bool flagLinked, const ASpace* space)
    : ADrift(space),
      _flagLinked(flagLinked),
      _driftEquationNumber(0),
      _coefDrift(),
      _drifts(),
      _filtered()
{
  _setDriftEquationNumber();
}

ADriftList::ADriftList(const ADriftList &r)
    : ADrift(r),
      _flagLinked(r._flagLinked),
      _driftEquationNumber(r._driftEquationNumber),
      _coefDrift(r._coefDrift),
      _drifts(),
      _filtered(r._filtered)
{
  _setDriftEquationNumber();
  for (auto e: r._drifts)
  {
    _drifts.push_back(dynamic_cast<ADriftElem*>(e->clone()));
  }
}

ADriftList& ADriftList::operator=(const ADriftList &r)
{
  if (this != &r)
  {
    ADrift::operator=(r);
    _flagLinked = r._flagLinked;
    _driftEquationNumber = r._driftEquationNumber;
    _coefDrift  = r._coefDrift;
    for (auto e: r._drifts)
    {
      _drifts.push_back(dynamic_cast<ADriftElem*>(e->clone()));
    }
    _filtered = r._filtered;
  }
  return *this;
}

ADriftList::~ADriftList()
{
  delAllDrift();
}

IClonable* ADriftList::clone() const
{
  return new ADriftList(*this);
}

bool ADriftList::isConsistent(const ASpace* /*space*/) const
{
  return true;
}

String ADriftList::toString(int /*level*/) const
{
  std::stringstream sstr;
  for (int i = 0; i < (int) getDriftNumber(); i++)
  {
    sstr << _drifts[i]->toString();
    if (_filtered[i])
      sstr << " (filtered)";
    sstr << std::endl;
  }
  return sstr.str();
}

double ADriftList::eval(const Db* /*db*/, int /*iech1*/) const
{
  double drift = 0;
  return drift;
}

void ADriftList::addDrift(const ADriftElem* drift)
{
  _drifts.push_back(dynamic_cast<ADriftElem*>(drift->clone()));
  _filtered.push_back(false);
  _updateCoefDrift();
}

void ADriftList::delDrift(unsigned int i)
{
  if (! _isDriftIndexValid(i)) return;
  _drifts.erase(_drifts.begin() + i);
  _filtered.erase(_filtered.begin() + i);
  _updateCoefDrift();
}

void ADriftList::delAllDrift()
{
  for (auto e: _drifts)
  {
    delete e;
  }
  _drifts.clear();
  _filtered.clear();
  _coefDrift.clear();
}

bool ADriftList::isFiltered(int i) const
{
  if (! _isDriftIndexValid(i)) return false;
  return _filtered[i];
}

void ADriftList::setFiltered(int i, bool filter)
{
  if (! _isDriftIndexValid(i)) return;
  _filtered[i] = filter;
}

const ADriftElem* ADriftList::getDrift(int il) const
{
  if (! _isDriftIndexValid(il)) return nullptr;
  return _drifts[il];
}

ADriftElem* ADriftList::getDrift(int il)
{
  if (! _isDriftIndexValid(il)) return nullptr;
  return _drifts[il];
}

const EDrift& ADriftList::getType(int il) const
{
  if (! _isDriftIndexValid(il)) return EDrift::UNKNOWN;
  return _drifts[il]->getType();
}

int ADriftList::getRankFex(int il) const
{
  if (! _isDriftIndexValid(il)) return 0;
  return _drifts[il]->getRankFex();
}

String ADriftList::getDriftName(int il) const
{
  if (! _isDriftIndexValid(il)) return String();
  return _drifts[il]->getDriftName();
}

void ADriftList::setType(int il, const EDrift& type)
{
  if (! _isDriftIndexValid(il)) return;
  _drifts[il]->setType(type);
}

void ADriftList::_setDriftEquationNumber()
{
  int ndrift = getDriftNumber();
  _driftEquationNumber = (_flagLinked) ? ndrift : ndrift * getNVariables();
}

int ADriftList::getNVariables() const
{
  if (getDriftNumber() > 0)
    return _drifts[0]->getNVariables();
  return 0;
}

bool ADriftList::_isDriftIndexValid(int i) const
{
  if (i < 0 || i > getDriftNumber())
  {
    mesArg("Drift Rank",i,getDriftNumber());
    return false;
  }
  return true;
}

bool ADriftList::_isDriftEquationValid(int ib) const
{
  if (ib < 0 || ib > getDriftEquationNumber())
  {
    mesArg("Drift Equation",ib,getDriftEquationNumber());
    return false;
  }
  return true;
}

void ADriftList::_updateCoefDrift()
{
  _setDriftEquationNumber();
  int nvar = getNVariables();
  int nfeq = getDriftEquationNumber();
  int nbfl = getDriftNumber();
  _coefDrift.resize(nvar * nfeq * nbfl);

  /* Copy the coefficients from the old to the new structure */

  if (_flagLinked)
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int ib = 0; ib < nfeq; ib++)
        for (int il = 0; il < nbfl; il++)
        {
          setCoefDrift(ivar, il, ib, (ib == il));
        }
  }
  else
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
        for (int jl = 0; jl < nbfl; jl++)
          for (int il = 0; il < nbfl; il++)
          {
            int ib = jvar + nvar * jl;
            setCoefDrift(ivar, il, ib, (ivar == jvar && il == jl));
          }
  }
}

VectorDouble ADriftList::getDrift(const Db* db, int ib, bool useSel)
{
  if (! _isDriftIndexValid(ib)) return VectorDouble();

  int nech = db->getActiveSampleNumber();
  VectorDouble vec(nech);

  int ecr = 0;
  for (int iech=0; iech<db->getSampleNumber(); iech++)
  {
    if (useSel && ! db->isActive(iech)) continue;
    vec[ecr++] = _drifts[ib]->eval(db, iech);
  }
  return vec;
}

VectorVectorDouble ADriftList::getDrifts(const Db* db, bool useSel)
{
  int ndrift = static_cast<int>(_drifts.size());

  VectorVectorDouble vec;
  for (int ib=0; ib<ndrift; ib++)
  {
    vec.push_back(getDrift(db, ib, useSel));
  }
  return vec;
}
