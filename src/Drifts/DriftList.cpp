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
#include "Drifts/DriftList.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"

DriftList::DriftList(bool flagLinked, const ASpace* space)
    : ADrift(space),
      _flagLinked(flagLinked),
      _coefDrift(),
      _drifts(),
      _filtered()
{
}

DriftList::DriftList(const DriftList &r)
    : ADrift(r),
      _flagLinked(r._flagLinked),
      _coefDrift(r._coefDrift),
      _drifts(),
      _filtered(r._filtered)
{
  for (auto e: r._drifts)
  {
    _drifts.push_back(dynamic_cast<ADriftElem*>(e->clone()));
  }
}

DriftList& DriftList::operator=(const DriftList &r)
{
  if (this != &r)
  {
    ADrift::operator=(r);
    _flagLinked = r._flagLinked;
    _coefDrift  = r._coefDrift;
    for (auto e: r._drifts)
    {
      _drifts.push_back(dynamic_cast<ADriftElem*>(e->clone()));
    }
    _filtered = r._filtered;
  }
  return *this;
}

DriftList::~DriftList()
{
  delAllDrift();
}

bool DriftList::isConsistent(const ASpace* /*space*/) const
{
  return true;
}

String DriftList::toString(int /*level*/) const
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

double DriftList::eval(const Db* /*db*/, int /*iech1*/) const
{
  double drift = 0;
  return drift;
}

void DriftList::addDriftList(const DriftList* drifts)
{
  int ndrift = drifts->getDriftNumber();
  _flagLinked = drifts->isFlagLinked();
  for (int idrift = 0; idrift < ndrift; idrift++)
    addDrift(drifts->getDrift(idrift));
}
void DriftList::addDrift(const ADriftElem* drift)
{
  _drifts.push_back(dynamic_cast<ADriftElem*>(drift->clone()));
  _filtered.push_back(false);
  _updateCoefDrift();
}

void DriftList::delDrift(unsigned int i)
{
  if (! _isDriftIndexValid(i)) return;
  _drifts.erase(_drifts.begin() + i);
  _filtered.erase(_filtered.begin() + i);
  _updateCoefDrift();
}

void DriftList::delAllDrift()
{
  for (auto e: _drifts)
  {
    delete e;
  }
  _drifts.clear();
  _filtered.clear();
  _coefDrift.clear();
}

bool DriftList::isFiltered(int i) const
{
  if (! _isDriftIndexValid(i)) return false;
  return _filtered[i];
}

void DriftList::setFiltered(int i, bool filter)
{
  if (! _isDriftIndexValid(i)) return;
  _filtered[i] = filter;
}

const ADriftElem* DriftList::getDrift(int il) const
{
  if (! _isDriftIndexValid(il)) return nullptr;
  return _drifts[il];
}

ADriftElem* DriftList::getDrift(int il)
{
  if (! _isDriftIndexValid(il)) return nullptr;
  return _drifts[il];
}

const EDrift& DriftList::getType(int il) const
{
  if (! _isDriftIndexValid(il)) return EDrift::UNKNOWN;
  return _drifts[il]->getType();
}

int DriftList::getRankFex(int il) const
{
  if (! _isDriftIndexValid(il)) return 0;
  return _drifts[il]->getRankFex();
}

String DriftList::getDriftName(int il) const
{
  if (! _isDriftIndexValid(il)) return String();
  return _drifts[il]->getDriftName();
}

void DriftList::setType(int il, const EDrift& type)
{
  if (! _isDriftIndexValid(il)) return;
  _drifts[il]->setType(type);
}

int DriftList::getDriftEquationNumber() const
{
  int ndrift = getDriftNumber();
  int ndriftEquationNumber = (_flagLinked) ? ndrift : ndrift * getNVariables();
  return ndriftEquationNumber;
}

int DriftList::getNVariables() const
{
  if (getDriftNumber() > 0)
    return _drifts[0]->getNVariables();
  return 0;
}

bool DriftList::_isDriftIndexValid(int i) const
{
  if (i < 0 || i > getDriftNumber())
  {
    mesArg("Drift Rank",i,getDriftNumber());
    return false;
  }
  return true;
}

bool DriftList::_isDriftEquationValid(int ib) const
{
  if (ib < 0 || ib > getDriftEquationNumber())
  {
    mesArg("Drift Equation",ib,getDriftEquationNumber());
    return false;
  }
  return true;
}

void DriftList::_updateCoefDrift()
{
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

VectorDouble DriftList::getDrift(const Db* db, int ib, bool useSel)
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

VectorVectorDouble DriftList::getDrifts(const Db* db, bool useSel)
{
  int ndrift = static_cast<int>(_drifts.size());

  VectorVectorDouble vec;
  for (int ib=0; ib<ndrift; ib++)
  {
    vec.push_back(getDrift(db, ib, useSel));
  }
  return vec;
}
