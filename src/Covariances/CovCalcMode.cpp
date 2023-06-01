/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "Covariances/CovCalcMode.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"

CovCalcMode::CovCalcMode(const ECalcMember& member,
                         bool asVario,
                         bool normalized,
                         bool filterNugget,
                         int keepOnlyCovIdx,
                         bool unitary,
                         int envelop,
                         int orderVario)
: AStringable(),
  _factorySettings(true),
  _member(member),
  _asVario(asVario),
  _normalized(normalized),
  _filterNugget(filterNugget),
  _keepOnlyCovIdx(keepOnlyCovIdx),
  _unitary(unitary),
  _envelop(envelop),
  _orderVario(orderVario),
  _indexClass(0),
  _covFiltered()
{
}

CovCalcMode::CovCalcMode(const CovCalcMode &r)
    : AStringable(r),
      _factorySettings(r._factorySettings),
      _member(r._member),
      _asVario(r._asVario),
      _normalized(r._normalized),
      _filterNugget(r._filterNugget),
      _keepOnlyCovIdx(r._keepOnlyCovIdx),
      _unitary(r._unitary),
      _envelop(r._envelop),
      _orderVario(r._orderVario),
      _indexClass(r._indexClass),
      _covFiltered(r._covFiltered)
{

}

CovCalcMode& CovCalcMode::operator=(const CovCalcMode &r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _factorySettings = r._factorySettings;
    _member = r._member;
    _asVario = r._asVario;
    _normalized = r._normalized;
    _filterNugget = r._filterNugget;
    _keepOnlyCovIdx = r._keepOnlyCovIdx;
    _unitary = r._unitary;
    _envelop = r._envelop;
    _orderVario = r._orderVario;
    _indexClass = r._indexClass;
    _covFiltered = r._covFiltered;
  }
  return *this;
}
CovCalcMode::~CovCalcMode()
{
}

CovCalcMode* CovCalcMode::create(const ECalcMember &member,
                                 bool asVario,
                                 bool normalized,
                                 bool filterNugget,
                                 int keepOnlyCovIdx,
                                 bool unitary,
                                 int envelop,
                                 int orderVario)
{
  return new CovCalcMode(member, asVario, normalized, filterNugget, keepOnlyCovIdx, unitary, envelop, orderVario);
}


void CovCalcMode::_checkFactorySettings(const ECalcMember& member,
                                        bool asVario,
                                        bool normalized,
                                        bool filterNugget,
                                        int keepOnlyCovIdx,
                                        bool unitary,
                                        int envelop,
                                        int orderVario)
{
  _factorySettings = false;
  if (_member != member) return;
  if (_asVario != asVario) return;
  if (_normalized != normalized) return;
  if (_filterNugget != filterNugget) return;
  if (_keepOnlyCovIdx != (int) keepOnlyCovIdx) return;
  if (_unitary != unitary) return;
  if (_envelop != envelop) return;
  if (_orderVario != orderVario) return;
  _factorySettings = true;
}

bool CovCalcMode::getCovFiltered(int i) const
{
  if (_covFiltered.empty()) return false;
  if (i < 0 || i >= (int) _covFiltered.size()) return false;
  return _covFiltered[i];
}

void CovCalcMode::setCovFiltered(int i, bool status)
{
  if (_covFiltered.empty()) return;
  if (i < 0 || i >= (int) _covFiltered.size()) return;
  _covFiltered[i] = status;
  _checkFactorySettings();
}

void CovCalcMode::setAllCovFiltered(int ncov, bool status)
{
  if (_covFiltered.empty())
    _covFiltered.resize(ncov, status);
  else
  {
    for (int i=0; i<(int) _covFiltered.size(); i++)
      _covFiltered[i] = status;
  }
  _checkFactorySettings();
}
