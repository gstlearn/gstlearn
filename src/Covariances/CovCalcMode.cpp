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
                         int keepOnlyCovIdx,
                         bool unitary,
                         int orderVario)
: AStringable(),
  _factorySettings(true),
  _member(member),
  _asVario(asVario),
  _keepOnlyCovIdx(keepOnlyCovIdx),
  _unitary(unitary),
  _orderVario(orderVario),
  _covFiltered()
{
  _checkFactorySettings();
}

CovCalcMode::CovCalcMode(const CovCalcMode &r)
    : AStringable(r),
      _factorySettings(r._factorySettings),
      _member(r._member),
      _asVario(r._asVario),
      _keepOnlyCovIdx(r._keepOnlyCovIdx),
      _unitary(r._unitary),
      _orderVario(r._orderVario),
      _covFiltered(r._covFiltered)
{
  _checkFactorySettings();
}

CovCalcMode& CovCalcMode::operator=(const CovCalcMode &r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _factorySettings = r._factorySettings;
    _member = r._member;
    _asVario = r._asVario;
    _keepOnlyCovIdx = r._keepOnlyCovIdx;
    _unitary = r._unitary;
    _orderVario = r._orderVario;
    _covFiltered = r._covFiltered;
    _checkFactorySettings();
  }
  return *this;
}
CovCalcMode::~CovCalcMode()
{
}

CovCalcMode* CovCalcMode::create(const ECalcMember &member,
                                 bool asVario,
                                 int keepOnlyCovIdx,
                                 bool unitary,
                                 int orderVario)
{
  return new CovCalcMode(member, asVario, keepOnlyCovIdx, unitary, orderVario);
}


void CovCalcMode::_checkFactorySettings(const ECalcMember& member,
                                        bool asVario,
                                        int keepOnlyCovIdx,
                                        bool unitary,
                                        int orderVario)
{
  _factorySettings = false;
  if (_member != member) return;
  if (_asVario != asVario) return;
  if (_keepOnlyCovIdx != (int) keepOnlyCovIdx) return;
  if (_unitary != unitary) return;
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
