/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "Covariances/CovCalcMode.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"

CovCalcMode::CovCalcMode(const ECalcMember &member)
    : AStringable(),
      _member(member),
      _asVario(false),
      _unitary(false),
      _orderVario(0),
      _activeCovList()
{
}

CovCalcMode::CovCalcMode(const CovCalcMode &r)
    : AStringable(r),
      _member(r._member),
      _asVario(r._asVario),
      _unitary(r._unitary),
      _orderVario(r._orderVario),
      _activeCovList(r._activeCovList)
{
}

CovCalcMode& CovCalcMode::operator=(const CovCalcMode &r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _member = r._member;
    _asVario = r._asVario;
    _unitary = r._unitary;
    _orderVario = r._orderVario;
    _activeCovList = r._activeCovList;

  }
  return *this;
}
CovCalcMode::~CovCalcMode()
{
}

CovCalcMode* CovCalcMode::create(const ECalcMember &member)
{
  return new CovCalcMode(member);
}

void CovCalcMode::setActiveCovListFromOne(int keepOnlyCovIdx)
{
  _activeCovList.clear();
  if (keepOnlyCovIdx >= 0) _activeCovList.push_back(keepOnlyCovIdx);
}

/**
 * Set the list of active covariances from an intervam
 * @param inddeb Lower bound of the interval (included)
 * @param indto  Upper bound of the interval (excluded)
 */
void CovCalcMode::setActiveCovListFromInterval(int inddeb, int indto)
{
  _activeCovList.clear();
  for (int i = inddeb; i < indto; i++)
    _activeCovList.push_back(i);
}
