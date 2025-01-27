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
#include "Covariances/CovCalcMode.hpp"

CovCalcMode::CovCalcMode(const ECalcMember &member,
                         bool asVario,
                         bool unitary,
                         int orderVario,
                         bool allActiveCov,
                         const VectorInt &activeCovList)
    : AStringable(),
      _member(member),
      _asVario(asVario),
      _unitary(unitary),
      _orderVario(orderVario),
      _allActiveCov(allActiveCov),
      _activeCovList(activeCovList)
{
}

CovCalcMode::CovCalcMode(const CovCalcMode &r)
    : AStringable(r),
      _member(r._member),
      _asVario(r._asVario),
      _unitary(r._unitary),
      _orderVario(r._orderVario),
      _allActiveCov(r._allActiveCov),
      _activeCovList(r._activeCovList)
{
}

CovCalcMode& CovCalcMode::operator=(const CovCalcMode &r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _member        = r._member;
    _asVario       = r._asVario;
    _unitary       = r._unitary;
    _orderVario    = r._orderVario;
    _allActiveCov  = r._allActiveCov;
    _activeCovList = r._activeCovList;
  }
  return *this;
}

CovCalcMode::~CovCalcMode() {}

CovCalcMode* CovCalcMode::create(const ECalcMember &member,
                                 bool asVario,
                                 bool unitary,
                                 int orderVario,
                                 bool allActiveCov,
                                 const VectorInt& activeCovList)
{
  return new CovCalcMode(member, asVario, unitary, orderVario, allActiveCov, activeCovList);
}

void CovCalcMode::setActiveCovListFromOne(int keepOnlyCovIdx)
{
  _allActiveCov = true;
  _activeCovList.clear();
  if (keepOnlyCovIdx >= 0)
  {
    _activeCovList.push_back(keepOnlyCovIdx);
    _allActiveCov = false;
  }
}

/**
 * Set the list of active covariances from an interval
 * @param inddeb Lower bound of the interval (included)
 * @param indto  Upper bound of the interval (excluded)
 */
void CovCalcMode::setActiveCovListFromInterval(int inddeb, int indto)
{
  _activeCovList.clear();
  for (int i = inddeb; i < indto; i++) _activeCovList.push_back(i);
  _allActiveCov = false;
}

void CovCalcMode::setActiveCovList(const VectorInt &activeCovList, bool allActiveCov)
{
  _activeCovList = activeCovList;
  _allActiveCov = allActiveCov;
}
