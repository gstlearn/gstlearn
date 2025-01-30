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

CovCalcMode::CovCalcMode(const ECalcMember& member,
                         bool asVario,
                         bool unitary,
                         int orderVario)
  : AStringable()
  , _member(member)
  , _asVario(asVario)
  , _unitary(unitary)
  , _orderVario(orderVario)
{
}

CovCalcMode::CovCalcMode(const CovCalcMode& r)
  : AStringable(r)
  , _member(r._member)
  , _asVario(r._asVario)
  , _unitary(r._unitary)
  , _orderVario(r._orderVario)
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
  }
  return *this;
}

CovCalcMode::~CovCalcMode() {}

CovCalcMode* CovCalcMode::create(const ECalcMember &member,
                                 bool asVario,
                                 bool unitary,
                                 int orderVario)
{
  return new CovCalcMode(member, asVario, unitary, orderVario);
}
