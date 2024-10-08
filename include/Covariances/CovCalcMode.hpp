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
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/ECalcMember.hpp"

#include "Basic/AStringable.hpp"

class GSTLEARN_EXPORT CovCalcMode : public AStringable
{
public:
  CovCalcMode(const ECalcMember &member = ECalcMember::fromKey("LHS"),
              bool asVario = false,
              bool unitary = false,
              int orderVario = 0,
              bool allActiveCov = true,
              const VectorInt &activeCovList = VectorInt());
  CovCalcMode(const CovCalcMode &r);
  CovCalcMode& operator= (const CovCalcMode &r);
  virtual ~CovCalcMode();

  static CovCalcMode* create(const ECalcMember &member = ECalcMember::fromKey("LHS"), 
                             bool asVario = false,
                             bool unitary = false,
                             int orderVario = 0,
                             bool allActiveCov = true,
                             const VectorInt& activeCovList = VectorInt());

  const ECalcMember&    getMember()             const { return _member; }
  bool                  getAsVario()            const { return _asVario; }
  bool                  getUnitary()            const { return _unitary; }
  int                   getOrderVario()         const { return _orderVario; }
  const VectorInt&      getActiveCovList()      const { return _activeCovList; }
  int                   getActiveCovList(int i) const { return _activeCovList[i]; }
  bool                  isAllActiveCov()        const { return _allActiveCov; }

  void setAsVario(bool asVario) { _asVario = asVario; }
  void setMember(const ECalcMember& member) { _member = member; }
  void setUnitary(bool unitary) { _unitary = unitary; }
  void setOrderVario(int orderVario) { _orderVario = orderVario; }

  void setActiveCovListFromOne(int keepOnlyCovIdx);
  void setActiveCovListFromInterval(int inddeb, int indto);
  void setActiveCovList(const VectorInt &activeCovList, bool allActiveCov);

private:
  ECalcMember   _member;         /*! LHS (default), RHS or VAR(IANCE) */
  bool          _asVario;        /*! True to calculate variogram instead of covariance */
  bool          _unitary;        /*! True to calculate covariance without sill (in Goulard) */
  int           _orderVario;     /*! Higher Variogram Order (0: standard) */
  bool          _allActiveCov;   /*! True if all covariances are active */
  VectorInt     _activeCovList;  /*! List of indices of the active covariances */
};
