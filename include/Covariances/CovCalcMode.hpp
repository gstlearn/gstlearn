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
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/ECalcMember.hpp"

#include "Basic/AStringable.hpp"

class GSTLEARN_EXPORT CovCalcMode : public AStringable
{
public:
  CovCalcMode(const ECalcMember& member = ECalcMember::fromKey("LHS"),
              bool asVario = false,
              bool filterNugget = false,
              int keepOnlyCovIdx = -1,
              bool unitary = false,
              int orderVario = 0);
  CovCalcMode(const CovCalcMode &r);
  CovCalcMode& operator= (const CovCalcMode &r);
  virtual ~CovCalcMode();

  static CovCalcMode* create(const ECalcMember &member = ECalcMember::fromKey("LHS"),
                             bool asVario = false,
                             bool filterNugget = false,
                             int keepOnlyCovIdx = -1,
                             bool unitary = false,
                             int orderVario = 0);

  bool                  isFactorySettings()   const { return _factorySettings; }
  const ECalcMember&    getMember()           const { return _member; }
  bool                  getAsVario()          const { return _asVario; }
  bool                  isFilterNugget()      const { return _filterNugget; }
  int                   getKeepOnlyCovIdx()   const { return _keepOnlyCovIdx; }
  bool                  getUnitary()          const { return _unitary; }
  int                   getOrderVario()       const { return _orderVario; }
  const VectorBool&     getCovFiltered()      const { return _covFiltered; }
  bool                  getCovFiltered(int i) const;

  void setAsVario(bool asVario) { _asVario = asVario; _checkFactorySettings(); }
  void setMember(const ECalcMember& member) { _member = member; _checkFactorySettings(); }
  void setFilterNugget(bool filterNugget) { _filterNugget = filterNugget; _checkFactorySettings(); }
  void setKeepOnlyCovIdx(int keepOnlyCovIdx) { _keepOnlyCovIdx = keepOnlyCovIdx; _checkFactorySettings(); }
  void setUnitary(bool unitary) { _unitary = unitary; _checkFactorySettings(); }
  void setOrderVario(int orderVario) { _orderVario = orderVario; _checkFactorySettings(); }
  void setCovFiltered(const VectorBool& covFiltered) { _covFiltered = covFiltered; _checkFactorySettings(); }
  void setCovFiltered(int i, bool status);
  void setAllCovFiltered(int ncov, bool status);

private:
  void _checkFactorySettings(const ECalcMember& member = ECalcMember::fromKey("LHS"),
                             bool asVario = false,
                             bool filterNugget = false,
                             int keepOnlyCovIdx = -1,
                             bool unitary = false,
                             int orderVario = 0);

private:
  bool          _factorySettings;
  ECalcMember   _member;         /*! LHS (default), RHS or VAR(IANCE) */
  bool          _asVario;        /*! True to calculate variogram and not covariance (default = false) */
  bool          _filterNugget;   /*! True to filter nugget structure (default = false) */
  int           _keepOnlyCovIdx; /*! Index of the covariance to be kept (default is -1) */
  bool          _unitary;        /*! True to calculate covariance without sill (in Goulard) */
  int           _orderVario;     /*! Higher Variogram Order (0: standard) */
  VectorBool    _covFiltered;    /*! Array of Covariance filtering flags (optional) */
};
