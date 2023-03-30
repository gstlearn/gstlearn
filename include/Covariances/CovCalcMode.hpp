/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
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
              bool normalized = false,
              bool filterNugget = false,
              unsigned int keepOnlyCovIdx = -1,
              bool unitary = false,
              int envelop = 0,
              int orderVario = 0);
  CovCalcMode(const CovCalcMode &r);
  CovCalcMode& operator= (const CovCalcMode &r);
  virtual ~CovCalcMode();

  bool isEqual(const CovCalcMode &r) const;

  void update(const ECalcMember& member     = ECalcMember::fromKey("LHS"),
              int                nugget_opt = 0,
              int                nostd      = 0,
              int                icov_r     = -1,
              int                flag_norm  = 0,
              int                flag_cov   = 1);

  const ECalcMember&    getMember()           const { return _member; }
  bool                  getAsVario()          const { return _asVario; }
  bool                  getNormalized()       const { return _normalized; }
  bool                  isFilterNugget()      const { return _filterNugget; }
  int                   getKeepOnlyCovIdx()   const { return _keepOnlyCovIdx; }
  bool                  getUnitary()          const { return _unitary; }
  int                   getOrderVario()       const { return _orderVario; }
  int                   getEnvelop()          const { return _envelop; }
  const VectorBool&     getCovFiltered()      const { return _covFiltered; }
  bool                  getCovFiltered(int i) const;
  int                   getIndexClass()       const { return _indexClass; }

  void setAsVario(bool asVario) { _asVario = asVario; }
  void setMember(const ECalcMember& member) { _member = member; }
  void setFilterNugget(bool filterNugget) { _filterNugget = filterNugget; }
  void setKeepOnlyCovIdx(int keepOnlyCovIdx) { _keepOnlyCovIdx = keepOnlyCovIdx; }
  void setUnitary(bool unitary) { _unitary = unitary; }
  void setNormalized(bool normalized) { _normalized = normalized; }
  void setEnvelop(int envelop) { _envelop = envelop; }
  void setOrderVario(int orderVario) { _orderVario = orderVario; }
  void setCovFiltered(const VectorBool& covFiltered) { _covFiltered = covFiltered; }
  void setCovFiltered(int i, bool status);
  void setAllCovFiltered(int ncov, bool status);
  void setIndexClass(int indexClass) { _indexClass = indexClass; }

private:
  ECalcMember   _member;         /*! LHS (default), RHS or VAR(IANCE) */
  bool          _asVario;        /*! True to calculate variogram and not covariance (default = false)*/
  bool          _normalized;     /*! Normalized variogram */
  bool          _filterNugget;   /*! True to filter nugget structure (default = false) */
  int           _keepOnlyCovIdx; /*! Index of the covariance to be kept (default is -1) */
  bool          _unitary;        /*! True to calculate covariance without sill (in Goulard) */
  int           _envelop;        /*! Envelop of Multivariate model: 1(upper) or -1(lower) */
  int           _orderVario;     /*! Higher Variogram Order (0: standard) */
  int           _indexClass;     /*! Index of the working class (used when calculating Covariance) */
  VectorBool    _covFiltered;    /*! Array of Covariance filtering flags (optional) */
};
