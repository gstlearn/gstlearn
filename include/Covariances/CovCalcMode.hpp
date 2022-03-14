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
#pragma once

#include "gstlearn_export.hpp"

// WARNING: Make this include list as small as possible!
#include "Basic/AStringable.hpp"
#include "Covariances/ECalcMember.hpp"

class GSTLEARN_EXPORT CovCalcMode : public AStringable
{
public:
  CovCalcMode(const ECalcMember& member = ECalcMember::LHS,
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

  const ECalcMember& getMember()         const { return _member; }
  bool                  getAsVario()        const { return _asVario; }
  bool                  getNormalized()     const { return _normalized; }
  bool                  isFilterNugget()    const { return _filterNugget; }
  unsigned int          getKeepOnlyCovIdx() const { return _keepOnlyCovIdx; }
  bool                  getUnitary()        const { return _unitary; }
  int                   getOrderVario()     const { return _orderVario; }
  int                   getEnvelop()        const { return _envelop; }

  void setAsVario(bool asVario) { _asVario = asVario; }
  void setMember(const ECalcMember& member) { _member = member; }
  void setFilterNugget(bool filterNugget) { _filterNugget = filterNugget; }
  void setKeepOnlyCovIdx(unsigned int keepOnlyCovIdx) { _keepOnlyCovIdx = keepOnlyCovIdx; }
  void setUnitary(bool unitary) { _unitary = unitary; }
  void setNormalized(bool normalized) { _normalized = normalized; }
  void setEnvelop(int envelop) { _envelop = envelop; }
  void setOrderVario(int orderVario) { _orderVario = orderVario; }

  void update(int                nugget_opt = 0,
              int                nostd      = 0,
              const ECalcMember& member     = ECalcMember::LHS,
              int                icov_r     = -1,
              int                flag_norm  = 0,
              int                flag_cov   = 1);

private:
  ECalcMember   _member;         /*! LHS (default), RHS or VAR(IANCE) */
  bool          _asVario;        /*! True to calculate variogram and not covariance (default = false)*/
  bool          _normalized;     /*! Normalized variogram */
  bool          _filterNugget;   /*! True to filter nugget structure (default = false) */
  unsigned int  _keepOnlyCovIdx; /*! Index of the covariance to be kept (default is -1) */
  bool          _unitary;        /*! True to calculate covariance without sill (in Goulard) */
  int           _envelop;        /*! Envelop of Multivariate model: 1(upper) or -1(lower) */
  int           _orderVario;     /*! Higher Variogram Order (0: standard) */
};
