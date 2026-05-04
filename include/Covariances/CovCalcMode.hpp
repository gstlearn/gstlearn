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

namespace gstlrn
{
  /**
   * @brief This class contains the information about the mode of covariance calculation, such as:
   * - whether to calculate variogram instead of covariance,
   * - whether to calculate unitary covariance (without sill, needed in Goulard's method),
   * - whether to calculate sill for coregionalization envelop,
   * - the higher variogram order (0 for standard variogram).
   */
  class GSTLEARN_EXPORT CovCalcMode: public AStringable
  {
  public:
    CovCalcMode(
      const ECalcMember& member = ECalcMember::fromKey("LHS"),
      bool asVario = false,
      bool unitary = false,
      Id orderVario = 0);
    CovCalcMode(const CovCalcMode& r);
    CovCalcMode& operator=(const CovCalcMode& r);
    virtual ~CovCalcMode();

    static CovCalcMode* create(
      const ECalcMember& member = ECalcMember::fromKey("LHS"),
      bool asVario = false,
      bool unitary = false,
      Id orderVario = 0);

    const ECalcMember& getMember() const { return _member; }

    bool getAsVario() const { return _asVario; }

    bool getUnitary() const { return _unitary; }

    bool getEnvelop() const { return _envelop; }

    Id getOrderVario() const { return _orderVario; }

    void setAsVario(bool asVario) { _asVario = asVario; }

    void setMember(const ECalcMember& member) { _member = member; }

    void setUnitary(bool unitary) { _unitary = unitary; }

    void setEnvelop(bool envelop) { _envelop = envelop; }

    void setOrderVario(Id orderVario) { _orderVario = orderVario; }

  private:
    ECalcMember _member; /*! LHS (default), RHS or VAR(IANCE) */
    bool _asVario; /*! True to calculate variogram instead of covariance */
    bool _unitary; /*! True to calculate covariance without sill (in Goulard) */
    bool _envelop; /*! True to calculate sill for coregionalization envelop */
    Id _orderVario; /*! Higher Variogram Order (0: standard) */
  };
} // namespace gstlrn
