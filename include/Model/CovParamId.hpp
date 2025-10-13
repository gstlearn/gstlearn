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

#include "Enum/EConsElem.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"

namespace gstlrn
{

/**
 * This class is used:
 * - to define the constraints for the Model Automatic Fitting (authAssign true)
 * - to define the non-stationary parameters of a Model
 */
class GSTLEARN_EXPORT CovParamId : public AStringable, public ICloneable
{
public:
  CovParamId(Id igrf = 0,
             Id icov = 0,
             const EConsElem& elem = EConsElem::fromKey("UNKNOWN"),
             Id iv1 = 0,
             Id iv2 = 0);
  CovParamId(const CovParamId &m);
  CovParamId& operator=(const CovParamId &m);
  virtual ~CovParamId();

  /// ICloneable interface
  IMPLEMENT_CLONING(CovParamId)

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  static CovParamId* create(Id igrf = 0,
                            Id icov = 0,
                            const EConsElem &elem = EConsElem::fromKey("UNKNOWN"),
                            Id iv1 = 0,
                            Id iv2 = 0);

  Id init(Id igrf, Id icov, const EConsElem& type, Id iv1, Id iv2);

  const EConsElem& getType() const { return _elemType; }
  Id getIGrf()  const { return _igrf; }
  Id getICov()  const { return _icov; }
  Id getIV1()   const { return _iv1; }
  Id getIV2()   const { return _iv2; }

  void setType(const EConsElem& type) { _elemType = type; }

  bool matchIGrf(Id igrf0)              const { return (igrf0 < 0 || _igrf == igrf0); }
  bool matchICov(Id icov0)              const { return (icov0 < 0 || _icov == icov0); }
  bool matchType(const EConsElem& type0) const { return (type0 == EConsElem::fromKey("UNKNOWN") || _elemType == type0); }
  bool matchIV1(Id iv10)                const { return (iv10 < 0 || _iv1 == iv10); }
  bool matchIV2(Id iv20)                const { return (iv20 < 0 || _iv2 == iv20); }

private:
  Id       _igrf;       /* Rank of the Gaussian Random Function */
  Id       _icov;       /* Structure rank */
  EConsElem _elemType;   /* Type of element */
  Id       _iv1;        /* Rank of the first variable */
  Id       _iv2;        /* Rank of the second variable */
};
}
