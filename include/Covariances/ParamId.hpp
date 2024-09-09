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

#include "Enum/EConsElem.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"

/**
 * This class is used:
 * - to define the constraints for the Model Automatic Fitting (authAssign true)
 * - to define the non-stationary parameters of a Model
 */
class GSTLEARN_EXPORT ParamId : public AStringable, public ICloneable
{
public:
  ParamId(const EConsElem& elem = EConsElem::fromKey("UNKNOWN"),
          int iv1 = 0,
          int iv2 = 0);
  ParamId(const ParamId &m);
  ParamId& operator=(const ParamId &m);
  virtual ~ParamId();

  /// ICloneable interface
  IMPLEMENT_CLONING(ParamId)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static ParamId* create(   const EConsElem &elem = EConsElem::fromKey("UNKNOWN"),
                            int iv1 = 0,
                            int iv2 = 0);

  int init(const EConsElem& type, int iv1, int iv2);

  const EConsElem& getType() const { return _elemType; }
  int getIV1()   const { return _iv1; }
  int getIV2()   const { return _iv2; }

  void setType(const EConsElem& type) { _elemType = type; }

  bool matchType(const EConsElem& type0) const { return (type0 == EConsElem::fromKey("UNKNOWN") || _elemType == type0); }
  bool matchIV1(int iv10)                const { return (iv10 < 0 || _iv1 == iv10); }
  bool matchIV2(int iv20)                const { return (iv20 < 0 || _iv2 == iv20); }

private:
  EConsElem _elemType;   /* Type of element */
  int       _iv1;        /* Rank of the first variable */
  int       _iv2;        /* Rank of the second variable */
};
