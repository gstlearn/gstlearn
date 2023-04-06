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

#include "Enum/EConsElem.hpp"
#include "Enum/EConsType.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"
#include "Model/CovParamId.hpp"

/**
 * This class is used:
 * - to define the constraints for the Model Automatic Fitting (authAssign true)
 * - to define the non-stationary parameters of a Model
 */
class GSTLEARN_EXPORT ConsItem : public AStringable, public ICloneable
{
public:
  ConsItem(const CovParamId& paramid,
           const EConsType& type = EConsType::fromKey("DEFAULT"),
           double value = 0.);
  ConsItem(const ConsItem &m);
  ConsItem& operator= (const ConsItem &m);
  virtual ~ConsItem();

  /// ICloneable interface
  IMPLEMENT_CLONING(ConsItem)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static ConsItem* create(const CovParamId &paramid,
                          const EConsType &type = EConsType::fromKey("DEFAULT"),
                          double value = 0.);
  static ConsItem* createFromParamId(int icov = 0,
                                     const EConsElem& elem = EConsElem::fromKey("UNKNOWN"),
                                     const EConsType &type = EConsType::fromKey("DEFAULT"),
                                     double value = 0.,
                                     int igrf = 0,
                                     int iv1 = 0,
                                     int iv2 = 0);

  // Pipe to the CovParamId class
  const EConsElem& getType() const { return _paramId.getType(); }
  int getIGrf()  const { return _paramId.getIGrf(); }
  int getICov()  const { return _paramId.getICov(); }
  int getIV1()   const { return _paramId.getIV1(); }
  int getIV2()   const { return _paramId.getIV2(); }

  void setValue(double value)          { _value = value; }
  void setIcase(const EConsType& type) { _type = type;   }
  const EConsType& getIcase() const { return _type; }
  double getValue() const { return _value; }

  bool matchIGrf(int igrf0)              const { return _paramId.matchIGrf(igrf0); }
  bool matchICov(int icov0)              const { return _paramId.matchICov(icov0); }
  bool matchType(const EConsElem& type0) const { return _paramId.matchType(type0); }
  bool matchIV1(int iv10)                const { return _paramId.matchIV1(iv10); }
  bool matchIV2(int iv20)                const { return _paramId.matchIV2(iv20); }

  const CovParamId& getParamId() const { return _paramId; }

  /**
   * This function creates a constraints on a parameter
   * This constraint will be used subsequently during the variogram fitting
   * @param elem The type of item on which the constraints applies (EConsElem.hpp)
   * @param icov The rank of the covariance
   * @param iv1  The rank of the first variable
   * @param iv2  The rank of the second variable
   * @param type The type of constraints (EConsType.hpp)
   * @param value The value assigned to the constraint
   *
   * @remark Do not forget to delete object after usage
   * @return
   */
  static ConsItem define(const EConsElem& elem = EConsElem::fromKey("UNKNOWN"),
                         int icov = 0,
                         int iv1 = 0,
                         int iv2 = 0,
                         const EConsType& type = EConsType::fromKey("DEFAULT"),
                         double value = 0.);

private:
  int _init(const CovParamId &paramid,
            const EConsType &type,
            double value = TEST);

private:
  CovParamId _paramId;
  EConsType  _type;       /* 0: Parameter; -1: Lower; 1: Upper; 2: Equal */
  double     _value;      /* Assigned value */
};
