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
#include "Basic/AStringable.hpp"
#include "Basic/IClonable.hpp"
#include "Model/CovParamId.hpp"
#include "Model/EConsElem.hpp"
#include "Model/EConsType.hpp"

/**
 * This class is used:
 * - to define the constraints for the Model Automatic Fitting (authAssign true)
 * - to define the non-stationary parameters of a Model
 */
class GSTLEARN_EXPORT ConsItem : public AStringable, public IClonable
{
public:
  ConsItem(const CovParamId& paramid,
           const EConsType& type = EConsType::DEFAULT,
           double value = 0.);
  ConsItem(const ConsItem &m);
  ConsItem& operator= (const ConsItem &m);
  virtual ~ConsItem();

  int init(const CovParamId& paramid,
           const EConsType& type,
           double value = TEST);

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;
  virtual IClonable* clone() const override { return new ConsItem(*this); };

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
  static ConsItem define(const EConsElem& elem = EConsElem::UNKNOWN,
                         int icov = 0,
                         int iv1 = 0,
                         int iv2 = 0,
                         const EConsType& type = EConsType::DEFAULT,
                         double value = 0.);


private:
  CovParamId _paramId;
  EConsType  _type;       /* 0: Parameter; -1: Lower; 1: Upper; 2: Equal */
  double     _value;      /* Assigned value */
};
