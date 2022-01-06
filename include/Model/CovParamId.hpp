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
#include "Model/EConsElem.hpp"

/**
 * This class is used:
 * - to define the constraints for the Model Automatic Fitting (authAssign true)
 * - to define the non-stationary parameters of a Model
 */
class GSTLEARN_EXPORT CovParamId : public AStringable, public IClonable
{
public:
  CovParamId(int igrf = 0,
             int icov = 0,
             const EConsElem& elem = EConsElem::UNKNOWN,
             int iv1 = 0,
             int iv2 = 0);
  CovParamId(const CovParamId &m);
  CovParamId& operator=(const CovParamId &m);
  virtual ~CovParamId();

  int init(int igrf, int icov, const EConsElem& type, int iv1, int iv2);

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;
  virtual IClonable* clone() const override { return new CovParamId(*this); };

  const EConsElem& getType() const { return _elemType; }
  int getIGrf()  const { return _igrf; }
  int getICov()  const { return _icov; }
  int getIV1()   const { return _iv1; }
  int getIV2()   const { return _iv2; }

  void setType(const EConsElem& type) { _elemType = type; }

  bool matchIGrf(int igrf0)              const { return (igrf0 < 0 || _igrf == igrf0); }
  bool matchICov(int icov0)              const { return (icov0 < 0 || _icov == icov0); }
  bool matchType(const EConsElem& type0) const { return (type0 == EConsElem::UNKNOWN || _elemType == type0); }
  bool matchIV1(int iv10)                const { return (iv10 < 0 || _iv1 == iv10); }
  bool matchIV2(int iv20)                const { return (iv20 < 0 || _iv2 == iv20); }

private:
  int       _igrf;       /* Rank of the Gaussian Random Function */
  int       _icov;       /* Structure rank */
  EConsElem _elemType;   /* Type of element */
  int       _iv1;        /* Rank of the first variable */
  int       _iv2;        /* Rank of the second variable */
};
