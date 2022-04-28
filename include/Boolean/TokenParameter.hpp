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

#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"
#include "Boolean/ETLaw.hpp"

// TODO Will be replaced by future class"Law" or "Distribution" which does not
// actually exist
class GSTLEARN_EXPORT TokenParameter: public AStringable
{
public:
  TokenParameter(ETLaw law = ETLaw::CONSTANT, double value = 0.);
  TokenParameter(const TokenParameter &r);
  TokenParameter& operator=(const TokenParameter &r);
  virtual ~TokenParameter();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  ETLaw getLaw() const { return _law; }
  const VectorDouble& getValarg() const { return _valarg; }
  double getValarg(int iarg) const;

  void setLaw(ETLaw law) { _law = law; }
  void setValarg(int iarg, double value);

  double generateValue() const;

private:
  bool _isValidArgIndex(int iarg) const;

private:
  ETLaw _law; /* Type of law */
  VectorDouble _valarg; /* Randomization arguments */
};
