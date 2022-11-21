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

#include "Enum/ETLaw.hpp"

#include "Basic/AStringable.hpp"

// TODO Will be replaced by future class"Law" or "Distribution" which does not
// actually exist
class GSTLEARN_EXPORT ShapeParameter: public AStringable
{
public:
  ShapeParameter(ETLaw law = ETLaw::fromKey("CONSTANT"), double value = 0.);
  ShapeParameter(const ShapeParameter &r);
  ShapeParameter& operator=(const ShapeParameter &r);
  virtual ~ShapeParameter();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  ETLaw getLaw() const { return _law; }
  const VectorDouble& getValarg() const { return _valarg; }
  double getValarg(int iarg) const;
  int getNbValarg() const { return (int) _valarg.size(); }

  void setLaw(ETLaw law) { _law = law; }
  void setValarg(int iarg, double value);

  double generateValue() const;

private:
  bool _isValidArgIndex(int iarg) const;

private:
  ETLaw _law; /* Type of law */
  VectorDouble _valarg; /* Randomization arguments */
};
