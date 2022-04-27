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

class GSTLEARN_EXPORT TokenParameter: public AStringable
{
public:
  TokenParameter();
  TokenParameter(const TokenParameter &r);
  TokenParameter& operator=(const TokenParameter &r);
  virtual ~TokenParameter();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  ETLaw getLaw() const { return _law; }
  void setLaw(ETLaw law) { _law = law; }
  const VectorDouble& getValarg() const { return _valarg; }
  double getValarg(int i) const { return _valarg[i]; }
  void setValarg(int i, double value) { _valarg[i] = value; }

  double getValue() const;

private:
  ETLaw _law; /* Type of law */
  VectorDouble _valarg; /* Randomization arguments */
};
