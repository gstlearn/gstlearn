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
#include "Basic/Vector.hpp"
#include "Boolean/TokenParameter.hpp"
#include "Boolean/ETShape.hpp"

class GSTLEARN_EXPORT AToken: public AStringable
{
public:
  AToken();
  AToken(const AToken &r);
  AToken& operator=(const AToken &r);
  virtual ~AToken();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface for AToken
  virtual ETShape getType() const = 0;
  virtual int getNArgs() const = 0;
  virtual int getFlagSymZ() const = 0;
  virtual String getParamName(int i) const = 0;
  virtual const TokenParameter& getParam(int i) const = 0;

  double getFactorX2Y() const { return _factorX2Y; }
  void setFactorX2Y(double factorX2Y) { _factorX2Y = factorX2Y; }
  double getFactorX2Z() const { return _factorX2Z; }
  void setFactorX2Z(double factorX2Z) { _factorX2Z = factorX2Z; }
  double getFactorY2Z() const { return _factorY2Z; }
  void setFactorY2Z(double factorY2Z) { _factorY2Z = factorY2Z; }
  double getProportion() const { return _proportion; }
  void setProportion(double proportion) { _proportion = proportion; }

private:
  double _factorX2Y; /* Link factor for the geometry from x to y */
  double _factorX2Z; /* Link factor for the geometry from x to z */
  double _factorY2Z; /* Link factor for the geometry from y to z */
  double _proportion; /* Token Proportion */
};
