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
#include "Basic/Vector.hpp"
#include "Boolean/ShapeParameter.hpp"
#include "Boolean/ETShape.hpp"
#include "Boolean/ETLaw.hpp"

class BooleanObject;

/**
 * Class defining the generic shape of the objects for Boolean Model
 */
class GSTLEARN_EXPORT AShape: public AStringable, public IClonable
{
public:
  AShape();
  AShape(const AShape &r);
  AShape& operator=(const AShape &r);
  virtual ~AShape();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface for IClonable
  virtual IClonable* clone() const override = 0;

  /// Interface for AShape
  virtual ETShape getType() const = 0;
  virtual int  getNParams() const = 0;
  virtual bool getFlagCutZ() const = 0;
  virtual BooleanObject* generateObject(int ndim = 3) = 0;
  virtual bool belongObject(const VectorDouble& coor, const BooleanObject* object) const = 0;

  double getFactorX2Y() const { return _factorX2Y; }
  double getFactorX2Z() const { return _factorX2Z; }
  double getFactorY2Z() const { return _factorY2Z; }
  double getProportion() const { return _proportion; }
  String getParamName(int ipar) const;
  double getParam(int ipar, int iarg) const;
  const ShapeParameter& getParam(int ipar) const;

  void setFactorX2Y(double factorX2Y) { _factorX2Y = factorX2Y; }
  void setFactorX2Z(double factorX2Z) { _factorX2Z = factorX2Z; }
  void setFactorY2Z(double factorY2Z) { _factorY2Z = factorY2Z; }
  void setProportion(double proportion) { _proportion = proportion; }
  void setParamName(int ipar, const String& name);
  void setParam(int ipar, int iarg, double value);
  void setParamDefault(int ipar, const String& name, double value);
  void setLaw(int ipar, ETLaw law);

  void initParams(int count);

  double generateParam(int ipar) const;

private:
  bool _isValidParamIndex(int ipar) const;

private:
  double _factorX2Y; /* Link factor for the geometry from x to y */
  double _factorX2Z; /* Link factor for the geometry from x to z */
  double _factorY2Z; /* Link factor for the geometry from y to z */
  double _proportion; /* Proportion for object of current shape */
  VectorString _paramNames;
  std::vector<ShapeParameter> _params; // TODO map (regrouping the two last lines)
};
