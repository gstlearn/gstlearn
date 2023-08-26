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

#include "Enum/EShape.hpp"
#include "Enum/ELaw.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"
#include "Boolean/ShapeParameter.hpp"

class BooleanObject;

/**
 * Class defining the generic shape of the objects for Boolean Model
 */
class GSTLEARN_EXPORT AShape: public AStringable, public ICloneable
{
public:
  AShape();
  AShape(const AShape &r);
  AShape& operator=(const AShape &r);
  virtual ~AShape();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface for AShape
  virtual EShape getType() const = 0;
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
  void setLaw(int ipar, ELaw law);

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
