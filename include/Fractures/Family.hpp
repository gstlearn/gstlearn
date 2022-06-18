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

class GSTLEARN_EXPORT Family: public AStringable
{
public:
  Family(double orient = 0.,
         double dorient = 0.,
         double theta0 = 0.,
         double alpha = 0.,
         double ratcst = 0.,
         double prop1 = 0.,
         double prop2 = 0.,
         double aterm = 0.,
         double bterm = 0.,
         double range = 0.);
  Family(const Family& r);
  Family& operator=(const Family& r);
  virtual ~Family();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  double getAlpha() const { return _alpha; }
  void setAlpha(double alpha) { _alpha = alpha; }
  double getAterm() const { return _aterm; }
  void setAterm(double aterm) { _aterm = aterm; }
  double getBterm() const { return _bterm; }
  void setBterm(double bterm) { _bterm = bterm; }
  double getDorient() const { return _dorient; }
  void setDorient(double dorient) { _dorient = dorient; }
  double getOrient() const { return _orient; }
  void setOrient(double orient) { _orient = orient; }
  double getProp1() const { return _prop1; }
  void setProp1(double prop1) { _prop1 = prop1; }
  double getProp2() const { return _prop2; }
  void setProp2(double prop2) { _prop2 = prop2; }
  double getRange() const { return _range; }
  void setRange(double range) { _range = range; }
  double getRatcst() const { return _ratcst; }
  void setRatcst(double ratcst) { _ratcst = ratcst; }
  double getTheta0() const { return _theta0; }
  void setTheta0(double theta0) { _theta0 = theta0; }

private:
  double _orient;              //!< Mean orientation
  double _dorient;             //!< Standard deviation for orientation
  double _theta0;              //!< Reference Poisson intensity
  double _alpha;               //!< Power dependency between layer & intensity
  double _ratcst;              //!< Ratio of Constant vs. shaped intensity
  double _prop1;               //!< Survival probability (constant term)
  double _prop2;               //!< Survival probability (length dependent term)
  double _aterm;               //!< Survival probability (cumulative length term)
  double _bterm;               //!< Survival probability (layer thickness term)
  double _range;               //!< Range of fracture repulsion area
};
