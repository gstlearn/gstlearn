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

class GSTLEARN_EXPORT Fault: public AStringable
{
public:
  Fault(double coord = 0., double orient = 0.);
  Fault(const Fault& r);
  Fault& operator=(const Fault& r);
  virtual ~Fault();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  double getCoord() const { return _coord; }
  double getOrient() const { return _orient; }
  VectorDouble getRangel() const { return _rangel; }
  VectorDouble getRanger() const { return _ranger; }
  VectorDouble getThetal() const { return _thetal; }
  VectorDouble getThetar() const { return _thetar; }
  double getRangel(int i) const { return _rangel[i]; }
  double getRanger(int i) const { return _ranger[i]; }
  double getThetal(int i) const { return _thetal[i]; }
  double getThetar(int i) const { return _thetar[i]; }

  int getNFamilies() const { return  (int) _thetal.size(); }
  double faultAbscissae(double cote) const;

  void addFaultPerFamily(double thetal,
                         double thetar,
                         double rangel,
                         double ranger);

  void setRangel(const VectorDouble& rangel) { _rangel = rangel; }
  void setRanger(const VectorDouble& ranger) { _ranger = ranger; }
  void setThetal(const VectorDouble& thetal) { _thetal = thetal; }
  void setThetar(const VectorDouble& thetar) { _thetar = thetar; }

private:
  double _coord;                //!< Abscissas of the first Fault point
  double _orient;               //!< Fault orientation
  VectorDouble _thetal;         //!< Maximum density on left
  VectorDouble _thetar;         //!< Maximum density on right
  VectorDouble _rangel;         //!< Decrease range on left
  VectorDouble _ranger;         //!< Decrease range on right
};
