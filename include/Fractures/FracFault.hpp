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

#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"

namespace gstlrn
{
class GSTLEARN_EXPORT FracFault: public AStringable, public ASerializable
{
public:
  FracFault(double coord = 0., double orient = 0.);
  FracFault(const FracFault& r);
  FracFault& operator=(const FracFault& r);
  virtual ~FracFault();

  /// Interface for AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  double getCoord() const { return _coord; }
  double getOrient() const { return _orient; }
  VectorDouble getRangel() const { return _rangel; }
  VectorDouble getRanger() const { return _ranger; }
  VectorDouble getThetal() const { return _thetal; }
  VectorDouble getThetar() const { return _thetar; }
  double getRangel(Id i) const { return _rangel[i]; }
  double getRanger(Id i) const { return _ranger[i]; }
  double getThetal(Id i) const { return _thetal[i]; }
  double getThetar(Id i) const { return _thetar[i]; }

  Id getNFamilies() const { return static_cast<Id>(_thetal.size()); }
  double faultAbscissae(double cote) const;

  void addFaultPerFamily(double thetal,
                         double thetar,
                         double rangel,
                         double ranger);

  void setRangel(const VectorDouble& rangel) { _rangel = rangel; }
  void setRanger(const VectorDouble& ranger) { _ranger = ranger; }
  void setThetal(const VectorDouble& thetal) { _thetal = thetal; }
  void setThetar(const VectorDouble& thetar) { _thetar = thetar; }

protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "FracFault"; }

private:
  double _coord;        //!< Abscissas of the first Fault point
  double _orient;       //!< Fault orientation
  VectorDouble _thetal; //!< Maximum density on left
  VectorDouble _thetar; //!< Maximum density on right
  VectorDouble _rangel; //!< Decrease range on left
  VectorDouble _ranger; //!< Decrease range on right

  friend class FracEnviron;
};
} // namespace gstlrn