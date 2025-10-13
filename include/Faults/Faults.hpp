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
#include "Basic/PolyLine2D.hpp"

namespace gstlrn
{
class SpacePoint;

class GSTLEARN_EXPORT Faults: public AStringable, public ASerializable
{
public:
  Faults();
  Faults(const Faults& r);
  Faults& operator=(const Faults& r);
  virtual ~Faults();

  /// Interface for AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  static Faults* createFromNF(const String& NFFilename, bool verbose = true);

  Id getNFaults() const { return static_cast<Id>(_faults.size()); }
  void addFault(const PolyLine2D& fault);

  inline const std::vector<PolyLine2D>& getFaults() const { return _faults; }
  inline const PolyLine2D& getFault(Id ifault) const { return _faults[ifault]; }

  bool isSplitByFault(double xt1, double yt1, double xt2, double yt2) const;
  bool isSplitByFaultSP(const SpacePoint& P1, const SpacePoint& P2) const;

protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "Faults"; }

private:
  std::vector<PolyLine2D> _faults;
};
} // namespace gstlrn