/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/PolyLine2D.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"

class GSTLEARN_EXPORT Faults: public AStringable, public ASerializable
{
public:
  Faults();
  Faults(const Faults& r);
  Faults& operator=(const Faults& r);
  virtual ~Faults();

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static Faults* createFromNF(const String& neutralFilename, bool verbose = true);

  int getNFaults() const { return (int) _faults.size(); }
  void addFault(const PolyLine2D& fault);

  std::vector<PolyLine2D> getFaults() const { return _faults; }
  const PolyLine2D getFault(int ifault) const { return _faults[ifault]; }

  bool isSplitByFault(double xt1,double yt1, double xt2, double yt2) const;

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "Faults"; }

private:
  std::vector<PolyLine2D> _faults;
};
