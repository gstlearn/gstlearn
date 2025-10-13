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

#include "Basic/ICloneable.hpp"
#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Enum/ENeigh.hpp"

#include "Neigh/ANeigh.hpp"
#include "Space/ASpace.hpp"

namespace gstlrn
{
class Db;

/**
 * \brief
 * Unique Neighborhood definition.
 *
 * The Neighborhood is usually meant to select a sub-population from the input Data Base,
 * containing the active samples close to the target.
 *
 * The Unique Neighborhood selects all the active samples. Nevertheless, it offers
 * the possibility to suppress any sample which would be too close to (or coincide with)
 * the target: this is the cross-validation option.
 */
class GSTLEARN_EXPORT NeighUnique: public ANeigh
{
public:
  NeighUnique(bool flag_xvalid = false, const ASpaceSharedPtr& space = ASpaceSharedPtr());
  NeighUnique(const NeighUnique& r);
  NeighUnique& operator=(const NeighUnique& r);
  virtual ~NeighUnique();

  IMPLEMENT_CLONING(NeighUnique)
  /// Interface for ANeigh
  void getNeigh(Id iech_out, VectorInt& ranks) override;
  Id getNSampleMax(const Db* db) const override;
  bool hasChanged(Id iech_out) const override;
  ENeigh getType() const override { return ENeigh::fromKey("UNIQUE"); }

  /// Interface for AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  static NeighUnique* create(bool flag_xvalid = false, const ASpaceSharedPtr& space = ASpaceSharedPtr());
  static NeighUnique* createFromNF(const String& NFFilename, bool verbose = true);

protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "NeighUnique"; }

private:
  void _unique(Id iech_out, VectorInt& ranks);
};
}