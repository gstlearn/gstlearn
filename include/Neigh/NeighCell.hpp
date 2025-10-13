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

#include <Geometry/BiTargetCheckCell.hpp>
#include "Basic/ICloneable.hpp"
#include "Space/ASpace.hpp"
#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Enum/ENeigh.hpp"

#include "Neigh/ANeigh.hpp"
#include "Space/SpaceTarget.hpp"

namespace gstlrn
{
class Db;

/**
 * \brief
 * Neighborhood definition by Cell.
 *
 * The Neighborhood is usually meant to select a sub-population from the input Data Base,
 * containing the active samples close to the target.
 *
 * The selected samples belong to the same 'cell' as the target. This obviously requires
 * that the target belongs to a Grid.
 *
 * The neighborhood also offers the possibility to suppress any sample which would be too close to (or coincide with)
 * the target: this is the cross-validation option.
 */
class GSTLEARN_EXPORT NeighCell: public ANeigh
{
public:
  NeighCell(bool flag_xvalid             = false,
            Id nmini                    = 1,
            bool useBallTree             = false,
            Id leaf_size                = 10,
            const ASpaceSharedPtr& space = ASpaceSharedPtr());
  NeighCell(const NeighCell& r);
  NeighCell& operator=(const NeighCell& r);
  virtual ~NeighCell();

  IMPLEMENT_CLONING(NeighCell)
  /// Interface for ANeigh
  Id attach(const Db *dbin, const Db *dbout = nullptr) override;
  void getNeigh(Id iech_out, VectorInt& ranks) override;
  bool hasChanged(Id iech_out) const override;
  Id getNSampleMax(const Db* db) const override {
    DECLARE_UNUSED(db);
    return 0;
  }
  ENeigh getType() const override { return ENeigh::fromKey("CELL"); }

  /// Interface for AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  static NeighCell* create(bool flag_xvalid             = false,
                           Id nmini                    = 1,
                           bool useBallTree             = false,
                           Id leaf_size                = 10,
                           const ASpaceSharedPtr& space = ASpaceSharedPtr());
  static NeighCell* createFromNF(const String& NFFilename, bool verbose = true);

  Id getNMini() const { return _nMini; }

private:
  Id _cell(Id iech_out, VectorInt& ranks);

protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "NeighCell"; }

private:
  Id _nMini;
  BiTargetCheckCell* _biPtCell;

  mutable SpaceTarget _T1;
  mutable SpaceTarget _T2;
};
}