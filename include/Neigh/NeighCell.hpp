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
#include "Space/ASpace.hpp"
#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Enum/ENeigh.hpp"

#include "Neigh/ANeigh.hpp"
#include "Space/SpaceTarget.hpp"

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
  NeighCell(bool flag_xvalid = false, int nmini = 1, const ASpaceSharedPtr& space = ASpaceSharedPtr());
  NeighCell(const NeighCell& r);
  NeighCell& operator=(const NeighCell& r);
  virtual ~NeighCell();

  /// Interface for ANeigh
  virtual int attach(const Db *dbin, const Db *dbout = nullptr) override;
  virtual void getNeigh(int iech_out, VectorInt& ranks) override;
  virtual bool hasChanged(int iech_out) const override;
  virtual int getNSampleMax(const Db* db) const override {
    DECLARE_UNUSED(db);
    return 0;
  }
  virtual ENeigh getType() const override { return ENeigh::fromKey("CELL"); }

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static NeighCell* create(bool flag_xvalid = false,
                           int nmini = 1,
                           const ASpaceSharedPtr& space = ASpaceSharedPtr());
  static NeighCell* createFromNF(const String& neutralFilename, bool verbose = true);

  int getNMini() const { return _nMini; }

private:
  int _cell(int iech_out, VectorInt& ranks);

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "NeighCell"; }

private:
  int _nMini;
  BiTargetCheckCell* _biPtCell;

  mutable SpaceTarget _T1;
  mutable SpaceTarget _T2;
};
