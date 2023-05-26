/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include <Geometry/BiTargetCheckCell.hpp>
#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Enum/ENeigh.hpp"

#include "Neigh/ANeigh.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Space/SpaceTarget.hpp"

class Db;

class GSTLEARN_EXPORT NeighCell: public ANeigh
{
public:
  NeighCell(bool flag_xvalid = false, int nmini = 1, const ASpace* space = nullptr);
  NeighCell(const NeighCell& r);
  NeighCell& operator=(const NeighCell& r);
  virtual ~NeighCell();

  /// Interface for ANeigh
  virtual int attach(const Db *dbin, const Db *dbout = nullptr) override;
  virtual VectorInt getNeigh(int iech_out) override;
  virtual bool hasChanged(int iech_out) const override;
  virtual int getMaxSampleNumber(const Db* db) const override { return 0; }
  virtual ENeigh getType() const override { return ENeigh::fromKey("CELL"); }

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static NeighCell* create(bool flag_xvalid = false,
                           int nmini = 1,
                           const ASpace *space = nullptr);
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

  mutable const DbGrid* _dbgrid;
  mutable SpaceTarget _T1;
  mutable SpaceTarget _T2;
};
