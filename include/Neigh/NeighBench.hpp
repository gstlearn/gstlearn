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

#include <Geometry/BiTargetCheckBench.hpp>
#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Enum/ENeigh.hpp"

#include "Neigh/ANeigh.hpp"
#include "Space/SpaceTarget.hpp"

class Db;

/**
 * \brief
 * Neighborhood definition by Bench.
 *
 * The Neighborhood is usually meant to select a sub-population from the input Data Base,
 * containing the active samples close to the target.
 *
 * The selected samples belong to the same 'bench' as the target: the distance (according
 * to the last space coordinate, e.g. the elevation in the 3-D case) between a selected sample
 * and the target is smaller than the bench width.
 *
 * The neighborhood also offers the possibility to suppress any sample which would be too close to (or coincide with)
 * the target: this is the cross-validation option.
 */
class GSTLEARN_EXPORT NeighBench: public ANeigh
{
public:
  NeighBench(bool flag_xvalid = false, double width = 0., const std::shared_ptr<const ASpace> &space = nullptr);
  NeighBench(const NeighBench& r);
  NeighBench& operator=(const NeighBench& r);
  virtual ~NeighBench();

  /// Interface for ANeigh
  virtual int attach(const Db *dbin, const Db *dbout = nullptr) override;
  virtual void getNeigh(int iech_out, VectorInt& ranks) override;
  virtual bool hasChanged(int iech_out) const override;
  virtual int getMaxSampleNumber(const Db* db) const override;
  virtual ENeigh getType() const override { return ENeigh::fromKey("BENCH"); }

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static NeighBench* create(bool flag_xvalid = false,
                            double width = 0,
                            const std::shared_ptr<const ASpace>& space = nullptr);
  static NeighBench* createFromNF(const String& neutralFilename, bool verbose = true);

  double getWidth() const { return _width; }

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "NeighBench"; }

private:
  bool _isSameTargetBench(int iech_out) const;
  void _bench(int iech_out, VectorInt& ranks);

private:
  double _width;
  BiTargetCheckBench* _biPtBench;

  mutable SpaceTarget _T1;
  mutable SpaceTarget _T2;
};
