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

#include <Geometry/BiTargetCheckBench.hpp>
#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Enum/ENeigh.hpp"

#include "Neigh/ANeigh.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Space/SpaceTarget.hpp"

class Db;

class GSTLEARN_EXPORT NeighBench: public ANeigh
{
public:
  NeighBench(bool flag_xvalid = false, double width = 0., const ASpace* space = nullptr);
  NeighBench(const NeighBench& r);
  NeighBench& operator=(const NeighBench& r);
  virtual ~NeighBench();

  /// Interface for ANeigh
  virtual int attach(const Db *dbin, const Db *dbout = nullptr) override;
  virtual VectorInt getNeigh(int iech_out) override;
  virtual bool hasChanged(int iech_out) const override;
  virtual int getMaxSampleNumber(const Db* db) const override;
  virtual ENeigh getType() const override { return ENeigh::fromKey("BENCH"); }

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static NeighBench* create(bool flag_xvalid = false,
                            double width = 0,
                            const ASpace *space = nullptr);
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
