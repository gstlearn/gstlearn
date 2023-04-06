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

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Enum/ENeigh.hpp"

#include "Neigh/ANeighParam.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"

class Db;

// TODO : inherits from ASpaceObject (see _init)
class GSTLEARN_EXPORT NeighBench: public ANeighParam
{
public:
  NeighBench(bool flag_xvalid = false, double width = 0., const ASpace* space = nullptr);
  NeighBench(const NeighBench& r);
  NeighBench& operator=(const NeighBench& r);
  virtual ~NeighBench();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static NeighBench* create(bool flag_xvalid = false,
                            double width = 0,
                            const ASpace *space = nullptr);
  static NeighBench* createFromNF(const String& neutralFilename, bool verbose = true);

  virtual int getMaxSampleNumber(const Db* db) const override;
  virtual ENeigh getType() const override { return ENeigh::fromKey("BENCH"); }

  double getWidth() const { return _width; }
  void setWidth(double width) { _width = width; }

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "NeighBench"; }

private:
  double _width;                 /* Width of the slice - bench */
};
