/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Enum/ENeigh.hpp"

#include "Neigh/ANeighParam.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"

class Db;

class GSTLEARN_EXPORT NeighUnique: public ANeighParam
{
public:
  NeighUnique(bool flag_xvalid = false, const ASpace* space = nullptr);
  NeighUnique(const NeighUnique& r);
  NeighUnique& operator=(const NeighUnique& r);
  virtual ~NeighUnique();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  virtual ENeigh getType() const override { return ENeigh::fromKey("UNIQUE"); }
  virtual int getMaxSampleNumber(const Db* db) const override;

  static NeighUnique* create(bool flag_xvalid = false, const ASpace* space = nullptr);
  static NeighUnique* createFromNF(const String& neutralFilename, bool verbose = true);

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "NeighUnique"; }
};
