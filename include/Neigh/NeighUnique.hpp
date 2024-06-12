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
#include "geoslib_define.h"

#include "Enum/ENeigh.hpp"

#include "Neigh/ANeigh.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"

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
  NeighUnique(bool flag_xvalid = false, const ASpace* space = nullptr);
  NeighUnique(const NeighUnique& r);
  NeighUnique& operator=(const NeighUnique& r);
  virtual ~NeighUnique();

  /// Interface for ANeigh
  virtual void getNeigh(int iech_out, VectorInt& ranks) override;
  virtual int getMaxSampleNumber(const Db* db) const override;
  virtual bool hasChanged(int iech_out) const override;
  virtual ENeigh getType() const override { return ENeigh::fromKey("UNIQUE"); }

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static NeighUnique* create(bool flag_xvalid = false, const ASpace* space = nullptr);
  static NeighUnique* createFromNF(const String& neutralFilename, bool verbose = true);

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "NeighUnique"; }

private:
  void _unique(int iech_out, VectorInt& ranks);
};
