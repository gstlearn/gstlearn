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

#include "Db/Db.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/ICloneable.hpp"

/**
 * \brief
 * Class containing the Data Information organized as a set of Lines
 *
 * This class is derived from the Db class, with a specific decoration: samples
 * within the Db are organized sequenially along lines.
 *
 * Note that this particular Db does not allow the modification of the sample
 * number by addition or deletion.
 *
 * The Line decoration provides a vector giving the number of samples per line.
 * Within one line, samples are ordered sequentailly (the order must not be modified).
 */
class GSTLEARN_EXPORT DbLine: public Db
{
public:
  DbLine();
  DbLine(const DbLine& r);
  DbLine& operator=(const DbLine& r);
  virtual ~DbLine();

public:
  /// ICloneable interface
  IMPLEMENT_CLONING(DbLine)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Db Interface
  inline bool isGrid() const override { return false; }
  bool mayChangeSampleNumber() const override { return false; }

  static DbLine* createFromNF(const String& neutralFilename,
                              bool verbose = true);

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "DbLine"; }
};
