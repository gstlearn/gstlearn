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
#include "geoslib_old_f.h"
#include "geoslib_define.h"

#include "Db/Db.hpp"
#include "Db/DbLine.hpp"
#include "Db/DbStringFormat.hpp"
#include "Polygon/Polygons.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/VectorNumT.hpp"
#include "Stats/Classical.hpp"

#include <math.h>

DbLine::DbLine()
    : Db()
{
  _clear();
}

DbLine::DbLine(const DbLine& r)
    : Db(r)
{
}

DbLine& DbLine::operator=(const DbLine& r)
{
  if (this != &r)
  {
    Db::operator=(r);
  }
  return *this;
}

DbLine::~DbLine()
{
}

String DbLine::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  const DbStringFormat* dbfmt = dynamic_cast<const DbStringFormat*>(strfmt);
  DbStringFormat dsf;
  if (dbfmt != nullptr) dsf = *dbfmt;

  sstr << toTitle(0, "Data Base Line Characteristics");

  sstr << _toStringCommon(&dsf);

  return sstr.str();
}

bool DbLine::_deserialize(std::istream& is, bool verbose)
{
  int ndim = 0;
  VectorString locators;
  VectorString names;
  VectorDouble values;
  VectorDouble allvalues;

  /* Initializations */

  bool ret = true;
  ret = ret && _recordRead<int>(is, "Space Dimension", ndim);

  ret = ret && Db::_deserialize(is, verbose);

  return ret;
}

bool DbLine::_serialize(std::ostream& os, bool verbose) const
{
  bool ret = true;

  /* Writing the header */

  ret = ret && _recordWrite<int>(os, "Space Dimension", getNDim());

  /* Writing the tail of the file */

  ret && Db::_serialize(os, verbose);

  return ret;
}

/**
 * Create a Db by loading the contents of a Neutral File
 *
 * @param neutralFilename Name of the Neutral File (Db format)
 * @param verbose         Verbose
 *
 * @remarks The name does not need to be completed in particular when defined by absolute path
 * @remarks or read from the Data Directory (in the gstlearn distribution)
 */
DbLine* DbLine::createFromNF(const String& neutralFilename, bool verbose)
{
  DbLine* dbLine = nullptr;
  std::ifstream is;
  dbLine = new DbLine;
  bool success = false;
  if (dbLine->_fileOpenRead(neutralFilename, is, verbose))
  {
    success = dbLine->deserialize(is, verbose);
  }
  if (! success)
  {
    delete dbLine;
    dbLine = nullptr;
  }
  return dbLine;
}
