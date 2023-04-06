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
#include "geoslib_old_f.h"

#include "Neigh/NeighUnique.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Db/Db.hpp"

NeighUnique::NeighUnique(bool flag_xvalid, const ASpace* space)
    : ANeighParam(flag_xvalid, space)
{
}

NeighUnique::NeighUnique(const NeighUnique& r)
    : ANeighParam(r)
{
}

NeighUnique& NeighUnique::operator=(const NeighUnique& r)
{
  if (this != &r)
  {
    ANeighParam::operator=(r);
   }
  return *this;
}

NeighUnique::~NeighUnique()
{
}

String NeighUnique::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  sstr << toTitle(0,"Unique Neighborhood");
  sstr << ANeighParam::toString(strfmt);

  return sstr.str();
}

bool NeighUnique::_deserialize(std::istream& is, bool verbose)
{
  bool ret = true;
  ret = ret && ANeighParam::_deserialize(is, verbose);
  return ret;
}

bool NeighUnique::_serialize(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret = ret && ANeighParam::_serialize(os, verbose);
  return ret;
}

NeighUnique* NeighUnique::create(bool flag_xvalid, const ASpace* space)
{
  return new NeighUnique(flag_xvalid, space);
}

/**
 * Create a NeighUniqueborhood by loading the contents of a Neutral File
 * @param neutralFilename Name of the Neutral File
 * @param verbose         Verbose flag
 * @return
 */
NeighUnique* NeighUnique::createFromNF(const String& neutralFilename, bool verbose)
{
  NeighUnique* neigh = nullptr;
  std::ifstream is;
  neigh = new NeighUnique;
  bool success = false;
  if (neigh->_fileOpenRead(neutralFilename, is, verbose))
  {
    success =  neigh->deserialize(is, verbose);
  }
  if (! success)
  {
    delete neigh;
    neigh = nullptr;
  }
  return neigh;
}

/**
 * Given a Db, returns the maximum number of samples per NeighUniqueborhood
 * @param db Pointer to the taregt Db
 * @return
 */
int NeighUnique::getMaxSampleNumber(const Db* db) const
{
  return db->getSampleNumber(true);
}
