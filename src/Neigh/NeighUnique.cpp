/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "Neigh/NeighUnique.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Db/Db.hpp"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

NeighUnique::NeighUnique(int ndim, bool flag_xvalid)
    : ANeighParam(ndim, flag_xvalid)
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

/**
 * Constructor of a Unique NeighUniqueborhood
 */
int NeighUnique::reset(int ndim, bool flag_xvalid)
{
  setNDim(ndim);
  setFlagXvalid(flag_xvalid);

  return 0;
}

String NeighUnique::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  sstr << ANeighParam::toString(strfmt);

  return sstr.str();
}

int NeighUnique::_deserialize(FILE* file, bool verbose)
{
  if (ANeighParam::_deserialize(file, verbose))
  {
    if (verbose)
      messerr("Problem reading from the Neutral File.");
    return 1;
  }
  return 0;
}

int NeighUnique::_serialize(FILE* file, bool verbose) const
{
  if (ANeighParam::_serialize(file, verbose))
  {
    if (verbose) messerr("Problem writing in the Neutral File.");
    return 1;
  }
  return 0;
}

NeighUnique* NeighUnique::create(int ndim, bool flag_xvalid)
{
  NeighUnique* neighU = new NeighUnique;
  if (neighU->reset(ndim, flag_xvalid))
  {
    messerr("Problem when creating Unique NeighUniqueborhood");
    delete neighU;
    neighU = nullptr;
  }
  return neighU;
}

int NeighUnique::dumpToNF(const String& neutralFilename, bool verbose) const
{
  FILE* file = _fileOpen(neutralFilename, "NeighUnique", "w", verbose);
  if (file == nullptr) return 1;

  if (_serialize(file, verbose))
  {
    if (verbose) messerr("Problem writing in the Neutral File.");
    _fileClose(file, verbose);
    return 1;
  }
  _fileClose(file, verbose);
  return 0;
}

/**
 * Create a NeighUniqueborhood by loading the contents of a Neutral File
 * @param neutralFilename Name of the Neutral File
 * @param verbose         Verbose flag
 * @return
 */
NeighUnique* NeighUnique::createFromNF(const String& neutralFilename, bool verbose)
{
  FILE* file = _fileOpen(neutralFilename, "NeighUnique", "r", verbose);
  if (file == nullptr) return nullptr;

  NeighUnique* neigh = new NeighUnique;
  if (neigh->_deserialize(file, verbose))
  {
    if (verbose) messerr("Problem reading the Neutral File.");
    delete neigh;
    neigh = nullptr;
  }
  _fileClose(file, verbose);
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
