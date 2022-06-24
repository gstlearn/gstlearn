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
