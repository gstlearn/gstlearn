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

int NeighUnique::_deserialize(std::istream& is, bool verbose)
{
  if (ANeighParam::_deserialize(is, verbose))
  {
    if (verbose)
      messerr("Problem reading from the Neutral File.");
    return 1;
  }
  return 0;
}

int NeighUnique::_serialize(std::ostream& os, bool verbose) const
{
  if (ANeighParam::_serialize(os, verbose))
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
  std::ofstream os;
  int ret = 1;
  if (_fileOpenWrite(neutralFilename, "NeighUnique", os, verbose))
  {
    ret = _serialize(os, verbose);
    if (ret && verbose) messerr("Problem writing in the Neutral File.");
    os.close();
  }
  return ret;
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
  if (_fileOpenRead(neutralFilename, "NeighUnique", is, verbose))
  {
    neigh = new NeighUnique;
    if (neigh->_deserialize(is, verbose))
    {
      if (verbose) messerr("Problem reading the Neutral File.");
      delete neigh;
      neigh = nullptr;
    }
    is.close();
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
