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
#include "Neigh/NeighBench.hpp"
#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Db/Db.hpp"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

NeighBench::NeighBench(int ndim, bool flag_xvalid, double width)
    : ANeighParam(ndim, flag_xvalid),
      _width(width)
{
}

NeighBench::NeighBench(const NeighBench& r)
    : ANeighParam(r),
      _width(r._width)
{
}

NeighBench& NeighBench::operator=(const NeighBench& r)
{
  if (this != &r)
  {
    ANeighParam::operator=(r);
    _width = r._width;
   }
  return *this;
}

NeighBench::~NeighBench()
{
}

int NeighBench::reset(int ndim, bool flag_xvalid, double width)
{
  setNDim(ndim);
  setFlagXvalid(flag_xvalid);

  _width = width;
  return 0;
}

String NeighBench::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  sstr << ANeighParam::toString(strfmt);

  sstr << "Bench width     = " << _width << std::endl;

  return sstr.str();
}

int NeighBench::_deserialize(std::istream& is, bool verbose)
{
  if (ANeighParam::_deserialize(is, verbose))
  {
    if (verbose)
      messerr("Problem reading from the Neutral File.");
    return 1;
  }

  bool ret = _recordRead<double>(is, "Bench Width", _width);

  if (! ret) return 1;
  return 0;
}

int NeighBench::_serialize(std::ostream& os, bool verbose) const
{
  if (ANeighParam::_serialize(os, verbose))
   {
     if (verbose) messerr("Problem writing in the Neutral File.");
     return 1;
   }

  bool ret = _recordWrite<double>(os, "Bench Width", getWidth());

  return ret ? 0 : 1;
}

NeighBench* NeighBench::create(int ndim, bool flag_xvalid, double width)
{
  NeighBench* neighB = new NeighBench;
  if (neighB->reset(ndim, flag_xvalid, width))
  {
    messerr("Problem when creating Moving NeighBenchborhood");
    delete neighB;
    neighB =  nullptr;
  }
  return neighB;
}

int NeighBench::dumpToNF(const String& neutralFilename, bool verbose) const
{
  std::ofstream os;
  int ret = 1;
  if (_fileOpenWrite(neutralFilename, "NeighBench", os, verbose))
  {
    ret = _serialize(os, verbose);
    if (ret && verbose) messerr("Problem writing in the Neutral File.");
    os.close();
  }
  return ret;
}

/**
 * Create a Neighborhood by loading the contents of a Neutral File
 * @param neutralFilename Name of the Neutral File
 * @param verbose         Verbose flag
 * @return
 */
NeighBench* NeighBench::createFromNF(const String& neutralFilename, bool verbose)
{
  NeighBench* neigh = nullptr;
  std::ifstream is;
  if (_fileOpenRead(neutralFilename, "NeighBench", is, verbose))
  {
    neigh = new NeighBench();
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
 * Given a Db, returns the maximum number of samples per NeighBenchborhood
 * @param db Pointer to the taregt Db
 * @return
 */
int NeighBench::getMaxSampleNumber(const Db* db) const
{
  bool useSel = false;
  int nech = db->getSampleNumber();
  int nmax = nech;
  int ndim = db->getNDim();
  if (db->getNDim() <= 2) return nech;

  /* Read the vector of the last coordinates */
  VectorDouble vec = db->getCoordinates(ndim-1, useSel);

  /* Sort the third coordinate vector */
  VectorDouble tab = ut_vector_sort(vec, true);

  /* Loop on the first point */
  nmax = 0;
  for (int iech = 0; iech < nech - 1; iech++)
  {

    /* Loop on the second point */
    int nloc = 1;
    for (int jech = iech + 1; jech < nech; jech++)
    {
      if (ABS(tab[jech] - tab[iech]) > 2. * _width) break;
      nloc++;
    }

    /* Store the maximum number of samples */
    if (nloc > nmax) nmax = nloc;
  }
  return nmax;
}
