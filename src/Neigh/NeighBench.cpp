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
#include "geoslib_old_f.h"

#include "Neigh/NeighBench.hpp"
#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Db/Db.hpp"

NeighBench::NeighBench(bool flag_xvalid, double width, const ASpace* space)
    : ANeighParam(flag_xvalid, space),
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

String NeighBench::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  sstr << toTitle(0,"Bench Neighborhood");
  sstr << ANeighParam::toString(strfmt);

  sstr << "Bench width     = " << _width << std::endl;

  return sstr.str();
}

bool NeighBench::_deserialize(std::istream& is, bool verbose)
{
  bool ret = true;
  ret = ret && ANeighParam::_deserialize(is, verbose);
  ret = ret && _recordRead<double>(is, "Bench Width", _width);
  return ret;
}

bool NeighBench::_serialize(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret = ret && ANeighParam::_serialize(os, verbose);
  ret = ret && _recordWrite<double>(os, "Bench Width", getWidth());
  return ret;
}

NeighBench* NeighBench::create(bool flag_xvalid, double width, const ASpace* space)
{
  return new NeighBench(flag_xvalid, width, space);
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
  neigh = new NeighBench();
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
