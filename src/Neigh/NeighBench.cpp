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

#include "Neigh/NeighBench.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"

NeighBench::NeighBench(bool flag_xvalid, double width, const ASpace *space)
    : ANeigh(space),
      _width(width),
      _biPtBench(),
      _T1(space),
      _T2(space)
{
  setFlagXvalid (flag_xvalid);
  _biPtBench = BiTargetCheckBench::create(-1, _width);
}

NeighBench::NeighBench(const NeighBench& r)
    : ANeigh(r),
      _width(r._width),
      _biPtBench(r._biPtBench),
      _T1(r._T1),
      _T2(r._T2)
{
}

NeighBench& NeighBench::operator=(const NeighBench& r)
{
  if (this != &r)
  {
    ANeigh::operator=(r);
    _width = r._width;
    _biPtBench = r._biPtBench;
    _T1 = r._T1;
    _T2 = r._T2;
   }
  return *this;
}

NeighBench::~NeighBench()
{
}

int NeighBench::attach(const Db *dbin, const Db *dbout)
{
  if (ANeigh::attach(dbin, dbout)) return 1;

  if (! _biPtBench->isValid(dbin, dbout)) return 1;

  return 0;
}

String NeighBench::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  sstr << toTitle(0,"Bench Neighborhood");

  sstr << _biPtBench->toString();

  return sstr.str();
}

bool NeighBench::_deserialize(std::istream& is, bool verbose)
{
  double width;
  bool ret = true;
  ret = ret && ANeigh::_deserialize(is, verbose);
  ret = ret && _recordRead<double>(is, "Bench Width", width);

  _biPtBench = BiTargetCheckBench::create(-1, width); // idim_bench will be updated in 'attach'

  return ret;
}

bool NeighBench::_serialize(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret = ret && ANeigh::_serialize(os, verbose);
  ret = ret && _recordWrite<double>(os, "Bench Width", _biPtBench->getWidth());
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
 * @param db Pointer to the target Db
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
  VectorDouble tab = VH::sort(vec, true);

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

bool NeighBench::hasChanged(int iech_out) const
{
  if (_iechMemo < 0 || _isNbghMemoEmpty()) return true;

  return _isSameTargetBench(iech_out);
}

bool NeighBench::_isSameTargetBench(int iech_out) const
{
  // Check if current target and previous target belong to the same bench

  int ndim = _dbout->getNDim();
  if (_dbgrid != nullptr)
  {
    int nval = 1;
    for (int idim = 0; idim < ndim - 1; idim++)
      nval *= _dbgrid->getNX(idim);
    if ((iech_out / nval) != (_iechMemo / nval)) return false;
  }
  else
  {
    if (_dbout->getCoordinate(iech_out, ndim - 1) !=
        _dbout->getCoordinate(_iechMemo, ndim - 1)) return false;
  }
  return true;
}

/****************************************************************************/
/*!
 **  Select the neighborhood
 **
 ** \return  Vector of sample ranks in neighborhood (empty when error)
 **
 ** \param[in]  iech_out      Valid Rank of the sample in the output Db
 **
 *****************************************************************************/
VectorInt NeighBench::getNeigh(int iech_out)
{
  int nech = _dbin->getSampleNumber();
  VectorInt ranks(nech, -1);

  // Select the neighborhood samples as the target sample has changed
  _bench(iech_out, ranks);

  // In case of debug option, dump out neighborhood characteristics
  if (OptDbg::query(EDbg::NBGH)) _display(ranks);

  // Compress the vector of returned sample ranks
  _neighCompress(ranks);

  return ranks;
}

/****************************************************************************/
/*!
 **  Search for the bench neighborhood, according to the last
 **  coordinate
 **
 ** \param[in]  iech_out  rank of the output sample
 **
 ** \param[out]  ranks    Vector of samples elected in the Neighborhood
 **
 *****************************************************************************/
void NeighBench::_bench(int iech_out, VectorInt& ranks)
{
  int nech = _dbin->getSampleNumber();

  // Load the target sample as a Space Target
  _dbout->getSampleAsST(iech_out, _T1);

  /* Loop on samples */

  for (int iech = 0; iech < nech; iech++)
  {
    /* Discard the masked input sample */

    if (! _dbin->isActive(iech)) continue;

    /* Discard samples where all variables are undefined */

    if (_discardUndefined(iech)) continue;

    /* Discard the target sample for the cross-validation option */

    if (getFlagXvalid())
    {
      if (_xvalid(iech, iech_out)) continue;
    }

    _dbin->getSampleAsST(iech, _T2);

    /* Discard sample located outside the bench */

    if (! _biPtBench->isOK(_T1, _T2)) continue;

    ranks[iech] = 0;
  }
}

