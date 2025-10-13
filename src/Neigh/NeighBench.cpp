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
#include "Neigh/NeighBench.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"

namespace gstlrn
{
NeighBench::NeighBench(bool flag_xvalid,
                       double width,
                       bool useBallTree,
                       Id leaf_size,
                       const ASpaceSharedPtr& space)
  : ANeigh(space)
  , _width(width)
  , _biPtBench()
  , _T1(space)
  , _T2(space)
{
  setFlagXvalid(flag_xvalid);

  setBallSearch(useBallTree, leaf_size);

  _biPtBench = BiTargetCheckBench::create(-1, _width);
}

NeighBench::NeighBench(const NeighBench& r)
  : ANeigh(r)
  , _width(r._width)
  , _biPtBench(r._biPtBench->clone())
  , _T1(r._T1)
  , _T2(r._T2)
{
}

NeighBench& NeighBench::operator=(const NeighBench& r)
{
  if (this != &r)
  {
    ANeigh::operator=(r);
    _width     = r._width;
    _biPtBench = r._biPtBench->clone();
    _T1        = r._T1;
    _T2        = r._T2;
  }
  return *this;
}

NeighBench::~NeighBench()
{
  delete _biPtBench;
}

Id NeighBench::attach(const Db* dbin, const Db* dbout)
{
  if (ANeigh::attach(dbin, dbout)) return 1;

  if (!_biPtBench->isValid(dbin, dbout)) return 1;

  return 0;
}

String NeighBench::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);

  std::stringstream sstr;

  sstr << toTitle(0, "Bench Neighborhood");

  sstr << _biPtBench->toString();

  return sstr.str();
}

bool NeighBench::_deserializeAscii(std::istream& is, bool verbose)
{
  double width = 0.;
  bool ret     = true;
  ret          = ret && ANeigh::_deserializeAscii(is, verbose);
  ret          = ret && _recordRead<double>(is, "Bench Width", width);

  _biPtBench = BiTargetCheckBench::create(-1, width); // idim_bench will be updated in 'attach'

  return ret;
}

bool NeighBench::_serializeAscii(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret      = ret && ANeigh::_serializeAscii(os, verbose);
  ret      = ret && _recordWrite<double>(os, "Bench Width", _biPtBench->getWidth());
  return ret;
}

NeighBench* NeighBench::create(bool flag_xvalid,
                               double width,
                               bool useBallTree,
                               Id leaf_size,
                               const ASpaceSharedPtr& space)
{
  return new NeighBench(flag_xvalid, width, useBallTree, leaf_size, space);
}

/**
 * Create a Neighborhood by loading the contents of a Neutral File
 * @param NFFilename Name of the Neutral File
 * @param verbose    Verbose flag
 * @return
 */
NeighBench* NeighBench::createFromNF(const String& NFFilename, bool verbose)
{
  auto* neigh = new NeighBench();
  if (neigh->_fileOpenAndDeserialize(NFFilename, verbose)) return neigh;
  delete neigh;
  return nullptr;
}

/**
 * Given a Db, returns the maximum number of samples per NeighBenchborhood
 * @param db Pointer to the target Db
 * @return
 */
Id NeighBench::getNSampleMax(const Db* db) const
{
  bool useSel = false;
  Id nech    = db->getNSample();
  Id ndim    = db->getNDim();
  if (db->getNDim() <= 2) return nech;

  /* Read the vector of the last coordinates */
  VectorDouble vec = db->getOneCoordinate(ndim - 1, useSel);

  /* Sort the third coordinate vector */
  VectorDouble tab = VH::sort(vec, true);

  /* Loop on the first point */
  Id nmax = 0;
  for (Id iech = 0; iech < nech - 1; iech++)
  {

    /* Loop on the second point */
    Id nloc = 1;
    for (Id jech = iech + 1; jech < nech; jech++)
    {
      if (ABS(tab[jech] - tab[iech]) > 2. * _width) break;
      nloc++;
    }

    /* Store the maximum number of samples */
    if (nloc > nmax) nmax = nloc;
  }
  return nmax;
}

bool NeighBench::hasChanged(Id iech_out) const
{
  if (_iechMemo < 0 || _isNbghMemoEmpty()) return true;

  return _isSameTargetBench(iech_out);
}

bool NeighBench::_isSameTargetBench(Id iech_out) const
{
  // Check if current target and previous target belong to the same bench

  Id ndim = _dbout->getNDim();
  if (_dbgrid != nullptr)
  {
    Id nval = 1;
    for (Id idim = 0; idim < ndim - 1; idim++)
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

/**
 * Select the neighborhood
 * @param iech_out Valid Rank of the sample in the output Db
 * @param ranks Vector of sample ranks in neighborhood (empty when error)
 */
void NeighBench::getNeigh(Id iech_out, VectorInt& ranks)
{
  // Select the neighborhood samples as the target sample has changed
  _bench(iech_out, ranks);

  // In case of debug option, dump out neighborhood characteristics
  if (OptDbg::query(EDbg::NBGH)) _display(ranks);

  // Compress the vector of returned sample ranks
  _neighCompress(ranks);
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
void NeighBench::_bench(Id iech_out, VectorInt& ranks)
{
  Id nech = _dbin->getNSample();
  ranks.resize(nech);
  ranks.fill(-1);

  // Load the target sample as a Space Target
  _dbout->getSampleAsSTInPlace(iech_out, _T1);

  /* Loop on samples */

  for (Id iech = 0; iech < nech; iech++)
  {
    /* Discard the masked input sample */

    if (!_dbin->isActive(iech)) continue;

    /* Discard samples where all variables are undefined */

    if (_discardUndefined(iech)) continue;

    /* Discard the target sample for the cross-validation option */

    if (getFlagXvalid())
    {
      if (_xvalid(iech, iech_out)) continue;
    }

    _dbin->getSampleAsSTInPlace(iech, _T2);

    /* Discard sample located outside the bench */

    if (!_biPtBench->isOK(_T1, _T2)) continue;

    ranks[iech] = 0;
  }
}

#ifdef HDF5
bool NeighBench::_deserializeH5(H5::Group& grp, [[maybe_unused]] bool verbose)
{
  auto neighG = SerializeHDF5::getGroup(grp, "NeighBench");
  if (!neighG)
  {
    return false;
  }

  /* Read the grid characteristics */
  bool ret     = true;
  double width = 0.;

  ret = ret && SerializeHDF5::readValue(*neighG, "Bench", width);

  ret = ret && ANeigh::_deserializeH5(*neighG, verbose);

  _biPtBench = BiTargetCheckBench::create(-1, width); // idim_bench will be updated in 'attach'

  return ret;
}

bool NeighBench::_serializeH5(H5::Group& grp, [[maybe_unused]] bool verbose) const
{
  auto neighG = grp.createGroup("NeighBench");

  bool ret = true;

  ret = ret && SerializeHDF5::writeValue(neighG, "Bench", _biPtBench->getWidth());
  ret = ret && ANeigh::_serializeH5(neighG, verbose);

  return ret;
}
#endif
} // namespace gstlrn