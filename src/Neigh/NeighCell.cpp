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
#include "Neigh/NeighCell.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"

namespace gstlrn
{
NeighCell::NeighCell(bool flag_xvalid, Id nmini, bool useBallTree, Id leaf_size, const ASpaceSharedPtr& space)
  : ANeigh(space)
  , _nMini(nmini)
  , _biPtCell()
  , _T1(space)
  , _T2(space)
{
  setFlagXvalid(flag_xvalid);

  setBallSearch(useBallTree, leaf_size);

  _biPtCell = BiTargetCheckCell::create();
}

NeighCell::NeighCell(const NeighCell& r)
  : ANeigh(r)
  , _nMini(r._nMini)
  , _biPtCell(r._biPtCell)
  , _T1(r._T1)
  , _T2(r._T2)
{
}

NeighCell& NeighCell::operator=(const NeighCell& r)
{
  if (this != &r)
  {
    ANeigh::operator=(r);
    _nMini    = r._nMini;
    _biPtCell = r._biPtCell;
    _T1       = r._T1;
    _T2       = r._T2;
  }
  return *this;
}

NeighCell::~NeighCell()
{
}

Id NeighCell::attach(const Db* dbin, const Db* dbout)
{
  if (ANeigh::attach(dbin, dbout)) return 1;
  if (!_biPtCell->isValid(dbin, dbout)) return 1;

  _dbgrid = dynamic_cast<const DbGrid*>(dbout);
  return 0;
}

String NeighCell::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);
  std::stringstream sstr;

  sstr << toTitle(0, "Cell Neighborhood");

  if (_biPtCell != nullptr)
    sstr << _biPtCell->toString();

  return sstr.str();
}

bool NeighCell::_deserializeAscii(std::istream& is, bool verbose)
{
  bool ret = true;
  ret      = ret && ANeigh::_deserializeAscii(is, verbose);
  ret      = ret && _recordRead<Id>(is, "Minimum Number of samples", _nMini);

  return ret;
}

bool NeighCell::_serializeAscii(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret      = ret && ANeigh::_serializeAscii(os, verbose);
  ret      = ret && _recordWrite<Id>(os, "", getNMini());
  return ret;
}

NeighCell* NeighCell::create(bool flag_xvalid,
                             Id nmini,
                             bool useBallTree,
                             Id leaf_size,
                             const ASpaceSharedPtr& space)
{
  return new NeighCell(flag_xvalid, nmini, useBallTree, leaf_size, space);
}

/**
 * Create a Neighborhood by loading the contents of a Neutral File
 * @param NFFilename Name of the Neutral File
 * @param verbose    Verbose flag
 * @return
 */
NeighCell* NeighCell::createFromNF(const String& NFFilename, bool verbose)
{
  auto* neigh = new NeighCell();
  if (neigh->_fileOpenAndDeserialize(NFFilename, verbose)) return neigh;
  delete neigh;
  return nullptr;
}

bool NeighCell::hasChanged(Id iech_out) const
{
  DECLARE_UNUSED(iech_out);

  if (_iechMemo < 0 || _isNbghMemoEmpty()) return true;

  return true;
}

/**
 * Select the neighborhood
 * @param iech_out Valid Rank of the sample in the output Db
 * @param ranks Vector of sample ranks in neighborhood (empty when error)
 */
void NeighCell::getNeigh(Id iech_out, VectorInt& ranks)
{
  // Select the neighborhood samples as the target sample has changed
  if (_cell(iech_out, ranks))
  {
    ranks.clear();
    return;
  }

  // In case of debug option, dump out neighborhood characteristics
  if (OptDbg::query(EDbg::NBGH)) _display(ranks);

  // Compress the vector of returned sample ranks
  _neighCompress(ranks);
}

/****************************************************************************/
/*!
 **  Search for the cell neighborhood
 **
 ** \return Error returned code
 **
 ** \param[in]  iech_out  rank of the output sample
 **
 ** \param[out]  ranks    Vector of samples elected in the Neighborhood
 **
 *****************************************************************************/
Id NeighCell::_cell(Id iech_out, VectorInt& ranks)
{
  Id nech = _dbin->getNSample();
  ranks.resize(nech);
  ranks.fill(-1);

  // Load the target sample as a Space Target
  _dbgrid->getSampleAsSTInPlace(iech_out, _T1);

  /* Loop on samples */

  Id nsel = 0;
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

    if (!_biPtCell->isOK(_T1, _T2)) continue;

    ranks[iech] = 0;
    nsel++;
  }

  return (nsel < getNMini());
}

#ifdef HDF5
bool NeighCell::_deserializeH5(H5::Group& grp, [[maybe_unused]] bool verbose)
{
  auto neighG = SerializeHDF5::getGroup(grp, "NeighCell");
  if (!neighG)
  {
    return false;
  }

  /* Read the grid characteristics */
  bool ret = true;

  ret = ret && SerializeHDF5::readValue(*neighG, "NMini", _nMini);

  ret = ret && ANeigh::_deserializeH5(*neighG, verbose);

  return ret;
}

bool NeighCell::_serializeH5(H5::Group& grp, [[maybe_unused]] bool verbose) const
{
  auto neighG = grp.createGroup("NeighCell");

  bool ret = true;
  ret      = ret && SerializeHDF5::writeValue(neighG, "NMini", getNMini());

  ret = ret && ANeigh::_serializeH5(neighG, verbose);

  return ret;
}
#endif
} // namespace gstlrn