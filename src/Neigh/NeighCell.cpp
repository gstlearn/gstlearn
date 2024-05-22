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

#include "Neigh/NeighCell.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"

NeighCell::NeighCell(bool flag_xvalid, int nmini, const ASpace *space)
    : ANeigh(space),
      _nMini(nmini),
      _biPtCell(),
      _dbgrid(),
      _T1(space),
      _T2(space)
{
  setFlagXvalid (flag_xvalid);

  _biPtCell = BiTargetCheckCell::create();
}

NeighCell::NeighCell(const NeighCell& r)
    : ANeigh(r),
      _nMini(r._nMini),
      _biPtCell(r._biPtCell),
      _dbgrid(r._dbgrid),
      _T1(r._T1),
      _T2(r._T2)
{
}

NeighCell& NeighCell::operator=(const NeighCell& r)
{
  if (this != &r)
  {
    ANeigh::operator=(r);
    _nMini = r._nMini;
    _biPtCell = r._biPtCell;
    _dbgrid = r._dbgrid;
    _T1 = r._T1;
    _T2 = r._T2;
   }
  return *this;
}

NeighCell::~NeighCell()
{
}

int NeighCell::attach(const Db *dbin, const Db *dbout)
{
  if (ANeigh::attach(dbin, dbout)) return 1;
  if (! _biPtCell->isValid(dbin, dbout)) return 1;

  _dbgrid = dynamic_cast<const DbGrid*>(dbout);
  return 0;
}

String NeighCell::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);
  std::stringstream sstr;

  sstr << toTitle(0,"Cell Neighborhood");

  if (_biPtCell != nullptr)
    sstr << _biPtCell->toString();

  return sstr.str();
}

bool NeighCell::_deserialize(std::istream& is, bool verbose)
{
  bool ret = true;
  ret = ret && ANeigh::_deserialize(is, verbose);
  ret = ret && _recordRead<int>(is, "Minimum Number of samples", _nMini);

  return ret;
}

bool NeighCell::_serialize(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret = ret && ANeigh::_serialize(os, verbose);
  ret = ret && _recordWrite<int>(os, "", getNMini());
  return ret;
}

NeighCell* NeighCell::create(bool flag_xvalid, int nmini, const ASpace* space)
{
  return new NeighCell(flag_xvalid, nmini, space);
}

/**
 * Create a Neighborhood by loading the contents of a Neutral File
 * @param neutralFilename Name of the Neutral File
 * @param verbose         Verbose flag
 * @return
 */
NeighCell* NeighCell::createFromNF(const String& neutralFilename, bool verbose)
{
  NeighCell* neigh = nullptr;
  std::ifstream is;
  neigh = new NeighCell();
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

bool NeighCell::hasChanged(int iech_out) const
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
void NeighCell::getNeigh(int iech_out, VectorInt& ranks)
{
  int nech = _dbin->getSampleNumber();
  ranks.resize(nech);
  ranks.fill(-1);

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
int NeighCell::_cell(int iech_out, VectorInt& ranks)
{
  int nech = _dbin->getSampleNumber();

  // Load the target sample as a Space Target
  _dbgrid->getSampleAsST(iech_out, _T1);

  /* Loop on samples */

  int nsel = 0;
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

    if (! _biPtCell->isOK(_T1, _T2)) continue;

    ranks[iech] = 0;
    nsel++;
  }

  return (nsel < getNMini());
}

