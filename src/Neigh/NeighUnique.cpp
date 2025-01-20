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
#include "Neigh/NeighUnique.hpp"
#include "Mesh/AMesh.hpp"
#include "Basic/OptDbg.hpp"
#include "Db/Db.hpp"

NeighUnique::NeighUnique(bool flag_xvalid,  std::shared_ptr<const ASpace> space)
    : ANeigh(space)
{
  setFlagXvalid(flag_xvalid);
}

NeighUnique::NeighUnique(const NeighUnique& r)
    : ANeigh(r)
{
}

NeighUnique& NeighUnique::operator=(const NeighUnique& r)
{
  if (this != &r)
  {
    ANeigh::operator=(r);
   }
  return *this;
}

NeighUnique::~NeighUnique()
{
}

String NeighUnique::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);
  std::stringstream sstr;

  sstr << toTitle(0,"Unique Neighborhood");

  return sstr.str();
}

bool NeighUnique::_deserialize(std::istream& is, bool verbose)
{
  bool ret = true;
  ret = ret && ANeigh::_deserialize(is, verbose);
  return ret;
}

bool NeighUnique::_serialize(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret = ret && ANeigh::_serialize(os, verbose);
  return ret;
}

NeighUnique* NeighUnique::create(bool flag_xvalid, std::shared_ptr<const ASpace> space)
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
 * @param db Pointer to the target Db
 * @return
 */
int NeighUnique::getMaxSampleNumber(const Db* db) const
{
  return db->getSampleNumber(true);
}

bool NeighUnique::hasChanged(int iech_out) const
{
  DECLARE_UNUSED(iech_out);
  return (_iechMemo < 0 || _isNbghMemoEmpty());
}

/**
 * Select the neighborhood
 * @param iech_out Valid Rank of the sample in the output Db
 * @param ranks Vector of sample ranks in neighborhood (empty when error)
 */
void NeighUnique::getNeigh(int iech_out, VectorInt& ranks)
{
  int nech = _dbin->getSampleNumber();
  ranks.resize(nech);
  ranks.fill(-1);

  // Select the neighborhood samples as the target sample has changed
  _unique(iech_out, ranks);

  // In case of debug option, dump out neighborhood characteristics
  if (OptDbg::query(EDbg::NBGH)) _display(ranks);

  // Compress the vector of returned sample ranks
  _neighCompress(ranks);
}

/****************************************************************************/
/*!
 **  Select the unique neighborhood (or Image Neighborhood)
 **
 ** \param[in]  iech_out  rank of the output sample
 **
 ** \param[out]  ranks   Vector of samples elected in the Neighborhood
 **
 *****************************************************************************/
void NeighUnique::_unique(int iech_out, VectorInt& ranks)
{
  int nech = _dbin->getSampleNumber();

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
    ranks[iech] = 0;
  }
}

