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
#include "Basic/SerializeHDF5.hpp"
#include "Db/Db.hpp"
#include "Space/ASpace.hpp"

namespace gstlrn
{
NeighUnique::NeighUnique(bool flag_xvalid,  const ASpaceSharedPtr& space)
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

bool NeighUnique::_deserializeAscii(std::istream& is, bool verbose)
{
  bool ret = true;
  ret = ret && ANeigh::_deserializeAscii(is, verbose);
  return ret;
}

bool NeighUnique::_serializeAscii(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret = ret && ANeigh::_serializeAscii(os, verbose);
  return ret;
}

NeighUnique* NeighUnique::create(bool flag_xvalid, const ASpaceSharedPtr& space)
{
  return new NeighUnique(flag_xvalid, space);
}

/**
 * Create a NeighUniqueborhood by loading the contents of a Neutral File
 * @param NFFilename Name of the Neutral File
 * @param verbose    Verbose flag
 * @return
 */
NeighUnique* NeighUnique::createFromNF(const String& NFFilename, bool verbose)
{
  auto* neigh = new NeighUnique;
  if (neigh->_fileOpenAndDeserialize(NFFilename, verbose)) return neigh;
  delete neigh;
  return nullptr;
}

/**
 * Given a Db, returns the maximum number of samples per NeighUniqueborhood
 * @param db Pointer to the target Db
 * @return
 */
Id NeighUnique::getNSampleMax(const Db* db) const
{
  return db->getNSample(true);
}

bool NeighUnique::hasChanged(Id iech_out) const
{
  DECLARE_UNUSED(iech_out);
  return (_iechMemo < 0 || _isNbghMemoEmpty());
}

/**
 * Select the neighborhood
 * @param iech_out Valid Rank of the sample in the output Db
 * @param ranks Vector of sample ranks in neighborhood (empty when error)
 */
void NeighUnique::getNeigh(Id iech_out, VectorInt& ranks)
{
  Id nech = _dbin->getNSample();
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
void NeighUnique::_unique(Id iech_out, VectorInt& ranks)
{
  Id nech = _dbin->getNSample();

  /* Loop on samples */

  for (Id iech = 0; iech < nech; iech++)
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

#ifdef HDF5
bool NeighUnique::_deserializeH5(H5::Group& grp, [[maybe_unused]] bool verbose)
{
  auto neighG = SerializeHDF5::getGroup(grp, "NeighUnique");
  if (!neighG)
  {
    return false;
  }

  /* Read the grid characteristics */
  bool ret = true;

  ret = ret && ANeigh::_deserializeH5(*neighG, verbose);
  
  return ret;
}

bool NeighUnique::_serializeH5(H5::Group& grp, [[maybe_unused]] bool verbose) const
{
  auto neighG = grp.createGroup("NeighUnique");

  bool ret = true;

  /* Writing the tail of the file */

  ret = ret && ANeigh::_serializeH5(neighG, verbose);

  return ret;
}
#endif
}