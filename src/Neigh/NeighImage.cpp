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
#include "Neigh/NeighImage.hpp"
#include "Basic/Law.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Mesh/AMesh.hpp"

#include "geoslib_old_f.h"

namespace gstlrn
{
NeighImage::NeighImage(const VectorInt& radius, Id skip, const ASpaceSharedPtr& space)
  : ANeigh(space)
  , _skip(skip)
  , _imageRadius(radius)
{
}

NeighImage::NeighImage(const NeighImage& r)
  : ANeigh(r)
  , _skip(r._skip)
  , _imageRadius(r._imageRadius)
{
}

NeighImage& NeighImage::operator=(const NeighImage& r)
{
  if (this != &r)
  {
    ANeigh::operator=(r);
    _skip        = r._skip;
    _imageRadius = r._imageRadius;
  }
  return *this;
}

NeighImage::~NeighImage()
{
}

String NeighImage::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);
  std::stringstream sstr;

  sstr << toTitle(0, "Image Neighborhood");

  sstr << "Skipping factor = " << _skip << std::endl;
  sstr << toMatrix("Image radius :", VectorString(), VectorString(), true,
                   1, static_cast<Id>(getNDim()), _imageRadius);

  return sstr.str();
}

bool NeighImage::_deserializeAscii(std::istream& is, bool verbose)
{
  bool ret = true;

  ret = ret && ANeigh::_deserializeAscii(is, verbose);
  ret = ret && _recordRead<Id>(is, "Skipping factor", _skip);
  for (Id idim = 0; ret && idim < static_cast<Id>(getNDim()); idim++)
  {
    double loc_radius  = 0.;
    ret                = ret && _recordRead<double>(is, "Image NeighImageborhood Radius",
                                                    loc_radius);
    _imageRadius[idim] = static_cast<Id>(loc_radius);
  }

  return ret;
}

bool NeighImage::_serializeAscii(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret      = ret && ANeigh::_serializeAscii(os, verbose);
  ret      = ret && _recordWrite<Id>(os, "", getSkip());
  for (Id idim = 0; ret && idim < static_cast<Id>(getNDim()); idim++)
    ret = ret && _recordWrite<double>(os, "", static_cast<double>(getImageRadius(idim)));
  ret = ret && _commentWrite(os, "Image NeighImageborhood parameters");
  return ret;
}

NeighImage* NeighImage::create(const VectorInt& radius, Id skip, const ASpaceSharedPtr& space)
{
  return new NeighImage(radius, skip, space);
}

/**
 * Create a NeighImageborhood by loading the contents of a Neutral File
 * @param NFFilename Name of the Neutral File
 * @param verbose    Verbose flag
 * @return
 */
NeighImage* NeighImage::createFromNF(const String& NFFilename, bool verbose)
{
  auto* neigh = new NeighImage();
  if (neigh->_fileOpenAndDeserialize(NFFilename, verbose)) return neigh;
  delete neigh;
  return nullptr;
}

/**
 * Given a Db, returns the maximum number of samples per NeighImageborhood
 * @param db Pointer to the target Db
 * @return
 */
Id NeighImage::getNSampleMax(const Db* /*db*/) const
{
  Id nmax = 1;
  for (Id idim = 0; idim < static_cast<Id>(getNDim()); idim++)
    nmax *= (2 * _imageRadius[idim] + 1);
  return nmax;
}

bool NeighImage::hasChanged(Id iech_out) const
{
  DECLARE_UNUSED(iech_out);
  return (_iechMemo < 0 || _isNbghMemoEmpty());
}

/**
 * Select the neighborhood
 * @param iech_out Valid Rank of the sample in the output Db
 * @param ranks Vector of sample ranks in neighborhood (empty when error)
 */
void NeighImage::getNeigh(Id iech_out, VectorInt& ranks)
{
  Id nech = _dbin->getNSample();
  ranks.resize(nech);
  ranks.fill(-1);

  // Select the neighborhood samples as the target sample has changed
  _uimage(iech_out, ranks);

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
void NeighImage::_uimage(Id iech_out, VectorInt& ranks)
{
  Id nech = _dbin->getNSample();

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
    ranks[iech] = 0;
  }
}

/**
 * @brief Create a subgrid containing the minimum pattern for Image Neighborhood.
 * The output subgrid is "parallel" to the input 'dbgrid'.
 *
 * @param dbgrid Input DbGrid
 * @param seed   Seed used for random number generation
 * @return Pointer to the newly created DbGrid
 *
 * @remark When a sample is skipped ('using 'skip' Neighborhood parameter)
 * the value of the corresponding variable is set to 'TEST'.
 * @remark The center point can never be skipped.
 */
DbGrid* NeighImage::buildImageGrid(const DbGrid* dbgrid, Id seed) const
{
  DbGrid* dbsub = nullptr;

  double seuil = 1. / (1. + _skip);
  Id ndim      = dbgrid->getNDim();
  Id nvar      = dbgrid->getNLoc(ELoc::Z);

  /* Core allocation */

  VectorInt nx(ndim);
  Id nech = 1;
  for (Id i = 0; i < ndim; i++)
  {
    nx[i] = 2 * _imageRadius[i] + 1;
    nech *= nx[i];
  }

  law_set_random_seed(seed);
  VectorBool sel(nech);
  for (Id iech = 0; iech < nech; iech++) sel[iech] = (law_uniform(0., 1.) < seuil);
  sel[nech / 2] = 1.;

  VectorDouble tab(nech * nvar);
  Id iecr = 0;
  for (Id ivar = 0; ivar < nvar; ivar++)
    for (Id iech = 0; iech < nech; iech++) tab[iecr++] = (sel[iech]) ? 0. : TEST;

  /* Create the grid */

  dbsub = DbGrid::create(nx, dbgrid->getDXs(), dbgrid->getX0s(), dbgrid->getAngles());
  dbsub->addColumns(tab, "Test", ELoc::Z);

  /* Shift the origin */

  VectorDouble coor(ndim);
  dbsub->rankToCoordinatesInPlace(nech / 2, coor);
  for (Id i = 0; i < ndim; i++) dbsub->setX0(i, dbsub->getX0(i) - coor[i]);
  if (db_grid_define_coordinates(dbsub)) return dbsub;

  return dbsub;
}
#ifdef HDF5
bool NeighImage::_deserializeH5(H5::Group& grp, [[maybe_unused]] bool verbose)
{
  auto neighG = SerializeHDF5::getGroup(grp, "NeighImage");
  if (!neighG)
  {
    return false;
  }

  /* Read the grid characteristics */
  bool ret = true;
  Id skip  = 0;

  ret = ret && SerializeHDF5::readValue(*neighG, "Skip", skip);

  ret = ret && SerializeHDF5::readVec(*neighG, "Radius", _imageRadius);

  ret = ret && ANeigh::_deserializeH5(*neighG, verbose);

  return ret;
}

bool NeighImage::_serializeH5(H5::Group& grp, [[maybe_unused]] bool verbose) const
{
  auto neighG = grp.createGroup("NeighImage");

  bool ret = true;

  ret = ret && SerializeHDF5::writeValue(neighG, "Skip", getSkip());
  ret = ret && SerializeHDF5::writeVec(neighG, "Radius", getImageRadius());

  ret = ret && ANeigh::_serializeH5(neighG, verbose);

  return ret;
}
#endif
} // namespace gstlrn