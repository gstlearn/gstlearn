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
#include "Basic/OptDbg.hpp"
#include "Basic/Law.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Mesh/AMesh.hpp"

#include "geoslib_old_f.h"

NeighImage::NeighImage(const VectorInt& radius, int skip, const ASpaceSharedPtr& space)
    : ANeigh(space),
      _skip(skip),
      _imageRadius(radius)
{
}

NeighImage::NeighImage(const NeighImage& r)
    : ANeigh(r),
      _skip(r._skip),
      _imageRadius(r._imageRadius)
{
}

NeighImage& NeighImage::operator=(const NeighImage& r)
{
  if (this != &r)
  {
    ANeigh::operator=(r);
    _skip = r._skip;
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

  sstr << toTitle(0,"Image Neighborhood");

  sstr << "Skipping factor = " << _skip << std::endl;
  sstr << toMatrix("Image radius :", VectorString(), VectorString(), true,
                   1, getNDim(), _imageRadius);

  return sstr.str();
}

bool NeighImage::_deserialize(std::istream& is, bool verbose)
{
  bool ret = true;

  ret = ret && ANeigh::_deserialize(is, verbose);
  ret = ret && _recordRead<int>(is, "Skipping factor", _skip);
  for (int idim = 0; ret && idim < (int) getNDim(); idim++)
  {
    double loc_radius = 0.;
    ret = ret && _recordRead<double>(is, "Image NeighImageborhood Radius",
                                     loc_radius);
    _imageRadius[idim] = static_cast<int> (loc_radius);
  }

  return ret;
}

bool NeighImage::_serialize(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret = ret && ANeigh::_serialize(os, verbose);
  ret = ret && _recordWrite<int>(os, "", getSkip());
  for (int idim = 0; ret && idim < (int) getNDim(); idim++)
    ret = ret && _recordWrite<double>(os, "", (double) getImageRadius(idim));
  ret = ret && _commentWrite(os, "Image NeighImageborhood parameters");
  return ret;
}

NeighImage* NeighImage::create(const VectorInt& radius, int skip, const ASpaceSharedPtr& space)
{
  return new NeighImage(radius, skip, space);
}

/**
 * Create a NeighImageborhood by loading the contents of a Neutral File
 * @param neutralFilename Name of the Neutral File
 * @param verbose         Verbose flag
 * @return
 */
NeighImage* NeighImage::createFromNF(const String& neutralFilename, bool verbose)
{
  NeighImage* neigh = nullptr;
  std::ifstream is;
  neigh = new NeighImage();
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
 * Given a Db, returns the maximum number of samples per NeighImageborhood
 * @param db Pointer to the target Db
 * @return
 */
int NeighImage::getNSampleMax(const Db* /*db*/) const
{
  int nmax = 1;
  for (int idim = 0; idim < (int) getNDim(); idim++)
    nmax *= (2 * _imageRadius[idim] + 1);
  return nmax;
}

bool NeighImage::hasChanged(int iech_out) const
{
  DECLARE_UNUSED(iech_out);
  return (_iechMemo < 0 || _isNbghMemoEmpty());
}

/**
 * Select the neighborhood
 * @param iech_out Valid Rank of the sample in the output Db
 * @param ranks Vector of sample ranks in neighborhood (empty when error)
 */
void NeighImage::getNeigh(int iech_out, VectorInt& ranks)
{
  int nech = _dbin->getNSample();
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
void NeighImage::_uimage(int iech_out, VectorInt& ranks)
{
  int nech = _dbin->getNSample();

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
DbGrid* NeighImage::buildImageGrid(const DbGrid* dbgrid, int seed) const
{
  DbGrid* dbsub = nullptr;

  double seuil = 1. / (1. + _skip);
  int ndim     = dbgrid->getNDim();
  int nvar     = dbgrid->getNLoc(ELoc::Z);

  /* Core allocation */

  VectorInt nx(ndim);
  int nech = 1;
  for (int i = 0; i < ndim; i++)
  {
    nx[i] = 2 * _imageRadius[i] + 1;
    nech *= nx[i];
  }

  law_set_random_seed(seed);
  VectorBool sel(nech);
  for (int iech = 0; iech < nech; iech++) sel[iech] = (law_uniform(0., 1.) < seuil);
  sel[nech / 2] = 1.;

  VectorDouble tab(nech * nvar);
  int iecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int iech = 0; iech < nech; iech++) tab[iecr++] = (sel[iech]) ? 0. : TEST;

  /* Create the grid */

  dbsub = DbGrid::create(nx, dbgrid->getDXs(), dbgrid->getX0s(), dbgrid->getAngles());
  dbsub->addColumns(tab, "Test", ELoc::Z);

  /* Shift the origin */

  VectorDouble coor(ndim);
  dbsub->rankToCoordinatesInPlace(nech / 2, coor);
  for (int i = 0; i < ndim; i++) dbsub->setX0(i, dbsub->getX0(i) - coor[i]);
  if (db_grid_define_coordinates(dbsub)) return dbsub;

  return dbsub;
}
