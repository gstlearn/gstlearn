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

#include "Neigh/NeighImage.hpp"
#include "Basic/OptDbg.hpp"
#include "Db/Db.hpp"

NeighImage::NeighImage(const VectorInt& radius, int skip, const ASpace* space)
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

NeighImage* NeighImage::create(const VectorInt& image, int skip, const ASpace* space)
{
  return new NeighImage(image, skip, space);
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
int NeighImage::getMaxSampleNumber(const Db* /*db*/) const
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
  int nech = _dbin->getSampleNumber();
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

