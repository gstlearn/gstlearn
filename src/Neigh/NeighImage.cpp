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
#include "Neigh/NeighImage.hpp"
#include "Morpho/Morpho.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Db/Db.hpp"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

NeighImage::NeighImage(int ndim, VectorInt radius, int skip)
    : ANeighParam(ndim),
      _skip(skip),
      _imageRadius(radius)
{
}

NeighImage::NeighImage(const NeighImage& r)
    : ANeighParam(r),
      _skip(r._skip),
      _imageRadius(r._imageRadius)
{
}

NeighImage& NeighImage::operator=(const NeighImage& r)
{
  if (this != &r)
  {
    ANeighParam::operator=(r);
    _skip = r._skip;
    _imageRadius = r._imageRadius;
   }
  return *this;
}

NeighImage::~NeighImage()
{
}

int NeighImage::reset(int ndim, const VectorInt& image, int skip)
{
  setNDim(ndim);
  _skip = skip;

  int radius = 1;
  if (! image.empty()) radius = image[0];
  if ((int) image.size() != ndim)
    _imageRadius.resize(ndim,radius);
  else
    _imageRadius = image;
  return 0;
}

String NeighImage::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  sstr << toTitle(0,"Image Neighborhood");
  sstr << ANeighParam::toString(strfmt);

  sstr << "Skipping factor = " << _skip << std::endl;
  sstr << toMatrix("Image radius :", VectorString(), VectorString(), true, getNDim(),
                  1, _imageRadius);

  return sstr.str();
}

int NeighImage::_deserialize(std::istream& is, bool verbose)
{
  if (ANeighParam::_deserialize(is, verbose))
  {
    if (verbose)
      messerr("Problem reading from the Neutral File.");
    return 1;
  }

  bool ret = _recordRead<int>(is, "Skipping factor", _skip);
  for (int idim = 0; idim < getNDim(); idim++)
  {
    double loc_radius = 0.;
    ret = ret && _recordRead<double>(is, "Image NeighImageborhood Radius", loc_radius);
    _imageRadius[idim] = static_cast<int> (loc_radius);
  }

  if (! ret) return 1;
  return 0;
}

int NeighImage::_serialize(std::ostream& os, bool verbose) const
{
  if (ANeighParam::_serialize(os, verbose))
  {
    if (verbose) messerr("Problem writing in the Neutral File.");
    return 1;
  }

  bool ret = _recordWrite<int>(os, "", getSkip());
  for (int idim = 0; idim < getNDim(); idim++)
    ret = ret && _recordWrite<double>(os, "", (double) getImageRadius(idim));
  ret = ret && _commentWrite(os, "Image NeighImageborhood parameters");
  return ret ? 0 : 1;
}

NeighImage* NeighImage::create(int ndim, const VectorInt& image, int skip)
{
  NeighImage* neighI = new NeighImage;
  if (neighI->reset(ndim, image, skip))
  {
    messerr("Problem when creating Image NeighImageborhood");
    delete neighI;
    neighI = nullptr;
  }
  return neighI;
}

int NeighImage::dumpToNF(const String& neutralFilename, bool verbose) const
{
  std::ofstream os;
  int ret = 1;
  if (_fileOpenWrite(neutralFilename, "NeighImage", os, verbose))
  {
    ret = _serialize(os, verbose);
    if (ret && verbose) messerr("Problem writing in the Neutral File.");
    os.close();
  }
  return ret;
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
  if (_fileOpenRead(neutralFilename, "NeighImage", is, verbose))
  {
    neigh = new NeighImage();
    if (neigh->_deserialize(is, verbose))
    {
      if (verbose) messerr("Problem reading the Neutral File.");
      delete neigh;
      neigh = nullptr;
    }
    is.close();
  }
  return neigh;
}

/**
 * Given a Db, returns the maximum number of samples per NeighImageborhood
 * @param db Pointer to the taregt Db
 * @return
 */
int NeighImage::getMaxSampleNumber(const Db* /*db*/) const
{
  int nmax = 1;
  for (int idim = 0; idim < getNDim(); idim++)
    nmax *= (2 * _imageRadius[idim] + 1);
  return nmax;
}
