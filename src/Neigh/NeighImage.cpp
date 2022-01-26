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

NeighImage::NeighImage(int ndim, int skip, VectorInt image)
    : ANeighParam(ndim),
      _skip(skip),
      _imageRadius(image)
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

int NeighImage::reset(int ndim, int skip, const VectorInt& image)
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

  sstr << ANeighParam::toString(strfmt);

  sstr << "Skipping factor = " << _skip << std::endl;
  sstr << toMatrix("Image radius :", VectorString(), VectorString(), true, getNDim(),
                  1, _imageRadius);

  return sstr.str();
}

int NeighImage::_deserialize(FILE* file, bool verbose)
{
  if (ANeighParam::_deserialize(file, verbose))
  {
    if (verbose)
      messerr("Problem reading from the Neutral File.");
    return 1;
  }

  if (_recordRead(file, "Skipping factor", "%d", &_skip)) return 1;
  for (int idim = 0; idim < getNDim(); idim++)
  {
    double loc_radius;
    if (_recordRead(file, "Image NeighImageborhood Radius", "%lf", &loc_radius)) return 1;
    _imageRadius[idim] = static_cast<int> (loc_radius);
  }
  return 0;
}

int NeighImage::_serialize(FILE* file, bool verbose) const
{
  if (ANeighParam::_serialize(file, verbose))
  {
    if (verbose) messerr("Problem writing in the Neutral File.");
    return 1;
  }

  _recordWrite(file, "%d", getSkip());
  for (int idim = 0; idim < getNDim(); idim++)
    _recordWrite(file, "%lf", (double) getImageRadius(idim));
  _recordWrite(file, "#", "Image NeighImageborhood parameters");
  return 0;
}

NeighImage* NeighImage::create(int ndim, int skip, const VectorInt& image)
{
  NeighImage* neighI = new NeighImage;
  if (neighI->reset(ndim, skip, image))
  {
    messerr("Problem when creating Image NeighImageborhood");
    delete neighI;
    neighI = nullptr;
  }
  return neighI;
}

int NeighImage::dumpToNF(const String& neutralFilename, bool verbose) const
{
  FILE* file = _fileOpen(neutralFilename, "NeighImage", "w", verbose);
  if (file == nullptr) return 1;

  if (_serialize(file, verbose))
  {
    if (verbose) messerr("Problem writing in the Neutral File.");
    _fileClose(file, verbose);
    return 1;
  }
  _fileClose(file, verbose);
  return 0;
}

/**
 * Create a NeighImageborhood by loading the contents of a Neutral File
 * @param neutralFilename Name of the Neutral File
 * @param verbose         Verbose flag
 * @return
 */
NeighImage* NeighImage::createFromNF(const String& neutralFilename, bool verbose)
{
  FILE* file = _fileOpen(neutralFilename, "NeighImage", "r", verbose);
  if (file == nullptr) return nullptr;

  NeighImage* neigh = new NeighImage();
  if (neigh->_deserialize(file, verbose))
  {
    if (verbose) messerr("Problem reading the Neutral File.");
    delete neigh;
    neigh = nullptr;
  }
  _fileClose(file, verbose);
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
