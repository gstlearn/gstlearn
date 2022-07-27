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

NeighImage::NeighImage(int ndim, const VectorInt& radius, int skip)
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

bool NeighImage::_deserialize(std::istream& is, bool verbose)
{
  bool ret = true;

  ret = ret && ANeighParam::_deserialize(is, verbose);
  ret = ret && _recordRead<int>(is, "Skipping factor", _skip);
  for (int idim = 0; ret && idim < getNDim(); idim++)
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
  ret = ret && ANeighParam::_serialize(os, verbose);
  ret = ret && _recordWrite<int>(os, "", getSkip());
  for (int idim = 0; ret && idim < getNDim(); idim++)
    ret = ret && _recordWrite<double>(os, "", (double) getImageRadius(idim));
  ret = ret && _commentWrite(os, "Image NeighImageborhood parameters");
  return ret;
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
