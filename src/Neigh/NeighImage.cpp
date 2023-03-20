/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "geoslib_old_f.h"

#include "Neigh/NeighImage.hpp"
#include "Morpho/Morpho.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Db/Db.hpp"

NeighImage::NeighImage(const VectorInt& radius, int skip, const ASpace* space)
    : ANeighParam(false, space),
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
  ret = ret && ANeighParam::_serialize(os, verbose);
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
 * @param db Pointer to the taregt Db
 * @return
 */
int NeighImage::getMaxSampleNumber(const Db* /*db*/) const
{
  int nmax = 1;
  for (int idim = 0; idim < (int) getNDim(); idim++)
    nmax *= (2 * _imageRadius[idim] + 1);
  return nmax;
}
