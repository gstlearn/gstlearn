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
#include "Arrays/BImage.hpp"
#include "Arrays/BImageStringFormat.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/Utilities.hpp"

namespace gstlrn
{
BImage::BImage(const VectorInt& ndims)
    : AArray(ndims),
      _values()
{
  _update();
}

BImage::BImage(const BImage &r)
    : AArray(r),
      _values(r._values)
{

}

BImage& BImage::operator=(const BImage &r)
{
  if (this != &r)
  {
    AArray::operator=(r);
    _values = r._values;
  }
  return *this;
}

BImage::~BImage()
{

}

void BImage::init(const VectorInt& ndims)
{
  AArray::init(ndims);
  _update();
}

void BImage::_update()
{
  auto nchar = getAllocSize();
  _values.resize(nchar, 0);
}

Id BImage::getAllocSize() const
{
  auto npixels = getNPixels();
  if (npixels <= 0) return 0;
  Id nred = ((npixels - 1) / 8 + 1);
  return nred;
}

Id BImage::getAddress(Id i, Id j, Id k) const
{
  return ((i)+(getNDims(0)*((j)+getNDims(1)*(k))));
}

bool BImage::isInside(Id i, Id j, Id k) const
{
  if (i < 0 || i >= getNDims(0)) return false;
  if (j < 0 || j >= getNDims(1)) return false;
  if (k < 0 || k >= getNDims(2)) return false;
  return true;
}

bool BImage::getValue(Id i, Id j, Id k) const
{
  return (getBImage(i, j, k) & getOffset(i, j, k));
}

void BImage::setMaskoff(Id i, Id j, Id k)
{
  _values[_divide(i, j, k)] &= getMaskoff(i, j, k);
}

void BImage::setOffset(Id i, Id j, Id k)
{
  _values[_divide(i, j, k)] |= getOffset(i, j, k);
}

String BImage::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  sstr << AArray::toString(strfmt);

  const auto* bstrfmt = dynamic_cast<const BImageStringFormat*>(strfmt);
  if (getNDim() <= 3)
  {
    // Default values

    Id izmin = 0;
    auto izmax = getNDims(2);
    Id iymin = 0;
    auto iymax = getNDims(1);
    Id ixmin = 0;
    auto ixmax = getNDims(0);
    char zero = '0';
    char one  = '1';
    if (strfmt != nullptr)
    {
      izmin = bstrfmt->getIndMin(2);
      iymin = bstrfmt->getIndMin(1);
      ixmin = bstrfmt->getIndMin(0);
      izmax = bstrfmt->getIndMax(2);
      if (IFFFF(izmax)) izmax = getNDims(2);
      iymax = bstrfmt->getIndMax(1);
      if (IFFFF(iymax)) iymax = getNDims(1);
      ixmax = bstrfmt->getIndMax(0);
      if (IFFFF(ixmax)) ixmax = getNDims(0);

      zero = bstrfmt->getCharZero();
      one  = bstrfmt->getCharOne();
    }

    /* Loop on the levels */

    for (Id iz = izmin; iz < izmax; iz++)
    {
      if (getNDims(2) > 1)
        sstr << toTitle(2, "Level %d/%d", iz + 1, getNDims(2));
      else
        sstr << std::endl;

      /* Loop on the cells of the layer */

      sstr << "  ";
      for (Id ix = ixmin; ix < ixmax; ix++)
      {
        Id val = (ix + 1) % 10;
        sstr << val;
      }
      sstr << std::endl;

      for (Id iy = iymin; iy < iymax; iy++)
      {
        auto jy = getNDims(1) - iy - 1;
        sstr << (iy + 1) % 10 << " ";
        for (Id ix = ixmin; ix < ixmax; ix++)
        {
          auto val = getValue(ix, jy, iz);
          if (val)
            sstr << one;
          else
            sstr << zero;
        }
        sstr << std::endl;
      }
    }
  }
  return sstr.str();
}

unsigned char BImage::getOffset (Id i, Id j, Id k) const
{
  static unsigned char COffset[] = { 128, 64, 32, 16, 8, 4, 2, 1 };
  return COffset[_residu(i,j,k)];
}

unsigned char BImage::getMaskoff(Id i, Id j, Id k) const
{
  static unsigned char CMaskoff[] = { 127, 191, 223, 239, 247, 251, 253, 254 };
  return CMaskoff[_residu(i,j,k)];
}
}