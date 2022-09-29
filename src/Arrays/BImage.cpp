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
#include "Arrays/BImage.hpp"
#include "Arrays/BImageStringFormat.hpp"
#include "Basic/Vector.hpp"
#include "Basic/Utilities.hpp"

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
  int nchar = getAllocSize();
  _values.resize(nchar, 0);
}

int BImage::getAllocSize() const
{
  int npixels = getNPixels();
  if (npixels <= 0) return 0;
  int nred = ((npixels - 1) / 8 + 1);
  return nred;
}

int BImage::getAddress(int i, int j, int k) const
{
  return ((i)+(getNDims(0)*((j)+getNDims(1)*(k))));
}

bool BImage::isInside(int i, int j, int k) const
{
  if (i < 0 || i >= getNDims(0)) return false;
  if (j < 0 || j >= getNDims(1)) return false;
  if (k < 0 || k >= getNDims(2)) return false;
  return true;
}

bool BImage::getValue(int i, int j, int k) const
{
  return (getBImage(i, j, k) & getOffset(i, j, k));
}

void BImage::setMaskoff(int i, int j, int k)
{
  _values[_divide(i, j, k)] &= getMaskoff(i, j, k);
}

void BImage::setOffset(int i, int j, int k)
{
  _values[_divide(i, j, k)] |= getOffset(i, j, k);
}

String BImage::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  sstr << AArray::toString(strfmt);

  const BImageStringFormat* bstrfmt = dynamic_cast<const BImageStringFormat*>(strfmt);
  if (getNDim() <= 3)
  {
    // Default values

    int izmin = 0;
    int izmax = getNDims(2);
    int iymin = 0;
    int iymax = getNDims(1);
    int ixmin = 0;
    int ixmax = getNDims(0);
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

    for (int iz = izmin; iz < izmax; iz++)
    {
      if (getNDims(2) > 1)
        sstr << toTitle(2, "Level %d/%d", iz + 1, getNDims(2));
      else
        sstr << std::endl;

      /* Loop on the cells of the layer */

      sstr << "  ";
      for (int ix = ixmin; ix < ixmax; ix++)
      {
        int val = (ix + 1) % 10;
        sstr << val;
      }
      sstr << std::endl;

      for (int iy = iymin; iy < iymax; iy++)
      {
        int jy = getNDims(1) - iy - 1;
        sstr << (iy + 1) % 10 << " ";
        for (int ix = ixmin; ix < ixmax; ix++)
        {
          int val = getValue(ix, jy, iz);
          if (val > 0)
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

unsigned char BImage::getOffset (int i, int j, int k) const
{
  static unsigned char COffset[] = { 128, 64, 32, 16, 8, 4, 2, 1 };
  return COffset[_residu(i,j,k)];
}

unsigned char BImage::getMaskoff(int i, int j, int k) const
{
  static unsigned char CMaskoff[] = { 127, 191, 223, 239, 247, 251, 253, 254 };
  return CMaskoff[_residu(i,j,k)];
}
