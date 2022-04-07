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
#include "geoslib_old_f.h"

#include "OutputFormat/GridBmp.hpp"
#include "OutputFormat/AOF.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"

#include <string.h>
#include <stdio.h>

#define BF_TYPE 0x4D42             /* "MB" */
#define COLOR_MASK   -1
#define COLOR_FFFF   -2
#define COLOR_LOWER  -3
#define COLOR_UPPER  -4
#define N_SAMPLE(nx,nsample) ((int) ((nx-1) / nsample) + 1)

GridBmp::GridBmp(const char* filename, const Db* db)
  : AOF(filename, db)
  , _nsamplex(1)
  , _nsampley(1)
  , _nmult(1)
  , _ncolor(0)
  , _flag_low(true)
  , _flag_high(true)
  , _mask_red(0)
  , _mask_green(0)
  , _mask_blue(0)
  , _ffff_red(0)
  , _ffff_green(0)
  , _ffff_blue(0)
  , _low_red(0)
  , _low_green(0)
  , _low_blue(0)
  , _high_red(0)
  , _high_green(0)
  , _high_blue(0)
  , _valmin(TEST)
  , _valmax(TEST)
  , _reds()
  , _greens()
  , _blues()
{
}

GridBmp::GridBmp(const GridBmp& r)
    : AOF(r),
      _nsamplex(r._nsamplex),
      _nsampley(r._nsampley),
      _nmult(r._nmult),
      _ncolor(r._ncolor),
      _flag_low(r._flag_low),
      _flag_high(r._flag_high),
      _mask_red(r._mask_red),
      _mask_green(r._mask_green),
      _mask_blue(r._mask_blue),
      _ffff_red(r._ffff_red),
      _ffff_green(r._ffff_green),
      _ffff_blue(r._ffff_blue),
      _low_red(r._low_red),
      _low_green(r._low_green),
      _low_blue(r._low_blue),
      _high_red(r._high_red),
      _high_green(r._high_green),
      _high_blue(r._high_blue),
      _valmin(r._valmin),
      _valmax(r._valmax),
      _reds(r._reds),
      _greens(r._greens),
      _blues(r._blues)
{
}

GridBmp& GridBmp::operator=(const GridBmp& r)
{
  if (this != &r)
  {
    AOF::operator=(r);
    _nsamplex = r._nsamplex;
    _nsampley = r._nsampley;
    _nmult = r._nmult;
    _ncolor = r._ncolor;
    _flag_low = r._flag_low;
    _flag_high = r._flag_high;
    _mask_red = r._mask_red;
    _mask_green = r._mask_green;
    _mask_blue = r._mask_blue;
    _ffff_red = r._ffff_red;
    _ffff_green = r._ffff_green;
    _ffff_blue = r._ffff_blue;
    _low_red = r._low_red;
    _low_green = r._low_green;
    _low_blue = r._low_blue;
    _high_red = r._high_red;
    _high_green = r._high_green;
    _high_blue = r._high_blue;
    _valmin = r._valmin;
    _valmax = r._valmax;
    _reds = r._reds;
    _greens = r._greens;
    _blues = r._blues;
  }
  return *this;
}

GridBmp::~GridBmp()
{
}

void GridBmp::setColors(const VectorInt& reds, const VectorInt& greens, const VectorInt& blues)
{
  _reds = reds;
  _greens = greens;
  _blues = blues;
}

void GridBmp::setFFFF(int red, int green, int blue)
{
  _ffff_red = red;
  _ffff_green = green;
  _ffff_blue = blue;
}

void GridBmp::setHigh(int red, int green, int blue)
{
  _high_red = red;
  _high_green = green;
  _high_blue = blue;
}
void GridBmp::setLow(int red, int green, int blue)
{
  _low_red = red;
  _low_green = green;
  _low_blue = blue;
}
void GridBmp::setMask(int red, int green, int blue)
{
  _mask_red = red;
  _mask_green = green;
  _mask_blue = blue;
}

int GridBmp::writeInFile()
{
  VectorInt indg(2);
  unsigned char ired, igreen, iblue;

  /* Open the file */

  if (_fileWriteOpen()) return 1;

  /* Preliminary checks */

  int ncolor = _ncolor;
  bool flag_color_scale = (_ncolor > 0 && ! _reds.empty() && ! _greens.empty()
      && ! _blues.empty());
  if (! flag_color_scale) ncolor = 256;

  /* Initializations */

  int nx = _dbgrid->getNX(0);
  int ny = _dbgrid->getNX(1);

  /* Calculate the statistics */

  double vmin =  1.e30;
  double vmax = -1.e30;
  for (int i = 0; i < nx * ny; i++)
  {
    if (! _dbgrid->isActive(i)) continue;
    double value = _dbgrid->getArray(i, _cols[0]);
    if (FFFF(value)) continue;
    if (value < vmin) vmin = value;
    if (value > vmax) vmax = value;
  }
  double delta = (vmax - vmin);
  vmin -= delta / 100.;
  vmax += delta / 100.;
  if (! FFFF(_valmin)) vmin = _valmin;
  if (! FFFF(_valmax)) vmax = _valmax;

  /* Figure out the constants */

  int infosize = 40;
  int headersize = 14;
  int width  = _nmult * N_SAMPLE(nx, _nsamplex);
  int height = _nmult * N_SAMPLE(ny, _nsampley);
  int imagesize = 3 * width * height;

  /* Write the file header, bitmap information, and bitmap pixel data... */
  _writeOut( 0, BF_TYPE);
  _writeOut( 1, headersize); /* Size of File Header */
  _writeOut( 0, 0); /* Reserved */
  _writeOut( 0, 0); /* Reserved */
  _writeOut( 1, headersize + infosize); /* Offset */

  _writeOut( 1, infosize); /* Size of Information Block */
  _writeOut( 2, width); /* Width */
  _writeOut( 2, height); /* Height */
  _writeOut( 0, 1); /* Number of planes */
  _writeOut( 0, 24); /* 24-bits per pixel */
  _writeOut( 1, 0); /* No compression */
  _writeOut( 1, imagesize); /* Image size */
  _writeOut( 2, 0);
  _writeOut( 2, 0);
  _writeOut( 1, 0);
  _writeOut( 1, 0);

  /* Writing the pixels */

  indg[0] = 0;
  indg[1] = 0;
  int ipad = nx * _nmult;
  ipad = ipad - 4 * ((int) (ipad / 4));
  for (int iy = 0; iy < ny; iy++)
  {
    if (iy % _nsampley != 0) continue;
    for (int jmult = 0; jmult < _nmult; jmult++)
    {
      for (int ix = 0; ix < nx; ix++)
      {
        if (ix % _nsamplex != 0) continue;
        indg[0] = ix;
        indg[1] = iy;
        int iech = _dbgrid->indiceToRank(indg);
        int rank = _colorRank(iech, ncolor, vmin, vmax);
        _colorInRGB(rank, flag_color_scale, &ired, &igreen, &iblue);

        for (int imult = 0; imult < _nmult; imult++)
        {
          (void) fwrite(&iblue, 1, 1, _file);
          (void) fwrite(&igreen, 1, 1, _file);
          (void) fwrite(&ired, 1, 1, _file);
        }
      }

      /* Write the padding */

      int color = 0;
      for (int i = 0; i < ipad; i++)
        (void) fwrite(&color, 1, 1, _file);
    }
  }

  _fileClose();
  return 0;
}

/****************************************************************************/
/*!
 **   Print an integer
 **
 ** \param[in]  mode  Type of writing
 ** \li                0 : 16-bit unsigned integer
 ** \li                1 : 32-bit unsigned integer
 ** \li                2 : 32-bit signed integer
 ** \param[in]  ival  Integer value to be written
 **
 *****************************************************************************/
void GridBmp::_writeOut(int mode, unsigned int ival)
{
  switch (mode)
  {
    case 0: /* Unsigned 16-bit */
      putc(ival, _file);
      putc(ival >> 8, _file);
      break;

    case 1: /* Unsigned 32-bit */
      putc(ival, _file);
      putc(ival >> 8, _file);
      putc(ival >> 16, _file);
      putc(ival >> 24, _file);
      break;

    case 2: /* Signed 32-bit */
      int jval = (int) ival;
      putc(jval, _file);
      putc(jval >> 8, _file);
      putc(jval >> 16, _file);
      putc(jval >> 24, _file);
      break;
  }
}

/****************************************************************************/
/*!
 **   Color rank for the sample of the Db
 **
 ** \return  Color rank
 ** \return  COLOR_MASK  : The sample is masked off
 ** \return  COLOR_FFFF  : The value if FFFF
 ** \return  COLOR_LOWER : The value if below vmin and must be acknowledged
 ** \return  COLOR_UPPER : The value if above vmax and must be acknowledged
 **
 ** \param[in]  iech       Rank of the sample
 ** \param[in]  ncolor     Number of colors
 ** \param[in]  vmin       Minimum value to be represented
 ** \param[in]  vmax       Maximum value to be represented
 **
 *****************************************************************************/
int GridBmp::_colorRank(int iech, int ncolor, double vmin, double vmax)
{
  /* Check if the sample is masked off */
  if (! _dbgrid->isActive(iech)) return COLOR_MASK;

  /* Read the value */
  double value = _dbgrid->getArray(iech, _cols[0]);

  /* Check if the value is defined */
  if (FFFF(value)) return COLOR_FFFF;

  /* Find the color */
  int ival = (int) (ncolor * (value - vmin) / (vmax - vmin));

  /* Value lower than vmin */
  if (ival < 0)
  {
    if (_flag_low)
      return COLOR_LOWER;
    else
      return (ival);
  }

  /* Value larger than vmax */
  if (ival >= ncolor)
  {
    if (_flag_high)
      return COLOR_UPPER;
    else
      return (ncolor - 1);
  }

  /* Return the rank of the color */
  return (ival);
}

/****************************************************************************/
/*!
 ** Convert a color rank into the Red, Green, Blue color decomposition
 **
 ** \param[in]  rank        Rank of the color
 ** \param[in]  flag_color_scale 1 if the color scale must be used
 **                              0 use the grey scale instead
 **
 ** \param[out] ired        Value for the red beam
 ** \param[out] igreen      Value for the green beam
 ** \param[out] iblue       Value for the blue beam
 **
 *****************************************************************************/
void GridBmp::_colorInRGB(int rank,
                          bool flag_color_scale,
                          unsigned char *ired,
                          unsigned char *igreen,
                          unsigned char *iblue)
{
  switch (rank)
  {
    case COLOR_MASK:
      *ired = (unsigned char) _mask_red;
      *igreen = (unsigned char) _mask_green;
      *iblue = (unsigned char) _mask_blue;
      break;

    case COLOR_FFFF:
      *ired = (unsigned char) _ffff_red;
      *igreen = (unsigned char) _ffff_green;
      *iblue = (unsigned char) _ffff_blue;
      break;

    case COLOR_LOWER:
      *ired = (unsigned char) _low_red;
      *igreen = (unsigned char) _low_green;
      *iblue = (unsigned char) _low_blue;
      break;

    case COLOR_UPPER:
      *ired = (unsigned char) _high_red;
      *igreen = (unsigned char) _high_green;
      *iblue = (unsigned char) _high_blue;
      break;

    default:
      if (flag_color_scale)
      {
        *ired = (unsigned char) _reds[rank];
        *igreen = (unsigned char) _greens[rank];
        *iblue = (unsigned char) _blues[rank];
      }
      else
      {
        *ired = (unsigned char) rank;
        *igreen = (unsigned char) rank;
        *iblue = (unsigned char) rank;
      }
  }
}

DbGrid*  GridBmp::readGridFromFile()
{
  DbGrid* dbgrid = nullptr;
  int ir[256], ig[256], ib[256];
  VectorDouble x0(2);
  VectorDouble dx(2);
  VectorInt nx(2);

  /* Open the file */

  if (_fileReadOpen()) return dbgrid;

  /* Initializations */

  x0[0] = 0.;
  x0[1] = 0.;
  dx[0] = 1.;
  dx[1] = 1.;

  /* Reading the file header */

  (void) _compose( 2);
  (void) _compose( 4);
  (void) _compose( 4);
  (void) _compose( 4);

  /* Reading the bitmap information header */

  (void) _compose( 4);
  nx[0] = _compose( 4);
  nx[1] = _compose( 4);
  (void) _compose( 2);
  int nbits = _compose( 2);
  int compress = _compose( 4);
  (void) _compose( 4);
  int ndx = _compose( 4);
  int ndy = _compose( 4);
  int ncol = _compose( 4);
  (void) _compose( 4);
  if (ncol > 256)
  {
    messerr("Your file seems to contain more than 256 colors");
    messerr("It is probably not a valid BMP file");
    return dbgrid;
  }

  /* Read the Color Scale (optional) */

  for (int icol = 0; icol < ncol; icol++)
  {
    ir[icol] = _readIn();
    ig[icol] = _readIn();
    ib[icol] = _readIn();
    (void) _readIn();
  }

  /* Final checks */

  if (compress != 0)
  {
    messerr("Error : only uncompressed images are treated here");
    return dbgrid;
  }

  /* Final results */

  if (ndx > 0) dx[0] = 100. / (double) ndx;
  if (ndy > 0) dx[1] = 100. / (double) ndy;

  /* Reading the image (from bottom to up) */

  int noct = nx[0] * nbits / 8;
  int npad = 4 - noct % 4;
  if (npad == 4) npad = 0;
  int size = nx[0] * nx[1];
  VectorDouble tab(size);

  int ecr = 0;
  for (int jy = 0; jy < nx[1]; jy++)
  {
    for (int ix = 0; ix < nx[0]; ix++)
    {
      unsigned char ctab = 0;
      if (nbits == 8)
      {
        unsigned char c = _readIn();
        _rgb2num(ir[c], ig[c], ib[c], 0, &ctab);
      }
      else if (nbits == 24)
        _rgb2num(_readIn(), _readIn(), _readIn(), 0, &ctab);
      else if (nbits == 32)
        _rgb2num(_readIn(), _readIn(), _readIn(), _readIn(), &ctab);
      else
      {
        messerr("Error: Number of bits/pixel (%d) is not treated", nbits);
        return dbgrid;
      }
      tab[ecr++] = ctab;
    }

    /* Reading the padding */

    for (int ix = 0; ix < npad; ix++) (void) _readIn();
  }

   dbgrid = new DbGrid();
   dbgrid->reset(nx,dx,x0,VectorDouble(),ELoadBy::SAMPLE,tab);

  // Close the file

  _fileClose();

  return dbgrid;
}

/****************************************************************************/
/*!
 **   Compose several bytes into an integer
 **
 ** \return  Returned integer
 **
 ** \param[in]  nb    Number of bytes to be considered
 **
 *****************************************************************************/
int GridBmp::_compose(int nb)

{
  int i, value, factor;
  unsigned char c;

  /* Initializations */

  value = 0;
  factor = 1;

  for (i = 0; i < nb; i++)
  {
    c = _readIn();
    value += c * factor;
    factor *= 0x100;
  }
  return (value);
}

/****************************************************************************/
/*!
 ** Read a byte from the binary file
 **
 *****************************************************************************/
unsigned char GridBmp::_readIn()

{
  unsigned char c;
  c = (unsigned char) fgetc(_file);
  if (!feof(_file)) return c;
  return (c);
}

/****************************************************************************/
/*!
 **  Convert numeric to RGB
 **
 ** \param[in]   value   Input value
 **
 ** \param[out]  r     Red index
 ** \param[out]  g     Green index
 ** \param[out]  b     Blue index
 ** \param[out]  a     Transparency index
 **
 *****************************************************************************/
void GridBmp::_num2rgb(unsigned char value, int *r, int *g, int *b, int *a)
{
  *r = static_cast<int>((value >> 24) & 0xff);
  *g = static_cast<int>((value >> 16) & 0xff);
  *b = static_cast<int>((value >> 8) & 0xff);
  *a = static_cast<int>((value) & 0xff);
}

/****************************************************************************/
/*!
 **   Convert RGB into numeric
 **
 ** \param[in]  red     Red index
 ** \param[in]  green   Green index
 ** \param[in]  blue    Blue index
 **
 ** \param[out] c       Numeric value
 **
 *****************************************************************************/
void GridBmp::_rgb2num(int red, int green, int blue, int /*a*/, unsigned char *c)
{
  double value;

  value = (double) (red + green + blue) / 3.;

  if (value < 0.) value = 0.;
  if (value > 255.) value = 255.;
  *c = (unsigned char) value;

  return;
}

