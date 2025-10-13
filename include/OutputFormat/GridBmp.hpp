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
#pragma once

#include "gstlearn_export.hpp"
#include "OutputFormat/AOF.hpp"

namespace gstlrn
{
class Db;

class GSTLEARN_EXPORT GridBmp: public AOF
{
public:
  GridBmp(const char* filename, const Db* db = nullptr);
  GridBmp(const GridBmp& r);
  GridBmp& operator=(const GridBmp& r);
  virtual ~GridBmp();

  bool mustBeGrid() const override { return true; }
  bool mustBeOneVariable() const override { return true; }
  bool mustBeForNDim(Id ndim) const override { return ndim == 2; }
  bool mustBeForRotation(Id mode) const override { return mode == 0; }
  Id  writeInFile() override;
  DbGrid* readGridFromFile() override;

  void setColors(const VectorInt& reds, const VectorInt& greens, const VectorInt& blues);
  void setFFFF(Id red, Id green, Id blue);
  void setHigh(Id red, Id green, Id blue);
  void setLow(Id red, Id green, Id blue);
  void setMask(Id red, Id green, Id blue);
  void setFlagHigh(bool flagHigh) { _flag_high = flagHigh; }
  void setFlagLow(bool flagLow)   { _flag_low = flagLow; }

  void setNcolor(Id ncolor) { _ncolor = ncolor; }
  void setNmult(Id nmult) { _nmult = nmult; }
  void setNsamplex(Id nsamplex) { _nsamplex = nsamplex; }
  void setNsampley(Id nsampley) { _nsampley = nsampley; }
  void setValmax(double valmax) { _valmax = valmax; }
  void setValmin(double valmin) { _valmin = valmin; }

private:
  void _writeOut(Id mode, I32 ival);
  Id  _colorRank(Id iech, Id ncolor, double vmin, double vmax);
  void _colorInRGB(Id rank,
                   bool flag_color_scale,
                   unsigned char *ired,
                   unsigned char *igreen,
                   unsigned char *iblue);
  Id _compose(Id nb);
  unsigned char _readIn();
  //void _num2rgb(unsigned char value, Id *r, Id *g, Id *b, Id *a);
  static void _rgb2num(Id red, Id green, Id blue, Id a, unsigned char *c);

private:
  Id _nsamplex;
  Id _nsampley;
  Id _nmult;
  Id _ncolor;
  bool _flag_low;
  bool _flag_high;
  Id _mask_red;
  Id _mask_green;
  Id _mask_blue;
  Id _ffff_red;
  Id _ffff_green;
  Id _ffff_blue;
  Id _low_red;
  Id _low_green;
  Id _low_blue;
  Id _high_red;
  Id _high_green;
  Id _high_blue;
  double _valmin;
  double _valmax;
  VectorInt _reds;
  VectorInt _greens;
  VectorInt _blues;
};
} // namespace gstlrn
