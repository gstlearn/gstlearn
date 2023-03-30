/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "OutputFormat/AOF.hpp"

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
  bool mustBeForNDim(int ndim) const override { return ndim == 2; }
  bool mustBeForRotation(int mode) const override { return mode == 0; }
  int  writeInFile() override;
  DbGrid* readGridFromFile() override;

  void setColors(const VectorInt& reds, const VectorInt& greens, const VectorInt& blues);
  void setFFFF(int red, int green, int blue);
  void setHigh(int red, int green, int blue);
  void setLow(int red, int green, int blue);
  void setMask(int red, int green, int blue);
  void setFlagHigh(bool flagHigh) { _flag_high = flagHigh; }
  void setFlagLow(bool flagLow)   { _flag_low = flagLow; }

  void setNcolor(int ncolor) { _ncolor = ncolor; }
  void setNmult(int nmult) { _nmult = nmult; }
  void setNsamplex(int nsamplex) { _nsamplex = nsamplex; }
  void setNsampley(int nsampley) { _nsampley = nsampley; }
  void setValmax(double valmax) { _valmax = valmax; }
  void setValmin(double valmin) { _valmin = valmin; }

private:
  void _writeOut(int mode, unsigned int ival);
  int  _colorRank(int iech, int ncolor, double vmin, double vmax);
  void _colorInRGB(int rank,
                   bool flag_color_scale,
                   unsigned char *ired,
                   unsigned char *igreen,
                   unsigned char *iblue);
  int _compose(int nb);
  unsigned char _readIn();
  //void _num2rgb(unsigned char value, int *r, int *g, int *b, int *a);
  void _rgb2num(int red, int green, int blue, int a, unsigned char *c);

private:
  int _nsamplex;
  int _nsampley;
  int _nmult;
  int _ncolor;
  bool _flag_low;
  bool _flag_high;
  int _mask_red;
  int _mask_green;
  int _mask_blue;
  int _ffff_red;
  int _ffff_green;
  int _ffff_blue;
  int _low_red;
  int _low_green;
  int _low_blue;
  int _high_red;
  int _high_green;
  int _high_blue;
  double _valmin;
  double _valmax;
  VectorInt _reds;
  VectorInt _greens;
  VectorInt _blues;
};
