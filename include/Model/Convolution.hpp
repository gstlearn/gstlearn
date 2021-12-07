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
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"

typedef struct
{
  std::string convName;
  double convScale;
  double (*convFunc)(int, double);
} Def_Conv;

/* Prototyping the internal covariance functions */
double _conv_uniform(int number, double v);
double _conv_exponential(int number, double v);
double _conv_gaussian(int number, double v);
double _conv_sincard(int number, double v);

GSTLEARN_EXPORT Def_Conv& D_CONV(int rank);

class GSTLEARN_EXPORT Convolution : public AStringable
{
public:
  Convolution();
  Convolution(const Convolution &m);
  Convolution& operator= (const Convolution &m);
  virtual ~Convolution();

  virtual String toString(int level = 0) const override;

  int getConvNumber();

  const VectorDouble& getDelX() const { return _delX; }
  const VectorDouble& getDelY() const { return _delY; }
  const VectorDouble& getDelZ() const { return _delZ; }
  int getDir() const         { return _dir; }
  int getDiscNumber() const  { return _discNumber; }
  int getCount() const       { return _count; }
  double getRange() const    { return _range; }
  int getType() const        { return _type; }
  const VectorDouble& getWeights() const { return _weight; }

  void setType(int type)          { _type = type; }
  void setDir(int dir)            { _dir = dir; }
  void setDiscNumber(int ndisc)   { _discNumber = ndisc; }
  void setCount(int count)        { _count = count; }
  void setRange(double range)     { _range = range; }

  double evaluateConv(int number, double v)
  {
    return D_CONV(_type).convFunc(number, v);
  }

  double (*_convFunc)(int, double);
  std::string _convName;

  int init(int conv_type, int conv_idir, int conv_ndisc, double conv_range);

private:
  int _type; /* Convolution type */
  int _dir; /* Convolution direction */
  int _discNumber; /* Number of discretization per direction */
  int _count; /* Number of discretization points */
  double _range; /* Convolution Range */
  VectorDouble _delX; /* Discretization along X */
  VectorDouble _delY; /* Discretization along Y */
  VectorDouble _delZ; /* Discretization along Z */
  VectorDouble _weight; /* Weights for convolution */
};
