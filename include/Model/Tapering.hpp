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

#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"

typedef struct
{
  std::string tapeName;
  int tapeType; /* Rank in the covariance list */
  int maxNDim;  /* Maximum dimension for validity */
  double (*tapeFunc)(double);
} Def_Tapering;

/* Prototyping the internal covariance functions */
double _tape_spherical(double);
double _tape_cubic(double);
double _tape_triangle(double);
double _tape_penta(double);
double _tape_storkey(double);
double _tape_wendland1(double);
double _tape_wendland2(double);

Def_Tapering& D_TAPE(int rank);

class Tapering : public AStringable
{
public:
  Tapering();
  Tapering(const Tapering &m);
  Tapering& operator= (const Tapering &m);
  virtual ~Tapering();

  virtual std::string toString(int level = 0) const override;

  int getTapeNumber();

  double getRange() const         { return _range; }
  int getType() const             { return _type; }
  int getMaxNDim() const          { return _maxNDim; }
  void setRange(double taperange) { _range = taperange; }
  void setType(int tapetype)      { _type = tapetype; }

  double evaluate(double h) const
  {
    return D_TAPE(_type).tapeFunc(h);
  }

  int init(int tape_type,double tape_range);

private:
  int    _type;
  int    _maxNDim;
  double _range;
  String _name;
};
