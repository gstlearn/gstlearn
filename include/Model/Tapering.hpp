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

#include "Basic/AStringable.hpp"

namespace gstlrn
{
typedef struct
{
  std::string tapeName;
  Id tapeType; /* Rank in the covariance list */
  Id maxNDim;  /* Maximum dimension for validity */
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

GSTLEARN_EXPORT Def_Tapering& D_TAPE(Id rank);

class GSTLEARN_EXPORT Tapering : public AStringable
{
public:
  Tapering();
  Tapering(const Tapering &m);
  Tapering& operator= (const Tapering &m);
  virtual ~Tapering();

  String toString(const AStringFormat* strfmt = nullptr) const override;

  static Id getNTape();

  double getRange() const         { return _range; }
  Id getType() const             { return _type; }
  Id getMaxNDim() const          { return _maxNDim; }
  void setRange(double taperange) { _range = taperange; }
  void setType(Id tapetype)      { _type = tapetype; }

  double evaluate(double h) const { return D_TAPE(_type).tapeFunc(h); }

  Id init(Id tape_type,double tape_range);

private:
  Id    _type;
  Id    _maxNDim;
  double _range;
  String _name;
};
}