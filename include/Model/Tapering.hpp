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

GSTLEARN_EXPORT Def_Tapering& D_TAPE(int rank);

class GSTLEARN_EXPORT Tapering : public AStringable
{
public:
  Tapering();
  Tapering(const Tapering &m);
  Tapering& operator= (const Tapering &m);
  virtual ~Tapering();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static int getTapeNumber();

  double getRange() const         { return _range; }
  int getType() const             { return _type; }
  int getMaxNDim() const          { return _maxNDim; }
  void setRange(double taperange) { _range = taperange; }
  void setType(int tapetype)      { _type = tapetype; }

  double evaluate(double h) const { return D_TAPE(_type).tapeFunc(h); }

  int init(int tape_type,double tape_range);

private:
  int    _type;
  int    _maxNDim;
  double _range;
  String _name;
};
