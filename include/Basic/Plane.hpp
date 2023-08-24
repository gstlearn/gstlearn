/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include"Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"

class DbGrid;

class GSTLEARN_EXPORT Plane : public AStringable
{
public:
  Plane();
  Plane(const Plane &m);
  Plane& operator=(const Plane &m);
  virtual ~Plane();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  const VectorDouble& getCoor() const { return _coor; }
  void setCoor(const VectorDouble& coor) { _coor = coor; }
  double getIntercept() const { return _intercept; }
  void setIntercept(double intercept) { _intercept = intercept; }
  double getRndval() const { return _rndval; }
  void setRndval(double rndval) { _rndval = rndval; }
  double getValue() const { return _value; }
  void setValue(double value) { _value = value; }
  void setCoor(int idim, double value);
  double getCoor(int idim) const;

  static std::vector<Plane> poissonPlanesGenerate(DbGrid *dbgrid, int np);

private:
  VectorDouble _coor;
  double _intercept;
  double _value;
  double _rndval;
};
