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
#include "LinearOp/IProj.hpp"

class GSTLEARN_EXPORT ProjZero : public IProj
{
public:
  ProjZero(int npoint, int napex) : IProj(), _npoint(npoint), _napex(napex) {}
  ProjZero(const ProjZero&) = default;
  ProjZero& operator=(const ProjZero&) = default;
  ~ProjZero() override = default;

#ifndef SWIG
protected:
  int _addPoint2mesh(const constvect, vect) const override { return 0; }
  int _addMesh2point(const constvect, vect) const override { return 0; }
#endif
public:
  int getNApex() const override { return _napex; }
  int getNPoint() const override { return _npoint; }

private:
  int _npoint, _napex;
};
