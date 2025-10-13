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

namespace gstlrn
{

class GSTLEARN_EXPORT ProjZero : public IProj
{
public:
  ProjZero(Id npoint, Id napex) : IProj(), _npoint(npoint), _napex(napex) {}
  ProjZero(const ProjZero&) = default;
  ProjZero& operator=(const ProjZero&) = default;
  ~ProjZero() override = default;

#ifndef SWIG
protected:
  Id _addPoint2mesh(const constvect inv, vect outv) const override { DECLARE_UNUSED(inv, outv); return 0; }
  Id _addMesh2point(const constvect inv, vect outv) const override { DECLARE_UNUSED(inv, outv); return 0; }
#endif
public:
  Id getNApex() const override { return _napex; }
  Id getNPoint() const override { return _npoint; }

private:
  Id _npoint, _napex;
};
}