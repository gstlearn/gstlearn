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

class GSTLEARN_EXPORT ProjComposition : public IProj
{
public:
  ProjComposition(std::vector<const IProj*> projs);
  ProjComposition(const ProjComposition&) = delete;
  ProjComposition& operator=(const ProjComposition&) = delete;
  ~ProjComposition() override = default;

#ifndef SWIG
protected:
  int _addPoint2mesh(const constvect in, vect out) const override;
  int _addMesh2point(const constvect in, vect out) const override;
#endif

public:
  int getNApex() const override { return (_projs.size() == 0 ? 0 : _projs.front()->getNApex()); }
  int getNPoint() const override { return (_projs.size() == 0 ? 0 : _projs.back()->getNPoint()); }

private:
  std::vector<std::unique_ptr<const IProj>> _projs;
  mutable std::vector<std::vector<double>> _works;
};
