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

#include "LinearOp/IProj.hpp"
#include "gstlearn_export.hpp"

#include <functional>

namespace gstlrn
{
  class GSTLEARN_EXPORT ProjComposition: public IProj
  {
  public:
    ProjComposition(std::vector<const IProj*> projs);
    ProjComposition(const ProjComposition&) = delete;
    ProjComposition& operator=(const ProjComposition&) = delete;
    ~ProjComposition() override = default;

#ifndef SWIG

  protected:
    Id _addPoint2mesh(const constvect in, vect out) const override;
    Id _addMesh2point(const constvect in, vect out) const override;
#endif

  public:
    Id getNApex() const override
    {
      return (_projs.size() == 0 ? 0 : _projs.front().get().getNApex());
    }

    Id getNPoint() const override
    {
      return (_projs.size() == 0 ? 0 : _projs.back().get().getNPoint());
    }

  private:
    using ProjVect = std::vector<std::reference_wrapper<const IProj>>;

    ProjVect _projs;
    mutable VectorVectorDouble _works;
  };
} // namespace gstlrn
