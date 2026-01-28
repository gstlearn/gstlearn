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

#include "Basic/AStringable.hpp"
#include "LinearOp/IProj.hpp"
#include "gstlearn_export.hpp"

#include <functional>
#include <optional>
#include <vector>

namespace gstlrn
{
class GSTLEARN_EXPORT ProjMulti: public IProj, public AStringable
{
public:
  using ProjVect = std::vector<std::vector<std::optional<std::reference_wrapper<const IProj>>>>;

  ProjMulti(const std::vector<std::vector<const IProj*>>& projs, bool silent = false);

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  Id getNApex() const override;
  Id getNPoint() const override;
  Id getNVariable() const { return _nvariable; }
  Id getNLatent() const { return _nlatent; }
  ~ProjMulti() override = default;
  bool empty() const { return _projs.empty(); }

#ifndef SWIG

protected:
  Id _addPoint2mesh(const constvect inv, vect outv) const override;
  Id _addMesh2point(const constvect inv, vect outv) const override;
#endif

private:
  bool _checkArg(const ProjVect& projs) const;
  void _init();

protected:
  Id findFirstNoNullOnRow(Id j) const;
  Id findFirstNoNullOnCol(Id j) const;
  const VectorInt& getNPoints() const { return _pointNumbers; }
  const VectorInt& getNApexs() const { return _apexNumbers; }

protected:
  ProjVect _projs;

private:
  Id _pointNumber;
  Id _apexNumber;
  Id _nlatent;
  Id _nvariable;
  VectorInt _pointNumbers;
  VectorInt _apexNumbers;
  bool _silent;
  mutable VectorDouble _work;
  mutable VectorDouble _workmesh;
};
} // namespace gstlrn
