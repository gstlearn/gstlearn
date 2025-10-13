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

#include "Geometry/ABiTargetCheck.hpp"

namespace gstlrn
{
class GSTLEARN_EXPORT BiTargetCheckGeometry: public ABiTargetCheck
{
public:
  BiTargetCheckGeometry(Id ndim,
                        const VectorDouble &codir = VectorDouble(),
                        double tolang = 90.,
                        double bench = 0.,
                        double cylrad = 0.,
                        bool flagasym = false);
  BiTargetCheckGeometry(const BiTargetCheckGeometry& r);
  BiTargetCheckGeometry& operator=(const BiTargetCheckGeometry& r);
  virtual ~BiTargetCheckGeometry();

  /// ICloneable Interface
  IMPLEMENT_CLONING(BiTargetCheckGeometry)

  bool isOK(const SpaceTarget &T1, const SpaceTarget &T2) const override;

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  static BiTargetCheckGeometry* create(Id ndim,
                                       const VectorDouble &codir = VectorDouble(),
                                       double tolang = 90.,
                                       double bench = 0.,
                                       double cylrad = 0.,
                                       bool flagasym = false);

  double getDist() const { return _dist; }

private:
  Id _ndim;
  VectorDouble _codir;
  double _tolAng;
  double _bench;
  double _cylrad;
  bool   _flagAsym;

  mutable VectorDouble _delta;
  mutable double _psmin;
  mutable double _dist;
};
}