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

class GSTLEARN_EXPORT BiTargetCheckGeometry: public ABiTargetCheck
{
public:
  BiTargetCheckGeometry(int ndim,
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

  virtual bool isOK(const SpaceTarget &T1, const SpaceTarget &T2) const override;

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static BiTargetCheckGeometry* create(int ndim,
                                       const VectorDouble &codir = VectorDouble(),
                                       double tolang = 90.,
                                       double bench = 0.,
                                       double cylrad = 0.,
                                       bool flagasym = false);

  double getDist() const { return _dist; }

private:
  int _ndim;
  VectorDouble _codir;
  double _tolAng;
  double _bench;
  double _cylrad;
  bool   _flagAsym;

  mutable double _psmin;
  mutable double _dist;
};
