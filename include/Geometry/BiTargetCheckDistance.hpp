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
class GSTLEARN_EXPORT BiTargetCheckDistance: public ABiTargetCheck
{
public:
  BiTargetCheckDistance(double radius = TEST,
                       const VectorDouble& coeffs = VectorDouble(),
                       const VectorDouble& angles = VectorDouble());
  BiTargetCheckDistance(const BiTargetCheckDistance& r);
  BiTargetCheckDistance& operator=(const BiTargetCheckDistance& r);
  virtual ~BiTargetCheckDistance();

  /// ICloneable Interface
  IMPLEMENT_CLONING(BiTargetCheckDistance)

  bool isOK(const SpaceTarget &T1, const SpaceTarget &T2) const override;

  static BiTargetCheckDistance*
  create(double radius              = TEST,
         const VectorDouble& coeffs = VectorDouble(),
         const VectorDouble& angles = VectorDouble());

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  Id getNDim() const { return _ndim; }
  const VectorDouble& getAnisoCoeffs() const { return _anisoCoeffs; }
  double getAnisoCoeff(Id i) const { return _anisoCoeffs[i]; }
  const VectorDouble& getAnisoRotMats() const { return _anisoRotMat; }
  double getAnisoRotMat(Id i) const { return _anisoRotMat[i]; }
  Id getFlagAniso() const { return _flagAniso; }
  Id getFlagRotation() const { return _flagRotation; }
  double getRadius() const { return _radius; }

  void setNDim(Id ndim) { _ndim = ndim; }
  void setAnisoCoeffs(const VectorDouble& anisoCoeffs) { _anisoCoeffs = anisoCoeffs; }
  void setAnisoRotMat(const VectorDouble& anisoRotMat) { _anisoRotMat = anisoRotMat; }
  void setFlagAniso(Id flagAniso) { _flagAniso = flagAniso; }
  void setFlagRotation(Id flagRotation) { _flagRotation = flagRotation; }
  void setRadius(double radius) { _radius = radius; }

  double getDistance() const { return _dist; }
  const VectorDouble& getIncr() const { return _movingIncr; }
  double getNormalizedDistance(const VectorDouble& dd) const;

private:
  void _calculateDistance() const;

private:
  Id    _ndim;                  /* Space dimension (used for array dimensioning) */
  bool   _flagAniso;             /* 1 if the MOVING neigh. is anisotropic */
  bool   _flagRotation;          /* 1 if the anisotropy is rotated */
  double _radius;                /* Maximum isotropic distance */
  VectorDouble _anisoCoeffs;     /* Anisotropy ratio for moving neighborhood */
  VectorDouble _anisoRotMat;     /* Anisotropy rotation matrix */

  mutable double _dist;
  mutable VectorDouble _movingIncr;
  mutable VectorDouble _movingAux;
};
}