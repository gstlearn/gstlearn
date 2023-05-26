/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Geometry/ABiPointCheck.hpp"

class GSTLEARN_EXPORT BiPointCheckDistance: public ABiPointCheck
{
public:
  BiPointCheckDistance(double radius = TEST,
                       const VectorDouble coeffs = VectorDouble(),
                       const VectorDouble angles = VectorDouble());
  BiPointCheckDistance(const BiPointCheckDistance& r);
  BiPointCheckDistance& operator=(const BiPointCheckDistance& r);
  virtual ~BiPointCheckDistance();

  /// ICloneable Interface
  IMPLEMENT_CLONING(BiPointCheckDistance)

  virtual bool isOK(const SpacePoint &P1,
                    const SpacePoint &P2,
                    int iech1 = -1,
                    int iech2 = -1) const override;

  static BiPointCheckDistance* create(double radius = TEST,
                                      const VectorDouble coeffs = VectorDouble(),
                                      const VectorDouble angles = VectorDouble());

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int getNDim() const { return _ndim; }
  const VectorDouble& getAnisoCoeffs() const { return _anisoCoeffs; }
  double getAnisoCoeff(int i) const { return _anisoCoeffs[i]; }
  const VectorDouble& getAnisoRotMats() const { return _anisoRotMat; }
  double getAnisoRotMat(int i) const { return _anisoRotMat[i]; }
  int getFlagAniso() const { return _flagAniso; }
  int getFlagRotation() const { return _flagRotation; }
  double getRadius() const { return _radius; }

  void setAnisoCoeffs(const VectorDouble& anisoCoeffs) { _anisoCoeffs = anisoCoeffs; }
  void setAnisoCoeff(int idim, double value);
  void setAnisoRotMat(const VectorDouble& anisoRotMat) { _anisoRotMat = anisoRotMat; }
  void setFlagAniso(int flagAniso) { _flagAniso = flagAniso; }
  void setFlagRotation(int flagRotation) { _flagRotation = flagRotation; }
  void setRadius(double radius) { _radius = radius; }

  double getDistance() const { return _dist; }
  VectorDouble getIncr() const { return _movingIncr; }
  double getNormalizedDistance(const VectorDouble& dd) const;

private:
  void   _calculateDistance() const;

private:
  int    _ndim;                  /* Space dimension (used for array dimensioning) */
  bool   _flagAniso;             /* 1 if the MOVING neigh. is anisotropic */
  bool   _flagRotation;          /* 1 if the anisotropy is rotated */
  double _radius;                /* Maximum isotropic distance */
  VectorDouble _anisoCoeffs;     /* Anisotropy ratio for moving neighborhood */
  VectorDouble _anisoRotMat;     /* Anisotropy rotation matrix */

  mutable double _dist;
  mutable VectorDouble _movingIncr;
  mutable VectorDouble _movingAux;
};
