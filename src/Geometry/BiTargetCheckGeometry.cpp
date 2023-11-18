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
#include "geoslib_f.h"

#include "Geometry/BiTargetCheckGeometry.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "Space/SpaceTarget.hpp"

BiTargetCheckGeometry::BiTargetCheckGeometry(int ndim,
                                             const VectorDouble &codir,
                                             double tolang,
                                             double bench,
                                             double cylrad,
                                             bool flagasym)
    : ABiTargetCheck(),
      _ndim(ndim),
      _codir(codir),
      _tolAng(tolang),
      _bench(bench),
      _cylrad(cylrad),
      _flagAsym(flagasym),
      _psmin(0.),
      _dist(0.)
{
  _psmin = GH::getCosineAngularTolerance(tolang);
}

BiTargetCheckGeometry::BiTargetCheckGeometry(const BiTargetCheckGeometry &r)
    : ABiTargetCheck(r),
      _ndim(r._ndim),
      _codir(r._codir),
      _tolAng(r._tolAng),
      _bench(r._bench),
      _cylrad(r._cylrad),
      _flagAsym(r._flagAsym),
      _psmin(r._psmin),
      _dist(r._dist)
{
}

BiTargetCheckGeometry& BiTargetCheckGeometry::operator=(const BiTargetCheckGeometry &r)
{
  if (this != &r)
  {
    ABiTargetCheck::operator=(r);
    _ndim = r._ndim;
    _codir = r._codir;
    _tolAng = r._tolAng;
    _bench = r._bench;
    _cylrad = r._cylrad;
    _flagAsym = r._flagAsym;
    _psmin = r._psmin;
    _dist = r._dist;
  }
  return *this;
}

BiTargetCheckGeometry::~BiTargetCheckGeometry()
{
}

BiTargetCheckGeometry* BiTargetCheckGeometry::create(int ndim,
                                                     const VectorDouble &codir,
                                                     double tolang,
                                                     double bench,
                                                     double cylrad,
                                                     bool flagasym)
{
  return new BiTargetCheckGeometry(ndim, codir, tolang, bench, cylrad, flagasym);
}

String BiTargetCheckGeometry::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "- Direction" << VH::toStringAsVD(_codir) << std::endl;
  sstr << "- Tolerance angular" << _tolAng << std::endl;
  if (!FFFF(_bench) && _bench > 0.)
    sstr << "Bench (%lf)" << _bench << std::endl;
  if (!FFFF(_cylrad) && _cylrad > 0.)
    sstr << "Cylinder check (%lf)" << _cylrad << std::endl;

  return sstr.str();
}

bool BiTargetCheckGeometry::isOK(const SpaceTarget &T1, const SpaceTarget &T2) const
{
  // Calculate the distance between the two samples
  _dist = T1.getDistance(T2);

  // When the distance is zero, the pair is always accepted
  if (_dist <= 0.) return true;

  // Increment between two samples
  VectorDouble delta = T1.getIncrement(T2);

  // Check if the angle of the pair matches the Calculation direction (up to angular tolerance)
  double dproj = 0.;
  double dn1 = 0.;
  double dn2 = 0.;
  for (int idim = 0; idim < _ndim; idim++)
  {
    dproj +=  delta[idim] * _codir[idim];
    dn1   +=  delta[idim] *  delta[idim];
    dn2   += _codir[idim] * _codir[idim];
  }
  double prod = dn1 * dn2;
  double ps = 1.;
  if (prod > 0.) ps = dproj / sqrt(prod);
  if (ABS(ps) < _psmin) return false;

  // Check for cylinder test
  if (!FFFF(_cylrad) && _cylrad > 0.)
  {
    double dortho = 0.;
    if (prod > 0.) dortho = sqrt(dn1 * (1. - ps * ps));
    if (dortho > _cylrad) return false;
  }

  // Check for vertical slicing test
  if (!FFFF(_bench) && _bench > 0.)
  {
    double dvect = ABS(delta[_ndim-1]);
    if (dvect > _bench) return false;
  }

  // Calculate the oriented distance
  if (_flagAsym && ps < _psmin) _dist = -_dist;

  return true;
}
