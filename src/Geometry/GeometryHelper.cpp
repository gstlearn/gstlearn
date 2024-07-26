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
#include "Geometry/GeometryHelper.hpp"

#include "Enum/ERotation.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Law.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpaceSN.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"

#include <math.h>

void GeometryHelper::_decodeConvRot(const ERotation &convrot,
                                    int *firstaxis,
                                    int *parity,
                                    int *repetition,
                                    int *frame)
{
  VectorInt ret(4);

  if (convrot == ERotation::SXYZ) ret = { 0, 0, 0, 0 };
  if (convrot == ERotation::SXYX) ret = { 0, 0, 1, 0 };
  if (convrot == ERotation::SXZY) ret = { 0, 1, 0, 0 };
  if (convrot == ERotation::SXZX) ret = { 0, 1, 1, 0 };
  if (convrot == ERotation::SYZX) ret = { 1, 0, 0, 0 };
  if (convrot == ERotation::SYZY) ret = { 1, 0, 1, 0 };
  if (convrot == ERotation::SYXZ) ret = { 1, 1, 0, 0 };
  if (convrot == ERotation::SYXY) ret = { 1, 1, 1, 0 };
  if (convrot == ERotation::SZXY) ret = { 2, 0, 0, 0 };
  if (convrot == ERotation::SZXZ) ret = { 2, 0, 1, 0 };
  if (convrot == ERotation::SZYX) ret = { 2, 1, 0, 0 };
  if (convrot == ERotation::SZYZ) ret = { 2, 1, 1, 0 };
  if (convrot == ERotation::RZYX) ret = { 0, 0, 0, 1 };
  if (convrot == ERotation::RXYX) ret = { 0, 0, 1, 1 };
  if (convrot == ERotation::RYZX) ret = { 0, 1, 0, 1 };
  if (convrot == ERotation::RXZX) ret = { 0, 1, 1, 1 };
  if (convrot == ERotation::RXZY) ret = { 1, 0, 0, 1 };
  if (convrot == ERotation::RYZY) ret = { 1, 0, 1, 1 };
  if (convrot == ERotation::RZXY) ret = { 1, 1, 0, 1 };
  if (convrot == ERotation::RYXY) ret = { 1, 1, 1, 1 };
  if (convrot == ERotation::RYXZ) ret = { 2, 0, 0, 1 };
  if (convrot == ERotation::RZXZ) ret = { 2, 0, 1, 1 };
  if (convrot == ERotation::RXYZ) ret = { 2, 1, 0, 1 };
  if (convrot == ERotation::RZYZ) ret = { 2, 1, 1, 1 };

  *firstaxis = ret[0];
  *parity = ret[1];
  *repetition = ret[2];
  *frame = ret[3];
}

/****************************************************************************/
/*!
 **  Calculates the trigonometric features
 **
 ** \param[in]  angle input angle (in degrees)
 **
 ** \param[out]  cosa  cosine function
 ** \param[out]  sina  sine function
 **
 *****************************************************************************/
void GeometryHelper::rotationGetSinCos(double angle, double *cosa, double *sina)
{
  double value;

  if (angle == 0.)
  {
    *cosa = 1.;
    *sina = 0.;
    return;
  }
  if (angle == 90.)
  {
    *cosa = 0.;
    *sina = 1.;
    return;
  }
  if (angle == 180.)
  {
    *cosa = -1.;
    *sina = 0.;
    return;
  }
  if (angle == 270.)
  {
    *cosa = 0.;
    *sina = -1.;
    return;
  }
  value = ut_deg2rad(angle);
  *cosa = cos(value);
  *sina = sin(value);
}

/****************************************************************************/
/*!
 **  Initialize a rotation matrix
 **
 ** \param[in]  ndim      Space dimension
 **
 ** \param[out] rot       Rotation matrix
 **
 *****************************************************************************/
void GeometryHelper::rotationMatrixIdentityInPlace(int ndim, VectorDouble &rot)
{
  int i, j, ecr;

  for (i = ecr = 0; i < ndim; i++)
    for (j = 0; j < ndim; j++, ecr++)
      rot[ecr] = (i == j);
}

/****************************************************************************/
/*!
 **  Calculates the 2-D rotation matrix
 **
 ** \param[in]  angle Rotation angle (in degrees)
 **
 ** \param[out]  rot   Rotation matrix (Dimension = 4)
 **
 *****************************************************************************/
void GeometryHelper::rotation2DMatrixInPlace(double angle, VectorDouble &rot)
{
  double ca, sa;

  GH::rotationGetSinCos(angle, &ca, &sa);

  /* Define the 2-D rotation matrix */

  rot[0] =  ca;
  rot[1] =  sa;
  rot[2] = -sa;
  rot[3] =  ca;
}

/*****************************************************************************/
/*!
 **  Calculates the 3-D rotation matrix
 **
 ** \param[in]  alpha angle (in degrees) / oz
 ** \param[in]  beta  angle (in degrees) / oy'
 ** \param[in]  gamma angle (in degrees) / ox''
 **
 ** \param[out] rot   direct rotation matrix (Dimension = 9)
 **
 *****************************************************************************/
void GeometryHelper::rotation3DMatrixInPlace(double alpha,
                                             double beta,
                                             double gamma,
                                             VectorDouble &rot)
{
  double ca[3], sa[3];

  /* Initializations */

  GH::rotationGetSinCos(alpha, &ca[0], &sa[0]);
  GH::rotationGetSinCos(beta, &ca[1], &sa[1]);
  GH::rotationGetSinCos(gamma, &ca[2], &sa[2]);

  /* Define the 3-D rotation matrix */

  rot[0] = ca[0] * ca[1];
  rot[1] = sa[0] * ca[1];
  rot[2] = -sa[1];
  rot[3] = -sa[0] * ca[2] + ca[0] * sa[1] * sa[2];
  rot[4] = ca[0] * ca[2] + sa[0] * sa[1] * sa[2];
  rot[5] = ca[1] * sa[2];
  rot[6] = sa[0] * sa[2] + ca[0] * sa[1] * ca[2];
  rot[7] = -ca[0] * sa[2] + sa[0] * sa[1] * ca[2];
  rot[8] = ca[1] * ca[2];
}

/*****************************************************************************/
/*!
 **  Calculates the rotation matrix
 **
 ** \param[in]  ndim   Space dimension
 ** \param[in]  angles Array of angles
 **
 ** \param[out] rot   direct rotation matrix (dimensionned to ndim*ndim)
 **
 *****************************************************************************/
void GeometryHelper::rotationMatrixInPlace(int ndim,
                                           const VectorDouble &angles,
                                           VectorDouble &rot)
{
  if (ndim == 2)
    GH::rotation2DMatrixInPlace(angles[0], rot);
  else if (ndim == 3)
    GH::rotation3DMatrixInPlace(angles[0], angles[1], angles[2], rot);
  else
    GH::rotationMatrixIdentityInPlace(ndim, rot);
}

/*****************************************************************************/
/*!
 **  Calculates the rotation matrix.
 **  Returns the rotation matrix as a VectorDouble
 **
 ** \param[in]  ndim   Space dimension
 ** \param[in]  angles Array of angles
 **
 *****************************************************************************/
VectorDouble GeometryHelper::rotationMatrix(int ndim,
                                            const VectorDouble &angles)
{
  VectorDouble rot;

  rot.resize(ndim * ndim);
  if (ndim == 2)
    GH::rotation2DMatrixInPlace(angles[0], rot);
  else if (ndim == 3)
    GH::rotation3DMatrixInPlace(angles[0], angles[1], angles[2], rot);
  else
    GH::rotationMatrixIdentityInPlace(ndim, rot);

  return rot;
}

/*****************************************************************************/
/*!
 **  Copy a rotation matrix
 **
 ** \param[in]  ndim   Space dimension
 ** \param[in]  rotin  Input rotation matrix
 **
 ** \param[out] rotout Output rotation matrix (already allocated)
 **
 *****************************************************************************/
void GeometryHelper::rotationCopy(int ndim, const double *rotin, double *rotout)
{
  for (int i = 0; i < ndim * ndim; i++)
    rotout[i] = rotin[i];
}

void GeometryHelper::rotationGetAnglesInPlace(const VectorDouble &rot,
                                              VectorDouble &angles)
{
  int ndim = sqrt((int) rot.size());
  GH::rotationGetAnglesInPlace(ndim, rot.data(), angles.data());
}

/****************************************************************************/
/*!
 **  Calculates the rotation angles from the rotation matrix
 **
 ** \param[in]  ndim  Space dimension
 ** \param[in]  rot   Rotation matrix (Dimension = 9)
 **
 ** \param[out]  angles Rotation angles (Dimension = ndim)
 **
 *****************************************************************************/
void GeometryHelper::rotationGetAnglesInPlace(int ndim,
                                              const double *rot,
                                              double *angles)
{
  int nval;

  /* Initializations */

  for (int i = 0; i < ndim; i++)
    angles[i] = 0.;
  if (rot == nullptr) return;

  /* Dispatch */

  if (ndim == 1)
  {
    nval = 1;
  }
  else if (ndim == 2)
  {
    angles[0] = atan2(rot[1], rot[0]);
    nval = 1;
  }
  else if (ndim == 3)
  {
    nval = 3;
    angles[0] = atan2(rot[1], rot[0]);
    angles[1] = atan2(-rot[2], sqrt(rot[5] * rot[5] + rot[8] * rot[8]));
    angles[2] = atan2(rot[5], rot[8]);

//    double s1 = -rot[2];
//    double c1 = sqrt(rot[0] * rot[0] + rot[1] * rot[1]);
//    if (isZero(c1))
//    {
//      if (s1 > 0.)
//      {
//        angles[0] = 0.;
//        angles[1] = GV_PI / 2.;
//        angles[2] = atan2(rot[3], rot[6]);
//      }
//      else
//      {
//        angles[0] = 0.;
//        angles[1] = -GV_PI / 2.;
//        angles[2] = atan2(-rot[3], -rot[6]);
//      }
//    }
//    else
//    {
//      double c0 = rot[0] / c1;
//      double s0 = rot[1] / c1;
//      double s2 = rot[5] / c1;
//      double c2 = rot[8] / c1;
//      angles[0] = atan2(s0, c0);
//      angles[1] = atan2(s1, c1);
//      angles[2] = atan2(s2, c2);
//    }
  }
  else
  {
    nval = 0;
  }

  /* Convert into degrees */

  for (int i = 0; i < nval; i++)
    angles[i] = ut_rad2deg(angles[i]);
}

/****************************************************************************/
/*!
 **  Calculates the rotation angle from the direction coefficient
 **
 ** \param[in]  codir  Direction vector (Dimension = ndim)
 **
 ** \param[out]  angles Rotation angles (Dimension = ndim)
 **
 *****************************************************************************/
void GeometryHelper::rotationGetAnglesFromCodirInPlace(const VectorDouble &codir,
                                                       VectorDouble &angles)
{
  int nval;
  int ndim = (int) codir.size();
  angles.resize(ndim);

  /* Initializations */

  for (int i = 0; i < ndim; i++)
    angles[i] = 0.;

  /* Dispatch */

  if (ndim == 1)
  {
    return;
  }
  if (ndim == 2)
  {
    angles[0] = atan2(codir[1], codir[0]);
    angles[1] = 0.;
    nval = 1;
  }
  else if (ndim == 3)
  {
    double norme = codir[0] * codir[0] + codir[1] * codir[1];
    if (norme > 0.)
    {
      norme = sqrt(norme);
      angles[0] = atan2(codir[1] / norme, codir[0] / norme);
      angles[1] = atan2(codir[2], norme);
    }
    else
    {
      angles[2] = GV_PI / 2.;
    }
    nval = 3;
  }
  else
  {
    nval = 0;
  }

  /* Convert into degrees */

  for (int i = 0; i < nval; i++)
    angles[i] = ut_rad2deg(angles[i]);
}

/**
 * From the vector of direction coefficients (codir) returns the vector of angles
 *
 * @param codir Input vector giving the direction coefficients
 * @param flagResize When TRUE (and if in 2-D) the returned vector is resized to 1
 * @return
 */
VectorDouble GeometryHelper::rotationGetAngles(const VectorDouble &codir,
                                               bool flagResize)
{
  int ndim = (int) codir.size();
  VectorDouble angles(ndim);
  GH::rotationGetAnglesFromCodirInPlace(codir, angles);
  if (ndim == 2 && flagResize) angles.resize(1);
  return angles;
}

/****************************************************************************/
/*!
 **   Convert angles to a set of Directions
 **
 ** \param[in]  angles   Vector giving the angles characteristics (in degrees)
 **                      As this provides rotation in 2D, 'angles' is dimensioned
 **                      to 'ndir' (one angle requires a single value)
 **
 ** \param[out] codir    Vector of the direction (Dim: ndir * ndim)
 **
 ** \remarks 'ndir' is given by the dimension of 'angles', 'ndim' is given by 'codir'
 **
 *****************************************************************************/
void GeometryHelper::rotationGetDirection2D(const VectorDouble &angles,
                                            VectorDouble &codir)
{
  int ndir = (int) angles.size();
  int ndim = 2;
  if (codir.size() > 0) ndim = (int) codir.size() / ndir;
  codir.resize(ndim * ndir);
  for (int i = 0; i < ndim * ndir; i++)
    codir[i] = 0.;

  for (int idir = 0; idir < ndir; idir++)
  {
    codir[idir * ndim + 0] = cos(angles[idir] * GV_PI / 180.);
    codir[idir * ndim + 1] = sin(angles[idir] * GV_PI / 180.);
  }
}

/****************************************************************************/
/*!
 **   Create a Direction (used as default)
 **
 ** \param[in]  ndim     Number of space dimensions
 **
 ** \param[out] codir    Vector of the direction (Dim: ndim)
 **
 *****************************************************************************/
void GeometryHelper::rotationGetDirectionDefault(int ndim, VectorDouble &codir)
{
  codir.resize(ndim, 0.);
  codir[0] = 1.;
}

/****************************************************************************/
/*!
 **  Starting from a rotation matrix, check it is different from the Identity
 **
 ** \return  False if a rotation is defined; True if it is an Identity
 **
 ** \param[in]  ndim     Space dimension
 ** \param[in]  rot      Rotation matrix
 ** \param[in]  eps      Tolerance
 **
 *****************************************************************************/
bool GeometryHelper::rotationIsIdentity(int ndim, const double *rot, double eps)
{
  for (int i = 0; i < ndim; i++)
    for (int j = 0; j < ndim; j++)
    {
      double rotval = rot[(i) * ndim + (j)];
      if (i == j)
      {
        if (ABS(rotval - 1.) > eps) return false;
      }
      else
      {
        if (ABS(rotval) > eps) return false;
      }
    }
  return true;
}

/****************************************************************************/
/*!
 **  Find the shortest distance between the point (x0,y0) and the segment
 **  with the two end points (x1,y1) and (x2,y2)
 **
 ** \return Minimum algebraic distance (positive or negative)
 **
 ** \param[in]  x0,y0   Coordinates of the target point
 ** \param[in]  x1,y1   Coordinate of the first end-point of the segment
 ** \param[in]  x2,y2   Coordinate of the second end-point of the segment
 **
 ** \param[out] xd,yd   Coordinates of the closest point
 ** \param[out] nint    =1 if the projection belongs to the segment
 **                     =0 if it is set to one of the segment vertices
 **
 *****************************************************************************/
double GeometryHelper::distancePointToSegment(double x0,
                                              double y0,
                                              double x1,
                                              double y1,
                                              double x2,
                                              double y2,
                                              double *xd,
                                              double *yd,
                                              int *nint)
{
  double dx = x2 - x1;
  double dy = y2 - y1;

  double ratio = (dx * (x0 - x1) + dy * (y0 - y1)) / (dx * dx + dy * dy);

  if (ratio < 0)
  {
    *xd = x1;
    *yd = y1;
    *nint = 0;
  }
  else if (ratio > 1)
  {
    *xd = x2;
    *yd = y2;
    *nint = 0;
  }
  else
  {
    *xd = x1 + ratio * dx;
    *yd = y1 + ratio * dy;
    *nint = 1;
  }

  double dxp = x0 - (*xd);
  double dyp = y0 - (*yd);
  double dist = sqrt(dxp * dxp + dyp * dyp);
  double signe = (dy * (x0 - x1) - dx * (y0 - y1) > 0.) ? 1 :
                                                          -1;

  return (dist * signe);
}

/****************************************************************************/
/*!
 **  Calculate the geodetic angular distance between two points on the sphere
 **
 ** \return Angular distance
 **
 ** \param[in]  long1  Longitude of the first point (in degrees)
 ** \param[in]  lat1   Latitude of the first point (in degrees)
 ** \param[in]  long2  Longitude of the second point (in degrees)
 ** \param[in]  lat2   Latitude of the second point (in degrees)
 ** \param[in]  radius Radius of the sphere
 **
 *****************************************************************************/
double GeometryHelper::geodeticAngularDistance(double long1,
                                               double lat1,
                                               double long2,
                                               double lat2,
                                               double radius)
{
  double rlon1 = ut_deg2rad(long1);
  double rlat1 = ut_deg2rad(lat1);
  double rlon2 = ut_deg2rad(long2);
  double rlat2 = ut_deg2rad(lat2);
  double dlong = rlon2 - rlon1;
  double angdst = acos(
      sin(rlat1) * sin(rlat2) + cos(rlat1) * cos(rlat2) * cos(dlong));
  return (radius * angdst);
}

/****************************************************************************/
/*!
 **  Extract A from a,b,c
 **
 ** \param[in]  cosa   Cosine of first angle
 ** \param[in]  sinb   Sine of second angle
 ** \param[in]  cosb   Cosine of second angle
 ** \param[in]  sinc   Sine of third angle
 ** \param[in]  cosc   Cosine of third angle
 **
 *****************************************************************************/
static double st_convert_geodetic_angle(double /*sina*/,
                                        double cosa,
                                        double sinb,
                                        double cosb,
                                        double sinc,
                                        double cosc)
{
  double prod = sinb * sinc;
  double cosA = (prod == 0.) ? 0. :
                               (cosa - cosb * cosc) / prod;
  if (cosA < -1) cosA = -1.;
  if (cosA > +1) cosA = +1.;
  return (acos(cosA));
}

/****************************************************************************/
/*!
 **  Calculate all geodetic angles from a spherical triangle
 **
 ** \param[in]  long1  Longitude of the first point (in degrees)
 ** \param[in]  lat1   Latitude of the first point (in degrees)
 ** \param[in]  long2  Longitude of the second point (in degrees)
 ** \param[in]  lat2   Latitude of the second point (in degrees)
 ** \param[in]  long3  Longitude of the third point (in degrees)
 ** \param[in]  lat3   Latitude of the third point (in degrees)
 **
 ** \param[out] a      Angle (P2,O,P3)
 ** \param[out] b      Angle (P3,O,P1)
 ** \param[out] c      Angle (P1,O,P2)
 ** \param[out] A      Angle (P2,P1,P3)
 ** \param[out] B      Angle (P3,P2,P1)
 ** \param[out] C      Angle (P1,P3,P2)
 **
 *****************************************************************************/
void GeometryHelper::geodeticAngles(double long1,
                                    double lat1,
                                    double long2,
                                    double lat2,
                                    double long3,
                                    double lat3,
                                    double *a,
                                    double *b,
                                    double *c,
                                    double *A,
                                    double *B,
                                    double *C)
{
  *a = GH::geodeticAngularDistance(long2, lat2, long3, lat3);
  *b = GH::geodeticAngularDistance(long1, lat1, long3, lat3);
  *c = GH::geodeticAngularDistance(long1, lat1, long2, lat2);

  double cosa = cos(*a);
  double cosb = cos(*b);
  double cosc = cos(*c);
  double sina = sin(*a);
  double sinb = sin(*b);
  double sinc = sin(*c);

  *A = st_convert_geodetic_angle(sina, cosa, sinb, cosb, sinc, cosc);
  *B = st_convert_geodetic_angle(sinb, cosb, sinc, cosc, sina, cosa);
  *C = st_convert_geodetic_angle(sinc, cosc, sina, cosa, sinb, cosb);
}

/****************************************************************************/
/*!
 **  Calculate the perimeter of the spherical triangle
 **
 ** \return The Perimeter of the spherical triangle
 **
 ** \param[in]  long1  Longitude of the first point (in degrees)
 ** \param[in]  lat1   Latitude of the first point (in degrees)
 ** \param[in]  long2  Longitude of the second point (in degrees)
 ** \param[in]  lat2   Latitude of the second point (in degrees)
 ** \param[in]  long3  Longitude of the third point (in degrees)
 ** \param[in]  lat3   Latitude of the third point (in degrees)
 **
 *****************************************************************************/
double GeometryHelper::geodeticTrianglePerimeter(double long1,
                                                 double lat1,
                                                 double long2,
                                                 double lat2,
                                                 double long3,
                                                 double lat3)
{
  double a, b, c, ga, gb, gc;
  GH::geodeticAngles(long1, lat1, long2, lat2, long3, lat3, &a, &b, &c, &ga,
                     &gb, &gc);
  return (a + b + c);
}

/****************************************************************************/
/*!
 **  Calculate the surface of the spherical triangle
 **
 ** \return The Surface of the spherical triangle (with unit radius)
 **
 ** \param[in]  long1  Longitude of the first point (in degrees)
 ** \param[in]  lat1   Latitude of the first point (in degrees)
 ** \param[in]  long2  Longitude of the second point (in degrees)
 ** \param[in]  lat2   Latitude of the second point (in degrees)
 ** \param[in]  long3  Longitude of the third point (in degrees)
 ** \param[in]  lat3   Latitude of the third point (in degrees)
 **
 *****************************************************************************/
double GeometryHelper::geodeticTriangleSurface(double long1,
                                               double lat1,
                                               double long2,
                                               double lat2,
                                               double long3,
                                               double lat3)
{
  double a, b, c, A, B, C;

  GH::geodeticAngles(long1, lat1, long2, lat2, long3, lat3, &a, &b, &c, &A, &B,
                     &C);
  return (A + B + C - GV_PI);
}

/****************************************************************************/
/*!
 **  Is a point inside a spherical triangle
 **
 ** \return True if the point belongs to the spherical triangle; False otherwise
 **
 ** \param[in]  coor    Coordinates of the target point (long,lat)
 ** \param[in]  ptsa    Coordinates of the first point of the triangle
 ** \param[in]  ptsb    Coordinates of the second point of the triangle
 ** \param[in]  ptsc    Coordinates of the third point of the triangle
 ** \param[in]  eps     Tolerance
 **
 ** \param[out] wgts    Array of weights
 **
 *****************************************************************************/
bool GeometryHelper::isInSphericalTriangleOptimized(const double *coor,
                                                    double *ptsa,
                                                    double *ptsb,
                                                    double *ptsc,
                                                    double *wgts,
                                                    double eps)
{
  double s[3];
  double total = 0.;

  double dab = GH::geodeticAngularDistance(ptsa[0], ptsa[1], ptsb[0], ptsb[1]);
  double dbc = GH::geodeticAngularDistance(ptsb[0], ptsb[1], ptsc[0], ptsc[1]);
  double dac = GH::geodeticAngularDistance(ptsa[0], ptsa[1], ptsc[0], ptsc[1]);
  double d0a = GH::geodeticAngularDistance(coor[0], coor[1], ptsa[0], ptsa[1]);
  double d0b = GH::geodeticAngularDistance(coor[0], coor[1], ptsb[0], ptsb[1]);
  double d0c = GH::geodeticAngularDistance(coor[0], coor[1], ptsc[0], ptsc[1]);

  double sinab = sin(dab);
  double cosab = cos(dab);
  double sinbc = sin(dbc);
  double cosbc = cos(dbc);
  double sinac = sin(dac);
  double cosac = cos(dac);
  double sin0a = sin(d0a);
  double cos0a = cos(d0a);
  double sin0b = sin(d0b);
  double cos0b = cos(d0b);
  double sin0c = sin(d0c);
  double cos0c = cos(d0c);

  double A = st_convert_geodetic_angle(sinbc, cosbc, sinac, cosac, sinab,
                                       cosab);
  double B = st_convert_geodetic_angle(sinac, cosac, sinab, cosab, sinbc,
                                       cosbc);
  double C = st_convert_geodetic_angle(sinab, cosab, sinbc, cosbc, sinac,
                                       cosac);
  double stot = (A + B + C - GV_PI);

  double OA = st_convert_geodetic_angle(sinbc, cosbc, sin0c, cos0c, sin0b,
                                        cos0b);
  double BA = st_convert_geodetic_angle(sin0c, cos0c, sin0b, cos0b, sinbc,
                                        cosbc);
  double CA = st_convert_geodetic_angle(sin0b, cos0b, sinbc, cosbc, sin0c,
                                        cos0c);
  s[0] = (OA + BA + CA - GV_PI);
  total += s[0];
  if (total > stot + eps) return false;

  double AB = st_convert_geodetic_angle(sin0c, cos0c, sinac, cosac, sin0a,
                                        cos0a);
  double OB = st_convert_geodetic_angle(sinac, cosac, sin0a, cos0a, sin0c,
                                        cos0c);
  double CB = st_convert_geodetic_angle(sin0a, cos0a, sin0c, cos0c, sinac,
                                        cosac);
  s[1] = (AB + OB + CB - GV_PI);
  total += s[1];
  if (total > stot + eps) return false;

  double AC = st_convert_geodetic_angle(sin0b, cos0b, sin0a, cos0a, sinab,
                                        cosab);
  double BC = st_convert_geodetic_angle(sin0a, cos0a, sinab, cosab, sin0b,
                                        cos0b);
  double OC = st_convert_geodetic_angle(sinab, cosab, sin0b, cos0b, sin0a,
                                        cos0a);
  s[2] = (AC + BC + OC - GV_PI);
  total += s[2];
  if (ABS(total - stot) > eps) return false;

  for (int i = 0; i < 3; i++)
    wgts[i] = s[i] / total;
  return true;
}

/****************************************************************************/
/*!
 **  Check if two 2-D segments intersect
 **
 ** \return True if there is an intersection; False otherwise
 **
 ** \param[in]  xd1,yd1     Starting point for the first segment
 ** \param[in]  xe1,ye1     Ending point for the first segment
 ** \param[in]  xd2,yd2     Starting point for the second segment
 ** \param[in]  xe2,ye2     Ending point for the second segment
 **
 ** \param[out]   xint,yint  Coordinates of the intersection
 **
 *****************************************************************************/
bool GeometryHelper::segmentIntersect(double xd1,
                                      double yd1,
                                      double xe1,
                                      double ye1,
                                      double xd2,
                                      double yd2,
                                      double xe2,
                                      double ye2,
                                      double *xint,
                                      double *yint)
{
  *xint = TEST;
  *yint = TEST;

  if (MAX(xd2, xe2) < MIN(xd1, xe1)) return false;
  if (MAX(xd1, xe1) < MIN(xd2, xe2)) return false;
  if (MAX(yd2, ye2) < MIN(yd1, ye1)) return false;
  if (MAX(yd1, ye1) < MIN(yd2, ye2)) return false;

  const double b1 = ye1 - yd1;
  const double b2 = ye2 - yd2;

  const bool b1_is_not_zero = b1 * b1 >= EPSILON20;
  const bool b2_is_not_zero = b2 * b2 >= EPSILON20;

  if (b1_is_not_zero && b2_is_not_zero)
  {

    /* This operation is safe as end-point coordinates cannot be equal */

    const double a1 = (xe1 - xd1) / b1;
    const double a2 = (xe2 - xd2) / b2;

    /* Skip the case of parallel segments */

    const double delta_a = a1 - a2;
    if (delta_a * delta_a < EPSILON20) return false;

    /* Discard intersection if located outside the segment */

    const double y = (xd2 - xd1 + a1 * yd1 - a2 * yd2) / delta_a;
    if ((y - yd1) * (y - ye1) > 0) return false;
    if ((y - yd2) * (y - ye2) > 0) return false;
    (*xint) = xd1 + a1 * (y - yd1);
    (*yint) = y;
    return true;
  }

  /* Case of two horizontal segments */

  if (!b1_is_not_zero && !b2_is_not_zero)
  {
    const double delta_y = ye1 - ye2;
    if (delta_y * delta_y > EPSILON20) return false;
    double x1m = MIN(xd1, xe1);
    double x1M = MAX(xd1, xe1);
    double x2m = MIN(xd2, xe2);
    double x2M = MAX(xd2, xe2);
    if (x1m > x2M || x2m > x1M) return false;
    (*xint) = MAX(x1m, x2m);
    (*yint) = ye1;
    return true;
  }

  /* Case of the horizontal first segment */

  if (b2_is_not_zero)
  {
    const double y = ye1;
    if ((y - yd2) * (y - ye2) > 0) return false;
    const double x = xe2 + (y - ye2) * (xe2 - xd2) / b2;
    if ((x - xd1) * (x - xe1) > 0) return false;
    if ((x - xd2) * (x - xe2) > 0) return false;
    (*xint) = x;
    (*yint) = y;
    return true;
  }

  /* Case of horizontal second segment */

  const double y = ye2;
  if ((y - yd1) * (y - ye1) > 0) return false;
  const double x = xe1 + (y - ye1) * (xe1 - xd1) / b1;
  if ((x - xd1) * (x - xe1) > 0) return false;
  if ((x - xd2) * (x - xe2) > 0) return false;
  (*xint) = x;
  (*yint) = y;
  return true;
}

/****************************************************************************/
/*!
 **  Is a point inside a spherical triangle
 **
 ** \return 1 if the point belongs to the spherical triangle; 0 otherwise
 **
 ** \param[in]  coor    Coordinates of the target point (long,lat)
 ** \param[in]  surface Surface of the spherical triangle
 ** \param[in]  pts1    Coordinates of the first point of the triangle
 ** \param[in]  pts2    Coordinates of the second point of the triangle
 ** \param[in]  pts3    Coordinates of the third point of the triangle
 ** \param[in]  eps     Tolerance
 **
 ** \param[out] wgts    Array of weights
 **
 *****************************************************************************/
bool GeometryHelper::isInSphericalTriangle(double *coor,
                                           double surface,
                                           double *pts1,
                                           double *pts2,
                                           double *pts3,
                                           double *wgts,
                                           double eps)
{
  double s[3];

  double total = 0.;
  s[0] = GH::geodeticTriangleSurface(coor[0], coor[1], pts2[0], pts2[1],
                                     pts3[0], pts3[1]);
  total += s[0];
  if (total > surface + eps) return false;
  s[1] = GH::geodeticTriangleSurface(pts1[0], pts1[1], coor[0], coor[1],
                                     pts3[0], pts3[1]);
  total += s[1];
  if (total > surface + eps) return false;
  s[2] = GH::geodeticTriangleSurface(pts1[0], pts1[1], pts2[0], pts2[1],
                                     coor[0], coor[1]);
  total += s[2];
  if (ABS(total - surface) > eps) return false;
  for (int i = 0; i < 3; i++)
    wgts[i] = s[i] / total;
  return true;
}

/****************************************************************************/
/*!
 **  Rotation of a Direction in 3-D
 **
 ** \param[in]  ct,st     Cosine and Sine of the rotation angle
 ** \param[in]  a         Random direction
 ** \param[in,out] codir  Direction to be rotated
 **
 *****************************************************************************/
void GeometryHelper::rotationGetRandomDirection(double ct,
                                                double st,
                                                const double *a,
                                                double *codir)
{
  double b[3], c[3], p[3];

  double rd = 0.;
  for (int k = 0; k < 3; k++)
    rd += codir[k] * a[k];
  for (int k = 0; k < 3; k++)
    p[k] = rd * a[k];
  for (int k = 0; k < 3; k++)
    b[k] = codir[k] - p[k];

  rd = 0.;
  for (int k = 0; k < 3; k++)
    rd += b[k] * b[k];

  rd = sqrt(rd);
  for (int k = 0; k < 3; k++)
    b[k] /= rd;

  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];

  for (int k = 0; k < 3; k++)
    codir[k] = p[k] + rd * (ct * b[k] + st * c[k]);
}

/**
 * Returns the Vector of Sample coordinates in 3-D from Longitude-Latitude
 * @param longitude Array of longitude values
 * @param latitude  Array of latitude values
 * @param dilate    Dilation applied to radius
 * @param radius_arg    Radius (if note defined, taken from variety definition)
 * @return
 */
VectorVectorDouble GeometryHelper::convertLongLat(const VectorDouble &longitude,
                                                  const VectorDouble &latitude,
                                                  double dilate,
                                                  double radius_arg)
{
  double radius = radius_arg;
  if (FFFF(radius))
  {
    const ASpace *space = getDefaultSpace();
    const SpaceSN *spaceSn = dynamic_cast<const SpaceSN*>(space);
    if (spaceSn != nullptr) radius = spaceSn->getRadius();
  }
  radius *= dilate;

  VectorVectorDouble tab;
  int number = (int) longitude.size();
  if (number != (int) latitude.size()) return tab;

  // Dimension the returned argument
  tab.resize(3);
  for (int idim = 0; idim < 3; idim++)
    tab[idim].resize(number, 0.);

  // Load the returned argument

  for (int ip = 0; ip < number; ip++)
  {
    double lon = longitude[ip];
    double lat = latitude[ip];
    if (FFFF(lon) || FFFF(lat))
    {
      tab[0][ip] = TEST;
      tab[1][ip] = TEST;
      tab[2][ip] = TEST;
    }
    else
    {
      lon = ut_deg2rad(lon);
      lat = ut_deg2rad(lat);
      tab[0][ip] = radius * cos(lon) * cos(lat);
      tab[1][ip] = radius * sin(lon) * cos(lat);
      tab[2][ip] = radius * sin(lat);
    }
  }
  return tab;
}

/****************************************************************************/
/*!
 **  Convert the cartesian coordinates into spherical coordinates
 **
 ** \param[in]  x     First cartesian coordinate
 ** \param[in]  y     Second cartesian coordinate
 ** \param[in]  z     Third cartesian coordinate
 ** \param[in]  radius_arg Radius of the sphere (Earth if TEST)
 **
 ** \param[out] rlong Longitude (in degrees)
 ** \param[out] rlat  Latitude (in degrees)
 **
 *****************************************************************************/
void GeometryHelper::convertCart2Sph(double x,
                                     double y,
                                     double z,
                                     double *rlong,
                                     double *rlat,
                                     double radius_arg)
{
  double radius = radius_arg;
  if (FFFF(radius))
  {
    const ASpace *space = getDefaultSpace();
    const SpaceSN *spaceSn = dynamic_cast<const SpaceSN*>(space);
    if (spaceSn != nullptr) radius = spaceSn->getRadius();
  }

  double loc_long, loc_lat;

  x /= radius;
  y /= radius;
  z /= radius;

  loc_long = ut_rad2deg(atan2(y, x));
  loc_lat = ut_rad2deg(asin(z));

  if (loc_long < 0.)
    loc_long += 360.;
  else if (loc_long > 360.) loc_long -= 360.;
  if (loc_lat < -90.)
    loc_lat += 180.;
  else if (loc_lat > 90.) loc_lat -= 180.;

  *rlong = loc_long;
  *rlat = loc_lat;
}

/****************************************************************************/
/*!
 **  Convert the spherical coordinates into cartesian coordinates
 **
 ** \param[in]  rlong Longitude (in degrees)
 ** \param[in]  rlat  Latitude (in degrees)
 ** \param[in]  radius_arg radius of the sphere (Earth if TEST)
 **
 ** \param[out] x     First cartesian coordinate
 ** \param[out] y     Second cartesian coordinate
 ** \param[out] z     Third cartesian coordinate
 **
 *****************************************************************************/
void GeometryHelper::convertSph2Cart(double rlong,
                                     double rlat,
                                     double *x,
                                     double *y,
                                     double *z,
                                     double radius_arg)
{
  double radius = radius_arg;
  if (FFFF(radius))
  {
    const ASpace *space = getDefaultSpace();
    const SpaceSN *spaceSn = dynamic_cast<const SpaceSN*>(space);
    if (spaceSn != nullptr) radius = spaceSn->getRadius();
  }

  double phi, theta, sinphi, cosphi, sinthe, costhe;

  phi = ut_deg2rad(rlat);
  theta = ut_deg2rad(rlong);
  sinphi = sin(phi);
  cosphi = cos(phi);
  sinthe = sin(theta);
  costhe = cos(theta);

  *x = radius * cosphi * costhe;
  *y = radius * cosphi * sinthe;
  *z = radius * sinphi;
}

/**
 * Calculate the angle to translate the gradients of a tilted plane
 * (defined by the partial derivatives along X and Y) into the axis of the rotation angle
 * @param dzoverdx Partial derivative along X
 * @param dzoverdy Partial derivative along Y
 */
double util_rotation_gradXYToAngle(double dzoverdx, double dzoverdy)
{
  int ndim = 3;

  // Vector orthogonal to the horizontal plane/ i.e. the vertical axis
  VectorDouble vert = VectorDouble(3, 0.);
  vert[ndim - 1] = -1.;

  // Vector orthogonal to the tilted plane
  VectorDouble vort = VectorDouble(ndim);
  vort[0] = dzoverdx;
  vort[1] = dzoverdy;
  vort[2] = -1.;
  double norme = VH::norm(vort);
  for (int idim = 0; idim < ndim; idim++)
    vort[idim] /= norme;

  // Cross product
  VectorDouble axis(ndim, 0.);
  axis[0] = vort[1];
  axis[1] = -vort[0];
  axis[2] = 0.;

  // Norm of the cross product and dot product between ref and vort
  double normcross = sqrt(VH::innerProduct(axis, axis));
  double dot = VH::innerProduct(vert, vort);

  // Rotation angle
  double angle = atan2(normcross, dot);
  return angle;
}

/**
 * Calculate the rotation matrix starting from the partial derivatives along X and Y
 * of a tilted plane
 * @param dzoverdx Partial derivative along X
 * @param dzoverdy Partial derivative along Y
 */
MatrixSquareGeneral GeometryHelper::gradXYToRotmat(double dzoverdx,
                                                   double dzoverdy)
{
  int ndim = 3;
  VectorDouble axis(ndim, 0.);

  // Vector orthogonal to the horizontal plane/ i.e. the vertical axis
  VectorDouble vert = VectorDouble(3, 0.);
  vert[ndim - 1] = -1.;

  // Vector orthogonal to the tilted plane
  VectorDouble vort = VectorDouble(ndim);
  vort[0] = dzoverdx;
  vort[1] = dzoverdy;
  vort[2] = -1.;
  double norme = VH::norm(vort);
  for (int idim = 0; idim < ndim; idim++)
    vort[idim] /= norme;

  // Cross product
  axis[0] = vort[1];
  axis[1] = -vort[0];
  axis[2] = 0.;

  // Norm of the cross product and dot product between ref and vort
  double normcross = sqrt(VH::innerProduct(axis, axis));
  double dot = VH::innerProduct(vert, vort);

  // Rotation axis (normalized
  for (int idim = 0; idim < ndim; idim++)
    axis[idim] /= normcross;

  // Rotation angle
  double angle = atan2(normcross, dot);
  double c = cos(angle);
  double s = sin(angle);

  MatrixSquareGeneral M(ndim);
  M.setValue(0, 0, axis[0] * axis[0] * (1 - c) + c);
  M.setValue(0, 1, axis[0] * axis[1] * (1 - c) - axis[2] * s);
  M.setValue(0, 2, axis[0] * axis[2] * (1 - c) + axis[1] * s);
  M.setValue(1, 0, axis[1] * axis[0] * (1 - c) + axis[2] * s);
  M.setValue(1, 1, axis[1] * axis[1] * (1 - c) + c);
  M.setValue(1, 2, axis[1] * axis[2] * (1 - c) - axis[0] * s);
  M.setValue(2, 0, axis[2] * axis[0] * (1 - c) - axis[1] * s);
  M.setValue(2, 1, axis[2] * axis[1] * (1 - c) + axis[0] * s);
  M.setValue(2, 2, axis[2] * axis[2] * (1 - c) + c);
  return M;
}

/**
 * Returns the Euler angles, starting from a rotation matrix
 * @param M Input matrix
 * @param convrot Rotation convention
 * @param eps Tolerance
 * @return
 *
 * @remark The code is coming from the following reference (BSD license)
 * @remark https://github.com/matthew-brett/transforms3d/blob/master/transforms3d/euler.py
 */
VectorDouble GeometryHelper::rotationToEuler(const MatrixSquareGeneral &M,
                                             const ERotation &convrot,
                                             double eps)
{
  int firstaxis, parity, repetition, frame;
  _decodeConvRot(convrot, &firstaxis, &parity, &repetition, &frame);

  VectorInt next_axis = { 1, 2, 0, 1 };
  int i = firstaxis;
  int j = next_axis[i + parity];
  int k = next_axis[i - parity + 1];

  double ax, ay, az;

  if (repetition)
  {
    double sy = sqrt(
        M.getValue(i, j) * M.getValue(i, j) + M.getValue(i, k)
            * M.getValue(i, k));
    if (sy > eps)
    {
      ax = atan2(M.getValue(i, j), M.getValue(i, k));
      ay = atan2(sy, M.getValue(i, i));
      az = atan2(M.getValue(j, i), -M.getValue(k, i));
    }
    else
    {
      ax = atan2(-M.getValue(j, k), M.getValue(j, j));
      ay = atan2(sy, M.getValue(i, i));
      az = 0.0;
    }
  }
  else
  {
    double cy = sqrt(
        M.getValue(i, i) * M.getValue(i, i) + M.getValue(j, i)
            * M.getValue(j, i));
    if (cy > eps)
    {
      ax = atan2(M.getValue(k, j), M.getValue(k, k));
      ay = atan2(-M.getValue(k, i), cy);
      az = atan2(M.getValue(j, i), M.getValue(i, i));
    }
    else
    {
      ax = atan2(-M.getValue(j, k), M.getValue(j, j));
      ay = atan2(-M.getValue(k, i), cy);
      az = 0.0;
    }
  }

  if (parity)
  {
    ax = -ax;
    ay = -ay;
    az = -az;
  }

  if (frame)
  {
    ax = az;
    az = ax;
  }

  VectorDouble angles = { ax, ay, az };
  return angles;
}

/**
 * Returns the Rotation matrix, starting from the Euler angles
 * @param angles Ordered list of Euler angles
 * @param convrot Rotation convention
 * @return
 *
 * @remark The code is coming from the following reference (BSD license)
 * @remark https://github.com/matthew-brett/transforms3d/blob/master/transforms3d/euler.py
 */
MatrixSquareGeneral GeometryHelper::EulerToRotation(const VectorDouble &angles,
                                                    const ERotation &convrot)
{
  int firstaxis, parity, repetition, frame;
  _decodeConvRot(convrot, &firstaxis, &parity, &repetition, &frame);

  VectorInt next_axis = { 1, 2, 0, 1 };
  int i = firstaxis;
  int j = next_axis[i + parity];
  int k = next_axis[i - parity + 1];

  int ndim = 3;
  MatrixSquareGeneral M(ndim);

  double ai = angles[0];
  double aj = angles[1];
  double ak = angles[2];
  if (frame)
  {
    ai = ak;
    ak = ai;
  }
  if (parity)
  {
    ai = -ai;
    aj = -aj;
    ak = -ak;
  }

  double si = sin(ai);
  double sj = sin(aj);
  double sk = sin(ak);
  double ci = cos(ai);
  double cj = cos(aj);
  double ck = cos(ak);

  double cc = ci * ck;
  double cs = ci * sk;
  double sc = si * ck;
  double ss = si * sk;

  if (repetition)
  {
    M.setValue(i, i, cj);
    M.setValue(i, j, sj * si);
    M.setValue(i, k, sj * ci);
    M.setValue(j, i, sj * sk);
    M.setValue(j, j, -cj * ss + cc);
    M.setValue(j, k, -cj * cs - sc);
    M.setValue(k, i, -sj * ck);
    M.setValue(k, j, cj * sc + cs);
    M.setValue(k, k, cj * cc - ss);
  }
  else
  {
    M.setValue(i, i, cj * ck);
    M.setValue(i, j, sj * sc - cs);
    M.setValue(i, k, sj * cc + ss);
    M.setValue(j, i, cj * sk);
    M.setValue(j, j, sj * ss + cc);
    M.setValue(j, k, sj * cs - sc);
    M.setValue(k, i, -sj);
    M.setValue(k, j, cj * si);
    M.setValue(k, k, cj * ci);
  }
  return M;
}

MatrixRectangular* GeometryHelper::getDirectionsInR3(const MatrixRectangular *U)
{
  int np = U->getNRows();
  int nd = U->getNCols();
  if (np <= 0)
  {
    messerr("The number of samples must be positive");
    return nullptr;
  }
  if (nd < 2)
  {
    messerr("This method requires at least 2 columns in U");
    return nullptr;
  }

  MatrixRectangular *X = new MatrixRectangular(np, 3);
  for (int ip = 0; ip < np; ip++)
  {
    double u1 = U->getValue(ip, 0);
    double u2 = U->getValue(ip, 1);
    X->setValue(ip, 0, cos(2 * GV_PI * u1) * sqrt(1 - u2 * u2));
    X->setValue(ip, 1, sin(2 * GV_PI * u1) * sqrt(1 - u2 * u2));
    X->setValue(ip, 2, u2);
  }
  return X;
}

/**
 * Function to compute directions in Rd
 * @param U a matrix [n, d] of uniform values between [0,1].
 * 'n' is the number of direction and 'd' is the dimension of the space
 * @return A vector of a matrix [n, d] of the coordinates of the 'n' directions
 * in the Euclidean space Rd.
 */
MatrixRectangular* GeometryHelper::getDirectionsInRn(const MatrixRectangular *U)
{
  int np = U->getNRows();
  int nd = U->getNCols();
  if (np <= 0)
  {
    messerr("The number of samples must be positive");
    return nullptr;
  }
  if (nd <= 0)
  {
    messerr("This method requires several columns in U");
    return nullptr;
  }

  // Check that all values in U lie in [0,1] interval
  if (U->getMinimum() < 0. || U->getMaximum() > 1.)
  {
    messerr("The argument 'U' must contain values lying within [0,1]");
    return nullptr;
  }

  MatrixRectangular *Y = new MatrixRectangular(np, nd);
  for (int ip = 0; ip < np; ip++)
  {
    double total = 0.;
    for (int id = 0; id < nd; id++)
    {
      double value = law_invcdf_gaussian(U->getValue(ip, id));
      Y->setValue(ip, id, value);
      total += value * value;
    }
    total = sqrt(total);
    for (int id = 0; id < nd; id++)
    {
      double value = Y->getValue(ip, id);
      value /= total;
      if (id == nd - 1) value = ABS(value);
      Y->setValue(ip, id, value);
    }
  }
  return Y;
}

/**
 * Format an Angle to be lying in [0, basis]
 * @param anglein Input angle value
 * @param basis   Basis (should be 360 by default; could be 180 for variogram or covariance)
 * @return
 */
double GeometryHelper::formatAngle(double anglein, double basis)
{
  double angle = anglein;
  if (angle < 0)
  {
    while (angle < 0.)
      angle += basis;
  }
  else
  {
    while (angle > basis)
      angle -= basis;
  }
  return angle;
}

VectorDouble GeometryHelper::formatAngles(const VectorDouble &anglesin, double basis)
{
  VectorDouble angles = anglesin;
  for (int idim = 0; idim < (int) angles.size(); idim++)
    angles[idim] = formatAngle(angles[idim], basis);
  return angles;
}

VectorDouble GeometryHelper::rayTriangleIntersect(const VectorDouble &dir,
                                                  const VectorDouble &v0,
                                                  const VectorDouble &v1,
                                                  const VectorDouble &v2)
{
  VectorDouble res(3, -1.);
  VectorDouble v20 = VH::subtract(v0, v2);
  VectorDouble v10 = VH::subtract(v0, v1);

  VectorDouble pvec = VH::crossProduct3D(dir, v20);
  double det = VH::innerProduct(pvec, v10);

  if (det == 0.) return res;

  double invDet = 1. / det;

  double u = -VH::innerProduct(pvec, v0) * invDet;
  if (u < 0 || u > 1) return res;

  VectorDouble vm0 = v0;
  VH::multiplyConstant(vm0, -1.);
  VectorDouble qvec = VH::crossProduct3D(vm0, v10);
  double v = VH::innerProduct(dir, qvec) * invDet;
  if (v < 0 || u + v > 1) return res;

  res[0] = VH::innerProduct(qvec, v20) * invDet;
  res[1] = u;
  res[2] = v;

  return res;
}

/**
 * Calculate the Barycenter coordinates of points with in the Spherical Meshing
 * @param sphPts Coordinates of target samples (dim = np * 3)
 * @param apices Coordinates of the apices of the Meshing
 * @param meshes Mesh information of the Meshing
 * @return A vector of barycenters (dimension = 4 * np) where
 * - 1: rank of the triangle
 * - 2, 3, 4: barycentric coordinates
 */
VectorVectorDouble GeometryHelper::sphBarCoord(const VectorVectorDouble &sphPts,
                                               const MatrixRectangular &apices,
                                               const MatrixInt &meshes)
{
  int np = (int) sphPts.size();
  int nmeshes = meshes.getNRows();

  // Dimension the output storage
  VectorVectorDouble res(4);
  for (int i = 0; i < 4; i++)
    res[i].resize(np, 0.);

  int iv0, iv1, iv2;
  VectorDouble w;

  for (int k = 0; k < np; k++)
  {
    bool notFound = true;
    int i = 0;

    while (i < nmeshes && notFound)
    {
      iv0 = meshes.getValue(i, 0);
      iv1 = meshes.getValue(i, 1);
      iv2 = meshes.getValue(i, 2);

      w = rayTriangleIntersect(sphPts[k], apices.getRow(iv0),
                               apices.getRow(iv1), apices.getRow(iv2));

      if (w[0] >= 0)
      {
        res[k][0] = i;
        res[k][1] = w[0];
        res[k][2] = w[1];
        res[k][3] = w[2];
        notFound = false;
      }
      i++;
    }
  }
  return res;
}

/****************************************************************************/
/*!
 **  Returns the cosine of the angular tolerance
 **
 ** \param[in]  tolang Angular tolerance
 **
 *****************************************************************************/
double GeometryHelper::getCosineAngularTolerance(double tolang)

{
  if (FFFF(tolang)) return 0.;
  if (tolang == 00.) return 1.;
  if (tolang == 90.) return 0.;
  return ABS(cos(ut_deg2rad(tolang)));
}

VectorVectorDouble GeometryHelper::getEllipse(const VectorDouble &center,
                                              double rx,
                                              double ry,
                                              double theta,
                                              int count)
{
  double angref = theta * GV_PI / 180.;
  double cosp = cos(angref);
  double sinp = sin(angref);

  VectorVectorDouble coords(2);
  coords[0].resize(count+1,0.);
  coords[1].resize(count+1,0.);

  for (int i = 0; i < count; i++)
  {
    double angle = 2. * i * GV_PI / count;
    double cosa = cos(angle);
    double sina = sin(angle);
    coords[0][i] = center[0] + rx * cosa * cosp - ry * sina * sinp;
    coords[1][i] = center[1] + rx * cosa * sinp + ry * sina * cosp;
  }

  coords[0][count] = coords[0][0];
  coords[1][count] = coords[1][0];
  return coords;
}
