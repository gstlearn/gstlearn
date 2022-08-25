/*
 * GlobalEnvironment.cpp
 *
 *  Created on: 22 juil. 2021
 *      Author: drenard
 */

#include "geoslib_f.h"
#include "Basic/Utilities.hpp"
#include "Basic/Geometry.hpp"

#include <math.h>

#define ROT(i,j)     (rot[(i) * ndim + (j)])

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
void ut_rotation_sincos(double angle, double *cosa, double *sina)
{
  double value;

  if (angle == 0.)
  {
    *cosa = 1.;
    *sina = 0.;
    return;
  }
  else if (angle == 90.)
  {
    *cosa = 0.;
    *sina = 1.;
    return;
  }
  else if (angle == 180.)
  {
    *cosa = -1.;
    *sina = 0.;
    return;
  }
  else if (angle == 270.)
  {
    *cosa = 0.;
    *sina = -1.;
    return;
  }
  else
  {
    value = ut_deg2rad(angle);
    *cosa = cos(value);
    *sina = sin(value);
  }
  return;
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
void ut_rotation_matrix_2D(double angle, double *rot)
{
  double ca, sa;

  ut_rotation_sincos(angle, &ca, &sa);

  /* Define the 2-D rotation matrix */

  rot[0] = ca;
  rot[1] = sa;
  rot[2] = -sa;
  rot[3] = ca;

  return;
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
void ut_rotation_matrix_3D(double alpha, double beta, double gamma, double *rot)
{
  double ca[3], sa[3];

  /* Initializations */

  ut_rotation_sincos(alpha, &ca[0], &sa[0]);
  ut_rotation_sincos(beta, &ca[1], &sa[1]);
  ut_rotation_sincos(gamma, &ca[2], &sa[2]);

  /* Define the 3-D rotation matrix */

  rot[0] = ca[0] * ca[1];
  rot[3] = -sa[0] * ca[2] + ca[0] * sa[1] * sa[2];
  rot[6] = sa[0] * sa[2] + ca[0] * sa[1] * ca[2];
  rot[1] = sa[0] * ca[1];
  rot[4] = ca[0] * ca[2] + sa[0] * sa[1] * sa[2];
  rot[7] = -ca[0] * sa[2] + sa[0] * sa[1] * ca[2];
  rot[2] = -sa[1];
  rot[5] = ca[1] * sa[2];
  rot[8] = ca[1] * ca[2];

  return;
}

/*****************************************************************************/
/*!
 **  Calculates the rotation matrix
 **
 ** \param[in]  ndim   Space dimension
 ** \param[in]  angles Array of angles
 **
 ** \param[out] rot   direct rotation matrix (Dimension = 9)
 **
 *****************************************************************************/
void ut_rotation_matrix(int ndim, const double *angles, double *rot)
{
  if (ndim == 2)
    ut_rotation_matrix_2D(angles[0], rot);
  else if (ndim == 3)
    ut_rotation_matrix_3D(angles[0], angles[1], angles[2], rot);
  else
    ut_rotation_init(ndim, rot);
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
VectorDouble ut_rotation_matrix_VD(int ndim, const VectorDouble &angles)
{
  VectorDouble rot;

  rot.resize(ndim * ndim);
  if (ndim == 2)
    ut_rotation_matrix_2D(angles[0], rot.data());
  else if (ndim == 3)
    ut_rotation_matrix_3D(angles[0], angles[1], angles[2], rot.data());
  else
    ut_rotation_init(ndim, rot.data());

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
void ut_rotation_copy(int ndim, const double *rotin, double *rotout)
{
  int i;

  for (i = 0; i < ndim * ndim; i++)
    rotout[i] = rotin[i];
}

/****************************************************************************/
/*!
 **  Merge the extensions of the boxes (parallel to main axes)
 **
 ** \param[in]  ndim     Space dimension
 ** \param[in]  mini1    Input array containing the minimum along each axis
 ** \param[in]  maxi1    Input array containing the maximum along each axis
 ** \param[in]  mini2    Output array containing the minimum along each axis
 ** \param[in]  maxi2    Output array containing the maximum along each axis
 **
 *****************************************************************************/
void merge_boxes(int ndim,
                 VectorDouble& mini1,
                 VectorDouble& maxi1,
                 VectorDouble& mini2,
                 VectorDouble& maxi2)
{
  for (int idim = 0; idim < ndim; idim++)
  {
    double mini = 1.e30;
    if (! mini1.empty()) mini = MIN(mini, mini1[idim]);
    if (! mini2.empty()) mini = MIN(mini, mini2[idim]);

    double maxi = -1.e30;
    if (! maxi1.empty()) maxi = MAX(maxi, maxi1[idim]);
    if (! maxi2.empty()) maxi = MAX(maxi, maxi2[idim]);
  }
}

/****************************************************************************/
/*!
 **  Calculates the rotation angles from the rotation matrix
 **
 ** \param[in]  rot   Rotation matrix (Dimension = 9)
 ** \param[in]  ndim  Space dimension
 **
 ** \param[out]  angles Rotation angles (Dimension = ndim)
 **
 *****************************************************************************/
int ut_angles_from_rotation_matrix(const double *rot, int ndim, double *angles)
{
  double s0, c0, s1, c1, s2, c2;
  int i, nval;

  /* Initializations */

  for (i = 0; i < ndim; i++)
    angles[i] = 0.;
  if (rot == nullptr) return (0);

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
    s1 = -rot[2];
    c1 = sqrt(rot[0] * rot[0] + rot[1] * rot[1]);
    if (ABS(c1) < EPSILON10)
    {
      if (s1 > 0.)
      {
        angles[0] = 0.;
        angles[1] = GV_PI / 2.;
        angles[2] = atan2(rot[3], rot[6]);
      }
      else
      {
        angles[0] = 0.;
        angles[1] = -GV_PI / 2.;
        angles[2] = atan2(-rot[3], -rot[6]);
      }
    }
    else
    {
      c0 = rot[0] / c1;
      s0 = rot[1] / c1;
      s2 = rot[5] / c1;
      c2 = rot[8] / c1;
      angles[0] = atan2(s0, c0);
      angles[1] = atan2(s1, c1);
      angles[2] = atan2(s2, c2);
    }
  }
  else
  {
    nval = 0;
  }

  /* Convert into degrees */

  for (i = 0; i < nval; i++)
    angles[i] = ut_rad2deg(angles[i]);

  return (nval);
}

/****************************************************************************/
/*!
 **  Calculates the rotation angle from the direction coefficient
 **
 ** \param[in]  ndim   Space dimension
 ** \param[in]  codir  Direction vector (Dimension = ndim)
 **
 ** \param[out]  angles Rotation angles (Dimension = ndim)
 **
 *****************************************************************************/
void ut_angles_from_codir(int ndim,
                          const VectorDouble &codir,
                          VectorDouble &angles)
{
  double norme;
  int i, nval;

  /* Initializations */

  for (i = 0; i < ndim; i++)
    angles[i] = 0.;

  /* Dispatch */

  if (ndim == 1)
  {
    return;
  }
  else if (ndim == 2)
  {
    angles[0] = atan2(codir[1], codir[0]);
    angles[1] = 0.;
    nval = 1;
  }
  else if (ndim == 3)
  {
    norme = codir[0] * codir[0] + codir[1] * codir[1];
    if (norme > 0.)
    {
      norme = sqrt(norme);
      angles[0] = atan2(codir[1] / norme, codir[0] / norme);
      angles[1] = atan2(codir[2], norme);
    }
    nval = 2;
  }
  else
  {
    nval = 0;
  }

  /* Convert into degrees */

  for (i = 0; i < nval; i++)
    angles[i] = ut_rad2deg(angles[i]);

  return;
}

/****************************************************************************/
/*!
 **   Convert angles to rotation matrices
 **
 ** \param[in]  ndim     Number of space dimensions
 ** \param[in]  ndir     Number of directions
 ** \param[in]  angles   Vector giving the angles characteristics (in degrees)
 **
 ** \param[out] codir    Vector of the direction (Dim: ndir * ndim)
 **
 ** \remarks If angles is not provided:
 ** \remarks - if ndir == ndim: return basic direction of space
 ** \remarks - if ndim < ndim: return 'ndir' directions regular in the 2-D
 **
 *****************************************************************************/
void ut_angles_to_codir(int ndim,
                        int ndir,
                        const VectorDouble &angles,
                        VectorDouble &codir)
{
  if (ndim <= 1) return;

  codir.resize(ndim * ndir);
  for (int i = 0; i < ndim * ndir; i++)
    codir[i] = 0.;

  if (angles.size() <= 0)
  {
    if (ndir == ndim)
    {
      int ecr = 0;
      for (int idir = 0; idir < ndir; idir++)
        for (int idim = 0; idim < ndim; idim++)
          codir[ecr++] = (idir == idim) ? 1. :
                                          0.;
    }
    else
    {
      {
        for (int idir = 0; idir < ndir; idir++)
        {
          double angref = 180. * idir / ndir;
          codir[idir * ndim + 0] = cos(angref * GV_PI / 180.);
          codir[idir * ndim + 1] = sin(angref * GV_PI / 180.);
        }
      }
    }
  }
  else
  {
    if ((int) angles.size() == ndir)
    {
      for (int idir = 0; idir < ndir; idir++)
      {
        codir[idir * ndim + 0] = cos(angles[idir] * GV_PI / 180.);
        codir[idir * ndim + 1] = sin(angles[idir] * GV_PI / 180.);
      }
    }
  }
  return;
}

/****************************************************************************/
/*!
 **  Starting from a rotation matrix, check it is different from the Identity
 **
 ** \return  1 if a rotation is defined; 0 otherwise
 **
 ** \param[in]  rot      Rotation matrix
 ** \param[in]  ndim     Space dimension
 **
 *****************************************************************************/
int ut_rotation_check(double *rot, int ndim)
{
  int i, j;

  for (i = 0; i < ndim; i++)
    for (j = 0; j < ndim; j++)
    {
      if (i == j)
      {
        if (ABS(ROT(i,j) - 1.) > EPSILON10) return (1);
      }
      else
      {
        if (ABS(ROT(i,j)) > EPSILON10) return (1);
      }
    }
  return (0);
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
double distance_point_to_segment(double x0,
                                 double y0,
                                 double x1,
                                 double y1,
                                 double x2,
                                 double y2,
                                 double *xd,
                                 double *yd,
                                 int *nint)
{
  double dx, dy, dxp, dyp, ratio, dist, signe;

  dx = x2 - x1;
  dy = y2 - y1;

  ratio = (dx * (x0 - x1) + dy * (y0 - y1)) / (dx * dx + dy * dy);

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

  dxp = x0 - (*xd);
  dyp = y0 - (*yd);
  dist = sqrt(dxp * dxp + dyp * dyp);
  signe = (dy * (x0 - x1) - dx * (y0 - y1) > 0.) ? 1 : -1;

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
double ut_geodetic_angular_distance(double long1,
                                    double lat1,
                                    double long2,
                                    double lat2,
                                    double radius)
{
  double rlon1, rlat1, rlon2, rlat2, dlong, angdst;

  rlon1 = ut_deg2rad(long1);
  rlat1 = ut_deg2rad(lat1);
  rlon2 = ut_deg2rad(long2);
  rlat2 = ut_deg2rad(lat2);
  dlong = rlon2 - rlon1;
  angdst = acos(sin(rlat1) * sin(rlat2) + cos(rlat1) * cos(rlat2) * cos(dlong));
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
  double prod, cosA;

  prod = sinb * sinc;
  cosA = (prod == 0.) ? 0. : (cosa - cosb * cosc) / prod;
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
void ut_geodetic_angles(double long1,
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
  double cosa, cosb, cosc, sina, sinb, sinc;

  *a = ut_geodetic_angular_distance(long2, lat2, long3, lat3);
  *b = ut_geodetic_angular_distance(long1, lat1, long3, lat3);
  *c = ut_geodetic_angular_distance(long1, lat1, long2, lat2);

  cosa = cos(*a);
  cosb = cos(*b);
  cosc = cos(*c);
  sina = sin(*a);
  sinb = sin(*b);
  sinc = sin(*c);

  *A = st_convert_geodetic_angle(sina, cosa, sinb, cosb, sinc, cosc);
  *B = st_convert_geodetic_angle(sinb, cosb, sinc, cosc, sina, cosa);
  *C = st_convert_geodetic_angle(sinc, cosc, sina, cosa, sinb, cosb);

  return;
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
double ut_geodetic_triangle_perimeter(double long1,
                                      double lat1,
                                      double long2,
                                      double lat2,
                                      double long3,
                                      double lat3)
{
  double a, b, c, ga, gb, gc, perimeter;

  ut_geodetic_angles(long1, lat1, long2, lat2, long3, lat3, &a, &b, &c, &ga,
                     &gb, &gc);
  perimeter = a + b + c;
  return (perimeter);
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
double ut_geodetic_triangle_surface(double long1,
                                    double lat1,
                                    double long2,
                                    double lat2,
                                    double long3,
                                    double lat3)
{
  double a, b, c, A, B, C, surface;

  ut_geodetic_angles(long1, lat1, long2, lat2, long3, lat3, &a, &b, &c, &A, &B,
                     &C);
  surface = (A + B + C - GV_PI);
  return (surface);
}

/****************************************************************************/
/*!
 **  Is a point inside a spherical triangle
 **
 ** \return 1 if the point belongs to the spherical triangle; 0 otherwise
 **
 ** \param[in]  coor    Coordinates of the target point (long,lat)
 ** \param[in]  ptsa    Coordinates of the first point of the triangle
 ** \param[in]  ptsb    Coordinates of the second point of the triangle
 ** \param[in]  ptsc    Coordinates of the third point of the triangle
 **
 ** \param[out] wgts    Array of weights
 **
 *****************************************************************************/
int is_in_spherical_triangle_optimized(const double *coor,
                                       double *ptsa,
                                       double *ptsb,
                                       double *ptsc,
                                       double *wgts)
{
  double total, s[3], stot, eps;
  double A, B, C, AB, AC, BA, BC, CA, CB, OA, OB, OC;
  double dab, dbc, dac, d0a, d0b, d0c;
  double sinab, cosab, sinbc, cosbc, sinac, cosac;
  double sin0a, cos0a, cos0b, sin0b, cos0c, sin0c;

  eps = 1.e-6;
  total = 0.;

  dab = ut_geodetic_angular_distance(ptsa[0], ptsa[1], ptsb[0], ptsb[1]);
  dbc = ut_geodetic_angular_distance(ptsb[0], ptsb[1], ptsc[0], ptsc[1]);
  dac = ut_geodetic_angular_distance(ptsa[0], ptsa[1], ptsc[0], ptsc[1]);
  d0a = ut_geodetic_angular_distance(coor[0], coor[1], ptsa[0], ptsa[1]);
  d0b = ut_geodetic_angular_distance(coor[0], coor[1], ptsb[0], ptsb[1]);
  d0c = ut_geodetic_angular_distance(coor[0], coor[1], ptsc[0], ptsc[1]);

  sinab = sin(dab);
  cosab = cos(dab);
  sinbc = sin(dbc);
  cosbc = cos(dbc);
  sinac = sin(dac);
  cosac = cos(dac);
  sin0a = sin(d0a);
  cos0a = cos(d0a);
  sin0b = sin(d0b);
  cos0b = cos(d0b);
  sin0c = sin(d0c);
  cos0c = cos(d0c);

  A = st_convert_geodetic_angle(sinbc, cosbc, sinac, cosac, sinab, cosab);
  B = st_convert_geodetic_angle(sinac, cosac, sinab, cosab, sinbc, cosbc);
  C = st_convert_geodetic_angle(sinab, cosab, sinbc, cosbc, sinac, cosac);
  stot = (A + B + C - GV_PI);

  OA = st_convert_geodetic_angle(sinbc, cosbc, sin0c, cos0c, sin0b, cos0b);
  BA = st_convert_geodetic_angle(sin0c, cos0c, sin0b, cos0b, sinbc, cosbc);
  CA = st_convert_geodetic_angle(sin0b, cos0b, sinbc, cosbc, sin0c, cos0c);
  s[0] = (OA + BA + CA - GV_PI);
  total += s[0];
  if (total > stot + eps) return (0);

  AB = st_convert_geodetic_angle(sin0c, cos0c, sinac, cosac, sin0a, cos0a);
  OB = st_convert_geodetic_angle(sinac, cosac, sin0a, cos0a, sin0c, cos0c);
  CB = st_convert_geodetic_angle(sin0a, cos0a, sin0c, cos0c, sinac, cosac);
  s[1] = (AB + OB + CB - GV_PI);
  total += s[1];
  if (total > stot + eps) return (0);

  AC = st_convert_geodetic_angle(sin0b, cos0b, sin0a, cos0a, sinab, cosab);
  BC = st_convert_geodetic_angle(sin0a, cos0a, sinab, cosab, sin0b, cos0b);
  OC = st_convert_geodetic_angle(sinab, cosab, sin0b, cos0b, sin0a, cos0a);
  s[2] = (AC + BC + OC - GV_PI);
  total += s[2];
  if (ABS(total - stot) > eps) return (0);

  for (int i = 0; i < 3; i++)
    wgts[i] = s[i] / total;
  return (1);
}

/****************************************************************************/
/*!
 **  Calculate the intersection between two segments
 **
 ** \return 0 there is an intersection; 1 if there is no intersection
 **
 ** \param[in]  xd1,yd1     Starting point for the first segment
 ** \param[in]  xe1,ye1     Ending point for the first segment
 ** \param[in]  xd2,yd2     Starting point for the second segment
 ** \param[in]  xe2,ye2     Ending point for the second segment
 **
 ** \param[out]   xint,yint  Coordinates of the intersection
 **
 *****************************************************************************/
int segment_intersect(double xd1,
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
  double a1, a2, b1, b2, x, y, x1m, x1M, x2m, x2M, testval;

  /* Preliminary check */

  b1 = ye1 - yd1;
  b2 = ye2 - yd2;

  /* Case of two horizontal segments */

  if (ABS(b1) < EPSILON10 && ABS(b2) < EPSILON10)
  {
    if (ABS(ye1 - ye2) > EPSILON10) return (1);
    x1m = MIN(xd1, xe1);
    x1M = MAX(xd1, xe1);
    x2m = MIN(xd2, xe2);
    x2M = MAX(xd2, xe2);
    if (x1m > x2M || x2m > x1M) return (1);
    (*xint) = MAX(x1m, x2m);
    (*yint) = ye1;
    return (0);
  }

  /* Case of the horizontal first segment */

  if (ABS(b1) < EPSILON10)
  {
    y = ye1;
    x = xe2 + (y - ye2) * (xe2 - xd2) / b2;
    if ((x - xd1) * (x - xe1) > 0) return (1);
    if ((y - yd1) * (y - ye1) > 0) return (1);
    if ((x - xd2) * (x - xe2) > 0) return (1);
    if ((y - yd2) * (y - ye2) > 0) return (1);
    (*xint) = x;
    (*yint) = y;
    return (0);
  }

  /* Case of horizontal second segment */

  if (ABS(b2) < EPSILON10)
  {
    y = ye2;
    x = xe1 + (y - ye1) * (xe1 - xd1) / b1;
    if ((x - xd1) * (x - xe1) > 0) return (1);
    if ((y - yd1) * (y - ye1) > 0) return (1);
    if ((x - xd2) * (x - xe2) > 0) return (1);
    if ((y - yd2) * (y - ye2) > 0) return (1);
    (*xint) = x;
    (*yint) = y;
    return (0);
  }

  /* This operation is safe as end-point ordinates cannot be equal */

  a1 = (xe1 - xd1) / b1;
  a2 = (xe2 - xd2) / b2;

  /* Skip the case of parallel fractures */

  if (ABS(a1 - a2) < EPSILON10) return (1);
  y = (xd2 - xd1 + a1 * yd1 - a2 * yd2) / (a1 - a2);

  /* Discard intersection if located outside the segment */

  if (ABS(b1) > 0)
  {
    testval = (y - yd1) * (y - ye1);
    if (testval > 0) return (1);
  }
  if (ABS(b2) > 0)
  {
    testval = (y - yd2) * (y - ye2);
    if (testval > 0) return (1);
  }

  /* Update the endpoint in case of intersection */

  (*xint) = xd1 + a1 * (y - yd1);
  (*yint) = y;
  return (0);
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
 **
 ** \param[out] wgts    Array of weights
 **
 *****************************************************************************/
int is_in_spherical_triangle(double *coor,
                             double surface,
                             double *pts1,
                             double *pts2,
                             double *pts3,
                             double *wgts)
{
  double total, s[3], eps;

  eps = 1.e-6;
  total = 0.;
  s[0] = ut_geodetic_triangle_surface(coor[0], coor[1], pts2[0], pts2[1],
                                      pts3[0], pts3[1]);
  total += s[0];
  if (total > surface + eps) return (0);
  s[1] = ut_geodetic_triangle_surface(pts1[0], pts1[1], coor[0], coor[1],
                                      pts3[0], pts3[1]);
  total += s[1];
  if (total > surface + eps) return (0);
  s[2] = ut_geodetic_triangle_surface(pts1[0], pts1[1], pts2[0], pts2[1],
                                      coor[0], coor[1]);
  total += s[2];
  if (ABS(total - surface) > eps) return (0);
  for (int i = 0; i < 3; i++)
    wgts[i] = s[i] / total;
  return (1);
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
void ut_rotation_direction(double ct, double st, double *a, double *codir)
{
  double rd, b[3], c[3], p[3];

  rd = 0.;
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
VectorVectorDouble util_convert_longlat(const VectorDouble& longitude,
                                        const VectorDouble& latitude,
                                        double dilate,
                                        double radius_arg)
{
  double radius = radius_arg;
  if (FFFF(radius)) variety_get_characteristics(&radius);
  radius *= dilate;

  VectorVectorDouble tab;
  int number = (int) longitude.size();
  if (number != (int) latitude.size()) return tab;

  // Dimension the returned argument
  tab.resize(3);
  for (int idim = 0; idim < 3; idim++)
    tab[idim].resize(number,0.);

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
void util_convert_cart2sph(double x,
                           double y,
                           double z,
                           double *rlong,
                           double *rlat,
                           double radius_arg)
{
  double radius = radius_arg;
  if (FFFF(radius)) variety_get_characteristics(&radius);

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
void util_convert_sph2cart(double rlong,
                           double rlat,
                           double *x,
                           double *y,
                           double *z,
                           double radius_arg)
{
  double radius = radius_arg;
  if (FFFF(radius)) variety_get_characteristics(&radius);

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
  vert[ndim-1] = -1.;

  // Vector orthogonal to the tilted plane
  VectorDouble vort = VectorDouble(ndim);
  vort[0] = dzoverdx;
  vort[1] = dzoverdy;
  vort[2] = -1.;
  double norme = ut_vector_norm(vort);
  for (int idim = 0; idim < ndim; idim++) vort[idim] /= norme;

  // Cross product
  VectorDouble axis(ndim,0.);
  axis[0] =  vort[1];
  axis[1] = -vort[0];
  axis[2] = 0.;

  // Norm of the cross product and dot product between ref and vort
  double normcross = sqrt(ut_vector_inner_product(axis, axis));
  double dot = ut_vector_inner_product(vert, vort);

  // Rotation angle
  double angle = atan2(normcross, dot);
  return angle;
}

/**
 * Calculate the rotation axis to translate the gradients of a tilted plane
 * (defined by the partial derivatives along X and Y) into the axis of the normal
 * @param dzoverdx Partial derivative along X
 * @param dzoverdy Partial derivative along Y
 */
VectorDouble util_rotation_gradXYToAxes(double dzoverdx, double dzoverdy)
{
  int ndim = 3;
  VectorDouble axis(ndim,0.);

  // Vector orthogonal to the horizontal plane/ i.e. the vertical axis
  VectorDouble vert = VectorDouble(3, 0.);
  vert[ndim-1] = -1.;

  // Vector orthogonal to the tilted plane
  VectorDouble vort = VectorDouble(ndim);
  vort[0] = dzoverdx;
  vort[1] = dzoverdy;
  vort[2] = -1.;
  double norme = ut_vector_norm(vort);
  for (int idim = 0; idim < ndim; idim++) vort[idim] /= norme;

  // Cross product
  axis[0] =  vort[1];
  axis[1] = -vort[0];
  axis[2] = 0.;

  // Norm of the cross product and dot product between ref and vort
  double normcross = sqrt(ut_vector_inner_product(axis, axis));

  // Rotation axis (normalized
  for (int idim = 0; idim < ndim; idim++) axis[idim] /= normcross;
  return axis;
}

MatrixSquareGeneral util_rotation_AxesAndAngleToMatrix(const VectorDouble &axis,
                                                       double angle)
{
  int ndim = 3;

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
 * @param mat Input matrix
 * @return
 *
 * @remark The code is coming from the following reference (BSD license)
 * @remark https://github.com/matthew-brett/transforms3d/blob/master/transforms3d/euler.py
 */
VectorDouble util_rotmatToEuler(const MatrixSquareGeneral &M,
                                const VectorInt &hyp,
                                double eps)
{
  VectorInt hyp_local = hyp;
  hyp_local.resize(4, 0);
  int firstaxis  = hyp_local[0];
  int parity     = hyp_local[0];
  int repetition = hyp_local[0];
  int frame      = hyp_local[0];

  VectorInt next_axis = {1, 2, 0, 1};
  int i = firstaxis;
  int j = next_axis[i+parity];
  int k = next_axis[i-parity+1];

  double ax, ay, az;

  if (repetition)
  {
    double sy = sqrt(M.getValue(i, j)*M.getValue(i, j) + M.getValue(i, k)*M.getValue(i, k));
    if (sy > eps)
    {
      ax = atan2( M.getValue(i, j),  M.getValue(i, k));
      ay = atan2( sy,                M.getValue(i, i));
      az = atan2( M.getValue(j, i), -M.getValue(k, i));
    }
    else
    {
      ax = atan2(-M.getValue(j, k),  M.getValue(j, j));
      ay = atan2( sy,                M.getValue(i, i));
      az = 0.0;
    }
  }
  else
  {
    double cy = sqrt(M.getValue(i, i)*M.getValue(i, i) + M.getValue(j, i)*M.getValue(j, i));
    if (cy > eps)
    {
      ax = atan2( M.getValue(k, j),  M.getValue(k, k));
      ay = atan2(-M.getValue(k, i),  cy);
      az = atan2( M.getValue(j, i),  M.getValue(i, i));
    }
    else
    {
      ax = atan2(-M.getValue(j, k),  M.getValue(j, j));
      ay = atan2(-M.getValue(k, i),  cy);
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

  VectorDouble angles = {ax, ay, az};
  return angles;
}

MatrixSquareGeneral util_EulerToRotmat(const VectorDouble &angles,
                                       const VectorInt &hyp)
{
  VectorInt hyp_local = hyp;
  hyp_local.resize(4, 0);
  int firstaxis  = hyp_local[0];
  int parity     = hyp_local[0];
  int repetition = hyp_local[0];
  int frame      = hyp_local[0];

  VectorInt next_axis = {1, 2, 0, 1};
  int i = firstaxis;
  int j = next_axis[i+parity];
  int k = next_axis[i-parity+1];

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

VectorInt util_Convention(const String& conv)
{
  VectorInt ret(4);

  if (conv == "sxyz") ret = {0, 0, 0, 0};
  if (conv == "sxyx") ret = {0, 0, 1, 0};
  if (conv == "sxzy") ret = {0, 1, 0, 0};
  if (conv == "sxzx") ret = {0, 1, 1, 0};
  if (conv == "syzx") ret = {1, 0, 0, 0};
  if (conv == "syzy") ret = {1, 0, 1, 0};
  if (conv == "syxz") ret = {1, 1, 0, 0};
  if (conv == "syxy") ret = {1, 1, 1, 0};
  if (conv == "szxy") ret = {2, 0, 0, 0};
  if (conv == "szxz") ret = {2, 0, 1, 0};
  if (conv == "szyx") ret = {2, 1, 0, 0};
  if (conv == "szyz") ret = {2, 1, 1, 0};
  if (conv == "rzyx") ret = {0, 0, 0, 1};
  if (conv == "rxyx") ret = {0, 0, 1, 1};
  if (conv == "ryzx") ret = {0, 1, 0, 1};
  if (conv == "rxzx") ret = {0, 1, 1, 1};
  if (conv == "rxzy") ret = {1, 0, 0, 1};
  if (conv == "ryzy") ret = {1, 0, 1, 1};
  if (conv == "rzxy") ret = {1, 1, 0, 1};
  if (conv == "ryxy") ret = {1, 1, 1, 1};
  if (conv == "ryxz") ret = {2, 0, 0, 1};
  if (conv == "rzxz") ret = {2, 0, 1, 1};
  if (conv == "rxyz") ret = {2, 1, 0, 1};
  if (conv == "rzyz") ret = {2, 1, 1, 1};
  return ret;
}
