/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Basic/ERotation.hpp"

GSTLEARN_EXPORT void ut_rotation_init(int ndim, double *rot);
GSTLEARN_EXPORT void ut_rotation_sincos(double angle,
                                        double *cosa,
                                        double *sina);
GSTLEARN_EXPORT void ut_rotation_matrix_2D(double angle, double *rot);
GSTLEARN_EXPORT void ut_rotation_matrix_3D(double alpha,
                                           double beta,
                                           double gamma,
                                           double *rot);
GSTLEARN_EXPORT void ut_rotation_matrix(int ndim,
                                        const double *angles,
                                        double *rot);
GSTLEARN_EXPORT VectorDouble ut_rotation_matrix_VD(int ndim,
                                                   const VectorDouble &angles);
GSTLEARN_EXPORT void ut_rotation_direction(double ct, double st, double *a, double *codir);
GSTLEARN_EXPORT void ut_rotation_copy(int ndim,
                                      const double *rotin,
                                      double *rotout);
GSTLEARN_EXPORT void merge_boxes(int ndim,
                                 VectorDouble& mini1,
                                 VectorDouble& maxi1,
                                 VectorDouble& mini2,
                                 VectorDouble& maxi2);
GSTLEARN_EXPORT int ut_angles_from_rotation_matrix(const double *rot,
                                                   int ndim,
                                                   double *angles);
GSTLEARN_EXPORT void ut_angles_from_codir(int ndim,
                                          const VectorDouble &codir,
                                          VectorDouble &angles);
GSTLEARN_EXPORT void ut_angles_to_codir(int ndim,
                                        int ndir,
                                        const VectorDouble &angles,
                                        VectorDouble &codir);
GSTLEARN_EXPORT int ut_rotation_check(double *rot, int ndim);
GSTLEARN_EXPORT double distance_point_to_segment(double x0,
                                                 double y0,
                                                 double x1,
                                                 double y1,
                                                 double x2,
                                                 double y2,
                                                 double *xd,
                                                 double *yd,
                                                 int *nint);
GSTLEARN_EXPORT double ut_geodetic_angular_distance(double long1,
                                                    double lat1,
                                                    double long2,
                                                    double lat2,
                                                    double radius = 1.);
GSTLEARN_EXPORT void ut_geodetic_angles(double long1,
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
                                        double *C);
GSTLEARN_EXPORT double ut_geodetic_triangle_perimeter(double long1,
                                                      double lat1,
                                                      double long2,
                                                      double lat2,
                                                      double long3,
                                                      double lat3);
GSTLEARN_EXPORT double ut_geodetic_triangle_surface(double long1,
                                                    double lat1,
                                                    double long2,
                                                    double lat2,
                                                    double long3,
                                                    double lat3);
GSTLEARN_EXPORT int segment_intersect(double xd1,
                                      double yd1,
                                      double xe1,
                                      double ye1,
                                      double xd2,
                                      double yd2,
                                      double xe2,
                                      double ye2,
                                      double *xint,
                                      double *yint);
GSTLEARN_EXPORT int is_in_spherical_triangle(double *coor,
                                             double surface,
                                             double *pts1,
                                             double *pts2,
                                             double *pts3,
                                             double *wgts);
GSTLEARN_EXPORT int is_in_spherical_triangle_optimized(const double *coor,
                                                       double *ptsa,
                                                       double *ptsb,
                                                       double *ptsc,
                                                       double *wgts);
GSTLEARN_EXPORT VectorVectorDouble util_convert_longlat(const VectorDouble& longitude,
                                                        const VectorDouble& latitude,
                                                        double dilate = 1.,
                                                        double radius_arg = 1.);
GSTLEARN_EXPORT void util_convert_cart2sph(double x,
                                           double y,
                                           double z,
                                           double *rlong,
                                           double *rlat,
                                           double radius_arg = 1.);
GSTLEARN_EXPORT void util_convert_sph2cart(double rlong,
                                           double rlat,
                                           double *x,
                                           double *y,
                                           double *z,
                                           double radius_arg = 1.);

GSTLEARN_EXPORT MatrixSquareGeneral util_gradXYToRotmat(double dzoverdx, double dzoverdy);

GSTLEARN_EXPORT VectorDouble util_rotmatToEuler(const MatrixSquareGeneral &mat,
                                                const ERotation& convrot = ERotation::SXYZ,
                                                double eps = EPSILON10);

GSTLEARN_EXPORT MatrixSquareGeneral util_EulerToRotmat(const VectorDouble& angles,
                                                       const ERotation& convrot = ERotation::SXYZ);
