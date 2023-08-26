/*
                                      gmt_sph

Original Author: Paul Wessel
Website: https://github.com/GenericMappingTools/gmt
License: LGPL v3
*/

/*--------------------------------------------------------------------
 *
 *  Copyright (c) 2008-2023 by the GMT Team (https://www.generic-mapping-tools.org/team.html)
 *  See LICENSE.TXT file for copying and redistribution conditions.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; version 3 or any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  Contact info: www.generic-mapping-tools.org
 *--------------------------------------------------------------------*/
/*
 * Spherical triangulation - Delaunay or Voronoi options.
 * Relies on STRIPACK Fortran F77 library (Renka, 1997). Reference:
 * Renka, R, J,, 1997, Algorithm 772: STRIPACK: Delaunay Triangulation
 *    and Voronoi Diagram on the Surface of a Sphere, AMC Trans. Math.
 *    Software, 23 (3), 416-434.
 * Spherical interpolation - tension or smoothing.
 * Relies on SSRFPACK Fortran F77 library (Renka, 1997). Reference:
 * Renka, R, J,, 1997, Algorithm 773: SSRFPACK: Interpolation of
 *    Scattered Data on the Surface of a Sphere with a Surface under Tension,
 *    AMC Trans. Math. Software, 23 (3), 435-442.
 * We translated both to C using f2c and removed/rewrite statements that needed -lf2c
 *
 * Author:      Paul Wessel
 * Date:        1-AUG-2011
 * Version: API 5 64-bit
 *
 */

/*
Modified by MINES Paris / ARMINES (2023)
Authors: gstlearn Team
Website: https://gstlearn.org
*/
#pragma once

class SphTriangle
{
public:
  int n_nodes; /* Number of nodes */
  int sph_size; /* Size of arrays sph_list and sph_lptr */
  double *sph_x; /* Array of X-coordinates for nodes */
  double *sph_y; /* Array of Y-coordinates for nodes */
  double *sph_z; /* Array of Z-coordinates for nodes */
  int *sph_list; /* Set of nodal indexes */
  int *sph_lptr; /* Set of pointers (sph_list indexes) */
  int *sph_lend; /* Set of pointers to adjacency lists */
};

int trmesh_(int *n,
            double *x,
            double *y,
            double *z__,
            int *list,
            int *lptr,
            int *lend,
            int *lnew,
            int *near__,
            int *next,
            double *dist,
            int *ier);
int trlist_(int *n,
            int *list,
            int *lptr,
            int *lend,
            int *nrow,
            int *nt,
            int *ltri,
            int *ier);

