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
#include "gmtsph.hpp"

#include "Basic/Law.hpp"
#include "Basic/MathFunc.hpp"
#include "Basic/VectorNumT.hpp"
#include "Core/Keypair.hpp"
#include "Db/Db.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "Mesh/AMesh.hpp"
#include "Mesh/LinkSphTriangle.hpp"

#include <cmath>

// External library

#define COORD(i, j) (coord[3 * (j) + (i)])

#define DEBUG 0

namespace gstlrn
{
/*****************************************************************************/
/*!
 **  Initialize the structure for triangles on a sphere
 **
 ** \param[in]  t      Pointer to the SphTriangle structure to be initialized
 **
 *****************************************************************************/
void meshes_2D_sph_init(SphTriangle* t)
{
  t->n_nodes  = 0;
  t->sph_size = 0;
  t->sph_x.clear();
  t->sph_y.clear();
  t->sph_z.clear();
  t->sph_list.clear();
  t->sph_lptr.clear();
  t->sph_lend.clear();
}

/*****************************************************************************/
/*!
 **  Free the structure for triangles on a sphere
 **
 ** \param[in]  t      Pointer to the SphTriangle structure to be freed
 ** \param[in]  mode   1 for partial deallocation
 **                    0 for total deallocation
 **
 *****************************************************************************/
void meshes_2D_sph_free(SphTriangle* t, Id mode)
{
  if (t == nullptr) return;
  if (mode == 0)
  {
    t->sph_x.clear();
    t->sph_y.clear();
    t->sph_z.clear();
    t->n_nodes = 0;
  }
  t->sph_list.clear();
  t->sph_lptr.clear();
  t->sph_lend.clear();
  t->sph_size = 0;
}

/*****************************************************************************/
/*!
 **  Define the memory areas by reading information from Db
 **
 ** \param[in]  db        Db structure where data are located
 ** \param[in]  t         Pointer to the SphTriangle structure to be loaded
 **
 ** \remarks This function adds vertices to an existing SphTriangle structure
 ** \remarks A conversion is launched to convert the 2-D information
 ** \remarks (longitude,latitude) into 3-D coordinates
 **
 *****************************************************************************/
Id meshes_2D_sph_from_db(Db* db, SphTriangle* t)
{
  Id error, nech, ndim, neff, nold, ecr;
  double xx, yy, zz;

  /* Initializations */

  if (db == nullptr) return (0);
  error = 1;
  nech  = db->getNSample();
  ndim  = db->getNDim();
  if (ndim != 2)
  {
    messerr(
      "In Spherical System, the Space Dimension of the data base Db must be 2 (%d)\n",
      ndim);
    return (1);
  }

  /* Count the number of active samples */

  neff = db->getNSample(true);

  /* Core allocation */

  nold = t->n_nodes;
  ecr  = nold;
  t->sph_x.resize(nold + neff);
  t->sph_y.resize(nold + neff);
  t->sph_z.resize(nold + neff);

  /* Load the points */

  for (Id iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    GH::convertSph2Cart(db->getCoordinate(iech, 0), db->getCoordinate(iech, 1),
                        &xx, &yy, &zz);
    t->sph_x[ecr] = xx;
    t->sph_y[ecr] = yy;
    t->sph_z[ecr] = zz;
    ecr++;
  }
  t->n_nodes = static_cast<I32>(nold + neff);

  /* Set the error return code */

  error = 0;

  if (error) meshes_2D_sph_free(t, 0);
  return (error);
}

/*****************************************************************************/
/*!
 **  Add fixed points to modify the SphTriangle structure
 **
 ** \param[in]  nech      Number of added samples
 ** \param[in]  x,y       Array containing the coordinates of added points
 ** \param[in]  t         Pointer to the SphTriangle structure to be loaded
 **
 ** \remarks This function adds the vertices to an existing
 ** \remarks SphTriangle structure
 **
 *****************************************************************************/
Id meshes_2D_sph_from_points(Id nech, double* x, double* y, SphTriangle* t)
{
  Id error, ecr, nold;
  double xx, yy, zz;

  /* Initializations */

  error = 1;

  /* List of points */

  nold = t->n_nodes;
  ecr  = nold;
  t->sph_x.resize(nold + nech);
  t->sph_y.resize(nold + nech);
  t->sph_z.resize(nold + nech);

  /* Load the information */

  for (Id iech = 0; iech < nech; iech++)
  {
    GH::convertSph2Cart(x[iech], y[iech], &xx, &yy, &zz);
    t->sph_x[ecr] = xx;
    t->sph_y[ecr] = yy;
    t->sph_z[ecr] = zz;
    ecr++;
  }
  t->n_nodes = static_cast<I32>(nold + nech);

  /* Set the error return code */

  error = 0;

  if (error) meshes_2D_sph_free(t, 0);
  return (error);
}

/*****************************************************************************/
/*!
 **  Add auxiliary random points
 **
 ** \param[in]  triswitch  Triangulation option
 ** \param[in]  t          SphTriangle structure
 **
 ** \remarks This function adds the vertices to an existing SphTriangle structure
 **
 *****************************************************************************/
Id meshes_2D_sph_from_auxiliary(const String& triswitch, SphTriangle* t)
{
  Id error, npoint, ecr, found_close, nech, nold, ndecode, flag_reg, flag_vdc;
  double c1[3], c2[3], dist;
  static double eps = 1.e-3;

  /* Initializations */

  error = 1;
  VectorDouble coord;
  nold    = t->n_nodes;
  ndecode = flag_vdc = flag_reg = npoint = 0;

  // We set the random seed (for reproductible exemples)

  law_set_random_seed(43241);

  /* Decode the 'triswitch' criterion */

  if (triswitch[0] == '-' && triswitch[1] == 'n')
  {
    flag_vdc = 1;
    ndecode  = static_cast<Id>(strtod(&triswitch[2], nullptr));
    if (ndecode <= 0) return (0);
  }
  if (triswitch[0] == '-' && triswitch[1] == 'r')
  {
    flag_reg = 1;
    ndecode  = static_cast<Id>(strtod(&triswitch[2], nullptr));
    if (ndecode <= 0) return (0);
  }
  if (triswitch[0] == '-' && triswitch[1] == 'h')
  {
    message("  usage [-n|r|h]\n");
    message("  -n  Use Van der Corput algorithm to generate N points.\n");
    message("  -r  Generate points from N iterated sphere discretization.\n");
    message("  -h  Help:  A brief instruction.\n");
    return (0);
  }

  /* Generate the random points */

  if (flag_vdc)
  {
    ut_vandercorput(ndecode, 1, 1, &npoint, coord);
  }
  else if (flag_reg)
  {
    if (ut_icosphere(ndecode, 1, &npoint, coord)) goto label_end;
  }

  /* Reallocate to maximum size */

  t->sph_x.resize(nold + npoint);
  t->sph_y.resize(nold + npoint);
  t->sph_z.resize(nold + npoint);

  /* Check that random points are not too close from hard nodes */

  ecr = nold;
  for (Id lec = 0; lec < npoint; lec++)
  {
    c1[0] = COORD(0, lec);
    c1[1] = COORD(1, lec);
    c1[2] = COORD(2, lec);

    /* Loop on the old points */

    found_close = -1;
    for (Id i = 0; i < t->n_nodes && found_close < 0; i++)
    {
      c2[0] = t->sph_x[i];
      c2[1] = t->sph_y[i];
      c2[2] = t->sph_z[i];

      /* We calculate the distance in the plane (we cannot use ut_distance */
      /* which requires long-lat information */

      dist = 0.;
      for (Id j = 0; j < 3; j++)
        dist += (c1[j] - c2[j]) * (c1[j] - c2[j]);
      if (dist < eps) found_close = i;
    }
    if (found_close >= 0) continue;
    t->sph_x[ecr] = c1[0];
    t->sph_y[ecr] = c1[1];
    t->sph_z[ecr] = c1[2];
    ecr++;
  }

  /* Final resize */

  nech = ecr;
  t->sph_x.resize(nech);
  t->sph_y.resize(nech);
  t->sph_z.resize(nech);
  t->n_nodes = static_cast<I32>(nech);

  /* Set the error return code */

  error = 0;

label_end:
  if (error) meshes_2D_sph_free(t, 0);
  return (error);
}

/*****************************************************************************/
/*!
 **  Print the contents of the SphTriangle structure
 **
 ** \param[in]  t         Pointer to the SphTriangle structure to be printed
 ** \param[in]  brief     1 for a brief output; 0 otherwise
 **
 *****************************************************************************/
void meshes_2D_sph_print(SphTriangle* t, Id brief)
{
  double rlong, rlat;

  message("- Number of nodes   = %d\n", t->n_nodes);

  if (!brief && t->n_nodes > 0 && !t->sph_x.empty() && !t->sph_y.empty() && !t->sph_z.empty())
  {
    message("\nCoordinates in Cartesian (R=1); then in Longitude - Latitude\n");
    for (Id i = 0; i < t->n_nodes; i++)
    {
      message("%3d", i + 1);
      GH::convertCart2Sph(t->sph_x[i], t->sph_y[i], t->sph_z[i], &rlong, &rlat);
      message(" Cartesian=%8.3lf %8.3lf %8.3lf - Long-Lat=%8.3lf %8.3lf\n",
              t->sph_x[i], t->sph_y[i], t->sph_z[i], rlong, rlat);
    }
    message("\n");
  }
}

/*****************************************************************************/
/*!
 **  Perform the spherical triangulation
 **
 ** \return Error return code
 **
 ** \param[in]  verbose    Verbose option
 ** \param[in]  t          SphTriangle structure
 **
 *****************************************************************************/
Id meshes_2D_sph_create(Id verbose, SphTriangle* t)
{
  std::vector<I32> loc_near, loc_next, loc_lnew;
  VectorDouble loc_dist;
  I32 error, skip_rnd, seed_memo;
  double memo[3][3], ampli, value, cste;

  /* Initializations */

  error    = 1;
  skip_rnd = static_cast<I32>(get_keypone("Skip_Random", 0));
  if (t == nullptr || t->n_nodes < 3) return (1);

  /* Re-allocate the arrays within the SphTriangle structure */

  meshes_2D_sph_free(t, 1);
  t->sph_size = 6 * t->n_nodes - 12;
  t->sph_list.resize(t->sph_size, ITEST);
  t->sph_lptr.resize(t->sph_size, ITEST);
  t->sph_lend.resize(t->sph_size, ITEST);

  /* Allocate local arrays */

  loc_lnew.resize(t->sph_size);
  loc_near.resize(t->n_nodes);
  loc_next.resize(t->n_nodes);
  loc_dist.resize(t->n_nodes);

  /* Avoid having three first points colinear */

  if (!skip_rnd)
  {
    seed_memo = static_cast<I32>(law_get_random_seed());
    law_set_random_seed(132231);
    ampli = 1.e-2;
    for (I32 i = 0; i < 3; i++)
    {
      memo[i][0] = t->sph_x[i];
      memo[i][1] = t->sph_y[i];
      memo[i][2] = t->sph_z[i];
      cste       = (t->sph_x[i] * t->sph_x[i] + t->sph_y[i] * t->sph_y[i] + t->sph_z[i] * t->sph_z[i]);
      t->sph_x[i] += law_uniform(-ampli, ampli);
      t->sph_y[i] += law_uniform(-ampli, ampli);
      value = sqrt(
        cste - t->sph_x[i] * t->sph_x[i] - t->sph_y[i] * t->sph_y[i]);
      t->sph_z[i] = (t->sph_z[i] > 0.) ? value : -value;
    }
    law_set_random_seed(seed_memo);
  }

  (void)trmesh_(&t->n_nodes, t->sph_x.data(), t->sph_y.data(), t->sph_z.data(), t->sph_list.data(),
                t->sph_lptr.data(), t->sph_lend.data(),
                loc_lnew.data(), loc_near.data(), loc_next.data(),
                loc_dist.data(), &error);

  /* Restore the initial coordinates */

  for (I32 i = 0; i < 3; i++)
  {
    t->sph_x[i] = memo[i][0];
    t->sph_y[i] = memo[i][1];
    t->sph_z[i] = memo[i][2];
  }

  switch (error)
  {
    case 0:
      if (verbose)
      {
        message("Spherical triangulation successfull\n");
        meshes_2D_sph_print(t, !DEBUG);
      }
      goto label_end;

    case -1:
      messerr("The initial number of nodes must be larger than 3");
      break;

    case -2:
      messerr("The first three nodes may not be colinear");
      break;

    default:
      meshes_2D_sph_print(t, 0);
      messerr("The node %d coincides with a node of higher rank.", error);
      break;
  }

  /* Free the result of the triangulation */

  meshes_2D_sph_free(t, 1);

label_end:
  return (error);
}

} // namespace gstlrn
