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
#include "Mesh/Delaunay.hpp"
#include "Basic/Utilities.hpp"
#include "Core/io.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Mesh/AMesh.hpp"
#include "Mesh/MeshEStandard.hpp"

#include <cmath>
#include <cstdio>

/*! \cond */
#define TRIANGLES(itri, j) (triangles[(itri) * 3 + (j)] - 1)
#define MESHES(imesh, j)   (meshes[(imesh) * ncorner + (j)] - 1)
#define POINTS(ip, idim)   (points[(ip) * ndim + (idim)])
#define EXT(idim, ip)      (ext[(idim) * number + (ip)])

/*! \endcond */

namespace gstlrn
{
/****************************************************************************/
/*!
 **  Returns the shift value for the apex('icorn')
 **
 ** \return Shift value
 **
 ** \param[in]  ndim       Space dimension (1, 2 or 3)
 ** \param[in]  ipol       Rank of the polarization (starting from 0)
 ** \param[in]  icas       Rank of the case (starting from 0)
 ** \param[in]  icorn      Rank of the corner (starting from 0)
 ** \param[in]  idim       Rank of the coordinate (starting from 0)
 **
 ** \remarks This function returns the shift value to be applied to each
 ** \remarks coordinate index, for each apex for each mesh which constitutes
 ** \remarks the partition of a cell into meshes (using Turbo facility)
 **
 ** \remarks When 'ndim' is provided as negative, a special case is programmed
 **
 *****************************************************************************/
Id MSS(Id ndim, Id ipol, Id icas, Id icorn, Id idim)
{
  // Arrays S*D are provided for [ipol][icas][icorn][idim]

  constexpr Id S1D[1][1][2][1] = {{{{0}, {1}}}};
  constexpr Id S2D[2][2][3][2] = {{{{0, 0}, {1, 0}, {0, 1}},
                                    {{0, 1}, {1, 0}, {1, 1}}},
                                   {{{0, 0}, {1, 0}, {1, 1}},
                                    {{0, 0}, {0, 1}, {1, 1}}}};
  constexpr Id S3D[1][6][4][3] = {{{{0, 0, 0},
                                     {1, 0, 0},
                                     {1, 0, 1},
                                     {1, 1, 1}},
                                    {{0, 0, 0},
                                     {0, 0, 1},
                                     {1, 0, 1},
                                     {1, 1, 1}},
                                    {{0, 0, 0},
                                     {0, 0, 1},
                                     {0, 1, 1},
                                     {1, 1, 1}},
                                    {{0, 0, 0},
                                     {0, 1, 0},
                                     {0, 1, 1},
                                     {1, 1, 1}},
                                    {{0, 0, 0},
                                     {1, 0, 0},
                                     {1, 1, 0},
                                     {1, 1, 1}},
                                    {{0, 0, 0},
                                     {0, 1, 0},
                                     {1, 1, 0},
                                     {1, 1, 1}}}};

  Id ival = 0;
  if (ipol < 0 || icorn < 0 || icas < 0 || idim < 0) return ival;

  if (ndim == 1)
  {
    ival = S1D[0][0][icorn][0];
  }
  else if (ndim == 2)
  {
    ival = S2D[ipol][icas][icorn][idim];
  }
  else
  {
    ival = S3D[0][icas][icorn][idim];
  }
  return (ival);
}

/****************************************************************************/
/*!
 **  Extend a Grid by gext
 **
 ** \return The coordinates of the extension points (Dimension: number * ndim)
 **
 ** \param[in]  db         Output Db grid structure
 ** \param[in]  gext       Array of domain dilation
 **
 ** \param[out] nout       Number of extension points
 **
 ** \remarks The returned array 'ext' must be freed by the calling function
 **
 *****************************************************************************/
VectorDouble extend_grid(DbGrid* db, const VectorDouble& gext, Id* nout)
{
  Id ndim, number, ndiv, ndiv0, rank, ival, delta;

  /* Initializations */

  ndim   = db->getNDim();
  number = static_cast<Id>(pow(2., ndim));
  ndiv0  = static_cast<Id>(pow(2., ndim - 1));
  *nout  = 0;

  /* Core allocation */

  VectorInt indg(ndim, 0);
  VectorDouble coor(ndim, 0.);
  VectorDouble ext(ndim * number);

  /* Generate the corner points */

  for (Id corner = 0; corner < number; corner++)
  {
    rank = corner;
    ndiv = ndiv0;
    for (Id idim = ndim - 1; idim >= 0; idim--)
    {
      delta = static_cast<Id>((ceil)(gext[idim] / db->getDX(idim)));
      ival  = rank / ndiv;
      rank  = rank - ndiv * ival;
      ndiv /= 2;
      indg[idim] = (ival == 0) ? -delta : db->getNX(idim) + delta;
    }
    db->indicesToCoordinateInPlace(indg, coor);

    for (Id idim = 0; idim < ndim; idim++)
      EXT(idim, corner) = coor[idim];
  }

  *nout = number;
  return ext;
}

/****************************************************************************/
/*!
 **  Extend a Point Domain by gext
 **
 ** \return The coordinates of the extension points (Dimension: number * ndim)
 **
 ** \param[in]  db         Output Db grid structure
 ** \param[in]  gext       Array of domain dilation
 **
 ** \param[out] nout       Number of extension points
 **
 ** \remarks The returned array 'ext' must be freed by the calling function
 **
 *****************************************************************************/
VectorDouble extend_point(Db* db, const VectorDouble& gext, Id* nout)
{
  Id ndim   = db->getNDim();
  Id number = static_cast<Id>(pow(2., ndim));
  Id ndiv0  = static_cast<Id>(pow(2., ndim - 1));

  /* Core allocation */

  *nout = 0;
  VectorDouble coor(ndim);
  VectorDouble mini(ndim);
  VectorDouble maxi(ndim);
  VectorDouble ext(ndim * number);

  /* Calculate the extension of the domain */

  db->getExtensionInPlace(mini, maxi);

  /* Generate the corner points */

  for (Id corner = 0; corner < number; corner++)
  {
    Id rank = corner;
    Id ndiv = ndiv0;
    for (Id idim = ndim - 1; idim >= 0; idim--)
    {
      Id ival = rank / ndiv;
      rank     = rank - ndiv * ival;
      ndiv /= 2;
      coor[idim] = (ival == 0) ? mini[idim] - gext[idim] : maxi[idim] + gext[idim];
    }

    for (Id idim = 0; idim < ndim; idim++)
      EXT(idim, corner) = coor[idim];
  }

  *nout = number;
  return ext;
}

/*****************************************************************************/
/*!
 **  Define the maximum extension between dbin and/or dbout
 **
 ** \return Array of extensions
 **
 ** \param[in]  dbin      Db input (optional)
 ** \param[in]  dbout     Db output (optional)
 **
 ** \param[out] nout       Number of extension points
 **
 ** \remarks The calling function must free the returned array 'ext'
 ** \remarks with its dimension: ndim * 2^ndim
 **
 *****************************************************************************/
VectorDouble get_db_extension(Db* dbin, Db* dbout, Id* nout)
{
  Id ndim = 0;
  if (dbin != nullptr) ndim = dbin->getNDim();
  if (dbout != nullptr) ndim = dbout->getNDim();
  Id number = static_cast<Id>(pow(2., ndim));
  Id ndiv0  = static_cast<Id>(pow(2., ndim - 1));

  /* Core allocation */

  VectorDouble ext(ndim * number);

  VectorDouble coor(ndim);
  VectorDouble mini_abs;
  VectorDouble maxi_abs;
  mini_abs.resize(ndim, TEST);
  maxi_abs.resize(ndim, TEST);

  /* Get the extension */

  if (dbin != nullptr)
  {
    dbin->getExtensionInPlace(mini_abs, maxi_abs, true);
  }
  if (dbout != nullptr)
  {
    dbout->getExtensionInPlace(mini_abs, maxi_abs, true);
  }

  /* Generate the corner points */

  for (Id corner = 0; corner < number; corner++)
  {
    Id rank = corner;
    Id ndiv = ndiv0;
    for (Id idim = ndim - 1; idim >= 0; idim--)
    {
      Id ival = rank / ndiv;
      rank     = rank - ndiv * ival;
      ndiv /= 2;
      coor[idim] = (ival == 0) ? mini_abs[idim] : maxi_abs[idim];
    }
    for (Id idim = 0; idim < ndim; idim++)
      EXT(idim, corner) = coor[idim];
  }

  *nout = number;
  return ext;
}

/*****************************************************************************/
/*!
 **  Load the vertices in a segment and check if the segment is masked
 **
 ** \return 1 If the segment is valid because at least one vertex is active
 **
 ** \param[in] dbgrid      Db structure
 ** \param[in] ipos        Position of newly created mesh information
 ** \param[in] ix1         Grid index along X for the vertex #1
 ** \param[in] ix2         Grid index along X for the vertex #2
 **
 ** \param[out] mesh       Array of triangle ranks (dimension = 3)
 ** \param[out] order      Array of relative ranks
 **
 ** \remarks The values in 'order' are the absolute indices (starting from 1),
 ** \remarks negative if the grid node is masked off
 **
 *****************************************************************************/
static Id st_load_segment(DbGrid* dbgrid,
                           VectorInt& mesh,
                           VectorInt& order,
                           Id ipos,
                           Id ix1,
                           Id ix2)
{
  Id iech1, iech2, imask1, imask2;
  VectorInt indg(1);

  Id nactive = 0;

  indg[0]        = ix1;
  iech1          = dbgrid->indiceToRank(indg);
  mesh[ipos + 0] = iech1;
  imask1         = dbgrid->isActive(iech1);
  nactive += imask1;

  indg[0]        = ix2;
  iech2          = dbgrid->indiceToRank(indg);
  mesh[ipos + 1] = iech2;
  imask2         = dbgrid->isActive(iech2);
  nactive += imask2;

  if (nactive <= 0) return (0);

  order[iech1] = (imask1) ? 1 : -1;
  order[iech2] = (imask2) ? 1 : -1;

  return (1);
}

/*****************************************************************************/
/*!
 **  Load the vertices in a triangle and check if the triangle is masked
 **
 ** \return 1 If the triangle is valid because at least one vertex is active
 **
 ** \param[in] dbgrid      Db structure
 ** \param[in] ipos        Position of newly created mesh information
 ** \param[in] ix1         Grid index along X for the vertex #1
 ** \param[in] iy1         Grid index along Y for the vertex #1
 ** \param[in] ix2         Grid index along X for the vertex #2
 ** \param[in] iy2         Grid index along Y for the vertex #2
 ** \param[in] ix3         Grid index along X for the vertex #3
 ** \param[in] iy3         Grid index along Y for the vertex #3
 **
 ** \param[out] mesh       Array of triangle ranks (dimension = 3)
 ** \param[out] order      Array of relative ranks
 **
 ** \remarks The values in 'order' are the absolute indices (starting from 1),
 ** \remarks negative if the grid node is masked off
 **
 *****************************************************************************/
static Id st_load_triangle(DbGrid* dbgrid,
                            VectorInt& mesh,
                            VectorInt& order,
                            Id ipos,
                            Id ix1,
                            Id iy1,
                            Id ix2,
                            Id iy2,
                            Id ix3,
                            Id iy3)
{
  Id iech1, iech2, iech3, imask1, imask2, imask3;
  VectorInt indg(2);

  Id nactive = 0;

  indg[0]        = ix1;
  indg[1]        = iy1;
  iech1          = dbgrid->indiceToRank(indg);
  mesh[ipos + 0] = iech1;
  imask1         = dbgrid->isActive(iech1);
  nactive += imask1;

  indg[0]        = ix2;
  indg[1]        = iy2;
  iech2          = dbgrid->indiceToRank(indg);
  mesh[ipos + 1] = iech2;
  imask2         = dbgrid->isActive(iech2);
  nactive += imask2;

  indg[0]        = ix3;
  indg[1]        = iy3;
  iech3          = dbgrid->indiceToRank(indg);
  mesh[ipos + 2] = iech3;
  imask3         = dbgrid->isActive(iech3);
  nactive += imask3;

  if (nactive <= 0) return (0);

  order[iech1] = (imask1) ? 1 : -1;
  order[iech2] = (imask2) ? 1 : -1;
  order[iech3] = (imask3) ? 1 : -1;

  return (1);
}

/*****************************************************************************/
/*!
 **  Load the vertices in a tetrahedron and check if the tetrahedron is masked
 **
 ** \return 1 If the tetrahedron is valid because at least one vertex is active
 **
 ** \param[in] dbgrid      Db structure
 ** \param[in] mesh        Array of triangle ranks (dimension = 4)
 ** \param[in] order       Array of relative ranks
 ** \param[in] ipos        Position of newly created mesh information
 ** \param[in] ix1         Grid index along X for the vertex #1
 ** \param[in] iy1         Grid index along Y for the vertex #1
 ** \param[in] iz1         Grid index along Z for the vertex #1
 ** \param[in] ix2         Grid index along X for the vertex #2
 ** \param[in] iy2         Grid index along Y for the vertex #2
 ** \param[in] iz2         Grid index along Z for the vertex #2
 ** \param[in] ix3         Grid index along X for the vertex #3
 ** \param[in] iy3         Grid index along Y for the vertex #3
 ** \param[in] iz3         Grid index along Z for the vertex #3
 ** \param[in] ix4         Grid index along X for the vertex #4
 ** \param[in] iy4         Grid index along Y for the vertex #4
 ** \param[in] iz4         Grid index along Z for the vertex #4
 **
 ** \remarks The values in 'order' are the absolute indices (starting from 1),
 ** \remarks negative if the grid node is masked off
 **
 *****************************************************************************/
static Id st_load_tetra(DbGrid* dbgrid,
                         VectorInt& mesh,
                         VectorInt& order,
                         Id ipos,
                         Id ix1,
                         Id iy1,
                         Id iz1,
                         Id ix2,
                         Id iy2,
                         Id iz2,
                         Id ix3,
                         Id iy3,
                         Id iz3,
                         Id ix4,
                         Id iy4,
                         Id iz4)
{
  Id iech1, iech2, iech3, iech4, imask1, imask2, imask3, imask4;
  VectorInt indg(3);

  Id nactive = 0;

  indg[0]        = ix1;
  indg[1]        = iy1;
  indg[2]        = iz1;
  iech1          = dbgrid->indiceToRank(indg);
  mesh[ipos + 0] = iech1;
  imask1         = dbgrid->isActive(iech1);
  nactive += imask1;

  indg[0]        = ix2;
  indg[1]        = iy2;
  indg[2]        = iz2;
  iech2          = dbgrid->indiceToRank(indg);
  mesh[ipos + 1] = iech2;
  imask2         = dbgrid->isActive(iech2);
  nactive += imask1;

  indg[0]        = ix3;
  indg[1]        = iy3;
  indg[2]        = iz3;
  iech3          = dbgrid->indiceToRank(indg);
  mesh[ipos + 2] = iech3;
  imask3         = dbgrid->isActive(iech3);
  nactive += imask1;

  indg[0]        = ix4;
  indg[1]        = iy4;
  indg[2]        = iz4;
  iech4          = dbgrid->indiceToRank(indg);
  mesh[ipos + 3] = iech4;
  imask4         = dbgrid->isActive(iech4);
  nactive += imask1;

  if (nactive <= 0) return (0);

  order[iech1] = (imask1) ? 1 : -1;
  order[iech2] = (imask2) ? 1 : -1;
  order[iech3] = (imask3) ? 1 : -1;
  order[iech4] = (imask4) ? 1 : -1;

  return (1);
}

/*****************************************************************************/
/*!
 **  Perform the ultimate task for the regular grid
 **
 ** \return Pointer to the newly created MeshEStandard structure
 **
 ** \param[in] dbgrid   Db structure
 ** \param[in] ndim     Space dimension
 ** \param[in] nmesh    Current number of meshes
 ** \param[in] ncorner  Number of corner per mesh
 ** \param[in] meshes   Initial array of meshes
 ** \param[in] order    Array of relative ranks
 **
 *****************************************************************************/
static MeshEStandard* st_ultimate_regular_grid(Db* dbgrid,
                                               Id ndim,
                                               Id nmesh,
                                               Id ncorner,
                                               VectorInt& meshes,
                                               VectorInt& order)
{
  /* Count the number of active vertices */

  Id number  = dbgrid->getNSample();
  Id nvertex = 0;
  Id nin     = 0;
  for (Id iech = 0; iech < number; iech++)
  {
    Id local = order[iech];
    if (IFFFF(local)) continue;
    if (local > 0) nin++;
    nvertex++;
  }

  /* Define the addresses assigned to information */

  Id rank_in  = 0;
  Id rank_out = nin;
  VectorInt ranks(number);
  for (Id iech = 0; iech < number; iech++)
  {
    Id local = order[iech];
    if (IFFFF(local)) continue;
    if (local > 0)
      ranks[iech] = rank_in++;
    else
      ranks[iech] = rank_out++;
  }

  /* Store the (active) vertices */

  VectorDouble points(nvertex * ndim, 0);
  for (Id iech = 0; iech < number; iech++)
  {
    Id local = order[iech];
    if (IFFFF(local)) continue;
    Id jech = ranks[iech];
    for (Id idim = 0; idim < ndim; idim++)
      points[jech * ndim + idim] = dbgrid->getCoordinate(iech, idim);
  }

  // Update the point ranks in the mesh */

  for (Id i = 0; i < nmesh * ncorner; i++)
    meshes[i] = ranks[meshes[i]];

  // Store the information in the returned AMesh structure

  auto* amesh = new MeshEStandard();
  amesh->reset(ndim, ncorner, points, meshes, false);

  return amesh;
}

/*****************************************************************************/
/*!
 **  Build the regular 2-D grid meshing
 **
 ** \return Pointer to the newly created AMesh structure
 **
 ** \param[in]  dbgrid    Db structure
 **
 *****************************************************************************/
AMesh* meshes_turbo_2D_grid_build(DbGrid* dbgrid)
{
  Id ndim    = 2;
  Id ncorner = 3;
  Id nx      = dbgrid->getNX(0);
  Id ny      = dbgrid->getNX(1);
  Id number  = nx * ny;

  /* Core allocation */

  VectorInt meshes(number * ncorner * 2, 0);
  VectorInt order(number, ITEST);

  /* Store the indices per mesh */

  Id nmesh = 0;
  for (Id ix = 0; ix < nx - 1; ix++)
    for (Id iy = 0; iy < ny - 1; iy++)
    {
      Id ipol = ((ix + iy) % 2 == 1) ? 0 : 1;
      for (Id i = 0; i < 2; i++)
        if (st_load_triangle(dbgrid, meshes, order, nmesh * ncorner,
                             ix + MSS(2, ipol, i, 0, 0),
                             iy + MSS(2, ipol, i, 0, 1),
                             ix + MSS(2, ipol, i, 1, 0),
                             iy + MSS(2, ipol, i, 1, 1),
                             ix + MSS(2, ipol, i, 2, 0),
                             iy + MSS(2, ipol, i, 2, 1))) nmesh++;
    }

  // Shrink the array (optional)

  meshes.resize(nmesh * ncorner);

  // Perform the ultimate tasks

  AMesh* amesh = st_ultimate_regular_grid(dbgrid, ndim, nmesh, ncorner, meshes, order);

  return amesh;
}

/*****************************************************************************/
/*!
 **  Dump the contents of a triangulation in an ASCII file according to the
 **  STL format
 **
 ** \param[in]  file_name Nmae of the created ASCII file
 ** \param[in]  obj_name  Name assigned to the object
 ** \param[in]  verbose   Verbose flag
 ** \param[in]  ndim      Space dimension
 ** \param[in]  ncode     Number of different codes
 ** \param[in]  ntri      Number of triangles (expected)
 ** \param[in]  npoints   Number of poins (expected)
 ** \param[in]  ntcode    Array of number of triangles per code
 ** \param[in]  triangles Array of vertex indices for each triangle
 ** \param[in]  points    Array of 3-D coordinates for triangle vertices
 **
 *****************************************************************************/
Id meshes_2D_write(const char* file_name,
                    const char* obj_name,
                    Id verbose,
                    Id ndim,
                    Id ncode,
                    Id ntri,
                    Id npoints,
                    const VectorInt& ntcode,
                    const VectorInt& triangles,
                    const VectorDouble& points)
{
  FILE* file;
  Id i, itri, ntriloc;
  double normal[3];

  /* Opening the ASCII File */

  file = _file_open(file_name, 1);
  if (file == nullptr)
  {
    messerr("Error when opening the file %s", file_name);
    return (1);
  }

  /* Print the statistics */

  if (verbose)
  {
    message("Number of codes     = %d\n", ncode);
    message("Number of triangles = %d\n", ntri);
    message("Number of vertices  = %d\n", npoints);
  }

  /* Write the Object name */

  fprintf(file, "solid %s\n", obj_name);

  /* Define the normal */

  normal[0] = -1.;
  normal[1] = 0.;
  normal[2] = 0.;

  /* Loop on the different codes */

  itri = 0;
  for (Id icode = 0; icode < ncode; icode++)
  {
    ntriloc = ntcode[icode];

    /* Loop on the triangles */

    for (Id jtri = 0; jtri < ntriloc; jtri++, itri++)
    {
      fprintf(file, " facet normal %lf %lf %lf\n", normal[0], normal[1], normal[2]);
      fprintf(file, "   outer loop\n");

      /* Loop on the vertices */

      for (Id icorn = 0; icorn < 3; icorn++)
      {
        i = TRIANGLES(itri, icorn);
        fprintf(file, "    vertex ");

        /* Loop on the space dimension */

        for (Id idim = 0; idim < ndim; idim++)
          fprintf(file, "%lf ", POINTS(i, idim));
        fprintf(file, "\n");
      }
      fprintf(file, "  endloop\n");
      fprintf(file, " endfacet\n");
    }
  }
  fprintf(file, "endsolid %s\n", obj_name);

  fclose(file);
  return (0);
}

/*****************************************************************************/
/*!
 **  Calculate the statistics on a set of meshes
 **
 ** \param[in]  ndim      Space dimension
 ** \param[in]  ncorner   Number of corners for each mesh
 ** \param[in]  nmesh     Number of meshes
 ** \param[in]  meshes    Array of vertex indices for each mesh
 ** \param[in]  points    Array of 'ndim' coordinates for mesh vertex
 **
 *****************************************************************************/
void mesh_stats(Id ndim, Id ncorner, Id nmesh, const I32* meshes, const double* points)
{
  VectorDouble mini(ndim, 0.);
  VectorDouble maxi(ndim, 0.);
  for (Id i = 0; i < ndim; i++)
  {
    mini[i] = MAXIMUM_BIG;
    maxi[i] = MINIMUM_BIG;
  }
  Id imin = 10000000;
  Id imax = -1;

  /* Loop on the meshes */

  for (Id imesh = 0; imesh < nmesh; imesh++)
  {
    for (Id icorn = 0; icorn < ncorner; icorn++)
    {
      Id i = MESHES(imesh, icorn);
      if (i < imin) imin = i;
      if (i > imax) imax = i;
      for (Id idim = 0; idim < ndim; idim++)
      {
        double value = POINTS(i, idim);
        if (value < mini[idim]) mini[idim] = value;
        if (value > maxi[idim]) maxi[idim] = value;
      }
    }
  }

  /* Print the statistics */

  bool flag_print = true;
  for (Id idim = 0; idim < ndim; idim++)
    if (mini[idim] > maxi[idim]) flag_print = false;
  if (imin > imax) flag_print = false;

  if (flag_print)
  {
    message("Statistics on the Meshes:\n");
    message("- Apex rank: from %d to %d\n", imin, imax);
    for (Id idim = 0; idim < ndim; idim++)
      message("- Coord#%d: from %lf to %lf\n", idim + 1, mini[idim], maxi[idim]);
  }
}

/*****************************************************************************/
/*!
 **  Build the regular 3-D grid meshing
 **
 ** \return Pointer to the newly created AMesh structure
 **
 ** \param[in]  dbgrid    Db structure
 **
 *****************************************************************************/
AMesh* meshes_turbo_3D_grid_build(DbGrid* dbgrid)
{
  Id ndim    = 3;
  Id ncorner = 4;
  Id nx      = dbgrid->getNX(0);
  Id ny      = dbgrid->getNX(1);
  Id nz      = dbgrid->getNX(2);
  Id number  = nx * ny * nz;

  /* Core allocation */

  VectorInt meshes(number * ncorner * 6, 0);
  VectorInt order(number, ITEST);

  /* Store the indices per mesh */

  Id nmesh = 0;
  for (Id ix = 0; ix < nx - 1; ix++)
    for (Id iy = 0; iy < ny - 1; iy++)
      for (Id iz = 0; iz < nz - 1; iz++)
        for (Id i = 0; i < 6; i++)
        {
          if (st_load_tetra(dbgrid, meshes, order, nmesh * ncorner,
                            ix + MSS(3, 0, i, 0, 0), iy + MSS(3, 0, i, 0, 1),
                            iz + MSS(3, 0, i, 0, 2), ix + MSS(3, 0, i, 1, 0),
                            iy + MSS(3, 0, i, 1, 1), iz + MSS(3, 0, i, 1, 2),
                            ix + MSS(3, 0, i, 2, 0), iy + MSS(3, 0, i, 2, 1),
                            iz + MSS(3, 0, i, 2, 2), ix + MSS(3, 0, i, 3, 0),
                            iy + MSS(3, 0, i, 3, 1), iz + MSS(3, 0, i, 3, 2)))
            nmesh++;
        }

  // Shrink the array (optional)

  meshes.resize(nmesh * ncorner);

  // Perform the ultimate tasks

  AMesh* amesh = st_ultimate_regular_grid(dbgrid, ndim, nmesh, ncorner, meshes, order);

  return amesh;
}

/*****************************************************************************/
/*!
 **  Build the regular meshing from a 1-D grid
 **
 ** \return The newly created AMesh structure
 **
 ** \param[in]  dbgrid    Db structure
 **
 *****************************************************************************/
AMesh* meshes_turbo_1D_grid_build(DbGrid* dbgrid)
{
  Id ndim    = 1;
  Id ncorner = 2;
  Id nx      = dbgrid->getNX(0);
  Id number  = nx - 1;

  /* Core allocation */

  VectorInt meshes(number * ncorner, 0);
  VectorInt order(number, ITEST);

  /* Store the indices per mesh */

  Id nmesh = 0;
  for (Id ix = 0; ix < nx - 1; ix++)
  {
    if (st_load_segment(dbgrid, meshes, order, nmesh * ncorner,
                        ix + MSS(1, 0, 1, 0, 0), ix + MSS(1, 0, 1, 1, 0)))
      nmesh++;
  }

  // Shrink the array (optional)

  meshes.resize(nmesh * ncorner);

  // Perform the ultimate tasks

  AMesh* amesh = st_ultimate_regular_grid(dbgrid, ndim, nmesh, ncorner, meshes, order);

  return amesh;
}

} // namespace gstlrn
