/******************************************************************************/
/* COPYRIGHT CARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include <Geometry/GeometryHelper.hpp>
#include "geoslib_old_f.h"
#include "geoslib_f_private.h"

#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Basic/MathFunc.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Mesh/tetgen.h"
#include "Mesh/MeshEStandard.hpp"

#include <math.h>
#include <string.h>
#include <stdio.h>

/*! \cond */
#define TRIANGLES(itri,j) (triangles[(itri) * 3 + (j)] - 1)
#define MESHES(imesh,j)   (meshes[(imesh) * ncorner + (j)] - 1)
#define POINTS(ip,idim)   (points[(ip) * ndim + (idim)])
#define FAULTS(ip,idim)   (faults[nfaults * (idim) + (ip)])
#define COORD(i,j)        (coord[3 * (j) + (i)])
#define EXT(idim,ip)      (ext[(idim) * number + (ip)])

#define DEBUG 0

/*! \endcond */

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
int MSS(int ndim, int ipol, int icas, int icorn, int idim)
{
  // Arrays S*D are provided for [ipol][icas][icorn][idim]

  constexpr int S1D[1][1][2][1] = { { { { 0 }, { 1 } } } };
  constexpr int S2D[2][2][3][2] = { { { { 0, 0 }, { 1, 0 }, { 0, 1 } },
                                      { { 0, 1 }, { 1, 0 }, { 1, 1 } } },
                                    { { { 0, 0 }, { 1, 0 }, { 1, 1 } },
                                      { { 0, 0 }, { 0, 1 }, { 1, 1 } } } };
  constexpr int S3D[1][6][4][3] = { { { { 0, 0, 0 },
                                        { 1, 0, 0 },
                                        { 1, 0, 1 },
                                        { 1, 1, 1 } },
                                      { { 0, 0, 0 },
                                        { 0, 0, 1 },
                                        { 1, 0, 1 },
                                        { 1, 1, 1 } },
                                      { { 0, 0, 0 },
                                        { 0, 0, 1 },
                                        { 0, 1, 1 },
                                        { 1, 1, 1 } },
                                      { { 0, 0, 0 },
                                        { 0, 1, 0 },
                                        { 0, 1, 1 },
                                        { 1, 1, 1 } },
                                      { { 0, 0, 0 },
                                        { 1, 0, 0 },
                                        { 1, 1, 0 },
                                        { 1, 1, 1 } },
                                      { { 0, 0, 0 },
                                        { 0, 1, 0 },
                                        { 1, 1, 0 },
                                        { 1, 1, 1 } } } };

  int ival = 0;
  if (ipol < 0 || icorn < 0 || icas < 0 || idim < 0) return ival;

  if (ndim == 1)
  {
    ival = S1D[ipol][icas][icorn][idim];
  }
  else if (ndim == 2)
  {
    ival = S2D[ipol][icas][icorn][idim];
  }
  else
  {
    ival = S3D[ipol][icas][icorn][idim];
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
static double* st_extend_grid(DbGrid *db, const double *gext, int *nout)
{
  int *indg, ndim, number, ndiv, ndiv0, rank, ival, error, delta;
  double *coor, *ext;

  /* Initializations */

  error = 1;
  ndim = db->getNDim();
  number = (int) pow(2., ndim);
  ndiv0 = (int) pow(2., ndim - 1);
  indg = nullptr;
  coor = ext = nullptr;
  *nout = 0;

  /* Core allocation */

  indg = db_indg_alloc(db);
  if (indg == nullptr) goto label_end;
  coor = (double*) mem_alloc(sizeof(double) * ndim, 0);
  if (coor == nullptr) goto label_end;
  ext = (double*) mem_alloc(sizeof(double) * ndim * number, 0);
  if (ext == nullptr) goto label_end;

  /* Generate the corner points */

  for (int corner = 0; corner < number; corner++)
  {
    rank = corner;
    ndiv = ndiv0;
    for (int idim = ndim - 1; idim >= 0; idim--)
    {
      delta = static_cast<int>((ceil)(gext[idim] / db->getDX(idim)));
      ival = rank / ndiv;
      rank = rank - ndiv * ival;
      ndiv /= 2;
      indg[idim] = (ival == 0) ? -delta :
                                 db->getNX(idim) + delta;
    }
    grid_to_point(db, indg, NULL, coor);

    for (int idim = 0; idim < ndim; idim++)
      EXT(idim,corner) = coor[idim];
  }

  // Set the error returned code

  error = 0;
  *nout = number;

  label_end: indg = db_indg_free(indg);
  coor = (double*) mem_free((char* ) coor);
  if (error) ext = (double*) mem_free((char* ) ext);
  return (ext);
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
static double* st_extend_point(Db *db, const double *gext, int *nout)
{
  double *ext;

  /* Initializations */

  int ndim = db->getNDim();
  int number = (int) pow(2., ndim);
  int ndiv0 = (int) pow(2., ndim - 1);

  /* Core allocation */

  *nout = 0;
  VectorDouble coor(ndim);
  VectorDouble mini(ndim);
  VectorDouble maxi(ndim);
  ext = (double*) mem_alloc(sizeof(double) * ndim * number, 0);
  if (ext == nullptr) return ext;

  /* Calculate the extension of the domain */

  db_extension(db, mini, maxi);

  /* Generate the corner points */

  for (int corner = 0; corner < number; corner++)
  {
    int rank = corner;
    int ndiv = ndiv0;
    for (int idim = ndim - 1; idim >= 0; idim--)
    {
      int ival = rank / ndiv;
      rank = rank - ndiv * ival;
      ndiv /= 2;
      coor[idim] = (ival == 0) ? mini[idim] - gext[idim] :
                                 maxi[idim] + gext[idim];
    }

    for (int idim = 0; idim < ndim; idim++)
      EXT(idim,corner) = coor[idim];
  }

  // Set the error returned code

  *nout = number;
  return (ext);
}

/*****************************************************************************/
/*!
 **  Free the triangulateio structure
 **
 ** \param[in]  t           Pointer to the triangulateio structure to be freed
 **
 ** \remark  When mode==1, the pointers which correspond to copies are
 ** \remark  simply set to NULL (not actually freed)
 **
 *****************************************************************************/
void meshes_2D_free(triangulateio *t, int /*mode*/)
{
  if (t == (triangulateio*) NULL) return;
  t->pointlist = (double*) mem_free((char* ) t->pointlist);
  t->pointattributelist = (double*) mem_free((char* ) t->pointattributelist);
  t->pointmarkerlist = (int*) mem_free((char* ) t->pointmarkerlist);
  t->numberofpoints = 0;
  t->numberofpointattributes = 0;

  t->trianglelist = (int*) mem_free((char* ) t->trianglelist);
  t->triangleattributelist = (double*) mem_free(
      (char* ) t->triangleattributelist);
  t->trianglearealist = (double*) mem_free((char* ) t->trianglearealist);
  t->neighborlist = (int*) mem_free((char* ) t->neighborlist);
  t->numberoftriangles = 0;
  t->numberofcorners = 0;
  t->numberoftriangleattributes = 0;

  t->segmentlist = (int*) mem_free((char* ) t->segmentlist);
  t->segmentmarkerlist = (int*) mem_free((char* ) t->segmentmarkerlist);
  t->numberofsegments = 0;

  t->holelist = (double*) mem_free((char* ) t->holelist);
  t->numberofholes = 0;

  t->regionlist = (double*) mem_free((char* ) t->regionlist);
  t->numberofregions = 0;

  t->edgelist = (int*) mem_free((char* ) t->edgelist);
  t->edgemarkerlist = (int*) mem_free((char* ) t->edgemarkerlist);
  t->normlist = (double*) mem_free((char* ) t->normlist);
  t->numberofedges = 0;
}

/*****************************************************************************/
/*!
 **  Initialize the triangulateio structure
 **
 ** \param[in]  t      Pointer to the triangulateio structure to be initialized
 **
 *****************************************************************************/
void meshes_2D_init(triangulateio *t)
{
  t->pointlist = nullptr;
  t->pointattributelist = nullptr;
  t->pointmarkerlist = nullptr;
  t->numberofpoints = 0;
  t->numberofpointattributes = 0;

  t->trianglelist = nullptr;
  t->triangleattributelist = nullptr;
  t->trianglearealist = nullptr;
  t->neighborlist = nullptr;
  t->numberoftriangles = 0;
  t->numberofcorners = 0;
  t->numberoftriangleattributes = 0;

  t->segmentlist = nullptr;
  t->segmentmarkerlist = nullptr;
  t->numberofsegments = 0;

  t->holelist = nullptr;
  t->numberofholes = 0;

  t->regionlist = nullptr;
  t->numberofregions = 0;

  t->edgelist = nullptr;
  t->edgemarkerlist = nullptr;
  t->normlist = nullptr;
  t->numberofedges = 0;
}

/*****************************************************************************/
/*!
 **  Define the memory areas by reading information from Db
 **
 ** \return Error return code
 **
 ** \param[in]  db        Db structure where data are located
 ** \param[in]  use_code  1 if the code (if available) is used to distinguish
 **                       vertices and holes
 ** \param[in]  t         Pointer to the triangulateio structure to be loaded
 **
 ** \remarks This function adds the vertices (and holes) to an existing
 ** \remarks triangulateio structure
 ** \remarks Only the first two coordinates are considered
 **
 *****************************************************************************/
int meshes_2D_from_db(Db *db, int use_code, triangulateio *t)
{
  int iech, nech, error, ecr, neff, ndim, nhole, ncode, nold;

  /* Initializations */

  if (db == nullptr) return (0);
  error = 1;
  nech = db->getSampleNumber();
  ndim = db->getNDim();
  ncode = db->getCodeNumber();
  if (ndim > 2) ndim = 2;

  /* Count the number of active samples */

  neff = nhole = 0;
  for (iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    if (use_code && ncode > 0 && db->getCode(iech) < 0)
      nhole++;
    else
      neff++;
  }

  /* List of points */

  nold = t->numberofpoints;
  ecr = nold * ndim;
  t->pointlist = (double*) mem_realloc((char* ) t->pointlist,
                                       sizeof(double) * (nold + neff) * ndim,
                                       0);
  if (t->pointlist == nullptr) goto label_end;

  for (iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    if (use_code && ncode > 0 && db->getCode(iech) < 0) continue;
    for (int idim = 0; idim < ndim; idim++)
      t->pointlist[ecr++] = db->getCoordinate(iech, idim);
  }
  t->numberofpoints = nold + neff;

  /* List of holes */

  if (nhole > 0)
  {
    nold = t->numberofholes;
    ecr = nold * ndim;
    t->holelist = (double*) mem_realloc((char* ) t->holelist,
                                        sizeof(double) * (nold + nhole) * ndim,
                                        0);
    if (t->holelist == nullptr) goto label_end;
    for (iech = 0; iech < nech; iech++)
    {
      if (!db->isActive(iech)) continue;
      if (!(use_code && ncode > 0 && db->getCode(iech) < 0)) continue;
      for (int idim = 0; idim < ndim; idim++)
        t->holelist[ecr++] = db->getCoordinate(iech, idim);
    }
    t->numberofholes = nold + nhole;
  }

  /* Set the error return code */

  error = 0;

  label_end: if (error)
  {
    t->pointlist = (double*) mem_free((char* ) t->pointlist);
    t->holelist = (double*) mem_free((char* ) t->holelist);
    t->numberofpoints = 0;
    t->numberofholes = 0;
  }
  return (error);
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
static double* st_get_db_extension(Db *dbin, Db *dbout, int *nout)
{
  int ndim = 0;
  if (dbin != nullptr) ndim = dbin->getNDim();
  if (dbout != nullptr) ndim = dbout->getNDim();
  int number = (int) pow(2., ndim);
  int ndiv0 = (int) pow(2., ndim - 1);

  /* Core allocation */

  double* ext = (double*) mem_alloc(sizeof(double) * ndim * number, 0);
  if (ext == nullptr) return ext;

  VectorDouble coor(ndim);
  VectorDouble mini_abs;
  VectorDouble maxi_abs;
  mini_abs.resize(ndim,TEST);
  maxi_abs.resize(ndim,TEST);

  /* Get the extension */

  if (dbin != nullptr)
  {
    db_extension(dbin, mini_abs, maxi_abs, true);
  }
  if (dbout != nullptr)
  {
    db_extension(dbout, mini_abs, maxi_abs, true);
  }

  /* Generate the corner points */

  for (int corner = 0; corner < number; corner++)
  {
    int rank = corner;
    int ndiv = ndiv0;
    for (int idim = ndim - 1; idim >= 0; idim--)
    {
      int ival = rank / ndiv;
      rank = rank - ndiv * ival;
      ndiv /= 2;
      coor[idim] = (ival == 0) ? mini_abs[idim] : maxi_abs[idim];
    }
    for (int idim = 0; idim < ndim; idim++)
      EXT(idim,corner) = coor[idim];
  }

  *nout = number;
  return (ext);
}

/*****************************************************************************/
/*!
 **  Define the triangulation from the frame around dbin and/or dbout
 **
 ** \param[in]  dbin      Db structure for input file (optional)
 ** \param[in]  dbout     Db structure for output file (optional)
 ** \param[in]  t         Pointer to the triangulateio structure to be loaded
 **
 *****************************************************************************/
void meshes_2D_default(Db *dbin, Db *dbout, triangulateio *t)
{
  double *ext;
  int number;

  /* Initializations */

  number = 0;

  /* Get the extension of the Dbs */

  ext = st_get_db_extension(dbin, dbout, &number);

  /* Extend the triangulation */

  (void) meshes_2D_from_points(number, &EXT(0, 0), &EXT(1, 0), t);

  /* Core deallocation */

  ext = (double*) mem_free((char* ) ext);

  return;
}

/*****************************************************************************/
/*!
 **  Add fixed points to modify the triangulateio basis
 **
 ** \param[in]  nech      Number of added samples
 ** \param[in]  x,y       Array containing the coordinates of added points
 ** \param[in]  t         Pointer to the triangulateio structure to be loaded
 **
 ** \remarks This function adds the vertices (and holes) to an existing
 ** \remarks triangulateio structure
 **
 *****************************************************************************/
int meshes_2D_from_points(int nech, double *x, double *y, triangulateio *t)
{
  int iech, error, ecr, ndim, nold;

  /* Initializations */

  error = 1;
  ndim = 2;

  /* List of points */

  nold = t->numberofpoints;
  ecr = nold * ndim;
  t->pointlist = (double*) mem_realloc((char* ) t->pointlist,
                                       sizeof(double) * (nold + nech) * ndim,
                                       0);
  if (t->pointlist == nullptr) goto label_end;
  for (iech = 0; iech < nech; iech++)
  {
    t->pointlist[ecr++] = x[iech];
    t->pointlist[ecr++] = y[iech];
  }
  t->numberofpoints = nold + nech;

  /* Set the error return code */

  error = 0;

  label_end: if (error)
  {
    t->pointlist = (double*) mem_free((char* ) t->pointlist);
    t->holelist = (double*) mem_free((char* ) t->holelist);
    t->numberofpoints = 0;
    t->numberofholes = 0;
  }
  return (error);
}

/*****************************************************************************/
/*!
 **  Define the segments by reading information from memory
 **
 ** \param[in]  nseg      Number of segments
 ** \param[in]  ncol      Number of columns (2 or 3 if marker)
 ** \param[in]  segments  Array containing the segment information
 ** \param[in]  t         Pointer to the triangulateio structure to be loaded
 **
 ** \remark The intput array 'segments' is sorted by column. The segments
 ** \remark are stored in memory by row
 **
 *****************************************************************************/
int meshes_2D_from_mem(int nseg, int ncol, int *segments, triangulateio *t)
{
  int i, j, error, ndim;

  /* Initializations */

  error = 1;
  ndim = 2;

  /* List of segments */

  t->segmentlist = (int*) mem_alloc(sizeof(int) * nseg * ndim, 0);
  if (t->segmentlist == nullptr) goto label_end;
  if (ncol > ndim)
  {
    t->segmentmarkerlist = (int*) mem_alloc(sizeof(int) * nseg, 0);
    if (t->segmentmarkerlist == nullptr) goto label_end;
  }

  for (i = 0; i < nseg; i++)
  {
    for (j = 0; j < ndim; j++)
      t->segmentlist[ndim * i + j] = segments[j * nseg + i];
    if (ncol > ndim) t->segmentmarkerlist[i] = segments[2 * nseg + i];
  }
  t->numberofsegments = nseg;

  /* Set the error return code */

  error = 0;

  label_end: if (error)
  {
    t->segmentlist = (int*) mem_free((char* ) t->segmentlist);
    t->numberofsegments = 0;
  }
  return (error);
}

/*****************************************************************************/
/*!
 **  Print the contents of the triangulateio structure
 **
 ** \param[in]  t         Pointer to the triangulateio structure to be printed
 ** \param[in]  brief     1 for a brief output; 0 otherwise
 **
 *****************************************************************************/
void meshes_2D_print(triangulateio *t, int brief)
{
  int ndim, i, j, lecp, leca, lecs, lech, lect, lecta, lecn;

  ndim = 2;

  // List of vertices 

  message("- Number of vertices   = %d\n", t->numberofpoints);
  message("- Number of attributes = %d\n", t->numberofpointattributes);

  if (!brief && t->numberofpoints > 0 && t->pointlist != nullptr)
  {
    lecp = leca = 0;
    for (i = 0; i < t->numberofpoints; i++)
    {
      message("%3d", i + 1);
      for (j = 0; j < ndim; j++)
        message(" %lf", t->pointlist[lecp++]);
      if (t->pointmarkerlist != nullptr) message(" %d", t->pointmarkerlist[i]);
      for (j = 0; j < t->numberofpointattributes; j++)
        message(" %lf", t->pointattributelist[leca++]);
      message("\n");
    }
    message("\n");
  }

  // List of segments

  message("- Number of segments   = %d\n", t->numberofsegments);

  if (!brief && t->numberofsegments > 0 && t->segmentlist != nullptr)
  {
    lecs = 0;
    for (i = 0; i < t->numberofsegments; i++)
    {
      message("%3d", i + 1);
      for (j = 0; j < ndim; j++)
        message(" %d", t->segmentlist[lecs++]);
      if (t->segmentmarkerlist != nullptr)
        message(" %d", t->segmentmarkerlist[i]);
      message("\n");
    }
    message("\n");
  }

  // List of holes

  message("- Number of holes      = %d\n", t->numberofholes);

  if (!brief && t->numberofholes > 0 && t->holelist != nullptr)
  {
    lech = 0;
    for (i = 0; i < t->numberofholes; i++)
    {
      message("%3d", i + 1);
      for (j = 0; j < ndim; j++)
        message(" %lf", t->holelist[lech++]);
      message("\n");
    }
    message("\n");
  }

  // List of triangles 

  message("- Number of triangles  = %d\n", t->numberoftriangles);
  message("- Number of corners    = %d\n", t->numberofcorners);

  if (!brief && t->numberoftriangles > 0 && t->trianglelist != nullptr)
  {
    lect = lecta = lecn = 0;
    for (i = 0; i < t->numberoftriangles; i++)
    {
      message("%3d", i + 1);
      for (j = 0; j < t->numberofcorners; j++)
        message(" %d", t->trianglelist[lect++]);
      if (t->triangleattributelist != nullptr)
        for (j = 0; j < t->numberoftriangleattributes; j++)
          message("%lf", t->triangleattributelist[lecta++]);
      if (t->trianglearealist != nullptr)
        message(" %lf", t->trianglearealist[i]);
      if (t->neighborlist != nullptr) for (j = 0; j < 3; j++)
        message(" %d", t->neighborlist[lecn++]);
      message("\n");
    }
    message("\n");
  }
}

/*****************************************************************************/
/*!
 **  Load the contents of the triangulateio structure into arrays
 **
 ** \param[in]  t      Pointer to the triangulateio structure to be downloaded
 **
 *****************************************************************************/
MeshEStandard* meshes_2D_load_vertices(triangulateio *t)
{
  MeshEStandard* amesh = new MeshEStandard();
  amesh->reset(2, 3, t->numberofpoints, t->numberofsegments,
                       t->pointlist, t->segmentlist);
  return amesh;
}

/*****************************************************************************/
/*!
 **  Load the vertices in a segment and check if the segment is masked
 **
 ** \return 1 If the segment is valid because at least one vertex is active
 **
 ** \param[in] dbgrid      Db structure
 ** \param[in] ix1         Grid index along X for the vertex #1
 ** \param[in] ix2         Grid index along X for the vertex #2
 **
 ** \param[out] mesh       Array of triangle ranks (dimension = 3)
 ** \param[out] order      Array of relative ranks
 ** \param[out] indg       Array for grid indexation
 **
 ** \remarks The values in 'order' are the absolute indices (starting from 1),
 ** \remarks negative if the grid node is masked off
 **
 *****************************************************************************/
static int st_load_segment(DbGrid *dbgrid,
                           VectorInt& mesh,
                           VectorInt& order,
                           int ipos,
                           int ix1,
                           int ix2)
{
  int nactive, iech1, iech2, imask1, imask2;
  VectorInt indg(1);

  nactive = 0;

  indg[0] = ix1;
  iech1 = dbgrid->indiceToRank(indg);
  mesh[ipos + 0] = iech1;
  imask1 = dbgrid->isActive(iech1);
  nactive += imask1;

  indg[0] = ix2;
  iech2 = dbgrid->indiceToRank(indg);
  mesh[ipos + 1] = iech2;
  imask2 = dbgrid->isActive(iech2);
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
static int st_load_triangle(DbGrid *dbgrid,
                            VectorInt& mesh,
                            VectorInt& order,
                            int ipos,
                            int ix1,
                            int iy1,
                            int ix2,
                            int iy2,
                            int ix3,
                            int iy3)
{
  int nactive, iech1, iech2, iech3, imask1, imask2, imask3;
  VectorInt indg(2);

  nactive = 0;

  indg[0] = ix1;
  indg[1] = iy1;
  iech1 = dbgrid->indiceToRank(indg);
  mesh[ipos + 0] = iech1;
  imask1 = dbgrid->isActive(iech1);
  nactive += imask1;

  indg[0] = ix2;
  indg[1] = iy2;
  iech2 = dbgrid->indiceToRank(indg);
  mesh[ipos + 1] = iech2;
  imask2 = dbgrid->isActive(iech2);
  nactive += imask2;

  indg[0] = ix3;
  indg[1] = iy3;
  iech3 = dbgrid->indiceToRank(indg);
  mesh[ipos + 2] = iech3;
  imask3 = dbgrid->isActive(iech3);
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
static int st_load_tetra(DbGrid *dbgrid,
                         VectorInt& mesh,
                         VectorInt& order,
                         int ipos,
                         int ix1,
                         int iy1,
                         int iz1,
                         int ix2,
                         int iy2,
                         int iz2,
                         int ix3,
                         int iy3,
                         int iz3,
                         int ix4,
                         int iy4,
                         int iz4)
{
  int nactive, iech1, iech2, iech3, iech4, imask1, imask2, imask3, imask4;
  VectorInt indg(3);

  nactive = 0;

  indg[0] = ix1;
  indg[1] = iy1;
  indg[2] = iz1;
  iech1 = dbgrid->indiceToRank(indg);
  mesh[ipos + 0] = iech1;
  imask1 = dbgrid->isActive(iech1);
  nactive += imask1;

  indg[0] = ix2;
  indg[1] = iy2;
  indg[2] = iz2;
  iech2 = dbgrid->indiceToRank(indg);
  mesh[ipos + 1] = iech2;
  imask2 = dbgrid->isActive(iech2);
  nactive += imask1;

  indg[0] = ix3;
  indg[1] = iy3;
  indg[2] = iz3;
  iech3 = dbgrid->indiceToRank(indg);
  mesh[ipos + 2] = iech3;
  imask3 = dbgrid->isActive(iech3);
  nactive += imask1;

  indg[0] = ix4;
  indg[1] = iy4;
  indg[2] = iz4;
  iech4 = dbgrid->indiceToRank(indg);
  mesh[ipos + 3] = iech4;
  imask4 = dbgrid->isActive(iech4);
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
static MeshEStandard* st_ultimate_regular_grid(Db *dbgrid,
                                               int ndim,
                                               int nmesh,
                                               int ncorner,
                                               VectorInt& meshes,
                                               VectorInt& order)
{
  /* Count the number of active vertices */

  int number = dbgrid->getSampleNumber();
  int nvertex = 0;
  int nin = 0;
  for (int iech = 0; iech < number; iech++)
  {
    int local = order[iech];
    if (IFFFF(local)) continue;
    if (local > 0) nin++;
    nvertex++;
  }

  /* Define the addresses assigned to information */

  int rank_in = 0;
  int rank_out = nin;
  VectorInt ranks(number);
  for (int iech = 0; iech < number; iech++)
  {
    int local = order[iech];
    if (IFFFF(local)) continue;
    if (local > 0)
      ranks[iech] = rank_in++;
    else
      ranks[iech] = rank_out++;
  }

  /* Store the (active) vertices */

  VectorDouble points(nvertex * ndim, 0);
  for (int iech = 0; iech < number; iech++)
  {
    int local = order[iech];
    if (IFFFF(local)) continue;
    int jech = ranks[iech];
    for (int idim = 0; idim < ndim; idim++)
      points[jech * ndim + idim] = dbgrid->getCoordinate(iech, idim);
  }

  // Update the point ranks in the mesh */

  for (int i = 0; i < nmesh * ncorner; i++)
    meshes[i] = ranks[meshes[i]];

  // Store the information in the returned AMesh structure

  MeshEStandard* amesh = new MeshEStandard();
  amesh->reset(ndim, ncorner, points, meshes, false);

  return amesh;
}

/*****************************************************************************/
/*!
 **  Build the regular 2-D grid meshing
 **
 ** \return Pointer to the newly created AMesh structure
 **
 ** \param[in]  verbose   Verbose flag
 ** \param[in]  dbgrid    Db structure
 **
 *****************************************************************************/
AMesh* meshes_turbo_2D_grid_build(int verbose, DbGrid *dbgrid)
{
  int ndim = 2;
  int ncorner = 3;
  int nx = dbgrid->getNX(0);
  int ny = dbgrid->getNX(1);
  int number = nx * ny;

  /* Core allocation */

  VectorInt meshes(number * ncorner * 2, 0);
  VectorInt order(number, ITEST);

  /* Store the indices per mesh */

  int nmesh = 0;
  for (int ix = 0; ix < nx - 1; ix++)
    for (int iy = 0; iy < ny - 1; iy++)
    {
      int ipol = ((ix + iy) % 2 == 1) ? 0 : 1;
      for (int i = 0; i < 2; i++)
        if (st_load_triangle(dbgrid, meshes, order, nmesh*ncorner,
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
 **  Strip the triangles as soon as one of its vertex is intercepted by
 **  a fault segment
 **
 ** \param[in]  t          triangulateio structure
 ** \param[in]  verbose    Verbose flag
 ** \param[in]  nfaults    Number of faults
 ** \param[in]  faults     Array giving the fault segments
 **
 *****************************************************************************/
static void st_strip_triangles_intercepted_faults(triangulateio *t,
                                                  int verbose,
                                                  int nfaults,
                                                  double *faults)
{
  int *meshes, *ranks, nmesh, ncorner, ndim, skip, i1, i2, j1, jcorn, ecr, lec;
  double *points, xx, yy;

  /* Initializations */

  if (nfaults < 2) return;
  nmesh = t->numberoftriangles;
  meshes = t->trianglelist;
  points = t->pointlist;
  ncorner = 3;
  ndim = 2;
  skip = 0;

  /* Core allocation */

  ranks = (int*) mem_alloc(sizeof(int) * nmesh, 1);

  /* Loop on the triangles */

  for (int imesh = 0; imesh < nmesh; imesh++)
  {
    skip = 0;

    /* Loop on the vertices */

    for (int icorn = 0; icorn < ncorner && !skip; icorn++)
    {
      i1 = MESHES(imesh, icorn);
      jcorn = (icorn + 1 == ncorner) ? 0 :
                                       icorn + 1;
      i2 = MESHES(imesh, jcorn);

      /* Loop on the faults vertices */

      for (int j2 = 1; j2 < nfaults && !skip; j2++)
      {
        j1 = j2 - 1;
        skip += (GH::segmentIntersect(POINTS(i1, 0), POINTS(i1, 1),
                                      POINTS(i2, 0), POINTS(i2, 1),
                                      FAULTS(j1, 0), FAULTS(j1, 1),
                                      FAULTS(j2, 0), FAULTS(j2, 1), &xx, &yy)
                 != 1);
      }
    }
    ranks[imesh] = skip;
  }

  /* Suppress the discarded triangles */

  ecr = 0;
  for (lec = 0; lec < nmesh; lec++)
  {
    if (ranks[lec]) continue;
    for (int icorn = 0; icorn < ncorner; icorn++)
      meshes[ecr * ncorner + icorn] = meshes[lec * ncorner + icorn];
    ecr++;
  }

  if (ecr < nmesh)
  {
    if (verbose)
    {
      message("Discarding Triangles due to Fault interception:\n");
      message("- Number of fault segments   = %d\n", nfaults);
      message("- Initial count of triangles = %d\n", nmesh);
      message("- Final count of triangles   = %d\n", ecr);
    }

    t->trianglelist = (int*) mem_realloc((char* ) t->trianglelist,
                                         sizeof(int) * ncorner * ecr, 1);
    t->numberoftriangles = ecr;
  }

  ranks = (int*) mem_free((char* ) ranks);
}

/*****************************************************************************/
/*!
 **  Perform the triangulation
 **
 ** \param[in]  verbose    Verbose option
 ** \param[in]  triswitch  Triangulation option
 ** \param[in]  in         triangulateio structure for input
 ** \param[in]  out        triangulateio structure for output
 ** \param[in]  vorout     triangulateio structure for voronoi
 **
 *****************************************************************************/
void meshes_2D_create(int verbose,
                      const String &triswitch,
                      triangulateio *in,
                      triangulateio *out,
                      triangulateio *vorout)
{
  int ndim, ncorner, ncol, nrow;
  double *faults;

  ndim = 2;
  ncorner = 3;
  faults = nullptr;
  triangulate(triswitch.c_str(), in, out, vorout);

  /* Strip off the triangles which intercept faults */

  if (!get_keypair("Intercept_Faults", &nrow, &ncol, &faults))
  {
    if (ncol == 2)
      st_strip_triangles_intercepted_faults(out, verbose, nrow, faults);
    faults = (double*) mem_free((char* ) faults);
  }

  /* Optional printout */

  if (verbose)
  {
    message("Triangulation performed (%s):\n", triswitch.c_str());
    meshes_2D_print(out, ! DEBUG);
    mesh_stats(ndim, ncorner, out->numberoftriangles, out->trianglelist,
               out->pointlist);
  }
}

/****************************************************************************/
/*!
 **  Create the extended domain
 **
 ** \param[in]  dbout      Output Db grid structure
 ** \param[in]  gext       Array of domain dilation
 ** \param[in]  t          Triangulation environment
 **
 *****************************************************************************/
void meshes_2D_extended_domain(Db *dbout,
                               const double *gext,
                               triangulateio *t)
{
  int number, flag_extend;
  double *ext;

  /* Initializations */

  if (dbout == nullptr) return;
  if (gext == nullptr) return;
  ext = nullptr;

  /* Check that the extension parameters are correctly defined */

  flag_extend = 0;
  for (int idim = 0; idim < dbout->getNDim(); idim++)
    if (gext[idim] > 0) flag_extend = 1;
  if (!flag_extend) return;

  /* Dilate the grid */

  if (dbout->isGrid())
  {
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(dbout);
    ext = st_extend_grid(dbgrid, gext, &number);
  }
  else
    ext = st_extend_point(dbout, gext, &number);
  if (meshes_2D_from_points(number, &EXT(0, 0), &EXT(1, 0), t)) goto label_end;

  /* Core deallocation */

  label_end: ext = (double*) mem_free((char* ) ext);
  return;
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
int meshes_2D_write(const char *file_name,
                    const char *obj_name,
                    int verbose,
                    int ndim,
                    int ncode,
                    int ntri,
                    int npoints,
                    int *ntcode,
                    int *triangles,
                    double *points)
{
  FILE *file;
  int i, itri, ntriloc;
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
  for (int icode = 0; icode < ncode; icode++)
  {
    ntriloc = ntcode[icode];

    /* Loop on the triangles */

    for (int jtri = 0; jtri < ntriloc; jtri++, itri++)
    {
      fprintf(file, " facet normal %lf %lf %lf\n", normal[0], normal[1],
              normal[2]);
      fprintf(file, "   outer loop\n");

      /* Loop on the vertices */

      for (int icorn = 0; icorn < 3; icorn++)
      {
        i = TRIANGLES(itri, icorn);
        fprintf(file, "    vertex ");

        /* Loop on the space dimension */

        for (int idim = 0; idim < ndim; idim++)
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
void mesh_stats(int ndim, int ncorner, int nmesh, int *meshes, double *points)
{
  double *mini, *maxi, value;
  int flag_print, imin, imax;

  /* Core allocation */

  mini = (double*) mem_alloc(sizeof(double) * ndim, 1);
  maxi = (double*) mem_alloc(sizeof(double) * ndim, 1);
  for (int i = 0; i < ndim; i++)
  {
    mini[i] = +1.e30;
    maxi[i] = -1.e30;
  }
  imin = 10000000;
  imax = -1;

  /* Loop on the meshes */

  for (int imesh = 0; imesh < nmesh; imesh++)
  {
    for (int icorn = 0; icorn < ncorner; icorn++)
    {
      int i = MESHES(imesh, icorn);
      if (i < imin) imin = i;
      if (i > imax) imax = i;
      for (int idim = 0; idim < ndim; idim++)
      {
        value = POINTS(i, idim);
        if (value < mini[idim]) mini[idim] = value;
        if (value > maxi[idim]) maxi[idim] = value;
      }
    }
  }

  /* Print the statistics */

  flag_print = 1;
  for (int idim = 0; idim < ndim; idim++)
    if (mini[idim] > maxi[idim]) flag_print = 0;
  if (imin > imax) flag_print = 0;

  if (flag_print)
  {
    message("Statistics on the Meshes:\n");
    message("- Apex rk: from %d to %d\n", imin, imax);
    for (int idim = 0; idim < ndim; idim++)
      message("- Coord#%d: from %lf to %lf\n", idim + 1, mini[idim],
              maxi[idim]);
  }

  /* Core deallocation */

  mini = (double*) mem_free((char* ) mini);
  maxi = (double*) mem_free((char* ) maxi);
  return;
}

/*****************************************************************************/
/*!
 **  Initialize the structure for triangles on a sphere
 **
 ** \param[in]  t      Pointer to the SphTriangle structure to be initialized
 **
 *****************************************************************************/
void meshes_2D_sph_init(SphTriangle *t)
{
  t->n_nodes = 0;
  t->sph_size = 0;
  t->sph_x = nullptr;
  t->sph_y = nullptr;
  t->sph_z = nullptr;
  t->sph_list = nullptr;
  t->sph_lptr = nullptr;
  t->sph_lend = nullptr;
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
void meshes_2D_sph_free(SphTriangle *t, int mode)
{
  if (t == (SphTriangle*) NULL) return;
  if (mode == 0)
  {
    t->sph_x = (double*) mem_free((char* ) t->sph_x);
    t->sph_y = (double*) mem_free((char* ) t->sph_y);
    t->sph_z = (double*) mem_free((char* ) t->sph_z);
    t->n_nodes = 0;
  }
  t->sph_list = (int*) mem_free((char* ) t->sph_list);
  t->sph_lptr = (int*) mem_free((char* ) t->sph_lptr);
  t->sph_lend = (int*) mem_free((char* ) t->sph_lend);
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
int meshes_2D_sph_from_db(Db *db, SphTriangle *t)
{
  int error, nech, ndim, neff, nold, ecr;
  double xx, yy, zz;

  /* Initializations */

  if (db == nullptr) return (0);
  error = 1;
  nech = db->getSampleNumber();
  ndim = db->getNDim();
  if (ndim != 2)
  {
    messerr(
        "In Spherical System, the Space Dimension of the data base Db must be 2 (%d)\n",
        ndim);
    return (1);
  }

  /* Count the number of active samples */

  neff = db->getActiveSampleNumber();

  /* Core allocation */

  nold = t->n_nodes;
  ecr = nold;
  t->sph_x = (double*) mem_realloc((char* ) t->sph_x,
                                   sizeof(double) * (nold + neff), 0);
  if (t->sph_x == nullptr) goto label_end;
  t->sph_y = (double*) mem_realloc((char* ) t->sph_y,
                                   sizeof(double) * (nold + neff), 0);
  if (t->sph_y == nullptr) goto label_end;
  t->sph_z = (double*) mem_realloc((char* ) t->sph_z,
                                   sizeof(double) * (nold + neff), 0);
  if (t->sph_z == nullptr) goto label_end;

  /* Load the points */

  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    GH::convertSph2Cart(db->getCoordinate(iech, 0), db->getCoordinate(iech, 1),
                        &xx, &yy, &zz);
    t->sph_x[ecr] = xx;
    t->sph_y[ecr] = yy;
    t->sph_z[ecr] = zz;
    ecr++;
  }
  t->n_nodes = nold + neff;

  /* Set the error return code */

  error = 0;

  label_end: if (error) meshes_2D_sph_free(t, 0);
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
int meshes_2D_sph_from_points(int nech,
                                              double *x,
                                              double *y,
                                              SphTriangle *t)
{
  int error, ecr, nold;
  double xx, yy, zz;

  /* Initializations */

  error = 1;

  /* List of points */

  nold = t->n_nodes;
  ecr = nold;
  t->sph_x = (double*) mem_realloc((char* ) t->sph_x,
                                   sizeof(double) * (nold + nech), 0);
  if (t->sph_x == nullptr) goto label_end;
  t->sph_y = (double*) mem_realloc((char* ) t->sph_y,
                                   sizeof(double) * (nold + nech), 0);
  if (t->sph_y == nullptr) goto label_end;
  t->sph_z = (double*) mem_realloc((char* ) t->sph_z,
                                   sizeof(double) * (nold + nech), 0);
  if (t->sph_z == nullptr) goto label_end;

  /* Load the information */

  for (int iech = 0; iech < nech; iech++)
  {
    GH::convertSph2Cart(x[iech], y[iech], &xx, &yy, &zz);
    t->sph_x[ecr] = xx;
    t->sph_y[ecr] = yy;
    t->sph_z[ecr] = zz;
    ecr++;
  }
  t->n_nodes = nold + nech;

  /* Set the error return code */

  error = 0;

  label_end: if (error) meshes_2D_sph_free(t, 0);
  return (error);
}

/*****************************************************************************/
/*!
 **  Add auxiliary random points
 **
 ** \param[in]  triswitch  Triangulation option
 ** \param[in]  t          SphTriangle structure
 **
 ** \remarks This function adds the vertices to an existing
 ** \remarks SphTriangle structure
 **
 *****************************************************************************/
int meshes_2D_sph_from_auxiliary(const String &triswitch,
                                                 SphTriangle *t)
{
  int error, npoint, ecr, found_close, nech, nold, ndecode, flag_reg, flag_vdc;
  double *coord, c1[3], c2[3], dist;
  static double eps = 1.e-3;

  /* Initializations */

  error = 1;
  coord = nullptr;
  nold = t->n_nodes;
  ndecode = flag_vdc = flag_reg = npoint = 0;

  // We set the random seed (for reproductible exemples)

  law_set_random_seed(43241);

  /* Decode the 'triswitch' criterion */

  if (triswitch[0] == '-' && triswitch[1] == 'n')
  {
    flag_vdc = 1;
    ndecode = (int) strtod(&triswitch[2], nullptr);
    if (ndecode <= 0) return (0);
  }
  if (triswitch[0] == '-' && triswitch[1] == 'r')
  {
    flag_reg = 1;
    ndecode = (int) strtod(&triswitch[2], nullptr);
    if (ndecode <= 0) return (0);
  }
  if (triswitch[0] == '-' && triswitch[1] == 'h')
  {
    message("  usage [-nrh]\n");
    message("  -n  Use Van der Corput algorithm to generate N points.\n");
    message("  -r  Generate points from N iterated sphere discretization.\n");
    message("  -h  Help:  A brief instruction.\n");
    return (0);
  }

  /* Generate the random points */

  if (flag_vdc)
  {
    ut_vandercorput(ndecode, 1, 1, &npoint, &coord);
  }
  else if (flag_reg)
  {
    if (ut_icosphere(ndecode, 1, &npoint, &coord)) goto label_end;
  }

  /* Reallocate to maximum size */

  t->sph_x = (double*) mem_realloc((char* ) t->sph_x,
                                   sizeof(double) * (nold + npoint), 0);
  if (t->sph_x == nullptr) goto label_end;
  t->sph_y = (double*) mem_realloc((char* ) t->sph_y,
                                   sizeof(double) * (nold + npoint), 0);
  if (t->sph_y == nullptr) goto label_end;
  t->sph_z = (double*) mem_realloc((char* ) t->sph_z,
                                   sizeof(double) * (nold + npoint), 0);
  if (t->sph_z == nullptr) goto label_end;

  /* Check that random points are not too close from hard nodes */

  ecr = nold;
  for (int lec = 0; lec < npoint; lec++)
  {
    c1[0] = COORD(0, lec);
    c1[1] = COORD(1, lec);
    c1[2] = COORD(2, lec);

    /* Loop on the old points */

    found_close = -1;
    for (int i = 0; i < t->n_nodes && found_close < 0; i++)
    {
      c2[0] = t->sph_x[i];
      c2[1] = t->sph_y[i];
      c2[2] = t->sph_z[i];

      /* We calculate the distance in the plane (we cannot use ut_distance */
      /* which requires long-lat information */

      dist = 0.;
      for (int j = 0; j < 3; j++)
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
  t->sph_x = (double*) mem_realloc((char* ) t->sph_x, sizeof(double) * nech, 0);
  if (t->sph_x == nullptr) goto label_end;
  t->sph_y = (double*) mem_realloc((char* ) t->sph_y, sizeof(double) * nech, 0);
  if (t->sph_y == nullptr) goto label_end;
  t->sph_z = (double*) mem_realloc((char* ) t->sph_z, sizeof(double) * nech, 0);
  if (t->sph_z == nullptr) goto label_end;
  t->n_nodes = nech;

  /* Set the error return code */

  error = 0;

  label_end: coord = (double*) mem_free((char* ) coord);
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
void meshes_2D_sph_print(SphTriangle *t, int brief)
{
  double rlong, rlat;

  message("- Number of nodes   = %d\n", t->n_nodes);

  if (!brief && t->n_nodes > 0 && t->sph_x != nullptr && t->sph_y != nullptr
      && t->sph_z != nullptr)
  {
    message("\nCoordinates in Cartesian (R=1); then in Longitude - Latitude\n");
    for (int i = 0; i < t->n_nodes; i++)
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
int meshes_2D_sph_create(int verbose, SphTriangle *t)
{
  int *loc_near, *loc_next, *loc_lnew, error, skip_rnd, seed_memo;
  double *loc_dist, memo[3][3], ampli, value, cste;

  /* Initializations */

  error = 1;
  skip_rnd = (int) get_keypone("Skip_Random", 0);
  loc_near = loc_next = loc_lnew = nullptr;
  loc_dist = nullptr;
  if (t == (SphTriangle*) NULL || t->n_nodes < 3) return (1);

  /* Re-allocate the arrays within the SphTriangle structure */

  meshes_2D_sph_free(t, 1);
  t->sph_size = 6 * t->n_nodes - 12;
  t->sph_list = (int*) mem_alloc(sizeof(int) * t->sph_size, 0);
  if (t->sph_list == nullptr) goto label_end;
  for (int i = 0; i < t->sph_size; i++)
    t->sph_list[i] = ITEST;
  t->sph_lptr = (int*) mem_alloc(sizeof(int) * t->sph_size, 0);
  if (t->sph_lptr == nullptr) goto label_end;
  for (int i = 0; i < t->sph_size; i++)
    t->sph_lptr[i] = ITEST;
  t->sph_lend = (int*) mem_alloc(sizeof(int) * t->n_nodes, 0);
  if (t->sph_lend == nullptr) goto label_end;
  for (int i = 0; i < t->n_nodes; i++)
    t->sph_lend[i] = ITEST;

  /* Allocate local arrays */

  loc_lnew = (int*) mem_alloc(sizeof(int) * t->sph_size, 0);
  if (loc_lnew == nullptr) goto label_end;
  loc_near = (int*) mem_alloc(sizeof(int) * t->n_nodes, 0);
  if (loc_near == nullptr) goto label_end;
  loc_next = (int*) mem_alloc(sizeof(int) * t->n_nodes, 0);
  if (loc_next == nullptr) goto label_end;
  loc_dist = (double*) mem_alloc(sizeof(double) * t->n_nodes, 0);
  if (loc_dist == nullptr) goto label_end;

  /* Avoid having three first points colinear */

  if (!skip_rnd)
  {
    seed_memo = law_get_random_seed();
    law_set_random_seed(132231);
    ampli = 1.e-2;
    for (int i = 0; i < 3; i++)
    {
      memo[i][0] = t->sph_x[i];
      memo[i][1] = t->sph_y[i];
      memo[i][2] = t->sph_z[i];
      cste = (t->sph_x[i] * t->sph_x[i] + t->sph_y[i] * t->sph_y[i]
              + t->sph_z[i] * t->sph_z[i]);
      t->sph_x[i] += law_uniform(-ampli, ampli);
      t->sph_y[i] += law_uniform(-ampli, ampli);
      value = sqrt(
          cste - t->sph_x[i] * t->sph_x[i] - t->sph_y[i] * t->sph_y[i]);
      t->sph_z[i] = (t->sph_z[i] > 0.) ? value :
                                         -value;
    }
    law_set_random_seed(seed_memo);
  }

  (void) trmesh_(&t->n_nodes, t->sph_x, t->sph_y, t->sph_z, t->sph_list,
                 t->sph_lptr, t->sph_lend, loc_lnew, loc_near, loc_next,
                 loc_dist, &error);

  /* Restore the initial coordinates */

  for (int i = 0; i < 3; i++)
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
        meshes_2D_sph_print(t, ! DEBUG);
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

  label_end: loc_near = (int*) mem_free((char* ) loc_near);
  loc_next = (int*) mem_free((char* ) loc_next);
  loc_lnew = (int*) mem_free((char* ) loc_lnew);
  loc_dist = (double*) mem_free((char* ) loc_dist);
  return (error);
}

/*****************************************************************************/
/*!
 **  Load the contents of the SphTriangle structure into arrays
 **
 ** \return A Pointer to a MeshEStandard newly created structure
 **
 ** \param[in]  t      Pointer to the SphTriangle structure to be downloaded
 **
 *****************************************************************************/
MeshEStandard* meshes_2D_sph_load_vertices(SphTriangle *t)
{
  int ecr, lec, nt, error;
  double rlong, rlat;

  int natt = 2;
  int ntab = t->n_nodes;

  ecr = 0;
  VectorDouble rtab(ntab * natt, 0);
  for (int i = 0; i < ntab; i++)
  {
    GH::convertCart2Sph(t->sph_x[i], t->sph_y[i], t->sph_z[i], &rlong, &rlat);
    rtab[ecr++] = rlong;
    rtab[ecr++] = rlat;
  }

  ecr = 0;
  lec = 0;
  int nrow = 6;
  VectorInt ltri(2 * nrow * t->n_nodes, 0);
  (void) trlist_(&t->n_nodes, t->sph_list, t->sph_lptr, t->sph_lend, &nrow,
                 &nt, ltri.data(), &error);
  if (error) return nullptr;

  natt = 3;
  ntab = nt;
  VectorInt itab(ntab * natt, 0);
  for (int i = 0; i < ntab; i++)
  {
    for (int j = 0; j < natt; j++)
      itab[ecr + j] = ltri[lec + j];
    ecr += natt;
    lec += nrow;
  }

  MeshEStandard* amesh = new MeshEStandard();
  amesh->reset(2, 3, t->n_nodes, nt, rtab.data(), itab.data());
  return amesh;
}

/*****************************************************************************/
/*!
 **  Delete one facet
 **
 ** \param[in]  f       facet structure to be deleted
 **
 *****************************************************************************/
static void st_del_facet(tetgenio::facet *f)

{
  tetgenio::polygon *p;

  for (int ip = 0; ip < f->numberofpolygons; ip++)
  {
    p = &f->polygonlist[ip];
    p->vertexlist = (int*) mem_free((char* ) p->vertexlist);
  }
  f->polygonlist = (tetgenio::polygon*) mem_free((char* ) f->polygonlist);
}

/*****************************************************************************/
/*!
 **  Add one facet
 **
 ** \param[in]  f       facet structure
 ** \param[in]  npref   Number of facets
 ** \param[in]  ir1     Rank of the first sample
 ** \param[in]  ir2     Rank of the second sample
 ** \param[in]  ir3     Rank of the third sample
 ** \param[in]  ir4     Rank of the fourth sample
 **
 *****************************************************************************/
static void st_add_facet(tetgenio::facet *f,
                         int npref,
                         int ir1,
                         int ir2,
                         int ir3,
                         int ir4)
{
  tetgenio::polygon *p;

  npref--;
  f->numberofpolygons = 1;
  f->polygonlist = (tetgenio::polygon*) mem_alloc(sizeof(tetgenio::polygon), 1);
  f->numberofholes = 0;
  f->holelist = NULL;
  p = &f->polygonlist[0];
  p->numberofvertices = 4;
  p->vertexlist = (int*) mem_alloc(sizeof(int) * 4, 1);
  p->vertexlist[0] = npref + ir1;
  p->vertexlist[1] = npref + ir2;
  p->vertexlist[2] = npref + ir3;
  p->vertexlist[3] = npref + ir4;
}

/*****************************************************************************/
/*!
 **  Free the tetgenio structure
 **
 ** \param[in]  t           Pointer to the tetgenio structure to be freed
 **
 *****************************************************************************/
void meshes_3D_free(tetgenio *t)

{
  if (t == nullptr) return;
  t->pointlist = (double*) mem_free((char* ) t->pointlist);
  t->pointattributelist = (double*) mem_free((char* ) t->pointattributelist);
  t->pointmtrlist = (double*) mem_free((char* ) t->pointmtrlist);
  t->numberofpoints = 0;
  t->numberofpointattributes = 0;
  t->numberofpointmtrs = 0;

  t->tetrahedronlist = (int*) mem_free((char* ) t->tetrahedronlist);
  t->tetrahedronattributelist = (double*) mem_free(
      (char* ) t->tetrahedronattributelist);
  t->tetrahedronvolumelist = (double*) mem_free(
      (char* ) t->tetrahedronvolumelist);
  t->neighborlist = (int*) mem_free((char* ) t->neighborlist);
  t->numberoftetrahedra = 0;
  t->numberofcorners = 0;
  t->numberoftetrahedronattributes = 0;

  if (t->numberoffacets > 0) for (int i = 0; i < t->numberoffacets; i++)
    st_del_facet(&t->facetlist[i]);
  t->facetlist = (tetgenio::facet*) mem_free((char* ) t->facetlist);
  t->facetmarkerlist = (int*) mem_free((char* ) t->facetmarkerlist);
  t->numberoffacets = 0;

  t->holelist = (double*) mem_free((char* ) t->holelist);
  t->numberofholes = 0;

  t->regionlist = (double*) mem_free((char* ) t->regionlist);
  t->numberofregions = 0;

  t->facetconstraintlist = (double*) mem_free((char* ) t->facetconstraintlist);
  t->numberoffacetconstraints = 0;

  t->segmentconstraintlist = (double*) mem_free(
      (char* ) t->segmentconstraintlist);
  t->numberofsegmentconstraints = 0;

  t->trifacelist = (int*) mem_free((char* ) t->trifacelist);
  t->trifacemarkerlist = (int*) mem_free((char* ) t->trifacemarkerlist);
  t->numberoftrifaces = 0;

  t->edgelist = (int*) mem_free((char* ) t->edgelist);
  t->edgemarkerlist = (int*) mem_free((char* ) t->edgemarkerlist);
  t->numberofedges = 0;
}

/*****************************************************************************/
/*!
 **  Print the tetrahedrization
 **
 ** \param[in]  t         tetgenio structure
 ** \param[in]  brief     1 for a brief output; 0 otherwise
 **
 *****************************************************************************/
void meshes_3D_print(tetgenio *t, int brief)
{
  int ndim, i, j, lecp, leca, lech, lect, lecta, lecn;

  /* Initializations */

  ndim = 3;

  message("- Number of vertices   = %d\n", t->numberofpoints);
  message("- Number of attributes = %d\n", t->numberofpointattributes);

  /* List of vertices */

  if (!brief && t->numberofpoints > 0 && t->pointlist != nullptr)
  {
    lecp = leca = 0;
    message("List of Vertices\n");
    for (i = 0; i < t->numberofpoints; i++)
    {
      message("%3d", i + 1);
      for (j = 0; j < ndim; j++)
        message(" %lf", t->pointlist[lecp++]);
      if (t->pointmarkerlist != nullptr) message(" %d", t->pointmarkerlist[i]);
      for (j = 0; j < t->numberofpointattributes; j++)
        message(" %lf", t->pointattributelist[leca++]);
      message("\n");
    }
    message("\n");
  }

  // List of holes

  message("- Number of holes      = %d\n", t->numberofholes);

  if (!brief && t->numberofholes > 0 && t->holelist != nullptr)
  {
    lech = 0;
    message("List of Holes\n");
    for (i = 0; i < t->numberofholes; i++)
    {
      message("%3d", i + 1);
      for (j = 0; j < ndim; j++)
        message(" %lf", t->holelist[lech++]);
      message("\n");
    }
    message("\n");
  }

  // List of tetrahedra

  message("- Number of tetrahedra = %d\n", t->numberoftetrahedra);
  message("- Number of corners    = %d\n", t->numberofcorners);

  if (!brief && t->numberoftetrahedra > 0 && t->tetrahedronlist != nullptr)
  {
    lect = lecta = lecn = 0;
    message("List of Tetrahedra\n");
    for (i = 0; i < t->numberoftetrahedra; i++)
    {
      message("%3d", i + 1);
      for (j = 0; j < t->numberofcorners; j++)
        message(" %d", t->tetrahedronlist[lect++]);
      if (t->tetrahedronattributelist != nullptr)
        for (j = 0; j < t->numberoftetrahedronattributes; j++)
          message("%lf", t->tetrahedronattributelist[lecta++]);
      if (t->tetrahedronvolumelist != nullptr)
        message(" %lf", t->tetrahedronvolumelist[i]);
      message("\n");
    }
    message("\n");
  }
}

/*****************************************************************************/
/*!
 **  Add the bounding facets
 **
 ** \param[in]  t          tetgenio structure
 **
 *****************************************************************************/
static void tetgen_bounding(tetgenio *t)

{
  double *points, mini[3], maxi[3], value, delta;
  int ndim, npref;

  /* Initializations */

  ndim = 3;
  npref = t->numberofpoints;
  points = t->pointlist;

  /* Delete any previously existing facet */
  if (t->numberoffacets > 0)
  {
    t->facetlist = (tetgenio::facet*) mem_free((char* ) t->facetlist);
    t->facetmarkerlist = (int*) mem_free((char* ) t->facetmarkerlist);
    t->numberoffacets = 0;
  }

  /* Calculate the extrema */
  for (int i = 0; i < ndim; i++)
  {
    mini[i] = +1.e30;
    maxi[i] = -1.e30;
  }

  /* Loop on the meshes */

  for (int ip = 0; ip < npref; ip++)
  {
    for (int idim = 0; idim < ndim; idim++)
    {
      value = POINTS(ip, idim);
      if (value < mini[idim]) mini[idim] = value;
      if (value > maxi[idim]) maxi[idim] = value;
    }
  }

  /* Shift the extrema */

  for (int idim = 0; idim < ndim; idim++)
  {
    delta = maxi[idim] - mini[idim];
    mini[idim] -= delta / 100.;
    maxi[idim] += delta / 100.;
  }

  /* Create the new vertices */

  meshes_3D_from_points(1, &mini[0], &mini[1], &mini[2], t);
  meshes_3D_from_points(1, &maxi[0], &mini[1], &mini[2], t);
  meshes_3D_from_points(1, &maxi[0], &maxi[1], &mini[2], t);
  meshes_3D_from_points(1, &mini[0], &maxi[1], &mini[2], t);
  meshes_3D_from_points(1, &mini[0], &mini[1], &maxi[2], t);
  meshes_3D_from_points(1, &maxi[0], &mini[1], &maxi[2], t);
  meshes_3D_from_points(1, &maxi[0], &maxi[1], &maxi[2], t);
  meshes_3D_from_points(1, &mini[0], &maxi[1], &maxi[2], t);

  /* Add the facets */

  t->numberoffacets = 6;
  t->facetlist = (tetgenio::facet*) mem_alloc(sizeof(tetgenio::facet) * 6, 1);
  t->facetmarkerlist = (int*) mem_alloc(sizeof(int) * 6, 1);

  st_add_facet(&t->facetlist[0], npref, 1, 2, 3, 4);
  st_add_facet(&t->facetlist[1], npref, 5, 6, 7, 8);
  st_add_facet(&t->facetlist[2], npref, 1, 5, 6, 2);
  st_add_facet(&t->facetlist[3], npref, 2, 6, 7, 3);
  st_add_facet(&t->facetlist[4], npref, 3, 7, 8, 4);
  st_add_facet(&t->facetlist[5], npref, 4, 8, 5, 1);

  // Set 'in.facetmarkerlist'

  t->facetmarkerlist[0] = -1;
  t->facetmarkerlist[1] = -2;
  t->facetmarkerlist[2] = 0;
  t->facetmarkerlist[3] = 0;
  t->facetmarkerlist[4] = 0;
  t->facetmarkerlist[5] = 0;
}

/*****************************************************************************/
/*!
 **  Perform the tetrahedrization
 **
 ** \param[in]  verbose    Verbose option
 ** \param[in]  triswitch  Tetrahedralization option
 ** \param[in]  in         tetgenio structure for input
 ** \param[in]  out        tetgenio structure for output
 **
 *****************************************************************************/
void meshes_3D_create(int verbose,
                                      const String &triswitch,
                                      tetgenio *in,
                                      tetgenio *out)
{
  tetgenio *addin, *bgmin;
  int ndim, ncorner, ecr;

  /* Initializations */

  ndim = 3;
  ncorner = 4;
  addin = NULL;
  bgmin = NULL;

  // Create the bounding facets
  tetgen_bounding(in);

  // Perform the tetrahedrization
  tetrahedralize(triswitch.c_str(), in, out, addin, bgmin);

  // Shift the mesh list by 1 (to be compatible with Triangles) 
  if (out->numberoftetrahedra > 0)
  {
    ecr = 0;
    for (int i = 0; i < out->numberoftetrahedra; i++)
      for (int j = 0; j < out->numberofcorners; j++, ecr++)
        out->tetrahedronlist[ecr]++;
  }

  if (verbose)
  {
    message("Tetrahedralization performed (%s):\n", triswitch.c_str());
    meshes_3D_print(out, ! DEBUG);
    mesh_stats(ndim, ncorner, out->numberoftetrahedra, out->tetrahedronlist,
               out->pointlist);
  }
}

/*****************************************************************************/
/*!
 **  Define the memory areas by reading information from Db
 **
 ** \param[in]  db        Db structure where data are located
 ** \param[in]  t         Pointer to the tetgenio structure to be loaded
 **
 ** \remarks This function adds the vertices to an existing tetgenio structure
 ** \remarks Only the first three coordinates are considered
 **
 *****************************************************************************/
int meshes_3D_from_db(Db *db, tetgenio *t)
{
  int iech, nech, error, ecr, neff, ndim, nold;

  /* Initializations */

  if (db == nullptr) return (0);
  error = 1;
  nech = db->getSampleNumber();
  ndim = db->getNDim();
  if (ndim > 3) ndim = 3;

  /* Count the number of active samples */

  neff = db->getActiveSampleNumber();

  /* List of points */

  nold = t->numberofpoints;
  ecr = nold * ndim;
  t->pointlist = (double*) mem_realloc((char* ) t->pointlist,
                                       sizeof(double) * (nold + neff) * ndim, 0);
  if (t->pointlist == nullptr) goto label_end;

  for (iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    for (int idim = 0; idim < ndim; idim++)
      t->pointlist[ecr++] = db->getCoordinate(iech, idim);
  }
  t->numberofpoints = nold + neff;

  /* Set the error return code */

  error = 0;

  label_end: if (error)
  {
    t->pointlist = (double*) mem_free((char* ) t->pointlist);
    t->numberofpoints = 0;
  }
  return (error);
}

/*****************************************************************************/
/*!
 **  Define the tetrahedrization from the frame around dbin and/or dbout
 **
 ** \param[in]  dbin      Db structure for input file (optional)
 ** \param[in]  dbout     Db structure for output file (optional)
 ** \param[in]  t         Pointer to the tetgenio structure to be loaded
 **
 *****************************************************************************/
void meshes_3D_default(Db *dbin, Db *dbout, tetgenio *t)
{
  double *ext;
  int number;

  /* Initializations */

  number = 0;

  /* Get the extension of the Dbs */

  ext = st_get_db_extension(dbin, dbout, &number);

  /* Extend the triangulation */

  (void) meshes_3D_from_points(number, &EXT(0, 0), &EXT(1, 0), &EXT(2, 0), t);

  /* Core deallocation */

  ext = (double*) mem_free((char* ) ext);

  return;
}

/*****************************************************************************/
/*!
 **  Add fixed points to modify the tetrahedrization basis
 **
 ** \param[in]  nech      Number of added samples
 ** \param[in]  x,y,z     Array containing the coordinates of added points
 ** \param[in]  t         Pointer to the tetgenio structure to be loaded
 **
 ** \remarks This function adds the vertices to an existing
 ** \remarks tetgenio structure
 **
 *****************************************************************************/
int meshes_3D_from_points(int nech,
                                          double *x,
                                          double *y,
                                          double *z,
                                          tetgenio *t)
{
  int iech, error, ecr, ndim, nold;

  /* Initializations */

  error = 1;
  ndim = 3;

  /* List of points */

  nold = t->numberofpoints;
  ecr = nold * ndim;
  t->pointlist = (double*) mem_realloc((char* ) t->pointlist,
                                       sizeof(double) * (nold + nech) * ndim,
                                       0);
  if (t->pointlist == nullptr) goto label_end;
  for (iech = 0; iech < nech; iech++)
  {
    t->pointlist[ecr++] = x[iech];
    t->pointlist[ecr++] = y[iech];
    t->pointlist[ecr++] = z[iech];
  }
  t->numberofpoints = nold + nech;

  /* Set the error return code */

  error = 0;

  label_end: if (error)
  {
    t->pointlist = (double*) mem_free((char* ) t->pointlist);
    t->numberofpoints = 0;
  }
  return (error);
}

/****************************************************************************/
/*!
 **  Create the extended domain
 **
 ** \param[in]  dbout      Output Db grid structure
 ** \param[in]  gext       Array of domain dilation
 ** \param[in]  t          Tetrahedrization environment
 **
 *****************************************************************************/
void meshes_3D_extended_domain(Db *dbout, const double *gext, tetgenio *t)
{
  int number, flag_extend;
  double *ext;

  if (dbout == nullptr) return;
  if (gext == nullptr) return;
  ext = nullptr;

  /* Check that the extension parameters are correctly defined */

  flag_extend = 0;
  for (int idim = 0; idim < dbout->getNDim(); idim++)
    if (gext[idim] > 0) flag_extend = 1;
  if (!flag_extend) return;

  /* Dilate the grid */

  if (dbout->isGrid())
  {
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(dbout);
    ext = st_extend_grid(dbgrid, gext, &number);
  }
  else
    ext = st_extend_point(dbout, gext, &number);
  if (meshes_3D_from_points(number, &EXT(0, 0), &EXT(1, 0), &EXT(2, 0), t))
    goto label_end;

  /* Core deallocation */

  label_end: ext = (double*) mem_free((char* ) ext);
  return;
}

/*****************************************************************************/
/*!
 **  Load the contents of the tetgenio structure into arrays
 **
 ** \return A pointer to the newly created AMesh
 **
 ** \param[in]  t      Pointer to the tetgenio structure to be downloaded
 **
 *****************************************************************************/
MeshEStandard* meshes_3D_load_vertices(tetgenio *t)
{
  MeshEStandard* amesh = new MeshEStandard();
  amesh->reset(3, 4, t->numberofpoints, t->numberoftetrahedra,
                       t->pointlist, t->tetrahedronlist);
  return amesh;
}

/*****************************************************************************/
/*!
 **  Build the regular 3-D grid meshing
 **
 ** \return Pointer to the newly created AMesh structure
 **
 ** \param[in]  verbose   Verbose flag
 ** \param[in]  dbgrid    Db structure
 **
 *****************************************************************************/
AMesh* meshes_turbo_3D_grid_build(int verbose, DbGrid *dbgrid)
{
  int ndim = 3;
  int ncorner = 4;
  int nx = dbgrid->getNX(0);
  int ny = dbgrid->getNX(1);
  int nz = dbgrid->getNX(2);
  int number = nx * ny * nz;

  /* Core allocation */

  VectorInt meshes(number * ncorner * 6, 0);
  VectorInt order(number, ITEST);

  /* Store the indices per mesh */

  int nmesh = 0;
  for (int ix = 0; ix < nx - 1; ix++)
    for (int iy = 0; iy < ny - 1; iy++)
      for (int iz = 0; iz < nz - 1; iz++)
        for (int i = 0; i < 6; i++)
        {
          if (st_load_tetra(dbgrid, meshes, order, nmesh*ncorner,
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
 **  Free the segmentio structure
 **
 ** \param[in]  t           Pointer to the segmentio structure to be freed
 **
 ** \remark  When mode==1, the pointers which correspond to copies are
 ** \remark  simply set to NULL (not actually freed)
 **
 *****************************************************************************/
void meshes_1D_free(segmentio *t, int /*mode*/)
{
  if (t == (segmentio*) NULL) return;
  t->pointlist = (double*) mem_free((char* ) t->pointlist);
  t->pointattributelist = (double*) mem_free((char* ) t->pointattributelist);
  t->numberofpoints = 0;
  t->numberofpointattributes = 0;

  t->segmentlist = (int*) mem_free((char* ) t->segmentlist);
  t->numberofsegments = 0;
  t->numberofcorners = 0;
}

/*****************************************************************************/
/*!
 **  Initialize the segmentio structure
 **
 ** \param[in]  t      Pointer to the segmentio structure to be initialized
 **
 *****************************************************************************/
void meshes_1D_init(segmentio *t)
{
  t->pointlist = nullptr;
  t->pointattributelist = nullptr;
  t->numberofpoints = 0;
  t->numberofpointattributes = 0;

  t->segmentlist = nullptr;
  t->numberofsegments = 0;
  t->numberofcorners = 0;
}

/*****************************************************************************/
/*!
 **  Define the memory areas by reading information from Db
 **
 ** \return Error return code
 **
 ** \param[in]  db        Db structure where data are located
 ** \param[in]  t         Pointer to the triangulateio structure to be loaded
 **
 ** \remarks This function adds the vertices to an existing
 ** \remarks segmentio structure
 ** \remarks Only the first coordinate is considered
 **
 *****************************************************************************/
int meshes_1D_from_db(Db *db, segmentio *t)
{
  int iech, nech, error, ecr, neff, ndim, nold;

  /* Initializations */

  if (db == nullptr) return (0);
  error = 1;
  nech = db->getSampleNumber();
  ndim = db->getNDim();
  if (ndim > 1) ndim = 1;

  /* Count the number of active samples */

  neff = db->getActiveSampleNumber();

  /* List of points */

  nold = t->numberofpoints;
  ecr = nold * ndim;
  t->pointlist = (double*) mem_realloc((char* ) t->pointlist,
                                       sizeof(double) * (nold + neff) * ndim,
                                       0);
  if (t->pointlist == nullptr) goto label_end;

  for (iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    for (int idim = 0; idim < ndim; idim++)
      t->pointlist[ecr++] = db->getCoordinate(iech, idim);
  }
  t->numberofpoints = nold + neff;

  /* Set the error return code */

  error = 0;

  label_end: if (error)
  {
    t->pointlist = (double*) mem_free((char* ) t->pointlist);
    t->numberofpoints = 0;
  }
  return (error);
}

/*****************************************************************************/
/*!
 **  Define the triangulation from the frame around dbin and/or dbout
 **
 ** \param[in]  dbin      Db structure for input file (optional)
 ** \param[in]  dbout     Db structure for output file (optional)
 ** \param[in]  t         Pointer to the segmentio structure to be loaded
 **
 *****************************************************************************/
void meshes_1D_default(Db *dbin, Db *dbout, segmentio *t)
{
  double *ext;
  int number;

  /* Initializations */

  number = 0;

  /* Get the extension of the Dbs */

  ext = st_get_db_extension(dbin, dbout, &number);

  /* Extend the triangulation */

  (void) meshes_1D_from_points(number, ext, t);

  /* Core deallocation */

  ext = (double*) mem_free((char* ) ext);

  return;
}

/*****************************************************************************/
/*!
 **  Add fixed points to modify the segmentio basis
 **
 ** \param[in]  nech      Number of added samples
 ** \param[in]  x         Array containing the coordinates of added points
 ** \param[in]  t         Pointer to the segmentio structure to be loaded
 **
 ** \remarks This function adds the vertices (and holes) to an existing
 ** \remarks segmentio structure
 **
 *****************************************************************************/
int meshes_1D_from_points(int nech, double *x, segmentio *t)
{
  int iech, error, ecr, ndim, nold;

  /* Initializations */

  error = 1;
  ndim = 1;

  /* List of points */

  nold = t->numberofpoints;
  ecr = nold * ndim;
  t->pointlist = (double*) mem_realloc((char* ) t->pointlist,
                                       sizeof(double) * (nold + nech) * ndim,
                                       0);
  if (t->pointlist == nullptr) goto label_end;
  for (iech = 0; iech < nech; iech++)
  {
    t->pointlist[ecr++] = x[iech];
  }
  t->numberofpoints = nold + nech;

  /* Set the error return code */

  error = 0;

  label_end: if (error)
  {
    t->pointlist = (double*) mem_free((char* ) t->pointlist);
    t->numberofpoints = 0;
  }
  return (error);
}

/*****************************************************************************/
/*!
 **  Print the contents of the segmentio structure
 **
 ** \param[in]  t         Pointer to the segmentio structure to be printed
 ** \param[in]  brief     1 for a brief output; 0 otherwise
 **
 *****************************************************************************/
void meshes_1D_print(segmentio *t, int brief)
{
  int ndim, i, j, lecp, leca, lect, lecta, lecn;

  ndim = 2;

  // List of vertices 

  message("- Number of vertices   = %d\n", t->numberofpoints);
  message("- Number of attributes = %d\n", t->numberofpointattributes);

  if (!brief && t->numberofpoints > 0 && t->pointlist != nullptr)
  {
    lecp = leca = 0;
    for (i = 0; i < t->numberofpoints; i++)
    {
      message("%3d", i + 1);
      for (j = 0; j < ndim; j++)
        message(" %lf", t->pointlist[lecp++]);
      for (j = 0; j < t->numberofpointattributes; j++)
        message(" %lf", t->pointattributelist[leca++]);
      message("\n");
    }
    message("\n");
  }

  // List of segments

  message("- Number of segments   = %d\n", t->numberofsegments);
  message("- Number of corners    = %d\n", t->numberofcorners);

  if (!brief && t->numberofsegments > 0 && t->segmentlist != nullptr)
  {
    lect = lecta = lecn = 0;
    for (i = 0; i < t->numberofsegments; i++)
    {
      message("%3d", i + 1);
      for (j = 0; j < t->numberofcorners; j++)
        message(" %d", t->segmentlist[lect++]);
      message("\n");
    }
    message("\n");
  }
}

/*****************************************************************************/
/*!
 **  Load the contents of the segmentio structure into arrays
 **
 ** \param[in]  t      Pointer to the segmentio structure to be downloaded
 **
 *****************************************************************************/
MeshEStandard* meshes_1D_load_vertices(segmentio *t)
{
  MeshEStandard* amesh = new MeshEStandard();
  amesh->reset(1, 2, t->numberofpoints, t->numberofsegments,
                       t->pointlist, t->segmentlist);
  return amesh;
}

/*****************************************************************************/
/*!
 **  Build the regular meshing from a 1-D grid
 **
 ** \return The newly created AMesh structure
 **
 ** \param[in]  verbose   Verbose flag
 ** \param[in]  dbgrid    Db structure
 **
 *****************************************************************************/
AMesh* meshes_turbo_1D_grid_build(int verbose, DbGrid *dbgrid)
{
  int ndim = 1;
  int ncorner = 2;
  int nx = dbgrid->getNX(0);
  int number = nx - 1;

  /* Core allocation */

  VectorInt meshes(number * ncorner, 0);
  VectorInt order(number, ITEST);

  /* Store the indices per mesh */

  int nmesh = 0;
  for (int ix = 0; ix < nx - 1; ix++)
  {
    if (st_load_segment(dbgrid, meshes, order, nmesh*ncorner,
                        ix + MSS(1, 0, 1, 0, 0), ix + MSS(1, 0, 1, 1, 0)))
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
 **  Perform the meshing in 1-D
 **
 ** \param[in]  verbose    Verbose option
 ** \param[in]  in         segmentio structure for input
 ** \param[in]  out        segmentio structure for output
 **
 *****************************************************************************/
void meshes_1D_create(int verbose, segmentio *in, segmentio *out)
{
  int *rank, ndim, ncorner, number, ecr;

  ndim = 1;
  ncorner = 2;

  /* Processing */

  meshes_1D_init(out);

  /* Create pointlist */

  number = in->numberofpoints * in->numberofpointattributes;
  out->numberofpoints = in->numberofpoints;
  out->pointlist = (double*) mem_alloc(sizeof(double) * number, 1);
  for (int i = 0; i < number; i++)
    out->pointlist[i] = in->pointlist[i];

  /* Create pointattributelist */

  number = in->numberofpointattributes;
  out->numberofpointattributes = number;
  out->pointattributelist = (double*) mem_alloc(sizeof(double) * number, 1);
  for (int i = 0; i < number; i++)
    out->pointattributelist[i] = in->pointattributelist[i];

  /* Sorting the vertices */

  rank = (int*) mem_alloc(sizeof(int) * number, 1);
  for (int i = 0; i < number; i++)
    rank[i] = i;
  ut_sort_double(1, number, rank, out->pointlist);

  /* Create the list of segments */

  out->numberofsegments = number - 1;
  out->numberofcorners = 2;
  number = out->numberofsegments * out->numberofcorners;
  out->segmentlist = (int*) mem_alloc(sizeof(int) * number, 1);
  ecr = 0;
  for (int i = 0; i < out->numberofsegments; i++)
  {
    out->segmentlist[ecr++] = rank[i];
    out->segmentlist[ecr++] = rank[i + 1];
  }

  /* Optional printout */

  if (verbose)
  {
    message("Meshing performed\n");
    meshes_1D_print(out, ! DEBUG);
    mesh_stats(ndim, ncorner, out->numberofsegments, out->segmentlist,
               out->pointlist);
  }
}

/****************************************************************************/
/*!
 **  Create the extended domain
 **
 ** \param[in]  dbout      Output Db grid structure
 ** \param[in]  gext       Array of domain dilation
 ** \param[in]  t          segmentio structure
 **
 *****************************************************************************/
void meshes_1D_extended_domain(Db *dbout, const double *gext, segmentio *t)
{
  int number, flag_extend;
  double *ext;

  /* Initializations */

  if (dbout == nullptr) return;
  if (gext == nullptr) return;
  ext = nullptr;

  /* Check that the extension parameters are correctly defined */

  flag_extend = 0;
  for (int idim = 0; idim < dbout->getNDim(); idim++)
    if (gext[idim] > 0) flag_extend = 1;
  if (!flag_extend) return;

  /* Dilate the grid */

  if (dbout->isGrid())
  {
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(dbout);
    ext = st_extend_grid(dbgrid, gext, &number);
  }
  else
    ext = st_extend_point(dbout, gext, &number);
  if (meshes_1D_from_points(number, ext, t)) goto label_end;

  /* Core deallocation */

  label_end: ext = (double*) mem_free((char* ) ext);
  return;
}

