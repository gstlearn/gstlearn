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
#include "geoslib_old_f.h"

#include "ExternalTools/LinkTriangle.hpp"

#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Mesh/MeshEStandard.hpp"

#include <math.h>

#define POINTS(ip,idim)   (points[(ip) * ndim + (idim)])
#define EXT(idim,ip)      (ext[(idim) * number + (ip)])
#define MESHES(imesh,j)   (meshes[(imesh) * ncorner + (j)] - 1)
#define FAULTS(ip,idim)   (faults[nfaults * (idim) + (ip)])

#define DEBUG 0

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
      jcorn = (icorn + 1 == ncorner) ? 0 : icorn + 1;
      i2 = MESHES(imesh, jcorn);

      /* Loop on the faults vertices */

      for (int j2 = 1; j2 < nfaults && !skip; j2++)
      {
        j1 = j2 - 1;
        skip += (GH::segmentIntersect(POINTS(i1, 0), POINTS(i1, 1),
                                      POINTS(i2, 0), POINTS(i2, 1),
                                      FAULTS(j1, 0), FAULTS(j1, 1),
                                      FAULTS(j2, 0), FAULTS(j2, 1), &xx, &yy) != 1);
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

  ext = get_db_extension(dbin, dbout, &number);

  /* Extend the triangulation */

  (void) meshes_1D_from_points(number, ext, t);

  /* Core deallocation */

  ext = (double*) mem_free((char* ) ext);

  return;
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
  amesh->reset(1, 2, t->numberofpoints, t->numberofsegments, t->pointlist,
               t->segmentlist, false);
  return amesh;
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
    ext = extend_grid(dbgrid, gext, &number);
  }
  else
    ext = extend_point(dbout, gext, &number);
  if (meshes_1D_from_points(number, ext, t)) goto label_end;

  /* Core deallocation */

  label_end: ext = (double*) mem_free((char* ) ext);
  return;
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
  ncode = db->getLocNumber(ELoc::C);
  if (ndim > 2) ndim = 2;

  /* Count the number of active samples */

  neff = nhole = 0;
  for (iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    if (use_code && ncode > 0 && db->getLocVariable(ELoc::C,iech,0) < 0)
      nhole++;
    else
      neff++;
  }

  /* List of points */

  nold = t->numberofpoints;
  ecr = nold * ndim;
  t->pointlist = (double*) mem_realloc((char* ) t->pointlist,
                                       sizeof(double) * (nold + neff) * ndim, 0);
  if (t->pointlist == nullptr) goto label_end;

  for (iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    if (use_code && ncode > 0 && db->getLocVariable(ELoc::C,iech,0) < 0) continue;
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
                                        sizeof(double) * (nold + nhole) * ndim, 0);
    if (t->holelist == nullptr) goto label_end;
    for (iech = 0; iech < nech; iech++)
    {
      if (!db->isActive(iech)) continue;
      if (!(use_code && ncode > 0 && db->getLocVariable(ELoc::C,iech,0) < 0)) continue;
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
  int iech, ecr, nold;

  /* Initializations */

  int error = 1;
  int ndim = 2;

  /* List of points */

  nold = t->numberofpoints;
  ecr = nold * ndim;
  t->pointlist = (double*) mem_realloc((char* ) t->pointlist,
                                       sizeof(double) * (nold + nech) * ndim, 0);
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
 **  Define the triangulation from the frame around dbin and/or dbout
 **
 ** \param[in]  dbin      Db structure for input file (optional)
 ** \param[in]  dbout     Db structure for output file (optional)
 ** \param[in]  t         Pointer to the triangulateio structure to be loaded
 **
 *****************************************************************************/
void meshes_2D_default(Db *dbin, Db *dbout, triangulateio *t)
{
  int number = 0;

  /* Get the extension of the Dbs */

  double* ext = get_db_extension(dbin, dbout, &number);

  /* Extend the triangulation */

  (void) meshes_2D_from_points(number, &EXT(0, 0), &EXT(1, 0), t);

  /* Core deallocation */

  ext = (double*) mem_free((char* ) ext);

  return;
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
  int error = 1;
  int ndim = 2;

  /* List of segments */

  t->segmentlist = (int*) mem_alloc(sizeof(int) * nseg * ndim, 0);
  if (t->segmentlist == nullptr) goto label_end;
  if (ncol > ndim)
  {
    t->segmentmarkerlist = (int*) mem_alloc(sizeof(int) * nseg, 0);
    if (t->segmentmarkerlist == nullptr) goto label_end;
  }

  for (int i = 0; i < nseg; i++)
  {
    for (int j = 0; j < ndim; j++)
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
  int i, j, lecp, leca, lecs, lech, lect, lecta, lecn;

  int ndim = 2;

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
  amesh->reset(2, 3, t->numberofpoints, t->numberoftriangles, t->pointlist,
               t->trianglelist, false);
  return amesh;
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
void meshes_2D_extended_domain(Db *dbout, const double *gext, triangulateio *t)
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
    ext = extend_grid(dbgrid, gext, &number);
  }
  else
    ext = extend_point(dbout, gext, &number);
  if (meshes_2D_from_points(number, &EXT(0, 0), &EXT(1, 0), t)) goto label_end;

  /* Core deallocation */

  label_end: ext = (double*) mem_free((char* ) ext);
  return;
}

