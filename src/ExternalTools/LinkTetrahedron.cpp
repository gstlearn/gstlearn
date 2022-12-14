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

#include "ExternalTools/LinkTetrahedron.hpp"

#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Mesh/MeshEStandard.hpp"

#include <math.h>

#define POINTS(ip,idim)   (points[(ip) * ndim + (idim)])
#define EXT(idim,ip)      (ext[(idim) * number + (ip)])

#define DEBUG 0

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
 **  Print the tetrahedrization
 **
 ** \param[in]  t         tetgenio structure
 ** \param[in]  brief     1 for a brief output; 0 otherwise
 **
 *****************************************************************************/
void meshes_3D_print(tetgenio *t, int brief)
{
  int i, j, lecp, leca, lech, lect, lecta, lecn;

  /* Initializations */

  int ndim = 3;

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

  ext = get_db_extension(dbin, dbout, &number);

  /* Extend the triangulation */

  (void) meshes_3D_from_points(number, &EXT(0, 0), &EXT(1, 0), &EXT(2, 0), t);

  /* Core deallocation */

  ext = (double*) mem_free((char* ) ext);

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
  amesh->reset(3, 4, t->numberofpoints, t->numberoftetrahedra, t->pointlist,
               t->tetrahedronlist, false);
  return amesh;
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
    ext = extend_grid(dbgrid, gext, &number);
  }
  else
    ext = extend_point(dbout, gext, &number);
  if (meshes_3D_from_points(number, &EXT(0, 0), &EXT(1, 0), &EXT(2, 0), t))
    goto label_end;

  /* Core deallocation */

  label_end: ext = (double*) mem_free((char* ) ext);
  return;
}
