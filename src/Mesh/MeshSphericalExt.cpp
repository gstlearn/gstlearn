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
#include "Mesh/MeshSphericalExt.hpp"
#include "Mesh/LinkSphTriangle.hpp"
#include "geoslib_old_f.h"

#include "Matrix/MatrixRectangular.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorHelper.hpp"
#include "Mesh/MeshEStandard.hpp"
#include "Db/Db.hpp"

// External library /// TODO : Dependency to gmtsph to be removed

MeshSphericalExt::MeshSphericalExt()
  : MeshSpherical()
{
}

MeshSphericalExt::MeshSphericalExt(const MeshSphericalExt &m)
  : MeshSpherical(m)
{
}

MeshSphericalExt& MeshSphericalExt::operator= (const MeshSphericalExt &m)
{
  if (this != &m)
  {
    MeshSpherical::operator=(m);
  }
  return *this;
}

MeshSphericalExt::~MeshSphericalExt()
{
}

/****************************************************************************/
/*!
** Create the meshing
**
** \param[in]  dbin            Pointer to the input Db (optional)
** \param[in]  dbout           Pointer to the output Db (optional)
** \param[in]  triswitch       Construction switch
** \param[in]  verbose         Verbose flag
**
*****************************************************************************/
int MeshSphericalExt::resetFromDb(Db *dbin,
                                  Db *dbout,
                                  const String &triswitch,
                                  bool verbose)
{
  SphTriangle in;

  /* Initializations */

  int ndim_ref = 2;

  /* Define the environment */

  _setNDim(ndim_ref);

  /* Initialize the Meshing output structure */

  meshes_2D_sph_init(&in);

  /* Set the control points for the triangulation */

  if (dbout != nullptr)
  {
    if (meshes_2D_sph_from_db(dbout,&in)) return 1;
  }
  if (dbin != nullptr)
  {
    if (meshes_2D_sph_from_db(dbin,&in)) return 1;
  }

  /* Add auxiliary random points */

  if (meshes_2D_sph_from_auxiliary(triswitch.c_str(),&in)) return 1;

  /* Perform the triangulation */

  if (meshes_2D_sph_create(verbose,&in)) return 1;

  /* Final meshing */

  _meshesSphLoadVertices(&in);

  meshes_2D_sph_free(&in,0);
  return 0;
}

/*****************************************************************************/
/*!
 **  Load the contents of the SphTriangle structure into arrays
 **
 ** \param[in]  t      Pointer to the SphTriangle structure to be downloaded
 **
 *****************************************************************************/
void MeshSphericalExt::_meshesSphLoadVertices(SphTriangle *t)
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
  if (error) return;

  natt = 3;
  ntab = nt;
  VectorInt itab(ntab * natt, 0);
  for (int i = 0; i < ntab; i++)
  {
    for (int j = 0; j < natt; j++)
      itab[ecr + j] = ltri[lec + j] - 1;
    ecr += natt;
    lec += nrow;
  }

  reset(2, 3, rtab, itab, false);
}

/****************************************************************************/
/*!
 **  Load the AMesh structure
 **
 ** \return  Pointer on the newly allocated AMesh
 **
 ** \param[in]  dbin       Db structure for the conditioning data
 ** \param[in]  dbout      Db structure of the grid
 ** \param[in]  gext       Array of domain dilation
 ** \param[in]  triswitch  Triswitch option
 ** \param[in]  verbose    Verbose option
 **
 *****************************************************************************/
AMesh* MeshSphericalExt::spde_mesh_load(Db *dbin,
                                        Db *dbout,
                                        const VectorDouble &gext,
                                        const String &triswitch,
                                        bool verbose)
{
  DECLARE_UNUSED(gext);

  int ndim_loc = 0;
  if (dbin != nullptr) ndim_loc = MAX(ndim_loc, dbin->getNDim());
  if (dbout != nullptr) ndim_loc = MAX(ndim_loc, dbout->getNDim());
  bool flag_sphere = isDefaultSpaceSphere();

  // Processing

  if (verbose)  message("Generating the meshes\n");

  if (flag_sphere)
  {
    if (verbose) message("Using Regular Meshing on Sphere\n");
    return _load2DSph(verbose, dbin, dbout, triswitch);
  }
  else
  {
    messerr("This method cannot be used for non Spherical Meshing");
    return nullptr;
  }
}

/****************************************************************************/
/*!
 **  Load the spherical triangles and vertex coordinates
 **
 ** \return  The newly created AMesh structure
 **
 ** \param[in]  verbose    Verbose option
 ** \param[in]  dbin       Db structure for the conditioning data
 ** \param[in]  dbout      Db structure of the grid
 ** \param[in]  triswitch  Triangulation option
 **
 *****************************************************************************/
AMesh* MeshSphericalExt::_load2DSph(bool verbose,
                                    Db *dbin,
                                    Db *dbout,
                                    const String &triswitch)
{
  DECLARE_UNUSED(verbose);
  SphTriangle in;

  /* Initialize the Meshing output structure */

  meshes_2D_sph_init(&in);

  /* Set the control points for the triangulation */

  if (dbout != nullptr)
  {
    if (meshes_2D_sph_from_db(dbout, &in)) return nullptr;
  }
  if (dbin != nullptr)
  {
    if (meshes_2D_sph_from_db(dbin, &in)) return nullptr;
  }

  /* Add auxiliary random points */

  if (meshes_2D_sph_from_auxiliary(triswitch, &in)) return nullptr;

  /* Perform the triangulation */

  if (meshes_2D_sph_create(verbose, &in)) return nullptr;

  /* Final meshing */

  MeshSphericalExt* amesh = new MeshSphericalExt();
  amesh->_meshesSphLoadVertices(&in);

  meshes_2D_sph_free(&in, 0);

  return amesh;
}
