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
#include "geoslib_old_f.h"

#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "ExternalTools/LinkTriangle.hpp"
#include "ExternalTools/LinkTetrahedron.hpp"
#include "ExternalTools/LinkSphTriangle.hpp"
#include "ExternalTools/MeshEStandardExt.hpp"
#include "Basic/AException.hpp"
#include "Db/Db.hpp"
#include "Basic/VectorHelper.hpp"

#include "csparse_f.h"

MeshEStandardExt::MeshEStandardExt()
  : MeshEStandard()
{
}

MeshEStandardExt::MeshEStandardExt(const MeshEStandardExt &m)
  : MeshEStandard(m)
{
}

MeshEStandardExt& MeshEStandardExt::operator= (const MeshEStandardExt &m)
{
  if (this != &m)
  {
    MeshEStandard::operator=(m);
  }
  return *this;
}

MeshEStandardExt::~MeshEStandardExt()
{
}

/****************************************************************************/
/*!
** Create the meshing
**
** \param[in]  dbin            Pointer to the input Db (optional)
** \param[in]  dbout           Pointer to the output Db (optional)
** \param[in]  dilate          Dilation of the Bounding box (optional)
** \param[in]  triswitch       Construction switch
** \param[in]  verbose         Verbose flag
**
*****************************************************************************/
int MeshEStandardExt::resetFromDb(Db *dbin,
                                  Db *dbout,
                                  const VectorDouble &dilate,
                                  const String &triswitch,
                                  bool verbose)
{
  int error = 1;

  // Discover the space dimension
  int ndim = 0;
  if (dbin != NULL)
    ndim = dbin->getNDim();
  if (dbout != NULL)
  {
    if (ndim > 0 && ndim != dbout->getNDim())
    {
      messerr("Space dimension should be the same in 'dbin'(%d) and 'dbout'(%d)",
               dbin->getNDim(), dbout->getNDim());
      return 1;
    }
    ndim = dbout->getNDim();
  }
  _setNDim(ndim);

  // Preliminary checks

  if (! dilate.empty())
  {
    for (int idim=0; idim<ndim; idim++)
      if (dilate[idim] < 0)
      {
        messerr("The dilation[%d] is negative (%d)",idim+1,dilate[idim]);
        return 1;
      }
  }

  // Dispatch according to Space Dimension

  if (ndim == 1)
    error = _create1D(verbose,dbin,dbout,dilate);
  else if (ndim == 2)
    error = _create2D(verbose,dbin,dbout,dilate,triswitch.c_str());
  else if (ndim == 3)
    error = _create3D(verbose,dbin,dbout,dilate,triswitch.c_str());
  else
  {
    messerr("Meshing is only provided for Space Dimension 1, 2 or 3");
  }
  if (error) return(error);

  // Define and store the Bounding Box extension

  _defineBoundingBox();

  // Optional printout

  if (verbose) messageFlush(toString());

  return(error);
}

/****************************************************************************/
/*!
** Create the meshing in 1D
**
** \param[in]  verbose         Verbose flag
** \param[in]  dbin            Pointer to the input Db (optional)
** \param[in]  dbout           Pointer to the output Db (optional)
** \param[in]  dilate          Dilation of the Bounding box (optional)
**
*****************************************************************************/
int MeshEStandardExt::_create1D(int verbose,
                                Db *dbin,
                                Db *dbout,
                                const VectorDouble &dilate)
{
  segmentio in, out;

  /* Initialize the Meshing structure */

  meshes_1D_init(&in);
  meshes_1D_init(&out);

  /* Set the control points for the triangulation */

  bool flag_defined = false;
  if (dbout != nullptr)
  {
    if (meshes_1D_from_db(dbout,&in)) return 1;
    flag_defined = true;
  }
  if (dbin != nullptr)
  {
    if (meshes_1D_from_db(dbin,&in)) return 1;
    flag_defined = true;
  }
  if (! flag_defined)
  {
    meshes_1D_default(dbin,dbout,&in);
    flag_defined = 1;
  }

  /* Extend the domain if 'dilate' is specified */

  if (dbout != nullptr && ! dilate.empty())
    meshes_1D_extended_domain(dbout,dilate.data(),&in);

  /* Perform the triangulation */

  meshes_1D_create(verbose,&in,&out);

  /* Create the final meshing */

  reset(1, 2, out.numberofpoints, out.numberofsegments, out.pointlist,
        out.segmentlist);

  meshes_1D_free(&in,1);
  meshes_1D_free(&out,0);
  return 0;
}

/****************************************************************************/
/*!
** Create the meshing in 2D
**
** \param[in]  verbose         Verbose flag
** \param[in]  dbin            Pointer to the input Db (optional)
** \param[in]  dbout           Pointer to the output Db (optional)
** \param[in]  dilate          Dilation of the Bounding box (optional)
** \param[in]  triswitch       Construction switch
**
*****************************************************************************/
int MeshEStandardExt::_create2D(bool verbose,
                                Db *dbin,
                                Db *dbout,
                                const VectorDouble &dilate,
                                const String& triswitch)
{
  triangulateio in, out, vorout;

  /* Initialize the Meshing structure */

  meshes_2D_init(&in);
  meshes_2D_init(&out);
  meshes_2D_init(&vorout);

  /* Set the control points for the triangulation */

  bool flag_defined = false;
  if (dbout != nullptr)
  {
    if (meshes_2D_from_db(dbout,1,&in)) return 1;
    flag_defined = true;
  }
  if (dbin != nullptr)
  {
    if (meshes_2D_from_db(dbin,1, &in)) return 1;
    flag_defined = true;
  }
  if (! flag_defined)
  {
    meshes_2D_default(dbin,dbout,&in);
    flag_defined = 1;
  }

  /* Extend the domain if 'dilate' is specified */

  if (dbout != nullptr && ! dilate.empty())
    meshes_2D_extended_domain(dbout,dilate.data(),&in);

  /* Perform the triangulation */

  meshes_2D_create(verbose,triswitch,&in,&out,&vorout);

  /* Coordinates of the triangle vertices */

  reset(2, 3, out.numberofpoints, out.numberoftriangles, out.pointlist,
        out.trianglelist, false);

  meshes_2D_free(&in,1);
  meshes_2D_free(&out,0);
  meshes_2D_free(&vorout,0);
  return 0;
}

/****************************************************************************/
/*!
** Create the meshing in 3D
**
** \param[in]  verbose         Verbose flag
** \param[in]  dbin            Pointer to the input Db (optional)
** \param[in]  dbout           Pointer to the output Db (optional)
** \param[in]  dilate          Dilation of the Bounding box (optional)
** \param[in]  triswitch       Construction switch
**
*****************************************************************************/
int MeshEStandardExt::_create3D(bool verbose,
                                Db *dbin,
                                Db *dbout,
                                const VectorDouble &dilate,
                                const String& triswitch)
{
  tetgenio in, out;

  /* Set the control points for the tetrahedralization */

  bool flag_defined = false;
  if (dbout != nullptr)
  {
    if (meshes_3D_from_db(dbout,&in)) return 1;
    flag_defined = true;
  }
  if (dbin != nullptr)
  {
    if (meshes_3D_from_db(dbin,&in)) return 1;
    flag_defined = true;
  }
  if (! flag_defined)
    meshes_3D_default(dbin,dbout,&in);

  /* Extend the domain if 'dilate' is specified */

  if (dbout != nullptr && ! dilate.empty())
    meshes_3D_extended_domain(dbout,dilate.data(),&in);

  /* Perform the tetrahedralization */

  meshes_3D_create(verbose,triswitch,&in,&out);

  /* Coordinates of the vertices */

  reset(3, 4, out.numberofpoints, out.numberoftetrahedra, out.pointlist,
        out.tetrahedronlist);

  meshes_3D_free(&in);
  meshes_3D_free(&out);
  return 0;
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
AMesh* MeshEStandardExt::spde_mesh_load(Db *dbin,
                                        Db *dbout,
                                        const VectorDouble &gext,
                                        const String &triswitch,
                                        bool verbose)
{
  int ndim_loc = 0;
  if (dbin != nullptr) ndim_loc = MAX(ndim_loc, dbin->getNDim());
  if (dbout != nullptr) ndim_loc = MAX(ndim_loc, dbout->getNDim());
  bool flag_sphere = isDefaultSpaceSphere();

  // Processing

  if (verbose)  message("Generating the meshes\n");

  if (flag_sphere)
  {
    messerr("This method cannot be run on a Spherical space");
    return nullptr;
  }
  else
  {
    if (verbose) message("Using Regular Meshing\n");

    if (ndim_loc == 1)
    {
      return _load1D(verbose, dbin, dbout, gext);
    }
    else if (ndim_loc == 2)
    {
      return _load2D(verbose, dbin, dbout, gext, triswitch);
    }
    else if (ndim_loc == 3)
    {
      return _load3D(verbose, dbin, dbout, gext, triswitch);
    }
  }
  return nullptr;
}

/****************************************************************************/
/*!
 **  Load the segments and vertex coordinates
 **
 ** \return  AMesh structure
 **
 ** \param[in]  verbose    Verbose option
 ** \param[in]  dbin       Db structure for the conditioning data
 ** \param[in]  dbout      Db structure of the grid
 ** \param[in]  gext       Array of domain dilation
 **
 *****************************************************************************/
AMesh* MeshEStandardExt::_load1D(bool verbose,
                                 Db *dbin,
                                 Db *dbout,
                                 const VectorDouble &gext)
{
  SYMBOL_UNUSED(verbose);
  segmentio in, out;

  /* Initialize the Meshing output structure */

  meshes_1D_init(&in);
  meshes_1D_init(&out);

  /* Set the control points for the triangulation */

  if (dbout != nullptr)
  {
    if (meshes_1D_from_db(dbout, &in)) return nullptr;
  }
  if (dbin != nullptr)
  {
    if (meshes_1D_from_db(dbin, &in)) return nullptr;
  }
  if (dbin == nullptr && dbout == nullptr)
  {
    meshes_1D_default(dbin, dbout, &in);
  }

  /* Extend the domain if gext is specified */

  if (!gext.empty())
  {
    meshes_1D_extended_domain(dbout, gext.data(), &in);
  }

  /* Perform the triangulation */

  meshes_1D_create(verbose, &in, &out);

  /* Create the final meshing */

  AMesh* amesh = meshes_1D_load_vertices(&out);

  /* Set the error return code */

  meshes_1D_free(&in, 1);
  meshes_1D_free(&out, 0);
  return amesh;
}

/****************************************************************************/
/*!
 **  Load the triangles and vertex coordinates
 **
 ** \return  AMesh structure
 **
 ** \param[in]  verbose    Verbose flag
 ** \param[in]  dbin       Db structure for the conditioning data
 ** \param[in]  dbout      Db structure of the grid
 ** \param[in]  gext       Array of domain dilation
 ** \param[in]  triswitch  triswitch option
 **
 *****************************************************************************/
AMesh* MeshEStandardExt::_load2D(bool verbose,
                                 Db *dbin,
                                 Db *dbout,
                                 const VectorDouble &gext,
                                 const String &triswitch)
{
  SYMBOL_UNUSED(verbose);
  triangulateio in, out, vorout;

  /* Initialize the Meshing output structure */

  meshes_2D_init(&in);
  meshes_2D_init(&out);
  meshes_2D_init(&vorout);

  /* Set the control points for the triangulation */

  if (dbout != nullptr)
  {
    if (meshes_2D_from_db(dbout, 1, &in)) return nullptr;
  }
  if (dbin != nullptr)
  {
    if (meshes_2D_from_db(dbin, 1, &in)) return nullptr;
  }
  if (dbin == nullptr && dbout == nullptr)
  {
    meshes_2D_default(dbin, dbout, &in);
  }

  /* Extend the domain if gext is specified */

  if (!gext.empty())
  {
    meshes_2D_extended_domain(dbout, gext.data(), &in);
  }

  /* Perform the triangulation */

  meshes_2D_create(verbose, triswitch, &in, &out, &vorout);

  /* Create the final meshing */

  AMesh* amesh = meshes_2D_load_vertices(&out);

  meshes_2D_free(&in, 1);
  meshes_2D_free(&out, 0);
  meshes_2D_free(&vorout, 0);
  return amesh;
}

/****************************************************************************/
/*!
 **  Load the tetrahedra and vertex coordinates
 **
 ** \return  Pointer to the newly created AMesh structure
 **
 ** \param[in]  verbose    Verbose flag
 ** \param[in]  dbin       Db structure for the conditioning data
 ** \param[in]  dbout      Db structure of the grid
 ** \param[in]  gext       Array of domain dilation
 ** \param[in]  triswitch  triswitch option
 **
 *****************************************************************************/
AMesh* MeshEStandardExt::_load3D(bool verbose,
                                 Db *dbin,
                                 Db *dbout,
                                 const VectorDouble &gext,
                                 const String &triswitch)
{
  tetgenio in, out;

  /* Set the control points for the tetrahedralization */

  if (dbout != nullptr)
  {
    if (meshes_3D_from_db(dbout, &in)) return nullptr;
  }
  if (dbin != nullptr)
  {
    if (meshes_3D_from_db(dbin, &in)) return nullptr;
  }
  if (dbin == nullptr && dbout == nullptr)
  {
    meshes_3D_default(dbin, dbout, &in);
  }

  /* Extend the domain if gext is specified */

  if (!gext.empty())
  {
    meshes_3D_extended_domain(dbout, gext.data(), &in);
  }

  /* Perform the tetrahedralization */

  meshes_3D_create(verbose, triswitch, &in, &out);

  /* Final meshing */

  AMesh* amesh = meshes_3D_load_vertices(&out);

  meshes_3D_free(&in);
  meshes_3D_free(&out);
  return amesh;
}

