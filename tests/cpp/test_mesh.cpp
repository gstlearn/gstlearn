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
#include "geoslib_d.h"
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "csparse_f.h"

#include "Model/Model.hpp"
#include "Mesh/MeshEStandard.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Mesh/MeshSpherical.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Mesh/MeshFactory.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "Db/ELoadBy.hpp"
#include "Space/ASpaceObject.hpp"
#include "Basic/String.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/File.hpp"
#include "Covariances/ECov.hpp"
#include "Covariances/CovContext.hpp"

#include <math.h>

/*********************/
/* Program principal */
/*********************/

int main(int /*argc*/, char */*argv*/[])
{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  AMesh *mesh,*meshb;
  Db    *dbin,*dbgrid;
  Model *model;
  CovContext ctxt;
  double  angle[3];
  double  range,param;
  MatrixRectangular apices;
  VectorInt meshes;
  char    triswitch[10];
  ShiftOpCs shiftop;
  VectorInt nx(3);
  VectorDouble rotmat;
  VectorDouble extendmin(3);
  VectorDouble extendmax(3);
  VectorDouble cellsize(3);
  VectorDouble dilate(3);
  double *loc_apices;
  int    *loc_meshes;

  /* Main options */

  int flag_grid = 0;
  int flag_mesh = 0;
  int ndim      = 3;
  int verbose   = 1;
  int variety   = 0;  // 0 for Euclidean; 1 for Spherical
  ASpaceObject::defineDefaultSpace(SPACE_RN, ndim);

  /* Cleverness of the options */

  if (variety == 1) flag_mesh = 0;
  if (flag_mesh)    flag_grid = 0;

  /* Initializations */

  extendmin[0] = 50.;
  extendmin[1] = 0.;
  extendmin[2] = 10.;
  extendmax[0] = 300.;
  extendmax[1] = 200.;
  extendmax[2] = 100.;
  cellsize[0]  = 10.;
  cellsize[1]  = 20.;
  cellsize[2]  = 10.;
  dilate[0]    = 0.;
  dilate[1]    = 0.;
  dilate[2]    = 0.;
  angle[0]     = 30.;
  angle[1]     = 0.;
  angle[2]     = 0.;
  dbin         = nullptr;
  dbgrid       = nullptr;
  loc_apices   = nullptr;
  loc_meshes   = nullptr;

  /* Triswitch option */

  if (variety == 0)
    (void) gslStrcpy((char *) triswitch,"Q");
  else
    (void) gslStrcpy((char *) triswitch,"-r1");
  
  /* Rotation definition */

  rotmat.resize(ndim * ndim);
  if (ndim == 2)
    ut_rotation_matrix_2D(angle[0],rotmat.data());
  else if (ndim == 3)
    ut_rotation_matrix_3D(angle[0],angle[1],angle[2],rotmat.data());

  /* Model definition */

  range        = 5.;
  param        = 1.;
  VectorDouble sill = { 2. };
  ctxt = CovContext(1, ndim);
  model = new Model(ctxt);
  if (model_add_cova(model,ECov::BESSEL_K,0,0,range,param,
                     VectorDouble(),VectorDouble(),sill,0.))
    messageAbort("Definition of the Model");

  /* Define the Grid Db */

  if (flag_grid)
  {
    nx[0] = 1+ceil((extendmax[0] - extendmin[0] - cellsize[0]/2.)/cellsize[0]);
    nx[1] = 1+ceil((extendmax[1] - extendmin[1] - cellsize[1]/2.)/cellsize[1]);
    dbgrid = db_create_grid(0,ndim,0,ELoadBy::COLUMN,1,nx,extendmin,cellsize);
  }

  /* Setup the license */

  if (setup_license("Demonstration")) return(0);

  /* Setup constants */

  OptDbg::reset();

  /* Instantiate the Meshing */

  mesh = MeshFactory::createMesh(variety,
                                 extendmin,extendmax,cellsize,rotmat,
                                 dilate,dbin,dbgrid,triswitch,
                                 apices,meshes,verbose);
  if (mesh == NULL) return(1);
  mesh->display();

  if (flag_mesh)
  {
    mesh->getElements(apices,meshes);
    meshb = MeshFactory::createMesh(variety,
                                    extendmin,extendmax,rotmat,cellsize,
                                    dilate,dbin,dbgrid,triswitch,
                                    apices,meshes,verbose);
    loc_apices = (double *) mem_free((char *) loc_apices);
    loc_meshes = (int    *) mem_free((char *) loc_meshes);
    if (meshb == NULL) return(1);
    meshb->display();
  }

  /* Instantiate the ShiftOp */

  if (! flag_mesh)
    shiftop.initFromMesh(mesh,model,NULL);
  else
    shiftop.initFromMesh(meshb,model,NULL,0,0,verbose);
  
  message("Test performed successfully\n");

  return(0);
}
