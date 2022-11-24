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

#include "Enum/ECov.hpp"
#include "Enum/ELoadBy.hpp"

#include "Model/Model.hpp"
#include "Mesh/AMesh.hpp"
#include "Mesh/MeshEStandard.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Mesh/MeshSpherical.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixInt.hpp"
#include "Mesh/MeshFactory.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "Space/ASpaceObject.hpp"
#include "Basic/String.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/File.hpp"
#include "Basic/ASerializable.hpp"
#include "Covariances/CovContext.hpp"

#include <math.h>


/*********************/
/* Program principal */
/*********************/

int main(int /*argc*/, char */*argv*/[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("TestMesh-");

  VectorInt nx(3);
  VectorDouble extendmin = {50., 0., 10.};
  VectorDouble extendmax = { 300., 200., 100.};
  VectorDouble cellsize = { 10., 20., 30.};
  VectorDouble dilate = {0., 0., 0.};
  VectorDouble angle = {30., 0., 0.};

  /* Main options */

  int flag_grid = 0;
  int flag_mesh = 0;
  int ndim      = 3;
  int verbose   = 1;
  int variety   = 0;  // 0 for Euclidean; 1 for Spherical
  if (variety == 0)
    defineDefaultSpace(ESpaceType::RN, ndim);
  else
    defineDefaultSpace(ESpaceType::SN, ndim, EARTH_RADIUS);

  /* Cleverness of the options */

  if (variety == 1) flag_mesh = 0;
  if (flag_mesh)    flag_grid = 0;

  /* Initializations */

  Db* dbin       = nullptr;
  DbGrid* dbgrid = nullptr;

  /* Triswitch option */

  char triswitch[10];
  if (variety == 0)
    (void) gslStrcpy((char *) triswitch,"Q");
  else
    (void) gslStrcpy((char *) triswitch,"-r1");
  
  /* Rotation definition */

  VectorDouble rotmat(ndim * ndim);
  if (ndim == 2)
    GH::rotationInit(angle[0],rotmat.data());
  else if (ndim == 3)
    GH::rotationInit(angle[0],angle[1],angle[2],rotmat.data());

  /* Model definition */

  double range = 5.;
  double param = 1.;
  double sill  = 2.;
  Model* model = Model::createFromParam(ECov::BESSEL_K, range, sill, param);

  /* Define the Grid Db */

  if (flag_grid)
  {
    nx[0] = 1+(int) ceil((extendmax[0] - extendmin[0] - cellsize[0]/2.)/cellsize[0]);
    nx[1] = 1+(int) ceil((extendmax[1] - extendmin[1] - cellsize[1]/2.)/cellsize[1]);
    dbgrid = DbGrid::create(nx, cellsize, extendmin, VectorDouble(),
                            ELoadBy::COLUMN, VectorDouble(),
                            VectorString(), VectorString(), 1);

  }

  /* Instantiate the Meshing */

  MatrixRectangular apices;
  MatrixInt meshes;
  AMesh* mesh = MeshFactory::createMesh(variety,
                                        extendmin,extendmax,cellsize,rotmat,
                                        dilate,dbin,dbgrid,triswitch,
                                        apices,meshes,verbose);
  if (mesh == nullptr) return(1);
  mesh->display();
  (void) mesh->dumpToNF("Standard.ascii");

  AMesh* meshb = nullptr;
  if (flag_mesh)
  {
    mesh->getElements(apices,meshes);
    meshb = MeshFactory::createMesh(variety,
                                    extendmin, extendmax, cellsize,rotmat,
                                    dilate, dbin, dbgrid, triswitch,
                                    apices, meshes, verbose);
    if (meshb == nullptr) return(1);
    meshb->display();
    (void) meshb->dumpToNF("Standard.bis.ascii");
  }

  /* Instantiate the ShiftOp */

  ShiftOpCs shiftop;
  if (! flag_mesh)
    shiftop.initFromMesh(mesh,model,NULL);
  else
    shiftop.initFromMesh(meshb,model,NULL,0,0,verbose);
  
  message("Test performed successfully\n");

  return(0);
}
