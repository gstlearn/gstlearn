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
#include "geoslib_f.h"
#include "csparse_f.h"

#include "Db/Db.hpp"
#include "Basic/File.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Variogram/DirParam.hpp"
#include "Variogram/VarioParam.hpp"
#include "Variogram/Vario.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Geometry/Geometry.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "Stats/Classical.hpp"

/**
 * This file is meant to perform any test that needs to be coded for a quick trial
 */
int main(int /*argc*/, char */*argv*/[])

{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
//  StdoutRedirect sr(sfn.str());

  ASpaceObject::defineDefaultSpace(ESpaceType::SPACE_RN, 2);

  Db* dat = Db::createFromNF("/home/drenard/a/dat.ascii");
  dat->display();

  Model* model = Model::createFromNF("/home/drenard/a/model.ascii");
  model->display();

  DbGrid* Result = DbGrid::createFromNF("/home/drenard/a/Result.ascii");
  Result->display();

  VectorDouble propGlob = dbStatisticsFacies(dat);
  int ncat = (int) propGlob.size();
  for (int i = 0; i < ncat; i++)
    message("Proportion of facies %d = %lf\n",i+1,propGlob[i]);

  (void) db_proportion_estimate(dat,Result,model);

//  MeshETurbo* mesh = MeshETurbo::createFromGrid(grid);
//  mesh->display();
//
//  ShiftOpCs S = ShiftOpCs(mesh, model, grid);
//  ut_vector_display_range("TildeC", S.getTildeC());
//  ut_vector_display_range("Lambda", S.getLambdas());
//  cs_print_dim("S", S.getS());
//  cs_print_range("S", S.getS());
//
//  PrecisionOp Qprop = PrecisionOp(&S, model->getCova(0), EPowerPT::ONE);
//  ProjMatrix Aproj = ProjMatrix(dat, mesh);
//  cs_print_dim("AProj", Aproj.getAproj());
//  cs_print_range("AProj", Aproj.getAproj());

  return (0);
}
