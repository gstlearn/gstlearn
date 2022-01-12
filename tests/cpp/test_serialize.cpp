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
#include "Db/Db.hpp"
#include "Variogram/VarioParam.hpp"
#include "Variogram/Vario.hpp"
#include "Model/Model.hpp"
#include "Basic/Table.hpp"
#include "Neigh/Neigh.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Polygon/Polygons.hpp"

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix algebra
**
*****************************************************************************/
int main(int /*argc*/, char */*argv*/[])

{
  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("TS-");

  // ===== Create the Db db1
  int nech = 20;
  int ndim = 2;
  bool verbose = false;

  Db db1(nech,VectorDouble(),VectorDouble(),ndim);
  VectorDouble vec1 = ut_vector_simulate_gaussian(nech);
  db1.addFields(vec1,"myvar1",ELoc::Z, 0);
  db1.display();
  
  // Serialize db1
  db1.serialize("Neutral.Db.ascii",verbose);

  // Deserialize db2
  Db db2("Neutral.Db.ascii",verbose);

  // ===== Create the Grid Db
  Db dbg1({12,10},{0.1,0.3},{0.2,0.4});
  vec1 = ut_vector_simulate_gaussian(dbg1.getSampleNumber());
  dbg1.addFields(vec1,"myvar1",ELoc::Z, 0);
  VectorDouble vec2 = ut_vector_simulate_gaussian(dbg1.getSampleNumber());
  vec2[2] = TEST;
  vec2[5] = TEST;
  dbg1.addFields(vec2,"myvar2",ELoc::Z, 1);
  dbg1.display();

  // Serialize dbg1
  dbg1.serialize("Neutral.Dbg.ascii",verbose);

  // Deserialize dbg2
  Db dbg2("Neutral.Dbg.ascii",verbose);
  dbg2.display();

  // ===== Create the Polygon poly1
  Polygons poly1;
  poly1.resetFromDb(&db1);
  Polygons polyb;
  polyb.resetFromDb(&dbg1);
  poly1.addPolySet(polyb.getPolySet(0));
  poly1.display();

  // Serialize poly1
  poly1.serialize("Neutral.Polygon.ascii",verbose);

  // Deserialize poly2
  Polygons* poly2 = Polygons::createFromNF("Neutral.Polygon.ascii",verbose);
  poly2->display();
  delete poly2;

  Polygons poly3;
  poly3.deSerialize("Neutral.Polygon.ascii", verbose);
  poly3.display();

  // ===== Compute an experimental variogram
  VarioParam varioparam1;
  DirParam dirparam(2, 10, 0.02);
  varioparam1.addDirs(dirparam);
  Vario vario1 = Vario(&varioparam1,&db1);
  vario1.compute("vg");
  vario1.display();

  // Serialize vario1
  vario1.serialize("Neutral.Vario.ascii",verbose);

  // Deserialize vario2
  Vario* vario2 = Vario::createFromNF("Neutral.Vario.ascii",verbose);
  vario2->display();

  // ===== Create a Model
  db1.display();
  Model model1(&db1);
  CovContext ctxt = model1.getContext();
  CovLMC covs(ctxt.getSpace());
  CovAniso cova(ECov::EXPONENTIAL, 0.3, 1., 0.2, ctxt);
  covs.addCov(&cova);
  model1.setCovList(&covs);
  model1.display();

  // Serialize model1
  model1.serialize("Neutral.Model.ascii",verbose);

  // Deserialize model2
  Model* model2 = Model::createFromNF("Neutral.Model.ascii",verbose);
  model2->display();

  // ===== Create a Table
  VectorVectorDouble table;
  int ncols = 3;
  int nrows = 10;
  table.resize(ncols);
  for (int icol = 0; icol < ncols; icol++)
    table[icol] = ut_vector_simulate_uniform(nrows);
  Table* table1 = Table::createFromArray(table);
  table1->display();

  // Serialize table
  table1->serialize("Neutral.Table.ascii",verbose);

  // Deserialize table1
  Table* table2 = Table::createFromNF("Neutral.Table.ascii",verbose);
  table2->display();

  delete vario2;
  delete model2;
  delete table1;
  delete table2;
  return(0);
}
