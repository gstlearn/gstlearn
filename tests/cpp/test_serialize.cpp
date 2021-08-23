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
#include "Db/Db.hpp"
#include "geoslib_d.h"
#include "geoslib_f.h"

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix algebra
**
*****************************************************************************/
int main(int argc, char *argv[])

{
  bool verbose = true;
  auto pygst = std::string(std::getenv("PYGSTLEARN_DIR"));

  // ===== Create the Db db1
  int nech = 20;
  int ndim = 2;
  Db db1(nech,VectorDouble(),VectorDouble(),ndim);
  VectorDouble vec = ut_vector_simulate_gaussian(nech);
  db1.addFields(vec,"myvar",LOC_Z);
  db1.display();
  
  // Serialize db1
  db1.serialize(pygst + "Neutral.Db.ascii",verbose);

  // Deserialize db2
  Db db2(pygst + "Neutral.Db.ascii",verbose);
  db2.display();

  // ===== Create the Grid Db
  Db dbg1({12,10},{0.1,0.3},{0.2,0.4});
  vec = ut_vector_simulate_gaussian(dbg1.getSampleNumber());
  dbg1.addFields(vec,"myvar",LOC_Z);
  dbg1.display();

  // Serialize dbg1
  dbg1.serialize(pygst + "Neutral.Dbg.ascii",verbose);

  // Deserialize dbg2
  Db dbg2(pygst + "Neutral.Dbg.ascii",verbose);
  dbg2.display();

  // ===== Create the Polygon poly1
  Polygons poly1(&db1);
  poly1.display();

  // Serialize poly1
  poly1.serialize(pygst + "Neutral.Polygon.ascii",verbose);

  // Deserialize poly2
  Polygons poly2(pygst + "Neutral.Polygon.ascii",verbose);
  poly2.display();

  // ===== Compute an experimental variogram
  Vario vario1 = Vario();
  Dir dir(2, 10, 0.02);
  vario1.addDirs(dir);
  vario1.compute(&db1,"vg");
  vario1.display();

  // Serialize vario1
  vario1.serialize(pygst + "Neutral.Vario.ascii",verbose);

  // Deserialize vario2
  Vario vario2(pygst + "Neutral.Vario.ascii",verbose);
  vario2.display();

  // ===== Compute a Model
  Model model1(&db1);
  CovAniso cova(COV_EXPONENTIAL, 0.3, 1., 0.2, model1.getContext());
  model1.addCova(&cova);
  model1.display();

  // Serialize model1
  model1.serialize(pygst + "Neutral.Model.ascii",verbose);

  // Deserialize model2
  Model model2(pygst + "Neutral.Model.ascii",verbose);
  model2.display();

  return(0);
}
