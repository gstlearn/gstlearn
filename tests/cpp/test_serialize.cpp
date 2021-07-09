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

  // ===== Create the Db db1
  int nech = 20;
  int ndim = 2;
  Db db1 = Db(nech,VectorDouble(),VectorDouble(),ndim);
  VectorDouble vec = ut_vector_simulate_gaussian(nech);
  db1.addFields(vec,"myvar",LOC_Z);
  db1.display();
  
  // Serialize db1
  db1.serialize("Db.Neutral",verbose);

  // Deserialize db2
  Db db2("Db.Neutral",verbose);
  db2.display();

  // ===== Create the Polygon poly1
  Polygons poly1 = Polygons(&db1);
  poly1.display();

  // Serialize poly1
  poly1.serialize("Polygon.Neutral",verbose);

  // Deserialize poly2
  Polygons poly2("Polygon.Neutral",verbose);
  poly2.display();

  // ===== Compute an experimental variogram
  Vario vario1 = Vario("vg");
  Dir dir = Dir(2, 10, 0.02);
  vario1.addDirs(dir);
  vario1.compute(&db1);
  vario1.display();

  // Serialize vario1
  vario1.serialize("Vario.Neutral",verbose);

  // Deserialize vario2
  Vario vario2("Vario.Neutral",verbose);
  vario2.display();

  // ===== Compute a Model
  Model model1 = Model(&db1);
  CovAniso cova = CovAniso(COV_EXPONENTIAL, 0.3, 1., 0.2, model1.getContext());
  model1.addCova(&cova);
  model1.display();

  // Serialize model1
  model1.serialize("Model.Neutral",verbose);

  // Unserialize model2
  Model model2("Model.Neutral",verbose);
  model2.display();

  return(0);
}
