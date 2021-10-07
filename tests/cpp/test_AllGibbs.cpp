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
#include "Basic/Law.hpp"
#include "Space/Space.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/CovAniso.hpp"
#include "geoslib_d.h"
#include "geoslib_f.h"

/****************************************************************************/
/*!
** Main Program
**
*****************************************************************************/
int main(int argc, char *argv[])

{
  int seed     = 5452;
  int ndim     = 2;
  int nvar     = 1;
  int nbsimu   = 1;
  int nlag     = 20;

  int niter    = 10000;
  int nburn    = 100;

  int nbgh_maxi = 20;
  double nbgh_radius = 1000.;

  VectorDouble ranges = { 10., 10.};
  double sill = 1.;

  bool verbose          = true;
  bool flag_moving      = true;
  bool flag_propagation = false;
  bool flag_multi_mono  = false;

  // Setup constants

  ASpaceObject::createGlobalSpace(SPACE_RN, ndim);
  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("AllGibbs-");
  law_set_random_seed(seed);
  
  // Data file

  VectorDouble dx = {1., 1.};
  Db* db = new Db({50,50},dx);
  db->addFields(1, -4, "Bounds", LOC_L);
  db->addFields(1,  4, "Bounds", LOC_U);
  
  // Model

  CovContext ctxt(nvar,2,1.);
  Model* model = new Model(ctxt);
  CovAniso cova(ECov::EXPONENTIAL,ctxt);
  cova.setRanges(ranges);
  cova.setSill({sill});
  model->addCova(&cova);
  model->display();

  // Neighborhood

  Neigh* neigh = nullptr;
  if (flag_moving)
  {
    neigh = new Neigh(ndim, nbgh_maxi, nbgh_radius);
    neigh->display();
  }

  // Gibbs

  gibbs_sampler(db, model, neigh, nbsimu, seed, nburn, niter, false,
                flag_multi_mono, flag_propagation, 2,
                5., EPSILON3, false, false, verbose);
  db->displayMore(FLAG_STATS);
  db->serialize("Result");

  // Calculate a variogram on the samples

  VarioParam varioparam;
  std::vector<DirParam> dirparams = generateMultipleGridDirs(ndim, nlag);
  varioparam.addDirs(dirparams);
  Vario vario(&varioparam,db);
  VectorString names = db->getNames("Gibbs*");
  for (int isimu=0; isimu<nbsimu; isimu++)
  {
    db->clearLocators(LOC_Z);
    db->setLocator(names[isimu],LOC_Z);
    vario.compute("vg",true);
    vario.serialize(incrementStringVersion("Vario",isimu+1));
  }

  // Cleaning structures

  db    = db_delete(db);
  model = model_free(model);
  neigh = neigh_free(neigh);
  return(0);
}
