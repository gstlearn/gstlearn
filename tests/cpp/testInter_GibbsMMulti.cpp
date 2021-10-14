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
#include "Basic/AException.hpp"
#include "Basic/ASerializable.hpp"
#include "Model/Model.hpp"
#include "Variogram/VarioParam.hpp"
#include "Variogram/Vario.hpp"
#include "Gibbs/GibbsMMulti.hpp"
#include "Db/Db.hpp"
#include "geoslib_d.h"
#include "geoslib_f.h"

/****************************************************************************/
/*!
** Main Program
**
*****************************************************************************/
int main(int argc, char *argv[])

{
  bool flag_inter = true;

  int nx        = 10;
  int niter     = 10000;
  int nburn     = 100;
  int nmaxi     = 4;
  double range  = 10.;
  double bound  = TEST;
  bool flag_sym_neigh = true;
  bool flag_sym_Q = true;

  if (flag_inter)
  {
    nx = askInt("Number of grid mesh [in each direction]", nx);
    niter = askInt("Number of Gibbs iterations",niter);
    nburn = askInt("Number of burning steps",nburn);
    nmaxi = askInt("Number of samples in Neighborhood",nmaxi);
    range = askDouble("Isotropic Range",range);
    bound = askDouble("Bounds [None: -10]",-10.);
    if (bound <= -10.) bound = TEST;
    flag_sym_neigh = askBool("Symmetrization of Neighborhood",flag_sym_neigh);
    flag_sym_Q = askBool("Symmetrization of Q",flag_sym_Q);
  }

  int seed     = 5452;
  int ndim     = 2;
  int nvar     = 1;
  int nbsimu   = 1;
  double sill  = 1.;
  int nlag     = 20;
  double nbgh_radius = 10. * range;
  VectorDouble ranges = { range, range};
  bool verbose          = true;

  // Setup constants

  ASpaceObject::defineDefaultSpace(SPACE_RN, ndim);
  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("AllGibbs-");
  law_set_random_seed(seed);
  
  // Data file

  VectorDouble dx = {1., 1.};
  Db* db = new Db({nx,nx},dx);
  if (! FFFF(bound))
  {
    db->addFields(1, -bound, "Bounds", ELoc::L);
    db->addFields(1, +bound, "Bounds", ELoc::U);
  }
  else
  {
    db->addFields(1, TEST, "Bounds", ELoc::L);
    db->addFields(1, TEST, "Bounds", ELoc::U);
  }

  // Model

  CovContext ctxt(nvar,2,1.); // use default space
  Model* model = new Model(ctxt);
  CovAniso cova(ECov::EXPONENTIAL,ctxt);
  cova.setRanges(ranges);
  cova.setSill({sill});
  model->addCova(&cova);
  model->display();

  // Neighborhood

  Neigh* neigh = new Neigh(ndim, nmaxi, nbgh_radius);
  neigh->display();

  // Initialize Gibbs

  GibbsMMulti gibbs(db, model, neigh);
  gibbs.setOptionStats(true);
  gibbs.setFlagSymNeigh(flag_sym_neigh);
  gibbs.setFlagSymQ(flag_sym_neigh);
  gibbs.init(1, nvar, nburn, niter,0, false, true);

  /* Allocate the covariance matrix inverted */

  if (gibbs.covmatAlloc(verbose)) return 1;

  // Invoke the Gibbs calculator

  if (gibbs.run(nbsimu, verbose)) return 1;
  db->serialize("Result");

  // Calculate a variogram on the samples

  VarioParam varioparam;
  std::vector<DirParam> dirparams = generateMultipleGridDirs(ndim, nlag);
  varioparam.addDirs(dirparams);
  Vario vario(&varioparam,db);
  VectorString names = db->getNames("Gibbs*");
  for (int isimu=0; isimu<nbsimu; isimu++)
  {
    db->clearLocators(ELoc::Z);
    db->setLocator(names[isimu],ELoc::Z);
    vario.compute("vg",true);
    vario.serialize(incrementStringVersion("Vario",isimu+1));
  }

  // Cleaning structures

  db    = db_delete(db);
  model = model_free(model);
  neigh = neigh_free(neigh);
  return(0);
}
