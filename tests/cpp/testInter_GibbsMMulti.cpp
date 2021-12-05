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
#include "Neigh/Neigh.hpp"
#include "geoslib_d.h"
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "geoslib_define.h"

/****************************************************************************/
/*!
** Main Program
**
*****************************************************************************/
int main(int /*argc*/, char */*argv*/[])


{
  int iptr;
  bool flag_inter = true;

  int nx        = 30;
  int niter     = 100;
  int nburn     = 20;
  double range  = 10.;
  double bound  = TEST;
  double eps    = 0.05;
  bool storeTables = true;
  bool storeInternal = true;
  bool storeVario = true;

  if (flag_inter)
  {
    nx    = askInt("Number of grid mesh [in each direction]", nx);
    niter = askInt("Number of Gibbs iterations",niter);
    nburn = askInt("Number of burning steps",nburn);
    eps   = askDouble("Epsilon",eps);
    range = askDouble("Isotropic Range",range);
    storeInternal = askBool("Store Internal", storeInternal);
  }

  int seed     = 5452;
  int ndim     = 2;
  int nvar     = 1;
  int nbsimu   = 1;
  double sill  = 1.;
  int nlag     = 20;

  VectorDouble ranges = { range, range};
  bool verbose          = true;

  // Setup constants

  ASpaceObject::defineDefaultSpace(SPACE_RN, ndim);
  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("AllGibbs-");
  law_set_random_seed(seed);
  
  // Data file

  Db* db = new Db({nx,nx},{1.,1.});
  if (! FFFF(bound))
  {
    db->addFields(1, -bound, "Lower", ELoc::L);
    db->addFields(1, +bound, "Upper", ELoc::U);
  }
  else
  {
    db->addFields(1, TEST, "Lower", ELoc::L);
    db->addFields(1, TEST, "Upper", ELoc::U);
  }
  if (db_locator_attribute_add(db,ELoc::GAUSFAC,nbsimu*nvar,0,0.,&iptr)) return 1;

  // Model

  CovContext ctxt(nvar,2,1.); // use default space
  Model* model = new Model(ctxt);
  CovAniso cova(ECov::SPHERICAL,ctxt);
  cova.setRanges(ranges);
  cova.setSill({sill});
  model->addCova(&cova);
  model->display();

  // Initialize Gibbs

  GibbsMMulti gibbs(db, model);
  gibbs.setOptionStats(2);
  gibbs.setEps(eps);
  gibbs.setFlagStoreInternal(storeInternal);
  gibbs.setStoreTables(storeTables);
  gibbs.init(1, nvar, nburn, niter,0, true);

  // Allocate the Gaussian vector

  VectorVectorDouble y = gibbs.allocY();

  /* Allocate the covariance matrix inverted */

  if (gibbs.covmatAlloc(verbose)) return 1;

  // Invoke the Gibbs calculator

  for (int isimu = 0; isimu < nbsimu; isimu++)
    if (gibbs.run(y, 0, isimu, verbose, false)) return 1;
  db->serialize("Result");

  // Calculate a variogram on the samples

  if (storeVario)
  {
    VarioParam varioparam;
    std::vector<DirParam> dirparams = generateMultipleGridDirs(ndim, nlag);
    varioparam.addDirs(dirparams);
    VectorString names = db->getNames("gausfac*");
    for (int isimu = 0; isimu < (int) names.size(); isimu++)
    {
      db->clearLocators(ELoc::Z);
      db->setLocator(names[isimu], ELoc::Z);
      Vario vario(&varioparam, db);
      vario.compute("vg", true);
      vario.serialize(incrementStringVersion("Vario", isimu + 1));
    }
  }

  // Cleaning structures

  gibbs.cleanup();
  db    = db_delete(db);
  model = model_free(model);
  return(0);
}
