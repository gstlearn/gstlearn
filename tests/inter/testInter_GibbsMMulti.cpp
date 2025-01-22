/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "geoslib_old_f.h"
#include "geoslib_define.h"

#include "Enum/ESpaceType.hpp"

#include "Basic/Law.hpp"
#include "Basic/ASerializable.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovAnisoList.hpp"
#include "Model/Model.hpp"
#include "Variogram/VarioParam.hpp"
#include "Variogram/Vario.hpp"
#include "Gibbs/GibbsMMulti.hpp"
#include "Gibbs/GibbsMulti.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"

/****************************************************************************/
/*!
** Main Program
**
*****************************************************************************/
int main()
{
  int iptr;
  bool flag_inter = false;

  int nx        = 50;
  int niter     = 1000;
  int nburn     = 20;
  double range  = 10.;
  double bound  = TEST;
  double eps    = EPSILON6;
  bool storeInternal = false; // No HDF5 by default
  bool storeVario = false;

  if (flag_inter)
  {
    nx    = askInt("Number of grid mesh [in each direction]", nx);
    niter = askInt("Number of Gibbs iterations",niter);
    nburn = askInt("Number of burning steps",nburn);
    eps   = askDouble("Epsilon",eps);
    range = askDouble("Isotropic Range",range);
#ifdef _USE_HDF5
    storeInternal = askBool("Store Internal", storeInternal);
#endif
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

  defineDefaultSpace(ESpaceType::RN, ndim);
  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("AllGibbs-");
  law_set_random_seed(seed);
  
  // Data file

  DbGrid* db = DbGrid::create({nx,nx},{1.,1.});
  if (! FFFF(bound))
  {
    db->addColumnsByConstant(1, -bound, "Lower", ELoc::L);
    db->addColumnsByConstant(1, +bound, "Upper", ELoc::U);
  }
  else
  {
    db->addColumnsByConstant(1, TEST, "Lower", ELoc::L);
    db->addColumnsByConstant(1, TEST, "Upper", ELoc::U);
  }
  if (db_locator_attribute_add(db, ELoc::GAUSFAC, nbsimu * nvar, 0, 0., &iptr)) return 1;

  // Model

  CovContext ctxt(nvar,2,1.); // use default space
  Model model(ctxt);
  CovAnisoList covs(ctxt.getSpace());
  CovAniso cova(ECov::SPHERICAL,ctxt);
  cova.setRanges(ranges);
  cova.setSill(sill);
  covs.addCov(&cova);
  model.setCovAnisoList(&covs);
  model.display();

  // Initialize Gibbs

  GibbsMMulti gibbs(db, &model);
  gibbs.setOptionStats(2);
  gibbs.setEps(eps);
  gibbs.setFlagStoreInternal(storeInternal);
  gibbs.init(1, nvar, nburn, niter,0, true);

  // Allocate the Gaussian vector

  VectorVectorDouble y = gibbs.allocY();

  /* Allocate the covariance matrix inverted */

  if (gibbs.covmatAlloc(verbose)) return 1;

  // Invoke the Gibbs calculator

  for (int isimu = 0; isimu < nbsimu; isimu++)
    if (gibbs.run(y, 0, isimu, verbose, false)) return 1;
  (void) db->dumpToNF("Result");

  // Calculate a variogram on the samples

  if (storeVario)
  {
    VarioParam varioparam;
    std::vector<DirParam> dirparams = DirParam::createMultipleInSpace(nlag);
    varioparam.addMultiDirs(dirparams);
    VectorString names = db->getName("gausfac*");
    for (int isimu = 0; isimu < (int) names.size(); isimu++)
    {
      db->clearLocators(ELoc::Z);
      db->setLocator(names[isimu], ELoc::Z, 0);
      Vario vario(varioparam);
      vario.compute(db, ECalcVario::VARIOGRAM);
      (void) vario.dumpToNF(incrementStringVersion("Vario", isimu + 1));
    }
  }

  // Cleaning structures

  gibbs.cleanup();
  delete db;
  return(0);
}
