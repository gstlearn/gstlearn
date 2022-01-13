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
#include "Basic/Law.hpp"
#include "Basic/String.hpp"
#include "Space/Space.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Basic/AException.hpp"
#include "Basic/ASerializable.hpp"
#include "Model/Model.hpp"
#include "Variogram/VarioParam.hpp"
#include "Variogram/Vario.hpp"
#include "Db/Db.hpp"
#include "Neigh/Neigh.hpp"

/****************************************************************************/
/*!
** Main Program
**
*****************************************************************************/
int main(int /*argc*/, char * /*argv*/[])

{
  int error = 1;

  int nx        = 10;
  int niter     = 10000;
  int nburn     = 100;
  int nmaxi     = 4;
  double range  = 10.;
  double bound  = TEST;
  bool flag_sym_neigh = true;

  int seed     = 5452;
  int ndim     = 2;
  int nvar     = 1;
  int nbsimu   = 1;
  double sill  = 1.;
  int nlag     = 20;
  double nbgh_radius = 10. * range;
  VectorDouble ranges = { range, range};
  bool verbose          = true;
  bool flag_moving      = true;
  bool flag_propagation = false;
  bool flag_multi_mono  = false;

  // Setup constants

  ASpaceObject::defineDefaultSpace(SPACE_RN, ndim);
  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("AllGibbs-");
  law_set_random_seed(seed);
  
  // Data file

  VectorDouble dx = {1., 1.};
  Db* db = Db::createFromGrid({nx,nx},dx);
  if (! FFFF(bound))
  {
    db->addFieldsByConstant(1, -bound, "Bounds", ELoc::L);
    db->addFieldsByConstant(1, +bound, "Bounds", ELoc::U);
  }
  else
  {
    db->addFieldsByConstant(1, TEST, "Bounds", ELoc::L);
    db->addFieldsByConstant(1, TEST, "Bounds", ELoc::U);
  }

  // Model

  CovContext ctxt(nvar,2,1.);
  Model* model = new Model(ctxt);
  CovLMC covs(ctxt.getSpace());
  CovAniso cova(ECov::EXPONENTIAL,ctxt);
  cova.setRanges(ranges);
  cova.setSill({sill});
  covs.addCov(&cova);
  model->setCovList(&covs);
  model->display();

  // Neighborhood

  Neigh* neigh = nullptr;
  if (flag_moving)
  {
    neigh = Neigh::createMoving(ndim, nmaxi, nbgh_radius);
    neigh->display();
  }

  // Gibbs

  error = gibbs_sampler(db, model, neigh, nbsimu, seed, nburn, niter, false,
                        flag_multi_mono, flag_propagation,
                        flag_sym_neigh, 2,
                        5., false, false, verbose);
  if (error) return 1;
  db->serialize("Result");

  // Calculate a variogram on the samples

  VarioParam varioparam;
  std::vector<DirParam> dirparams = generateMultipleGridDirs(ndim, nlag);
  varioparam.addMultiDirs(dirparams);
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
