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
#include "Basic/File.hpp"
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
#include "Db/DbStringFormat.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Neigh/NeighMoving.hpp"

/****************************************************************************/
/*!
** Main Program
**
*****************************************************************************/
int main(int /*argc*/, char * /*argv*/[])

{
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
  bool verbose          = false;
  bool flag_moving      = true;
  bool flag_propagation = false;
  bool flag_multi_mono  = false;

  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  // Setup constants

  ASpaceObject::defineDefaultSpace(SPACE_RN, ndim);
  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("AllGibbs-");
  law_set_random_seed(seed);
  
  // Data file

  VectorDouble dx = {1., 1.};
  DbGrid* db = DbGrid::create({nx,nx},dx);
  if (! FFFF(bound))
  {
    db->addColumnsByConstant(1, -bound, "Bounds", ELoc::L);
    db->addColumnsByConstant(1, +bound, "Bounds", ELoc::U);
  }
  else
  {
    db->addColumnsByConstant(1, TEST, "Bounds", ELoc::L);
    db->addColumnsByConstant(1, TEST, "Bounds", ELoc::U);
  }
  db->display();

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

  ANeighParam* neighparam = nullptr;
  if (flag_moving)
  {
    neighparam = NeighMoving::create(ndim, false, nmaxi, nbgh_radius);
    neighparam->display();
  }

  // Gibbs

  if (gibbs_sampler(db, model, neighparam, nbsimu, seed, nburn, niter, false,
                    flag_multi_mono, flag_propagation, flag_sym_neigh, 2, 5.,
                    false, false, verbose)) return 1;
  DbStringFormat dbfmt(FLAG_STATS,{"*Gibbs*"});
  db->display(&dbfmt);
  (void) db->dumpToNF("Result");

  // Calculate a variogram on the samples

  VarioParam varioparam;
  std::vector<DirParam> dirparams = DirParam::createMultipleFromGrid(ndim, nlag);
  varioparam.addMultiDirs(dirparams);
  Vario vario(&varioparam,db);
  VectorString names = db->getNames("Gibbs*");
  for (int isimu=0; isimu<nbsimu; isimu++)
  {
    db->clearLocators(ELoc::Z);
    db->setLocator(names[isimu],ELoc::Z);
    vario.computeByKey("vg");
    (void) vario.dumpToNF(incrementStringVersion("Vario",isimu+1));
  }

  // Cleaning structures

  db    = db_delete(db);
  model = model_free(model);
  delete neighparam;
  return(0);
}
