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
#include "geoslib_f.h"

#include "Basic/ASerializable.hpp"
#include "Basic/File.hpp"
#include "Basic/Law.hpp"
#include "Basic/String.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovAnisoList.hpp"
#include "Covariances/CovContext.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Variogram/Vario.hpp"
#include "Variogram/VarioParam.hpp"

using namespace gstlrn;

/****************************************************************************/
/*!
** Main Program
**
*****************************************************************************/
int main(int argc, char* argv[])
{
  Id nx               = 10;
  Id niter            = 10000;
  Id nburn            = 100;
  double range        = 10.;
  double bound        = TEST;
  bool flag_sym_neigh = true;

  Id seed               = 5452;
  Id ndim               = 2;
  Id nvar               = 1;
  Id nbsimu             = 1;
  double sill           = 1.;
  Id nlag               = 20;
  VectorDouble ranges   = {range, range};
  bool verbose          = true;
  bool flag_moving      = true;
  bool flag_propagation = false;
  bool flag_multi_mono  = false;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // Setup constants

  defineDefaultSpace(ESpaceType::RN, ndim);
  ASerializable::setPrefixName("test_AllGibbs-");
  law_set_random_seed(seed);

  // Data file

  VectorDouble dx = {1., 1.};
  DbGrid* db      = DbGrid::create({nx, nx}, dx);
  if (!FFFF(bound))
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

  CovContext ctxt(nvar, 2, 1.);
  auto* model = new Model(ctxt);
  CovAnisoList covs(ctxt);
  CovAniso cova(ECov::EXPONENTIAL, ctxt);
  cova.setRanges(ranges);
  cova.setSill(sill);
  covs.addCov(cova);
  model->setCovAnisoList(&covs);
  model->display();

  // Gibbs

  if (gibbs_sampler(db, model, nbsimu, seed, nburn, niter, flag_moving, false,
                    flag_multi_mono, flag_propagation, flag_sym_neigh, 2, 5.,
                    false, false, verbose)) return 1;
  DbStringFormat dbfmt(FLAG_STATS, {"*Gibbs*"});
  db->display(&dbfmt);
  (void)db->dumpToNF("Result");

  // Calculate a variogram on the samples

  VarioParam varioparam;
  std::vector<DirParam> dirparams = DirParam::createMultipleInSpace(nlag);
  varioparam.addMultiDirs(dirparams);
  Vario vario(varioparam);
  VectorString names = db->getName("Gibbs*");
  for (Id isimu = 0; isimu < nbsimu; isimu++)
  {
    db->clearLocators(ELoc::Z);
    db->setLocator(names[isimu], ELoc::Z, 0);
    vario.compute(db, ECalcVario::VARIOGRAM);
    (void)vario.dumpToNF(incrementStringVersion("Vario", isimu + 1));
  }

  // Cleaning structures

  delete db;
  delete model;
  return (0);
}
