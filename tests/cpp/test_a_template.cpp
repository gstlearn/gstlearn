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
/**
 * This file is meant to parametrized the ModelGeneric in terms of ParamInfo
 * and to fit the values of these parameters according to the Maximum LogLikelihood
 * method and using the Vecchia approximation.
 */
#include "Basic/OptDbg.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Enum/ESpaceType.hpp"
#include "Estimation/CalcGlobal.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Space/ASpaceObject.hpp"
#include "Variogram/DirParam.hpp"
#include "Variogram/Vario.hpp"
#include "geoslib_define.h"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  DECLARE_UNUSED(argc);
  DECLARE_UNUSED(argv);

  defineDefaultSpace(ESpaceType::RN, 2);

  String filename = "/home/drenard/Rtest/Bez/myDb.NF";
  bool hasDrift   = true;
  bool flagLocal  = false;
  double propSel  = 0.1;
  Id optDebug     = 0;

  Db* myDb = Db::createFromNF(filename);
  myDb->setLocators({"longitude", "latitude"}, ELoc::X);
  myDb->setLocators({"density", "depth"}, ELoc::Z);
  myDb->deleteColumnsByLocator(ELoc::SEL);
  if (propSel > 0)
    myDb->addSelectionRandom(propSel);
  else
    myDb->clearSelection();

  DbStringFormat* dbfmt = DbStringFormat::createFromFlags(false, false, false, true, true);
  myDb->display(dbfmt);
  message("Number of active samples = %d\n", myDb->getNSampleActive());

  DbGrid* myGrid = DbGrid::createCoveringDb(myDb, VectorInt(), {10, 10});
  message("Number of grid nodes = %d\n", myGrid->getNSample());

  NeighUnique myNeighU = NeighUnique();

  VarioParam* varioparam = VarioParam::createOmniDirection(10, 50);
  Vario* myVarioOmni     = Vario::computeFromDb(*varioparam, myDb);

  Model myModelBi;
  myModelBi.fit(myVarioOmni, {ECov::NUGGET, ECov::EXPONENTIAL, ECov::SPHERICAL});
  if (hasDrift) myModelBi.setDriftIRF(0, 0);

  OptDbg::setReference(optDebug);
  Id mode = 0;

  if (mode == 0 || mode == 2)
  {
    // CoKriging
    myDb->setLocators({"density", "depth"}, ELoc::Z);
    if (flagLocal)
    {
      (void)kriging(myDb, myGrid, &myModelBi, &myNeighU, true, true, false, KrigOpt(),
                    NamingConvention("CoK"));
    }
    else
    {
      Global_Result resBivariate = global_kriging(myDb, myGrid, &myModelBi);
      resBivariate.display();
    }
  }

  if (mode == 0 || mode == 1)
  {
    // Kriging of first variable

    myDb->setLocator("density", ELoc::Z, 0);
    Model* myModelUni = myModelBi.createReduce({0});
    if (flagLocal)
    {
      (void)kriging(myDb, myGrid, myModelUni, &myNeighU, true, true, false, KrigOpt(),
                    NamingConvention("K"));
    }
    else
    {
      Global_Result resUnivariate = global_kriging(myDb, myGrid, myModelUni);
      resUnivariate.display();
    }
    delete myModelUni;
  }
}
