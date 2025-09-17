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
#include "Basic/AStringable.hpp"
#include "Enum/EFormatNF.hpp"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

#include "Basic/File.hpp"
#include "Core/Ascii.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Estimation/CalcKrigingGradient.hpp"
#include "Neigh/ANeigh.hpp"
#include "Neigh/NeighBench.hpp"
#include "Neigh/NeighImage.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Space/ASpaceObject.hpp"
#include "Variogram/Vario.hpp"

using namespace gstlrn;

/*********************/
/* Program principal */
/*********************/

int main(int argc, char* argv[])
{
  String filename;
  Db* dbin;
  DbGrid* dbout;
  Vario* vario;
  Model* model;
  ANeigh* neigh;
  Constraints constraints;
  DbStringFormat dbfmt;
  Id nbsimu, seed, nbtuba;
  static int nboot    = 10;
  static int niter    = 10;
  static bool verbose = false;

  /* Initializations */

  dbin  = nullptr;
  dbout = nullptr;
  vario = nullptr;
  model = nullptr;
  neigh = nullptr;

  /* Standard output redirection to file */

  StdoutRedirect sr("Result.out", argc, argv, 2);

  /* Create the output name (for storage of dump files) */

  VectorString subparts = separateKeywords(argv[1]);

  int nargs      = static_cast<int>(subparts.size());
  String outname = concatenateStrings("", subparts[nargs - 2], subparts[nargs - 1], "-");
  // if (outname == "Jeu3-") verbose = true; // Pour voir le resultat de Jeu3 en particulier
  ASerializable::setPrefixName(outname);

  /* Getting the Study name */

  if (argc < 2) messageAbort("Wrong number of arguments");
  ascii_study_define(argv[1]);

  /* Define the environment */

  ascii_filename("Environ", 0, 0, filename);
  ascii_environ_read(filename, verbose);

  /* Define the data */

  ascii_filename("Data", 0, 0, filename);
  dbin = Db::createFromNF(filename, verbose);
  if (dbin == nullptr) goto label_end;
  dbfmt.setFlags(true, false, true, true, true);
  dbin->display(&dbfmt);

  /* Define the Default Space according to the Dimension of the Input Db */

  defineDefaultSpace(ESpaceType::RN, dbin->getNDim());

  /* Define the output grid file */

  ascii_filename("Grid", 0, 0, filename);
  dbout = DbGrid::createFromNF(filename, verbose);
  //  if (dbout != nullptr) dbout->display(&dbfmt);

  /* Define the variogram */

  ascii_filename("Vario", 0, 0, filename);
  vario = Vario::createFromNF(filename, verbose);
  if (vario != nullptr)
  {
    vario->compute(dbin, ECalcVario::VARIOGRAM);
    vario->display();
    ascii_filename("Vario", 0, 1, filename);
    if (!vario->dumpToNF(filename, EFormatNF::DEFAULT, verbose))
      messageAbort("ascii_vario_write");
  }

  /* Define the model */

  ascii_filename("Model", 0, 0, filename);
  model = Model::createFromNF(filename, verbose);
  if (model == nullptr) goto label_end;
  if (vario != nullptr)
  {
    if (model_fitting_sills(vario, model, constraints)) goto label_end;
    ascii_filename("Model", 0, 1, filename);
    if (!model->dumpToNF(filename, EFormatNF::DEFAULT, verbose))
      messageAbort("ascii_model_write");
  }

  /* Define the neighborhood */

  ascii_filename("Neigh", 0, 0, filename);
  neigh = NeighUnique::createFromNF(filename, verbose);
  if (neigh == nullptr)
    neigh = NeighImage::createFromNF(filename, verbose);
  if (neigh == nullptr)
    neigh = NeighBench::createFromNF(filename, verbose);
  if (neigh == nullptr)
    neigh = NeighMoving::createFromNF(filename, verbose);

  /* Look for simulations */

  ascii_filename("Simu", 0, 0, filename);
  ascii_simu_read(filename, verbose, &nbsimu, &nbtuba, &seed);

  /* Conditional expectation */
  if (dbin->getNInterval() > 0)
  {
    if (verbose) message("Performing Gibbs Sampler\n");
    dbin->clearLocators(ELoc::Z);
    if (gibbs_sampler(dbin, model,
                      1, seed, nboot, niter, false, false, true, false, false, 0,
                      5., true, true, true))
      messageAbort("gibbs_sampler");
    /* Set the current variable to the conditional expectation */
    dbin->setLocatorByUID(dbin->getNColumn() - 1, ELoc::Z, 0);
  }

  /* Perform the estimation */

  if (dbin != nullptr && model != nullptr && neigh != nullptr)
  {
    if (nbsimu > 0)
    {

      /* Simulation case */

      if (verbose) message("Performing Simulations");
      if (simtub(dbin, dbout, model, neigh, nbsimu, seed, nbtuba, 0))
        messageAbort("Simulations");
      dbfmt.setFlags(true, false, true, true, true);
      dbout->display(&dbfmt);
      dbout->dumpToNF("Simu.out", EFormatNF::DEFAULT, verbose);
    }
    else
    {
      if (dbout == nullptr)
      {

        /* Cross-validation */
        if (verbose) message("Performing Cross-validation\n");
        if (xvalid(dbin, model, neigh, false, 1, 0, 0)) messageAbort("xvalid");
        dbfmt.setFlags(true, false, true, true, true);
        dbin->display(&dbfmt);
      }
      else
      {

        /* Estimation case */
        if (verbose) message("Performing Kriging\n");

        if (dbin->getNLoc(ELoc::G) > 0)
        {
          double ballradius = 0.01;
          if (krigingGradient(dbin, dbout, model, neigh,
                              true, true, ballradius, true)) messageAbort("kriging");
        }
        else
        {
          if (kriging(dbin, dbout, model, neigh,
                      true, true, false)) messageAbort("kriging");
        }
        dbfmt.setFlags(true, false, true, true, true);
        dbout->display(&dbfmt);
        dbout->dumpToNF("Krige.out", EFormatNF::DEFAULT, verbose);
      }
    }
  }

  /* Core deallocation */

label_end:
  delete dbin;
  delete dbout;
  delete model;
  delete neigh;
  delete vario;
  return (0);
}
