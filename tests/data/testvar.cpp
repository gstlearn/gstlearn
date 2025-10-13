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
#include "Enum/EFormatNF.hpp"
#include "geoslib_f.h"

#include "Basic/File.hpp"
#include "Basic/Law.hpp"
#include "Basic/OptDbg.hpp"
#include "Core/Ascii.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Model/Model.hpp"
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
  DbGrid* dbout;
  Vario* vario;
  Model* model;
  Option_AutoFit mauto;
  Option_VarioFit options;
  Constraints constraints;
  Id nbsimu, nbtuba, seed, flag_norm_sill, flag_goulard_used, repval;
  double gof;
  static bool verbose = false;

  /* Initializations */

  dbout             = nullptr;
  vario             = nullptr;
  model             = nullptr;
  flag_norm_sill    = 0;
  flag_goulard_used = 1;

  /* Standard output redirection to file */

  StdoutRedirect sr("Result.out", argc, argv, 2);

  /* Create the output name (for storage of dump files) */

  VectorString subparts = separateKeywords(argv[1]);
  int nargs             = static_cast<int>(subparts.size());
  String outname        = concatenateStrings("", subparts[nargs - 2], subparts[nargs - 1], "-");
  ASerializable::setPrefixName(outname);

  /* Setup constants */

  OptDbg::reset();

  /* Getting the Study name */

  if (argc < 2) messageAbort("Wrong number of arguments");
  ascii_study_define(argv[1]);

  /* Define the environment */

  ascii_filename("Environ", 0, 0, filename);
  ascii_environ_read(filename, verbose);

  /* Define the options */

  ascii_filename("Option", 0, 0, filename);
  if (ascii_option_defined(filename, "Norm_sill", &repval)) flag_norm_sill = repval;
  if (ascii_option_defined(filename, "Goulard_used", &repval)) flag_goulard_used = repval;

  /* Define the output grid file */

  ascii_filename("Grid", 0, 0, filename);
  dbout = DbGrid::createFromNF(filename, verbose);

  /* Look for simulations */

  ascii_filename("Simu", 0, 0, filename);
  ascii_simu_read(filename, verbose, &nbsimu, &nbtuba, &seed);

  /* Define the model */

  ascii_filename("Model", 0, 0, filename);
  model = Model::createFromNF(filename, verbose);
  if (model == nullptr) goto label_end;

  // Define and store the Space

  defineDefaultSpace(ESpaceType::RN, model->getNDim());

  /* Perform the non-conditional Simulation */

  if (dbout != nullptr)
  {
    if (simtub(nullptr, dbout, model, nullptr, nbsimu, seed, nbtuba, 0))
      messageAbort("Simulations");
    /* Set the current variable to the conditional expectation */
    dbout->setLocatorByUID(dbout->getNColumn() - 1, ELoc::Z, 0);
  }
  seed = law_get_random_seed();

  /* Define the variogram */

  ascii_filename("Vario", 0, 0, filename);
  vario = Vario::createFromNF(filename, verbose);
  if (vario == nullptr) goto label_end;
  if (dbout != nullptr)
  {
    vario->compute(dbout, ECalcVario::VARIOGRAM);
    ascii_filename("Vario", 0, 1, filename);
  }
  vario->display();
  if (!vario->dumpToNF("Vario.dat", EFormatNF::DEFAULT, verbose))
    messageAbort("ascii_vario_write");

  /* Fit the model */

  if (flag_norm_sill) constraints.setConstantSillValue(1.);
  options.setFlagGoulardUsed(flag_goulard_used);
  //  OptDbg::define(EDbg::CONVERGE);
  //  verbose = true;
  // Warning: initially (before using NLOPT), we discarded the use of Eigen library
  // in order to prevent diffs across platforms. This is now impossible.
  if (model_auto_fit(vario, model, verbose, mauto, constraints, options))
    messageAbort("model_auto_fit");
  model->display();
  ascii_filename("Model", 0, 1, filename);
  if (!model->dumpToNF("Model.out", EFormatNF::DEFAULT, verbose))
    messageAbort("ascii_model_write");

  // Produce the Goodness-of-fit score

  gof = model->gofToVario(vario, false);
  Model::gofDisplay(gof, false);

  /* Core deallocation */

label_end:
  delete model;
  delete dbout;
  delete vario;
  return (0);
}
