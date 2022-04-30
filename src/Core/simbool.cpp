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
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "geoslib_f_private.h"
#include "geoslib_define.h"

#include "Basic/Law.hpp"
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Boolean/Tokens.hpp"
#include "Boolean/ObjectList.hpp"

/*****************************************************************************/
/*!
 **  Performs the boolean simulation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin          Db structure containing the data (optional)
 ** \param[in]  dbout         DbGrid structure containing the simulated grid
 ** \param[in]  tokens        Tokens structure
 ** \param[in]  flagStat      1 if the Intensity is constant
 ** \param[in]  thetaCst      Intensity constant value
 ** \param[in]  tmax          Maximum time
 ** \param[in]  seed          Seed for the random number generator
 ** \param[in]  flag_simu     Store the boolean simulation
 ** \param[in]  flag_rank     Store the object rank
 ** \param[in]  background    Value assigned to the background
 ** \param[in]  facies        Value of the facies assigned
 ** \param[in]  dilate        Array of dilation radius (optional)
 ** \param[in]  maxiter       Maximum number of iterations
 ** \param[in]  verbose       1 for a verbose output
 ** \param[in]  namconv       Naming convention
 **
 *****************************************************************************/
int simbool(Db* dbin,
            DbGrid* dbout,
            Tokens* tokens,
            bool flagStat,
            double thetaCst,
            double tmax,
            int seed,
            bool flag_simu,
            bool flag_rank,
            double background,
            double facies,
            const VectorDouble& dilate,
            int maxiter,
            bool verbose,
            const NamingConvention& namconv)
{
  int iptr_cover = -1;
  if (dbin != nullptr)
  {
    if (dbin->getVariableNumber() != 1)
    {
      messerr("Conditional Boolean simulation needs 1 variable");
      return 1;
    }
    iptr_cover = dbin->addColumnsByConstant(1, 0.,"Cover",ELoc::Z,1);
    if (iptr_cover < 0) return 1;
  }

  /* Add the attributes for storing the simulation */

  int iptr_simu = -1;
  if (flag_simu)
  {
    iptr_simu = dbout->addColumnsByConstant(1, background);
    if (iptr_simu < 0) return 1;
  }
  int iptr_rank = -1;
  if (flag_rank)
  {
    iptr_rank = dbout->addColumnsByConstant(1, TEST);
    if (iptr_rank < 0) return 1;
  }

  /* Define the global variables */

  law_set_random_seed(seed);

  /* Count the number of conditioning pores and grains */

  if (verbose) mestitle(0,"Boolean simulation");
  ObjectList objlist;
  int nbgrain = 0;
  int nbpore = 0;
  objlist.countConditioning(dbin, &nbgrain, &nbpore, verbose);

  /*******************************/
  /* Simulate the Initial grains */
  /*******************************/

  if (dbin != nullptr)
  {
    if (verbose)
    {
      message("- Conditioning option               = YES\n");
      mestitle(1, "Simulating the initial tokens");
      message("- Number of grains to be covered = %d\n", nbgrain);
      message("- Number of conditioning pores      = %d\n", nbpore);
    }

    if (objlist.generatePrimary(dbin, dbout, tokens, flagStat, thetaCst,
                                dilate, maxiter)) return 1;
    if (verbose)
      message("- Number of Initial Objects = %d\n",objlist.getNObjects(1));
  }
  else
  {
    if (verbose) message("- Conditioning option               = NO\n");
  }

  /*********************************/
  /* Simulate the Secondary grains */
  /*********************************/

  if (verbose)
  {
    mestitle(1, "Simulating the secondary tokens");
    message("- Maximum time available = %lf\n", tmax);
    if (flagStat)
      message("- Poisson Intensity                 = %g\n", thetaCst);
    else
      message("- Variable Poisson Intensity\n");
  }

  if (objlist.generateSecondary(dbin, dbout, tokens, flagStat, thetaCst,
                                tmax, dilate, maxiter)) return 1;

  /* Print the list of retained tokens */

  if (verbose) objlist.display();

  // Project the objects on the output grid

  objlist.projectToGrid(dbout, iptr_simu, iptr_rank, facies);

  // Final statistics

  if (verbose)
  {
    if (dbin != nullptr)
    {
      message("- Ending number of primary objects  = %d\n",objlist.getNObjects(1));
    }
    message("- Total number of objects           = %d\n",objlist.getNObjects());
  }

  namconv.setNamesAndLocators(dbin, ELoc::Z, 1, dbout, iptr_simu, "Facies", 1,
                              false);
  namconv.setNamesAndLocators(dbin, ELoc::Z, 1, dbout, iptr_rank, "Rank", 1,
                              false);
  return 0;
}

