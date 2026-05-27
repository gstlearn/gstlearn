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
#include "Basic/Law.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Core/Keypair.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "Gibbs/AGibbs.hpp"
#include "Gibbs/GibbsFactory.hpp"
#include "Gibbs/GibbsMMulti.hpp"
#include "Gibbs/GibbsUMultiMono.hpp"
#include "Gibbs/GibbsUPropMono.hpp"
#include "LithoRule/PropDef.hpp"
#include "LithoRule/Rule.hpp"
#include "LithoRule/RuleProp.hpp"
#include "LithoRule/RuleShadow.hpp"
#include "LithoRule/RuleShift.hpp"
#include "Mesh/MeshSpherical.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Simulation/CalcSimuBoolean.hpp"
#include "Simulation/CalcSimuEden.hpp"
#include "Simulation/CalcSimuFFT.hpp"
#include "Simulation/CalcSimuRefine.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Simulation/SimuSpherical.hpp"
#include "Simulation/SimuSphericalParam.hpp"
#include "Simulation/Simulations.hpp"
#include "geoslib_d.h"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

#include <cmath>
#include <cstring>

/*! \cond */
#define DATA 0
#define RESULT 1

// TODO : transform this to enum
#define TYPE_GAUS 0
#define TYPE_FACIES 1
#define TYPE_PROP 2

/*! \endcond */

namespace gstlrn
{

  static double GIBBS_RHO, GIBBS_SQR;

  /****************************************************************************/
  /*!
   **  Initialize the global values
   **
   *****************************************************************************/
  static void st_simulation_environment(void)
  {
    GIBBS_RHO = 0.;
    GIBBS_SQR = 0.;
  }

  /****************************************************************************/
  /*!
   **  Check for the presence of mandatory attributes
   **
   ** \param[in]  method  Name of the method
   ** \param[in]  db      Db structure
   ** \param[in]  locatorType  Mandatory attribute type
   **
   *****************************************************************************/
  void check_mandatory_attribute(
    const char* method,
    Db* db,
    const ELoc& locatorType)
  {
    if (get_LOCATOR_NITEM(db, locatorType) <= 0)
      messageAbort(
        "%s : Attributes %d are mandatory", method, locatorType.getValue());
  }

  /****************************************************************************/
  /*!
   **  Checks the environment for simulations by Turning Bands
   **
   ** \return  Error return code
   **
   ** \param[in]  dbin       input Db structure (optional if non conditional)
   ** \param[in]  dbout      output Db structure
   ** \param[in]  model      Model structure
   ** \param[in]  neigh      ANeigh structure (optional if non conditional)
   **
   *****************************************************************************/
  static Id st_check_simtub_environment(
    Db* dbin,
    Db* dbout,
    Model* model,
    ANeigh* neigh)
  {
    Id nvar = 0;
    Id nfex = 0;
    bool flag_cond = (dbin != nullptr);
    size_t ndim = dbout->getNDim();

    /**************************************************************/
    /* Check if the Space dimension is compatible with the method */
    /**************************************************************/

    if (ndim > 3)
    {
      messerr("The Turning Band Method is not a relevant simulation model");
      messerr("for this Space Dimension (%d)", ndim);
      return 1;
    }

    /*********************************/
    /* Compatibility between two Dbs */
    /*********************************/

    if (flag_cond && !dbin->hasSameDimension(dbout)) return 1;

    /**********************/
    /* Checking the model */
    /**********************/

    if (model != nullptr)
    {
      nvar = model->getNVar();
      if (nvar <= 0)
      {
        messerr(
          "The number of variables must be positive = %d", model->getNVar());
        return 1;
      }
      if (flag_cond && dbin->getNLoc(ELoc::Z) != nvar)
      {
        messerr(
          "The number of variables of the Data (%d)", dbin->getNLoc(ELoc::Z));
        messerr(
          "does not match the number of variables of the Model (%d)", nvar);
        return 1;
      }
      if (model->getNCov() <= 0)
      {
        messerr("The number of covariance must be positive");
        return 1;
      }

      if (model->getNDim() <= 0)
      {
        messerr("The Space Dimension must be positive = %d", model->getNDim());
        return 1;
      }
      if (model->getNDim() != ndim)
      {
        messerr("The Space Dimension of the Db structure (%d)", ndim);
        messerr(
          "Does not correspond to the Space Dimension of the model (%d)",
          model->getNDim());
        return 1;
      }

      nfex = model->getNExtDrift();
      if (flag_cond && nfex != 0 && !dbout->isGrid()
          && dbin->getNLoc(ELoc::F) != nfex)
      {
        messerr(
          "The Model requires %d external drift(s)", model->getNExtDrift());
        messerr(
          "but the input Db refers to %d external drift variables",
          dbin->getNLoc(ELoc::F));
        return 1;
      }
      if (nfex != 0 && dbout->getNLoc(ELoc::F) != nfex)
      {
        messerr(
          "The Model requires %d external drift(s)", model->getNExtDrift());
        messerr(
          "but the output Db refers to %d external drift variables",
          dbout->getNLoc(ELoc::F));
        return 1;
      }
    }

    /*********************************/
    /* Calculate the field extension */
    /*********************************/

    VectorDouble db_mini(ndim, TEST);
    VectorDouble db_maxi(ndim, TEST);

    dbout->getExtensionInPlace(db_mini, db_maxi, true);

    if (flag_cond) dbin->getExtensionInPlace(db_mini, db_maxi, true);

    if (model != nullptr)
      model->setField(VH::extensionDiagonal(db_mini, db_maxi));

    /*****************************/
    /* Checking the Neighborhood */
    /*****************************/

    if (flag_cond && neigh != nullptr)
    {
      if (ndim != neigh->getNDim())
      {
        messerr(
          "The Space Dimension of the Neighborhood (%d)",
          static_cast<Id>(neigh->getNDim()));
        messerr(
          "does not correspond to the Space Dimension of the first Db (%d)",
          ndim);
        return 1;
      }
      if (neigh->getFlagXvalid() && neigh->getType() != ENeigh::MOVING)
      {
        messerr(
          "The Cross-Validation can only be processed with Moving "
          "neighborhood");
        return 1;
      }
    }
    return 0;
  }

  /****************************************************************************/
  /*!
   **  Convert series of simulations to conditional expectation and variance
   **
   ** \return  Error return code
   **
   ** \param[in]  db        Db structure
   ** \param[in]  locatorType    Type of pointer containing the simulations
   ** \param[in]  nbsimu    Number of simulations
   ** \param[in]  nvar      Number of variables
   **
   ** \param[out] iptr_ce_arg   Pointer to the Conditional Expectation attributes
   ** \param[out] iptr_cstd_arg Pointer to the Conditional St. Dev. attributes
   **
   *****************************************************************************/
  Id db_simulations_to_ce(
    Db* db,
    const ELoc& locatorType,
    Id nbsimu,
    Id nvar,
    Id* iptr_ce_arg,
    Id* iptr_cstd_arg)
  {
    Id error, iptr_ce, iptr_cstd, iptr_nb, nech;
    double value, count, mean, var;

    // Initializations

    error = 1;
    iptr_ce = iptr_cstd = iptr_nb = -1;
    if (db == nullptr) goto label_end;
    nech = db->getNSample();
    if (nbsimu <= 0 || nvar <= 0 || nech <= 0) return (1);

    // Allocate the new attributes:

    iptr_ce = db->addColumnsByConstant(nvar, 0.);
    if (iptr_ce < 0) goto label_end;
    iptr_cstd = db->addColumnsByConstant(nvar, 0.);
    if (iptr_cstd < 0) goto label_end;
    iptr_nb = db->addColumnsByConstant(nvar, 0.);
    if (iptr_nb < 0) goto label_end;

    // Loop on the simulations

    for (Id isimu = 0; isimu < nbsimu; isimu++)
    {
      // Loop on the samples

      for (Id iech = 0; iech < nech; iech++)
      {
        if (!db->isActive(iech)) continue;

        // Loop on the variables

        for (Id ivar = 0; ivar < nvar; ivar++)
        {
          // Arguments 'simu' and 'nvar' are interchanged to keep correct order
          value =
            db->getSimvar(locatorType, iech, ivar, isimu, 0, nvar, nbsimu);
          if (FFFF(value)) continue;
          db->updArray(iech, iptr_ce + ivar, EOperator::ADD, value);
          db->updArray(iech, iptr_cstd + ivar, EOperator::ADD, value * value);
          db->updArray(iech, iptr_nb + ivar, EOperator::ADD, 1.);
        }
      }
    }

    // Scale the conditional expectation and variance

    for (Id iech = 0; iech < nech; iech++)
    {
      if (!db->isActive(iech)) continue;

      for (Id ivar = 0; ivar < nvar; ivar++)
      {
        count = db->getArray(iech, iptr_nb + ivar);
        if (count <= 0)
        {
          db->setArray(iech, iptr_ce + ivar, TEST);
          db->setArray(iech, iptr_cstd + ivar, TEST);
        }
        else
        {
          mean = db->getArray(iech, iptr_ce + ivar) / count;
          db->setArray(iech, iptr_ce + ivar, mean);
          var = db->getArray(iech, iptr_cstd + ivar) / count - mean * mean;
          var = (var > 0.) ? sqrt(var) : 0.;
          db->setArray(iech, iptr_cstd + ivar, var);
        }
      }
    }

    // Set the error return code

    error = 0;

  label_end:
    db->deleteColumnsByUIDRange(iptr_nb, nvar);
    if (error)
    {
      db->deleteColumnsByUIDRange(iptr_ce, nvar);
      db->deleteColumnsByUIDRange(iptr_cstd, nvar);
      *iptr_ce_arg = -1;
      *iptr_cstd_arg = -1;
    }
    else
    {
      *iptr_ce_arg = iptr_ce;
      *iptr_cstd_arg = iptr_cstd;
    }
    return (error);
  }

  /****************************************************************************/
  /*!
   **  Perform the Gibbs sampler
   **
   ** \return  Error return code
   **
   ** \param[in]  dbin        Db structure
   ** \param[in]  model       Model structure
   ** \param[in]  nbsimu      Number of simulations
   ** \param[in]  seed        Seed for random number generator
   ** \param[in]  gibbs_nburn Initial number of iterations for bootstrapping
   ** \param[in]  gibbs_niter Maximum number of iterations
   ** \param[in]  flag_moving      TRUE for Moving
   ** \param[in]  flag_norm        TRUE if the Model must be normalized
   ** \param[in]  flag_multi_mono  TRUE for the Multi_mono algorithm
   ** \param[in]  flag_propagation TRUE for the propagation algorithm
   ** \param[in]  flag_sym_neigh   Deprecated argument
   ** \param[in]  gibbs_optstats   0: No stats - 1: Print - 2: Save Neutral file
   ** \param[in]  percent     Amount of nugget effect added to too continuous
   **                         model (expressed in percentage of total variance)
   ** \param[in]  flag_ce     TRUE if the conditional expectation
   **                         should be returned instead of simulations
   ** \param[in]  flag_cstd   TRUE if the conditional standard deviation
   **                         should be returned instead of simulations
   ** \param[in]  verbose     Verbose flag
   ** \param[in]  namconv     Naming convention
   **
   ** \remark 'dbin' should contain the target variable(s) (Locator: Z)
   **         and the bounds (Locators: L and U).
   **
   *****************************************************************************/
  Id gibbs_sampler(
    Db* dbin,
    Model* model,
    Id nbsimu,
    Id seed,
    Id gibbs_nburn,
    Id gibbs_niter,
    bool flag_moving,
    bool flag_norm,
    bool flag_multi_mono,
    bool flag_propagation,
    bool flag_sym_neigh,
    Id gibbs_optstats,
    double percent,
    bool flag_ce,
    bool flag_cstd,
    bool verbose,
    const NamingConvention& namconv)
  {
    DECLARE_UNUSED(flag_sym_neigh);
    Id iptr;
    PropDef propdef;

    /* Initializations */

    Id error = 1;
    Id npgs = 1;
    Id nvar = 0;
    Id iptr_ce = -1;
    Id iptr_cstd = -1;
    AGibbs* gibbs = nullptr;
    std::vector<Model*> modvec;
    VectorVectorDouble y;

    /**********************/
    /* Preliminary checks */
    /**********************/

    /* Model */

    if (model == nullptr)
    {
      messerr("No Model is provided");
      return 1;
    }
    nvar = model->getNVar();
    if (!flag_propagation)
    {
      if (model->stabilize(percent, true)) return 1;
    }
    if (flag_norm)
    {
      if (model->standardize(true)) return 1;
    }

    /* Db */

    Id nbounds = dbin->getNBounds();
    if (flag_propagation)
    {
      if (nbounds > 0)
      {
        messerr("The propagation algorithm is incompatible with bounds");
        return 1;
      }
    }
    else
    {
      if (nbounds > 0)
      {
        if (dbin->getNLoc(ELoc::L) != nvar)
        {
          messerr("There must be as many Lower bound variables (%d)", nbounds);
          messerr("as there are variables defined in the Model (%d)", nvar);
          return 1;
        }
        if (dbin->getNLoc(ELoc::U) != nvar)
        {
          messerr("There must be as many Upper bound variables (%d)", nbounds);
          messerr("as there are variables defined in the Model (%d)", nvar);
          return 1;
        }

        // Check the consistency of the bounds
        // If Z-locator is defined, check the consistency with the bounds
        Id nvarDb = dbin->getNLoc(ELoc::Z);
        if (nvarDb > 0)
        {
          if (nvarDb != nvar)
          {
            messerr("Some Z-variables are defined");
            messerr(
              "Their count (%d) must match the number of variables defined in "
              "the Model (%d)",
              nvarDb, nvar);
            return 1;
          }

          // Convert the Z-values into bounds
          for (Id iech = 0; iech < dbin->getNSample(); iech++)
          {
            if (!dbin->isActive(iech)) continue;
            for (Id ivar = 0; ivar < nvar; ivar++)
            {
              double value = dbin->getLocVariable(ELoc::Z, iech, ivar);
              if (FFFF(value)) continue;

              // Set the bounds to the known exact value
              dbin->setLocVariable(ELoc::L, iech, ivar, value);
              dbin->setLocVariable(ELoc::U, iech, ivar, value);
            }
          }

          // Cancel the Z-locator for the rest of the Gibbs function
          dbin->clearLocators(ELoc::Z);
        }
      }
    }

    /*******************/
    /* Core allocation */
    /*******************/

    if (propdef.define(
          false, true, {1, 0}, {nvar, 0}, dbin, NULL, VectorDouble()))
      goto label_end;

    /**********************/
    /* Add the attributes */
    /**********************/

    if (db_locator_attribute_add(
          dbin, ELoc::GAUSFAC, nbsimu * nvar, 0, 0., &iptr))
      goto label_end;

    /*****************/
    /* Gibbs sampler */
    /*****************/

    if (!flag_multi_mono)
    {
      gibbs = GibbsFactory::createGibbs(dbin, model, flag_moving);
    }
    else
    {
      modvec.push_back(model);
      gibbs = GibbsFactory::createGibbs(dbin, modvec, 0., flag_propagation);
    }
    if (gibbs == nullptr) goto label_end;
    gibbs->setOptionStats(gibbs_optstats);
    gibbs->init(npgs, nvar, gibbs_nburn, gibbs_niter, seed);

    // Allocate the Gaussian vector

    y = gibbs->allocY();

    /* Allocate the covariance matrix inverted */

    if (gibbs->covmatAlloc(verbose)) goto label_end;

    // Invoke the Gibbs calculator

    for (Id isimu = 0; isimu < nbsimu; isimu++)
      if (gibbs->run(y, 0, isimu)) goto label_end;

    delete gibbs;

    /* Convert the simulations */

    if (flag_ce || flag_cstd)
    {
      if (db_simulations_to_ce(
            dbin, ELoc::GAUSFAC, nbsimu, nvar, &iptr_ce, &iptr_cstd))
        goto label_end;

      // We release the attributes dedicated to simulations on Dbout

      if (!flag_ce)
      {
        dbin->deleteColumnsByUIDRange(iptr_ce, nvar);
        iptr_ce = -1;
      }
      if (!flag_cstd)
      {
        dbin->deleteColumnsByUIDRange(iptr_cstd, nvar);
        iptr_cstd = -1;
      }
      dbin->deleteColumnsByLocator(ELoc::GAUSFAC);
    }

    /* Set the error return flag */

    error = 0;
    namconv.setNamesAndLocators(
      dbin, VectorString(), ELoc::UNDEFINED, nvar, dbin, iptr, String(),
      nbsimu);
    if (iptr_cstd >= 0)
      namconv.setNamesAndLocators(
        dbin, VectorString(), ELoc::UNDEFINED, nvar, dbin, iptr_cstd, "STD",
        nvar);
    if (iptr_ce >= 0)
      namconv.setNamesAndLocators(
        dbin, VectorString(), ELoc::UNDEFINED, nvar, dbin, iptr_ce, "CE", nvar);

  label_end:
    return (error);
  }

  /****************************************************************************/
  /*!
   **  Perform a set of valid conditional or non-conditional simulations
   **
   ** \return  Error return code
   **
   ** \param[in]  dbin         Input Db structure (optional)
   ** \param[in]  dbout        Output Db structure
   ** \param[in]  model        Model structure
   ** \param[in]  neigh        ANeigh structure (optional)
   ** \param[in]  seed         Seed for random number generator
   ** \param[in]  nbtuba       Number of turning bands
   ** \param[in]  nbsimu_min   Minimum number of simulations
   ** \param[in]  nbsimu_quant Additional quantum of simulations
   ** \param[in]  niter_max    Maximum number of iterations
   ** \param[in]  cols         Vector of column indices
   ** \param[in]  func_valid   Testing function
   **
   ** \remarks  This function calls a simulation outcome.
   ** \remarks  It provides a return code:
   ** \remarks  -1: the simulation outcome is not valid and the process must be
   ** \remarks      interrupted gently with no error (case of non-convergence)
   ** \remarks   0: the simulation outcome is not valid
   ** \remarks   1: the simulation outcome is valid and must be kept
   ** \remarks   2: the simulation outcome is valid and its returned version
   ** \remarks      must be kept (this is usefull if func_valid() has modified it).
   ** \remarks
   ** \remarks  The func_valid prototype has the following arguments:
   ** \li       ndim    Space dimension
   ** \li       nx      Array of number of grid meshes along all space direction
   ** \li       dx      Array of grid mesh along all space direction
   ** \li       x0      Array of grid origin along all space direction
   ** \li       nonval  Value corresponding to a missing grid value
   ** \li       percent Percentage of the simulations already validated
   ** \li       tab     Array containing the simulated outcome
   **
   ** \remarks  The following lines give an example of func_valid() which considers
   ** \remarks  a simulation outcome as valid if more than 50% of the valid samples
   ** \remarks  have a positive value.
   **
   ** \code
   **   Id func_valid(Id flag_grid,Id ndim,Id nech,
   **                  Id *nx,double *dx,double *x0,
   **                  double nonval, double percent, double *tab)
   **  {
   **    double ratio;
   **    Id i,npositive,nvalid;
   **
   **    for (i=0; i<nech; i++)
   **      {
   **         if (tab[i] == nonval) continue;
   **         nvalid++;
   **         if (tab[i] > 10.) npositive++;
   **      }
   **      ratio = (nvalid > 0) ? npositive / nvalid : 0.;
   **      return(ratio > 0.5);
   **   }
   ** \endcode
   **
   *****************************************************************************/
  Id simtub_constraints(
    Db* dbin,
    Db* dbout,
    Model* model,
    ANeigh* neigh,
    Id seed,
    Id nbtuba,
    Id nbsimu_min,
    Id nbsimu_quant,
    Id niter_max,
    VectorInt& cols,
    Id (*func_valid)(
      Id flag_grid,
      Id nDim,
      Id nech,
      Id* nx,
      double* dx,
      double* x0,
      double nonval,
      double percent,
      VectorDouble& tab))
  {
    Id iatt, retval, nbtest;
    Id error, nbsimu, nvalid, isimu, ndim, iter, nech, flag_grid, i;
    double percent;
    VectorDouble tab;
    VectorDouble dx;
    VectorDouble x0;
    VectorInt nx;

    /* Initializations */

    error = 1;
    law_set_random_seed(seed);
    cols.clear();

    /* Preliminary check */

    flag_grid = dbout->isGrid();
    ndim = dbout->getNDim();
    nech = dbout->getNSample();
    tab.resize(dbout->getNSample());
    if (flag_grid)
    {
      auto* dbgrid = dynamic_cast<DbGrid*>(dbout);
      nx.resize(ndim);
      dx.resize(ndim);
      x0.resize(ndim);

      for (i = 0; i < ndim; i++)
      {
        nx[i] = dbgrid->getNX(i);
        dx[i] = dbgrid->getDX(i);
        x0[i] = dbgrid->getX0(i);
      }
    }

    /* Implicit loop on the simulations */

    iatt = dbout->getNColumn();
    nvalid = iter = nbtest = 0;
    nbsimu = nbsimu_min + nbsimu_quant;
    while (nvalid < nbsimu_min && iter < niter_max)
    {

      /* Performing the simulations */

      iter++;
      nbtest += nbsimu;
      if (simtub(dbin, dbout, model, neigh, nbsimu, 0, nbtuba)) goto label_end;

      /* Check if the simulated outcomes are valid */

      for (isimu = 0; isimu < nbsimu; isimu++, iatt++)
      {

        /* Load the target simulation into the interface buffer */

        tab = dbout->getColumnByUID(iatt, true);

        /* Check if the simulation is valid */

        percent = 100. * nvalid / nbsimu_min;
        retval = func_valid(
          flag_grid, ndim, nech, nx.data(), dx.data(), x0.data(), TEST, percent,
          tab);
        if (retval == 0)
        {

          /* Delete the current simulation */

          dbout->deleteColumnByUID(iatt);

          /* Interrupt the loop (if requested) */

          if (retval < 0)
          {
            error = 0;
            goto label_end;
          }
        }
        else
        {

          /* The current simulation is accepted */

          if (retval > 1)
          {

            /* Update the vector (optional) */

            dbout->setColumnByUID(tab, iatt);
          }
          cols.push_back(iatt);
          nvalid++;
        }
      }

      /* Optional printout */

      if (OptDbg::query(EDbg::CONVERGE))
        message(
          "Iteration #%2d - Simulations %3d tested, %2d valid\n", iter, nbtest,
          nvalid);

      /* Define the number of simulations for the next batch */

      nbsimu = nbsimu_quant;
    }

    /* Set the error return code */

    error = 0;

  label_end:
    return (error);
  }

  /****************************************************************************/
  /*!
   **  Mask the grid nodes whose value is already too large
   **
   ** \return Number of cells left to be filled
   **
   ** \param[in]  dbout     Output Db structure
   ** \param[in]  seuil     Threshold
   ** \param[in]  scale     Scaling factor for the new simulation
   ** \param[in]  iptrv     Pointer to the max-stable outcome
   ** \param[in]  iptrs     Pointer to the current selection
   **
   *****************************************************************************/
  static Id
    st_maxstable_mask(Db* dbout, double seuil, double scale, Id iptrv, Id iptrs)
  {
    Id iech, number;
    double valsim;

    for (iech = number = 0; iech < dbout->getNSample(); iech++)
    {
      if (!dbout->isActive(iech)) continue;
      valsim = dbout->getArray(iech, iptrv);
      if (valsim > seuil / scale)
        dbout->setArray(iech, iptrs, 0.);
      else
        number++;
    }
    return (number);
  }

  /****************************************************************************/
  /*!
   **  Combine the simulations of the max-stable process
   **
   ** \param[in]  dbout     Output Db structure
   ** \param[in]  scale     Scaling factor for the new simulation
   ** \param[in]  iter0     Rank of the current iteration
   ** \param[in]  iptrg     Pointer to the newly simulated outcome
   ** \param[in]  iptrv     Pointer to the max-stable outcome
   ** \param[in]  iptrr     Pointer to the max-stable rank outcome
   **
   ** \param[in,out] last   Rank of Iteration where the last grid node is covered
   **
   *****************************************************************************/
  static void st_maxstable_combine(
    Db* dbout,
    double scale,
    Id iter0,
    Id iptrg,
    Id iptrv,
    Id iptrr,
    Id* last)
  {
    Id iech;
    double valsim, valold;

    for (iech = 0; iech < dbout->getNSample(); iech++)
    {
      if (!dbout->isActive(iech)) continue;
      valold = dbout->getArray(iech, iptrv);
      valsim = dbout->getArray(iech, iptrg) / scale;
      if (valsim > valold)
      {
        dbout->setArray(iech, iptrv, valsim);
        dbout->setArray(iech, iptrr, iter0);
        (*last) = iter0;
      }
    }
  }

  /****************************************************************************/
  /*!
   **  Perform the non-conditional simulation of the Max-Stable Model
   **
   ** \return  Error return code
   **
   ** \param[in]  dbout     Output Db structure
   ** \param[in]  model     Model structure
   ** \param[in]  ratio     Ratio modifying the range at each iteration
   ** \param[in]  seed      Seed for random number generator
   ** \param[in]  nbtuba    Number of turning bands
   ** \param[in]  flag_simu 1 if the simulation must be stored
   ** \param[in]  flag_rank 1 if the iteration rank must be stored
   ** \param[in]  verbose   Verbose flag
   **
   ** \remarks This function uses a threshold that can be defined using
   ** \remarks keypair mechanism with keyword "MaxStableThresh".
   **
   *****************************************************************************/
  Id simmaxstable(
    Db* dbout,
    Model* model,
    double ratio,
    Id seed,
    Id nbtuba,
    Id flag_simu,
    Id flag_rank,
    Id verbose)
  {
    double tpois, seuil;
    Id error, iptrg, iptrv, iptrr, iptrs, niter, nleft, icov, last;
    static double seuil_ref = 5.;

    /* Initializations */

    error = 1;
    iptrv = iptrg = iptrs = iptrr = -1;
    law_set_random_seed(seed);
    if (st_check_simtub_environment(NULL, dbout, model, NULL)) goto label_end;
    seuil = get_keypone("MaxStableThresh", seuil_ref);

    /* Preliminary checks */

    if (model->getNVar() != 1)
    {
      messerr("This feature is limited to the monovariate case");
      goto label_end;
    }
    if (!flag_simu && !flag_rank)
    {
      messerr("You must choose 'flag_simu' or  'flag_rank' or both");
      goto label_end;
    }

    /* Define the environment variables for printout */

    st_simulation_environment();

    /* Add the attributes for storing the results */

    iptrv = dbout->addColumnsByConstant(1, 0.);
    if (iptrv < 0) goto label_end;
    iptrr = dbout->addColumnsByConstant(1, 0.);
    if (iptrr < 0) goto label_end;
    if (db_locator_attribute_add(dbout, ELoc::SEL, 1, 0, 0., &iptrs))
      goto label_end;
    if (db_locator_attribute_add(dbout, ELoc::SIMU, 1, 0, 0., &iptrg))
      goto label_end;

    /* Implicit loop on the simulations */

    if (verbose)
    {
      message("Total number of cells = %d\n", dbout->getNSample());
      message("Maximum simulation value = %lf\n", seuil);
    }

    tpois = 0.;
    niter = last = 0;
    while (1)
    {
      niter++;
      tpois -= log(law_uniform(0., 1.));

      /* Mask the nodes that cannot be accessed anymore */

      nleft = st_maxstable_mask(dbout, seuil, tpois, iptrv, iptrs);
      if (nleft <= 0) break;

      /* Processing the Turning Bands algorithm */

      {
        CalcSimuTurningBands situba(1, nbtuba, seed);
        if (situba.simulate(nullptr, dbout, model, nullptr, 0)) goto label_end;
      }

      /* Combine the newly simulated outcome to the background */

      st_maxstable_combine(dbout, tpois, niter, iptrg, iptrv, iptrr, &last);

      if (verbose)
        message(
          "Iteration #%2d - Scale = %6.3lf - Nb. cells left = %d\n", niter,
          tpois, nleft);

      /* Update the model for next iteration */

      for (icov = 0; icov < model->getNCov(); icov++)
        model->setRangeIsotropic(icov, model->getRange(icov) * ratio);
    }

    if (verbose)
    {
      message("Number of iterations = %d\n", niter);
      message("Rank of the last covering iteration = %d\n", last);
    }

    /* Set the error return flag */

    error = 0;

  label_end:
    if (iptrs >= 0) dbout->deleteColumnByUID(iptrs);
    if (iptrg >= 0) dbout->deleteColumnByUID(iptrg);
    if (!flag_rank && iptrr >= 0) dbout->deleteColumnByUID(iptrr);
    if (!flag_simu && iptrv >= 0) dbout->deleteColumnByUID(iptrv);
    return (error);
  }

  /****************************************************************************/
  /*!
   **  Calculate the quantile for a given array
   **
   ** \return  Quantile value
   **
   ** \param[in]  dbout     Output Db structure
   ** \param[in]  proba     Probability
   **
   ** \param[out] sort      Sorting array
   **
   *****************************************************************************/
  static double st_quantile(Db* dbout, double proba, double* sort)
  {
    Id iech, nech, nval, rank;

    /* Initializations */

    nech = dbout->getNSample();

    /* Load the non-masked simulated values */

    for (iech = nval = 0; iech < nech; iech++)
    {
      if (!dbout->isActive(iech)) continue;
      sort[nval++] = dbout->getSimvar(ELoc::SIMU, iech, 0, 0, 0, 1, 1);
    }

    /* Sorting the array */

    ut_sort_double(0, nval, NULL, sort);

    /* Calculate the quantile */

    rank = static_cast<Id>(proba * static_cast<double>(nval));
    return (sort[rank]);
  }

  /****************************************************************************/
  /*!
   **  Perform the non-conditional simulation of the Orthogonal Residual Model
   **
   ** \return  Error return code
   **
   ** \param[in]  dbout     Output Db structure
   ** \param[in]  model     Model structure
   ** \param[in]  ncut      Number of cutoffs
   ** \param[in]  zcut      Array of cutoffs
   ** \param[in]  wcut      Array of weights
   ** \param[in]  seed      Seed for random number generator
   ** \param[in]  nbtuba    Number of turning bands
   ** \param[in]  verbose   Verbose flag
   **
   *****************************************************************************/
  Id simRI(
    Db* dbout,
    Model* model,
    Id ncut,
    double* zcut,
    double* wcut,
    Id seed,
    Id nbtuba,
    Id verbose)
  {
    double cumul, simval, proba, seuil;
    Id icut, error, iptrg, iptrs, nech, iech, count, total;
    VectorDouble sort;
    VectorDouble pton;
    VectorDouble pres;

    /* Initializations */

    error = 1;
    iptrg = iptrs = -1;
    nech = dbout->getNSample();
    law_set_random_seed(seed);
    if (st_check_simtub_environment(NULL, dbout, model, NULL)) goto label_end;

    /* Preliminary checks */

    if (model->getNVar() != 1)
    {
      messerr("This feature is limited to the monovariate case");
      goto label_end;
    }

    /* Define the environment variables for printout */

    st_simulation_environment();

    /* Add the attributes for storing the results */

    sort.resize(nech);
    pton.resize(ncut);
    pres.resize(ncut - 1);
    if (db_locator_attribute_add(dbout, ELoc::SEL, 1, 0, 0., &iptrs))
      goto label_end;
    if (db_locator_attribute_add(dbout, ELoc::SIMU, 1, 0, 0., &iptrg))
      goto label_end;

    /* Preliminary calculations */

    cumul = 0.;
    for (icut = 0; icut < ncut; icut++)
    {
      if (icut > 0 && zcut[icut] <= zcut[icut - 1])
      {
        messerr("The cutoff values must be ordered increasingly");
        goto label_end;
      }
      if (wcut[icut] < 0)
      {
        messerr("The weight of class (%d) cannot be negative", icut + 1);
        goto label_end;
      }
      cumul += wcut[icut];
    }
    if (cumul <= 0.)
    {
      messerr("The sum of weights cannot be negative or null");
      goto label_end;
    }
    for (icut = 0; icut < ncut; icut++) wcut[icut] /= cumul;
    pton[0] = 1.;
    for (icut = 1; icut < ncut; icut++)
      pton[icut] = pton[icut - 1] - wcut[icut];
    for (icut = 0; icut < ncut - 1; icut++)
      pres[icut] = pton[icut + 1] / pton[icut];

    /* Set the mask to the whole set of grid nodes */

    for (iech = 0; iech < nech; iech++)
      dbout->setLocVariable(ELoc::SEL, iech, 0, 1.);

    /* Loop on the cutoff classes */

    total = 0;
    for (icut = 0; icut < ncut; icut++)
    {

      /* Simulation in the non-masked part of the grid */

      {
        CalcSimuTurningBands situba(1, nbtuba, seed);
        if (situba.simulate(nullptr, dbout, model, nullptr, 0)) goto label_end;
      }

      /* Look for the quantile */

      proba = 1. - pres[icut];
      seuil = (icut < ncut - 1) ? st_quantile(dbout, proba, sort.data()) : TEST;

      /* Update the current selection */

      for (iech = count = 0; iech < nech; iech++)
      {
        if (!dbout->getSelection(iech)) continue;
        simval = dbout->getSimvar(ELoc::SIMU, iech, 0, 0, 0, 1, 1);
        if (!FFFF(seuil) && simval >= seuil) continue;
        dbout->setSimvar(
          ELoc::SIMU, iech, 0, 0, 0, 1, 1, static_cast<double>(icut + 1));
        dbout->setLocVariable(ELoc::SEL, iech, 0, 0.);
        count++;
      }
      total += count;
      if (verbose)
        message(
          "Level %3d - Proba=%lf - Affected=%7d - Total=%7d\n", icut + 1, proba,
          count, total);
    }

    /* Perform the final coding */

    for (iech = 0; iech < nech; iech++)
    {
      icut = static_cast<Id>(dbout->getSimvar(ELoc::SIMU, iech, 0, 0, 0, 1, 1));
      if (icut < 1 || icut > ncut)
        dbout->setSimvar(ELoc::SIMU, iech, 0, 0, 0, 1, 1, TEST);
      else
        dbout->setSimvar(ELoc::SIMU, iech, 0, 0, 0, 1, 1, zcut[icut - 1]);
    }

    /* Set the error return flag */

    error = 0;

  label_end:
    if (iptrs >= 0) dbout->deleteColumnByUID(iptrs);
    return (error);
  }

  /****************************************************************************/
  /*!
   **  Perform the conditional simulations under inequality constraints
   **
   ** \return  Error return code
   **
   ** \param[in]  dbin        Input Db structure (optional)
   ** \param[in]  dbout       Output Db structure
   ** \param[in]  model       Model structure
   ** \param[in]  seed        Seed for random number generator
   ** \param[in]  nbsimu      Number of simulations
   ** \param[in]  nbtuba      Number of turning bands
   ** \param[in]  gibbs_nburn Initial number of iterations for bootstrapping
   ** \param[in]  gibbs_niter Maximum number of iterations
   ** \param[in]  flag_ce     1 if the conditional expectation
   **                         should be returned instead of simulations
   ** \param[in]  flag_cstd   1 if the conditional standard deviation
   **                         should be returned instead of simulations
   ** \param[in]  verbose     Verbose flag
   **
   ** \remarks The Neighborhood does not have to be defined as this method
   ** \remarks only functions using a Unique Neighborhood. For consistency
   ** \remarks it is generated internally.
   **
   *****************************************************************************/
  Id simcond(
    Db* dbin,
    Db* dbout,
    Model* model,
    Id seed,
    Id nbsimu,
    Id nbtuba,
    Id gibbs_nburn,
    Id gibbs_niter,
    Id flag_ce,
    Id flag_cstd,
    Id verbose)
  {
    PropDef propdef;
    Id nvar, error, iptr, iptr_ce, iptr_cstd;

    /* Initializations */

    error = 1;
    bool flag_ext_created = false;
    bool flag_nostat_created = false;
    nvar = model->getNVar();
    iptr = -1;

    /* Preliminary checks */

    NeighUnique* neighU = NeighUnique::create(false);
    law_set_random_seed(seed);
    if (st_check_simtub_environment(dbin, dbout, model, NULL)) goto label_end;
    if (manageExternalInformation(1, ELoc::F, dbin, dbout, &flag_ext_created))
      goto label_end;
    if (manageExternalInformation(
          1, ELoc::NOSTAT, dbin, dbout, &flag_nostat_created))
      goto label_end;

    /* Limitations */

    if (nvar > 1)
    {
      messerr("This method is restricted to the monovariate case");
      goto label_end;
    }
    if (dbin->getNBounds() <= 0)
    {
      messerr("No bound is defined: use 'simtub' instead");
      goto label_end;
    }

    /* Core allocation */

    if (propdef.define(
          false, true, {1, 0}, {nvar, 0}, dbin, NULL, VectorDouble()))
      goto label_end;

    /* Define the environment variables for printout */

    st_simulation_environment();

    /* Add the attributes for storing the results */

    if (db_locator_attribute_add(dbin, ELoc::GAUSFAC, nbsimu, 0, 0., &iptr))
      goto label_end;
    if (db_locator_attribute_add(dbin, ELoc::SIMU, nvar * nbsimu, 0, 0., &iptr))
      goto label_end;
    if (db_locator_attribute_add(
          dbout, ELoc::SIMU, nvar * nbsimu, 0, 0., &iptr))
      goto label_end;

    /*****************/
    /* Gibbs sampler */
    /*****************/

    {
      AGibbs* gibbs = GibbsFactory::createGibbs(dbin, model, false);
      gibbs->init(1, 1, gibbs_nburn, gibbs_niter, seed);

      /* Allocate the covariance matrix inverted */

      if (gibbs->covmatAlloc(verbose)) goto label_end;

      /* Allocate the Gaussian Vector */

      VectorVectorDouble y = gibbs->allocY();

      /* Loop on the simulations */

      for (Id isimu = 0; isimu < nbsimu; isimu++)
      {
        if (gibbs->run(y, 0, isimu, verbose)) goto label_end;
      }

      delete gibbs;
    }

    /* Processing the Turning Bands algorithm */

    {
      CalcSimuTurningBands situba(nbsimu, nbtuba, seed);
      if (situba.simulate(dbin, dbout, model, neighU, 0, false, false, true))
        goto label_end;
    }

    /* Free the temporary variables not used anymore */

    dbin->deleteColumnsByLocator(ELoc::GAUSFAC);
    dbin->deleteColumnsByLocator(ELoc::SIMU);

    /* Convert the simulations */

    if (flag_ce || flag_cstd)
    {
      if (db_simulations_to_ce(
            dbout, ELoc::SIMU, nbsimu, nvar, &iptr_ce, &iptr_cstd))
        goto label_end;

      // We release the attributes dedicated to simulations on Dbout

      dbout->deleteColumnsByLocator(ELoc::SIMU);
      if (!flag_ce)
      {
        dbout->deleteColumnsByUIDRange(iptr_ce, nvar);
        iptr_ce = -1;
      }
      if (!flag_cstd)
      {
        dbout->deleteColumnsByUIDRange(iptr_cstd, nvar);
        iptr_cstd = -1;
      }
    }

    /* Set the error return flag */

    error = 0;

  label_end:
    delete neighU;
    (void)
      manageExternalInformation(-1, ELoc::F, dbin, dbout, &flag_ext_created);
    (void)manageExternalInformation(
      -1, ELoc::NOSTAT, dbin, dbout, &flag_nostat_created);
    return (error);
  }

  /*****************************************************************************/
  /*!
   **  Simulates the random function on the sphere
   **
   ** \param[in]  db          Data base containing the coordinates of target points
   **                         These coordinates must be expressed in long/lat
   ** \param[in]  model       Model (defined in Euclidean space) to be used
   ** \param[in]  sphepar     SimuSphericalParam structure
   ** \param[in]  seed        Seed for random number generation
   ** \param[in]  verbose     Verbose flag
   ** \param[in]  namconv     Naming convention
   **
   *****************************************************************************/
  Id simsph(
    DbGrid* db,
    Model* model,
    const SimuSphericalParam& sphepar,
    Id seed,
    bool verbose,
    const NamingConvention& namconv)
  {
    Id flag_sphere;

    /* Preliminary checks */

    flag_sphere = isDefaultSpaceSphere();
    if (!flag_sphere)
    {
      messerr(
        "The Spherical Simulation is restricted to Spherical coordinates");
      return 1;
    }
    if (db->getNDim() != 2)
    {
      messerr("The Simulation on Sphere is restricted to 2-D case");
      return 1;
    }
    for (Id icova = 0; icova < model->getNCov(); icova++)
    {
      if (model->getCovAniso(icova)->getFlagAniso())
      {
        messerr("Only Isotropic Models may be used for Spherical Simulations");
        return 1;
      }
    }

    /* Create the new variable in the Data base */

    Id iptr = db->addColumnsByConstant(1, 0., String(), ELoc::SIMU);

    SimuSpherical simsphe(1, seed);
    if (simsphe.simulate(db, model, sphepar, iptr, verbose)) return 1;

    namconv.setNamesAndLocators(
      db, VectorString(), ELoc::UNDEFINED, 1, db, iptr, "Simu");
    return 0;
  }

  /*****************************************************************************/
  /*!
   **  Simulates the random function on the sphere
   **
   ** \return The Vector simulated values
   **
   ** \param[in]  mesh        MeshSpherical object
   ** \param[in]  model       Model (defined in Euclidean space) to be used
   ** \param[in]  sphepar     SimuSphericalParam structure
   ** \param[in]  seed        Seed for random number generation
   ** \param[in]  verbose     Verbose flag
   **
   *****************************************************************************/
  VectorDouble simsph_mesh(
    MeshSpherical* mesh,
    Model* model,
    const SimuSphericalParam& sphepar,
    Id seed,
    Id verbose)
  {
    VectorDouble simu;

    bool flag_sphere = (model->getSpace()->getType() == ESpaceType::SN);
    if (!flag_sphere)
    {
      messerr(
        "The Spherical Simulation is restricted to Spherical coordinates");
      return simu;
    }
    for (Id icova = 0; icova < model->getNCov(); icova++)
    {
      if (model->getCovAniso(icova)->getFlagAniso())
      {
        messerr("Only Isotropic Models may be used for Spherical Simulations");
        return simu;
      }
    }

    SimuSpherical simsphe(1, seed);
    simu = simsphe.simulate_mesh(mesh, model, sphepar, verbose);

    return simu;
  }

  /****************************************************************************/
  /*!
   **  Check if the sample belongs to the time slice
   **
   ** \return  Rank of the time slice (or -1)
   **
   ** \param[in]  date     Date attached to a sample
   ** \param[in]  ntime    Number of time intervals
   ** \param[in]  time0    Origin of the first time interval
   ** \param[in]  dtime    Time interval
   **
   *****************************************************************************/
  static Id
    st_getTimeInterval(double date, Id ntime, double time0, double dtime)
  {
    for (Id itime = 0; itime < ntime; itime++)
    {
      double time_deb = time0 + dtime * itime;
      double time_fin = time0 + dtime * (itime + 1);
      if (date >= time_deb && date < time_fin) return (itime);
    }
    return (-1);
  }

  static Id
    st_getFACIES(const DbGrid* dbgrid, Id nfacies, Id indFacies, Id iech)
  {
    Id ifacies = static_cast<Id>(dbgrid->getArray(iech, indFacies));
    if (ifacies < 0 || ifacies > nfacies || isNA(ifacies)) ifacies = 0;
    return (ifacies);
  }

  static double st_getPORO(const DbGrid* dbgrid, Id indPoro, Id iech)
  {
    if (indPoro <= 0) return (1);
    double poro = dbgrid->getArray(iech, indPoro);
    if (FFFF(poro) || poro < 0.) poro = 0.;
    return (poro);
  }

  static double st_getDATE(const DbGrid* dbgrid, Id indDate, Id iech)
  {
    double date;

    if (indDate <= 0) return (0);
    date = dbgrid->getArray(iech, indDate);
    if (FFFF(date)) return (0);
    date = MAX(1., date);
    return (date);
  }

  static Id st_getFLUID(const DbGrid* dbgrid, Id nfluids, Id indFluid, Id iech)
  {
    Id ifluid = static_cast<Id>(dbgrid->getArray(iech, indFluid));
    if (ifluid < 0 || ifluid > nfluids || isNA(ifluid)) ifluid = 0;
    return (ifluid);
  }

  /*****************************************************************************/
  /*!
  **  Extract time charts from the fluid propagation block
  **
  ** \return  The returned matrix
  **
  ** \param[in]  dbgrid        Db grid structure
  ** \param[in]  name_facies   Name of variable containing Facies
  ** \param[in]  name_fluid    Name of variable containing Fluid
  ** \param[in]  name_poro     Name of variable containing Porosity (optional)
  ** \param[in]  name_date     Name of variable containing Date
  ** \param[in]  nfacies       number of facies (facies 0 excluded)
  ** \param[in]  nfluids       number of fluids
  ** \param[in]  facies0       Value of the target facies
  ** \param[in]  fluid0        Value of the target fluid
  ** \param[in]  ntime         Number of Time intervals
  ** \param[in]  time0         Starting time
  ** \param[in]  dtime         Time interval
  ** \param[in]  verbose       1 for a verbose option
  **
  *****************************************************************************/
  MatrixDense fluid_extract(
    DbGrid* dbgrid,
    const String& name_facies,
    const String& name_fluid,
    const String& name_poro,
    const String& name_date,
    Id nfacies,
    Id nfluids,
    Id facies0,
    Id fluid0,
    Id ntime,
    double time0,
    double dtime,
    bool verbose)
  {
    MatrixDense tab;
    if (!dbgrid->isGrid())
    {
      messerr("The Fluid Propagation is restricted to regular grid");
      return tab;
    }
    if (dbgrid->getNDim() > 3)
    {
      messerr("Fluid propagation is limited to 3-D space (maximum)");
      return tab;
    }
    if (ntime < 0 || time0 < 0 || dtime <= 0)
    {
      messerr("Error in Time Interval Definition");
      messerr("Origin=%lf - Step=%lf - Number=%d", time0, dtime, ntime);
      return tab;
    }

    /* Define global variables */

    Id ind_facies = dbgrid->getUID(name_facies);
    Id ind_fluid = dbgrid->getUID(name_fluid);
    Id ind_date = dbgrid->getUID(name_date);
    if (ind_facies < 0)
    {
      messerr("Variable 'Facies' must be provided");
      return 1;
    }
    if (ind_fluid < 0)
    {
      messerr("Variable 'Fluid' must be provided");
      return 1;
    }
    if (ind_date < 0)
    {
      messerr("Variable 'Date' must be provided");
      return 1;
    }
    Id ind_poro = dbgrid->getUID(name_poro);

    /* Initialize the array */

    tab = MatrixDense(ntime, 4);
    Id nxyz = dbgrid->getNSample();
    for (Id itime = 0; itime < ntime; itime++)
    {
      tab.setValue(itime, 0, time0 + dtime * itime);
      tab.setValue(itime, 1, time0 + dtime * (itime + 1));
      tab.setValue(itime, 2, 0.);
      tab.setValue(itime, 3, 0.);
    }

    /* Loop on the blocks */

    double totnum = 0.;
    double totvol = 0.;
    double locnum = 0.;
    double locvol = 0.;
    double datmax = 0;
    for (Id iech = 0; iech < nxyz; iech++)
    {

      if (st_getFACIES(dbgrid, nfacies, ind_facies, iech) != facies0) continue;
      if (st_getFLUID(dbgrid, nfluids, ind_fluid, iech) != fluid0) continue;
      double volume = st_getPORO(dbgrid, ind_poro, iech);
      double date = st_getDATE(dbgrid, ind_date, iech);
      if (date > datmax) datmax = date;

      totnum += 1;
      totvol += volume;
      auto itime = st_getTimeInterval(date, ntime, time0, dtime);
      if (itime < 0) continue;
      locnum += 1;
      locvol += volume;

      tab.setValue(itime, 2, tab.getValue(itime, 2) + 1);
      tab.setValue(itime, 3, tab.getValue(itime, 3) + volume);
    }

    /* Final printout */

    if (verbose)
    {
      mestitle(1, "Extraction for Fluid(%d) and Facies(%d)", facies0, fluid0);
      message(
        "Time slices: From %lf to %lf by step of %lf\n", time0,
        time0 + dtime * ntime, dtime);
      message("Total Number of Cells               = %d\n", nxyz);
      message("Maximum Date                        = %lf\n", datmax);
      message("Total Number of Invaded Cells       = %lf\n", totnum);
      message("Total Volume of Invaded Cells       = %lf\n", totvol);
      message("Total Number of Cells in Time Slice = %lf\n", locnum);
      message("Total Volume of Cells in Time Slice = %lf\n", locvol);
    }
    return tab;
  }
} // namespace gstlrn
