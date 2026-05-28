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
#include "Simulation/CalcSimuPGS.hpp"
#include "Basic/VectorHelper.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Db/Db.hpp"
#include "Enum/ENeigh.hpp"
#include "Enum/ERule.hpp"
#include "Gibbs/AGibbs.hpp"
#include "Gibbs/GibbsFactory.hpp"
#include "LithoRule/RuleProp.hpp"
#include "Model/Model.hpp"
#include "Neigh/ANeigh.hpp"
#include "Simulation/ACalcSimulation.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "geoslib_old_f.h"
#include <cmath>

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
  CalcSimuPGS::CalcSimuPGS(
    Id nbsimu,
    Id nbtuba,
    Id gibbsNBurn,
    Id gibbsNIter,
    double percent,
    Id seed,
    bool verbose)
    : ACalcSimulation(nbsimu, seed, verbose)
    , _numberPGS(1)
    , _numberGRF{}
    , _numberFacies{}
    , _ngrftot(0)
    , _nfactot(0)
    , _flagGaus(false)
    , _flagProp(false)
    , _flagCheck(false)
    , _flagShow(false)
    , _nbtuba(nbtuba)
    , _gibbsNBurn(gibbsNBurn)
    , _gibbsNIter(gibbsNIter)
    , _percent(percent)
    , _ruleprop(nullptr)
    , _models{{nullptr, nullptr}, {nullptr, nullptr}}
    , _flagStat(false)
    , _flagCond(false)
    , _iattZ(2, -1)
    , _iptrRP(-1)
    , _iptrRF(-1)
    , _iptrDF(-1)
    , _iptrDN(-1)
    , _iptrRN(-1)
    , _nechInitial(0)
    , _flagUsed(2, VectorBool(2, false))
    , _propCst()
  {
  }

  bool CalcSimuPGS::_check()
  {
    if (!ACalcSimulation::_check()) return false;

    if (!hasDbout()) return false;

    if (_numberPGS != 1 && _numberPGS != 2)
    {
      messerr("Argument 'numberPGS' should be 1 or 2");
      return false;
    }
    if (_ruleprop == nullptr)
    {
      messerr("RuleProp must be defined");
      return false;
    }

    // The simulation engine is Turning Bands
    if (_nbtuba <= 0)
    {
      messerr("You must define 'nbtuba' > 0");
      return false;
    }

    // Global assignments
    _flagCond = isConditional();
    _flagStat = _ruleprop->isFlagStat();
    _propCst = _ruleprop->getPropCst();

    _nfactot = 0;
    _numberPGS = _ruleprop->getNRule();
    _numberGRF = {};
    _numberFacies = {};
    for (Id ipgs = 0; ipgs < _numberPGS; ++ipgs)
    {
      if (_ruleprop->getRule(ipgs) == nullptr)
      {
        messerr("Rule #%d is not defined in the RuleProp", ipgs + 1);
        return false;
      }
      _numberGRF[ipgs] = _ruleprop->getRule(ipgs)->getNGRF();
      _numberFacies[ipgs] = _ruleprop->getRule(ipgs)->getNFacies();
      _ngrftot += _numberGRF[ipgs];
      _nfactot += _numberFacies[ipgs];
    }

    // Check on the Input Db
    if (_flagCond)
    {
      _nechInitial = getDbin()->getNSample();
      if (!getDbin()->isNVarComparedTo(_numberPGS)) return false;
      for (Id ipgs = 0; ipgs < _numberPGS; ++ipgs)
        _iattZ[ipgs] = getDbin()->getUIDByLocator(ELoc::Z, ipgs);
    }

    /* Output Db */
    if (_flagProp && _flagGaus)
    {
      messerr(
        "Calculating the facies proportions is incompatible with storing the "
        "Gaussian values");
      return false;
    }

    // Check the Models
    for (Id ipgs = 0; ipgs < _numberPGS; ipgs++)
    {
      for (Id igrf = 0; igrf < 2; igrf++)
      {
        _flagUsed[ipgs][igrf] = _ruleprop->getRule(ipgs)->isYUsed(igrf);
        if (!_flagUsed[ipgs][igrf]) continue;
        if (_models[ipgs][igrf] == nullptr)
        {
          messerr(
            "Variable #%d needs the underlying GRF #%d", ipgs + 1, igrf + 1);
          messerr("No corresponding Model is provided");
          return false;
        }
        if (_models[ipgs][igrf]->getNVar() != 1)
        {
          messerr(
            "The number of variables in Model #%d (%d) for Variable %d should "
            "be 1",
            igrf + 1, ipgs + 1, _models[ipgs][igrf]->getNVar());
          return false;
        }
        if (_models[ipgs][igrf]->stabilize(_percent, true)) return false;
        if (_models[ipgs][igrf]->standardize(true)) return false;
      }
    }

    // Check the Neighborhood
    if (getNeigh()->getType() != ENeigh::UNIQUE
        && getNeigh()->getType() != ENeigh::BENCH)
    {
      messerr("The only authorized Neighborhoods are UNIQUE or BENCH");
      return false;
    }

    // Check the Rules
    if (_numberPGS > 1)
      for (Id ipgs = 0; ipgs < _numberPGS; ipgs++)
      {
        // Check the Rules (only ERule::STD case is authorized)
        if (_ruleprop->getRule(ipgs)->getModeRule() != ERule::STD)
        {
          messerr(
            "In the Bi-PGS application, only Standard Rule is authorized");
        }
      }

    // Final checks
    for (Id ipgs = 0; ipgs < _numberPGS; ipgs++)
    {
      if (_ruleprop->getRule(ipgs)->particularities(
            getDbout(), _ruleprop->getDbprop(), _models[ipgs][0], 1, _flagStat))
        return false;

      if (_flagCond)
      {
        getDbin()->clearLocators(ELoc::Z);
        getDbin()->setLocatorByUID(_iattZ[ipgs], ELoc::Z, 0);
      }
      if (!_isEnvironmentValid(ipgs)) return false;
    }
    return true;
  }

  bool CalcSimuPGS::_isEnvironmentValid(Id ipgs)
  {
    Id nvar = 0;
    Id nfex = 0;
    size_t ndim = getDbout()->getNDim();

    /**************************************************************/
    /* Check if the Space dimension is compatible with the method */
    /**************************************************************/

    if (ndim > 3)
    {
      messerr("The Turning Band Method is not a relevant simulation model");
      messerr("for this Space Dimension (%d)", ndim);
      return false;
    }

    /*********************************/
    /* Compatibility between two Dbs */
    /*********************************/

    if (_flagCond && !getDbin()->hasSameDimension(getDbout())) return false;

    /**********************/
    /* Checking the model */
    /**********************/

    Model* modelLocal = _models[ipgs][0];
    if (modelLocal != nullptr)
    {
      nvar = modelLocal->getNVar();
      if (nvar <= 0)
      {
        messerr(
          "The number of variables must be positive = %d",
          modelLocal->getNVar());
        return false;
      }
      if (_flagCond && getDbin()->getNLoc(ELoc::Z) != nvar)
      {
        messerr(
          "The number of variables of the Data (%d)",
          getDbin()->getNLoc(ELoc::Z));
        messerr(
          "does not match the number of variables of the Model (%d)", nvar);
        return false;
      }
      if (modelLocal->getNCov() <= 0)
      {
        messerr("The number of covariance must be positive");
        return false;
      }

      if (modelLocal->getNDim() <= 0)
      {
        messerr(
          "The Space Dimension must be positive = %d", modelLocal->getNDim());
        return false;
      }
      if (modelLocal->getNDim() != ndim)
      {
        messerr("The Space Dimension of the Db structure (%d)", ndim);
        messerr(
          "Does not correspond to the Space Dimension of the model (%d)",
          modelLocal->getNDim());
        return false;
      }

      nfex = modelLocal->getNExtDrift();
      if (_flagCond && nfex != 0 && !getDbout()->isGrid()
          && getDbin()->getNLoc(ELoc::F) != nfex)
      {
        messerr(
          "The Model requires %d external drift(s)",
          modelLocal->getNExtDrift());
        messerr(
          "but the input Db refers to %d external drift variables",
          getDbin()->getNLoc(ELoc::F));
        return false;
      }
      if (nfex != 0 && getDbout()->getNLoc(ELoc::F) != nfex)
      {
        messerr(
          "The Model requires %d external drift(s)",
          modelLocal->getNExtDrift());
        messerr(
          "but the output Db refers to %d external drift variables",
          getDbout()->getNLoc(ELoc::F));
        return false;
      }
    }

    /*********************************/
    /* Calculate the field extension */
    /*********************************/

    VectorDouble db_mini(ndim, TEST);
    VectorDouble db_maxi(ndim, TEST);

    getDbout()->getExtensionInPlace(db_mini, db_maxi, true);

    if (_flagCond) getDbin()->getExtensionInPlace(db_mini, db_maxi, true);

    if (modelLocal != nullptr)
      modelLocal->setField(VH::extensionDiagonal(db_mini, db_maxi));

    /*****************************/
    /* Checking the Neighborhood */
    /*****************************/

    if (_flagCond && getNeigh() != nullptr)
    {
      if (ndim != getNeigh()->getNDim())
      {
        messerr(
          "The Space Dimension of the Neighborhood (%d)",
          static_cast<Id>(getNeigh()->getNDim()));
        messerr(
          "does not correspond to the Space Dimension of the first Db (%d)",
          ndim);
        return false;
      }
      if (getNeigh()->getFlagXvalid()
          && getNeigh()->getType() != ENeigh::MOVING)
      {
        messerr(
          "The Cross-Validation can only be processed with Moving "
          "neighborhood");
        return false;
      }
    }
    return true;
  }

  bool CalcSimuPGS::_preprocess()
  {
    if (!ACalcSimulation::_preprocess()) return false;

    Id iptr;

    /* Storage of the proportions */
    if (_flagProp)
    {
      if (db_locator_attribute_add(
            getDbout(), ELoc::P, _nfactot, 0, 0., &_iptrRP))
        return false;
    }

    /* Storage of the facies simulations in the output file */
    if (db_locator_attribute_add(
          getDbout(), ELoc::FACIES, _numberPGS * getNbSimu(), 0, 0., &_iptrRF))
      return false;

    /* Storage of the facies simulations in the input file */
    if (_flagCond)
    {
      if (db_locator_attribute_add(
            getDbin(), ELoc::FACIES, _numberPGS * getNbSimu(), 0, 0., &_iptrDF))
        return false;
    }

    /* Gaussian transform of the facies input data */
    if (_flagCond)
    {
      if (db_locator_attribute_add(
            getDbin(), ELoc::GAUSFAC, _ngrftot * getNbSimu(), 0, 0., &iptr))
        return false;
    }

    /* Non-conditional gaussian simulations at data points */
    if (_flagCond)
    {
      if (db_locator_attribute_add(
            getDbin(), ELoc::SIMU, _ngrftot * getNbSimu(), 0, 0., &_iptrDN))
        return false;
    }

    /* Non-conditional gaussian simulations at target points */
    if (db_locator_attribute_add(
          getDbout(), ELoc::SIMU, _ngrftot * getNbSimu(), 0, 0., &_iptrRN))
      return false;

    if (_flagCond)
    {
      /* Lower bound at input data points */
      if (db_locator_attribute_add(getDbin(), ELoc::L, _ngrftot, 0, 0., &iptr))
        return false;

      /* Upper bound at input data points */
      if (db_locator_attribute_add(getDbin(), ELoc::U, _ngrftot, 0, 0., &iptr))
        return false;
    }
    return true;
  }

  bool CalcSimuPGS::_run()
  {
    PropDef propdef;
    if (propdef.define(
          true, _flagStat, _numberGRF, _numberFacies, getDbin(),
          _ruleprop->getDbprop(), _propCst))
      return false;

    /************************/
    /* Main loop on the PGS */
    /************************/

    for (Id ipgs = 0; ipgs < _numberPGS; ipgs++)
    {
      if (_flagCond)
      {
        getDbin()->clearLocators(ELoc::Z);
        getDbin()->setLocatorByUID(_iattZ[ipgs], ELoc::Z, 0);
      }

      if (_numberPGS == 1)
        propdef.defineRuleMethod(EProcessOper::COPY);
      else if (ipgs == 0)
        propdef.defineRuleMethod(EProcessOper::MARGINAL);
      else
        propdef.defineRuleMethod(EProcessOper::CONDITIONAL);

      /****************************************/
      /* Convert facies into gaussian at data */
      /****************************************/

      if (_flagCond)
      {

        // Create the Gibbs sampler (multi_mono case)
        std::vector<Model*> modvec;
        for (Id igrf = 0; igrf < 2; igrf++)
        {
          if (!_flagUsed[ipgs][igrf]) continue;
          modvec.push_back(_models[ipgs][igrf]);
        }
        AGibbs* gibbs = GibbsFactory::createGibbs(
          getDbin(), modvec, _ruleprop->getRule(ipgs)->getRho(), false);
        gibbs->init(
          _numberPGS, _numberGRF[ipgs], _gibbsNBurn, _gibbsNIter, getSeed());

        /* Allocate the covariance matrix inverted */

        if (gibbs->covmatAlloc(getVerbose())) return false;

        // Core allocation

        VectorVectorDouble y = gibbs->allocY();

        /* Loop on the simulations */

        for (Id isimu = 0; isimu < getNbSimu(); isimu++)
        {
          /* Update the proportions */

          for (Id igrf = 0; igrf < _numberGRF[ipgs]; igrf++)
            if (_ruleprop->getRule(ipgs)->evaluateBounds(
                  propdef, getDbin(), getDbout(), isimu, igrf, ipgs,
                  getNbSimu()))
              return false;

          if (gibbs->run(y, ipgs, isimu, getVerbose())) return false;
        }

        /* Convert gaussian to facies on data point */

        if (_numberPGS > 1)
          for (Id isimu = 0; isimu < getNbSimu(); isimu++)
          {
            if (_ruleprop->getRule(ipgs)->gaus2facData(
                  propdef, getDbin(), getDbout(), _flagUsed[ipgs], ipgs, isimu,
                  getNbSimu()))
              return false;
          }

        delete gibbs;
      }

      /***************************************************/
      /* Perform the conditional simulation for each GRF */
      /***************************************************/

      /* Define the environment variables for printout */

      Id local_seed = getSeed();
      for (Id igrf = 0; igrf < 2; igrf++)
      {
        if (!_flagUsed[ipgs][igrf]) continue;
        Id icase = propdef.getRank(ipgs, igrf);
        CalcSimuTurningBands situba(getNbSimu(), _nbtuba, local_seed);
        situba.setFlagAllocationAlreadyDone(true);
        local_seed = 0;
        if (situba.simulate(
              getDbin(), getDbout(), _models[ipgs][igrf], getNeigh(), icase,
              false, true))
          return false;
      }

      /* Convert gaussian to facies at target point */

      if (!_flagGaus)
        for (Id isimu = 0; isimu < getNbSimu(); isimu++)
          propdef.transformCategorical(
            _ruleprop->getRule(ipgs), getDbout(), getVerbose(), _flagUsed[ipgs],
            ipgs, isimu, getNbSimu());

      /* Update facies proportions at target points */

      if (_flagProp)
      {
        for (Id isimu = 0; isimu < getNbSimu(); isimu++)
          propdef.updateCategorical(
            getDbout(), getVerbose(), ipgs, isimu, getNbSimu());
        propdef.scaleCategorical(getDbout(), getVerbose(), ipgs, getNbSimu());
      }

      /* Check/show facies at data against facies at the closest grid node */

      if (_flagCond && !_flagGaus && _flagShow) _checkFaciesDataToGrid(ipgs);
    }
    return true;
  }

  void CalcSimuPGS::_checkFaciesDataToGrid(Id ipgs)
  {
    Id iech, jech, isimu, facdat, facres, number;
    Id nfacies = _numberFacies[ipgs];
    Id nbsimu = getNbSimu();

    /* Initializations */

    if (!getDbout()->isGrid()) return;
    auto* dbgrid = dynamic_cast<DbGrid*>(getDbout());
    check_mandatory_attribute(
      "st_check_facies_data2grid", dbgrid, ELoc::FACIES);
    number = 0;
    if (_flagCheck)
      mestitle(
        1, "Checking facies of data against closest grid node (PGS=%d)",
        ipgs + 1);

    /* Core allocation */

    Id ndim = getDbin()->getNDim();
    VectorDouble coor(ndim);

    /* Loop on the data */

    for (iech = 0; iech < _nechInitial; iech++)
    {
      if (!getDbin()->isActive(iech)) continue;
      facdat = static_cast<Id>(getDbin()->getZVariable(iech, 0));
      if (facdat < 1 || facdat > nfacies) continue;
      jech = index_point_to_grid(getDbin(), iech, 0, dbgrid, coor.data());
      if (jech < 0) continue;

      for (isimu = 0; isimu < nbsimu; isimu++)
      {
        facres = static_cast<Id>(
          dbgrid->getSimvar(ELoc::FACIES, jech, isimu, 0, ipgs, nbsimu, 1));
        if (_flagShow)
        {
          if (facdat == facres)
            dbgrid->setSimvar(
              ELoc::FACIES, jech, isimu, 0, ipgs, nbsimu, 1, -facdat);
          else
            dbgrid->setSimvar(
              ELoc::FACIES, jech, isimu, 0, ipgs, nbsimu, 1, 0.);
        }

        if (facdat == facres) continue;
        number++;

        /* The data facies is different from the grid facies */

        if (_flagCheck)
        {
          message("Inconsistency for Simulation (%d) between :\n", isimu + 1);
          message("- Facies (%d) at Data (#%d)\n", facdat, iech + 1);
          message("- Facies (%d) at Grid (#%d)\n", facres, jech + 1);
        }
      }
    }

    if (_flagCheck && number <= 0) message("No problem found\n");
  }

  bool CalcSimuPGS::_postprocess()
  {
    /* Free the temporary variables */
    _cleanVariableDb(2);

    // _renameVariable(
    //   2, VectorString(), ELoc::Z, 1, _iattOut, String(), getNbSimu());

    auto namconv = getNamingConvention();
    if (getDbout() != nullptr)
    {
      if (!_keep(RESULT, TYPE_PROP) && _iptrRP >= 0)
        getDbout()->deleteColumnsByLocator(ELoc::P);
      else
        namconv.setNamesAndLocators(
          NULL, VectorString(), ELoc::Z, -1, getDbout(), _iptrRP, "Props",
          _nfactot, false);

      if (!_keep(RESULT, TYPE_GAUS) && _iptrRN >= 0)
        getDbout()->deleteColumnsByLocator(ELoc::SIMU);
      else
        namconv.setNamesAndLocators(
          NULL, VectorString(), ELoc::Z, -1, getDbout(), _iptrRN, "Gaus",
          _ngrftot * getNbSimu(), false);

      if (!_keep(RESULT, TYPE_FACIES) && _iptrRF >= 0)
        getDbout()->deleteColumnsByLocator(ELoc::FACIES);
      else
        namconv.setNamesAndLocators(
          NULL, VectorString(), ELoc::Z, -1, getDbout(), _iptrRF, String(),
          _numberPGS * getNbSimu());
    }

    if (isConditional())
    {
      if (!_keep(DATA, TYPE_GAUS) && _iptrDN >= 0)
        getDbin()->deleteColumnsByLocator(ELoc::SIMU);
      else
        namconv.setNamesAndLocators(
          NULL, VectorString(), ELoc::Z, -1, getDbin(), _iptrDN, "Gaus",
          _ngrftot * getNbSimu(), false);

      if (!_keep(DATA, TYPE_FACIES) && _iptrDF >= 0)
        getDbin()->deleteColumnsByLocator(ELoc::FACIES);
      else
        namconv.setNamesAndLocators(
          NULL, VectorString(), ELoc::Z, -1, getDbin(), _iptrDF, String(),
          _numberPGS * getNbSimu(), false);

      getDbin()->deleteColumnsByLocator(ELoc::GAUSFAC);
      getDbin()->deleteColumnsByLocator(ELoc::L);
      getDbin()->deleteColumnsByLocator(ELoc::U);
    }

    // Some samples may have been added to the input Db for PGS with Shift or Shadow options.
    // They must be removed from the final Db.
    _suppressAddedSamples();

    return true;
  }

  void CalcSimuPGS::_suppressAddedSamples()
  {
    Db* dbin = getDbin();
    if (!isConditional()) return;
    Id nfrom = dbin->getNSample();
    for (Id iech = nfrom - 1; iech >= _nechInitial; iech--)
      (void)dbin->deleteSample(iech);
  }

  bool CalcSimuPGS::_keep(Id file, Id type) const
  {
    if (file == DATA)
    {
      // Input Db
      return false;
    }

    // Output Db
    switch (type)
    {
      case 0: /* Gaussian */ return (_flagGaus && !_flagProp);
      case 1: /* Facies */ return (!_flagGaus && !_flagProp);
      case 2: /* Proportion */ return (_flagProp);
    }
    return false;
  }

  void CalcSimuPGS::_rollback()
  {
    _cleanVariableDb(1);
  }

  /****************************************************************************/
  /*!
   **  Perform conditional or non-conditional Pluri-gaussian simulations
   **
   ** \return  Error return code
   **
   ** \param[in]  dbin        Input Db structure (optional)
   ** \param[in]  dbout       Output Db structure
   ** \param[in]  ruleprop    RuleProp structure
   ** \param[in]  model1      First Model structure
   ** \param[in]  model2      Second Model structure (optional)
   ** \param[in]  neigh       ANeigh structure
   ** \param[in]  nbsimu      Number of simulations
   ** \param[in]  seed        Seed for random number generator
   ** \param[in]  flag_gaus   1 if results must be gaussian; otherwise facies
   ** \param[in]  flag_prop   1 for facies proportion
   ** \param[in]  flag_check  1 if the facies at data must be checked against
   **                         the closest simulated grid node
   ** \param[in]  flag_show   1 if the grid node which coincides with the data
   **                         should be represented with the data facies
   **                         (only if flag_cond && !flag_gaus)
   ** \param[in]  nbtuba      Number of turning bands
   ** \param[in]  gibbs_nburn Number of bootstrap iterations
   ** \param[in]  gibbs_niter Maximum number of iterations
   ** \param[in]  percent     Amount of nugget effect added to too much continous
   **                         model (expressed in percentage of the total variance)
   ** \param[in]  namconv     Naming convention
   **
   ** \remark  The argument 'dbin' is optional: it must be defined only for
   ** \remark  conditional simulations.
   ** \remark  When conditional, the unique variable in the input Db structure
   ** \remark  should correspond to the facies index (starting from 1)
   **
   *****************************************************************************/
  Id simpgs(
    Db* dbin,
    Db* dbout,
    RuleProp* ruleprop,
    Model* model1,
    Model* model2,
    ANeigh* neigh,
    Id nbsimu,
    Id seed,
    Id flag_gaus,
    Id flag_prop,
    Id flag_check,
    Id flag_show,
    Id nbtuba,
    Id gibbs_nburn,
    Id gibbs_niter,
    double percent,
    const NamingConvention& namconv)
  {
    CalcSimuPGS simu(nbsimu, nbtuba, gibbs_nburn, gibbs_niter, percent, seed);
    ANeigh* neighLocal = ANeigh::createDefaultNeighborhood(neigh, dbin, dbout);

    // Set the members of the Calculator
    simu.setDbin(dbin);
    simu.setDbout(dbout);
    simu.setRuleProp(ruleprop);
    simu.setModel(0, 0, model1);
    simu.setModel(0, 1, model2);
    simu.setFlagGaus(flag_gaus);
    simu.setFlagCheck(flag_check);
    simu.setFlagProp(flag_prop);
    simu.setFlagShow(flag_show);
    simu.setNeigh(neighLocal);
    simu.setNamingConvention(namconv);

    // Run the calculator
    Id error = (simu.run()) ? 0 : 1;
    if (neigh != neighLocal) delete neighLocal;
    return error;
  }

  /****************************************************************************/
  /*!
   **  Perform the conditional or non-conditional Bi Pluri-gaussian
   **  simulations
   **
   ** \return  Error return code
   **
   ** \param[in]  dbin        Input Db structure (optional)
   ** \param[in]  dbout       Output Db structure
   ** \param[in]  ruleprop    Ruleprop definition
   ** \param[in]  model11     First Model structure for First Lithotype Rule
   ** \param[in]  model12     Second Model structure for First Lithotype Rule
   ** \param[in]  model21     First Model structure for Second Lithotype Rule
   ** \param[in]  model22     Second Model structure for Second Lithotype Rule
   ** \param[in]  neigh       ANeigh structure
   ** \param[in]  nbsimu      Number of simulations
   ** \param[in]  seed        Seed for random number generator
   ** \param[in]  flag_gaus   1 gaussian results; otherwise facies
   ** \param[in]  flag_prop   1 for facies proportion
   ** \param[in]  flag_check  1 if the facies at data must be checked against
   **                         the closest simulated grid node
   ** \param[in]  flag_show   1 if the grid node which coincides with the data
   **                         should be represented with the data facies
   **                         (only if flag_cond && !flag_gaus)
   ** \param[in]  nbtuba      Number of turning bands
   ** \param[in]  gibbs_nburn Number of bootstrap iterations
   ** \param[in]  gibbs_niter Maximum number of iterations
   ** \param[in]  percent     Amount of nugget effect added to too continuous
   **                         model (expressed in percentage of the total variance)
   ** \param[in]  namconv     Naming convention
   **
   ** \remark  When conditional, the two first variables in the input Db
   ** \remark  should correspond to the two facies indices (starting from 1)
   ** \remark  The argument 'dbin' is optional: it must be defined only for
   ** \remark  conditional simulations
   ** \remark  The proportions (nfac1 * nfac2) must be ordered as follows:
   ** \remark  f1af2a, f1bf2a, f1cf2a, ..., f1bf2a, f1bf2b, ..., f1nf2m
   **
   *****************************************************************************/
  Id simbipgs(
    Db* dbin,
    Db* dbout,
    RuleProp* ruleprop,
    Model* model11,
    Model* model12,
    Model* model21,
    Model* model22,
    ANeigh* neigh,
    Id nbsimu,
    Id seed,
    Id flag_gaus,
    Id flag_prop,
    Id flag_check,
    Id flag_show,
    Id nbtuba,
    Id gibbs_nburn,
    Id gibbs_niter,
    double percent,
    const NamingConvention& namconv)
  {
    CalcSimuPGS simu(nbsimu, nbtuba, gibbs_nburn, gibbs_niter, percent, seed);
    ANeigh* neighLocal = ANeigh::createDefaultNeighborhood(neigh, dbin, dbout);

    // Set the members of the Calculator
    simu.setDbin(dbin);
    simu.setDbout(dbout);
    simu.setRuleProp(ruleprop);
    simu.setModel(0, 0, model11);
    simu.setModel(0, 1, model12);
    simu.setModel(1, 0, model21);
    simu.setModel(1, 1, model22);
    simu.setFlagGaus(flag_gaus);
    simu.setFlagCheck(flag_check);
    simu.setFlagProp(flag_prop);
    simu.setFlagShow(flag_show);
    simu.setNeigh(neighLocal);
    simu.setNamingConvention(namconv);

    // Run the calculator
    Id error = (simu.run()) ? 0 : 1;
    if (neigh != neighLocal) delete neighLocal;
    return error;
  }

} // namespace gstlrn
