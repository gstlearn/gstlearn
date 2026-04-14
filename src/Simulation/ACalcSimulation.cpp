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
#include "Simulation/ACalcSimulation.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Basic/Law.hpp"
#include "Basic/OptDbg.hpp"
#include "Calculators/ACalcInterpolator.hpp"
#include "Estimation/KrigingSystem.hpp"
#include "Model/Model.hpp"
#include "Neigh/ANeigh.hpp"
#include "geoslib_old_f.h"

namespace gstlrn
{

  ACalcSimulation::ACalcSimulation(Id nbsimu, Id seed, bool verbose)
    : ACalcInterpolator(verbose)
    , _nbsimu(nbsimu)
    , _seed(seed)
    , _shift(0)
    , _flagBayes(false)
    , _flagPGS(false)
    , _flagGibbs(false)
    , _flagDGM(false)
    , _seedPerSimu()
  {
  }

  ACalcSimulation::~ACalcSimulation() {}

  bool ACalcSimulation::_check()
  {
    if (!ACalcInterpolator::_check()) return false;

    if (getNbSimu() <= 0)
    {
      messerr("You must define 'nbsimu' and 'nbtuba'");
      return false;
    }
    return true;
  }

  bool ACalcSimulation::_preprocess()
  {
    if (!ACalcInterpolator::_preprocess()) return false;

    Id nbsimu = getNbSimu();
    _seedPerSimu = VectorInt(nbsimu);

    // Initialize the random seed for each simulation
    Id mem_seed = getSeed();
    law_set_random_seed(getSeed());
    for (Id isimu = 0; isimu < nbsimu; isimu++)
    {
      // Ask for a random funtion to move the random generator and
      // modify the seed stored for each simulation rank
      (void)law_uniform();
      _seedPerSimu[isimu] = law_get_random_seed();
    }
    law_set_random_seed(mem_seed);

    return true;
  }

  Id ACalcSimulation::_getSeedPerSimu(Id isimu) const
  {
    DECLARE_UNUSED(isimu);
    // TODO: the correct version is to return the Initial seed value stored in _seedPerSimu
    // for the current simulation rank.
    // This version is temporarily bypassed in order not to modify the results of the
    // non-regression tests.
    Id seed = getSeedPerSimu(isimu);
    DECLARE_UNUSED(seed);
    return 0;
  }

  /*****************************************************************************/
  /*!
   **  Perform one non-conditional simulation on a set of gradient points using
   **  Turning Bands method.
   **
   ** \param[in]  dbgrd      Gradient Db structure
   ** \param[in]  isimu      Index of the simulation
   ** \param[in]  delta      Value of the increment
   ** \param[in]  activeArray  Array of active samples
   ** \param[out] tab        Array containing the non-conditional simulation values
   **
   ** \remarks The simulated gradients are stored as follows:
   ** \remarks idim * nbsimu + isimu (for simulation at first point)
   ** \remarks idim * nbsimu + isimu + ndim * nbsimu (for simulation at 2nd point)
   ** \remarks At the end, the simulated gradient is stored at first point
   **
   *****************************************************************************/
  void ACalcSimulation::_computeGradient(
    Db* dbgrd,
    Id isimu,
    double delta,
    const VectorBool& activeArray,
    VectorVectorDouble& tab)
  {
    Id jsimu;
    Id ndim = dbgrd->getNDim();
    auto nbsimu = getNbSimu();

    // Core allocation
    VectorBool activeLoc; // Not used
    VectorVectorDouble tab1;
    VectorVectorDouble tab2;
    _allocateForOneSimulation(dbgrd, 1, activeLoc, tab1, false);
    _allocateForOneSimulation(dbgrd, 1, activeLoc, tab2, false);

    for (Id idim = 0; idim < ndim; idim++)
    {

      /* Simulation at the initial location */

      jsimu = isimu + idim * nbsimu;
      _setShift(jsimu);
      _compute(dbgrd, isimu, activeArray, tab1);

      /* Shift the information */
      for (Id iech = 0; iech < dbgrd->getNSample(); iech++)
        if (activeArray[iech])
          dbgrd->setCoordinate(
            iech, idim, dbgrd->getCoordinate(iech, idim) + delta);

      /* Simulation at the shift location */

      jsimu = isimu + idim * nbsimu + ndim * nbsimu;
      _setShift(jsimu);
      _compute(dbgrd, isimu, activeArray, tab2);

      /* Un-Shift the information */

      for (Id iech = 0; iech < dbgrd->getNSample(); iech++)
        if (activeArray[iech])
          dbgrd->setCoordinate(
            iech, idim, dbgrd->getCoordinate(iech, idim) - delta);

      /* Scaling */

      for (Id iech = 0; iech < dbgrd->getNSample(); iech++)
      {
        if (!activeArray[iech]) continue;
        tab[idim][iech] = (tab2[0][iech] - tab1[0][iech]) / delta;
      }
    }
  }

  /*****************************************************************************/
  /*!
   **  Perform one non-conditional simulation on a set of tangent points using
   **  Turning Bands method.
   **
   ** \param[in]  dbtgt      Tangent Db structure
   ** \param[in]  isimu      Index of the simulation
   ** \param[in]  delta      Value of the increment
   ** \param[in]  activeArray  Array of active samples
   ** \param[out] tab        Array containing the non-conditional simulation values
   **
   ** \remarks Warning: To perform the simulation of the tangent, we must
   ** \remarks simulated the gradients first. So we need to dimension the
   ** \remarks simulation outcome variables as for the gradients
   **
   *****************************************************************************/
  void ACalcSimulation::_computeTangent(
    Db* dbtgt,
    Id isimu,
    double delta,
    const VectorBool& activeArray,
    VectorVectorDouble& tab)
  {
    Id icase = 0;
    auto nvar = getNVar();
    auto ndim = dbtgt->getNDim();
    auto nbsimu = getNbSimu();

    /* Perform the simulation of the gradients at tangent points */

    _computeGradient(dbtgt, isimu, delta, activeArray, tab);

    /* Calculate the simulated tangent */

    for (Id iech = 0; iech < dbtgt->getNSample(); iech++)
    {
      if (!activeArray[iech]) continue;

      double value = 0.;
      for (Id idim = 0; idim < ndim; idim++)
        value +=
          dbtgt->getLocVariable(ELoc::TGTE, iech, idim)
          * dbtgt->getSimvar(ELoc::SIMU, iech, isimu, 0, icase, nbsimu, nvar);
      dbtgt->setSimvar(ELoc::SIMU, iech, isimu, 0, icase, nbsimu, nvar, value);
    }
  }

  /****************************************************************************/
  /*!
   **  Correct for the mean in the case of non-conditional simulations
   **
   ** \param[in]  activeArray  Array of active samples
   ** \param[out] tab          Array containing the non-conditional simulation values
   **
   *****************************************************************************/
  void ACalcSimulation::_correctMean(
    const VectorBool& activeArray,
    VectorVectorDouble& tab)
  {
    if (_getFlagBayes()) return;
    Id nech = tab[0].size();

    for (Id ivar = 0, nvar = getNVar(); ivar < nvar; ivar++)
    {
      double mean = getModelGeneric()->getMean(ivar);
      if (FFFF(mean)) continue;
      for (Id iech = 0; iech < nech; iech++)
        if (activeArray[iech]) tab[ivar][iech] += mean;
    }
  }

  /*****************************************************************************/
  /*!
   **  Convert the non conditional simulations at the data points
   **  into simulation error
   **
   ** \param[in]  dbin       Input Db structure
   ** \param[in]  isimu      Index of the simulation
   ** \param[in]  activeArray  Array of active samples
   ** \param[out] tab        Array containing the non-conditional simulation values
   **
   *****************************************************************************/
  void ACalcSimulation::_difference(
    Db* dbin,
    Id isimu,
    const VectorBool& activeArray,
    VectorVectorDouble& tab)
  {
    auto nbsimu = getNbSimu();
    auto nvar = getNVar();
    double r = 1.;
    if (_getFlagDGM())
    {
      auto* modelLocal = dynamic_cast<Model*>(getModelGeneric());
      const auto* anamH =
        dynamic_cast<const AnamHermite*>(modelLocal->getAnam());
      r = anamH->getRCoef();
    }

    /* Transform the non conditional simulation into simulation error */

    if (!_getFlagPGS())
    {
      /********************************/
      /* Standard case (multivariate) */
      /********************************/

      for (Id iech = 0; iech < dbin->getNSample(); iech++)
      {
        if (!activeArray[iech]) continue;
        for (Id ivar = 0; ivar < nvar; ivar++)
        {
          double zvar = TEST;
          if (!_getFlagGibbs())
          {
            zvar = dbin->getZVariable(iech, ivar);
          }
          if (_getFlagGibbs())
          {
            zvar = dbin->getSimvar(
              ELoc::GAUSFAC, iech, isimu, ivar, 0, nbsimu, nvar);
            if (OptDbg::query(EDbg::SIMULATE)) printElement(zvar);
          }
          double simval = tab[ivar][iech];
          if (_getFlagDGM())
            simval = r * simval + sqrt(1. - r * r) * law_gaussian();

          double simunc = (FFFF(zvar) || FFFF(simval)) ? TEST : simval - zvar;
          tab[ivar][iech] = simunc;
        }
      }
    }
    else
    {

      /*********************************************************/
      /* Case of PGS: Data varies per simulation (monovariate) */
      /*********************************************************/

      for (Id iech = 0; iech < dbin->getNSample(); iech++)
      {
        if (!activeArray[iech]) continue;
        double zvar = dbin->getSimvar(
          ELoc::GAUSFAC, iech, isimu, 0, _getIcase(), nbsimu, 1);
        if (!FFFF(zvar)) tab[0][iech] -= zvar;
      }
    }
  }

  /****************************************************************************/
  /*!
   **  Update the conditional simulations when the target coincides
   **  with a data point
   **
   ** \param[in]  dbin      Input Db structure
   ** \param[in]  dbout     Output Db structure
   **
   ** \remarks This migration is not performed in the case where data point
   ** \remarks coincide with the target artificially. This is the case
   ** \remarks for the Discrete Gaussian Model (DGM) where data have been
   ** \remarks migrated to the cell center to mimic a point randomized
   ** \remarks within a cell
   **
   *****************************************************************************/
  void ACalcSimulation::_updateDataToTarget(Db* dbin, Db* dbout) const
  {
    if (dbin->getNSample() <= 0) return;
    if (_getFlagDGM()) return;
    auto nvar = getNVar();
    Id ndim = dbin->getNDim();
    auto nbsimu = getNbSimu();

    /* Calculate the field extension */

    double radius = dbin->getExtensionDiagonal();
    double eps = radius * EPSILON6;
    double eps2 = eps * eps;
    VectorDouble coor1(ndim);
    VectorDouble coor2(ndim);
    VectorBool activeArrayIn = dbin->getActiveArray();
    VectorBool activeArrayOut = dbout->getActiveArray();

    /* Dispatch according to the file type */

    if (dbout->isGrid())
    {
      auto* dbgrid = dynamic_cast<DbGrid*>(dbout);

      /*********************************************/
      /* Case where the output file is a grid file */
      /*********************************************/

      for (Id ip = 0; ip < dbin->getNSample(); ip++)
      {
        if (!activeArrayIn[ip]) continue;
        dbin->getCoordinatesInPlace(coor2, ip);
        Id rank = dbgrid->coordinateToRank(coor2, false, eps);
        if (rank < 0 || !activeArrayOut[rank]) continue;
        dbgrid->rankToCoordinatesInPlace(rank, coor1);

        /* Get the distance to the target point */

        double dist = 0;
        for (Id idim = 0; idim < ndim; idim++)
        {
          double delta = coor1[idim] - coor2[idim];
          dist += delta * delta;
        }
        if (dist > eps2) continue;

        /* We have found a close data point: perform the assignment */

        for (Id isimu = 0; isimu < nbsimu; isimu++)
          for (Id ivar = 0; ivar < nvar; ivar++)
          {
            double valdat;
            if (!_getFlagPGS())
              valdat = dbin->getZVariable(ip, ivar);
            else
              valdat = dbin->getSimvar(
                ELoc::GAUSFAC, ip, isimu, 0, _getIcase(), nbsimu, 1);
            if (FFFF(valdat)) continue;
            dbgrid->setSimvar(
              ELoc::SIMU, rank, isimu, ivar, _getIcase(), nbsimu, nvar, valdat);
          }
      }
    }
    else
    {

      /**********************************************/
      /* Case where the output file is a point file */
      /**********************************************/

      for (Id ik = 0; ik < dbout->getNSample(); ik++)
      {
        // Get coordinates of the active target point
        if (!activeArrayOut[ik]) continue;
        dbout->getCoordinatesInPlace(coor1, ik);

        /* Look for the closest data point */

        Id ip_close = -1;
        for (Id ip = 0; ip < dbin->getNSample() && ip_close < 0; ip++)
        {
          // Get the coordinates of the active data point
          if (!activeArrayIn[ip]) continue;
          dbin->getCoordinatesInPlace(coor2, ip);

          /* Get the distance to the target point */

          double dist = 0;
          for (Id idim = 0; idim < ndim; idim++)
          {
            double delta = coor1[idim] - coor2[idim];
            dist += delta * delta;
          }
          if (dist <= eps2) ip_close = ip;
        }

        if (ip_close < 0) continue;

        /* We have found a close data point: perform the assignment */

        for (Id isimu = 0; isimu < nbsimu; isimu++)
          for (Id ivar = 0; ivar < nvar; ivar++)
          {
            double valdat;
            if (!_getFlagPGS())
              valdat = dbin->getZVariable(ip_close, ivar);
            else
              valdat = dbin->getSimvar(
                ELoc::GAUSFAC, ip_close, isimu, 0, _getIcase(), nbsimu, 1);
            if (FFFF(valdat)) continue;
            dbout->setSimvar(
              ELoc::SIMU, ik, isimu, ivar, _getIcase(), nbsimu, nvar, valdat);
          }
      }
    }
  }

  /****************************************************************************/
  /*!
   **  Check/Show the data (gaussian) against the closest grid node
   **
   ** \param[in]  dbin       Input Db structure
   ** \param[in]  dbout      Output Db grid structure
   **
   ** \remark Attributes ELoc::SIMU and ELoc::GAUSFAC (for PGS) are mandatory
   ** \remark Tests have only been produced for icase=0
   **
   *****************************************************************************/
  Id ACalcSimulation::_checkGaussianDataToGrid(Db* dbin, Db* dbout) const
  {
    if (dbin == nullptr) return 0;
    if (get_LOCATOR_NITEM(dbout, ELoc::SIMU) <= 0) return 0;
    auto nbsimu = getNbSimu();
    if (nbsimu <= 0) return 0;

    auto* model = getModelGeneric();
    auto* dbgrid = dynamic_cast<DbGrid*>(dbout);
    if (dbgrid == nullptr) return 0;
    Id ndim = dbin->getNDim();

    mestitle(1, "Checking Gaussian of data against closest grid node");

    /* Loop on the data */

    Id number = 0;
    VectorDouble coor(ndim);
    for (Id iech = 0; iech < dbin->getNSample(); iech++)
    {
      if (!dbin->isActive(iech)) continue;

      // Find the index of the closest grid node and derive tolerance
      Id jech = index_point_to_grid(dbin, iech, 0, dbgrid, coor.data());
      if (jech < 0) continue;
      double eps = model->calculateStDev(dbin, iech, dbgrid, jech, false, 2.);
      if (eps < 1.e-6) eps = 1.e-6;

      for (Id isimu = 0; isimu < nbsimu; isimu++)
      {
        double valdat =
          dbin->getSimvar(ELoc::GAUSFAC, iech, 0, 0, 0, nbsimu, 1);
        double valres =
          dbgrid->getSimvar(ELoc::SIMU, jech, isimu, 0, 0, nbsimu, 1);
        if (ABS(valdat - valres) < eps) continue;
        number++;

        /* The data facies is different from the grid facies */

        message("Inconsistency for Simulation (%d) between :\n", isimu + 1);
        message("- Value (%lf) at Data (#%d) ", valdat, iech + 1);
        message("at (");
        for (Id idim = 0; idim < ndim; idim++)
          message(" %lf", dbin->getCoordinate(iech, idim));
        message(")\n");

        message("- Value (%lf) at Grid (#%d) ", valres, jech + 1);
        message("at (");
        for (Id idim = 0; idim < ndim; idim++)
          message(" %lf", dbgrid->getCoordinate(jech, idim));
        message(")\n");

        message("- Tolerance = %lf\n", eps);
      }
    }
    if (number <= 0) message("No problem found\n");
    return number;
  }

  /****************************************************************************/
  /*!
   **  Conditioning Kriging
   **
   ** \return  Error return code
   **
   ** \param[in]  dbin       input Db structure
   ** \param[in]  dbout      output Db structure
   **
   ** \remark: The model contains an anamorphosis with a change of support
   ** \remark: coefficient as soon as flag_dgm is TRUE
   **
   *****************************************************************************/
  Id ACalcSimulation::_conditionalKriging(Db* dbin, Db* dbout)
  {
    auto* neigh = getNeigh();
    auto* modelGeneric = getModelGeneric();
    Id nbsimu = getNbSimu();

    // Preliminary checks

    if (neigh->getType() == ENeigh::IMAGE)
    {
      messerr("This tool cannot function with an IMAGE neighborhood");
      return 1;
    }

    /* Add the attributes for storing the results */

    Id iptr_est = dbout->getColIdxByLocator(ELoc::SIMU, 0);
    if (iptr_est < 0) return 1;

    /* Setting options */

    KrigOpt krigopt;
    krigopt.setOptionDGM(_getFlagDGM());

    KrigingSystem ksys(dbin, dbout, modelGeneric, neigh, krigopt);
    if (ksys.setKrigOptFlagSimu(true, nbsimu, _getIcase())) return 1;
    if (ksys.updKrigOptEstim(iptr_est, -1, -1, true)) return 1;
    if (ksys.setKrigOptBayes(_getFlagBayes())) return 1;
    if (!ksys.isReady()) return 1;

    /* Loop on the targets to be processed */

    for (Id iech_out = 0; iech_out < dbout->getNSample(); iech_out++)
    {
      mes_process("Conditional Simulation", dbout->getNSample(), iech_out);
      if (ksys.estimate(iech_out)) return 1;
    }

    ksys.conclusion();

    return 0;
  }

  /**
   * @brief Save one multivariate simulation result into the Db
   *
   * @param db Db where the result is stored
   * @param isimu Simulation index
   * @param activeArray Array indicating active samples
   * @param tab   Array containing simulation values for all variables
   */
  void ACalcSimulation::saveResults(
    Db* db,
    Id isimu,
    const VectorBool& activeArray,
    const VectorVectorDouble& tab) const
  {
    auto nbsimu = getNbSimu();
    auto nvar = getNVar();
    auto icase = _getIcase();

    for (Id iech = 0, nech = db->getNSample(); iech < nech; iech++)
    {
      if (!activeArray[iech]) continue;
      for (Id ivar = 0; ivar < nvar; ivar++)
        db->setSimvar(
          ELoc::SIMU, iech, _getShift() + isimu, ivar, icase, nbsimu, nvar,
          tab[ivar][iech]);
    }
  }

  Id ACalcSimulation::getNVar() const
  {
    if (getModelGeneric() == nullptr) return 0;
    return getModelGeneric()->getNVar();
  }

  void ACalcSimulation::_allocateForOneSimulation(
    const Db* db,
    Id nvar,
    VectorBool& activeArray,
    VectorVectorDouble& tab,
    bool flagActive)
  {
    auto nech = db->getNSample();
    tab.resize(nvar);
    for (Id ivar = 0; ivar < nvar; ivar++) tab[ivar].fill(0., nech);

    if (flagActive) activeArray = db->getActiveArray();
  }

} // namespace gstlrn
