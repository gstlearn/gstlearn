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
#include "LithoRule/PropDef.hpp"
#include "Basic/Message.hpp"
#include "Basic/OptDbg.hpp"
#include "LithoRule/Rule.hpp"
#include "LithoRule/RuleShadow.hpp"
#include "geoslib_old_f.h"

namespace gstlrn
{
  Id PropDef::_getINDLOC(Id ifac1, Id ifac2) const
  {
    return ((ifac2)*_nfac[0] + (ifac1));
  }

  double PropDef::_getPROPFIX(Id ifac, Id ifac2) const
  {
    return _propfix[_getINDLOC(ifac, ifac2)];
  }

  double PropDef::_getPROPWRK(Id ifac, Id ifac2) const
  {
    return _propwrk[_getINDLOC(ifac, ifac2)];
  }

  void PropDef::_setPROPWRK(Id ifac, Id ifac2, double prop) const
  {
    _propwrk[_getINDLOC(ifac, ifac2)] = prop;
  }

  Id PropDef::getRank(Id ipgs, Id igrf) const
  {
    if (ipgs <= 0) return (igrf);
    return (_ngrf[0] + igrf);
  }

  /****************************************************************************/
  /*!
   **  Set the method to compute Proportions
   **
   ** \param[in]  oper     Type of operation (EProcessOper)
   **
   ****************************************************************************/
  void PropDef::defineRuleMethod(const EProcessOper& oper)
  {
    /* Assignments */

    _mode = oper;

    // Assign the current value for the number of facies
    switch (_mode.toEnum())
    {
      case EProcessOper::E_COPY:
      case EProcessOper::E_MARGINAL: _nfaccur = _nfac[0]; break;

      case EProcessOper::E_CONDITIONAL: _nfaccur = _nfac[1]; break;

      default: break;
    }

    /* In the stationary case, transform the proportions (from CST to WRK) */

    if (_caseStat) _proportionTransform();
  }

  /****************************************************************************/
  /*!
   **  Set the (non-stationary) proportions
   **
   ** \return  Error return code
   ** \return  - the target point does not lie within the proportion grid
   ** \return  - in conditional processing, the reference facies does not exist
   **
   ** \param[in]  db         Db input structure
   ** \param[in]  iech       Rank of the data in the input Db
   ** \param[in]  isimu      Rank of the simulation (EProcessOper::CONDITIONAL)
   ** \param[in]  nbsimu     Number of simulations
   **
   ** \param[out] jech       Rank of the auxiliary data in the input Db
   **
   ** \remark  At the end of this function, the local proportions are stored
   ** \remark  in the array proploc of the structure PropDef
   ** \remark  The argument 'isimu' is only used for
   ** \remark            mode == EProcessOper::CONDITIONAL (simbipgs)
   **
   *****************************************************************************/
  Id PropDef::proportionDefine(
    const Db* db,
    Id iech,
    Id isimu,
    Id nbsimu,
    Id* jech) const
  {
    Id ifac, ifac_ref;

    /* Non-stationary case : Load the proportions in propcst */

    (*jech) = 0;
    if (!_caseStat)
    {
      if (_casePropInterp)
      {

        /* Case where the proportions must be interpolated */

        (*jech) = index_point_to_grid(db, iech, 1, _dbprop, _coor.data());
        if ((*jech) < 0)
        {
          messerr(
            "At the data #%d, the proportion matrix is undefined", iech + 1);
          return (1);
        }

        /* Load the proportions (into CST) */

        double total = 0.;
        for (ifac = 0; ifac < _nfacprod; ifac++)
        {
          _propfix[ifac] = _dbprop->getLocVariable(ELoc::P, *jech, ifac);
          total += _propfix[ifac];
        }
        // Discard case where the proportions are not available
        if (total <= 0.) return 1;
      }
      else
      {

        /* The proportions are already available from the dbin */

        for (ifac = 0; ifac < _nfacprod; ifac++)
          _propfix[ifac] = db->getLocVariable(ELoc::P, iech, ifac);
      }

      /* Transform proportions (from CST to WRK) */

      _proportionTransform();
    }

    /* Locate the current proportions (from WRK to LOC) */

    ifac_ref = -1;
    if (_mode == EProcessOper::CONDITIONAL)
    {
      ifac_ref = static_cast<Id>(
        db->getSimvar(ELoc::FACIES, iech, isimu, 0, 0, nbsimu, 1));
      if (ifac_ref < 1 || ifac_ref > _nfac[0]) return (1);
    }
    _proportionLocate(ifac_ref);

    return (0);
  }

  /****************************************************************************/
  /*!
   **  Locate the current proportions
   **
   ** \param[in]  ifac_ref   Conditional (first variable) facies
   **                        (Only used for EProcessOper::CONDITIONAL)
   **
   ****************************************************************************/
  Id PropDef::_proportionLocate(Id ifac_ref) const
  {
    switch (_mode.toEnum())
    {
      case EProcessOper::E_COPY:
      case EProcessOper::E_MARGINAL: _proploc = _propwrk; break;

      case EProcessOper::E_CONDITIONAL:
        for (Id ifac = 0, nfac = _nfaccur; ifac < nfac; ifac++)
          _proploc[ifac] = _getPROPWRK(ifac_ref - 1, ifac);
        break;

      default: messerr("Unknown process operation"); break;
    }
    return (0);
  }

  /****************************************************************************/
  /*!
   **  Transform the proportions (from CST to WRK)
   **
   ** \return  -1 if the proportion is not defined; 0 otherwise
   **
   ****************************************************************************/
  Id PropDef::_proportionTransform() const

  {
    double total, pp;

    /* Dispatch */

    switch (_mode.toEnum())
    {
      case EProcessOper::E_COPY:
        for (Id ifac1 = 0, nfac1 = _nfac[0]; ifac1 < nfac1; ifac1++)
        {
          pp = _getPROPFIX(ifac1, 0);
          if (FFFF(pp)) return (-1);
          _setPROPWRK(ifac1, 0, pp);
        }
        break;

      case EProcessOper::E_MARGINAL:
        for (Id ifac1 = 0, nfac1 = _nfac[0]; ifac1 < nfac1; ifac1++)
        {
          _setPROPWRK(ifac1, 0, 0.);
          for (Id ifac2 = 0, nfac2 = _nfac[1]; ifac2 < nfac2; ifac2++)
          {
            pp = _getPROPFIX(ifac1, ifac2);
            if (FFFF(pp)) return (-1);
            _setPROPWRK(ifac1, 0, _getPROPWRK(ifac1, 0) + pp);
          }
        }
        break;

      case EProcessOper::E_CONDITIONAL:
        for (Id ifac1 = 0, nfac1 = _nfac[0]; ifac1 < nfac1; ifac1++)
          for (Id ifac2 = 0, nfac2 = _nfac[1]; ifac2 < nfac2; ifac2++)
          {
            pp = _getPROPFIX(ifac1, ifac2);
            if (FFFF(pp)) return (-1);
            _setPROPWRK(ifac1, ifac2, pp);
          }

        for (Id ifac1 = 0, nfac1 = _nfac[0]; ifac1 < nfac1; ifac1++)
        {
          total = 0.;
          for (Id ifac2 = 0, nfac2 = _nfac[1]; ifac2 < nfac2; ifac2++)
            total += _getPROPWRK(ifac1, ifac2);

          for (Id ifac2 = 0, nfac2 = _nfac[1]; ifac2 < nfac2; ifac2++)
            _setPROPWRK(
              ifac1, ifac2,
              (total <= 0.) ? 1. / _nfac[1]
                            : _getPROPWRK(ifac1, ifac2) / total);
        }
        break;

      default:
        messageAbort("This should never happen in _proportionTransform");
        break;
    }
    return (0);
  }

  void PropDef::reset() const
  {
    if (!_propmem.empty()) _propmem.fill(-1);
  }

  void PropDef::printInfo() const
  {
    mestitle(0, "Proportions");

    printMatrix(_propfix, _nfac[1], _nfac[0], "Initial :", 0, 1);

    printMatrix(_propwrk, _nfac[1], _nfac[0], "Working :", 0, 1);

    printMatrix(_proploc, _nfaccur, 1, "Current :", 0, 1);
  }

  /****************************************************************************/
  /*!
   **  Set the (non-stationary) proportions and define thresholds
   **
   ** \return  Error return code
   **
   ** \param[in]  db         Db input structure
   ** \param[in]  rule       Rule structure
   ** \param[in]  facies     Facies of interest (or ITEST) starting from 1
   ** \param[in]  iech       Rank of the data in the input Db
   ** \param[in]  isimu      Rank of the simulation (EProcessOper::CONDITIONAL)
   ** \param[in]  nbsimu     Number of simulations
   ** \param[in]  flag_check 1 if the consistency check with the actual
   **                        proportion of the current facies must be done
   **
   ** \param[out] t1min      Minimum threshold for Y1
   ** \param[out] t1max      Maximum threshold for Y1
   ** \param[out] t2min      Minimum threshold for Y2
   ** \param[out] t2max      Maximum threshold for Y2
   **
   *****************************************************************************/
  Id PropDef::ruleThreshDefine(
    Db* db,
    const Rule* rule,
    Id facies,
    Id iech,
    Id isimu,
    Id nbsimu,
    Id flag_check,
    double* t1min,
    double* t1max,
    double* t2min,
    double* t2max) const
  {
    /* Set the debugging information */

    OptDbg::setCurrentIndex(iech + 1);

    /* Processing an "unknown" facies */

    if (!isNA(facies) && (facies < 1 || facies > _nfaccur))
    {
      *t1min = *t2min = get_rule_extreme(-1);
      *t1max = *t2max = get_rule_extreme(+1);
      return (0);
    }

    /* Define the proportions */

    Id jech;
    if (proportionDefine(db, iech, isimu, nbsimu, &jech))
    {
      *t1min = *t2min = get_rule_extreme(-1);
      *t1max = *t2max = get_rule_extreme(+1);
      return (0);
    }

    /* Check if the proportions have been changed */

    bool modify = _proportionChanged();

    /* Check that the facies is compatible with the proportions */

    if (flag_check && !isNA(facies) && rule->getModeRule() == ERule::STD)
    {
      if (_proploc[facies - 1] <= 0.)
      {
        messerr(
          "The presence of facies (%d) at sample (%d) is not consistent with "
          "the zero proportion",
          facies, iech + 1);
        if (!_caseStat)
          messerr(
            "Check the proportions in the cell (%d) of the Proportion Db",
            jech + 1);
        return (1);
      }
    }

    /* Set the proportions and translate proportions into thresholds */

    if (modify)
    {
      if (rule->setProportions(_proploc)) return (1);

      /* In the case of SHIFT, update the thresholds */

      if (rule->getModeRule() == ERule::SHIFT && 0) rule->updateShift();
    }

    /* Convert the proportions into thresholds */

    Id facloc = (isNA(facies)) ? 1 : facies;
    auto bounds = rule->getThresh(facloc);
    *t1min = bounds[0];
    *t1max = bounds[1];
    *t2min = bounds[2];
    *t2max = bounds[3];

    return (0);
  }

  /****************************************************************************/
  /*!
   **  Check if the proportion has changed since the previous usage
   **  and store the current proportions for future comparison
   **
   ** \return  True if the proportions have changed; false otherwise
   **
   ****************************************************************************/
  bool PropDef::_proportionChanged() const
  {
    if (_proploc.isEqual(_propmem)) return false;

    if (OptDbg::query(EDbg::PROPS)) printInfo();

    _propmem = _proploc;

    return true;
  }

  /****************************************************************************/
  /*!
   **  Set the (non-stationary) proportions and define thresholds (for shadow only)
   **
   ** \return  Error return code
   **
   ** \param[in]  db         Db input structure
   ** \param[in]  rule       Rule structure
   ** \param[in]  facies     Facies of interest (or ITEST)
   ** \param[in]  iech       Rank of the data in the input Db
   ** \param[in]  isimu      Rank of the simulation (EProcessOper::CONDITIONAL)
   ** \param[in]  nbsimu     Number of simulations (EProcessOper::CONDITIONAL)
   **
   ** \param[out] t1min      Minimum threshold for Y1
   ** \param[out] t1max      Maximum threshold for Y1
   ** \param[out] t2min      Minimum threshold for Y2
   ** \param[out] t2max      Maximum threshold for Y2
   ** \param[out] sh_dsup    Local or global upwards shift (shadow)
   ** \param[out] sh_down    Local or global downwards shift (shadow)
   **
   *****************************************************************************/
  Id PropDef::ruleThreshShadowDefine(
    Db* db,
    const RuleShadow* rule,
    Id facies,
    Id iech,
    Id isimu,
    Id nbsimu,
    double* t1min,
    double* t1max,
    double* t2min,
    double* t2max,
    double* sh_dsup,
    double* sh_down) const
  {
    Id facloc, jech;

    /* Set the debugging information */

    OptDbg::setCurrentIndex(iech + 1);
    auto* dbgrid = dynamic_cast<DbGrid*>(db);

    /* Processing an "unknown" facies */

    if (!isNA(facies) && (facies < 1 || facies > _nfaccur))
    {
      *t1min = *t2min = get_rule_extreme(-1);
      *t1max = *t2max = get_rule_extreme(+1);
      return (0);
    }

    /* Define the proportions */

    if (proportionDefine(dbgrid, iech, isimu, nbsimu, &jech))
    {
      *t1min = *t2min = get_rule_extreme(-1);
      *t1max = *t2max = get_rule_extreme(+1);
      return (0);
    }

    /* Check if the proportions have been changed */

    bool modify = _proportionChanged();

    /* In case of Shadow, return the upwards and downwards values */

    *sh_dsup = (_caseStat) ? rule->getShDsup() : _proploc[1];
    *sh_down = (_caseStat) ? rule->getShDown() : _proploc[2];

    /* In the special cases, only the first proportion is significant */

    _proploc[1] = (1 - _proploc[0]) / 2;
    _proploc[2] = (1 - _proploc[0]) / 2;

    /* Set the proportions and translate proportions into thresholds */

    if (modify)
    {
      if (rule->setProportions(_proploc)) return (1);
    }

    /* Convert the proportions into thresholds */

    facloc = (isNA(facies)) ? 1 : facies;
    auto bounds = rule->getThresh(facloc);
    *t1min = bounds[0];
    *t1max = bounds[1];
    *t2min = bounds[2];
    *t2max = bounds[3];

    return (0);
  }

  /****************************************************************************/
  /*!
   **  Allocate or deallocate a proportion array
   **
   ** \return  Pointer on the returned PropDef structure
   **
   ** \param[in]  flag_facies 1 if Gibbs is used for facies
   ** \param[in]  flag_stat   1 if the proportions are stationary
   ** \param[in]  ngrf        Number of GRFs for the first PGS
   ** \param[in]  nfac        Number of facies for the first PGS
   ** \param[in]  db          Db structure containing the data
   ** \param[in]  dbprop      Db structure containing the proportions
   **                         (only used in the non-stationary case)
   ** \param[in]  propcst     Constant set of proportions (used if flag_stat)
   **
   ****************************************************************************/
  Id PropDef::define(
    bool flag_facies,
    bool flag_stat,
    const std::array<Id, 2>& ngrf,
    const std::array<Id, 2>& nfac,
    Db* db,
    const Db* dbprop,
    const VectorDouble& propcst)
  {
    const Db* db_loc;
    Id nfactot = nfac[0];
    if (nfac[1] > 0) nfactot *= nfac[1];

    /* Dispatch */

    _caseFacies = flag_facies;
    _caseStat = flag_stat;
    _casePropInterp = (dbprop != nullptr && dbprop->isGrid());
    _ngrf = ngrf;
    _nfac = nfac;
    _nfaccur = nfac[0];
    _nfacprod = nfactot;
    _mode = EProcessOper::UNDEFINED;
    if (_nfaccur <= 0)
    {
      messerr(" The number of facies may not be zero");
      return 1;
    }
    _propfix.resize(_nfacprod, 0.);
    _propwrk.resize(_nfacprod, 0.);
    _proploc.resize(_nfacprod, 0.);
    _propmem.resize(_nfacprod, 0.);

    if (flag_facies)
    {

      // Case of facies: Use of the proportions

      if (!flag_stat)
      {
        // Non-stationary case

        db_loc = (_casePropInterp) ? dbprop : db;
        if (db_loc == nullptr)
        {
          messerr("You have requested Non-stationary proportions");
          messerr("No file is provided containing Proportion variables");
          messerr("Please provide variables with 'proportion' locators");
          messerr("either in the input 'Db' or in 'dbprop'");
          return 1;
        }
        if (!db_loc->isGrid())
        {
          messerr("The 'Db' used for Proportions must be a Grid");
          return 1;
        }
        if (db_loc->getNLoc(ELoc::P) != _nfacprod)
        {
          messerr(
            "In the non-stationary case, the number of proportion variables "
            "(%d)",
            db_loc->getNLoc(ELoc::P));
          messerr(
            "must be equal to the number of facies (%d) in the Lithotype "
            "Rule",
            _nfacprod);
          return 1;
        }
        _dbprop = dynamic_cast<const DbGrid*>(dbprop);
        _coor.resize(db_loc->getNDim());
      }
      else
      {

        // Stationary case

        double pref = 1. / static_cast<double>(_nfacprod);
        for (Id ifac = 0; ifac < _nfacprod; ifac++)
        {
          _propfix[ifac] = (propcst.empty()) ? pref : propcst[ifac];
          _propwrk[ifac] = (propcst.empty()) ? pref : propcst[ifac];
          _proploc[ifac] = (propcst.empty()) ? pref : propcst[ifac];
          _propmem[ifac] = (propcst.empty()) ? pref : propcst[ifac];
        }
      }

      /* Set memory proportion so as to provoke the update at first usage */

      reset();
    }
    return 0;
  }

  Id PropDef::_getFacies(Id ipgs, Id ifac) const
  {
    if (ipgs <= 0) return (ifac);
    return (_nfac[0] + ifac);
  }

  void PropDef::updateCategorical(
    Db* db,
    bool verbose,
    Id ipgs,
    Id isimu,
    Id nbsimu) const
  {
    Id iptr_simu = Db::getSimRank(isimu, 0, ipgs, nbsimu, 1);

    /* Loop on the grid cells */

    for (Id iech = 0; iech < db->getNSample(); iech++)
    {
      if (!db->isActive(iech)) continue;
      Id facies =
        static_cast<Id>(db->getFromLocator(ELoc::FACIES, iech, iptr_simu)) - 1;
      Id rank = _getFacies(ipgs, facies);
      double prop = db->getLocVariable(ELoc::P, iech, rank) + 1.;
      db->setLocVariable(ELoc::P, iech, rank, prop);
    }

    /* Optional printout */

    if (verbose)
      message("Simulation Categorical Update (%d/%d)\n", isimu + 1, nbsimu);
  }

  void PropDef::updateContinuous(Db* db, bool verbose, Id isimu, Id nbsimu)
  {
    Id iptr_simu = Db::getSimRank(isimu, 0, 0, nbsimu, 1);

    /* Loop on the grid cells */

    for (Id iech = 0; iech < db->getNSample(); iech++)
    {
      if (!db->isActive(iech)) continue;
      double simval = db->getFromLocator(ELoc::SIMU, iech, iptr_simu);
      db->updLocVariable(ELoc::Z, iech, 0, EOperator::ADD, simval);
      db->updLocVariable(ELoc::Z, iech, 1, EOperator::ADD, simval * simval);
    }

    /* Optional printout */

    if (verbose)
      message("Simulation Continuous Update (%d/%d)\n", isimu + 1, nbsimu);
  }

  void PropDef::scaleContinuous(Db* db, bool verbose, Id nbsimu)
  {

    /* Loop on the grid cells */

    double mean, stdv;
    for (Id iech = 0; iech < db->getNSample(); iech++)
    {
      if (!db->isActive(iech)) continue;
      mean = db->getZVariable(iech, 0) / nbsimu;
      db->setLocVariable(ELoc::Z, iech, 0, mean);
      stdv = db->getZVariable(iech, 1) / nbsimu - mean * mean;
      stdv = (stdv > 0) ? sqrt(stdv) : 0.;
      db->setLocVariable(ELoc::Z, iech, 1, stdv);
    }

    /* Optional printout */

    if (verbose) message("Simulation Continuous Scaling (%d)\n", nbsimu);
  }

  void PropDef::scaleCategorical(Db* db, bool verbose, Id ipgs, Id nbsimu) const
  {
    Id nfacies = _nfac[ipgs];

    /* Loop on the grid cells */

    for (Id iech = 0; iech < db->getNSample(); iech++)
    {
      if (!db->isActive(iech)) continue;
      for (Id ifac = 0; ifac < nfacies; ifac++)
      {
        Id rank = _getFacies(ipgs, ifac);
        double prop =
          db->getLocVariable(ELoc::P, iech, rank) / static_cast<double>(nbsimu);
        db->setLocVariable(ELoc::P, iech, rank, prop);
      }
    }

    /* Optional printout */

    if (verbose) message("Simulation Categorical Scaling (%d)\n", nbsimu);
  }

  void PropDef::transformCategorical(
    const Rule* rule,
    Db* db,
    bool verbose,
    const VectorBool& flag_used,
    Id ipgs,
    Id isimu,
    Id nbsimu) const
  {
    rule->gaus2facResult(*this, db, flag_used, ipgs, isimu, nbsimu);

    /* Optional printout */

    if (verbose)
      message(
        "Simulation Categorical Transformation (%d/%d)\n", isimu + 1, nbsimu);
  }

} // namespace gstlrn
