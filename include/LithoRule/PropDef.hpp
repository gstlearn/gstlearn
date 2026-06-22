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
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/EProcessOper.hpp"

#include <array>

namespace gstlrn
{
  class Db;
  class DbGrid;
  class Rule;
  class RuleShadow;

  class GSTLEARN_EXPORT PropDef
  {
  public:
    PropDef() = default;
    PropDef(const PropDef& m) = default;
    PropDef& operator=(const PropDef& m) = default;
    virtual ~PropDef() = default;

    void printInfo() const;
    Id getRank(Id ipgs, Id igrf) const;
    void defineRuleMethod(const EProcessOper& oper);
    void reset() const;

    Id proportionDefine(const Db* db, Id iech, Id isimu, Id nbsimu, Id* jech)
      const;
    Id define(
      bool flag_facies,
      bool flag_stat,
      const std::array<Id, 2>& ngrf,
      const std::array<Id, 2>& nfac,
      Db* db,
      const Db* dbprop,
      const VectorDouble& propcst);
    Id ruleThreshDefine(
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
      double* t2max) const;
    Id ruleThreshShadowDefine(
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
      double* sh_down) const;

    static void updateContinuous(Db* db, bool verbose, Id isimu, Id nbsimu);
    void updateCategorical(Db* db, bool verbose, Id ipgs, Id isimu, Id nbsimu)
      const;
    static void scaleContinuous(Db* db, bool verbose, Id nbsimu);
    void scaleCategorical(Db* db, bool verbose, Id ipgs, Id nbsimu) const;
    void transformCategorical(
      const Rule* rule,
      Db* db,
      bool verbose,
      const VectorBool& flag_used,
      Id ipgs,
      Id isimu,
      Id nbsimu) const;

    double getPropFix(Id ifac) const { return _propfix[ifac]; }

  private:
    Id _proportionTransform() const;
    Id _proportionLocate(Id ifac_ref) const;
    Id _getINDLOC(Id ifac1, Id ifac2) const;
    double _getPROPFIX(Id ifac, Id ifac2) const;
    double _getPROPWRK(Id ifac, Id ifac2) const;
    void _setPROPWRK(Id ifac, Id ifac2, double prop) const;
    bool _proportionChanged() const;
    Id _getFacies(Id ipgs, Id ifac) const;

    // TODO To be transformed in private URGENT

  public:
    bool _caseFacies; /* TRUE when Gibbs is used for Facies */
    bool _caseStat; /* TRUE if proportions are constant */
    bool _casePropInterp; /* TRUE when props are given in proportion file */
    std::array<Id, 2> _ngrf{}; /* Number of GRF for the PGSs */
    std::array<Id, 2> _nfac{}; /* Number of facies for the PGSs */
    Id _nfaccur; /* Number of facies for current PGS */
    Id _nfacprod; /* Product of the number of facies */
    EProcessOper _mode; /* Type of process */
    mutable VectorDouble _propfix;
    mutable VectorDouble _propmem;
    mutable VectorDouble _propwrk;
    mutable VectorDouble _proploc;
    mutable VectorDouble _coor;
    const DbGrid* _dbprop; /* Pointer not to be deleted */
  };

} // namespace gstlrn
