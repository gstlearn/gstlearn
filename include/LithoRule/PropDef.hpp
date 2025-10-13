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

namespace gstlrn
{
class Db;
class DbGrid;
class Rule;
class RuleShadow;

class GSTLEARN_EXPORT PropDef
{
  // TODO To be transformed in private URGENT
public:
  Id case_facies; /* TRUE when Gibbs used for Facies */
  Id case_stat; /* TRUE if proportions are constant */
  Id case_prop_interp; /* TRUE when props are given in proportion file */
  Id ngrf[2]; /* Number of GRF for the PGSs */
  Id nfac[2]; /* Number of facies for the PGSs */
  Id nfaccur; /* Number of facies for current PGS */
  Id nfacprod; /* Product of the number of facies */
  Id nfacmax; /* Maximum number of facies over all PGS */
  EProcessOper mode; /* Type of process */
  VectorDouble propfix;
  VectorDouble propmem;
  VectorDouble propwrk;
  VectorDouble proploc;
  VectorDouble coor;
  const DbGrid *dbprop; /* Pointer to the Proportion file */
};

GSTLEARN_EXPORT Id get_rank_from_propdef(PropDef* propdef, Id ipgs, Id igrf);
GSTLEARN_EXPORT Id rule_thresh_define_shadow(PropDef* propdef,
                                              Db* dbin,
                                              const RuleShadow* rule,
                                              Id facies,
                                              Id iech,
                                              Id isimu,
                                              Id nbsimu,
                                              double* t1min,
                                              double* t1max,
                                              double* t2min,
                                              double* t2max,
                                              double* dsup,
                                              double* down);
GSTLEARN_EXPORT Id rule_thresh_define(PropDef* propdef,
                                       Db* dbin,
                                       const Rule* rule,
                                       Id facies,
                                       Id iech,
                                       Id isimu,
                                       Id nbsimu,
                                       Id flag_check,
                                       double* t1min,
                                       double* t1max,
                                       double* t2min,
                                       double* t2max);
GSTLEARN_EXPORT void proportion_rule_process(PropDef* propdef,
                                             const EProcessOper& mode);
GSTLEARN_EXPORT PropDef* proportion_manage(Id mode,
                                           Id flag_facies,
                                           Id flag_stat,
                                           Id ngrf1,
                                           Id ngrf2,
                                           Id nfac1,
                                           Id nfac2,
                                           Db* db,
                                           const Db* dbprop,
                                           const VectorDouble& propcst,
                                           PropDef* proploc);
GSTLEARN_EXPORT void propdef_reset(PropDef* propdef);
GSTLEARN_EXPORT void proportion_print(PropDef* propdef);
}