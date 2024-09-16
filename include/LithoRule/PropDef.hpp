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

class Db;
class DbGrid;
class Rule;
class RuleShadow;

class GSTLEARN_EXPORT PropDef
{
  // TODO To be transformed in private URGENT
public:
  int case_facies; /* TRUE when Gibbs used for Facies */
  int case_stat; /* TRUE if proportions are constant */
  int case_prop_interp; /* TRUE when props are given in proportion file */
  int ngrf[2]; /* Number of GRF for the PGSs */
  int nfac[2]; /* Number of facies for the PGSs */
  int nfaccur; /* Number of facies for current PGS */
  int nfacprod; /* Product of the number of facies */
  int nfacmax; /* Maximum number of facies over all PGS */
  EProcessOper mode; /* Type of process */
  VectorDouble propfix;
  VectorDouble propmem;
  VectorDouble propwrk;
  VectorDouble proploc;
  VectorDouble coor;
  const DbGrid *dbprop; /* Pointer to the Proportion file */
};

GSTLEARN_EXPORT int get_rank_from_propdef(PropDef* propdef, int ipgs, int igrf);
GSTLEARN_EXPORT int rule_thresh_define_shadow(PropDef* propdef,
                                              Db* dbin,
                                              const RuleShadow* rule,
                                              int facies,
                                              int iech,
                                              int isimu,
                                              int nbsimu,
                                              double* t1min,
                                              double* t1max,
                                              double* t2min,
                                              double* t2max,
                                              double* dsup,
                                              double* down);
GSTLEARN_EXPORT int rule_thresh_define(PropDef* propdef,
                                       Db* dbin,
                                       const Rule* rule,
                                       int facies,
                                       int iech,
                                       int isimu,
                                       int nbsimu,
                                       int flag_check,
                                       double* t1min,
                                       double* t1max,
                                       double* t2min,
                                       double* t2max);
GSTLEARN_EXPORT void proportion_rule_process(PropDef* propdef,
                                             const EProcessOper& mode);
GSTLEARN_EXPORT PropDef* proportion_manage(int mode,
                                           int flag_facies,
                                           int flag_stat,
                                           int ngrf1,
                                           int ngrf2,
                                           int nfac1,
                                           int nfac2,
                                           Db* db,
                                           const Db* dbprop,
                                           const VectorDouble& propcst,
                                           PropDef* proploc);
GSTLEARN_EXPORT void propdef_reset(PropDef* propdef);
GSTLEARN_EXPORT void proportion_print(PropDef* propdef);
