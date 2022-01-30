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
#pragma once

#include "gstlearn_export.hpp"
#include "LithoRule/EProcessOper.hpp"

class Db;
class Dbgrid;

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
  const Dbgrid *dbprop; /* Pointer to the Proportion file */
};
