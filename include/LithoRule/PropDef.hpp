/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/EProcessOper.hpp"

class Db;
class DbGrid;

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
