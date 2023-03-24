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

// WARNING: Make this include list as small as possible!
#include "geoslib_define.h"

#include "Enum/EKrigOpt.hpp"
#include "Model/Option_VarioFit.hpp"
#include "Mesh/AMesh.hpp"

class Koption
{
public:
  EKrigOpt calcul; /* Type of calculation (EKrigOpt) */
  int ndim; /* Space dimension */
  int ntot; /* Number of discretization points */
  int *ndisc; /* Array of discretization counts */
  double *disc1; /* Discretization coordinates */
  double *disc2; /* Discretization randomized coordinates */
  int flag_data_disc; /* Discretization flag */
  double *dsize;
};

class Model;
typedef struct
{
  int norder;
  int nmodel;
  int npar_init;
  Model *models[2];
  Option_VarioFit optvar;
  void *user_data;
  VectorInt parid;
  VectorDouble covtab;
} StrMod;

class Db;
typedef struct
{
  int case_facies; /* TRUE when Gibbs used for Facies */
  int case_stat; /* TRUE if proportions are constant */
  int case_prop_interp; /* TRUE when props are given in proportion file */
  int ngrf[2]; /* Number of GRF for the PGSs */
  int nfac[2]; /* Number of facies for the PGSs */
  int nfaccur; /* Number of facies for current PGS */
  int nfacprod; /* Product of the number of facies */
  int nfacmax; /* Maximum number of facies over all PGS */
  int mode; /* Type of process */
  VectorDouble propfix;
  VectorDouble propmem;
  VectorDouble propwrk;
  VectorDouble proploc;
  VectorDouble coor;
  const Db *dbprop; /* Pointer to the Proportion file */
} Props;

class Rule;
class PropDef;
typedef struct
{
  int ipgs;
  int flag_used[2];
  const Rule  *rule;
  PropDef *propdef;
} Modif_Categorical;

class DbGrid;
typedef struct
{
  int nalloc;
  int npair;
  int size_aux;
  int flag_dist;
  int *tab_iech;
  int *tab_jech;
  int *tab_ipas;
  int *tab_sort;
  char *tab_aux_iech;
  char *tab_aux_jech;
  double *tab_dist;
} Vario_Order;

typedef struct
{
  double coor[3];
  double intercept;
  double value;
  double rndval;
} SubPlan;

typedef struct
{
  int nplan;
  std::vector<SubPlan> plans;
} SubPlanes;

class QChol;
typedef struct
{
  QChol *QCtt;
  QChol *QCtd;
} QSimu;

class Cheb_Elem
{
public:
  int ncoeffs; /* Number of coefficients */
  int ncmax; /* Maximum number of polynomials */
  int ndisc; /* Number of discretizations */
  double power; /* Power of the transform */
  double a;
  double b;
  double v1;
  double v2;
  double tol; /* Tolerance */
  double *coeffs; /* Array of coefficients */
};

#ifndef SWIG
class cs_MGS;
typedef struct
{
  VectorDouble Lambda;
  cs *S;
  cs *Aproj;
  QChol *QC;
  QChol **QCov;
  double *Isill;
  double *Csill;
  QSimu *qsimu;
  cs_MGS *mgs;
  Cheb_Elem *s_cheb;
  AMesh *amesh;
} SPDE_Matelem;
#endif

typedef struct
{
  bool mesh_dbin;
  bool mesh_dbout;
  String triswitch;
} SPDE_SS_Option;

typedef struct
{
  std::vector<SPDE_SS_Option> options;
} SPDE_Option;

typedef struct
{
  double *res;
} CTable;

typedef struct
{
  int nconf;                // Number of covariance configurations
  int ndisc;                // Number of discretization steps
  int flag_cumul;           // 1 if storing integer from -infinity to value
                            // 0 if storing the value per discretized class
  double cmin;              // Minimum correlation value
  double cmax;              // Maximum correlation value
  double dc;                // Covariance class interval
  double dp;                // Probability quantum for discretization
  double *v;                // Array of thresholds (Dim: ndisc+1)
  CTable** CT;
} CTables;

struct Local_Relem;

struct Local_Split
{
  int oper;                   // Rank of operator
  int nrule;                  // Number of generated rules
  int nbyrule;                // Number of symbols in the Rules
  int *Srules;                // List of rules (Dim: [nitem][NRULE])
  int *Sfipos;                // Position of facies (Dim: [nprod][NCOLOR])
  Local_Relem *old_relem;     // Not allocated
  std::vector<Local_Relem *> relems;
};

struct Local_Relem
{
  VectorInt facies;           // List of facies
  int nrule;                  // Number of generated rules
  int nbyrule;                // Number of symbols in the Rules
  int nsplit;                 // Number of splits
  int *Rrules;                // List of rules (Dim: [nitem][NRULE])
  int *Rfipos;                // Position of facies (Dim: [nprod][NCOLOR])
  Local_Split *old_split;     // Not allocated
  std::vector<Local_Split *> splits;
};

struct Global_Res
{
  int ntot; // Total Number of Data
  int np;   // Number of active Data
  int ng;   // Number of grid nodes discretizing Domain
  double surface; // Surface of Domain
  double zest;    // Estimate
  double sse;     // Standard deviation of estimation
  double cvgeo;   // Coefficient of Variation
  double cvv;     // Variance of Domain
  VectorDouble weights; // Weights attached to data
};

typedef struct Local_Relem Relem;
typedef struct Local_Split Split;
