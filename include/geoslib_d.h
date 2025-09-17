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

// WARNING: Make this include list as small as possible!
#include "geoslib_define.h"

#include "Enum/EKrigOpt.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Mesh/AMesh.hpp"
#include "Model/Option_VarioFit.hpp"

namespace gstlrn
{
class Koption
{
public:
  EKrigOpt calcul;    /* Type of calculation (EKrigOpt) */
  Id ndim;            /* Space dimension */
  Id ntot;            /* Number of discretization points */
  VectorInt ndisc;    /* Array of discretization counts */
  VectorDouble disc1; /* Discretization coordinates */
  VectorDouble disc2; /* Discretization randomized coordinates */
  Id flag_data_disc;  /* Discretization flag */
  VectorDouble dsize;
};

class Model;
typedef struct
{
  Id norder;
  Id nmodel;
  Id npar_init;
  Model* models[2];
  Option_VarioFit optvar;
  void* user_data;
  VectorInt parid;
  VectorDouble covtab;
} StrMod;

class Db;
typedef struct
{
  Id case_facies;      /* TRUE when Gibbs used for Facies */
  Id case_stat;        /* TRUE if proportions are constant */
  Id case_prop_interp; /* TRUE when props are given in proportion file */
  Id ngrf[2];          /* Number of GRF for the PGSs */
  Id nfac[2];          /* Number of facies for the PGSs */
  Id nfaccur;          /* Number of facies for current PGS */
  Id nfacprod;         /* Product of the number of facies */
  Id nfacmax;          /* Maximum number of facies over all PGS */
  Id mode;             /* Type of process */
  VectorDouble propfix;
  VectorDouble propmem;
  VectorDouble propwrk;
  VectorDouble proploc;
  VectorDouble coor;
  const Db* dbprop; /* Pointer to the Proportion file */
} Props;

class Rule;
class PropDef;
typedef struct
{
  Id ipgs;
  Id flag_used[2];
  const Rule* rule;
  PropDef* propdef;
} Modif_Categorical;

typedef struct
{
  double coor[3];
  double intercept;
  double value;
  double rndval;
} SubPlan;

typedef struct
{
  Id nplan;
  std::vector<SubPlan> plans;
} SubPlanes;

struct QChol;
typedef struct
{
  QChol* QCtt;
  QChol* QCtd;
} QSimu;

class Cheb_Elem
{
public:
  Id ncoeffs;   /* Number of coefficients */
  Id ncmax;     /* Maximum number of polynomials */
  Id ndisc;     /* Number of discretizations */
  double power; /* Power of the transform */
  double a;
  double b;
  double v1;
  double v2;
  double tol;          /* Tolerance */
  VectorDouble coeffs; /* Array of coefficients */
};

#ifndef SWIG
typedef struct
{
  VectorDouble Lambda;
  MatrixSparse* S;
  MatrixSparse* Aproj;
  QChol* QC;
  std::vector<QChol*> QCov;
  VectorDouble Isill;
  VectorDouble Csill;
  QSimu* qsimu;
  Cheb_Elem* s_cheb;
  AMesh* amesh;
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
  Id nconf;               // Number of covariance configurations
  Id ndisc;               // Number of discretization steps
  Id flag_cumul;          // 1 if storing integer from -infinity to value
                          // 0 if storing the value per discretized class
  double cmin;            // Minimum correlation value
  double cmax;            // Maximum correlation value
  double dc;              // Covariance class interval
  double dp;              // Probability quantum for discretization
  VectorDouble v;         // Vector of thresholds (Dim: ndisc+1)
  VectorVectorDouble res; // Dimension: [nconf][size]
} CTables;

struct Local_Relem;

struct Local_Split
{
  Id oper;                // Rank of operator
  Id nrule;               // Number of generated rules
  Id nbyrule;             // Number of symbols in the Rules
  VectorInt Srules;       // List of rules (Dim: [nitem][NRULE])
  VectorInt Sfipos;       // Position of facies (Dim: [nprod][NCOLOR])
  Local_Relem* old_relem; // Not allocated
  std::vector<Local_Relem*> relems;
};

struct Local_Relem
{
  VectorInt facies;       // List of facies
  Id nrule;               // Number of generated rules
  Id nbyrule;             // Number of symbols in the Rules
  Id nsplit;              // Number of splits
  VectorInt Rrules;       // List of rules (Dim: [nitem][NRULE])
  VectorInt Rfipos;       // Position of facies (Dim: [nprod][NCOLOR])
  Local_Split* old_split; // Not allocated
  std::vector<Local_Split*> splits;
};

typedef struct Local_Relem Relem;
typedef struct Local_Split Split;
} // namespace gstlrn
