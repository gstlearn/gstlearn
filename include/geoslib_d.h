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
#include "Matrix/MatrixSquare.hpp"
#include "Mesh/AMesh.hpp"
#include "Model/Option_VarioFit.hpp"

namespace gstlrn
{
  class Koption
  {
  public:
    EKrigOpt calcul; /* Type of calculation (EKrigOpt) */
    Id ndim; /* Space dimension */
    Id ntot; /* Number of discretization points */
    VectorInt ndisc; /* Array of discretization counts */
    VectorDouble disc1; /* Discretization coordinates */
    VectorDouble disc2; /* Discretization randomized coordinates */
    Id flag_data_disc; /* Discretization flag */
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
    Id ncoeffs; /* Number of coefficients */
    Id ncmax; /* Maximum number of polynomials */
    Id ndisc; /* Number of discretizations */
    double power; /* Power of the transform */
    double a;
    double b;
    double v1;
    double v2;
    double tol; /* Tolerance */
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
    MatrixSquare Isill;
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

} // namespace gstlrn
