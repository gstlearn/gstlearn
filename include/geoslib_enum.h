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

/* Different ENUM lists */

#ifndef SWIG
// Internal enums (currently not exported via SWIG)
typedef enum
{
  TYPE_DB = 0,        //!< Data Base
  TYPE_VARIO = 1,     //!< Experimental variogram
  TYPE_MODEL = 2,     //!< Variogram Model
  TYPE_NEIGH = 3,     //!< Neighborhood
  TYPE_RULE = 4,      //!< Lithotype Rule
  TYPE_ANAM = 5,      //!< Gaussian Anamorphosis
  TYPE_TOKENS = 6,    //!< Object Definition
  TYPE_POLYGON = 7,   //!< Polygon
  TYPE_FRACTURE = 8,  //!< Fracture
  TYPE_PCA = 9,       //!< PCA Transform
} ENUM_TYPES;

typedef enum
{
  CST_NTCAR = 0,             //!< Number of characters in printout
  CST_NTDEC = 1,             //!< Number of decimal digits in printout
  CST_NTROW = 2,             //!< Maximum number of rows in table printout
  CST_NTCOL = 3,             //!< Maximum number of columns in table printout
  CST_NPROC = 4,             //!< Display the Progress Bar
  CST_LOCMOD = 5,            //!< Update Locators Option
  CST_LOCNEW = 6,            //!< Delete all similar locators
  CST_RGL = 7,               //!< RGL graphic rendition
  CST_ASP = 8,               //!< Scaling Factor
  CST_TOLINV = 9,            //!< Tolerance for matrix inversion
  CST_TOLGEN = 10,           //!< Tolerance for matrix generalized inversion
  CST_EPSMAT = 11,           //!< Tolerance value for Matrix calculation
  CST_EPSSVD = 12,           //!< Tolerance value for SVD Matrix calculation
  CST_NUMBER = 13,           //!< Maximum number of CST Enums
} ENUM_CSTS;
#endif

// These functions have been removed from new_cova
//  COV_EXP2DFACT = 20,    //!< Factorized Factorized in 2-D
//  COV_EXPFACT = 21,      //!< Factorized Exponential

#ifndef SWIG
// Internal enums (currently not exported via SWIG)
typedef enum
{
  CONV_UNIFORM = 1,         //!< Uniform
  CONV_EXPONENTIAL = 2,     //!< Exponential
  CONV_GAUSSIAN = 3,        //!< Gaussian
  CONV_SINCARD = 4,         //!< Sine Cardinal
} ENUM_CONVS;

typedef enum
{
  CONV_DIRX = 1,         //!< Along X
  CONV_DIRY = 2,         //!< Along Y
  CONV_DIRZ = 3,         //!< Along Z
  CONV_DIRXY = 4,        //!< Along XY
  CONV_DIRXYZ = 5,       //!< Along XYZ
} ENUM_CONVDIRS;
#endif
/*
typedef enum
{
  MODEL_CALCUL_NATURAL = 0,  //!< Standard Calculation
  MODEL_CALCUL_Z_Z = 10,     //!< Calculation between Z and Z
  MODEL_CALCUL_Z_GX = 20,    //!< Calculation between Z and X-Grad
  MODEL_CALCUL_Z_GY = 21,    //!< Calculation between Z and Y-Grad
  MODEL_CALCUL_GX_GX = 30,   //!< Calculation between X-Grad and X-Grad
  MODEL_CALCUL_GY_GY = 31,   //!< Calculation between Y-Grad and Y-Grad
  MODEL_CALCUL_GX_GY = 32,   //!< Calculation between X-Grad and Y-Grad
} ENUM_MODEL_CALCULS; // No more needed
*/

/*
typedef enum
{
  MODEL_PROPERTY_NONE = 0,      //!< No specific property
  MODEL_PROPERTY_CONV = 1,      //!< Convolution mode
  MODEL_PROPERTY_ANAM = 2,      //!< Anamorphosis mode
  MODEL_PROPERTY_TAPE = 3,      //!< Tapering mode
  MODEL_PROPERTY_GRAD = 4,      //!< Gradient mode
} ENUM_MODEL_PROPERTIES; // Now see EModelProperty.hpp

typedef enum
{
  NEIGH_UNIQUE = 0,        //!< Unique Neighborhood
  NEIGH_BENCH = 1,         //!< Bench Neighborhood
  NEIGH_MOVING = 2,        //!< Moving Neighborhood
  NEIGH_IMAGE = 3,         //!< Image Neighborhood
} ENUM_NEIGHS; // Now see ENeigh.hpp

typedef enum
{
  CALCUL_UNDEFINED = -1,    //!< Undefined
  CALCUL_VARIOGRAM = 0,     //!< Variogram
  CALCUL_COVARIANCE = 1,    //!< Covariance
  CALCUL_COVARIOGRAM = 2,   //!< Transitive Covariogram
  CALCUL_MADOGRAM = 3,      //!< Madogram
  CALCUL_RODOGRAM = 4,      //!< Rodogram
  CALCUL_POISSON = 5,       //!< Poisson
  CALCUL_GENERAL1 = 6,      //!< Generalized Variogram of order 1
  CALCUL_GENERAL2 = 7,      //!< Generalized Variogram of order 2
  CALCUL_GENERAL3 = 8,      //!< Generalized Variogram of order 3
  CALCUL_COVARIANCE_NC = 9, //!< Non-centered Covariance
  CALCUL_ORDER4 = 10,       //!< Order-4 Variogram
  CALCUL_TRANS1 = 11,       //!< Transition probability G12/G1
  CALCUL_TRANS2 = 12,       //!< Transition probability G12/G2
  CALCUL_BINORMAL = 13,     //!< Binormal hypothesis G12/sqrt(G1 * G2)
} ENUM_CALCUL_VARIO; // Now see ECalcVario.hpp
*/
/*
typedef enum
{
  RULE_STD = 0,     //!< Standard Lithotype Rule
  RULE_SHIFT = 1,   //!< Shift Rule
  RULE_SHADOW = 2,  //!< Shadow Rule
} ENUM_RULES;   // Now see ERule.hpp
*/
#ifndef SWIG
// Internal enums (currently not exported via SWIG)
typedef enum
{
  SHADOW_IDLE = 0,      //!< No Shadow
  SHADOW_ISLAND = 1,    //!< Island for Shadow
  SHADOW_WATER = 2,     //!< Water for Shadow
  SHADOW_SHADOW = 3,    //!< Shadow
} ENUM_SHADOWS;

typedef enum
{
  SEISMIC_NOP = 0,        //!< No Operation
  SEISMIC_FABS = 1,       //!< Absolute value
  SEISMIC_SSQRT = 2,      //!< Signed square root
  SEISMIC_SQR = 3,        //!< Square
  SEISMIC_SSQR = 4,       //!< Signed square
  SEISMIC_SIGN = 5,       //!< Signum Function
  SEISMIC_EXP = 6,        //!< Exponentiate
  SEISMIC_SLOG = 7,       //!< Signed Natural Log
  SEISMIC_SLOG10 = 8,     //!< Signed Common Log
  SEISMIC_COS = 9,        //!< Cosine
  SEISMIC_SIN = 10,       //!< Sine
  SEISMIC_TAN = 11,       //!< Tangent
  SEISMIC_COSH = 12,      //!< Hyperbolic Cosine
  SEISMIC_SINH = 13,      //!< Hyperbolic Sine
  SEISMIC_TANH = 14,      //!< Hyperbolic Tangent
  SEISMIC_NORM = 15,      //!< Divide trace by Max. Value
  SEISMIC_DB = 16,        //!< 20 * slog10 (data)
  SEISMIC_NEG = 17,       //!< Negate value
  SEISMIC_ONLY_POS = 18,  //!< Pass only positive values
  SEISMIC_ONLY_NEG = 19,  //!< Pass only negative values
  SEISMIC_SUM = 20,       //!< Running sum trace integration
  SEISMIC_DIFF = 21,      //!< Running diff trace differentiation
  SEISMIC_REFL = 22,      //!< (v[i+1] - v[i]) / (v[i+1] + v[i])
  SEISMIC_MOD_2PI = 23,   //!< Modulo 2 PI
  SEISMIC_INV = 24,       //!< Inverse
  SEISMIC_AVG = 25,       //!< Remove average value
} ENUM_SEISMICS;

typedef enum
{
  WAVELET_NONE = 0,       //!< No Convolution
  WAVELET_RICKER1 = 1,    //!< Ricker wavelet with peak frequency "fpeak"
  WAVELET_RICKER2 = 2,    //!< Ricker wavelet
  WAVELET_AKB = 3,        //!< Wavelet Alford, Kelly, and Boore
  WAVELET_SPIKE = 4,      //!< Spike
  WAVELET_UNIT = 5,       //!< Constant unit shift
} ENUM_WAVELETS;
#endif
/*
typedef enum
{
  MEMBER_LHS = 0,        //!< Left-hand Side of the Kriging System
  MEMBER_RHS = 1,        //!< Right-hand Side of the Kriging System
  MEMBER_VAR = 2,        //!< Variance of the Kriging System
} ENUM_MEMBERS;  // Now see ECalcMember.hpp

typedef enum
{
  ANAM_UNDEFINED = -1,     //!< Undefined anamorphosis
  ANAM_EXTERNAL = 0,       //!< External anamorphosis
  ANAM_HERMITIAN = 1,      //!< Hermitian anamorphosis
  ANAM_EMPIRICAL = 2,      //!< Empirical anamorphosis
  ANAM_DISCRETE_DD = 3,    //!< Discrete anamorphosis
  ANAM_DISCRETE_IR = 4,    //!< Discrete Indicator Residuals anamorphosis
} ENUM_ANAMS;   /: Now see EAnam.hpp

typedef enum
{
  CONS_UNKNOWN = 0,
  CONS_RANGE = 1,    //!< Non-stationary range
  CONS_ANGLE = 2,    //!< Non-stationary anisotropy rotation angle (degree)
  CONS_PARAM = 3,    //!< Non-stationary auxiliary parameter
  CONS_SILL = 4,     //!< Non-stationary sill
  CONS_SCALE = 5,    //!< Non-stationary scale
  CONS_T_RANGE = 6,  //!< Non-stationary tapering range
  CONS_VELOCITY = 7, //!< Non-stationary velocity (advection)
  CONS_SPHEROT = 8,  //!< Non-stationary rotation angle for Sphere
  CONS_ROTMAT = 9,   //!< Non-stationary anisotropy matrix term
} ENUM_CONS; // Now see EConsElem.hpp

typedef enum
{
  CONS_TYPE_LOWER = -1,  //!< Lower Bound
  CONS_TYPE_DEFAULT = 0, //!< Default parameter
  CONS_TYPE_UPPER = 1,   //!< Upper Bound
  CONS_TYPE_EQUAL = 2,   //!< Equality
} ENUM_CONS_TYPE; // Now see EConsType.hpp
*/
#ifndef SWIG
// Internal enums (currently not exported via SWIG)
typedef enum
{
  ANAM_QT_Z = 0,
  ANAM_QT_T = 1,
  ANAM_QT_Q = 2,
  ANAM_QT_B = 3,
  ANAM_QT_M = 4,
  ANAM_QT_PROBA = 5,
  ANAM_QT_QUANT = 6,
  ANAM_N_QT = 7,
} ENUM_ANAM_QT;
#endif
/*
typedef enum
{
  GD_J_LEFT = -1,
  GD_J_CENTER = 0,
  GD_J_RIGHT = 1,
} ENUM_GD_J;  // Now see EJustify.hpp

typedef enum
{
  PROCESS_UNDEFINED = -1,
  PROCESS_COPY = 0,
  PROCESS_MARGINAL = 1,
  PROCESS_CONDITIONAL = 2,
} ENUM_PROCESS; // Now see EProcessOper.hpp

typedef enum
{
  POPT_UNDEFINED = -1,     //!< Undefined
  POPT_ONE = 0,            //!< Power is 1
  POPT_MINUSONE = 1,       //!< Power is -1
  POPT_MINUSHALF = 2,      //!< Power is -0.5
  POPT_HALF = 3,           //!< Power is 0.5
  POPT_LOG = 4,            //!< Logarithm
} ENUM_POPTS; // Now see EPowerPT.hpp
typedef enum
{
  CALCUL_KRIGING     = 0,    //!< Kriging
  CALCUL_SIMUCOND    = 1,   //!< Conditional simulations
  CALCUL_SIMUNONCOND = 2 //!< Non conditional simulations
} ENUM_CALCUL_MODE;  // Now see ESPDECalcMode.hpp
*/

