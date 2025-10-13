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

#include "Basic/NamingConvention.hpp"

namespace gstlrn
{
class MatrixSparse;
class Model;
class Vario;
class ANeigh;
class AMesh;
class MeshEStandard;
class RuleProp;
class Cheb_Elem;
class Rule;
class VarioParam;
class AAnam;
class AnamHermite;
class Selectivity;
class DbGrid;
class NeighImage;
class EMorpho;
class MatrixSymmetric;

/****************************************/
/* Prototyping the functions in krige.c */
/****************************************/

Id _krigsim(Db* dbin,
             Db* dbout,
             const Model* model,
             ANeigh* neigh,
             bool flag_bayes,
             const VectorDouble& dmean,
             const MatrixSymmetric& dcov,
             Id icase,
             Id nbsimu,
             bool flag_dgm);

/*******************************************/
/* Prototyping the functions in variopgs.c */
/*******************************************/

Rule* _rule_auto(Db* db,
                 const VarioParam* varioparam,
                 const RuleProp* ruleprop,
                 Id ngrfmax = 1,
                 Id verbose = false);

/*****************************************/
/* Prototyping the functions in thresh.c */
/*****************************************/

Id _db_rule(Db* db,
             const RuleProp* ruleprop,
             Model* model                    = nullptr,
             const NamingConvention& namconv = NamingConvention("Facies", true, true, true, ELoc::fromKey("FACIES")));
Id _db_bounds(Db* db,
               const RuleProp* ruleprop,
               Model* model                    = nullptr,
               const NamingConvention& namconv = NamingConvention("Bounds"));
Id _db_threshold(Db* db,
                  const RuleProp* ruleprop,
                  Model* model                    = nullptr,
                  const NamingConvention& namconv = NamingConvention("Thresh"));

} // namespace gstlrn