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
class MatrixSquareSymmetric;
class MatrixSparse;

/****************************************/
/* Prototyping the functions in krige.c */
/****************************************/

int _krigsim(Db* dbin,
             Db* dbout,
             const Model* model,
             ANeigh* neigh,
             bool flag_bayes,
             const VectorDouble& dmean,
             const MatrixSquareSymmetric& dcov,
             int icase,
             int nbsimu,
             bool flag_dgm);
void _image_smoother(DbGrid *dbgrid,
                     const NeighImage *neigh,
                     int type,
                     double range,
                     int iptr0);

/***************************************/
/* Prototyping the functions in spde.c */
/***************************************/

double* _spde_get_mesh_dimension(AMesh *amesh);
MatrixSparse* _spde_fill_S(AMesh *amesh, Model *model, const double *units);
VectorDouble _spde_fill_TildeC(AMesh* amesh, const double* units);
VectorDouble _spde_fill_Lambda(Model *model,
                               AMesh *amesh,
                               const VectorDouble &TildeC);
MatrixSparse* _spde_build_Q(MatrixSparse *S,
                            const VectorDouble &Lambda,
                            int nblin,
                            double *blin);
Cheb_Elem* _spde_cheb_duplicate(Cheb_Elem *cheb_in);

/*******************************************/
/* Prototyping the functions in variopgs.c */
/*******************************************/

Rule* _rule_auto(Db *db,
                 const VarioParam *varioparam,
                 const RuleProp *ruleprop,
                 int ngrfmax = 1,
                 int verbose = false);

/*****************************************/
/* Prototyping the functions in thresh.c */
/*****************************************/

int _db_rule(Db *db,
             const RuleProp *ruleprop,
             Model *model = nullptr,
             const NamingConvention& namconv = NamingConvention("Facies", true, true, true,
                                                         ELoc::fromKey("FACIES")));
int _db_bounds(Db *db,
               const RuleProp *ruleprop,
               Model *model = nullptr,
               const NamingConvention& namconv = NamingConvention("Bounds"));
int _db_threshold(Db *db,
                  const RuleProp *ruleprop,
                  Model *model = nullptr,
                  const NamingConvention& namconv = NamingConvention("Thresh"));

