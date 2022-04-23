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

//#include "geoslib_d.h"

#include "Basic/NamingConvention.hpp"

class Model;
class Vario;
class ANeighParam;
class MeshEStandard;
class RuleProp;
class cs;
class Cheb_Elem;
class Rule;
class VarioParam;

/*************************************/
/* Prototyping the functions in io.c */
/*************************************/
int _file_read(FILE *file, const char *format, va_list ap);
int _file_get_ncol(FILE *file);
void _file_delimitors(char del_com, char del_sep, char del_blk);
FILE* _file_open(const char *filename, int mode);
int _record_read(FILE *file, const char *format, ...);
int _buffer_read(char **buffer, const char *format, va_list ap);
void _file_write(FILE *file, const char *format, va_list ap);
void _buffer_write(char *buffer, const char *format, va_list ap);
void _lire_string(const char *question,
                  int flag_def,
                  const char *valdef,
                  char *answer);
int _lire_int(const char *question,
              int flag_def,
              int valdef,
              int valmin,
              int valmax);
double _lire_double(const char *question,
                    int flag_def,
                    double valdef,
                    double valmin,
                    double valmax);
int _lire_logical(const char *question, int flag_def, int valdef);
void _erase_current_string(void);

/****************************************/
/* Prototyping the functions in vario.c */
/****************************************/

double _variogram_convert_angular_tolerance(double tolang);
int _variogram_compute(Db *db,
                       Vario *vario,
                       int flag_grid = 0,
                       int flag_gen = 0,
                       int flag_sample = 0,
                       int verr_mode = 0,
                       Model *model = nullptr,
                       int verbose = 0);

/****************************************/
/* Prototyping the functions in krige.c */
/****************************************/

int _krigsim(Db* dbin,
             Db* dbout,
             const Model* model,
             ANeighParam* neighparam,
             int flag_bayes,
             const VectorDouble& dmean,
             const VectorDouble& dcov,
             int icase,
             int nbsimu,
             bool flag_dgm,
             double r_coeff);

/***************************************/
/* Prototyping the functions in spde.c */
/***************************************/

double* _spde_get_mesh_dimension(MeshEStandard *amesh);
cs* _spde_fill_S(MeshEStandard *amesh, Model *model, double *units);
VectorDouble _spde_fill_TildeC(MeshEStandard *amesh, double *units);
VectorDouble _spde_fill_Lambda(Model *model,
                               MeshEStandard *amesh,
                               const VectorDouble &TildeC);
cs* _spde_build_Q(cs *S, const VectorDouble &Lambda, int nblin, double *blin);
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
             const NamingConvention& namconv = NamingConvention("Facies",
                                                         ELoc::FACIES));
int _db_bounds(Db *db,
               const RuleProp *ruleprop,
               Model *model = nullptr,
               const NamingConvention& namconv = NamingConvention("Bounds"));
int _db_threshold(Db *db,
                  const RuleProp *ruleprop,
                  Model *model = nullptr,
                  const NamingConvention& namconv = NamingConvention("Thresh"));

/******************************************/
/* Prototyping the functions in dbtools.c */
/******************************************/

int _db_indicator(Db *db,
                  int iatt,
                  int flag_disc,
                  const VectorDouble &mini = VectorDouble(),
                  const VectorDouble &maxi = VectorDouble(),
                  const VectorBool &incmini = VectorBool(),
                  const VectorBool &incmaxi = VectorBool(),
                  const NamingConvention& namconv = NamingConvention("Indicator"));
int _db_category(Db *db,
                 int ivar,
                 const VectorDouble &mini = VectorDouble(),
                 const VectorDouble &maxi = VectorDouble(),
                 const VectorBool &incmini = VectorBool(),
                 const VectorBool &incmaxi = VectorBool(),
                 const NamingConvention& namconv = NamingConvention("Category"));

