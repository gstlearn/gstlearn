/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
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
class cs;
class Cheb_Elem;
class Rule;
class VarioParam;
class AAnam;
class AnamHermite;
class Selectivity;
class DbGrid;
class NeighImage;
class EMorpho;

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
             ANeigh* neigh,
             bool flag_bayes,
             const VectorDouble& dmean,
             const VectorDouble& dcov,
             int icase,
             int nbsimu,
             bool flag_dgm);
void _image_smoother(DbGrid *dbgrid,
                     const NeighImage *neigh,
                     int type,
                     double range,
                     int iptr0);
int _db_morpho_calc(DbGrid *dbgrid,
                    int iptr0,
                    const EMorpho &oper,
                    double vmin = 0.,
                    double vmax = 1.5,
                    int option = 0,
                    const VectorInt &radius = VectorInt(),
                    bool flagDistErode = false,
                    bool verbose = false);
void _morpho_angle2D(DbGrid *dbgrid, const VectorInt &radius, int iptr0);
void _morpho_gradients(DbGrid *dbgrid, int iptr0);

/***************************************/
/* Prototyping the functions in spde.c */
/***************************************/

double* _spde_get_mesh_dimension(AMesh *amesh);
cs* _spde_fill_S(AMesh *amesh, Model *model, double *units);
VectorDouble _spde_fill_TildeC(AMesh* amesh, double* units);
VectorDouble _spde_fill_Lambda(Model *model,
                               AMesh *amesh,
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
                  bool flagBelow = false,
                  bool flagAbove = false,
                  const NamingConvention& namconv = NamingConvention("Indicator"));
int _db_category(Db *db,
                 int ivar,
                 const VectorDouble &mini = VectorDouble(),
                 const VectorDouble &maxi = VectorDouble(),
                 const VectorBool &incmini = VectorBool(),
                 const VectorBool &incmaxi = VectorBool(),
                 const NamingConvention& namconv = NamingConvention("Category"));
VectorDouble _db_limits_statistics(Db *db,
                                   int iatt,
                                   const VectorDouble &mini,
                                   const VectorDouble &maxi,
                                   const VectorBool &incmini,
                                   const VectorBool &incmaxi,
                                   int optionStat,
                                   bool flagBelow,
                                   bool flagAbove);
int _migrate(Db *db1,
             Db *db2,
             int iatt1,
             int iatt2,
             int ldmax,
             const VectorDouble &dmax,
             bool flag_fill,
             bool flag_inter,
             bool flag_ball = false);

/***************************************/
/* Prototyping the functions in anam.c */
/***************************************/

int _conditionalExpectation(Db *db,
                            AAnam *anam,
                            const Selectivity *selectivity,
                            int iptr0,
                            int col_est,
                            int col_std,
                            bool flag_OK,
                            double proba,
                            int nbsimu);
int _uniformConditioning(Db *db,
                         AnamHermite *anam,
                         Selectivity *selectivity,
                         int iptr0,
                         int col_est,
                         int col_var);
