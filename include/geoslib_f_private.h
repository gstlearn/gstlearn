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
#ifndef GEOSLIB_F_PRIVATE_H
#define GEOSLIB_F_PRIVATE_H

#include "geoslib_d_private.h"
#include "geoslib_d.h"
#include "Neigh/Neigh.hpp"
#include "Model/Model.hpp"
#include "Variogram/Vario.hpp"
#include "Mesh/MeshEStandard.hpp"
#include "Model/ANoStat.hpp"
#include "LithoRule/RuleProp.hpp"

/***********************************************/
/* Prototyping the functions in dirent_win32.c */
/***********************************************/
GEOSLIB_API DIR   *opendir(const char *dirname);
GEOSLIB_API int    closedir(DIR *dirp);
GEOSLIB_API struct dirent *readdir(DIR *dirp);

/*************************************/
/* Prototyping the functions in io.c */
/*************************************/
GEOSLIB_API int   _file_read(FILE *file,const char *format, va_list ap);
GEOSLIB_API int   _file_get_ncol(FILE *file);
GEOSLIB_API void  _file_delimitors(char del_com,char del_sep,char del_blk);
GEOSLIB_API FILE *_file_open(const char *filename,int mode);
GEOSLIB_API int   _record_read(FILE *file,const char *format,...);
GEOSLIB_API int   _buffer_read(char **buffer,const char *format,va_list ap);
GEOSLIB_API void  _file_write(FILE *file,const char *format,va_list ap);
GEOSLIB_API void  _buffer_write(char *buffer,const char *format,va_list ap);
GEOSLIB_API int   _lire_key(const char  *question,int nkeys,const char **keys);
GEOSLIB_API void  _lire_string(const char *question,int flag_def,
                               const char *valdef,char *answer);
GEOSLIB_API int   _lire_int(const char *question,int flag_def,
                            int valdef,int valmin,int valmax);
GEOSLIB_API double _lire_double(const char *question,int flag_def,
                                double valdef,double valmin,double valmax);
GEOSLIB_API int   _lire_logical(const char *question,int flag_def,int valdef);
GEOSLIB_API char *_next_file(char *dirname,char *in_string,char *ex_string);
GEOSLIB_API void  _erase_current_string(void);

/****************************************/
/* Prototyping the functions in vario.c */
/****************************************/

GEOSLIB_API double _variogram_convert_angular_tolerance(double tolang);
GEOSLIB_API int _variogram_compute(Db* db,
                                   Vario* vario,
                                   const VectorDouble& means = VectorDouble(),
                                   const VectorDouble& vars = VectorDouble(),
                                   int flag_grid = 0,
                                   int flag_gen = 0,
                                   int flag_sample = 0,
                                   int verr_mode = 0,
                                   int flag_model = 0,
                                   Model* model = nullptr,
                                   int verbose = 0);

/****************************************/
/* Prototyping the functions in krige.c */
/****************************************/

GEOSLIB_API int krigsim(const char *string,
                        Db *dbin,
                        Db *dbout,
                        Model *model,
                        Neigh *neigh,
                        double *dmean,
                        double *dcov,
                        int icase,
                        int nbsimu,
                        int flag_dgm,
                        double rval);

/*****************************************/
/* Prototyping the functions in simtub.c */
/*****************************************/

GEOSLIB_API int _gibbs_init_monovariate(int flag_order,Props *propdef,
                                        Db *db,Model *model,int isimu,
                                        int ipgs,int igrf,int nbsimu,
                                        int verbose);
GEOSLIB_API int _gibbs_init_multivar(int flag_order,Props *propdef,Db *db,
                                     Model *model,int isimu,int ipgs,
                                     int nbsimu,int verbose);

/***************************************/
/* Prototyping the functions in spde.c */
/***************************************/

GEOSLIB_API double *spde_get_mesh_dimension(MeshEStandard *amesh);
GEOSLIB_API cs     *spde_fill_S(MeshEStandard *amesh,Model* model, double *units);
GEOSLIB_API VectorDouble spde_fill_TildeC(MeshEStandard *amesh,double *units);
GEOSLIB_API VectorDouble spde_fill_Lambda(Model* model, MeshEStandard* amesh,
                                          const VectorDouble& TildeC);
GEOSLIB_API cs *spde_build_Q(cs *S,
                             const VectorDouble& Lambda,
                             int nblin,double *blin);
GEOSLIB_API Cheb_Elem *spde_cheb_duplicate(Cheb_Elem *cheb_in);
GEOSLIB_API Cheb_Elem *spde_cheb_manage(int mode,int verbose,double power,
                                        int nblin,double *blin,
                                        cs *S,Cheb_Elem *cheb_old);
GEOSLIB_API int spde_chebychev_operate(cs *S,Cheb_Elem *cheb_elem,
                                       const VectorDouble& lambda,
                                       const double *x,
                                       double *y);
GEOSLIB_API Rule *rule_auto(Db *db,
                            Vario *vario,
                            RuleProp* ruleprop,
                            int ngrfmax = 1,
                            int verbose = false);
GEOSLIB_API int db_rule(Db* db,
                        RuleProp* ruleprop,
                        Model *model = nullptr,
                        NamingConvention namconv = NamingConvention("Facies",LOC_FACIES));
GEOSLIB_API int db_bounds(Db *db,
                          RuleProp* ruleprop,
                          Model *model = nullptr,
                          NamingConvention namconv = NamingConvention("Bounds"));
GEOSLIB_API int db_threshold(Db *db,
                             RuleProp* ruleprop,
                             Model *model = nullptr,
                             NamingConvention namconv = NamingConvention("Thresh"));


/*************************************/
/* Prototyping the functions in db.c */
/*************************************/

GEOSLIB_API int db_indicator(Db *db,
                             int iatt,
                             int flag_disc,
                             const VectorDouble& mini = VectorDouble(),
                             const VectorDouble& maxi = VectorDouble(),
                             const VectorBool& incmini = VectorBool(),
                             const VectorBool& incmaxi = VectorBool(),
                             NamingConvention namconv = NamingConvention("Indicator"));
GEOSLIB_API int db_category(Db *db,
                            int ivar,
                            const VectorDouble& mini = VectorDouble(),
                            const VectorDouble& maxi = VectorDouble(),
                            const VectorBool& incmini = VectorBool(),
                            const VectorBool& incmaxi = VectorBool(),
                            NamingConvention namconv = NamingConvention("Category"));

#endif
