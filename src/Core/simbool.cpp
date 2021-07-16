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
#include "geoslib_e.h"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"

static Bool_Object *Start_object_init,*Start_object;
static Token_Def   *Def;
static double       Origin[3],Field[3],Coor[3],Angle_Z,Theta_cste;
static int          Nb_object_init,Nb_object,Flag_stat;
static Db          *Dbout;
static int PHASE, NDIM;
static int ITER = 0;

/*! \cond */
#define EPS            1.e-3
#define MAXITER        100000
#define STOPITER      (ITER >= MAXITER)
#define NB_TOKEN_TYPES 9
#define NB_PARAM_LAWS  8
#define DEBUG          0
#define NB_FIELDS      8
#define SAVE_TAB(i,j)  (save_tab[(i) * NB_FIELDS + (j)])
/*! \endcond */

static void st_generate_type0(Bool_Object *object);
static void st_generate_type1(Bool_Object *object);
static void st_generate_type2(Bool_Object *object);
static void st_generate_type3(Bool_Object *object);
static void st_generate_type4(Bool_Object *object);
static void st_generate_type5(Bool_Object *object);
static void st_generate_type6(Bool_Object *object);
static void st_generate_type7(Bool_Object *object);
static void st_generate_type8(Bool_Object *object);
static int  st_check_type0(double dx,double dy,double dz,Bool_Object *object);
static int  st_check_type1(double dx,double dy,double dz,Bool_Object *object);
static int  st_check_type2(double dx,double dy,double dz,Bool_Object *object);
static int  st_check_type3(double dx,double dy,double dz,Bool_Object *object);
static int  st_check_type4(double dx,double dy,double dz,Bool_Object *object);
static int  st_check_type5(double dx,double dy,double dz,Bool_Object *object);
static int  st_check_type6(double dx,double dy,double dz,Bool_Object *object);
static int  st_check_type7(double dx,double dy,double dz,Bool_Object *object);
static int  st_check_type8(double dx,double dy,double dz,Bool_Object *object);
static double st_random_constant   (double *valarg);
static double st_random_uniform    (double *valarg);
static double st_random_gaussian   (double *valarg);
static double st_random_exponential(double *valarg);
static double st_random_gamma      (double *valarg);
static double st_random_stable     (double *valarg);
static double st_random_beta1      (double *valarg);
static double st_random_beta2      (double *valarg);

typedef struct {
  char name[STRING_LENGTH];
  int  npar;
  int  ind_extension_x;
  int  ind_extension_y;
  int  ind_extension_z;
  int  flag_symz;               /*  0 if the object if Z-symmetric */
  /*  1 if it is truncated upwards */
  /* -1 if it is truncated downwards */
  void (*generate_type)(Bool_Object *object);
  int  (*check_grain)(double, double, double, Bool_Object *object);
  char title[6][STRING_LENGTH];
} Def_Token;

typedef struct {
  int  law;
  char name[STRING_LENGTH];
  double (*random_generate)(double *valarg);
} Def_Rand;

static Def_Rand DEF_RAND[] = {
  {0, "Constant"    , st_random_constant},
  {1, "Uniform"     , st_random_uniform},
  {2, "Gaussian"    , st_random_gaussian},
  {3, "Exponential" , st_random_exponential},
  {4, "Gamma"       , st_random_gamma},
  {5, "Stable"      , st_random_stable},
  {6, "Beta1"       , st_random_beta1},
  {7, "Beta2"       , st_random_beta2}
};

static Def_Token DEF_TOKEN[] = { 
  {"Parallelepiped"       ,4, 0,1,2,0,st_generate_type0,st_check_type0,
   {"X-Extension","Y-Extension","Z-Extension","Orientation Angle","",""}},
  {"Lower-half Ellipsoid" ,4, 0,1,2,1,st_generate_type1,st_check_type1,
   {"X-Extension","Y-Extension","Z-Half Extension","Orientation Angle","",""}},
  {"Lower-half Sinusoid"  ,6, 3,4,-1,1,st_generate_type2,st_check_type2,
   {"Period","Amplitude","Thickness","X-Extension","Z-Half Extension",
    "Orientation Angle"}},
  {"Lower-half Paraboloid",4, 0,1,2,1,st_generate_type3,st_check_type3,
   {"X-Extension","Y-Extension","Z-Half Extension","Orientation Angle","",""}},
  {"Upper-half Ellipsoid" ,4,0,1,2,-1,st_generate_type4,st_check_type4,
   {"X-Extension","Y-Extension","Z-Half Extension","Orientation Angle","",""}},
  {"Upper-half Sinusoid"  ,6,3,4,-1,-1,st_generate_type5,st_check_type5,
   {"Period","Amplitude","Thickness","X-Extension","Z-Half Extension",
    "Orientation Angle"}},
  {"Upper-half Paraboloid",4,0,1,2,-1,st_generate_type6,st_check_type6,
   {"X-Extension","Y-Extension","Z-Half Extension","Orientation Angle","",""}},
  {"Full Ellipsoid"       ,4,0,1,2, 0,st_generate_type7,st_check_type7,
   {"X-Extension","Y-Extension","Z-Extension","Orientation Angle","",""}},
  {"Full Paraboloid"      ,4,0,1,2, 0,st_generate_type8,st_check_type8,
   {"X-Extension","Y-Extension","Z-Extension","Orientation Angle","",""}},
};

/****************************************************************************/
/*!
**  Clear the Object
**
** \return  Pointer to the newly freed Tokens structure
** 
** \param[in]  tokens Pointer to the Tokens structure to be freed
**
*****************************************************************************/
GEOSLIB_API Tokens *tokens_free(Tokens *tokens)

{
  if (tokens == (Tokens *) NULL) return(tokens);
  delete tokens;
  tokens = (Tokens *) nullptr;
  return(tokens);
}

/****************************************************************************/
/*!
**  Normalize the proportions
**
** \param[in]  tokens Tokens structure
** 
*****************************************************************************/
static void st_normalize_proportions(Tokens *tokens)

{
  int i;
  double total;

  /* Loop on the different tokens */
  
  total = 0.;
  for (i=0; i<tokens->nb_tokens; i++)
    total += tokens->defs[i].prop;

  if (ABS(total) <= 0.)
  {
    for (i=0; i<tokens->nb_tokens; i++)
      tokens->defs[i].prop = 1./ (double) tokens->nb_tokens;
  }
  else
  {
    for (i=0; i<tokens->nb_tokens; i++)
      tokens->defs[i].prop /= total;
  }
  return;
}

/****************************************************************************/
/*!
**  Create the Tokens Object
**
** \return  Pointer to the newly created Tokens structure
**
** \param[in]  nb_tokens  Number of tokens
** 
*****************************************************************************/
GEOSLIB_API Tokens *tokens_create(int nb_tokens)

{
  Tokens *tokens;

  /* Initializations */

  tokens = (Tokens *) NULL;
  if (nb_tokens <= 0) return(tokens);
  
  /* Allocate the main structure */

  tokens = new Tokens;
  tokens->nb_tokens = nb_tokens;
  tokens->defs.resize(nb_tokens);

  return tokens;
}

/****************************************************************************/
/*!
**  Create one Token Object
**
** \return  Error return code
**
** \param[in]  tokens     Tokens structure
** \param[in]  rank       Rank of the target Token
** \param[in]  type       Token type
** \param[in]  npar       Token number of parameters
** \param[in]  prop       Token proportion
** \param[in]  factor_x2y Factor of linkage between X and Y
** \param[in]  factor_x2z Factor of linkage between X and Z
** \param[in]  factor_y2z Factor of linkage between XY and Z
** \param[in]  law        Array of laws
** \param[in]  valarg     Array of randomization values
**                        (Dimension: 4 * npar)
** 
*****************************************************************************/
GEOSLIB_API int tokone_create(Tokens *tokens,
                              int     rank,
                              int     type,
                              int     npar,
                              double  prop,
                              double  factor_x2y,
                              double  factor_x2z,
                              double  factor_y2z,
                              int    *law,
                              double *valarg)
{
  if (tokens == (Tokens *) NULL) return 1;
  if (rank < 0 || rank >= tokens->nb_tokens) return 1;
  if (npar != DEF_TOKEN[type].npar) return 1;

  /* Load the Token_Def information */

  Token_Def& def = tokens->defs[rank];
  def.type       = type;
  def.prop        = prop;
  def.npar       = DEF_TOKEN[def.type].npar;
  def.factor_x2y = factor_x2y;
  def.factor_x2z = factor_x2z;
  def.factor_y2z = factor_y2z;

  /* Allocate the Token_Pars sub-structures */

  def.pars.resize(def.npar);

  /* Loop on the Token_Par sub-structures */

  for (int j=0; j<def.npar; j++)
  {
    Token_Par& param = def.pars[j];
    param.law = law[j];
    for (int k=0; k<4; k++) param.valarg[k] = valarg[k + 4 * j];
  }
  return(0);
}

/****************************************************************************/
/*!
**  Define the token type interactively
**
** \param[in] def   Token_Def structure
**
*****************************************************************************/
static void st_lire_token_type(Token_Def& def)

{
  int j,lrep;

  message("List of Token Set available\n");
  for (j=0; j<NB_TOKEN_TYPES; j++)
    message("Type %d : %s\n",j+1,DEF_TOKEN[j].name);
  lrep = _lire_int("Enter the token type",0,ITEST,1,NB_TOKEN_TYPES);
  lrep--;
  def.type = lrep;
}

/****************************************************************************/
/*!
**  Initialize a parameter
**
** \param[in,out] param   Token_Par structure
**
*****************************************************************************/
static void st_init_param(Token_Par& param)
{
  int i;

  param.law = -1;
  for (i=0; i<4; i++)  param.valarg[i] = 0.;
}

/****************************************************************************/
/*!
**  Define the parameter law interactively
**
** \param[in] param   Token_Par structure
**
*****************************************************************************/
static void st_lire_param_law(Token_Par& param)
{
  int j,lrep;

  message("List of Laws available:\n");
  for (j=0; j<NB_PARAM_LAWS; j++)
    message("Law %d : %s\n",j+1,DEF_RAND[j].name);
  lrep = _lire_int("Enter the Law for Parameter",1,1,1,NB_PARAM_LAWS);
  lrep--;
  param.law = lrep;
}

/****************************************************************************/
/*!
**  Define the randomization parameters interactively
**
** \param[in] param   Token_Par structure
**
*****************************************************************************/
static void st_lire_param_valarg(Token_Par& param)
{
  int k;
  for (k=0; k<4; k++) param.valarg[k] = 0.;

  switch (param.law)
  {
    case 0:
      param.valarg[0] = _lire_double("Constant Parameter Value",1,0.,0.,TEST);
      break;
    
    case 1:
      param.valarg[0] = _lire_double("Minimum Value",1,0.,0.,TEST);
      param.valarg[1] = _lire_double("Maximum Value",1,0.,0.,TEST);
      break;
    
    case 2:
      param.valarg[0] = _lire_double("Mean Value",1,0.,0.,TEST);
      param.valarg[1] = _lire_double("St. Dev. Value",1,0.,0.,TEST);
      break;
    
    case 3:
      param.valarg[0] = _lire_double("Mean Value",1,0.,0.,TEST);
      param.valarg[1] = _lire_double("Scale Factor",1,0.,0.,TEST);
      break;
    
    case 4:
      param.valarg[0] = _lire_double("Mean Value",1,0.,0.,TEST);
      param.valarg[1] = _lire_double("Parameter",1,0.,0.,TEST);
      break;
    
    case 5:
      param.valarg[0] = _lire_double("Alpha Parameter",1,0.,0.,TEST);
      param.valarg[1] = _lire_double("Beta Parameter",1,0.,0.,TEST);
      param.valarg[2] = _lire_double("Gamma Parameter",1,0.,0.,TEST);
      param.valarg[3] = _lire_double("Delta Parameter",1,0.,0.,TEST);
      break;
    
    case 6:
      param.valarg[0] = _lire_double("First Parameter",1,0.,0.,TEST);
      param.valarg[1] = _lire_double("Second Parameter",1,0.,0.,TEST);
      break;
    
    case 7:
      param.valarg[0] = _lire_double("First Parameter",1,0.,0.,TEST);
      param.valarg[1] = _lire_double("Second Parameter",1,0.,0.,TEST);
      break;
  }
}

/****************************************************************************/
/*!
**  Print the randomization parameters
**
** \param[in] param   Token_Par structure
**
*****************************************************************************/
static void st_print_param_valarg(Token_Par& param)
{
  switch (param.law)
  {
    case -1:
      message(" Obtained as X * %lf\n",param.valarg[0]);
      break;

    case -2:
      message(" Obtained as Y * %lf\n",param.valarg[0]);
      break;

    case 0:
      message(" Constant=%lf\n",param.valarg[0]);
      break;
    
    case 1:
      message(" Uniform within [%lf; %lf]\n",param.valarg[0],param.valarg[1]);
      break;
    
    case 2:
      message(" Gaussian - Mean=%lf - Stdv=%lf\n",param.valarg[0],param.valarg[1]);
      break;
    
    case 3:
      message(" Exponential - Mean=%lf - Scale=%lf\n",param.valarg[0],param.valarg[1]);
      break;
    
    case 4:
      message(" Gamma - Mean=%lf - Scale=%lf\n",param.valarg[0],param.valarg[1]);
      break;
    
    case 5:
      message(" Stable - Alpha=%lf - Beta=%lf - Gamma=%lf - Delta=%lf\n",
              param.valarg[0],param.valarg[1],param.valarg[2],param.valarg[3]);
      break;
    
    case 6:
      message(" Beta1 - Par1=%lf - Par2=%lf\n",param.valarg[0],param.valarg[1]);
      break;
    
    case 7:
      message(" Beta2 - Par1=%lf - Par2=%lf\n",param.valarg[0],param.valarg[1]);
      break;
  }
}

/****************************************************************************/
/*!
**  Check if the argument corresponds to an extension parameter
**  If positive, ask if there is a linkage with a previously defined extension
**
** \return  1 if the current parameter is a linked extension
**
** \param[in]  def     Token_Def structure
** \param[in]  rank    Rank of the current parameter
**
*****************************************************************************/
static int st_lire_extension_linkage(Token_Def& def,
                                     int rank)
{
  if (rank == DEF_TOKEN[def.type].ind_extension_y)
  {
    def.factor_x2y =
      _lire_double("Ratio Y/X extensions for object (0: no linkage)",1,
                   0.,0.,ITEST);
    if (def.factor_x2y > 0.)
    {
      def.pars[rank].law = -1;
      def.pars[rank].valarg[0] = def.factor_x2y;
      return(1);
    }
  }

  if (rank == DEF_TOKEN[def.type].ind_extension_z)
  {
    def.factor_x2z =
      _lire_double("Ratio Z/X extensions for object (0: no linkage)",1,
                   0.,0.,ITEST);
    if (def.factor_x2z > 0.)
    {
      def.pars[rank].law = -1;
      def.pars[rank].valarg[0] = def.factor_x2z;
      return(1);
    }
  }

  if (rank == DEF_TOKEN[def.type].ind_extension_z)
  {
    def.factor_y2z =
      _lire_double("Ratio Z/Y extensions for object (0: no linkage)",1,
                   0.,0.,ITEST);
    if (def.factor_y2z > 0.)
    {
      def.pars[rank].law = -2;
      def.pars[rank].valarg[0] = def.factor_y2z;
      return(1);
    }
  }
  return(0);
}

/****************************************************************************/
/*!
**  Create the Object interactively
**
** \return  Pointer to the newly created Tokens structure
**
*****************************************************************************/
GEOSLIB_API Tokens *tokens_input(void)

{
  Tokens    *tokens;
  int        i,j,nb_tokens;

  /* Initializations */

  tokens = (Tokens *) NULL;
  nb_tokens = _lire_int("Number of tokens",1,1,1,ITEST);
  
  /* Allocate the main structure */

  tokens = new Tokens;
  tokens->nb_tokens = nb_tokens;
  tokens->defs.resize(nb_tokens);

  /* Loop on the Token_Def structures */

  for (i=0; i<nb_tokens; i++)
  {
    Token_Def& def = tokens->defs[i];

    /* Load the Token_Def information */

    st_lire_token_type(def);
    message("\nDefinition of the parameters for the token '%s'\n",
            DEF_TOKEN[def.type].name);
    def.prop = _lire_double("- Proportion",1,1.,0.,TEST);
    def.npar = DEF_TOKEN[def.type].npar;
    def.factor_x2y = 0.;
    def.factor_x2z = 0.;
    def.factor_y2z = 0.;

    /* Allocate the Token_Pars sub-structures */

    def.pars.resize(def.npar);

    /* Loop on the Token_Par sub-structures */

    for (j=0; j<def.npar; j++)
    {
      message("\n- %s : \n",DEF_TOKEN[def.type].title[j]);
      Token_Par& param  = def.pars[j];
      st_init_param(param);
      if (st_lire_extension_linkage(def,j)) continue;
      st_lire_param_law(param);
      st_lire_param_valarg(param);
    }
  }

  /* Normalize the proportions */

  st_normalize_proportions(tokens);

  return(tokens);
}

/****************************************************************************/
/*!
**  Returns the number of parameters for the current Token
**
** \param[in]  tokens   Tokens structure
** \param[in]  rank     Rank of the target token
**
** \param[out]  type    Token type
** \param[out]  npar    Number of paramaters for the token
** \param[out]  prop    Token proportions
** 
*****************************************************************************/
GEOSLIB_API void tokone_get_nbparams(Tokens *tokens,
                                     int     rank,
                                     int    *type,
                                     int    *npar,
                                     double *prop)
{
  Token_Def& def = tokens->defs[rank];
  *type = def.type;
  *npar = def.npar;
  *prop = def.prop;

  return;
}

/****************************************************************************/
/*!
**  Returns the list of parameters of the target Token
**
** \param[in]  tokens   Tokens structure
** \param[in]  rank     Rank of the target token
**
** \param[out]  factor_x2y Y/X extension factor
** \param[out]  factor_x2z Z/X extension factor
** \param[out]  factor_y2z Z/Y extension factor
** \param[out]  law        Array of laws for arguments (Dimension: npar)
** \param[out]  valarg     Array of randomization parameters (Dimension: 4*npar)
** 
*****************************************************************************/
GEOSLIB_API void tokone_get_params(Tokens *tokens,
                                   int     rank,
                                   double *factor_x2y,
                                   double *factor_x2z,
                                   double *factor_y2z,
                                   int    *law,
                                   double *valarg)
{
  Token_Def& def = tokens->defs[rank];
  *factor_x2y = def.factor_x2y;
  *factor_x2z = def.factor_x2z;
  *factor_y2z = def.factor_y2z;

  int ecr = 0;
  for (int j=0; j<def.npar; j++)
  {
    Token_Par& param = def.pars[j];
    law[j] = param.law;
    for (int k=0; k<4; k++, ecr++)
      valarg[ecr] = param.valarg[k];
  }
  return;
}

/****************************************************************************/
/*!
**  Print information about one Token structure
**
** \param[in]  tokens   Tokens structure
** \param[in]  rank     Rank of the Token set
** 
*****************************************************************************/
GEOSLIB_API void tokone_print(Tokens *tokens,
                              int rank)
{
  if (tokens == (Tokens *) NULL) return;
  if (rank < 0 || rank >= tokens->nb_tokens) return;

  Token_Def& def = tokens->defs[rank];
  message("Token %d : %s (Nb. params=%d) - Proportion=%lf\n",rank+1,
          DEF_TOKEN[def.type].name,def.npar,def.prop);
  for (int j=0; j<def.npar; j++)
  {
    message("  %20s :",DEF_TOKEN[def.type].title[j]);
    st_print_param_valarg(def.pars[j]);
  }
  return;
}

/****************************************************************************/
/*!
**  Print information about the Tokens structure
**
** \param[in]  tokens   Tokens structure
** 
*****************************************************************************/
GEOSLIB_API void tokens_print(Tokens *tokens)

{
  int i;

  mestitle(0,"Tokens Definition");
  message("Number of tokens = %d\n",tokens->nb_tokens);

  for (i=0; i<tokens->nb_tokens; i++)
    tokone_print(tokens,i);
}

/*****************************************************************************/
/*!
**  Blank out a given object
**
** \param[in]  object characteristics of the object
**
*****************************************************************************/
static void st_blank_object(Bool_Object *object)
{
  int  i;

  /* Fixed size fields */

  for (i=0; i<3; i++)
  {
    object->center[i]    = 0.;
    object->values[i]    = 0.;
    object->extension[i] = 0.;
  }

  object->orientation = 0.;
  object->address     = (Bool_Object *) NULL;

  return;
}

/*****************************************************************************/
/*!
**  Generate a constant value
**
** \param[in] valarg Argument list
**
*****************************************************************************/
static double st_random_constant(double *valarg)
{
  return(valarg[0]);
}

/*****************************************************************************/
/*!
**  Generate a random value according to a UNIFORM distribution
**
** \param[in] valarg Argument list
**
*****************************************************************************/
static double st_random_uniform(double *valarg)
{
  double value;
  value = law_uniform(valarg[0],valarg[1]);
  return(value);
}

/*****************************************************************************/
/*!
**  Generate a random value according to a GAUSSIAN distribution
**
** \param[in] valarg Argument list
**
*****************************************************************************/
static double st_random_gaussian(double *valarg)
{
  double value;
  value = valarg[0] + valarg[1] * law_gaussian();
  return(value);
}

/*****************************************************************************/
/*!
**  Generate a random value according to an EXPONENTIAL distribution
**
** \param[in] valarg Argument list
**
*****************************************************************************/
static double st_random_exponential(double *valarg)
{
  double value;
  value = valarg[0] + valarg[1] * law_gaussian();
  return(value);
}

/*****************************************************************************/
/*!
**  Generate a random value according to an GAMMA distribution
**
** \param[in] valarg Argument list
**
*****************************************************************************/
static double st_random_gamma(double *valarg)
{
  double value;
  value = valarg[0] + law_gamma(valarg[1]);
  return(value);
}

/*****************************************************************************/
/*!
**  Generate a random value according to an STABLE distribution
**
** \param[in] valarg Argument list
**
*****************************************************************************/
static double st_random_stable(double *valarg)
{
  double value;
  value = law_stable(valarg[0],valarg[1],valarg[2],valarg[3]);
  return(value);
}

/*****************************************************************************/
/*!
**  Generate a random value according to an BETA1 distribution
**
** \param[in] valarg Argument list
**
*****************************************************************************/
static double st_random_beta1(double *valarg)
{
  double value;
  value = law_beta1(valarg[0],valarg[1]);
  return(value);
}

/*****************************************************************************/
/*!
**  Generate a random value according to an BETA2 distribution
**
** \param[in] valarg Argument list
**
*****************************************************************************/
static double st_random_beta2(double *valarg)
{
  double value;
  value = law_beta2(valarg[0],valarg[1]);
  return(value);
}

/*****************************************************************************/
/*!
**  Returns the value of the parameter according to its law 
**
** \param[in]  rank rank of the parameter 
**
*****************************************************************************/
static double st_generate_value(int rank)
{
  Token_Par& param  = Def->pars[rank];
  if (param.law < 0) return(0.);
  double value  = DEF_RAND[param.law].random_generate(param.valarg);

  return(value);
}

/****************************************************************************/
/*!
**  Check if the pixel (x,y,z) belongs to the object of type 0
**
** \return  1 if the pixel is in the grain, 0 if it is in the pore
**
** \param[in]  dx      location of the pixel along X 
** \param[in]  dy      location of the pixel along Y
** \param[in]  dz      location of the pixel along Z
** \param[in]  object  Bool_Object to be generated
**
*****************************************************************************/
static int st_check_type0(double       dx,
                          double       dy,
                          double       dz,
                          Bool_Object *object)
{
  return(1);
}

/*****************************************************************************/
/*!
**  Generate the geometry of the PARALLELEPIPED object
**
** \param[in]  object Bool_Object to be generated
**
*****************************************************************************/
static void st_generate_type0(Bool_Object *object)

{
  if (NDIM >= 1)
    object->extension[0] = st_generate_value(0);
  if (NDIM >= 2)
    object->extension[1] = st_generate_value(1);
  if (NDIM >= 3)
    object->extension[2] = st_generate_value(2);
  object->orientation  = st_generate_value(3) - Angle_Z;
  
  return;
}

/****************************************************************************/
/*!
**  Check if the pixel (x,y,z) belongs to the object of type 1
**
** \return  1 if the pixel is in the grain, 0 if it is in the pore
**
** \param[in]  dx      location of the pixel along X 
** \param[in]  dy      location of the pixel along Y
** \param[in]  dz      location of the pixel along Z
** \param[in]  object  Bool_Object to be generated
**
*****************************************************************************/
static int st_check_type1(double       dx,
                          double       dy,
                          double       dz,
                          Bool_Object *object)
{
  dx = (NDIM >= 1) ? dx / (object->extension[0] / 2.) : 0.;
  dy = (NDIM >= 2) ? dy / (object->extension[1] / 2.) : 0.;
  dz = (NDIM >= 3) ? dz / (object->extension[2]     ) : 0.;

  if (dx*dx + dy*dy + dz*dz > 1) return(0);
  return(1);
}

/*****************************************************************************/
/*!
**  Generate the geometry of the LOWER-HALF ELLIPSOID object
**
** \param[in]  object  Bool_Object to be generated
**
*****************************************************************************/
static void st_generate_type1(Bool_Object *object)

{
  if (NDIM >= 1)
    object->extension[0] = st_generate_value(0);
  if (NDIM >= 2)
    object->extension[1] = st_generate_value(1);
  if (NDIM >= 3)
    object->extension[2] = st_generate_value(2);
  object->orientation  = st_generate_value(3) - Angle_Z;
  
  return;
}

/****************************************************************************/
/*!
**  Check if the pixel (x,y,z) belongs to the object of type 2
**
** \return  1 if the pixel is in the grain, 0 if it is in the pore
**
** \param[in]  dx      location of the pixel along X 
** \param[in]  dy      location of the pixel along Y
** \param[in]  dz      location of the pixel along Z
** \param[in]  object  Bool_Object to be generated
**
*****************************************************************************/
static int st_check_type2(double       dx,
                          double       dy,
                          double       dz,
                          Bool_Object *object)
{
  double yloc;

  dx   = (NDIM >= 1) ? dx / object->values[0] : 0.;
  dz   = (NDIM >= 3) ? dz / object->extension[2] : 0.;
  yloc = object->values[1] * cos(2. * GV_PI * dx) / 2.;
  dy   = (NDIM >= 2) ? (dy - yloc) / (object->values[2] / 2.) : 0.;

  if (dz*dz + dy*dy > 1) return(0);
  return(1);
}

/*****************************************************************************/
/*!
**  Generate the geometry of the LOWER-HALF SINUSOID object
**
** \param[in]  object  Bool_Object to be generated
**
*****************************************************************************/
static void st_generate_type2(Bool_Object *object)

{
  if (NDIM >= 1)
    object->values[0]    = st_generate_value(0);
  if (NDIM >= 2)
    object->values[1]    = st_generate_value(1);
  if (NDIM >= 3)
    object->values[2]    = st_generate_value(2);
  if (NDIM >= 1)
    object->extension[0] = st_generate_value(3);
  if (NDIM >= 2)
    object->extension[1] = (object->values[1] + object->values[2]);
  if (NDIM >= 3)
    object->extension[2] = st_generate_value(4);
  object->orientation  = st_generate_value(5) - Angle_Z;
  
  return;
}

/****************************************************************************/
/*!
**  Check if the pixel (x,y,z) belongs to the object of type 3
**
** \return  1 if the pixel is in the grain, 0 if it is in the pore
**
** \param[in]  dx      location of the pixel along X 
** \param[in]  dy      location of the pixel along Y
** \param[in]  dz      location of the pixel along Z
** \param[in]  object  Bool_Object to be generated
**
*****************************************************************************/
static int st_check_type3(double       dx,
                          double       dy,
                          double       dz,
                          Bool_Object *object)
{
  dx = (NDIM >= 1) ? dx / (object->extension[0] / 2.) : 0.;
  dy = (NDIM >= 2) ? dy / (object->extension[1] / 2.) : 0.;
  dz = (NDIM >= 3) ? dz / (object->extension[2]     ) : 0.;

  if (dx*dx + dy*dy - dz > 1) return(0);
  return(1);
}

/*****************************************************************************/
/*!
**  Generate the geometry of the LOWER-HALF PARABOLOID object
**
** \param[in]  object  Bool_Object to be generated
**
*****************************************************************************/
static void st_generate_type3(Bool_Object *object)

{
  if (NDIM >= 1)
    object->extension[0] = st_generate_value(0);
  if (NDIM >= 2)
    object->extension[1] = st_generate_value(1);
  if (NDIM >= 3)
    object->extension[2] = st_generate_value(2);
  object->orientation  = st_generate_value(3) - Angle_Z;
  
  return;
}

/****************************************************************************/
/*!
**  Check if the pixel (x,y,z) belongs to the object of type 4
**
** \return  1 if the pixel is in the grain, 0 if it is in the pore
**
** \param[in]  dx      location of the pixel along X 
** \param[in]  dy      location of the pixel along Y
** \param[in]  dz      location of the pixel along Z
** \param[in]  object  Bool_Object to be generated
**
*****************************************************************************/
static int st_check_type4(double       dx,
                          double       dy,
                          double       dz,
                          Bool_Object *object)
{
  dx = (NDIM >= 1) ? dx / (object->extension[0] / 2.) : 0.;
  dy = (NDIM >= 2) ? dy / (object->extension[1] / 2.) : 0.;
  dz = (NDIM >= 3) ? dz / (object->extension[2]     ) : 0.;

  if (dx*dx + dy*dy + dz*dz > 1) return(0);
  return(1);
}

/*****************************************************************************/
/*!
**  Generate the geometry of the UPPER-HALF ELLIPSOID object
**
** \param[in]  object  Bool_Object to be generated
**
*****************************************************************************/
static void st_generate_type4(Bool_Object *object)

{
  if (NDIM >= 1)
    object->extension[0] = st_generate_value(0);
  if (NDIM >= 2)
    object->extension[1] = st_generate_value(1);
  if (NDIM >= 3)
    object->extension[2] = st_generate_value(2);
  object->orientation  = st_generate_value(3) - Angle_Z;
  
  return;
}

/****************************************************************************/
/*!
**  Check if the pixel (x,y,z) belongs to the object of type 5
**
** \return  1 if the pixel is in the grain, 0 if it is in the pore
**
** \param[in]  dx      location of the pixel along X 
** \param[in]  dy      location of the pixel along Y
** \param[in]  dz      location of the pixel along Z
** \param[in]  object  Bool_Object to be generated
**
*****************************************************************************/
static int st_check_type5(double       dx,
                          double       dy,
                          double       dz,
                          Bool_Object *object)
{
  double yloc;

  dx   = (NDIM >= 1) ? dx / object->values[0] : 0.;
  dz   = (NDIM >= 3) ? dz / object->extension[2] : 0.;
  yloc = object->values[1] * cos(2. * GV_PI * dx) / 2.;
  dy   = (NDIM >= 2) ? (dy - yloc) / (object->values[2] / 2.) : 0.;

  if (dz*dz + dy*dy > 1) return(0);
  return(1);
}

/*****************************************************************************/
/*!
**  Generate the geometry of the UPPER-HALF SINUSOID object
**
** \param[in]  object  Bool_Object to be generated
**
*****************************************************************************/
static void st_generate_type5(Bool_Object *object)

{
  if (NDIM >= 1)
    object->values[0]    = st_generate_value(0);
  if (NDIM >= 2)
    object->values[1]    = st_generate_value(1);
  if (NDIM >= 3)
    object->values[2]    = st_generate_value(2);
  if (NDIM >= 1)
    object->extension[0] = st_generate_value(3);
  if (NDIM >= 2)
    object->extension[1] = (object->values[1] + object->values[2]);
  if (NDIM >= 3)
    object->extension[2] = st_generate_value(4);
  object->orientation  = st_generate_value(5) - Angle_Z;
  
  return;
}

/****************************************************************************/
/*!
**  Check if the pixel (x,y,z) belongs to the object of type 6
**
** \return  1 if the pixel is in the grain, 0 if it is in the pore
**
** \param[in]  dx      location of the pixel along X 
** \param[in]  dy      location of the pixel along Y
** \param[in]  dz      location of the pixel along Z
** \param[in]  object  Bool_Object to be generated
**
*****************************************************************************/
static int st_check_type6(double       dx,
                          double       dy,
                          double       dz,
                          Bool_Object *object)
{
  dx = (NDIM >= 1) ? dx / (object->extension[0] / 2.) : 0.;
  dy = (NDIM >= 2) ? dy / (object->extension[1] / 2.) : 0.;
  dz = (NDIM >= 3) ? dz / (object->extension[2]     ) : 0.;

  if (dx*dx + dy*dy + dz > 1) return(0);
  return(1);
}

/*****************************************************************************/
/*!
**  Generate the geometry of the UPPER-HALF PARABOLOID object
**
** \param[in]  object  Bool_Object to be generated
**
*****************************************************************************/
static void st_generate_type6(Bool_Object *object)

{
  if (NDIM >= 1)
    object->extension[0] = st_generate_value(0);
  if (NDIM >= 2)
    object->extension[1] = st_generate_value(1);
  if (NDIM >= 3)
    object->extension[2] = st_generate_value(2);
  object->orientation  = st_generate_value(3) - Angle_Z;
  
  return;
}

/****************************************************************************/
/*!
**  Check if the pixel (x,y,z) belongs to the object of type 7
**
** \return  1 if the pixel is in the grain, 0 if it is in the pore
**
** \param[in]  dx      location of the pixel along X 
** \param[in]  dy      location of the pixel along Y
** \param[in]  dz      location of the pixel along Z
** \param[in]  object  Bool_Object to be generated
**
*****************************************************************************/
static int st_check_type7(double       dx,
                          double       dy,
                          double       dz,
                          Bool_Object *object)
{
  dx = (NDIM >= 1) ? dx / (object->extension[0] / 2.) : 0.;
  dy = (NDIM >= 2) ? dy / (object->extension[1] / 2.) : 0.;
  dz = (NDIM >= 3) ? dz / (object->extension[2] / 2.) : 0.;

  if (dx*dx + dy*dy + dz*dz > 1) return(0);
  return(1);
}

/*****************************************************************************/
/*!
**  Generate the geometry of the FULL ELLIPSOID object
**
** \param[in]  object  Bool_Object to be generated
**
*****************************************************************************/
static void st_generate_type7(Bool_Object *object)

{
  if (NDIM >= 1)
    object->extension[0] = st_generate_value(0);
  if (NDIM >= 2)
    object->extension[1] = st_generate_value(1);
  if (NDIM >= 3)
    object->extension[2] = st_generate_value(2);
  object->orientation  = st_generate_value(3) - Angle_Z;
  
  return;
}

/****************************************************************************/
/*!
**  Check if the pixel (x,y,z) belongs to the object of type 8
**
** \return  1 if the pixel is in the grain, 0 if it is in the pore
**
** \param[in]  dx      location of the pixel along X 
** \param[in]  dy      location of the pixel along Y
** \param[in]  dz      location of the pixel along Z
** \param[in]  object  Bool_Object to be generated
**
*****************************************************************************/
static int st_check_type8(double       dx,
                          double       dy,
                          double       dz,
                          Bool_Object *object)
{
  dx = (NDIM >= 1) ? dx / (object->extension[0] / 2.) : 0.;
  dy = (NDIM >= 2) ? dy / (object->extension[1] / 2.) : 0.;
  dz = (NDIM >= 3) ? dz / (object->extension[2] / 2.) : 0.;

  if (dx*dx + dy*dy - dz > 1) return(0);
  if (dx*dx + dy*dy + dz > 1) return(0);
  return(1);
}

/*****************************************************************************/
/*!
**  Generate the geometry of the FULL PARABOLOID object
**
** \param[in]  object  Bool_Object to be generated
**
*****************************************************************************/
static void st_generate_type8(Bool_Object *object)

{
  if (NDIM >= 1)
    object->extension[0] = st_generate_value(0);
  if (NDIM >= 2)
    object->extension[1] = st_generate_value(1);
  if (NDIM >= 3)
    object->extension[2] = st_generate_value(2);
  object->orientation  = st_generate_value(3) - Angle_Z;
  
  return;
}

/*****************************************************************************/
/*!
**  Check if the pixel (xyz) belongs to the grain.
**  Perform the rotation (if necessary) before calling the specific
**  grain-dependent procedure 
**
** \return  1 if the pixel is in the grain, 0 if it is in the pore 
**
** \param[in]  xyz     location of the pixel 
** \param[in]  object  characteristics of the object
**
*****************************************************************************/
static int st_check_object(double       xyz[3],
                           Bool_Object *object)
{
  double angle,sint,cost,dx,dy,dz,dxr,dyr;
  int    memo,value;

  dx = xyz[0] - object->center[0];
  dy = xyz[1] - object->center[1];
  dz = xyz[2] - object->center[2];

  if (object->orientation)
  {
    angle = ut_deg2rad(object->orientation);
    sint  = sin(angle);
    cost  = cos(angle);
    dxr   = dx*cost + dy*sint;
    dyr   = dy*cost - dx*sint;
    dx    = dxr;
    dy    = dyr;
  }

  /* Check if the grain is outside the box */

  if (ABS(dx) > object->extension[0] / 2.) return(0);
  if (ABS(dy) > object->extension[1] / 2.) return(0);

  if (NDIM > 2)
    switch (DEF_TOKEN[object->type].flag_symz)
    {
      case 0:
        if (ABS(dz) > object->extension[2] / 2.) return(0);
        break;
      
      case 1:
        if (dz > 0) return(0);
        if (ABS(dz) > object->extension[2]     ) return(0);
        break;
      
      case -1:
        if (dz < 0) return(0);
        if (ABS(dz) > object->extension[2]     ) return(0);
        break;
    }

  /* Check the pixel according to the grain definition */

  memo  = law_get_random_seed();
  value = DEF_TOKEN[object->type].check_grain(dx,dy,dz,object);
  law_set_random_seed(memo);

  return (value);
}

/*****************************************************************************/
/*!
**  Check if an object may be generated according to the value
**  of the Intensity 
**  This Intensity can be local or not (if Flag_stat)
**
** \return  Error return code: 
** \return  0 the token is created 
** \return  1 the token may not be created 
**
*****************************************************************************/
static int st_check_intensity(void)

{
  int    i,ind[3],iech;
  double theta;

  /* Initializations */

  if (Flag_stat)
    theta = Theta_cste;
  else
  {
    for (i=0; i<3; i++)
      ind[i] = (int) ((Coor[i] - Dbout->getX0(i)) / Dbout->getDX(i));
    iech  = db_index_grid_to_sample(Dbout,ind);
    theta = Dbout->getProportion(iech,0);
  }
  
  return (law_uniform(0.,1.) > theta);
}

/*****************************************************************************/
/*!
**  Function to link the geometries of an object
**
** \param[in,out] object     Bool_Object to be filled 
**
*****************************************************************************/
static void st_extension_linkage(Bool_Object *object)
{
  if (Def->factor_x2y > 0.)
    object->extension[1] = object->extension[0] * Def->factor_x2y;
  if (Def->factor_x2z > 0.)
    object->extension[2] = object->extension[0] * Def->factor_x2z;
  if (Def->factor_y2z > 0.)
    object->extension[2] = object->extension[1] * Def->factor_y2z;
}

/*****************************************************************************/
/*!
**  Print the characteristics of an object
**
** \param[in]  title    Title for the printout
** \param[in]  object   Bool_Object description of the object
**
*****************************************************************************/
static void st_print_object(const char *title,
                            Bool_Object *object)

{
  int idim;

  if (object == (Bool_Object *) NULL) return;

  message("Characteristics of the %s Object\n",title);
  message("- Type        = %d\n",object->type);
  message("- Seed        = %d\n",object->seed);
  message("- Center      =");
  for (idim=0; idim<NDIM; idim++) message(" %lf",object->center[idim]);
  message("\n");
  message("- Extension   =");
  for (idim=0; idim<NDIM; idim++) message(" %lf",object->extension[idim]);
  message("\n");
  message("- Orientation = %lf\n",object->orientation);
  message("  Current count of objects = %d\n",Nb_object);
}

/*****************************************************************************/
/*!
**  Print the list of all objects
**
*****************************************************************************/
static void st_print_all_objects(void)
{
  Bool_Object *object;

  mestitle(1,"List of retained tokens");
  object = Start_object;
  while (object)
  {
    st_print_object("current",object);
    object = object->address;
  }
}

/*****************************************************************************/
/*!
**  Function used to generate the geometry of an object
**
** \return  Error return code: 
** \return   0 : no problem 
** \return   1 : the grain cannot be generated because of Intensity 
** \return  -1 : this type of grain does not exist 
**
** \param[in]  cdgrain    characteristics of the conditioning grain 
** \li                    if NULL the object must be drawn at random within
**                        the field 
** \li                    otherwise, it must be randomized within the
**                        extension of the object
** \param[in]  tokens     Description of the tokens
**
** \param[out] object     Bool_Object to be filled 
**
*****************************************************************************/
static int st_generate_object(Bool_Cond   *cdgrain,
                              Tokens      *tokens,
                              Bool_Object *object)
{
  double value,cumul,dx,dy,dz,sint,cost,angle,valrand,total;
  int    i,itoken,flag_fix;

  /* Blank out all the parameters of the buffer object */

  flag_fix = (cdgrain != (Bool_Cond *) NULL);
  valrand  = 0.;
  st_blank_object(object);

  /* Define the (primary) location of the object */

  if (flag_fix)
  {
    for (i=0; i<NDIM; i++) Coor[i] = cdgrain->center[i];
  }
  else
  {
    do
    {
      ITER++;
      for (i=0; i<NDIM; i++)
        Coor[i] = Origin[i] + Field[i] * law_uniform(0.,1.);
    } while (st_check_intensity() && ! STOPITER);
    if (STOPITER) return(1);
  }

  /* Calculate the total probability */

  total = 0.;
  for (itoken=0; itoken<tokens->nb_tokens; itoken++)
    total += tokens->defs[itoken].prop;
  if (total <= 0.) return(-1);

  /* Find the type of token to be generated */

  value = total * law_uniform(0.,1.);
  cumul = 0.;
  for (itoken=0; itoken<tokens->nb_tokens; itoken++)
  {
    cumul += tokens->defs[itoken].prop;
    if (value < cumul) break;
  }
  Def = &tokens->defs[itoken];

  /* Fill the object characteristics */
  
  object->seed = law_get_random_seed();
  object->type = Def->type;
  DEF_TOKEN[object->type].generate_type(object);

  /* Operate the linkage */

  st_extension_linkage(object);

  /* Store the coordinates of the object center */

  if (flag_fix)
  {
    do {
      ITER++;
      for (i=0; i<NDIM; i++)
      {
        if (i < 2)
        {
          valrand = law_uniform(0.,1.) - 0.5;
        }
        else
        {
          switch (DEF_TOKEN[object->type].flag_symz)
          {
            case 0:
              valrand = law_uniform(0.,1.) - 0.5;
              break;
              
            case 1:
              valrand = law_uniform(0.,1.);
              break;
              
            case -1:
              valrand = -law_uniform(0.,1.);
              break;
          }
        }
        object->center[i] = Coor[i] + object->extension[i] * valrand;
      }
    } while (! st_check_object(cdgrain->center,object) && ! STOPITER);
    if (STOPITER) return(1);
  }
  else
    for (i=0; i<NDIM; i++) object->center[i] = Coor[i];

  /* Determine the inclusive box */

  if (ABS(object->orientation) < EPS)
  {
    dx = object->extension[0];
    dy = object->extension[1];
    dz = object->extension[2];
  }
  else

  {
    angle = ut_deg2rad(object->orientation);
    sint  = ABS(sin(angle));
    cost  = ABS(cos(angle));
    dx    = cost * object->extension[0] + sint * object->extension[1];
    dy    = sint * object->extension[0] + cost * object->extension[1];
    dz    = object->extension[2];
  }

  object->box[0][0] = object->center[0] - dx / 2;
  object->box[0][1] = object->center[0] + dx / 2;
  object->box[1][0] = object->center[1] - dy / 2;
  object->box[1][1] = object->center[1] + dy / 2;
  switch (DEF_TOKEN[object->type].flag_symz)
  {
    case 0:
      object->box[2][0] = object->center[2] - dz / 2;
      object->box[2][1] = object->center[2] + dz / 2;
      break;

    case 1:
      object->box[2][0] = object->center[2] - dz;
      object->box[2][1] = object->center[2];
      break;

    case -1:
      object->box[2][0] = object->center[2];
      object->box[2][1] = object->center[2] + dz;
      break;
  }

  if (PHASE == 2) ITER = 0;
  return(0);
}

/*****************************************************************************/
/*!
**  Check if the pixel (xyz) belongs to the object bounding box 
**
** \return  1 if the pixel is in the box, 0 otherwise 
**
** \param[in]  xyz     location of the pixel 
** \param[in]  object  characteristics of the object
**
*****************************************************************************/
static int st_check_box(double       xyz[3],
                        Bool_Object *object)

{
  int i;

  for (i=0; i<NDIM; i++)
  {
    if (xyz[i] < object->box[i][0]) return(0);
    if (xyz[i] > object->box[i][1]) return(0);
  }
  return(1);
}

/*****************************************************************************/
/*!
**  Check if the current object is compatible with the constraining
**  pores 
**
** \return  Error return code: 1 for incompatibility; 0 otherwise 
**
** \param[in]  object  Bool_Object structure describing the object
** \param[in]  nbpore  count of constraining pores 
** \param[in]  cdpore  Bool_Cond describing the constraining pores 
**
*****************************************************************************/
static int st_check_pore(Bool_Object *object,
                         int          nbpore,
                         Bool_Cond   *cdpore[])
{
  int  i;

  for (i=0; i<nbpore; i++)
  {
    if (! st_check_box(cdpore[i]->center,object)) continue;
    law_set_random_seed(object->seed);
    if (st_check_object(cdpore[i]->center,object)) return(1);
  }
  return(0);
}

/*****************************************************************************/
/*!
**  Check if an object can be added or deleted with regards to the
**  constraining grains 
**
** \return  Error return code: 1  if impossible; 0 otherwise 
**
** \param[in]  object   Bool_Object describing object to be operated 
** \param[in]  nbgrain  number of constraining grains 
** \param[in]  cdgrain  Bool_Cond describing the conditioning grains 
** \param[in]  val      type of the operation to be tested 
**                      1 for addition; -1 for deletion 
**
*****************************************************************************/
static int st_ask_grain(Bool_Object *object,
                        int          nbgrain,
                        Bool_Cond   *cdgrain[],
                        int          val)
{
  int i;

  /* Dispatch */

  if (val < 0)
  {

    /* Deleting the object */

    for (i=0; i<nbgrain; i++)
    {
      if (! st_check_box(cdgrain[i]->center,object)) continue;
      if (cdgrain[i]->nb_cover > 1) continue;
      if (st_check_object(cdgrain[i]->center,object)) return(1);
    }
  }
  else
  {

    /* Adding the object */

    for (i=0; i<nbgrain; i++)
    {
      if (! st_check_box(cdgrain[i]->center,object)) continue;
      if (! st_check_object(cdgrain[i]->center,object)) continue;
    }
  }

  return(0);
}

/*****************************************************************************/
/*!
**  Update the covering value of each constraining grain after a
**  deletion or an addition operation 
**
** \return  Count of grains not covered after the operation 
**
** \param[in]  object   Bool_Object describing token to be operated 
** \param[in]  nbgrain  number of constraining grains 
** \param[in]  cdgrain  Bool_Cond describing the conditioning grains 
** \param[in]  val      type of the operation to be tested 
**                      1 for addition; -1 for deletion 
**
*****************************************************************************/
static int st_cover_update(Bool_Object *object,
                           int          nbgrain,
                           Bool_Cond   *cdgrain[],
                           int          val)
{
  int i ,not_covered;

  /* Dispatch */

  for (i=not_covered=0; i<nbgrain; i++)
  {
    if (st_check_box(cdgrain[i]->center,object))
    {
      law_set_random_seed(object->seed);
      if (st_check_object(cdgrain[i]->center,object))
      {
        if (val < 0)
        {

          /* Deletion */

          cdgrain[i]->nb_cover--;
        }
        else
        {

          /* Addition */

          cdgrain[i]->nb_cover++;
        }
      }
    }
    if (cdgrain[i]->nb_cover <= 0) not_covered++;
  }
  return(not_covered);
}

/*****************************************************************************/
/*!
**  Add the object to the chain and update the covering value of each
**  grain 
**
** \return  Error return code
**
** \param[in]  mode     0 for primary objects; 1 for secondary objects 
** \param[in]  nbgrain  number of conditioning grains 
** \param[in]  cdgrain  Bool_Cond structures describing the grains 
** \param[in]  object   Bool_Object description of the object to be added 
**
** \param[out] not_covered count of token not covered 
**
*****************************************************************************/
static int st_add_object(int          mode,
                         int          nbgrain,
                         Bool_Cond   *cdgrain[],
                         Bool_Object *object,
                         int         *not_covered)
{
  Bool_Object *new_object;

  new_object = (Bool_Object *) mem_alloc(sizeof(Bool_Object),0);
  if (new_object == (Bool_Object *) NULL) return(1);

  (void) memcpy((char *) new_object,(char *) object,sizeof(Bool_Object));

  if (mode == 0)
  {
    new_object->address = Start_object_init;
    Start_object_init  = new_object;
    Nb_object_init    += 1;
  }
  else
    new_object->address = Start_object;

  Start_object  = new_object;
  Nb_object    += 1;

  /* Update the covering value */

  *not_covered = st_cover_update(object,nbgrain,cdgrain,1);

  /* Printout the object */
  
  if (DEBUG) st_print_object("added",object);

  /* Blank out the object (as it is now part of the chain) */

  st_blank_object(object);

  return(0);
}

/*****************************************************************************/
/*!
**  Attempts to delete an object (primary or secondary) 
**
** \return  Error return code: 1 if the object cannot be deleted; 0 otherwise 
**
** \param[in]  mode     0 for primary tokens; 1 for secondary tokens 
** \param[in]  nbgrain  number of conditioning grains 
** \param[in]  cdgrain  Bool_Cond structures describing the grains 
**
*****************************************************************************/
static int st_delete_object(int        mode,
                            int        nbgrain,
                            Bool_Cond *cdgrain[])
{
  Bool_Object *object,*last_object;
  int          i, count, rank;

  if (mode == 0)
  {
    object = Start_object_init;
    count  = Nb_object_init;
  }
  else
  {
    object = Start_object;
    count  = Nb_object;
  }
  if (count <= 0) return(1);

  /* Search for the object to be deleted */

  last_object = (Bool_Object *) NULL;
  rank = (int) (count * law_uniform(0.,1.));

  for (i=0; i<rank; i++)
  {
    last_object = object;
    object = object->address;
    if (object == (Bool_Object *) NULL) 
      messageAbort("st_delete_object; Error #1");
  }

  /* Check if the object can be deleted */

  if (st_ask_grain(object,nbgrain,cdgrain,-1)) return(1);

  /* Delete the object */

  if (st_cover_update(object,nbgrain,cdgrain,-1)) 
    messageAbort("st_delete_object : Error #2");

  /* Update the pointed address of the previous object */

  if (last_object == (Bool_Object *) NULL)
  {
    if (mode == 0 && Start_object != Start_object_init)
    {

      /* In the case of the last primary object deleted with mode = 0   */
      /* when some secondary objects have already been defined          */
      /* Search for the secondary object ancestor of the current object */

      last_object = Start_object;
      while (last_object->address != Start_object_init)
        last_object = last_object->address;
    }
  }
  if (last_object != (Bool_Object *) NULL)
    last_object->address = (rank == count - 1) ?
      (Bool_Object *) NULL : object->address;

  /* Update the global pointers (if necessary) */

  if (object == Start_object)      Start_object      = object->address;
  if (object == Start_object_init) Start_object_init = object->address;

  if (mode == 0)
    Nb_object_init--;
  else
    if (rank >= Nb_object - Nb_object_init) Nb_object_init--;
  Nb_object--;

  /* Printout the deleted object */

  if (DEBUG) st_print_object("deleted",object);
  
  /* Free the current object */

  object = (Bool_Object *) mem_free((char *) object);

  return(0);
}

/*****************************************************************************/
/*!
**  Project the objects on the output grid
**
** \param[in]  background    Value for the background
** \param[in]  facies        Value of the facies assigned
** \param[in]  iptr_simu     Pointer for storing the Boolean simulation (or <0)
** \param[in]  iptr_rank     Pointer for storing the object ranks (or <0)
**
*****************************************************************************/
static void st_project_objects(double background,
                               double facies,
                               int    iptr_simu,
                               int    iptr_rank)
{
  Bool_Object *object;
  int    i,ix,iy,iz,ix0,ix1,iy0,iy1,iz0,iz1,iad,nb_write,flag_save,ecr,rank;
  int    ind[3];
  double xyz[3],*save_tab;
    
  /* Initializations */
  
  save_tab  = (double *) NULL;
  flag_save = (int) get_keypone("Boolean_Save",0);
  if (flag_save)
    save_tab = (double *) mem_alloc(sizeof(double) * NB_FIELDS * Nb_object,1);

  /* Initialize the output result */

  for (i=0; i<get_NECH(Dbout); i++)
  {
    if (iptr_simu >= 0) Dbout->setArray(i,iptr_simu,background);
    if (iptr_rank >= 0) Dbout->setArray(i,iptr_rank,TEST);
  }

  /* Loop on the objects */

  rank   = 0;
  object = Start_object;
  while (object)
  {

    /* Look for the nodes in the box of influence of the object */

    ix0 = (int) ((object->box[0][0] - Dbout->getX0(0)) / Dbout->getDX(0) - 1);
    ix0 = MAX(ix0, 0);
    ix1 = (int) ((object->box[0][1] - Dbout->getX0(0)) / Dbout->getDX(0) + 1);
    ix1 = MIN(ix1, Dbout->getNX(0) - 1);
    iy0 = (int) ((object->box[1][0] - Dbout->getX0(1)) / Dbout->getDX(1) - 1);
    iy0 = MAX(iy0, 0);
    iy1 = (int) ((object->box[1][1] - Dbout->getX0(1)) / Dbout->getDX(1) + 1);
    iy1 = MIN(iy1, Dbout->getNX(1) - 1);
    iz0 = (int) ((object->box[2][0] - Dbout->getX0(2)) / Dbout->getDX(2) - 1);
    iz0 = MAX(iz0, 0);
    iz1 = (int) ((object->box[2][1] - Dbout->getX0(2)) / Dbout->getDX(2) + 1);
    iz1 = MIN(iz1, Dbout->getNX(2) - 1);

    /* Check the pixels within the box */

    nb_write = 0;
    for (ix=ix0; ix<=ix1; ix++)
      for (iy=iy0; iy<=iy1; iy++)
        for (iz=iz0; iz<=iz1; iz++)
        {
          xyz[0] = Dbout->getX0(0) + ix * Dbout->getDX(0);
          xyz[1] = Dbout->getX0(1) + iy * Dbout->getDX(1);
          xyz[2] = Dbout->getX0(2) + iz * Dbout->getDX(2);
          if (! st_check_object(xyz,object)) continue;
          ind[0] = ix;
          ind[1] = iy;
          ind[2] = iz;
          iad = db_index_grid_to_sample(Dbout,ind);

          /* Bypass writing if the cell is masked off */

          if (! Dbout->isActive(iad)) continue;
          
          /* Set the values */

          nb_write++;
          if (iptr_simu >= 0)
            Dbout->setArray(iad,iptr_simu,facies);
          if (iptr_rank >= 0)
          {
            if (FFFF(get_ARRAY(Dbout,iad,iptr_rank)))
              Dbout->setArray(iad,iptr_rank,(double) (rank+1));
          }
        }

    /* Save the object description (optional) */

    if (flag_save)
    {
      ecr = 0;
      SAVE_TAB(rank,ecr++) = (double) object->type;
      SAVE_TAB(rank,ecr++) = object->orientation;
      for (int idim=0; idim<3; idim++)
        SAVE_TAB(rank,ecr++) = object->center[idim];
      for (int idim=0; idim<3; idim++)
        SAVE_TAB(rank,ecr++) = object->extension[idim];
    }
    
    object = object->address;
    if (nb_write > 0) rank++;
  }

  /* Save the results in keypair mechanism */
  
  if (flag_save) set_keypair("Boolean_Tab",1,rank,NB_FIELDS,save_tab);

  save_tab = (double *) mem_free((char *) save_tab);
  return;
}

/*****************************************************************************/
/*!
**  Print the references of a grain
**
** \param[in]  cdgrain       Bool_Cond structure
**
*****************************************************************************/
static void st_print_grain(Bool_Cond *cdgrain)
{
  if (NDIM == 2)
    messerr("Grain : Center = (%lf %lf)",
            cdgrain->center[0],
            cdgrain->center[1]);
  else
    messerr("Grain : Center = (%lf %lf %lf)",
            cdgrain->center[0],
            cdgrain->center[1],
            cdgrain->center[2]);
}

/*****************************************************************************/
/*!
**  Performs the boolean simulation 
**
** \return  Error return code
**
** \param[in]  dbin          Db structure containing the data (optional)
** \param[in]  dbout         Db structure containing the simulated grid
** \param[in]  tokens        Tokens structure 
** \param[in]  seed          Seed for the random number generator
** \param[in]  nb_average    Average number of boolean objects
** \param[in]  flag_stat     1 if the Intensity is constant
** \param[in]  flag_simu     Store the boolean simulation
** \param[in]  flag_rank     Store the object rank
** \param[in]  background    Value assigned to the background
** \param[in]  facies        Value of the facies assigned
** \param[in]  dilate        Array of dilation radius (optional)
** \param[in]  theta_cste    Intensity constant value 
** \param[in]  tmax          Maximum time 
** \param[in]  verbose       1 for a verbose output
**
*****************************************************************************/
GEOSLIB_API int simbool_f(Db     *dbin,
                        Db     *dbout,
                        Tokens *tokens,
                        int     seed,
                        int     nb_average,
                        int     flag_stat,
                        int     flag_simu,
                        int     flag_rank,
                        double  background,
                        double  facies,
                        double *dilate,
                        double  theta_cste,
                        double  tmax,
                        int     verbose)
{
  Bool_Object  object;
  Bool_Cond  **cdgrain,**cdpore;
  int          iech,i,j,draw_more,iref,nbgrain,nbpore,memo_init,rank;
  int          error,status,iptr_simu,iptr_rank;
  double       tabtime,coor[3],data;

  /* Initializations */

  Start_object_init = Start_object = (Bool_Object *) NULL;
  Nb_object_init = Nb_object = memo_init = nbgrain = nbpore = 0;
  iptr_simu = iptr_rank = -1;
  tabtime = 0.;
  st_blank_object(&object);
  cdgrain = cdpore = (Bool_Cond **) NULL;

  /* Preliminary checks */

  error = 3;
  if (! is_grid(dbout))
  {
    messerr("The output Db file must be a grid");
    return(1);
  }
  NDIM = dbout->getNDim();

  /* Add the attributes for storing the simulation */

  if (flag_simu)
  {
    iptr_simu = dbout->addFields(1,background);
    if (iptr_simu < 0) goto label_end;
  }
  if (flag_rank)
  {
    iptr_rank = dbout->addFields(1,TEST);
    if (iptr_rank < 0) goto label_end;
  }

  /* Define the global variables */

  law_set_random_seed(seed);
  Dbout      = dbout;
  Flag_stat  = flag_stat;
  Theta_cste = theta_cste;
  Angle_Z    = 0.;
  for (i=0; i<NDIM; i++)
  {
    Origin[i] = dbout->getX0(i) - dbout->getDX(i) / 2.;
    if (dilate != (double *) NULL) Origin[i] -= dilate[i];
    Field[i]  = dbout->getDX(i) * dbout->getNX(i);
    if (dilate != (double *) NULL) Field[i] += 2. * dilate[i];
  }
  if (verbose)
  {
    if (dbin == (Db *) NULL)
      message("Boolean non conditional simulation. Average of %d objects\n",
              nb_average);
    else
      message("Boolean conditional simulation. Average of %d objects\n",
              nb_average);
  }

  /* Count the number of conditioning pores and grains */

  if (dbin != (Db *) NULL)
  {
    for (iech=0; iech<get_NECH(dbin); iech++)
    {
      for (i=0; i<NDIM; i++)
      {
        coor[i] = get_IDIM(dbout,iech,i);
        if (coor[i] < Origin[i] || coor[i] > Origin[i] + Field[i]) continue;
      }
      data = dbin->getVariable(iech,0);
      if (FFFF(data)) continue;
      if (data)
        nbgrain++;
      else
        nbpore++;
    }
    if (verbose)
      message("Conditioning data: %d grains and %d pores\n",nbgrain,nbpore);

    /* Store the constraining grain characteristics */
    
    if (nbgrain > 0)
    {
      cdgrain = (Bool_Cond **) mem_alloc(nbgrain * sizeof(Bool_Cond *),0);
      if (cdgrain == (Bool_Cond **) NULL) goto label_end;
      for (i=0; i<nbgrain; i++) cdgrain[i] = (Bool_Cond *) NULL;
      for (i=0; i<nbgrain; i++)
      {
        cdgrain[i] = (Bool_Cond *) mem_alloc(sizeof(Bool_Cond),0);
        if (cdgrain[i] == (Bool_Cond *) NULL) goto label_end;
      }
    }
    
    /* Store the constraining pore characteristics */
    
    if (nbpore > 0)
    {
      cdpore = (Bool_Cond **) mem_alloc(nbpore * sizeof(Bool_Cond *),0);
      if (cdpore == (Bool_Cond **) NULL) goto label_end;
      for (i=0; i<nbpore; i++) cdpore[i] = (Bool_Cond *) NULL;
      for (i=0; i<nbpore; i++)
      {
        cdpore[i] = (Bool_Cond *) mem_alloc(sizeof(Bool_Cond),0);
        if (cdpore[i] == (Bool_Cond *) NULL) goto label_end;
      }
    }
    
    for (iech=nbgrain=nbpore=0; iech<get_NECH(dbin); iech++)
    {
      for (i=0; i<NDIM; i++)
      {
        coor[i] = get_IDIM(dbin,iech,i);
        if (coor[i] < Origin[i] || coor[i] > Origin[i] + Field[i]) continue;
      }
      data = dbin->getVariable(iech,0);
      if (FFFF(data)) continue;
      if (data)
      {
        cdgrain[nbgrain]->nb_cover  = 0;
        for (i=0; i<NDIM; i++) cdgrain[nbgrain]->center[i] = coor[i];
        nbgrain++;
      }
      else
      {
        cdpore[nbpore]->nb_cover    = 0;
        for (i=0; i<NDIM; i++) cdpore[nbpore]->center[i] = coor[i];
        nbpore++;
      }
    }
  }

  /*******************************/
  /* Simulate the Initial grains */
  /*******************************/

  if (verbose)
  {
    mestitle(1,"Simulating the initial tokens");
    message("- Number of grains to be covered = %d\n",nbgrain);
  }
  draw_more = nbgrain;
  error = 1;
  ITER  = 0;
  PHASE = 1;
  while (draw_more)
  {
    ITER++;
    if (STOPITER)
    {
      messerr("Simulation of the initial objects failed after %d iterations",
              MAXITER);
      messerr("to cover %d of the %d grains",draw_more, nbgrain);
      for (i=0; i<nbgrain; i++)
        if (cdgrain[i]->nb_cover <= 0) st_print_grain(cdgrain[i]);
      messerr("Check the Token definition or the Intensity value(s)");
      goto label_end;
    }

    /* Debugging statement */

    if (debug_query("converge"))
      message("Initial grain iteration %d: Number of non covered grains = %d\n",
              ITER,draw_more);

    /* Look for a non-covered grain */

    rank = (int) (draw_more * law_uniform(0.,1.));
    for (iref=j=0; iref<nbgrain && j<=rank; iref++)
      if (cdgrain[iref]->nb_cover <= 0) j++;
    iref--;

    /* Generate an object covering the grain(x,y,z) */

    status = st_generate_object(cdgrain[iref],tokens,&object);
    if (status > 0) continue;
    if (status < 0)
    {
      messerr("Grain #%d cannot be covered",iref+1);
      st_print_grain(cdgrain[iref]);
      goto label_end;
    }

    /* Check if the object is compatible with the constraining pores */

    if (st_check_pore(&object,nbpore,cdpore)) continue;

    /* Check if object is compatible with the constraining grains */

    if (st_ask_grain(&object,nbgrain,cdgrain,1)) continue;

    /* Add the grain */

    if (st_add_object(0,nbgrain,cdgrain,&object,&draw_more)) goto label_end;
  }
  memo_init = Nb_object_init;
  if (verbose)
    message("- Number of Initial Tokens = %d\n",Nb_object_init);

  /*********************************/
  /* Simulate the Secondary grains */
  /*********************************/

  if (verbose)
  {
    mestitle(1,"Simulating the secondary tokens");
    message("- Maximum time available = %lf\n",tmax);
  }
  
  error = 2;
  ITER = 0;
  PHASE = 2;
  tabtime = 0.;
  while (tabtime < tmax)
  {
    if (STOPITER) break;
    /* The next line is not correct but is kept for compatibility.
       The correct version should be implemented on next case study 
       update.
    */
    tabtime += law_exponential()  / (nb_average + Nb_object);
    /* This should be the right version */
    /*    tabtime += law_exponential(); */

    if (law_uniform(0.,1.) <=
        ((double) nb_average / (double) (nb_average + Nb_object)))
    {

      /* Add an object */

      if (st_generate_object((Bool_Cond *) NULL,tokens,&object)) continue;

      /* Check if the object is compatible with the constraining pores */

      if (st_check_pore(&object,nbpore,cdpore)) continue;

      /* Check if the object is compatible with the constraining grains */

      if (st_ask_grain(&object,nbgrain,cdgrain,1)) continue;

      /* Add the object */

      if (st_add_object(1,nbgrain,cdgrain,
                        &object,&draw_more)) goto label_end;
    }
    else
    {

      /* Delete a primary object */

      if (! st_delete_object(0,nbgrain,cdgrain)) continue;

      /* Delete a secondary object */

      if (st_delete_object(1,nbgrain,cdgrain)) continue;
    }

    if (debug_query("converge"))
      message("At time=%lf (< tmax=%lf), %d objects have been simulated\n",
              tabtime,tmax,Nb_object);
  }

  /* Print the list of retained tokens */

  if (DEBUG) st_print_all_objects();
  
  /******************************************/
  /* Project the objects on the output grid */
  /******************************************/

  st_project_objects(background,facies,iptr_simu,iptr_rank);

  error = 0;

label_end:
  if (! error && verbose)
  {
    message("\n"
            "Boolean simulation\n"
            "==================\n");

    if (dbin != (Db *) NULL && get_NECH(dbin) > 0)
    {
      message("Conditioning option               = YES\n");
      message("Number of conditioning grains     = %d\n", nbgrain);
      message("Number of conditioning pores      = %d\n", nbpore);
    }
    else
      message("Conditioning option               = NO\n");

    message("Ending simulation time            = %g\n", tmax);
    message("Average number of objects         = %d\n", nb_average);

    if (dbin != (Db *) NULL && get_NECH(dbin) > 0)
    {
      message("Initial number of primary objects = %d\n", memo_init);
      message("Ending number of primary objects  = %d\n", Nb_object_init);
    }
    message("Total number of objects           = %d\n", Nb_object);
  }

  /* Free memory */

  st_blank_object(&object);

  if (nbgrain > 0)
    for (i=0; i<nbgrain; i++)
      cdgrain[i] = (Bool_Cond *) mem_free((char *) cdgrain[i]);
  cdgrain = (Bool_Cond **) mem_free((char *) cdgrain);
  if (nbpore > 0)
    for (i=0; i<nbpore; i++)
      cdpore[i] = (Bool_Cond *) mem_free((char *) cdpore[i]);
  cdpore = (Bool_Cond **) mem_free((char *) cdpore);

  return(error);
}

/****************************************************************************/
/*!
**  Returns the number of parameters for the current Token type
**
** \return Number of parameters
**
** \param[in]   type    Token type (starting from 0)
**
*****************************************************************************/
GEOSLIB_API int toktype_get_nbparams(int type)
{
  if (type < 0 || type >= NB_TOKEN_TYPES)
  {
    messerr("The token type must vary between 0 and %d",NB_TOKEN_TYPES-1);
    return(0);
  }
  return(DEF_TOKEN[type].npar);
}

