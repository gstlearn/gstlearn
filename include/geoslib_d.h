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

// WARNING: Make this include list as small as possible!
#include "geoslib_define.h"

#include "Enum/EKrigOpt.hpp"
#include "Model/Option_VarioFit.hpp"

// TODO : strcasecmp macro to be kept ?
#if defined(_WIN32) || defined(_WIN64)
#if !defined(strcasecmp)
#define strcasecmp _stricmp
#endif
#if !defined(strncasecmp)
#define strncasecmp _strnicmp
#endif
#endif

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

// TODO: Transform to a class and prevent calling malloc and memfree
typedef struct
{
  double tmin; /* Minimum abscissa along line */
  double tmax; /* Maximum abscissa along line */
  double scale; /* Scaling factor */
  double t00; /* Origin along the line */
  double dxp; /* Increment along X */
  double dyp; /* Increment along Y */
  double dzp; /* Increment along Z */
  double ang[3]; /* Angles for the line orientation */
} Direction;

typedef struct
{
  int nbands; /* Total number of bands */
  int nbtuba; /* Number of turning bands */
  int nbsimu; /* Number of simulations */
  int max_alloc; /* Maximum allocated core space */
  int nb_points_simu; /* Total number of target (point + grid) */
  double theta; /* Poisson intensity along the line */
  double field; /* Field extension */
  int *seeds; /* Array for seeds */
  int *senug; /* Array for seeds (nugget effect) */
  Direction **codir; /* Direction sub-structures */
} Situba;


typedef struct
{
  int law; /* Type of law */
  double valarg[4]; /* Randomization arguments */
} Token_Par;

typedef struct
{
  int type; /* Token type */
  int npar; /* Number of parameters */
  double factor_x2y; /* Link factor for the geometry from x to y */
  double factor_x2z; /* Link factor for the geometry from x to z */
  double factor_y2z; /* Link factor for the geometry from y to z */
  double prop; /* Token Proportion */
  std::vector<Token_Par> pars; /* Token parameter array */
} Token_Def;

typedef struct
{
  int nb_tokens; /* Number of tokens */
  std::vector<Token_Def> defs; /* Token Definition array */
} Tokens;

struct Local_Bool_Object
{
  int type; /* Type of the token */
  int seed; /* Seed for the token generation */
  struct Local_Bool_Object *address; /* Pointer to the next token */
  double center[3]; /* Coordinates of the center of the token */
  double extension[3]; /* Extension of the token */
  double orientation; /* Orientation angle for the token (radian) */
  double values[3]; /* List of additional arguments */
  double box[3][2]; /* Box containing the token */
};
typedef struct Local_Bool_Object Bool_Object;

typedef struct
{
  double center[3]; /* Location of the conditioning information */
  int nb_cover; /* Number of covering tokens */
} Bool_Cond;

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

class Dbgrid;
typedef struct
{
  int nxyz;
  int ndim;
  int nval;
  int size;
  int quant;
  int date;
  int nalloc;
  int nval_max;
  double  total;
  double  total_max;
  int    *address;
  double *energy;
  Dbgrid *db;
} Skin;

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
  double coord;                //!< Abscissas of the first Fault point
  double orient;               //!< Fault orientation
  VectorDouble thetal;               //!< Maximum density on left
  VectorDouble thetar;               //!< Maximum density on right
  VectorDouble rangel;               //!< Decrease range on left
  VectorDouble ranger;               //!< Decrease range on right
} Frac_Fault;

typedef struct
{
  double orient;              //!< Mean orientation
  double dorient;             //!< Standard deviation for orientation
  double theta0;              //!< Reference Poisson intensity
  double alpha;               //!< Power dependency between layer & intensity
  double ratcst;              //!< Ratio of Constant vs. shaped intensity
  double prop1;               //!< Survival probability (constant term)
  double prop2;               //!< Survival probability (length dependent term)
  double aterm;               //!< Survival probability (cumulative length term)
  double bterm;               //!< Survival probability (layer thickness term)
  double range;               //!< Range of fracture repulsion area
} Frac_Fam;

typedef struct
{
  int nfamilies;            //!< Number of families
  int nfaults;              //!< Number of main faults
  double xmax;                 //!< Maximum horizontal distance
  double ymax;                 //!< Maximum vertical distance
  double deltax;               //!< Dilation along the horizontal axis
  double deltay;               //!< Dilation along the vertical axis
  double xextend;              //!< Field extension along horizontal axis
  double mean;                 //!< Mean of thickness distribution
  double stdev;                //!< Standard deviation of thickness distribution
  std::vector<Frac_Fam> frac_fams; //!< Family definition (dim: nfamilies)
  std::vector<Frac_Fault> frac_faults; //!< Fault definition (dim: nfaults)
} Frac_Environ;

typedef struct
{
  int npoint;
  int family;
  double orient;
  VectorDouble xy;
} Frac_Desc;

typedef struct
{
  int nfracs;               //<! Number of fractures
  std::vector<Frac_Desc> frac_descs; //<! Array of fracture descriptions
} Frac_List;

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

typedef struct
{
  int ndim;
  int rank;
  double dist;
  double *coor;
} PL_Dist;

struct triangulateio
{
  double *pointlist; /* In / out */
  double *pointattributelist; /* In / out */
  int *pointmarkerlist; /* In / out */
  int numberofpoints; /* In / out */
  int numberofpointattributes; /* In / out */
  int *trianglelist; /* In / out */
  double *triangleattributelist; /* In / out */
  double *trianglearealist; /* In only  */
  int *neighborlist; /* Out only */
  int numberoftriangles; /* In / out */
  int numberofcorners; /* In / out */
  int numberoftriangleattributes; /* In / out */
  int *segmentlist; /* In / out */
  int *segmentmarkerlist; /* In / out */
  int numberofsegments; /* In / out */
  double *holelist; /* In / pointer to array copied out */
  int numberofholes; /* In / copied out */
  double *regionlist; /* In / pointer to array copied out */
  int numberofregions; /* In / copied out */
  int *edgelist; /* Out only */
  int *edgemarkerlist; /* Not used with Voronoi diagram; out only */
  double *normlist; /* Used only with Voronoi diagram; out only */
  int numberofedges; /* Out only */
};

struct segmentio
{
  double *pointlist; /* In / out */
  double *pointattributelist; /* In / out */
  int numberofpoints; /* In / out */
  int numberofpointattributes; /* In / out */
  int *segmentlist; /* In / out */
  int numberofsegments; /* In / out */
  int numberofcorners; /* In / out */
};

typedef struct
{
  int n_nodes; /* Number of nodes */
  int sph_size; /* Size of arrays sph_list and sph_lptr */
  double *sph_x; /* Array of X-coordinates for nodes */
  double *sph_y; /* Array of Y-coordinates for nodes */
  double *sph_z; /* Array of Z-coordinates for nodes */
  int *sph_list; /* Set of nodal indexes */
  int *sph_lptr; /* Set of pointers (sph_list indexes) */
  int *sph_lend; /* Set of pointers to adjacency lists */
} SphTriangle;

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

typedef struct
{
  int ndupl;
  int *dupl_data;
  int *dupl_dabs;
  int *dupl_grid;
} Vercoloc;

typedef struct
{
  int order;
  int nvertex;
  int ngibbs;
  int nb1;
  int nb2;
  int *vt;
  int *r_g;
  int *r_abs;
} Vertype;

class SPDE_Mesh
{
public:

  int ndim;
  int ncorner;
  int nmesh;
  int nvertex;
  int *meshes;
  double *points;
  Vercoloc *vercoloc;
  Vertype *vertype;
};

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
  SPDE_Mesh *s_mesh;
} SPDE_Matelem;

typedef struct
{
  int mesh_dbin;
  int mesh_dbout;
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

typedef struct Local_Relem Relem;
typedef struct Local_Split Split;
