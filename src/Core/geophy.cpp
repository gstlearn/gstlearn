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

/*---------------------------------------------------------------------------*/
/*  TIME_3DB: FINITE DIFFERENCE COMPUTATION OF FIRST ARRIVAL TIMES IN 3D.    */
/*---------------------------------------------------------------------------*/
/*  P.PODVIN, Geophysique, Ecole des Mines de Paris, Fontainebleau.          */
/*  e-mail: Pascal.Podvin@ensmp.fr          Tel: 33-(1) 64 69 49 25.         */
/*                                                                           */
/*  Version bidouillee (18Dec2007). Voir timegrad_3db.c                      */
/*                                                                           */
/*  First release: April 1991; Previous revision: 7 May 1993                 */
/*  A version with improved portability was released 19 January 2004         */
/*  Previous update was dated 16 February 2006 (with corrections in the      */
/*  initialization procedure, with improved consistency when source point    */
/*  is at the immediate vicinity of a velocity heterogeneity)                */
/*  This version is dated 16 October 2006 (two patches, seek '16102006')     */
/*                                                                           */
/* int time_3db(double *hs, double *t, int nx, int ny, int nz,               */
/*              int bx, int by, int bz,                                      */
/*              double xs, double ys, double zs, double eps_init,            */
/*              int messages);                                               */
/*                                                                           */
/* ARGUMENTS (C description; all FORTRAN arguments are pointers)             */
/*                                                                           */
/*       (int)     nx,ny,nz        : dimensions of the time field (number of */
/*                                   grid points). Cells are cubic. No       */
/*                                   dimension may be lower than 2.          */
/*                                                                           */
/*       (int)     bx,by,bz        : args controlling the form factor of the */
/*                                   computing box. Use 1,1,1 if you don't   */
/*                                   know what this is about !!!             */
/*                                                                           */
/*       (double *) hs,t            : 1D arrays of nx*ny*nz elements.         */
/*                                   t will contain the computed time field  */
/*                                   arranged as a succession of planes x    */
/*                                   (resp. z)=Cst, each plane consisting of */
/*                                   a suite of y=Cst vectors when this      */
/*                                   routine is called from a program in C   */
/*                                   (resp. in FORTRAN). hs contains model   */
/*                                   data (slowness*mesh spacing) organized  */
/*                                   in the same manner. Within the routine, */
/*                                   values stored in hs are implicitly      */
/*                                   ascribed to cell centers, whereas times */
/*                                   computed and stored in t refer to cell  */
/*                                   corners (grid points). Cells located at */
/*                                   coordinates x=nx-1, y=ny-1 or z=nz-1    */
/*                                   (x=nx, y=ny, z=nz in FORTRAN dialect)   */
/*                                   are thus dummy cells (out of the model).*/
/*                                   The corresponding values will be        */
/*                                   ignored and treated as "infinite".      */
/*                              Note:                                        */
/*                                   Values stored in hs may not be negative */
/*                                   and must not exceed 0.499e+19 (see      */
/*                                   macro FD_HUGE defined below).           */
/*                                                                           */
/*       (double)   xs,ys,zs      :  Point source coordinates referred to    */
/*                                   "upper left" corner of model (grid      */
/*                                   point with lowest indices), expressed   */
/*                                   in mesh spacing (h) unit.               */
/*                                   Licit ranges: [0.0,nx-1.][0.0,ny-1.]    */
/*                                   [0.,nz-1.]                              */
/*                                   If source point is found to be out      */
/*                                   of the licit range, the timefield is    */
/*                                   treated as already initialized in the   */
/*                                   calling program (extended source, e.g., */
/*                                   exploding reflector). In such case, the */
/*                                   timefield must not be uniform, and      */
/*                                   nodes where times are computed must be  */
/*                                   initialized with a value certainly      */
/*                                   higher than their final arrival time.   */
/*                              Note: although this is not required when the */
/*                                    source is licit, you should always     */
/*                                    initialize array t with whatever       */
/*                                    constant value (why not 0?), so that   */
/*                                    you'll be warned if a wrong (illicit)  */
/*                                    source location is entered.            */
/*                                                                           */
/*       (double)   eps_init      :  tolerance on relative inhomogeneity     */
/*                                   used (only) during initialization.      */
/*                                   Although licit values are in [0,1],     */
/*                                   relevant values should be <0.01.        */
/*                                   0.001 is a reasonable choice.           */
/*                                   (arg. ignored if source is illicit)     */
/*                                                                           */
/*       (int)     messages        : 0 makes the routine silent (except on   */
/*                                   diagnosed error); 1: small talk mode;   */
/*                                   >1: debug mode (verbose).               */
/*                                   (a negative value is treated as 1).     */
/*                                                                           */
/* VALUE :                                                                   */
/*                                                                           */
/*       time_3d() returns a nonzero value if an error was detected.         */
/*       An explicit error message is printed.                               */
/*                                                                           */
/* CALLING TIME_2D FROM A PROGRAM WRITTEN IN FORTRAN :                       */
/*                                                                           */
/*       The C-function time_3d_() is provided as the interface to Fortran.  */
/*       In this routine, all arguments are pointers as required by Fortran  */
/*       (where routine arguments are passed by address, not by value), and  */
/*       dimensions 'x','z' are swapped to mimic standard Fortran memory     */
/*       mapping (hs[i][j][k] in C "means" HS(K,J,I) in Fortran).            */
/*       Normally, calling TIME_3D (without the trailing underscore) in      */
/*       Fortran will automatically result in a call to time_3d_() in C.     */
/*    Compiling :                                                            */
/*       With Sun compilers, nothing particular is required to obtain this   */
/*       automatic behaviour.                                                */
/*       With the GNU compilers (gcc, g77), you will need to compile your    */
/*       Fortran program with option -fno-second-underscore to obtain this.  */ 
/*       Because the C/FORTRAN interface is not standardized, you might      */
/*       have to browse your documentation on other platforms.               */
/*       If you get into trouble :                                           */
/*       Compiling this program with the option -DDEBUG_ARGS allows to check */
/*       what C-function is actually called (time_3d or time_3d_) and how    */
/*       its arguments are interpreted.                                      */
/*     Declarations in the calling program :                                 */
/*       As seen from Fortran, TIME_3D is an INTEGER FUNCTION.               */
/*       It should thus be declared as such in the calling program.          */
/*       Please note that arguments passed to TIME_3D MUST be declared with  */
/*       their correct types (e.g., XS,YS,ZS MUST NOT be declared INTEGER,   */
/*       while NX,NY,NZ MUST NOT be declared REAL !)                         */
/*       Not conforming to this may result into incomprehensible situations  */
/*       (once again, compile option -DDEBUG_ARGS may help...)               */
/*       All scalar arguments are read-only and may thus be declared         */
/*       constant (using PARAMETER statements).                              */
/*     Program skeleton :                                                    */
/*       INTEGER NX,NY,NZ                                                    */
/*       REAL HS(NX,NY,NZ),T(NX,NY,NZ)                                       */
/* C or  REAL HS(NX*NY*NZ),T(NX*NY*NZ)                                       */
/*       REAL XS,YS,ZS,EPS_INIT                                              */
/*       INTEGER MESSAGES,TIME_3D,STATUS                                     */
/*       .....                                                               */
/*       STATUS=TIME_3D(HS,T,NX,NY,NZ,XS,YS,ZS,EPS_INIT,MESSAGES)            */
/*       IF(STATUS.NE.0)                                                     */
/*         STOP "time_3d diagnosed a (presumably fatal) error"               */
/*       ENDIF                                                               */
/*       .....                                                               */
/*       STOP                                                                */
/*       END                                                                 */
/*                                                                           */
/* RECENT UPDATES :                                                          */
/*                                                                           */
/*        Although time_3d has been used by dozens of people worldwide       */
/*        for many years, some changes were recently introduced (2003).      */
/*        These changes do not affect the routine's usage (except for        */
/*        some values returned on error, and the treatment of "infinite"     */
/*        slowness values, now forbidden in hs). Their only justification    */
/*        is improved portability.                                           */
/*        The changes are the following :                                    */
/*        - The program now conforms to ANSI-C (you must include fdtimes.h   */
/*          in your calling program, if it is written in C).                 */
/*        - I decided to drop Sun-specific calls to routines implementing    */
/*          the IEEE standards for infinity and related tests. This is non   */
/*          portable, and would even create problems on Sun platforms when   */
/*          the calling program is in Fortran !                              */
/*          As a consequence, a finite threshold is now used to test whether */
/*          a double is treated as infinite. No value in array hs is allowed */
/*          to exceed this threshold (see macro FD_HUGE below).              */
/*        - Unpredictible crashes were seldom encountered on Intel based PCs */
/*          due to the fact that the routine's behaviour strongly depended   */
/*          on tests comparing doubles in an "exact" manner (e.g., replacing */
/*          a test x<y by x<=y altered significantly the code behaviour).    */
/*          These tests are now "fuzzified" with no significant consequence  */
/*          on precision (using macro EPS_FUZZY below).                      */
/*                                                                           */
/*        In 2005, Ari Trygvasson and Bjorn Bergman signaled an unexpected   */
/*        asymmetry of results in a particular case. The initialization      */
/*        procedure was found to provide wrong results when the source was   */
/*        located at the immediate vicinity of a velocity heterogeneity      */
/*        (with significant but not severe impact on precision). A quite     */
/*        significantly revised version was released on 2 February 2006.     */
/*                                                                           */
/* REFERENCE : Podvin & Lecomte, Geophys.J.Intern. 105, 271-284, 1991.       */
/*----------------------------------------------------------------------------*/

/*! \cond */
#ifndef M_SQRT2
#define M_SQRT2     1.41421356237309504880
#endif
#ifndef M_SQRT3
#define M_SQRT3     1.732050807568877076
#endif
#define MIN3(x,y,z) (MIN(x,MIN(y,z)))
#define MIN4(x,y,z,t) (MIN(x,MIN(y,MIN(z,t))))
#define    NINT(x)     (int)floor((x)+0.5)
/* NINT is "lazy nearest integer", used here only with positive values */

#define    T3D_INF     0.500e+19
#define    FD_HUGE     0.499e+19
#define    ISINF(x)    ((x)>FD_HUGE)
/* Note: T3D_INF must be lower than the SQUARE ROOT of the highest   */
/* acceptable double value (machine dependent) to prevent overflow.  */
/* FD_HUGE must be chosen such that T3D_INF>FD_HUGE is true.         */
/* This last condition is actually tested by the program !           */

#define EPS_FUZZY	1.2e-07
/* slightly more than 1/2^23 (4-byte-double mantissas are encoded on 23 bits) */

#define SMALLTALK        messages
#define VERBOSE          messages>1

/*-------------------------------------Static functions-----------------------*/

static int
pre_init(void),
  init_point(void),                         /* routine modified in Dec 2005 */
  recursive_init(void),
  propagate_point(int),
  x_side(int,int,int,int,int,int),
  y_side(int,int,int,int,int,int),
  z_side(int,int,int,int,int,int),
  scan_x_ff(int,int,int,int,int,int,int),
  scan_x_fb(int,int,int,int,int,int,int),
  scan_x_bf(int,int,int,int,int,int,int),
  scan_x_bb(int,int,int,int,int,int,int),
  scan_y_ff(int,int,int,int,int,int,int),
  scan_y_fb(int,int,int,int,int,int,int),
  scan_y_bf(int,int,int,int,int,int,int),
  scan_y_bb(int,int,int,int,int,int,int),
  scan_z_ff(int,int,int,int,int,int,int),

  scan_z_fb(int,int,int,int,int,int,int),
  scan_z_bf(int,int,int,int,int,int,int),
  scan_z_bb(int,int,int,int,int,int,int);
/* the only fully commented "side" functions are x_side() and scan_x_ff() */

static void
error(int),
  init_nearest(void),                       /* routine modified in Dec 2005 */
  init_cell(double,double,double,int,int,int), /* routine modified in Dec 2005 */
  free_ptrs(int);

static double
/* new function init_cellh(): see patches 271205[1] and 271205[2] */
init_cellh(double vh, double vv, double hsc, double hsn),
  exact_delay(double,double,double,int,int,int);


static int
t_1d(int,int,int,double,double,double,double,double),
  t_2d(int,int,int,double,double,double,double),
  diff_2d(int,int,int,double,double,double),
  t_3d_(int,int,int,double,double,double,double,double,int),
  t_3d_part1(int,int,int,double,double,double,double),
  point_diff(int,int,int,double,double),
  edge_diff(int,int,int,double,double,double);

/*-------------------------------------Static variables-----------------------*/

/* MODEL */

static int
nmesh_x,nmesh_y,nmesh_z;          /* Model dimensions (cells) */
static double
***hs,*hs_buf,                    /* 1D and 3D arrays */
  *hs_keep=(double *)NULL;           /* to save boundary values */

/* TIMEFIELD */

static int
nx,ny,nz;                         /* Timefield dimensions (nodes) */
static double
***t,*t_buf;                      /* 1D and 3D arrays */
static double
timeshift;                        /* required by "fuzzy tests" */
/* for more comments, see init_point() */

/* SOURCE */

static double
fxs,fys,fzs;                      /* Point source coordinates */
static int
xs,ys,zs;                         /* Nearest node */

/* PARAMETERS */

#ifndef INIT_MIN
#define INIT_MIN    7
#endif /* INIT_MIN */
#define N_INIT_X    (4*INIT_MIN+3)
#define N_INIT      N_INIT_X*N_INIT_X*N_INIT_X
/* This conventional value defines  */
/* the maximum size of the box that */
/* will be initialized recursively. */
/* (ADJUSTABLE at compile time).    */
/* Default value is reasonable.     */
/* Cost is already high: 3584 to    */
/* 26416 more points according to   */
/* source position.                 */

#ifndef INIT_RECURS_LIMIT
#define INIT_RECURS_LIMIT    1
#endif /* INIT_RECURS_LIMIT */
/* This parameter defines the maximal */
/* level of recursivity during init.  */
/* (ADJUSTABLE at compile time).      */
/* Value zero would practically rest- */
/* rain initialization to the source  */
/* point in heterogeneous models.     */
/* Value 2 is advisable only when     */
/* VERY severe heterogeneities are    */
/* located close to the source point. */
/*! \endcond */

static int
messages,                         /* message flag (0:silent)              */
  source_at_node=0,                 /* are source coordinate int's ? (0/1)  */
  init_stage=0,                     /* level of recursivity during init.    */
  current_side_limit,               /* actual boundary of computations      */
  X0,X1,Y0,Y1,Z0,Z1,                /* inclusive boundaries of timed region */
  reverse_order,                    /* level of recursivity in FD scheme    */
  *longflags,                       /* local headwave flags.                */
  flag_fb,x_start_fb,y_start_fb,z_start_fb,
  flag_bf,x_start_bf,y_start_bf,z_start_bf,
  flag_ff,x_start_ff,y_start_ff,z_start_ff,
  flag_bb,x_start_bb,y_start_bb,z_start_bb;
/* control current side scanning.       */

static double
hs_eps_init;                      /* tolerance on homogeneity
                                     (fraction of slowness at source point) */

static int bx,by,bz; /* hack dated 18 Dec 2007 */

/*------------------------------------------------Error flags---------------*/

/*! \cond */
#define ERROR_FREE         0
#define ERR_INFBUG       (-1)
#define ERR_MULTUNI      (-2)
#define ERR_MALLOC       (-3)
#define ERR_RECURS       (-4)
#define ERR_EPS          (-5)
#define ERR_RANGE        (-6)
#define ERR_PHYS         (-7)
#define ERR_DIM          (-8)
#define ERR_BOX          (-9)
/*! \endcond */

static const  char *err_msg[]=
{
  "\ntime_3d: Computations terminated normally.\n",
  "\ntime_3d: [Bug] macros T3D_INF, FD_HUGE not properly set.\n",
  "\ntime_3d: Multiple source but input time map is uniform.\n",
  "\ntime_3d: Not enough virtual memory (malloc failed).\n",
  "\ntime_3d: Recursive init failed: see message(s) above.\n",
  "\ntime_3d: [Init] Illegal tolerance on inhomogeneity.\n",
  "\ntime_3d: Illegal ('infinite') value(s) in array hs.\n",
  "\ntime_3d: Non-physical negative value(s) in array hs.\n",
  "\ntime_3d: a dimension (nx,ny,nz) is too small or negative.\n",
  "\ntime_3d: args bx,by,bx not correctly set. use 1,1,1 !!!\n",
};

/*-------------------------------------------------Error()------------------*/

static void error(int flag)

{
  if(messages || flag) messerr("%s",err_msg[-flag]);
}

/*-------------------------------------------------Time_3d()----------------*/

GEOSLIB_API int time_3db(double *HS, double *T, int NX, int NY, int NZ,
                         int BX, int BY, int BZ,
                         double XS, double YS, double ZS, 
                         double HS_EPS_INIT, int MSG)
{
  int    signal;
   
#ifdef DEBUG_ARGS
  messerr("******** time_3db: Option DEBUG_ARGS is on.");
  if(init_stage) messerr("Recursively entering ");
  else messerr("Initially entering ");
  messerr("time_3db() in C-style, using 'time_3db'.");
  messerr("Array dimensions nx=%d ny=%d nz=%d",NX,NY,NZ);
  messerr("Args hs,t are seen as arrays[nx][ny][nz], \
i.e. arrays[%d][%d][%d])",NX,NY,NZ);
  messerr("Licit src coordinate ranges: xs in [0.,%d] ys in [0.,%d] \
zs in [0.,%d]",NX-1,NY-1,NZ-1);
  messerr("Requested source location xs=%g ys=%g zs=%g",XS,YS,ZS);
  messerr("Other Args: EPS_INIT=%g MESSAGES=%d",
          HS_EPS_INIT,MSG);
  messerr("First elements of input arrays: HS[0][0][0]=%g, \
T[0][0][0]=%g",*HS,*T);
  messerr("******** time_3db: Option DEBUG_ARGS done.");
#endif

  /* This section merely copies arguments to internal variables. */
  /* This is where you might change things in order to build     */
  /* a customized calling interface (with an alternate arglist). */
  /* If you choose to do so, you must keep time_3d() unchanged,  */
  /* as it is needed internally (recursive init). Design your    */
  /* customized interface as another function (having another    */
  /* name...) and prototype it in fdtimes.h                      */
  
  hs_buf=HS;
  t_buf=T;
  nx=NX;
  ny=NY;
  nz=NZ;
  bx=BX;
  by=BY;
  bz=BZ;
  if((bx*by*bz!=bx && bx*by*bz!=by && bx*by*bz!=bz) || bx*by*bz<=0){
    error(ERR_BOX);
    return ERR_BOX;
  }
  fxs=XS;
  fys=YS;
  fzs=ZS;
  hs_eps_init=HS_EPS_INIT;
  if(MSG<0) messages=1;
  else messages=MSG;
  /* You should change nothing below this  */

  if((signal=pre_init())==ERROR_FREE){
    signal=propagate_point(init_point());
    free_ptrs(nx);
  }
  if(init_stage==0 || signal!=ERROR_FREE) error(signal);
  return signal;
}

/*------------------------------------------------Pre_init()----------------*/

static int pre_init(void)

{
  int
    x,y,z,
    np,nt,
    n0,n1,errtest;
  double
    *pf;
  
  if(nx<2 || ny<2 || nz<2) return ERR_DIM;
  if(!(ISINF(T3D_INF))) return ERR_INFBUG;
  /* if you encounter this error, it probably means you played */
  /* around with the values of macros FD_HUGE and T3D_INF !   */
  
  nmesh_x=nx-1;
  nmesh_y=ny-1;
  nmesh_z=nz-1;
  np=ny*nz;
  nt=nx*np;
  n1=MAX(nx,ny);
  n0=MAX(nx,nz);
  if(n1==n0) n0=MAX(ny,nz);
  n1*=n0;
  
  /* allocate pointers */
  if(!(hs=(double ***)malloc((unsigned)nx*sizeof(double **))))
    return ERR_MALLOC;
  if(!(t =(double ***)malloc((unsigned)nx*sizeof(double **)))){
    free((char *)hs);
    return ERR_MALLOC;
  }
  if(!(longflags =(int *)malloc((unsigned)n1*sizeof(int)))){
    free((char *)t);
    free((char *)hs);
    return ERR_MALLOC;
  }/* size of the largest side of the model */
  for(x=0;x<nx;x++)
    if(   !(hs[x]=(double **)malloc((unsigned)ny*sizeof(double *))) ||
          !(t[x] =(double **)malloc((unsigned)ny*sizeof(double *)))) {
      timeshift=0.0; /* possibly uninitialized */
      free_ptrs(x);
      return ERR_MALLOC;
    }
  for(x=0;x<nx;x++)
    for(y=0;y<ny;y++){
      hs[x][y]=hs_buf+x*np+y*nz;
      t[x][y]=t_buf+x*np+y*nz;
    }
  
  /* stop here if recursive call */
  if(init_stage) return ERROR_FREE;
  
  /* initialize all times as T3D_INF if licit point source */
  if(fxs>=0.0 && fxs<=nx-1 && fys>=0 && fys<=ny-1 && fzs>=0 && fzs<=nz-1)
    for(x=0,pf=t_buf;x<nt;x++) *pf++=T3D_INF;
  
  /* assign T3D_INF to hs in dummy meshes (x=nmesh_x|y=nmesh_y|z=nmesh_z) */
  /* and keep masked values in hs_keep[] (will be restored in free_ptrs()) */
  x=((nx+1)*(ny+1)+(nx+1)*nz+nz*ny)*sizeof(double);
  if(!(hs_keep=(double *)malloc((unsigned)x))) {
    timeshift=0.0; /* possibly uninitialized */
    free_ptrs(nx);
    return ERR_MALLOC;
  }
  pf=hs_keep;
  for(x=0;x<nx;x++){
    for(y=0;y<ny;y++) {
      *pf++=hs[x][y][nmesh_z];
      hs[x][y][nmesh_z]=T3D_INF;
    }
    for(z=0;z<nmesh_z;z++) {
      *pf++=hs[x][nmesh_y][z];
      hs[x][nmesh_y][z]=T3D_INF;
    }
  }
  for(y=0;y<nmesh_y;y++)
    for(z=0;z<nmesh_z;z++) {
      *pf++=hs[nmesh_x][y][z];
      hs[nmesh_x][y][z]=T3D_INF;
    }
  
  /* test for non-masked negative or infinite slowness value */
  errtest=ERROR_FREE;
  for(x=0;x<nmesh_x;x++)
    for(y=0;y<nmesh_y;y++)
      for(z=0;z<nmesh_z && errtest==ERROR_FREE;z++){
        if(ISINF(hs[x][y][z])) errtest=ERR_RANGE;
        if(hs[x][y][z]<0.0) errtest=ERR_PHYS;
      }
  if(errtest!=ERROR_FREE){
    timeshift=0.0; /* possibly uninitialized */
    free_ptrs(nx);
  }
  
  return errtest;
}

/*------------------------------------------------Init_point()--------------*/

static int init_point(void)

{
  int
    signal=ERROR_FREE,
    x,y,z,
    xsc,ysc,zsc, /* three variables added : patches of December 2005 */
    test,
    t_X0,t_X1,t_Y0,t_Y1,t_Z0,t_Z1;
  double
    min_t, max_t,
    hs0=0.0,           /* initialization required by gcc -Wall, unused */
    allowed_delta_hs,
    dist;

  /* if illicit src location, locate minimum time source point and return */
  if(fxs<0.0 || fxs>nx-1 || fys<0 || fys>ny-1 || fzs<0 || fzs>nz-1){
    /* patch 16102006[2]: bug: xs,ys,zs remained uninitialized in this loop */
    for(x=0,min_t=max_t=t[0][0][0],xs=ys=zs=0;x<nx;x++)
      for(y=0;y<ny;y++)
        for(z=0;z<nz;z++){
          if(t[x][y][z]<min_t){
            min_t=t[x][y][z];
            xs=x;
            ys=y;
            zs=z;
          }
          if(t[x][y][z]>max_t)
            max_t=t[x][y][z];
        }
    if(min_t==max_t) return ERR_MULTUNI;
    source_at_node=1;
    /* FUZZIFIED COMPARISONS: fuzzy tests take it for granted that times are*/
    /* non-negative. This is why global variable timeshift was introduced...*/
    /* If needed, time-shifting must be done only once at init_stage zero,  */
    /* and finally undone at return time (see free_ptrs()).                 */
    if(init_stage==0){
      if(min_t<0.0){
        timeshift=min_t;
        for(x=0;x<nx;x++)
          for(y=0;y<ny;y++)
            for(z=0;z<nz;z++)
              t[x][y][z]-=timeshift;
      }/* shift all times to work with non-negative values */
      else timeshift=0.0;
    }
    if(SMALLTALK)
      message("\nMultiple source starting at node [%d,%d,%d] at time %g.",
              xs,ys,zs,min_t);
    X0=X1=xs;
    Y0=Y1=ys;
    Z0=Z1=zs;
    return ERROR_FREE; /* no-op: this is a preinitialized time field */
  }
    
  /* else, src location is licit: initialize properly */
  /* locate node closest to source and associated cell */
  xs=NINT(fxs);
  ys=NINT(fys);
  zs=NINT(fzs);
  if(xs==fxs && ys==fys && zs==fzs){
    source_at_node=1;
    if(SMALLTALK) message("\nPoint source at node [%d,%d,%d].",
                          xs,ys,zs);
    xsc=(xs==nmesh_x)? xs-1:xs;
    ysc=(ys==nmesh_y)? ys-1:ys;
    zsc=(zs==nmesh_z)? zs-1:zs;
  }
  else {
    /* patch 261205[1] : upper limits were not handled properly
     *      x=(fxs<xs) ? xs-1:xs;
     *      y=(fys<ys) ? ys-1:ys;
     *      z=(fzs<zs) ? zs-1:zs;
     */
    if(xs==nmesh_x) xsc=xs-1; else xsc=(fxs<xs) ? xs-1:xs;
    if(ys==nmesh_y) ysc=ys-1; else ysc=(fys<ys) ? ys-1:ys;
    if(zs==nmesh_z) zsc=zs-1; else zsc=(fzs<zs) ? zs-1:zs;
    /* end 261205[1] */
    if(SMALLTALK) message("\nPoint source at [%g,%g,%g]; \
Nearest node [%d,%d,%d].",fxs,fys,fzs,xs,ys,zs);
  }
  hs0=hs[xsc][ysc][zsc];
  timeshift=0.0; /* all times will be positive by construction */
    
  /* initialize inclusive boundaries of explored region */
  /* patch 070206[1] : many changes in the way the homogeneous region
   * is searched (do-loop and its initialization)
   *  X0=MAX(xs-1,0);
   *  Y0=MAX(ys-1,0);
   *  Z0=MAX(zs-1,0);
   *  X1=MIN(xs+1,nmesh_x-1);
   *  Y1=MIN(ys+1,nmesh_y-1);
   *  Z1=MIN(zs+1,nmesh_z-1);
   */
  /* search largest parallelepipedic homogeneous box centered on the source */
  /* Note: during the do-loop, X0,X1,Y0,Y1,Z0,Z1 are the inclusive boundaries
   * of the indices of cells with quasi-constant slowness. */
  X0=X1=xsc;
  Y0=Y1=ysc;
  Z0=Z1=zsc;
  t_X0=t_X1=t_Y0=t_Y1=t_Z0=t_Z1=0;
  /* these flags will signal that a heterogeneity has been reached */
  if(hs_eps_init<0.0 || hs_eps_init>1.0) {
    error(ERR_EPS);
    return ERR_EPS;
  }
  allowed_delta_hs=hs0*hs_eps_init;
  /* defines tolerated inhomogeneity for exact initialization */
  do{
    int count;
    test=0;
    for(count=0;count<bx;count++) if(X0 && !t_X0){
      test++;
      x=--X0;
      for(y=Y0;y<=Y1 && !t_X0;y++)
        for(z=Z0;z<=Z1 && !t_X0;z++)
          if(fabs(hs[x][y][z]-hs0)>allowed_delta_hs) t_X0=1;
      if(t_X0) X0++;
    }
    for(count=0;count<by;count++) if(Y0 && !t_Y0){
      test++;
      y=--Y0;
      for(x=X0;x<=X1 && !t_Y0;x++)
        for(z=Z0;z<=Z1 && !t_Y0;z++)
          if(fabs(hs[x][y][z]-hs0)>allowed_delta_hs) t_Y0=1;
      if(t_Y0) Y0++;
    }
    for(count=0;count<bz;count++) if(Z0 && !t_Z0){
      test++;
      z=--Z0;
      for(x=X0;x<=X1 && !t_Z0;x++)
        for(y=Y0;y<=Y1 && !t_Z0;y++)
          if(fabs(hs[x][y][z]-hs0)>allowed_delta_hs) t_Z0=1;
      if(t_Z0) Z0++;
    }
    for(count=0;count<bx;count++) if(X1<nmesh_x-1 && !t_X1){
      test++;
      x=++X1;
      for(y=Y0;y<=Y1 && !t_X1;y++)
        for(z=Z0;z<=Z1 && !t_X1;z++)
          if(fabs(hs[x][y][z]-hs0)>allowed_delta_hs) t_X1=1;
      if(t_X1) X1--;
    }
    for(count=0;count<by;count++) if(Y1<nmesh_y-1 && !t_Y1){
      test++;
      y=++Y1;
      for(x=X0;x<=X1 && !t_Y1;x++)
        for(z=Z0;z<=Z1 && !t_Y1;z++)
          if(fabs(hs[x][y][z]-hs0)>allowed_delta_hs) t_Y1=1;
      if(t_Y1) Y1--;
    }
    for(count=0;count<bz;count++) if(Z1<nmesh_z-1 && !t_Z1){
      test++;
      z=++Z1;
      for(x=X0;x<=X1 && !t_Z1;x++)
        for(y=Y0;y<=Y1 && !t_Z1;y++)
          if(fabs(hs[x][y][z]-hs0)>allowed_delta_hs) t_Z1=1;
      if(t_Z1) Z1--;
    }
    /* patch 16102006[1]: trick to enforce quasi-cubic initialization boxes
     * to avoid severe performance degradation due to the fact that at large
     * offsets, critical conditions will probably trigger numerous recursive
     * calls during propagation of computations. */
    if(test) test=(t_X0+t_X1+t_Y0+t_Y1+t_Z0+t_Z1)? 0:1;
    /* end patch 16102006[1] */
  } while(test);
  X1++;
  Y1++;
  Z1++;
  /* X1,Y1,Z1 are incremented so that X0,X1,Y0,Y1,Z0,Z1 are now the inclusive
   * boundaries of the indices of the grid-points that will be initialized
   * "exactly" */
  /* end 070206[1] */
    
  /* patch 261205[4] : homogeneous zone shrinkage only at heterogeneous
   *                   boundaries (more straightforward...) 
   *  if(X0) X0++;
   *  if(Y0) Y0++;
   *  if(Z0) Z0++;
   *  if(X1<nmesh_x) X1--;
   *  if(Y1<nmesh_y) Y1--;
   *  if(Z1<nmesh_z) Z1--;
   */
  if(t_X0) X0++;
  if(t_Y0) Y0++;
  if(t_Z0) Z0++;
  if(t_X1) X1--;
  if(t_Y1) Y1--;
  if(t_Z1) Z1--;
  /* end 261205[4] */
  /* limits are decremented so that interfaces where heterogeneities     */
  /* were detected are dealt with by finite differences (cf. headwaves). */
  /* (but this is not necessary when we reached model boundaries !)      */
    
  /* patch 261205[5] : once shrinked, the homogeneous region may endup not
   *                   containing the source point anymore (because it is
   *                   located at the immediate vicinity of a velocity
   *                   heterogeneity). 
   *                   In such case, only minimal initialization will be 
   *                   performed (via init_nearest()).
   * The following 8 lines were added :
   */
  if(X0>fxs || X1<fxs || Y0>fys || Y1<fxs || Z0>fzs || Z1<fzs){
    X0=xsc;
    Y0=ysc;
    Z0=zsc;
    X1=xsc+1;
    Y1=ysc+1;
    Z1=zsc+1;
  }
  /* end 261205[5] */
  if( init_stage>=INIT_RECURS_LIMIT ||
      (   (X0==0 || (xs-X0)>=INIT_MIN) &&
          (Y0==0 || (ys-Y0)>=INIT_MIN) &&
          (Z0==0 || (zs-Z0)>=INIT_MIN) &&
          (X1==nmesh_x || (X1-xs)>=INIT_MIN) &&
          (Y1==nmesh_y || (Y1-ys)>=INIT_MIN) &&
          (Z1==nmesh_z || (Z1-zs)>=INIT_MIN)     )   ) {
    /* patch 261205[6] : this was simply wrong !
     *      if((X1-X0+1)*(Y1-Y0+1)*(Z1-Z0+1)==1)
     */
    if((X1-X0)*(Y1-Y0)*(Z1-Z0)==1)
      /* end 261205[6] */
      init_nearest();
    else for(x=X0;x<=X1;x++)
      for(y=Y0;y<=Y1;y++)
        for(z=Z0;z<=Z1;z++){
          dist=(x-fxs)*(x-fxs)+(y-fys)*(y-fys)+(z-fzs)*(z-fzs);
          t[x][y][z]=hs0*sqrt(dist);
        }
    if(SMALLTALK)
      message("\nHomogeneous region: x[%d->%d] y[%d->%d] z[%d->%d]\n",
              X0,X1,Y0,Y1,Z0,Z1);
  } /* if smallest distance from source to boundaries of the homogeneous */
  /* box exceeds conventional limit INIT_MIN, OR if no further recursi-*/
  /* vity is allowed, then exact arrivals are computed in this region. */
    
  else {
    if((signal=recursive_init())!=ERROR_FREE) return signal;
    X0=MAX(xs-INIT_MIN,0);
    Y0=MAX(ys-INIT_MIN,0);
    Z0=MAX(zs-INIT_MIN,0);
    X1=MIN(xs+INIT_MIN,nmesh_x);
    Y1=MIN(ys+INIT_MIN,nmesh_y);
    Z1=MIN(zs+INIT_MIN,nmesh_z);
  } /* otherwise, time_3d() is used recursively   */
  /* on a re-discretized (2*INIT_MIN+1)^3 cube. */
    
  return signal;
}

/*------------------------------------------------Init_nearest()------------*/

static void init_nearest(void)

  /* initialize the 8|12|18 nearest neighbour nodes of the source    */
  /* according to source position (inside a mesh or at a boundary).  */
  /* WARNING: errors are maximal when the source is located close to */
  /* a grid-point. Best configurations are close to the centre of a  */
  /* mesh face, or of a mesh. Errors increase (anisotropically) when */
  /* the source gets close to a grid-point. Better use the grid-     */
  /* point itself as the source in such case...                      */
{
  int   x,y,z;
  double distx,disty,distz;

  if(source_at_node){
    /* formerly t[x][y][z]=0; now initializing the 26 nearest neighbours */
    if(xs<nmesh_x && ys<nmesh_y && zs<nmesh_z)
      init_cell(0.,0.,0.,xs,ys,zs);
    if(xs && ys<nmesh_y && zs<nmesh_z) init_cell(1.,0.,0.,xs-1,ys,zs);
    if(xs<nmesh_x && ys && zs<nmesh_z) init_cell(0.,1.,0.,xs,ys-1,zs);
    if(xs<nmesh_x && ys<nmesh_y && zs) init_cell(0.,0.,1.,xs,ys,zs-1);
    if(xs<nmesh_x && ys && zs)         init_cell(0.,1.,1.,xs,ys-1,zs-1);
    if(xs && ys<nmesh_y && zs)         init_cell(1.,0.,1.,xs-1,ys,zs-1);
    if(xs && ys && zs<nmesh_z)         init_cell(1.,1.,0.,xs-1,ys-1,zs);
    if(xs && ys && zs)                 init_cell(1.,1.,1.,xs-1,ys-1,zs-1);
    return;
  }
  x=(fxs<xs) ? xs-1:xs;
  y=(fys<ys) ? ys-1:ys;
  z=(fzs<zs) ? zs-1:zs;
  /* x,y,z : coordinates of current cell */
  /* patch 171205[1] : fabs unnecessary, because args here are always positive
   *  distx=fabs(fxs-x);
   *  disty=fabs(fys-y);
   *  distz=fabs(fzs-z);
   */
  distx=fxs-x;
  disty=fys-y;
  distz=fzs-z;
  /* end 171205[1] */
  /* dist* : distances from source to node minx,miny,minz of current cell */

  init_cell(distx,disty,distz,x,y,z);
  /* this is enough if the source is strictly located */
  /* within the current cell (init: 8 neighbours).    */

  if(fxs==xs){
    if(fys==ys){
      if(x) init_cell(1.,0.,distz,x-1,y,z);
      if(y) init_cell(0.,1.,distz,x,y-1,z);
      if(x && y) init_cell(1.,1.,distz,x-1,y-1,z);
    }/* source located on cell edge parallel to z (18 neighbours) */
    else
      if(fzs==zs){
        if(x) init_cell(1.,disty,0.,x-1,y,z);
        if(z) init_cell(0.,disty,1.,x,y,z-1);
        if(z && x) init_cell(1.,disty,1.,x-1,y,z-1);
      }/* source located on cell edge parallel to y (18 neighbours) */
      else{
        if(x) init_cell(1.,disty,distz,x-1,y,z);
      }/* source located on cell face perpendicular to x (12 neighbours) */
  }
  else
    if(fys==ys){
      if(fzs==zs){
        if(y) init_cell(distx,1.,0.,x,y-1,z);
        if(z) init_cell(distz,0.,1.,x,y,z-1);
        if(y && z) init_cell(distx,1.,1.,x,y-1,z-1);
      }/* source located on cell edge parallel to x (18 neighbours) */
      else {
        if(y) init_cell(distx,1.,distz,x,y-1,z);
      }/* source located on cell face perpendicular to y (12 neighbours) */
    }
    else
      if(fzs==zs){
        /* patch 171205[2] : bug!
         *      if(y) init_cell(distx,disty,1.,x,y,z-1);
         */
        if(z) init_cell(distx,disty,1.,x,y,z-1);
        /* end 171205[2] */
      }/* source located on cell face perpendicular to z (12 neighbours) */

}
/*------------------------------------------------Init_cell()---------------*/

static void init_cell(double vx, double vy, double vz, int xl, int yl, int zl)

  /* compute delays between doubleing source and nodes of current cell     */
  /* xl,yl,zl are current cell coordinates,                               */
  /* vx,vy,vz are distances from source to node xl,yl,zl (0<=vx<=1.0,...) */

{
  double est,hsc,hsn,vh,vv,hs3;
  est=exact_delay(vx,vy,vz,xl,yl,zl);
  if(est<t[xl][yl][zl]) t[xl][yl][zl]=est;
  est=exact_delay(1.0-vx,vy,vz,xl,yl,zl);
  if(est<t[xl+1][yl][zl]) t[xl+1][yl][zl]=est;
  est=exact_delay(vx,1.0-vy,vz,xl,yl,zl);
  if(est<t[xl][yl+1][zl]) t[xl][yl+1][zl]=est;
  est=exact_delay(vx,vy,1.0-vz,xl,yl,zl);
  if(est<t[xl][yl][zl+1]) t[xl][yl][zl+1]=est;
  est=exact_delay(1.0-vx,1.0-vy,vz,xl,yl,zl);
  if(est<t[xl+1][yl+1][zl]) t[xl+1][yl+1][zl]=est;
  est=exact_delay(1.0-vx,vy,1.0-vz,xl,yl,zl);
  if(est<t[xl+1][yl][zl+1]) t[xl+1][yl][zl+1]=est;
  est=exact_delay(vx,1.0-vy,1.0-vz,xl,yl,zl);
  if(est<t[xl][yl+1][zl+1]) t[xl][yl+1][zl+1]=est;
  est=exact_delay(1.0-vx,1.0-vy,1.0-vz,xl,yl,zl);
  if(est<t[xl+1][yl+1][zl+1]) t[xl+1][yl+1][zl+1]=est;

  /* patch 271205[1]: manage headwaves at cell boundaries (edges & vertices) */
  hsc=hs[xl][yl][zl];
  if(xl && (hsn=hs[xl-1][yl][zl])<hsc){
    vh=sqrt(vy*vy+vz*vz);
    if((est=init_cellh(vh,vx,hsc,hsn))<t[xl][yl][zl])
      t[xl][yl][zl]=est;
    vh=sqrt((1.0-vy)*(1.0-vy)+vz*vz);
    if((est=init_cellh(vh,vx,hsc,hsn))<t[xl][yl+1][zl])
      t[xl][yl+1][zl]=est;
    vh=sqrt(vy*vy+(1.0-vz)*(1.0-vz));
    if((est=init_cellh(vh,vx,hsc,hsn))<t[xl][yl][zl+1])
      t[xl][yl][zl+1]=est;
    vh=sqrt((1.0-vy)*(1.0-vy)+(1.0-vz)*(1.0-vz));
    if((est=init_cellh(vh,vx,hsc,hsn))<t[xl][yl+1][zl+1])
      t[xl][yl+1][zl+1]=est;
  }
  if(xl<nmesh_x-1 && (hsn=hs[xl+1][yl][zl])<hsc){
    vh=sqrt(vy*vy+vz*vz);
    if((est=init_cellh(vh,1.0-vx,hsc,hsn))<t[xl+1][yl][zl])
      t[xl+1][yl][zl]=est;
    vh=sqrt((1.0-vy)*(1.0-vy)+vz*vz);
    if((est=init_cellh(vh,1.0-vx,hsc,hsn))<t[xl+1][yl+1][zl])
      t[xl+1][yl+1][zl]=est;
    vh=sqrt(vy*vy+(1.0-vz)*(1.0-vz));
    if((est=init_cellh(vh,1.0-vx,hsc,hsn))<t[xl+1][yl][zl+1])
      t[xl+1][yl][zl+1]=est;
    vh=sqrt((1.0-vy)*(1.0-vy)+(1.0-vz)*(1.0-vz));
    if((est=init_cellh(vh,1.0-vx,hsc,hsn))<t[xl+1][yl+1][zl+1])
      t[xl+1][yl+1][zl+1]=est;
  }
  if(yl && (hsn=hs[xl][yl-1][zl])<hsc){
    vh=sqrt(vx*vx+vz*vz);
    if((est=init_cellh(vh,vy,hsc,hsn))<t[xl][yl][zl])
      t[xl][yl][zl]=est;
    vh=sqrt((1.0-vx)*(1.0-vx)+vz*vz);
    if((est=init_cellh(vh,vy,hsc,hsn))<t[xl+1][yl][zl])
      t[xl+1][yl][zl]=est;
    vh=sqrt(vx*vx+(1.0-vz)*(1.0-vz));
    if((est=init_cellh(vh,vy,hsc,hsn))<t[xl][yl][zl+1])
      t[xl][yl][zl+1]=est;
    vh=sqrt((1.0-vx)*(1.0-vx)+(1.0-vz)*(1.0-vz));
    if((est=init_cellh(vh,vy,hsc,hsn))<t[xl+1][yl][zl+1])
      t[xl+1][yl][zl+1]=est;
  }
  if(yl<nmesh_y-1 && (hsn=hs[xl][yl+1][zl])<hsc){
    vh=sqrt(vx*vx+vz*vz);
    if((est=init_cellh(vh,1.0-vy,hsc,hsn))<t[xl][yl+1][zl])
      t[xl][yl+1][zl]=est;
    vh=sqrt((1.0-vx)*(1.0-vx)+vz*vz);
    if((est=init_cellh(vh,1.0-vy,hsc,hsn))<t[xl+1][yl+1][zl])
      t[xl+1][yl+1][zl]=est;
    vh=sqrt(vx*vx+(1.0-vz)*(1.0-vz));
    if((est=init_cellh(vh,1.0-vy,hsc,hsn))<t[xl][yl+1][zl+1])
      t[xl][yl+1][zl+1]=est;
    vh=sqrt((1.0-vx)*(1.0-vx)+(1.0-vz)*(1.0-vz));
    if((est=init_cellh(vh,1.0-vy,hsc,hsn))<t[xl+1][yl+1][zl+1])
      t[xl+1][yl+1][zl+1]=est;
  }
  if(zl && (hsn=hs[xl][yl][zl-1])<hsc){
    vh=sqrt(vx*vx+vy*vy);
    if((est=init_cellh(vh,vz,hsc,hsn))<t[xl][yl][zl])
      t[xl][yl][zl]=est;
    vh=sqrt((1.0-vx)*(1.0-vx)+vy*vy);
    if((est=init_cellh(vh,vz,hsc,hsn))<t[xl+1][yl][zl])
      t[xl+1][yl][zl]=est;
    vh=sqrt(vx*vx+(1.0-vy)*(1.0-vy));
    if((est=init_cellh(vh,vz,hsc,hsn))<t[xl][yl+1][zl])
      t[xl][yl+1][zl]=est;
    vh=sqrt((1.0-vx)*(1.0-vx)+(1.0-vy)*(1.0-vy));
    if((est=init_cellh(vh,vz,hsc,hsn))<t[xl+1][yl+1][zl])
      t[xl+1][yl+1][zl]=est;
  }
  if(zl<nmesh_z-1 && (hsn=hs[xl][yl][zl+1])<hsc){
    vh=sqrt(vx*vx+vy*vy);
    if((est=init_cellh(vh,1.0-vz,hsc,hsn))<t[xl][yl][zl+1])
      t[xl][yl][zl+1]=est;
    vh=sqrt((1.0-vx)*(1.0-vx)+vy*vy);
    if((est=init_cellh(vh,1.0-vz,hsc,hsn))<t[xl+1][yl][zl+1])
      t[xl+1][yl][zl+1]=est;
    vh=sqrt(vx*vx+(1.0-vy)*(1.0-vy));
    if((est=init_cellh(vh,1.0-vz,hsc,hsn))<t[xl][yl+1][zl+1])
      t[xl][yl+1][zl+1]=est;
    vh=sqrt((1.0-vx)*(1.0-vx)+(1.0-vy)*(1.0-vy));
    if((est=init_cellh(vh,1.0-vz,hsc,hsn))<t[xl+1][yl+1][zl+1])
      t[xl+1][yl+1][zl+1]=est;
  }
  if(xl && yl){
    hs3=MIN3(hs[xl-1][yl][zl],hs[xl][yl-1][zl],hs[xl-1][yl-1][zl]);
    if(hs3<hsc){
      vv=sqrt(vx*vx+vy*vy);
      if((est=init_cellh(vz,vv,hsc,hs3))<t[xl][yl][zl])
        t[xl][yl][zl]=est;
      if((est=init_cellh(1.0-vz,vv,hsc,hs3))<t[xl][yl][zl+1])
        t[xl][yl][zl+1]=est;
    }
  }
  if(yl && zl){
    hs3=MIN3(hs[xl][yl-1][zl],hs[xl][yl][zl-1],hs[xl][yl-1][zl-1]);
    if(hs3<hsc){
      vv=sqrt(vy*vy+vz*vz);
      if((est=init_cellh(vx,vv,hsc,hs3))<t[xl][yl][zl])
        t[xl][yl][zl]=est;
      if((est=init_cellh(1.0-vx,vv,hsc,hs3))<t[xl+1][yl][zl])
        t[xl+1][yl][zl]=est;
    }
  }
  if(zl && xl){
    hs3=MIN3(hs[xl][yl][zl-1],hs[xl-1][yl][zl],hs[xl-1][yl][zl-1]);
    if(hs3<hsc){
      vv=sqrt(vz*vz+vx*vx);
      if((est=init_cellh(vy,vv,hsc,hs3))<t[xl][yl][zl])
        t[xl][yl][zl]=est;
      if((est=init_cellh(1.0-vy,vv,hsc,hs3))<t[xl][yl+1][zl])
        t[xl][yl+1][zl]=est;
    }
  }
  if(xl<nmesh_x-1 && yl<nmesh_y-1){
    hs3=MIN3(hs[xl+1][yl][zl],hs[xl][yl+1][zl],hs[xl+1][yl+1][zl]);
    if(hs3<hsc){
      vv=sqrt((1.0-vx)*(1.0-vx)+(1.0-vy)*(1.0-vy));
      if((est=init_cellh(vz,vv,hsc,hs3))<t[xl+1][yl+1][zl])
        t[xl+1][yl+1][zl]=est;
      if((est=init_cellh(1.0-vz,vv,hsc,hs3))<t[xl+1][yl+1][zl+1])
        t[xl+1][yl+1][zl+1]=est;
    }
  }
  if(yl<nmesh_y-1 && zl<nmesh_z-1){
    hs3=MIN3(hs[xl][yl+1][zl],hs[xl][yl][zl+1],hs[xl][yl+1][zl+1]);
    if(hs3<hsc){
      vv=sqrt((1.0-vy)*(1.0-vy)+(1.0-vz)*(1.0-vz));
      if((est=init_cellh(vx,vv,hsc,hs3))<t[xl][yl+1][zl+1])
        t[xl][yl+1][zl+1]=est;
      if((est=init_cellh(1.0-vx,vv,hsc,hs3))<t[xl+1][yl+1][zl+1])
        t[xl+1][yl+1][zl+1]=est;
    }
  }
  if(zl<nmesh_z-1 && xl<nmesh_x-1){
    hs3=MIN3(hs[xl][yl][zl+1],hs[xl+1][yl][zl],hs[xl+1][yl][zl+1]);
    if(hs3<hsc){
      vv=sqrt((1.0-vz)*(1.0-vz)+(1.0-vx)*(1.0-vx));
      if((est=init_cellh(vy,vv,hsc,hs3))<t[xl+1][yl][zl+1])
        t[xl+1][yl][zl+1]=est;
      if((est=init_cellh(1.0-vy,vv,hsc,hs3))<t[xl+1][yl+1][zl+1])
        t[xl+1][yl+1][zl+1]=est;
    }
  }
  if(xl && yl<nmesh_y-1){
    hs3=MIN3(hs[xl-1][yl][zl],hs[xl][yl+1][zl],hs[xl-1][yl+1][zl]);
    if(hs3<hsc){
      vv=sqrt(vx*vx+(1.0-vy)*(1.0-vy));
      if((est=init_cellh(vz,vv,hsc,hs3))<t[xl][yl+1][zl])
        t[xl][yl+1][zl]=est;
      if((est=init_cellh(1.0-vz,vv,hsc,hs3))<t[xl][yl+1][zl+1])
        t[xl][yl+1][zl+1]=est;
    }
  }
  if(yl && zl<nmesh_z-1){
    hs3=MIN3(hs[xl][yl-1][zl],hs[xl][yl][zl+1],hs[xl][yl-1][zl+1]);
    if(hs3<hsc){
      vv=sqrt(vy*vy+(1.0-vz)*(1.0-vz));
      if((est=init_cellh(vx,vv,hsc,hs3))<t[xl][yl][zl+1])
        t[xl][yl][zl+1]=est;
      if((est=init_cellh(1.0-vx,vv,hsc,hs3))<t[xl+1][yl][zl+1])
        t[xl+1][yl][zl+1]=est;
    }
  }
  if(zl && xl<nmesh_x-1){
    hs3=MIN3(hs[xl][yl][zl-1],hs[xl+1][yl][zl],hs[xl+1][yl][zl-1]);
    if(hs3<hsc){
      vv=sqrt(vz*vz+(1.0-vx)*(1.0-vx));
      if((est=init_cellh(vy,vv,hsc,hs3))<t[xl+1][yl][zl])
        t[xl+1][yl][zl]=est;
      if((est=init_cellh(1.0-vy,vv,hsc,hs3))<t[xl+1][yl+1][zl])
        t[xl+1][yl+1][zl]=est;
    }
  }
  if(xl<nmesh_x-1 && yl){
    hs3=MIN3(hs[xl+1][yl][zl],hs[xl][yl-1][zl],hs[xl+1][yl-1][zl]);
    if(hs3<hsc){
      vv=sqrt((1.0-vx)*(1.0-vx)+vy*vy);
      if((est=init_cellh(vz,vv,hsc,hs3))<t[xl+1][yl][zl])
        t[xl+1][yl][zl]=est;
      if((est=init_cellh(1.0-vz,vv,hsc,hs3))<t[xl+1][yl][zl+1])
        t[xl+1][yl][zl+1]=est;
    }
  }
  if(yl<nmesh_y-1 && zl){
    hs3=MIN3(hs[xl][yl+1][zl],hs[xl][yl][zl-1],hs[xl][yl+1][zl-1]);
    if(hs3<hsc){
      vv=sqrt((1.0-vy)*(1.0-vy)+vz*vz);
      if((est=init_cellh(vx,vv,hsc,hs3))<t[xl][yl+1][zl])
        t[xl][yl+1][zl]=est;
      if((est=init_cellh(1.0-vx,vv,hsc,hs3))<t[xl+1][yl+1][zl])
        t[xl+1][yl+1][zl]=est;
    }
  }
  if(zl<nmesh_z-1 && xl){
    hs3=MIN3(hs[xl][yl][zl+1],hs[xl-1][yl][zl],hs[xl-1][yl][zl+1]);
    if(hs3<hsc){
      vv=sqrt((1.0-vz)*(1.0-vz)+vx*vx);
      if((est=init_cellh(vy,vv,hsc,hs3))<t[xl][yl][zl+1])
        t[xl][yl][zl+1]=est;
      if((est=init_cellh(1.0-vy,vv,hsc,hs3))<t[xl][yl+1][zl+1])
        t[xl][yl+1][zl+1]=est;
    }
  }
  /* end 271205[1] */
    
}

/* patch 271205[2]: new function added */
static double init_cellh(double vh, double vv, double hsc, double hsn)

  /* hsc,hsn: hs in current,neighbour cell                                 */
  /* vv: euclidean distance from src to interface                          */
  /* vh: distance from src to target node projected onto interface         */
  /* returns headwave arrival time at target node if it exists, else INF   */
  /* The same function is used to compute 1D headwaves along vertices :    */
  /* hsn is the min value of hs in the four cells bound by the vertex      */
{
  double hsd;
  hsd=sqrt(hsc*hsc-hsn*hsn);
  if(vh*hsd>vv*hsn) return vh*hsn+vv*hsd; /* critical condition reached */
  else return T3D_INF;
}
/* end 271205[2] */

/*------------------------------------------------Recursive_init()----------*/

static int recursive_init(void)

{
  int
    signal,
    nx_,ny_,nz_,
    bx_,by_,bz_,
    xs_,ys_,zs_,
    X0_,X1_,Y0_,Y1_,Z0_,Z1_,
    n,d,
    i,ii,ihs,i0,
    j,jj,jhs,j0,
    k,kk,khs,k0;
  double
    *hs_buf_,*t_buf_,
    fxs_,fys_,fzs_,
    HS[N_INIT],T[N_INIT];

  /* increment count of recursivity level */
  init_stage++;
  if(SMALLTALK)
    message("\nRecursive initialization: level %d",init_stage);

  /* free locally allocated pointers (double ***) */
  free_ptrs(nx);

  /* save static parameters at this stage */
  nx_=nx;
  ny_=ny;
  nz_=nz;
  bx_=bx;
  by_=by;
  bz_=bz;
  hs_buf_=hs_buf;
  t_buf_=t_buf;
  xs_=xs;
  ys_=ys;
  zs_=zs;
  fxs_=fxs;
  fys_=fys;
  fzs_=fzs;
  X0_=X0;
  X1_=X1;
  Y0_=Y0;
  Y1_=Y1;
  Z0_=Z0;
  Z1_=Z1;

  /* build the re-discretized local model and the associated source position */
  for(i=0;i<N_INIT;i++) HS[i]=T[i]=T3D_INF;
  nx=ny=nz=N_INIT_X;
  xs=ys=zs=2*INIT_MIN+1;
  i0=j0=k0=1;
  ihs=xs_-INIT_MIN-1;
  if((d=INIT_MIN-xs_)>=0){
    ihs+=d+1;
    d=1+2*d;
    nx-=d;
    xs-=d;
    i0=0;
  }
  if((d=xs_+INIT_MIN-nx_+1)>=0) nx-=1+2*d;
  jhs=ys_-INIT_MIN-1;
  if((d=INIT_MIN-ys_)>=0){
    jhs+=d+1;
    d=1+2*d;
    ny-=d;
    ys-=d;
    j0=0;
  }
  if((d=ys_+INIT_MIN-ny_+1)>=0) ny-=1+2*d;
  khs=zs_-INIT_MIN-1;
  if((d=INIT_MIN-zs_)>=0){
    khs+=d+1;
    d=1+2*d;
    nz-=d;
    zs-=d;
    k0=0;
  }
  if((d=zs_+INIT_MIN-nz_+1)>=0) nz-=1+2*d;
  for(i=ihs,n=ii=0;ii<nx;ii++){
    for(j=jhs,jj=0;jj<ny;jj++){
      for(k=khs,kk=0;kk<nz;kk++,n++){
        HS[n]=0.5*hs_buf_[i*ny_*nz_+j*nz_+k];
        if(kk%2!=k0) k++;
      }
      if(jj%2!=j0) j++;
    }
    if(ii%2!=i0) i++;
  }/* No smoothing is associated with this re-discretization */
  fxs=xs+2.0*(fxs_-xs_);
  fys=ys+2.0*(fys_-ys_);
  fzs=zs+2.0*(fzs_-zs_);

  if(VERBOSE)
    message("\nRediscretized timefield dimensions: %d %d %d",nx,ny,nz);

  /* recursively compute times on this rediscretized model */
  signal=time_3db(HS,T,nx,ny,nz,bx,by,bz,fxs,fys,fzs,hs_eps_init,messages);

  /* assign relevant times to parent timefield */
  if(signal==ERROR_FREE){
    for(i=ihs+i0,ii=i0;ii<nx;ii+=2,i++)
      for(j=jhs+j0,jj=j0;jj<ny;jj+=2,j++)
        for(k=khs+k0,kk=k0;kk<nz;kk+=2,k++)
          t_buf_[i*ny_*nz_+j*nz_+k]=T[ii*ny*nz+jj*nz+kk];
  }
  else {
    error(signal);
    signal=ERR_RECURS;
  }
 
  /* retrieve initial static parameters */
  nx=nx_;
  ny=ny_;
  nz=nz_;
  bx=bx_;
  by=by_;
  bz=bz_;
  hs_buf=hs_buf_;
  t_buf=t_buf_;
  xs=xs_;
  ys=ys_;
  zs=zs_;
  fxs=fxs_;
  fys=fys_;
  fzs=fzs_;
  X0=X0_;
  X1=X1_;
  Y0=Y0_;
  Y1=Y1_;
  Z0=Z0_;
  Z1=Z1_;
 
  /* reallocate pointers (but do not re-initialize!) */
  if((i=pre_init())!=ERROR_FREE){
    error(i);
    signal=ERR_RECURS;
  }

  /* decrement count of recursivity level */
  init_stage--;
 
  return signal;
 
}

/*------------------------------------------------Propagate_point()---------*/

static int propagate_point(int start)

{
  int
    msg,test;

  if(start!=ERROR_FREE) return start; /* Initialization failed */

  /* Make recursive_init silent */
  if(SMALLTALK) message("\nStarting F.D. computation...");
  msg=messages;
  if(init_stage) messages=0;

  /* Increment boundaries of timed zone as long as necessary... */
  /* (Outwards propagation is adopted as an initial guess).     */
  do {
    int count;
    test=0;

    for(count=0;count<bx;count++) if(X0>0){
      X0--;
      if(VERBOSE) message("\nx_side %d->%d: ",X0+1,X0);
      x_side(Y0,Y1,Z0,Z1,X0,-1);
      test++;
    }

    for(count=0;count<by;count++) if(Y0>0){
      Y0--;
      if(VERBOSE) message("\ny_side %d->%d: ",Y0+1,Y0);
      y_side(X0,X1,Z0,Z1,Y0,-1);
      test++;
    }

    for(count=0;count<bz;count++) if(Z0>0){
      Z0--;
      if(VERBOSE) message("\nz_side %d->%d: ",Z0+1,Z0);
      z_side(X0,X1,Y0,Y1,Z0,-1);
      test++;
    }

    for(count=0;count<bx;count++) if(X1<nmesh_x){
      X1++;
      if(VERBOSE) message("\nx_side %d->%d: ",X1-1,X1);
      x_side(Y0,Y1,Z0,Z1,X1,1);
      test++;
    }

    for(count=0;count<by;count++) if(Y1<nmesh_y){
      Y1++;
      if(VERBOSE) message("\ny_side %d->%d: ",Y1-1,Y1);
      y_side(X0,X1,Z0,Z1,Y1,1);
      test++;
    }

    for(count=0;count<bz;count++) if(Z1<nmesh_z){
      Z1++;
      if(VERBOSE) message("\nz_side %d->%d: ",Z1-1,Z1);
      z_side(X0,X1,Y0,Y1,Z1,1);
      test++;
    }

  } while(test);

  messages=msg;

  return ERROR_FREE;

}

/*---------------------------------------------- Free_ptrs()------------------*/

static void free_ptrs(int max_x)

{
  int x,y,z;
  double *pf;

  /* if relevant, retrieve T3D_INF-masked hs values at model boundaries */
  if(init_stage==0 && hs_keep){
    pf=hs_keep;
    for(x=0;x<nx;x++){
      for(y=0;y<ny;y++) hs[x][y][nmesh_z]= *pf++;
      for(z=0;z<nmesh_z;z++) hs[x][nmesh_y][z]= *pf++;
    }
    for(y=0;y<nmesh_y;y++)
      for(z=0;z<nmesh_z;z++) hs[nmesh_x][y][z]= *pf++;
    free((char *)hs_keep);
  }

  /* if relevant, undo timeshift (see init_point() for more comments) */
  if(init_stage==0 && timeshift<0.0){
    for(x=0;x<nx;x++)
      for(y=0;y<ny;y++)
        for(z=0;z<nz;z++)
          t[x][y][z]-=timeshift;
  }

  /* free pointers */
  for(x=0;x<max_x;x++){
    free((char *)hs[x]);
    free((char *)t[x]);
  }
  free((char *)hs);
  free((char *)t);
  free((char *)longflags);

}
/****end mail1/3****/
/*--------------------LOCAL 3-D STENCILS (FUNCTIONS AND MACROS)---------------*/
/* Counts of stencils refer to the local inwards isotropic formulation of the */
/* local finite difference computation function (170 stencils for a given     */
/* point, taken as the center of a 2*2*2 cube, with 26 timed neighbours).     */
/* See Podvin and Lecomte, 1991.                                              */
/* In this implementation, this function is never fully computed, because we  */
/* sequentially select only a limited number of relevant directions of propa- */
/* gation, according to the recent history of the signal.                     */
/* As a consequence, only 28 stencils are generally tested at each point.     */
/* The corresponding count is indicated between brackets.                     */
/*----------------------------------------------------------------------------*/
/* Code was restructured in order to minimize the number of actually computed */
/* sqrt()s. Some redundancy was introduced in tests in order to cope properly */
/* with problems arising from numerical errors (e.g. sqrt(x*x) may be lower   */
/* than x, a fact leading to infinite recursions in some awkward cases).      */
/*----------------------------------------------------------------------------*/

/*----------------------------------------- exact_delay() ------------------- */

static double exact_delay(double vx, double vy, double vz, int xm, int ym, int zm)

{
  double estimate;

  if(xm<0 || xm>=nmesh_x || ym<0 || ym>=nmesh_y || zm<0 || zm>=nmesh_z)
    return T3D_INF;
  estimate=(vx*vx+vy*vy+vz*vz)*hs[xm][ym][zm]*hs[xm][ym][zm];
  return sqrt(estimate);
}

/*----------------------------------------- 1-D transmission : 6 stencils [3] */
/*------------------------------------- (Direct arrival from first neighbour) */

static int t_1d(int x, int y, int z,
                double t0, double hs0, double hs1, double hs2, double hs3)

{
  double estimate,dt;
  estimate=t0+MIN4(hs0,hs1,hs2,hs3);
  /* 12 July 2003: FUZZIFIED COMPARISON */
  dt=t[x][y][z]-estimate;
  if(dt>EPS_FUZZY*t[x][y][z]){
    t[x][y][z]=estimate;
    return 1;
  }
  return 0;
}

/*----------------------------------------- 2-D diffraction : 12 stencils [3] */
/*------------------------------------ (Direct arrival from second neighbour) */

static int diff_2d(int x, int y, int z, double t0, double hs0, double hs1)

{
  double estimate,dt;
  estimate=t0+M_SQRT2*MIN(hs0,hs1);
  /* 12 July 2003: FUZZIFIED COMPARISON */
  dt=t[x][y][z]-estimate;
  if(dt>EPS_FUZZY*t[x][y][z]){
    t[x][y][z]=estimate;
    return 1;
  }
  return 0;
}

/*------------------------------------ 3-D point diffraction : 8 stencils [1] */
/*------------------------------------- (Direct arrival from third neighbour) */

static int point_diff(int x, int y, int z, double t0, double hs0)

{
  double estimate;
  estimate=t0+hs0*M_SQRT3;
  /* not "fuzzified" as it cannot trigger recursive calls (not a headwave) */
  if(estimate<t[x][y][z]){
    t[x][y][z]=estimate;
    return 1;
  }
  return 0;
}

/*---------------------------------------- 2-D transmission : 24 stencils [6] */
/*----------------------------------------- (Arrival from coplanar mesh edge) */

static int t_2d(int x, int y, int z, double t0, double t1, double hs0, double hs1)

{   double estimate,dt,hsm,test2,u2;
  dt=t1-t0;
  test2=t[x][y][z]-t1;
  if(dt<0.0 || test2<0.0) return 0;
  test2*=test2;
  hsm=MIN(hs0,hs1);
  u2=hsm*hsm-dt*dt;
  if(dt<=hsm/M_SQRT2 && u2<=test2){
    estimate=t1+sqrt(u2);
    /* 12 July 2003: FUZZIFIED COMPARISON */
    dt=t[x][y][z]-estimate;
    if(dt>EPS_FUZZY*t[x][y][z]){
      t[x][y][z]=estimate;
      return 1;
    }
  }
  return 0;
}

/*------------------------------------ 3-D edge diffraction : 24 stencils [3] */
/*------------------------------------- (Arrival from non-coplanar mesh edge) */

static int edge_diff(int x, int y, int z, double t0, double t1, double hs0)

{   double estimate,u2,test2,dt;
  dt=t1-t0;
  test2=t[x][y][z]-t1;
  if(dt<0.0 || test2<0.0) return 0;
  test2*=test2;
  u2=hs0*hs0-dt*dt;
  if(dt<=hs0/M_SQRT3 && 2.0*u2<=test2){
    estimate=t1+M_SQRT2*sqrt(u2);
    /* not "fuzzified" as it cannot trigger recursive calls (not a headwave) */
    if(estimate<t[x][y][z]){
      t[x][y][z]=estimate;
      return 1;
    }
  }
  return 0;
}

/*--------------------------------------- 3-D transmission : 96 stencils [12] */
/*------------------------------------- (Arrival from non-coplanar interface) */
/* 4 stencils per function call or 1+3 using two function calls. */

/*! \cond */
#define t_3d(x,y,z,a,b,c,d,e)        t_3d_(x,y,z,a,b,c,d,e,0)
#define t_3d_part2(x,y,z,a,b,c,d,e)  t_3d_(x,y,z,a,b,c,d,e,1)
/*! \endcond */

static int t_3d_(int x, int y, int z, double t0, double tl, double tr, double td,
                 double hs0, int redundant)

  /* The current point is in diagonal position with respect to t0     */
  /* and it is a first neighbour of td. tl,tr are second neighbours.  */
  /* One of these estimators is redundant during first step of *_side */
  /* functions. See t_3d_part1() which also computes it.              */
  /* This function is always called through macros t_3d or t_3d_part2 */

  /* Note: no "fuzzification" is required because these estimates     */
  /* cannot trigger recursive calls (only due to headwave generation) */

{   double test2,r2,s2,t2,u2,dta,dtb,dta2,dtb2,estimate;
  int action;
  action=0;
  hs0*=hs0;

  dta=tl-t0;
  dtb=tr-t0;
  dta2=dta*dta;
  dtb2=dtb*dtb;
  if(dta>=0.0 && dtb>=0.0 && dta2+dtb2+dta*dtb>=0.5*hs0
     && 2.0*dta2+dtb2<=hs0 && 2.0*dtb2+dta2<=hs0){
    test2=t[x][y][z]-tr-tl+t0;
    if(test2>=0.0){
      test2*=test2;
      r2=hs0-dta2-dtb2;
      if(r2<test2){
        estimate=tr+tl-t0+sqrt(r2);
        if(estimate<t[x][y][z]){
          t[x][y][z]=estimate;
          action++;
        }
      }
    }
  }

  test2=t[x][y][z]-td;
  if(test2<0.0) return action;
  test2*=test2;
  s2=t2=u2=T3D_INF;

  dtb=td-tl;
  dtb2=dtb*dtb;
  if(dta>=0.0 && dtb>=dta && 2.0*dtb2+dta2<=hs0)
    s2=hs0-dta2-dtb2;

  dta=td-tr;
  dta2=dta*dta;
  if(!redundant && dta>=0.0 && dtb>=0.0 && dta2+dtb2+dta*dtb<=0.5*hs0)
    t2=hs0-dta2-dtb2;

  dtb=tr-t0;
  dtb2=dtb*dtb;
  if(dtb>=0.0 && dta>=dtb && 2.0*dta2+dtb2<=hs0)
    u2=hs0-dta2-dtb2;

  u2=MIN3(s2,t2,u2);
  if(u2<test2){
    estimate=td+sqrt(u2);
    if(estimate<t[x][y][z]){
      t[x][y][z]=estimate;
      action++;
    }
  }
  return action;
}

/* 3-D transmission: partial stencil, introduced because initial scheme */
/* failed to fullfill the exhaustivity condition requested by Fermat's  */
/* principle. (See "a-causal" step in *_side() functions; 18/07/91)     */
static int t_3d_part1(int x, int y, int z,
                      double t0, double tl, double tr, double hs0)

  /* The current point is a first neighbour of t0; tl,tr are two other */
  /* first neighbours of t0. Transmission through 0-l-r is tested.     */
{
  double dtl,dtr,s2,u2,estimate,test2;
  dtl=t0-tl;
  dtr=t0-tr;
  test2=t[x][y][z]-t0;
  if(test2<0.0 || dtl<0.0 || dtr<0.0) return 0;
  test2*=test2;
  hs0*=hs0;
  s2=dtl*dtl+dtr*dtr;
  if(s2+dtl*dtr>0.5*hs0) return 0;
  /* illumination condition */
  u2=hs0-s2;
  if(u2<test2){
    estimate=t0+sqrt(u2);
    if(estimate<t[x][y][z]){
      t[x][y][z]=estimate;
      return 1;
    }
  }
  return 0;
}

/*----------------------------------------------X_SIDE()--------------------*/

static int x_side(int y_begin, int y_end, int z_begin, int z_end,
                  int x, int future)

  /* Propagates computations from side x-future to current side x */
  /* between *_begin and *_end coordinates. Returns a nonzero     */
  /* integer if something actually happened (a time was lowered). */
  /* Extensions _bb, _fb etc... define simple orientation rules:  */
  /* _bf means backwards along y axis and forwards along z axis.  */
  /* So-called "longitudinal" headwaves refer to first arrivals   */
  /* due to a headwave propagating along the current side.        */

{
  int
    updated,                           /* counts adopted FD stencils */
    longhead,                          /* counts "longitudinal" headwaves */
    x0,                                /* past side coordinate */
    x_s,                               /* current meshes coordinate */
    y,z,                               /* current point coordinate */
    sign_ff,sign_bf,sign_bb,sign_fb,   /* sign flags for time differences */
    past,                              /* opposite to future ! */
    test;
  double
    hs_ff,hs_bf,hs_bb,hs_fb;           /* local slownesses */

  if(reverse_order==0)
    current_side_limit=x+future;
  updated=0;
  x0=x-future;
  if(future==1) x_s=x0;
  else          x_s=x;

  flag_fb=flag_bf=flag_ff=flag_bb=0;
  y_start_ff=y_start_fb=y_end;
  y_start_bf=y_start_bb=y_begin;
  z_start_ff=z_start_bf=z_end;
  z_start_fb=z_start_bb=z_begin;

  /* First,  Compute all stencils using only nodes of side x0.   */
  /* As only times on side x will be changed, these stencils     */
  /* are computed initially, once for all, in any order (no      */
  /* causality problem involved).                                */
  /* During this first pass, future directions of propagation    */
  /* are diagnosed, according to the time pattern on side x0.    */
  /* Borders of zones to be scanned are computed.                */
  /* This part may be seen as the explicit part of the FD scheme.*/

  for(y=y_begin;y<=y_end;y++){
    for(z=z_begin;z<=z_end;z++){

      hs_ff=hs[x_s][y][z];
      if(y>0) hs_bf=hs[x_s][y-1][z];
      else hs_bf=T3D_INF;
      if(z>0 && y>0) hs_bb=hs[x_s][y-1][z-1];
      else hs_bb=T3D_INF;
      if(z>0) hs_fb=hs[x_s][y][z-1];
      else hs_fb=T3D_INF;
      sign_fb=sign_bf=sign_ff=sign_bb=0;

      /* illuminate first neighbours */
      /* 1 1D transmission and 4 partial 3D transmission */
      updated+=t_1d(x,y,z,t[x0][y][z],hs_ff,hs_bf,hs_bb,hs_fb);
      if(y<y_end && z<z_end)
        updated+=t_3d_part1(x,y,z,
                            t[x0][y][z],t[x0][y+1][z],t[x0][y][z+1],hs_ff);
      if(y>y_begin && z<z_end)
        updated+=t_3d_part1(x,y,z,
                            t[x0][y][z],t[x0][y-1][z],t[x0][y][z+1],hs_bf);
      if(y>y_begin && z>z_begin)
        updated+=t_3d_part1(x,y,z,
                            t[x0][y][z],t[x0][y-1][z],t[x0][y][z-1],hs_bb);
      if(y<y_end && z>z_begin)
        updated+=t_3d_part1(x,y,z,
                            t[x0][y][z],t[x0][y+1][z],t[x0][y][z-1],hs_fb);

      /* illuminate second neighbours (if necessary)    */
      /* 4 2D diffraction and 4 2D transmission */
      if(y<y_end && t[x0][y][z]<=t[x0][y+1][z]){
        sign_fb++;
        sign_ff++;
        if(y<y_start_ff) y_start_ff=y;
        if(y<y_start_fb) y_start_fb=y;
        updated+=diff_2d(x,y+1,z,t[x0][y][z],hs_ff,hs_fb);
        updated+=t_2d(x,y+1,z,t[x0][y][z],t[x0][y+1][z],hs_ff,hs_fb);
      }
      if(y>y_begin && t[x0][y][z]<=t[x0][y-1][z]){
        sign_bb++; 
        sign_bf++;
        if(y>y_start_bf) y_start_bf=y;
        if(y>y_start_bb) y_start_bb=y;
        updated+=diff_2d(x,y-1,z,t[x0][y][z],hs_bf,hs_bb);
        updated+=t_2d(x,y-1,z,t[x0][y][z],t[x0][y-1][z],hs_bf,hs_bb);
      }
      if(z<z_end && t[x0][y][z]<=t[x0][y][z+1]){
        sign_bf++; 
        sign_ff++;
        if(z<z_start_ff) z_start_ff=z;
        if(z<z_start_bf) z_start_bf=z;
        updated+=diff_2d(x,y,z+1,t[x0][y][z],hs_ff,hs_bf);
        updated+=t_2d(x,y,z+1,t[x0][y][z],t[x0][y][z+1],hs_ff,hs_bf);
      }
      if(z>z_begin && t[x0][y][z]<=t[x0][y][z-1]){
        sign_bb++; 
        sign_fb++;
        if(z>z_start_fb) z_start_fb=z;
        if(z>z_start_bb) z_start_bb=z;
        updated+=diff_2d(x,y,z-1,t[x0][y][z],hs_bb,hs_fb);
        updated+=t_2d(x,y,z-1,t[x0][y][z],t[x0][y][z-1],hs_bb,hs_fb);
      }

      /* illuminate third neighbours (if necessary) */
      /* 4 3D point diffraction, 8 3D edge diffraction and 12 3D transmission */
      if(sign_ff==2){
        flag_ff=1;
        updated+=point_diff(x,y+1,z+1,t[x0][y][z],hs_ff);
        updated+=edge_diff(x,y+1,z+1,t[x0][y][z],t[x0][y+1][z],hs_ff);
        updated+=edge_diff(x,y+1,z+1,t[x0][y][z],t[x0][y][z+1],hs_ff);
        updated+=t_3d_part2(x,y+1,z+1,t[x0][y][z],
                            t[x0][y+1][z],t[x0][y][z+1],t[x0][y+1][z+1],hs_ff);
      }
      if(sign_bf==2){
        flag_bf=1;
        updated+=point_diff(x,y-1,z+1,t[x0][y][z],hs_bf);
        updated+=edge_diff(x,y-1,z+1,t[x0][y][z],t[x0][y-1][z],hs_bf);
        updated+=edge_diff(x,y-1,z+1,t[x0][y][z],t[x0][y][z+1],hs_bf);
        updated+=t_3d_part2(x,y-1,z+1,t[x0][y][z],
                            t[x0][y-1][z],t[x0][y][z+1],t[x0][y-1][z+1],hs_bf);
      }
      if(sign_bb==2){
        flag_bb=1;
        updated+=point_diff(x,y-1,z-1,t[x0][y][z],hs_bb);
        updated+=edge_diff(x,y-1,z-1,t[x0][y][z],t[x0][y-1][z],hs_bb);
        updated+=edge_diff(x,y-1,z-1,t[x0][y][z],t[x0][y][z-1],hs_bb);
        updated+=t_3d_part2(x,y-1,z-1,t[x0][y][z],
                            t[x0][y-1][z],t[x0][y][z-1],t[x0][y-1][z-1],hs_bb);
      }
      if(sign_fb==2){
        flag_fb=1;
        updated+=point_diff(x,y+1,z-1,t[x0][y][z],hs_fb);
        updated+=edge_diff(x,y+1,z-1,t[x0][y][z],t[x0][y+1][z],hs_fb);
        updated+=edge_diff(x,y+1,z-1,t[x0][y][z],t[x0][y][z-1],hs_fb);
        updated+=t_3d_part2(x,y+1,z-1,t[x0][y][z],
                            t[x0][y+1][z],t[x0][y][z-1],t[x0][y+1][z-1],hs_fb);
      }
    }
  }

  /* Now, all remaining stencils depend on nodes located on the current  */
  /* side. They must be propagated causally, and occurrences of critical */
  /* conditions must be diagnosed, in order to propagate associated head */
  /* waves exhaustively. Four independent scanning directions are succes-*/
  /* sively explored, and headwave flags are used to generate the supple-*/
  /* mentary scans requested by exhaustivity ("Reverse" Propagations).   */
  /* flag_* flag  is non-zero while the corresponding direction remains  */
  /* to be examined (or reexamined). Its examination may detect critical */
  /* conditions which in their turn will make another scan necessary.    */
  /* Initialization of this process was achieved during first step.      */
  /* This second step may be seen as the implicit part of the FD scheme. */

  /* initialize local headwave flags */
  for(y=0;y<ny*nz;y++) longflags[y]=0;

  /* enforce all scans if current side has a null surface */
  /* (This may only be encountered near the source point) */
  if(y_begin==y_end || z_begin==z_end){
    flag_ff=flag_fb=flag_bf=flag_bb=1;
    y_start_ff=y_start_fb=y_begin;
    y_start_bf=y_start_bb=y_end;
    z_start_ff=z_start_bf=z_begin;
    z_start_fb=z_start_bb=z_end;
  }

  /* Reexamine each direction, while necessary */
  do{
    test=0;
    if(flag_ff){
      test++;
      if(VERBOSE) message("ff ");
      updated+=scan_x_ff(y_start_ff,y_end,z_start_ff,z_end,x0,x,x_s);
    }
    if(flag_fb){
      test++;
      if(VERBOSE) message("fb ");
      updated+=scan_x_fb(y_start_fb,y_end,z_begin,z_start_fb,x0,x,x_s);
    }
    if(flag_bb){
      test++;
      if(VERBOSE) message("bb ");
      updated+=scan_x_bb(y_begin,y_start_bb,z_begin,z_start_bb,x0,x,x_s);
    }
    if(flag_bf){
      test++;
      if(VERBOSE) message("bf ");
      updated+=scan_x_bf(y_begin,y_start_bf,z_start_bf,z_end,x0,x,x_s);
    }
  } while(test);

  /* At this stage, all points of the current side have been timed.     */
  /* Now, Reverse propagation must be invoked if a headwave propagating */
  /* along the current side was generated, because a new branch of the  */
  /* timefield (conical wave) propagates towards the already timed box. */

  for(y=longhead=0;y<ny*nz;y++) longhead+=longflags[y];

  if(longhead){

    reverse_order++;

    if(VERBOSE) message("\nReverse#%d from x_side %d",reverse_order,x);
    past= -future;
    for(x=x0;x!=current_side_limit;x+=past){
      if(x<0 || x>=nx) break;
      if(VERBOSE) message("\nupdate side x=%d: ",x);
      if(x_side(y_begin,y_end,z_begin,z_end,x,past)==0) break;
      if(VERBOSE) message("x=%d <R#%d>updated.",x,reverse_order);
    }
    if(VERBOSE) message("\nEnd Reverse#%d\n",reverse_order);

    reverse_order--;

  }

  return updated;

}

/*--------------------------------------X_SIDE() : SCAN_X_EE()--------------*/

static int
scan_x_ff(int y_start, int y_end, int z_start, int z_end,
          int x0, int x, int x_s)

  /* scan x_side by increasing y and z ("ff"=forwards, forwards)      */
  /* propagating causal stencils with provisional a-priori that some  */
  /* significant wavefronts propagate in this direction in the zone   */
  /* defined by y_start,y_end,z_start,z_end.                          */
  /* Critical conditions on any interface are detected so that other  */
  /* relevant directions of propagation due to headwave generation    */
  /* may be exhaustively taken into account at a later stage.         */

{    int
    updated=0,
    alert0,alert1,
    y,z,
    x_sf;
  double
    hs_bf,hs_bb,hs_fb,
    hs_ube,hs_ubb,hs_ueb;

  x_sf=x_s+x-x0;

  /* We first propagate headwaves along the two borders of the current zone */
  /* These headwaves are usually relevant only when a local minimum valley  */
  /* is present along these borders on the preceding x_side (x=x0).         */
  /* This is analogous with timing local minima in the 2-D implementation.  */

  /* interface waves along y_side: 1 1D transmission and 1 2D transmission */
  hs_bb=hs_ubb=hs_ube=T3D_INF;
  for(y=y_start,z=z_start;y<y_end;y++){
    hs_bf=hs[x_s][y][z];
    if(z) hs_bb=hs[x_s][y][z-1];
    if(x_sf>=0 && x_sf<nmesh_x){
      hs_ube=hs[x_sf][y][z];
      if(z) hs_ubb=hs[x_sf][y][z-1];
    }
    alert1=t_1d(x,y+1,z,t[x][y][z],hs_bb,hs_bf,hs_ubb,hs_ube);
    alert0=t_2d(x,y+1,z,t[x0][y][z],t[x][y][z],hs_bb,hs_bf);
    updated+=alert0+alert1;
    if(alert1) longflags[y*nz+nz+z]=1;
    if(alert0) longflags[y*nz+nz+z]=0;
  }

  /* interface waves along z_side: 1 1D transmission and 1 2D transmission */
  hs_bb=hs_ubb=hs_ueb=T3D_INF;
  for(y=y_start,z=z_start;z<z_end;z++){
    hs_fb=hs[x_s][y][z];
    if(y) hs_bb=hs[x_s][y-1][z];
    if(x_sf>=0 && x_sf<nmesh_x){
      hs_ueb=hs[x_sf][y][z];
      if(y) hs_ubb=hs[x_sf][y-1][z];
    }
    alert1=t_1d(x,y,z+1,t[x][y][z],hs_bb,hs_fb,hs_ubb,hs_ueb);
    alert0=t_2d(x,y,z+1,t[x0][y][z],t[x][y][z],hs_bb,hs_fb);
    updated+=alert0+alert1;
    if(alert1) longflags[y*nz+z+1]=1;
    if(alert0) longflags[y*nz+z+1]=0;
  }

  /* We now propagate bulk and head waves into the region of interest. */

  for(y=y_start;y<y_end;y++){
    for(z=z_start;z<z_end;z++){

      hs_bb=hs[x_s][y][z];
      hs_bf=hs[x_s][y][z+1];
      hs_fb=hs[x_s][y+1][z];
      if(x_sf>=0 && x_sf<nmesh_x){
        hs_ubb=hs[x_sf][y][z];
        hs_ube=hs[x_sf][y][z+1];
        hs_ueb=hs[x_sf][y+1][z];
      }
      else hs_ubb=hs_ube=hs_ueb=T3D_INF;

      /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
      alert0=edge_diff(x,y+1,z+1,t[x0][y][z],t[x][y][z],hs_bb)
        +t_3d(x,y+1,z+1,t[x0][y][z],
              t[x0][y+1][z],t[x][y][z],t[x][y+1][z],hs_bb)
        +t_3d(x,y+1,z+1,t[x0][y][z],
              t[x0][y][z+1],t[x][y][z],t[x][y][z+1],hs_bb);
      if(alert0){
        updated+=alert0;
        longflags[y*nz+nz+z+1]=0;
      }

      /* interface waves along y_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x,y+1,z+1,t[x][y][z+1],hs_bb,hs_bf,hs_ubb,hs_ube);
      alert0=t_2d(x,y+1,z+1,t[x0][y][z+1],t[x][y][z+1],hs_bb,hs_bf);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_fb++; /* this scan must be (re-)examined */
        if(y_start_fb>y) y_start_fb=y;
        if(z_start_fb<z+1) z_start_fb=z+1;
        if(alert1) longflags[y*nz+nz+z+1]=1;
        else       longflags[y*nz+nz+z+1]=0;
      }

      /* interface waves along z_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x,y+1,z+1,t[x][y+1][z],hs_bb,hs_fb,hs_ubb,hs_ueb);
      alert0=t_2d(x,y+1,z+1,t[x0][y+1][z],t[x][y+1][z],hs_bb,hs_fb);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_bf++; /* this scan must be (re-)examined */
        if(y_start_bf<y+1) y_start_bf=y+1;
        if(z_start_bf>z) z_start_bf=z;
        if(alert1) longflags[y*nz+nz+z+1]=1;
        else       longflags[y*nz+nz+z+1]=0;
      }

      /* interface waves along x_side : 2 2D transmission and 1 2D diffraction */
      alert1=diff_2d(x,y+1,z+1,t[x][y][z],hs_bb,hs_ubb)
        +t_2d(x,y+1,z+1,t[x][y][z],t[x][y+1][z],hs_bb,hs_ubb)
        +t_2d(x,y+1,z+1,t[x][y][z],t[x][y][z+1],hs_bb,hs_ubb);
      if(alert1){
        updated+=alert1;
        longflags[y*nz+nz+z+1]=1;
      }
    }
  }

  flag_ff=0;
  y_start_ff=y_end;
  z_start_ff=z_end;
  /* this direction has been examined: unset corresponding flag */

  return updated;
}

/*--------------------------------------X_SIDE() : SCAN_X_BE()--------------*/

static int
scan_x_bf(int y_begin, int y_start, int z_start, int z_end,
          int x0, int x, int x_s)

{    int
    updated=0,
    alert0,alert1,
    y,z,
    x_sf;
  double
    hs_ff,hs_bb,hs_fb,
    hs_uee,hs_ubb,hs_ueb;

  x_sf=x_s+x-x0;

  hs_fb=hs_uee=hs_ueb=T3D_INF;
  for(y=y_start,z=z_start;y>y_begin;y--){
    hs_ff=hs[x_s][y-1][z];
    if(z) hs_fb=hs[x_s][y-1][z-1];
    if(x_sf>=0 && x_sf<nmesh_x){
      hs_uee=hs[x_sf][y-1][z];
      if(z) hs_ueb=hs[x_sf][y-1][z-1];
    }
    alert1=t_1d(x,y-1,z,t[x][y][z],hs_fb,hs_ff,hs_ueb,hs_uee);
    alert0=t_2d(x,y-1,z,t[x0][y][z],t[x][y][z],hs_fb,hs_ff);
    updated+=alert0+alert1;
    if(alert1) longflags[y*nz-nz+z]=1;
    if(alert0) longflags[y*nz-nz+z]=0;
  }

  hs_bb=hs_ubb=hs_ueb=T3D_INF;
  for(y=y_start,z=z_start;z<z_end;z++){
    if(y) hs_bb=hs[x_s][y-1][z];
    hs_fb=hs[x_s][y][z];
    if(x_sf>=0 && x_sf<nmesh_x){
      if(y) hs_ubb=hs[x_sf][y-1][z];
      hs_ueb=hs[x_sf][y][z];
    }
    alert1=t_1d(x,y,z+1,t[x][y][z],hs_bb,hs_fb,hs_ubb,hs_ueb);
    alert0=t_2d(x,y,z+1,t[x0][y][z],t[x][y][z],hs_bb,hs_fb);
    updated+=alert0+alert1;
    if(alert1) longflags[y*nz+z+1]=1;
    if(alert0) longflags[y*nz+z+1]=0;
  }

  for(y=y_start;y>y_begin;y--){
    for(z=z_start;z<z_end;z++){

      hs_ff=hs[x_s][y-1][z+1];
      if(y>1) hs_bb=hs[x_s][y-2][z];
      else hs_bb=T3D_INF;
      hs_fb=hs[x_s][y-1][z];
      if(x_sf>=0 && x_sf<nmesh_x){
        hs_uee=hs[x_sf][y-1][z+1];
        if(y>1) hs_ubb=hs[x_sf][y-2][z];
        else hs_ubb=T3D_INF;
        hs_ueb=hs[x_sf][y-1][z];
      }
      else hs_ubb=hs_uee=hs_ueb=T3D_INF;

      /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
      alert0=edge_diff(x,y-1,z+1,t[x0][y][z],t[x][y][z],hs_fb)
        +t_3d(x,y-1,z+1,t[x0][y][z],
              t[x0][y-1][z],t[x][y][z],t[x][y-1][z],hs_fb)
        +t_3d(x,y-1,z+1,t[x0][y][z],
              t[x0][y][z+1],t[x][y][z],t[x][y][z+1],hs_fb);
      if(alert0){
        updated+=alert0;
        longflags[y*nz-nz+z+1]=0;
      }

      /* interface waves along y_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x,y-1,z+1,t[x][y][z+1],hs_ff,hs_fb,hs_uee,hs_ueb);
      alert0=t_2d(x,y-1,z+1,t[x0][y][z+1],t[x][y][z+1],hs_ff,hs_fb);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_bb++; /* this scan must be (re-)examined */
        if(y_start_bb<y) y_start_bb=y;
        if(z_start_bb<z+1) z_start_bb=z+1;
        if(alert1) longflags[y*nz-nz+z+1]=1;
        else       longflags[y*nz-nz+z+1]=0;
      }

      /* interface waves along z_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x,y-1,z+1,t[x][y-1][z],hs_bb,hs_fb,hs_ubb,hs_ueb);
      alert0=t_2d(x,y-1,z+1,t[x0][y-1][z],t[x][y-1][z],hs_bb,hs_fb);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_ff++; /* this scan must be (re-)examined */
        if(y_start_ff>y-1) y_start_ff=y-1;
        if(z_start_ff>z) z_start_ff=z;
        if(alert1) longflags[y*nz-nz+z+1]=1;
        else       longflags[y*nz-nz+z+1]=0;
      }

      /* interface waves along x_side : 2 2D transmission and 1 2D diffraction */
      alert1=diff_2d(x,y-1,z+1,t[x][y][z],hs_fb,hs_ueb)
        +t_2d(x,y-1,z+1,t[x][y][z],t[x][y-1][z],hs_fb,hs_ueb)
        +t_2d(x,y-1,z+1,t[x][y][z],t[x][y][z+1],hs_fb,hs_ueb);
      if(alert1){
        updated+=alert1;
        longflags[y*nz-nz+z+1]=1;
      }
    }
  }

  flag_bf=0;
  y_start_bf=y_begin;
  z_start_bf=z_end;
  /* this direction has been examined: unset corresponding flag */

  return updated;
}

/*--------------------------------------X_SIDE() : SCAN_X_BB()--------------*/

static int scan_x_bb(int y_begin, int y_start, int z_begin, int z_start,
                     int x0, int x, int x_s)

{    int
    updated=0,
    alert0,alert1,
    y,z,
    x_sf;
  double
    hs_ff,hs_bf,hs_fb,
    hs_uee,hs_ube,hs_ueb;

  x_sf=x_s+x-x0;

  hs_ff=hs_uee=hs_ueb=T3D_INF;
  for(y=y_start,z=z_start;y>y_begin;y--){
    if(z) hs_ff=hs[x_s][y-1][z-1];
    hs_fb=hs[x_s][y-1][z];
    if(x_sf>=0 && x_sf<nmesh_x){
      if(z) hs_uee=hs[x_sf][y-1][z-1];
      hs_ueb=hs[x_sf][y-1][z];
    }
    alert1=t_1d(x,y-1,z,t[x][y][z],hs_fb,hs_ff,hs_ueb,hs_uee);
    alert0=t_2d(x,y-1,z,t[x0][y][z],t[x][y][z],hs_fb,hs_ff);
    updated+=alert0+alert1;
    if(alert1) longflags[y*nz-nz+z]=1;
    if(alert0) longflags[y*nz-nz+z]=0;
  }

  hs_ff=hs_uee=hs_ube=T3D_INF;
  for(y=y_start,z=z_start;z>z_begin;z--){
    if(y) hs_ff=hs[x_s][y-1][z-1];
    hs_bf=hs[x_s][y][z-1];
    if(x_sf>=0 && x_sf<nmesh_x){
      if(y) hs_uee=hs[x_sf][y-1][z-1];
      hs_ube=hs[x_sf][y][z-1];
    }
    alert1=t_1d(x,y,z-1,t[x][y][z],hs_bf,hs_ff,hs_ube,hs_uee);
    alert0=t_2d(x,y,z-1,t[x0][y][z],t[x][y][z],hs_bf,hs_ff);
    updated+=alert0+alert1;
    if(alert1) longflags[y*nz+z-1]=1;
    if(alert0) longflags[y*nz+z-1]=0;
  }

  for(y=y_start;y>y_begin;y--){
    for(z=z_start;z>z_begin;z--){

      hs_ff=hs[x_s][y-1][z-1];
      if(y>1) hs_bf=hs[x_s][y-2][z-1];
      else hs_bf=T3D_INF;
      if(z>1) hs_fb=hs[x_s][y-1][z-2];
      else hs_fb=T3D_INF;
      if(x_sf>=0 && x_sf<nmesh_x){
        hs_uee=hs[x_sf][y-1][z-1];
        if(y>1) hs_ube=hs[x_sf][y-2][z-1];
        else hs_ube=T3D_INF;
        if(z>1) hs_ueb=hs[x_sf][y-1][z-2];
        else hs_ueb=T3D_INF;
      }
      else hs_uee=hs_ube=hs_ueb=T3D_INF;

      /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
      alert0=edge_diff(x,y-1,z-1,t[x0][y][z],t[x][y][z],hs_ff)
        +t_3d(x,y-1,z-1,t[x0][y][z],
              t[x0][y-1][z],t[x][y][z],t[x][y-1][z],hs_ff)
        +t_3d(x,y-1,z-1,t[x0][y][z],
              t[x0][y][z-1],t[x][y][z],t[x][y][z-1],hs_ff);
      if(alert0){
        updated+=alert0;
        longflags[y*nz-nz+z-1]=0;
      }

      /* interface waves along y_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x,y-1,z-1,t[x][y][z-1],hs_ff,hs_fb,hs_uee,hs_ueb);
      alert0=t_2d(x,y-1,z-1,t[x0][y][z-1],t[x][y][z-1],hs_ff,hs_fb);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_bf++; /* this scan must be (re-)examined */
        if(y_start_bf<y) y_start_bf=y;
        if(z_start_bf>z-1) z_start_bf=z-1;
        if(alert1) longflags[y*nz-nz+z-1]=1;
        else       longflags[y*nz-nz+z-1]=0;
      }

      /* interface waves along z_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x,y-1,z-1,t[x][y-1][z],hs_ff,hs_bf,hs_uee,hs_ube);
      alert0=t_2d(x,y-1,z-1,t[x0][y-1][z],t[x][y-1][z],hs_ff,hs_bf);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_fb++; /* this scan must be (re-)examined */
        if(y_start_fb>y-1) y_start_fb=y-1;
        if(z_start_fb<z) z_start_fb=z;
        if(alert1) longflags[y*nz-nz+z-1]=1;
        else       longflags[y*nz-nz+z-1]=0;
      }

      /* interface waves along x_side : 2 2D transmission and 1 2D diffraction */
      alert1=diff_2d(x,y-1,z-1,t[x][y][z],hs_ff,hs_uee)
        +t_2d(x,y-1,z-1,t[x][y][z],t[x][y-1][z],hs_ff,hs_uee)
        +t_2d(x,y-1,z-1,t[x][y][z],t[x][y][z-1],hs_ff,hs_uee);
      if(alert1){
        updated+=alert1;
        longflags[y*nz-nz+z-1]=1;
      }
    }
  }

  flag_bb=0;
  y_start_bb=y_begin;
  z_start_bb=z_begin;
  /* this direction has been examined: unset corresponding flag */

  return updated;
}

/*--------------------------------------X_SIDE() : SCAN_X_EB()--------------*/

static int scan_x_fb(int y_start, int y_end, int z_begin, int z_start,
                     int x0, int x, int x_s)

{    int
    updated=0,
    alert0,alert1,
    y,z,
    x_sf;
  double
    hs_ff,hs_bb,hs_bf,
    hs_uee,hs_ubb,hs_ube;

  x_sf=x_s+x-x0;

  hs_bf=hs_ubb=hs_ube=T3D_INF;
  for(y=y_start,z=z_start;y<y_end;y++){
    if(z) hs_bf=hs[x_s][y][z-1];
    hs_bb=hs[x_s][y][z];
    if(x_sf>=0 && x_sf<nmesh_x){
      if(z) hs_ube=hs[x_sf][y][z-1];
      hs_ubb=hs[x_sf][y][z];
    }
    alert1=t_1d(x,y+1,z,t[x][y][z],hs_bf,hs_bb,hs_ube,hs_ubb);
    alert0=t_2d(x,y+1,z,t[x0][y][z],t[x][y][z],hs_bf,hs_bb);
    updated+=alert0+alert1;
    if(alert1) longflags[y*nz+nz+z]=1;
    if(alert0) longflags[y*nz+nz+z]=0;
  }

  hs_ff=hs_uee=hs_ube=T3D_INF;
  for(y=y_start,z=z_start;z>z_begin;z--){
    if(y) hs_ff=hs[x_s][y-1][z-1];
    hs_bf=hs[x_s][y][z-1];
    if(x_sf>=0 && x_sf<nmesh_x){
      if(y) hs_uee=hs[x_sf][y-1][z-1];
      hs_ube=hs[x_sf][y][z-1];
    }
    alert1=t_1d(x,y,z-1,t[x][y][z],hs_bf,hs_ff,hs_ube,hs_uee);
    alert0=t_2d(x,y,z-1,t[x0][y][z],t[x][y][z],hs_bf,hs_ff);
    updated+=alert0+alert1;
    if(alert1) longflags[y*nz+z-1]=1;
    if(alert0) longflags[y*nz+z-1]=0;
  }

  for(y=y_start;y<y_end;y++){
    for(z=z_start;z>z_begin;z--){

      hs_bf=hs[x_s][y][z-1];
      if(z>1) hs_bb=hs[x_s][y][z-2];
      else hs_bb=T3D_INF;
      hs_ff=hs[x_s][y+1][z-1];
      if(x_sf>=0 && x_sf<nmesh_x){
        hs_ube=hs[x_sf][y][z-1];
        if(z>1) hs_ubb=hs[x_sf][y][z-2];
        else hs_ubb=T3D_INF;
        hs_uee=hs[x_sf][y+1][z-1];
      }
      else hs_ubb=hs_ube=hs_uee=T3D_INF;

      /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
      alert0=edge_diff(x,y+1,z-1,t[x0][y][z],t[x][y][z],hs_bf)
        +t_3d(x,y+1,z-1,t[x0][y][z],
              t[x0][y+1][z],t[x][y][z],t[x][y+1][z],hs_bf)
        +t_3d(x,y+1,z-1,t[x0][y][z],
              t[x0][y][z-1],t[x][y][z],t[x][y][z-1],hs_bf);
      if(alert0){
        updated+=alert0;
        longflags[y*nz+nz+z-1]=0;
      }

      /* interface waves along y_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x,y+1,z-1,t[x][y][z-1],hs_bb,hs_bf,hs_ubb,hs_ube);
      alert0=t_2d(x,y+1,z-1,t[x0][y][z-1],t[x][y][z-1],hs_bb,hs_bf);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_ff++; /* this scan must be (re-)examined */
        if(y_start_ff>y) y_start_ff=y;
        if(z_start_ff>z-1) z_start_ff=z-1;
        if(alert1) longflags[y*nz+nz+z-1]=1;
        else       longflags[y*nz+nz+z-1]=0;
      }

      /* interface waves along z_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x,y+1,z-1,t[x][y+1][z],hs_ff,hs_bf,hs_uee,hs_ube);
      alert0=t_2d(x,y+1,z-1,t[x0][y+1][z],t[x][y+1][z],hs_ff,hs_bf);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_bb++; /* this scan must be (re-)examined */
        if(y_start_bb<y+1) y_start_bb=y+1;
        if(z_start_bb<z) z_start_bb=z;
        if(alert1) longflags[y*nz+nz+z-1]=1;
        else       longflags[y*nz+nz+z-1]=0;
      }

      /* interface waves along x_side : 2 2D transmission and 1 2D diffraction */
      alert1=diff_2d(x,y+1,z-1,t[x][y][z],hs_bf,hs_ube)
        +t_2d(x,y+1,z-1,t[x][y][z],t[x][y+1][z],hs_bf,hs_ube)
        +t_2d(x,y+1,z-1,t[x][y][z],t[x][y][z-1],hs_bf,hs_ube);
      if(alert1){
        updated++;
        longflags[y*nz+nz+z-1]=1;
      }
    }
  }

  flag_fb=0;
  y_start_fb=y_end;
  z_start_fb=z_begin;
  /* this direction has been examined: unset corresponding flag */

  return updated;
}

/****end mail2/3***-------------------------END_X_SIDE()--------------------*/
/*----------------------------------------------Y_SIDE()--------------------*/

static int y_side(int x_begin, int x_end, int z_begin, int z_end,
                  int y, int future)

  /* Propagates computations from side y-future to side y        */
  /* between *_begin and *_end coordinates. Returns a nonzero    */
  /* integer if something actually happened (a time was lowered).*/
  /* Extensions _bb, _fb etc... define simple orientation rules: */
  /* _bf means backwards along x axis and forwards along z axis. */
  /* See complete comments in function x_side().                 */

{
  int
    updated,
    longhead,
    y0,
    y_s,
    x,z,
    sign_ff,sign_bf,sign_bb,sign_fb,
    test,
    past;
  double
    hs_ff,hs_bf,hs_bb,hs_fb;

  if(reverse_order==0)
    current_side_limit=y+future;
  updated=0;
  y0=y-future;
  if(future==1) y_s=y0;
  else          y_s=y;

  flag_fb=flag_bf=flag_ff=flag_bb=0;
  x_start_ff=x_start_fb=x_end;
  x_start_bf=x_start_bb=x_begin;
  z_start_ff=z_start_bf=z_end;
  z_start_fb=z_start_bb=z_begin;

  /* First Step: "a-causal" stencils */

  for(x=x_begin;x<=x_end;x++){
    for(z=z_begin;z<=z_end;z++){
      hs_ff=hs[x][y_s][z];
      if(x>0) hs_bf=hs[x-1][y_s][z];
      else hs_bf=T3D_INF;
      if(x>0 && z>0) hs_bb=hs[x-1][y_s][z-1];
      else hs_bb=T3D_INF;
      if(z>0) hs_fb=hs[x][y_s][z-1];
      else hs_fb=T3D_INF;

      sign_ff=sign_bf=sign_bb=sign_fb=0;

      /* illuminate first neighbours */
      /* 1 1D transmission and 4 partial 3D transmission */
      updated+=t_1d(x,y,z,t[x][y0][z],hs_ff,hs_bf,hs_bb,hs_fb);
      if(x<x_end && z<z_end)
        updated+=t_3d_part1(x,y,z,
                            t[x][y0][z],t[x+1][y0][z],t[x][y0][z+1],hs_ff);
      if(x>x_begin && z<z_end)
        updated+=t_3d_part1(x,y,z,
                            t[x][y0][z],t[x-1][y0][z],t[x][y0][z+1],hs_bf);
      if(x>x_begin && z>z_begin)
        updated+=t_3d_part1(x,y,z,
                            t[x][y0][z],t[x-1][y0][z],t[x][y0][z-1],hs_bb);
      if(x<x_end && z>z_begin)
        updated+=t_3d_part1(x,y,z,
                            t[x][y0][z],t[x+1][y0][z],t[x][y0][z-1],hs_fb);

      /* illuminate second neighbours */
      /* 4 2D diffraction and 4 2D transmission */
      if(x<x_end && t[x][y0][z]<=t[x+1][y0][z]){
        sign_fb++;
        sign_ff++;
        if(x<x_start_ff) x_start_ff=x;
        if(x<x_start_fb) x_start_fb=x;
        updated+=diff_2d(x+1,y,z,t[x][y0][z],hs_ff,hs_fb);
        updated+=t_2d(x+1,y,z,t[x][y0][z],t[x+1][y0][z],hs_ff,hs_fb);
      }
      if(x>x_begin && t[x][y0][z]<=t[x-1][y0][z]){
        sign_bb++; 
        sign_bf++;
        if(x>x_start_bf) x_start_bf=x;
        if(x>x_start_bb) x_start_bb=x;
        updated+=diff_2d(x-1,y,z,t[x][y0][z],hs_bf,hs_bb);
        updated+=t_2d(x-1,y,z,t[x][y0][z],t[x-1][y0][z],hs_bf,hs_bb);
      }
      if(z<z_end && t[x][y0][z]<=t[x][y0][z+1]){
        sign_bf++; 
        sign_ff++;
        if(z<z_start_ff) z_start_ff=z;
        if(z<z_start_bf) z_start_bf=z;
        updated+=diff_2d(x,y,z+1,t[x][y0][z],hs_ff,hs_bf);
        updated+=t_2d(x,y,z+1,t[x][y0][z],t[x][y0][z+1],hs_ff,hs_bf);
      }
      if(z>z_begin && t[x][y0][z]<=t[x][y0][z-1]){
        sign_bb++; 
        sign_fb++;
        if(z>z_start_fb) z_start_fb=z;
        if(z>z_start_bb) z_start_bb=z;
        updated+=diff_2d(x,y,z-1,t[x][y0][z],hs_bb,hs_fb);
        updated+=t_2d(x,y,z-1,t[x][y0][z],t[x][y0][z-1],hs_bb,hs_fb);
      }

      /* illuminate third neighbours */
      /* 4 3D point diffraction, 8 3D edge diffraction and 12 3D transmission */
      if(sign_ff==2){
        flag_ff=1;
        updated+=point_diff(x+1,y,z+1,t[x][y0][z],hs_ff)
          +edge_diff(x+1,y,z+1,t[x][y0][z],t[x+1][y0][z],hs_ff)
          +edge_diff(x+1,y,z+1,t[x][y0][z],t[x][y0][z+1],hs_ff)
          +t_3d_part2(x+1,y,z+1,t[x][y0][z],t[x+1][y0][z],
                      t[x][y0][z+1],t[x+1][y0][z+1],hs_ff);
      }
      if(sign_bf==2){
        flag_bf=1;
        updated+=point_diff(x-1,y,z+1,t[x][y0][z],hs_bf)
          +edge_diff(x-1,y,z+1,t[x][y0][z],t[x-1][y0][z],hs_bf)
          +edge_diff(x-1,y,z+1,t[x][y0][z],t[x][y0][z+1],hs_bf)
          +t_3d_part2(x-1,y,z+1,t[x][y0][z],t[x-1][y0][z],
                      t[x][y0][z+1],t[x-1][y0][z+1],hs_bf);
      }
      if(sign_bb==2){
        flag_bb=1;
        updated+=point_diff(x-1,y,z-1,t[x][y0][z],hs_bb)
          +edge_diff(x-1,y,z-1,t[x][y0][z],t[x-1][y0][z],hs_bb)
          +edge_diff(x-1,y,z-1,t[x][y0][z],t[x][y0][z-1],hs_bb)
          +t_3d_part2(x-1,y,z-1,t[x][y0][z],t[x-1][y0][z],
                      t[x][y0][z-1],t[x-1][y0][z-1],hs_bb);
      }
      if(sign_fb==2){
        flag_fb=1;
        updated+=point_diff(x+1,y,z-1,t[x][y0][z],hs_fb)
          +edge_diff(x+1,y,z-1,t[x][y0][z],t[x+1][y0][z],hs_fb)
          +edge_diff(x+1,y,z-1,t[x][y0][z],t[x][y0][z-1],hs_fb)
          +t_3d_part2(x+1,y,z-1,t[x][y0][z],t[x+1][y0][z],
                      t[x][y0][z-1],t[x+1][y0][z-1],hs_fb);
      }
    }
  }

  /* Second Step: causal propagation */

  for(x=0;x<nx*nz;x++) longflags[x]=0;

  if(x_begin==x_end || z_begin==z_end){
    flag_ff=flag_fb=flag_bf=flag_bb=1;
    x_start_ff=x_start_fb=x_begin;
    x_start_bf=x_start_bb=x_end;
    z_start_ff=z_start_bf=z_begin;
    z_start_fb=z_start_bb=z_end;
  }

  do{
    test=0;
    if(flag_ff){
      test++;
      if(VERBOSE) message("ff ");
      updated+=scan_y_ff(x_start_ff,x_end,z_start_ff,z_end,y0,y,y_s);
    }
    if(flag_fb){
      test++;
      if(VERBOSE) message("fb ");
      updated+=scan_y_fb(x_start_fb,x_end,z_begin,z_start_fb,y0,y,y_s);
    }
    if(flag_bb){
      test++;
      if(VERBOSE) message("bb ");
      updated+=scan_y_bb(x_begin,x_start_bb,z_begin,z_start_bb,y0,y,y_s);
    }
    if(flag_bf){
      test++;
      if(VERBOSE) message("bf ");
      updated+=scan_y_bf(x_begin,x_start_bf,z_start_bf,z_end,y0,y,y_s);
    }
  } while(test);

  /* Third Step: Reverse propagation, if necessary */

  for(x=longhead=0;x<nx*nz;x++) longhead+=longflags[x];

  if(longhead){

    reverse_order++;

    if(VERBOSE) message("\nReverse#%d from y_side %d",reverse_order,y);
    past= -future;
    for(y=y0;y!=current_side_limit;y+=past){
      if(y<0 || y>=ny) break;
      if(VERBOSE) message("\nupdate side y=%d: ",y);
      if(y_side(x_begin,x_end,z_begin,z_end,y,past)==0) break;
      if(VERBOSE) message("y=%d <R#%d>updated.",y,reverse_order);
    }
    if(VERBOSE) message("\nEnd Reverse#%d\n",reverse_order);

    reverse_order--;

  }

  return updated;

}

/*--------------------------------------Y_SIDE() : SCAN_Y_EE()--------------*/

static int scan_y_ff(int x_start, int x_end, int z_start, int z_end,
                     int y0, int y, int y_s)

{    int
    updated=0,
    alert0,alert1,
    x,z,
    y_sf;
  double
    hs_bf,hs_bb,hs_fb,
    hs_ube,hs_ubb,hs_ueb;

  y_sf=y_s+y-y0;

  hs_bb=hs_ubb=hs_ube=T3D_INF;
  for(x=x_start,z=z_start;x<x_end;x++){
    hs_bf=hs[x][y_s][z];
    if(z) hs_bb=hs[x][y_s][z-1];
    if(y_sf>=0 && y_sf<nmesh_y){
      hs_ube=hs[x][y_sf][z];
      if(z) hs_ubb=hs[x][y_sf][z-1];
    }
    alert1=t_1d(x+1,y,z,t[x][y][z],hs_bb,hs_bf,hs_ubb,hs_ube);
    alert0=t_2d(x+1,y,z,t[x][y0][z],t[x][y][z],hs_bb,hs_bf);
    updated+=alert0+alert1;
    if(alert1) longflags[x*nz+nz+z]=1;
    if(alert0) longflags[x*nz+nz+z]=0;
  }
 
  hs_bb=hs_ubb=hs_ueb=T3D_INF;
  for(x=x_start,z=z_start;z<z_end;z++){
    hs_fb=hs[x][y_s][z];
    if(x) hs_bb=hs[x-1][y_s][z];
    if(y_sf>=0 && y_sf<nmesh_y){
      hs_ueb=hs[x][y_sf][z];
      if(x) hs_ubb=hs[x-1][y_sf][z];
    }
    alert1=t_1d(x,y,z+1,t[x][y][z],hs_bb,hs_fb,hs_ubb,hs_ueb);
    alert0=t_2d(x,y,z+1,t[x][y0][z],t[x][y][z],hs_bb,hs_fb);
    updated+=alert0+alert1;
    if(alert1) longflags[x*nz+z+1]=1;
    if(alert0) longflags[x*nz+z+1]=0;
  }

  for(x=x_start;x<x_end;x++){
    for(z=z_start;z<z_end;z++){
      hs_bb=hs[x][y_s][z];
      hs_bf=hs[x][y_s][z+1];
      hs_fb=hs[x+1][y_s][z];
      if(y_sf>=0 && y_sf<nmesh_y){
        hs_ubb=hs[x][y_sf][z];
        hs_ube=hs[x][y_sf][z+1];
        hs_ueb=hs[x+1][y_sf][z];
      }
      else hs_ubb=hs_ube=hs_ueb=T3D_INF;

      /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
      alert0=edge_diff(x+1,y,z+1,t[x][y0][z],t[x][y][z],hs_bb)
        +t_3d(x+1,y,z+1,t[x][y0][z],
              t[x+1][y0][z],t[x][y][z],t[x+1][y][z],hs_bb)
        +t_3d(x+1,y,z+1,t[x][y0][z],
              t[x][y0][z+1],t[x][y][z],t[x][y][z+1],hs_bb);
      if(alert0){
        updated++;
        longflags[x*nz+nz+z+1]=0;
      }

      /* interface waves along x_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x+1,y,z+1,t[x][y][z+1],hs_bb,hs_bf,hs_ubb,hs_ube);
      alert0=t_2d(x+1,y,z+1,t[x][y0][z+1],t[x][y][z+1],hs_bb,hs_bf);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_fb++; /* this scan must be (re-)examined */
        if(x_start_fb>x) x_start_fb=x;
        if(z_start_fb<z+1) z_start_fb=z+1;
        if(alert1) longflags[x*nz+nz+z+1]=1;
        else       longflags[x*nz+nz+z+1]=0;
      }

      /* interface waves along z_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x+1,y,z+1,t[x+1][y][z],hs_bb,hs_fb,hs_ubb,hs_ueb);
      alert0=t_2d(x+1,y,z+1,t[x+1][y0][z],t[x+1][y][z],hs_bb,hs_fb);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_bf++; /* this scan must be (re-)examined */
        if(x_start_bf<x+1) x_start_bf=x+1;
        if(z_start_bf>z) z_start_bf=z;
        if(alert1) longflags[x*nz+nz+z+1]=1;
        else       longflags[x*nz+nz+z+1]=0;
      }

      /* interface waves along y_side : 2 2D transmission and 1 2D diffraction */
      alert1=diff_2d(x+1,y,z+1,t[x][y][z],hs_bb,hs_ubb)
        +t_2d(x+1,y,z+1,t[x][y][z],t[x+1][y][z],hs_bb,hs_ubb)
        +t_2d(x+1,y,z+1,t[x][y][z],t[x][y][z+1],hs_bb,hs_ubb);
      if(alert1){
        updated+=alert1;
        longflags[x*nz+nz+z+1]=1;
      }
    }
  }

  flag_ff=0;
  x_start_ff=x_end;
  z_start_ff=z_end;

  return updated;
}

/*--------------------------------------Y_SIDE() : SCAN_Y_BE()--------------*/

static int scan_y_bf(int x_begin, int x_start, int z_start, int z_end,
                     int y0, int y, int y_s)

{    int
    updated=0,
    alert0,alert1,
    x,z,
    y_sf;
  double
    hs_ff,hs_bb,hs_fb,
    hs_uee,hs_ubb,hs_ueb;

  y_sf=y_s+y-y0;

  hs_fb=hs_uee=hs_ueb=T3D_INF;
  for(x=x_start,z=z_start;x>x_begin;x--){
    hs_ff=hs[x-1][y_s][z];
    if(z) hs_fb=hs[x-1][y_s][z-1];
    if(y_sf>=0 && y_sf<nmesh_y){
      hs_uee=hs[x-1][y_sf][z];
      if(z) hs_ueb=hs[x-1][y_sf][z-1];
    }   
    alert1=t_1d(x-1,y,z,t[x][y][z],hs_fb,hs_ff,hs_ueb,hs_uee);
    alert0=t_2d(x-1,y,z,t[x][y0][z],t[x][y][z],hs_fb,hs_ff);
    updated+=alert0+alert1;
    if(alert1) longflags[x*nz-nz+z]=1;
    if(alert0) longflags[x*nz-nz+z]=0;
  }
 
  hs_bb=hs_ubb=hs_ueb=T3D_INF;
  for(x=x_start,z=z_start;z<z_end;z++){
    if(x) hs_bb=hs[x-1][y_s][z];
    hs_fb=hs[x][y_s][z];
    if(y_sf>=0 && y_sf<nmesh_y){
      if(x) hs_ubb=hs[x-1][y_sf][z];
      hs_ueb=hs[x][y_sf][z];
    }
    alert1=t_1d(x,y,z+1,t[x][y][z],hs_bb,hs_fb,hs_ubb,hs_ueb);
    alert0=t_2d(x,y,z+1,t[x][y0][z],t[x][y][z],hs_bb,hs_fb);
    updated+=alert0+alert1;
    if(alert1) longflags[x*nz+z+1]=1;
    if(alert0) longflags[x*nz+z+1]=0;
  }
 
  for(x=x_start;x>x_begin;x--){
    for(z=z_start;z<z_end;z++){

      hs_ff=hs[x-1][y_s][z+1];
      if(x>1) hs_bb=hs[x-2][y_s][z];
      else hs_bb=T3D_INF;
      hs_fb=hs[x-1][y_s][z];
      if(y_sf>=0 && y_sf<nmesh_y){
        hs_uee=hs[x-1][y_sf][z+1];
        if(x>1) hs_ubb=hs[x-2][y_sf][z];
        else hs_ubb=T3D_INF;
        hs_ueb=hs[x-1][y_sf][z];
      }
      else hs_ubb=hs_uee=hs_ueb=T3D_INF;

      /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
      alert0=edge_diff(x-1,y,z+1,t[x][y0][z],t[x][y][z],hs_fb)
        +t_3d(x-1,y,z+1,t[x][y0][z],
              t[x-1][y0][z],t[x][y][z],t[x-1][y][z],hs_fb)
        +t_3d(x-1,y,z+1,t[x][y0][z],
              t[x][y0][z+1],t[x][y][z],t[x][y][z+1],hs_fb);
      if(alert0){
        updated+=alert0;
        longflags[x*nz-nz+z+1]=0;
      }

      /* interface waves along x_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x-1,y,z+1,t[x][y][z+1],hs_ff,hs_fb,hs_uee,hs_ueb);
      alert0=t_2d(x-1,y,z+1,t[x][y0][z+1],t[x][y][z+1],hs_ff,hs_fb);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_bb++; /* this scan must be (re-)examined */
        if(x_start_bb<x) x_start_bb=x;
        if(z_start_bb<z+1) z_start_bb=z+1;
        if(alert1) longflags[x*nz-nz+z+1]=1;
        else       longflags[x*nz-nz+z+1]=0;
      }

      /* interface waves along z_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x-1,y,z+1,t[x-1][y][z],hs_bb,hs_fb,hs_ubb,hs_ueb);
      alert0=t_2d(x-1,y,z+1,t[x-1][y0][z],t[x-1][y][z],hs_bb,hs_fb);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_ff++; /* this scan must be (re-)examined */
        if(x_start_ff>x-1) x_start_ff=x-1;
        if(z_start_ff>z) z_start_ff=z;
        if(alert1) longflags[x*nz-nz+z+1]=1;
        else       longflags[x*nz-nz+z+1]=0;
      }

      /* interface waves along y_side : 2 2D transmission and 1 2D diffraction */
      alert1=diff_2d(x-1,y,z+1,t[x][y][z],hs_fb,hs_ueb)
        +t_2d(x-1,y,z+1,t[x][y][z],t[x-1][y][z],hs_fb,hs_ueb)
        +t_2d(x-1,y,z+1,t[x][y][z],t[x][y][z+1],hs_fb,hs_ueb);
      if(alert1){
        updated+=alert1;
        longflags[x*nz-nz+z+1]=1;
      }
    }
  }

  flag_bf=0;
  x_start_bf=x_begin;
  z_start_bf=z_end;

  return updated;
}

/*--------------------------------------Y_SIDE() : SCAN_Y_BB()--------------*/

static int scan_y_bb(int x_begin, int x_start, int z_begin, int z_start,
                     int y0, int y, int y_s)

{    int
    updated=0,
    alert0,alert1,
    x,z,
    y_sf;
  double
    hs_ff,hs_bf,hs_fb,
    hs_uee,hs_ube,hs_ueb;

  y_sf=y_s+y-y0;

  hs_ff=hs_uee=hs_ueb=T3D_INF;
  for(x=x_start,z=z_start;x>x_begin;x--){
    if(z) hs_ff=hs[x-1][y_s][z-1];
    hs_fb=hs[x-1][y_s][z];
    if(y_sf>=0 && y_sf<nmesh_y){
      if(z) hs_uee=hs[x-1][y_sf][z-1];
      hs_ueb=hs[x-1][y_sf][z];
    }
    alert1=t_1d(x-1,y,z,t[x][y][z],hs_fb,hs_ff,hs_ueb,hs_uee);
    alert0=t_2d(x-1,y,z,t[x][y0][z],t[x][y][z],hs_fb,hs_ff);
    updated+=alert0+alert1;
    if(alert1) longflags[x*nz-nz+z]=1;
    if(alert0) longflags[x*nz-nz+z]=0;
  }   
 
  hs_ff=hs_uee=hs_ube=T3D_INF;
  for(x=x_start,z=z_start;z>z_begin;z--){
    if(x) hs_ff=hs[x-1][y_s][z-1];
    hs_bf=hs[x][y_s][z-1];
    if(y_sf>=0 && y_sf<nmesh_y){
      if(x) hs_uee=hs[x-1][y_sf][z-1];
      hs_ube=hs[x][y_sf][z-1];
    }   
    alert1=t_1d(x,y,z-1,t[x][y][z],hs_bf,hs_ff,hs_ube,hs_uee);
    alert0=t_2d(x,y,z-1,t[x][y0][z],t[x][y][z],hs_bf,hs_ff);
    updated+=alert0+alert1;
    if(alert1) longflags[x*nz+z-1]=1;
    if(alert0) longflags[x*nz+z-1]=0;
  }
 
  for(x=x_start;x>x_begin;x--){
    for(z=z_start;z>z_begin;z--){

      hs_ff=hs[x-1][y_s][z-1];
      if(x>1) hs_bf=hs[x-2][y_s][z-1];
      else hs_bf=T3D_INF;
      if(z>1) hs_fb=hs[x-1][y_s][z-2];
      else hs_fb=T3D_INF;
      if(y_sf>=0 && y_sf<nmesh_y){
        hs_uee=hs[x-1][y_sf][z-1];
        if(x>1) hs_ube=hs[x-2][y_sf][z-1];
        else hs_ube=T3D_INF;
        if(z>1) hs_ueb=hs[x-1][y_sf][z-2];
        else hs_ueb=T3D_INF;
      }
      else hs_uee=hs_ube=hs_ueb=T3D_INF;

      /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
      alert0=edge_diff(x-1,y,z-1,t[x][y0][z],t[x][y][z],hs_ff)
        +t_3d(x-1,y,z-1,t[x][y0][z],
              t[x-1][y0][z],t[x][y][z],t[x-1][y][z],hs_ff)
        +t_3d(x-1,y,z-1,t[x][y0][z],
              t[x][y0][z-1],t[x][y][z],t[x][y][z-1],hs_ff);
      if(alert0){
        updated+=alert0;
        longflags[x*nz-nz+z-1]=0;
      }

      /* interface waves along x_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x-1,y,z-1,t[x][y][z-1],hs_ff,hs_fb,hs_uee,hs_ueb);
      alert0=t_2d(x-1,y,z-1,t[x][y0][z-1],t[x][y][z-1],hs_ff,hs_fb);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_bf++; /* this scan must be (re-)examined */
        if(x_start_bf<x) x_start_bf=x;
        if(z_start_bf>z-1) z_start_bf=z-1;
        if(alert1) longflags[x*nz-nz+z-1]=1;
        else       longflags[x*nz-nz+z-1]=0;
      }

      /* interface waves along z_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x-1,y,z-1,t[x-1][y][z],hs_ff,hs_bf,hs_uee,hs_ube);
      alert0=t_2d(x-1,y,z-1,t[x-1][y0][z],t[x-1][y][z],hs_ff,hs_bf);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_fb++; /* this scan must be (re-)examined */
        if(x_start_fb>x-1) x_start_fb=x-1;
        if(z_start_fb<z) z_start_fb=z;
        if(alert1) longflags[x*nz-nz+z-1]=1;
        else       longflags[x*nz-nz+z-1]=0;
      }

      /* interface waves along y_side : 2 2D transmission and 1 2D diffraction */
      alert1=diff_2d(x-1,y,z-1,t[x][y][z],hs_ff,hs_uee)
        +t_2d(x-1,y,z-1,t[x][y][z],t[x-1][y][z],hs_ff,hs_uee)
        +t_2d(x-1,y,z-1,t[x][y][z],t[x][y][z-1],hs_ff,hs_uee);
      if(alert1){
        updated+=alert1;
        longflags[x*nz-nz+z-1]=1;
      }
    }
  }

  flag_bb=0;
  x_start_bb=x_begin;
  z_start_bb=z_begin;

  return updated;
}

/*--------------------------------------Y_SIDE() : SCAN_Y_EB()--------------*/

static int scan_y_fb(int x_start, int x_end, int z_begin, int z_start,
                     int y0, int y, int y_s)

{    int
    updated=0,
    alert0,alert1,
    x,z,
    y_sf;
  double
    hs_ff,hs_bb,hs_bf,
    hs_uee,hs_ubb,hs_ube;

  y_sf=y_s+y-y0;

  hs_bf=hs_ubb=hs_ube=T3D_INF;
  for(x=x_start,z=z_start;x<x_end;x++){
    if(z) hs_bf=hs[x][y_s][z-1];
    hs_bb=hs[x][y_s][z];
    if(y_sf>=0 && y_sf<nmesh_y){
      if(z) hs_ube=hs[x][y_sf][z-1];
      hs_ubb=hs[x][y_sf][z];
    }   
    alert1=t_1d(x+1,y,z,t[x][y][z],hs_bf,hs_bb,hs_ube,hs_ubb);
    alert0=t_2d(x+1,y,z,t[x][y0][z],t[x][y][z],hs_bf,hs_bb);
    updated+=alert0+alert1;
    if(alert1) longflags[x*nz+nz+z]=1;
    if(alert0) longflags[x*nz+nz+z]=0;
  }
 
  hs_ff=hs_uee=hs_ube=T3D_INF;
  for(x=x_start,z=z_start;z>z_begin;z--){
    if(x) hs_ff=hs[x-1][y_s][z-1];
    hs_bf=hs[x][y_s][z-1];
    if(y_sf>=0 && y_sf<nmesh_y){
      if(x) hs_uee=hs[x-1][y_sf][z-1];
      hs_ube=hs[x][y_sf][z-1];
    }
    alert1=t_1d(x,y,z-1,t[x][y][z],hs_bf,hs_ff,hs_ube,hs_uee);
    alert0=t_2d(x,y,z-1,t[x][y0][z],t[x][y][z],hs_bf,hs_ff);
    updated+=alert0+alert1;
    if(alert1) longflags[x*nz+z-1]=1;
    if(alert0) longflags[x*nz+z-1]=0;
  }
 
  for(x=x_start;x<x_end;x++){
    for(z=z_start;z>z_begin;z--){

      hs_bf=hs[x][y_s][z-1];
      if(z>1) hs_bb=hs[x][y_s][z-2];
      else hs_bb=T3D_INF;
      hs_ff=hs[x+1][y_s][z-1];
      if(y_sf>=0 && y_sf<nmesh_y){
        hs_ube=hs[x][y_sf][z-1];
        if(z>1) hs_ubb=hs[x][y_sf][z-2];
        else hs_ubb=T3D_INF;
        hs_uee=hs[x+1][y_sf][z-1];
      }
      else hs_ubb=hs_ube=hs_uee=T3D_INF;

      /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
      alert0=edge_diff(x+1,y,z-1,t[x][y0][z],t[x][y][z],hs_bf)
        +t_3d(x+1,y,z-1,t[x][y0][z],
              t[x+1][y0][z],t[x][y][z],t[x+1][y][z],hs_bf)
        +t_3d(x+1,y,z-1,t[x][y0][z],
              t[x][y0][z-1],t[x][y][z],t[x][y][z-1],hs_bf);
      if(alert0){
        updated+=alert0;
        longflags[x*nz+nz+z-1]=0;
      }

      /* interface waves along x_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x+1,y,z-1,t[x][y][z-1],hs_bb,hs_bf,hs_ubb,hs_ube);
      alert0=t_2d(x+1,y,z-1,t[x][y0][z-1],t[x][y][z-1],hs_bb,hs_bf);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_ff++; /* this scan must be (re-)examined */
        if(x_start_ff>x) x_start_ff=x;
        if(z_start_ff>z-1) z_start_ff=z-1;
        if(alert1) longflags[x*nz+nz+z-1]=1;
        else       longflags[x*nz+nz+z-1]=0;
      }

      /* interface waves along z_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x+1,y,z-1,t[x+1][y][z],hs_ff,hs_bf,hs_uee,hs_ube);
      alert0=t_2d(x+1,y,z-1,t[x+1][y0][z],t[x+1][y][z],hs_ff,hs_bf);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_bb++; /* this scan must be (re-)examined */
        if(x_start_bb<x+1) x_start_bb=x+1;
        if(z_start_bb<z) z_start_bb=z;
        if(alert1) longflags[x*nz+nz+z-1]=1;
        else       longflags[x*nz+nz+z-1]=0;
      }

      /* interface waves along y_side : 2 2D transmission and 1 2D diffraction */
      alert1=diff_2d(x+1,y,z-1,t[x][y][z],hs_bf,hs_ube)
        +t_2d(x+1,y,z-1,t[x][y][z],t[x+1][y][z],hs_bf,hs_ube)
        +t_2d(x+1,y,z-1,t[x][y][z],t[x][y][z-1],hs_bf,hs_ube);
      if(alert1){
        updated+=alert1;
        longflags[x*nz+nz+z-1]=1;
      }
    }
  }

  flag_fb=0;
  x_start_fb=x_end;
  z_start_fb=z_begin;
  /* this scan has been examined */

  return updated;
}

/*--------------------------------------------END_Y_SIDE()--------------------*/
/*----------------------------------------------Z_SIDE()--------------------*/

static int z_side(int x_begin, int x_end, int y_begin, int y_end,
                  int z, int future)

  /* Propagates computations from side z-future to side z.       */
  /* between *_begin and *_end coordinates. Returns a nonzero    */
  /* integer if something actually happened (a time was lowered).*/
  /* Extensions _bb, _fb etc... define simple orientation rules: */
  /* _bf means backwards along x axis and forwards along y axis. */
  /* See complete comments in function x_side().                 */

{
  int
    updated,
    longhead,
    z0,
    z_s,
    x,y,
    sign_ff,sign_bf,sign_bb,sign_fb,
    test,
    past;
  double
    hs_ff,hs_bf,hs_bb,hs_fb;

  if(reverse_order==0)
    current_side_limit=z+future;
  updated=0;
  z0=z-future;
  if(future==1) z_s=z0;
  else          z_s=z;

  flag_fb=flag_bf=flag_ff=flag_bb=0;
  x_start_ff=x_start_fb=x_end;
  x_start_bf=x_start_bb=x_begin;
  y_start_ff=y_start_bf=y_end;
  y_start_fb=y_start_bb=y_begin;

  /* First Step: "a-causal" stencils */

  for(x=x_begin;x<=x_end;x++){
    for(y=y_begin;y<=y_end;y++){

      hs_ff=hs[x][y][z_s];
      if(x>0) hs_bf=hs[x-1][y][z_s];
      else hs_bf=T3D_INF;
      if(x>0 && y>0) hs_bb=hs[x-1][y-1][z_s];
      else hs_bb=T3D_INF;
      if(y>0) hs_fb=hs[x][y-1][z_s];
      else hs_fb=T3D_INF;
      sign_ff=sign_bf=sign_bb=sign_fb=0;

      /* illuminate first neighbours */
      /* 1 1D transmission and 4 partial 3D transmission */
      updated+=t_1d(x,y,z,t[x][y][z0],hs_ff,hs_bf,hs_bb,hs_fb);
      if(x<x_end && y<y_end)
        updated+=t_3d_part1(x,y,z,t[x][y][z0],
                            t[x+1][y][z0],t[x][y+1][z0],hs_ff);
      if(x>x_begin && y<y_end)
        updated+=t_3d_part1(x,y,z,t[x][y][z0],
                            t[x-1][y][z0],t[x][y+1][z0],hs_bf);
      if(x>x_begin && y>y_begin)
        updated+=t_3d_part1(x,y,z,t[x][y][z0],
                            t[x-1][y][z0],t[x][y-1][z0],hs_bb);
      if(x<x_end && y>y_begin)
        updated+=t_3d_part1(x,y,z,t[x][y][z0],
                            t[x+1][y][z0],t[x][y-1][z0],hs_fb);

      /* illuminate second neighbours */
      /* 4 2D diffraction and 4 2D transmission */
      if(x<x_end && t[x][y][z0]<=t[x+1][y][z0]){
        sign_fb++;
        sign_ff++;
        if(x<x_start_ff) x_start_ff=x;
        if(x<x_start_fb) x_start_fb=x;
        updated+=diff_2d(x+1,y,z,t[x][y][z0],hs_ff,hs_fb);
        updated+=t_2d(x+1,y,z,t[x][y][z0],t[x+1][y][z0],hs_ff,hs_fb);
      }
      if(x>x_begin && t[x][y][z0]<=t[x-1][y][z0]){
        sign_bb++;
        sign_bf++;
        if(x>x_start_bf) x_start_bf=x;
        if(x>x_start_bb) x_start_bb=x;
        updated+=diff_2d(x-1,y,z,t[x][y][z0],hs_bf,hs_bb);
        updated+=t_2d(x-1,y,z,t[x][y][z0],t[x-1][y][z0],hs_bf,hs_bb);
      }
      if(y<y_end && t[x][y][z0]<=t[x][y+1][z0]){
        sign_bf++;
        sign_ff++;
        if(y<y_start_ff) y_start_ff=y;
        if(y<y_start_bf) y_start_bf=y;
        updated+=diff_2d(x,y+1,z,t[x][y][z0],hs_ff,hs_bf);
        updated+=t_2d(x,y+1,z,t[x][y][z0],t[x][y+1][z0],hs_ff,hs_bf);
      }
      if(y>y_begin && t[x][y][z0]<=t[x][y-1][z0]){
        sign_bb++;
        sign_fb++;
        if(y>y_start_fb) y_start_fb=y;
        if(y>y_start_bb) y_start_bb=y;
        updated+=diff_2d(x,y-1,z,t[x][y][z0],hs_bb,hs_fb);
        updated+=t_2d(x,y-1,z,t[x][y][z0],t[x][y-1][z0],hs_bb,hs_fb);
      }

      /* illuminate third neighbours */
      /* 4 3D point diffraction, 8 3D edge diffraction and 12 3D transmission */
      if(sign_ff==2){
        flag_ff=1;
        updated+=point_diff(x+1,y+1,z,t[x][y][z0],hs_ff)
          +edge_diff(x+1,y+1,z,t[x][y][z0],t[x+1][y][z0],hs_ff)
          +edge_diff(x+1,y+1,z,t[x][y][z0],t[x][y+1][z0],hs_ff)
          +t_3d_part2(x+1,y+1,z,t[x][y][z0],t[x+1][y][z0],
                      t[x][y+1][z0],t[x+1][y+1][z0],hs_ff);
      }
      if(sign_bf==2){
        flag_bf=1;
        updated+=point_diff(x-1,y+1,z,t[x][y][z0],hs_bf)
          +edge_diff(x-1,y+1,z,t[x][y][z0],t[x-1][y][z0],hs_bf)
          +edge_diff(x-1,y+1,z,t[x][y][z0],t[x][y+1][z0],hs_bf)
          +t_3d_part2(x-1,y+1,z,t[x][y][z0],t[x-1][y][z0],
                      t[x][y+1][z0],t[x-1][y+1][z0],hs_bf);
      }
      if(sign_bb==2){
        flag_bb=1;
        updated+=point_diff(x-1,y-1,z,t[x][y][z0],hs_bb)
          +edge_diff(x-1,y-1,z,t[x][y][z0],t[x-1][y][z0],hs_bb)
          +edge_diff(x-1,y-1,z,t[x][y][z0],t[x][y-1][z0],hs_bb)
          +t_3d_part2(x-1,y-1,z,t[x][y][z0],t[x-1][y][z0],
                      t[x][y-1][z0],t[x-1][y-1][z0],hs_bb);
      }
      if(sign_fb==2){
        flag_fb=1;
        updated+=point_diff(x+1,y-1,z,t[x][y][z0],hs_fb)
          +edge_diff(x+1,y-1,z,t[x][y][z0],t[x+1][y][z0],hs_fb)
          +edge_diff(x+1,y-1,z,t[x][y][z0],t[x][y-1][z0],hs_fb)
          +t_3d_part2(x+1,y-1,z,t[x][y][z0],t[x+1][y][z0],
                      t[x][y-1][z0],t[x+1][y-1][z0],hs_fb);
      }
    }
  }

  /* Second Step: causal propagation */

  for(x=0;x<nx*ny;x++) longflags[x]=0;

  if(x_begin==x_end || y_begin==y_end){
    flag_ff=flag_fb=flag_bf=flag_bb=1;
    x_start_ff=x_start_fb=x_begin;
    x_start_bf=x_start_bb=x_end;
    y_start_ff=y_start_bf=y_begin;
    y_start_fb=y_start_bb=y_end;
  }

  do{
    test=0;
    if(flag_ff){
      test++;
      if(VERBOSE) message("ff ");
      updated+=scan_z_ff(x_start_ff,x_end,y_start_ff,y_end,z0,z,z_s);
    }
    if(flag_fb){
      test++;
      if(VERBOSE) message("fb ");
      updated+=scan_z_fb(x_start_fb,x_end,y_begin,y_start_fb,z0,z,z_s);
    }
    if(flag_bb){
      test++;
      if(VERBOSE) message("bb ");
      updated+=scan_z_bb(x_begin,x_start_bb,y_begin,y_start_bb,z0,z,z_s);
    }
    if(flag_bf){
      test++;
      if(VERBOSE) message("bf ");
      updated+=scan_z_bf(x_begin,x_start_bf,y_start_bf,y_end,z0,z,z_s);
    }
  } while(test);

  /* Third Step: Reverse Propagation if necessary */
 
  for(x=longhead=0;x<nx*ny;x++) longhead+=longflags[x];

  if(longhead){

    reverse_order++;

    if(VERBOSE) message("\nReverse#%d from z_side %d",reverse_order,z);
    past= -future;
    for(z=z0;z!=current_side_limit;z+=past){
      if(z<0 || z>=nz) break;
      if(VERBOSE) message("\nupdate side z=%d: ",z);
      if(z_side(x_begin,x_end,y_begin,y_end,z,past)==0) break;
      if(VERBOSE) message("z=%d <R#%d>updated.",z,reverse_order);
    }
    if(VERBOSE) message("\nEnd Reverse#%d\n",reverse_order);

    reverse_order--;

  }

  return updated;

}

/*--------------------------------------Z_SIDE() : SCAN_Z_EE()--------------*/

static int scan_z_ff(int x_start, int x_end, int y_start, int y_end,
                     int z0, int z, int z_s)

{    int
    updated=0,
    alert0,alert1,
    x,y,
    z_sf;
  double
    hs_bf,hs_bb,hs_fb,hs_ube,hs_ubb,hs_ueb;

  z_sf=z_s+z-z0;

  hs_bb=hs_ubb=hs_ube=T3D_INF;
  for(x=x_start,y=y_start;x<x_end;x++){
    hs_bf=hs[x][y][z_s];
    if(y) hs_bb=hs[x][y-1][z_s];
    if(z_sf>=0 && z_sf<nmesh_z){
      hs_ube=hs[x][y][z_sf];
      if(y) hs_ubb=hs[x][y-1][z_sf];
    }   
    alert1=t_1d(x+1,y,z,t[x][y][z],hs_bb,hs_bf,hs_ubb,hs_ube);
    alert0=t_2d(x+1,y,z,t[x][y][z0],t[x][y][z],hs_bb,hs_bf);
    if(alert1) longflags[x*ny+ny+y]=1;
    if(alert0) longflags[x*ny+ny+y]=0;
  }    
 
  hs_bb=hs_ubb=hs_ueb=T3D_INF;
  for(x=x_start,y=y_start;y<y_end;y++){
    hs_fb=hs[x][y][z_s];
    if(x) hs_bb=hs[x-1][y][z_s];
    if(z_sf>=0 && z_sf<nmesh_z){
      hs_ueb=hs[x][y][z_sf];
      if(x) hs_ubb=hs[x-1][y][z_sf];
    }   
    alert1=t_1d(x,y+1,z,t[x][y][z],hs_bb,hs_fb,hs_ubb,hs_ueb);
    alert0=t_2d(x,y+1,z,t[x][y][z0],t[x][y][z],hs_bb,hs_fb);
    if(alert1) longflags[x*ny+y+1]=1;
    if(alert0) longflags[x*ny+y+1]=0;
  }    

  for(x=x_start;x<x_end;x++){
    for(y=y_start;y<y_end;y++){

      hs_bb=hs[x][y][z_s];
      hs_bf=hs[x][y+1][z_s];
      hs_fb=hs[x+1][y][z_s];
      if(z_sf>=0 && z_sf<nmesh_z){
        hs_ubb=hs[x][y][z_sf];
        hs_ube=hs[x][y+1][z_sf];
        hs_ueb=hs[x+1][y][z_sf];
      }
      else hs_ubb=hs_ube=hs_ueb=T3D_INF;

      /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
      alert0=edge_diff(x+1,y+1,z,t[x][y][z0],t[x][y][z],hs_bb)
        +t_3d(x+1,y+1,z,t[x][y][z0],
              t[x+1][y][z0],t[x][y][z],t[x+1][y][z],hs_bb)
        +t_3d(x+1,y+1,z,t[x][y][z0],
              t[x][y+1][z0],t[x][y][z],t[x][y+1][z],hs_bb);
      if(alert0){
        updated+=alert0;
        longflags[x*ny+ny+y+1]=0;
      }

      /* interface waves along x_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x+1,y+1,z,t[x][y+1][z],hs_bb,hs_bf,hs_ubb,hs_ube);
      alert0=t_2d(x+1,y+1,z,t[x][y+1][z0],t[x][y+1][z],hs_bb,hs_bf);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_fb++; /* this scan must be (re-)examined */
        if(x_start_fb>x) x_start_fb=x;
        if(y_start_fb<y+1) y_start_fb=y+1;
        if(alert1) longflags[x*ny+ny+y+1]=1;
        else       longflags[x*ny+ny+y+1]=0;
      }

      /* interface waves along y_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x+1,y+1,z,t[x+1][y][z],hs_bb,hs_fb,hs_ubb,hs_ueb);
      alert0=t_2d(x+1,y+1,z,t[x+1][y][z0],t[x+1][y][z],hs_bb,hs_fb);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_bf++; /* this scan must be (re-)examined */
        if(x_start_bf<x+1) x_start_bf=x+1;
        if(y_start_bf>y) y_start_bf=y;
        if(alert1) longflags[x*ny+ny+y+1]=1;
        else       longflags[x*ny+ny+y+1]=0;
      }

      /* interface waves along z_side : 2 2D transmission and 1 2D diffraction */
      alert1=diff_2d(x+1,y+1,z,t[x][y][z],hs_bb,hs_ubb)
        +t_2d(x+1,y+1,z,t[x][y][z],t[x+1][y][z],hs_bb,hs_ubb)
        +t_2d(x+1,y+1,z,t[x][y][z],t[x][y+1][z],hs_bb,hs_ubb);
      if(alert1){
        updated+=alert1;
        longflags[x*ny+ny+y+1]=1;
      }
    }
  }

  flag_ff=0;
  x_start_ff=x_end;
  y_start_ff=y_end;

  return updated;
}

/*--------------------------------------Z_SIDE() : SCAN_Z_BE()--------------*/

static int scan_z_bf(int x_begin, int x_start, int y_start, int y_end,
                     int z0, int z, int z_s)

{    int
    updated=0,
    alert0,alert1,
    x,y,
    z_sf;
  double
    hs_ff,hs_bb,hs_fb,hs_uee,hs_ubb,hs_ueb;

  z_sf=z_s+z-z0;

  hs_fb=hs_uee=hs_ueb=T3D_INF;
  for(x=x_start,y=y_start;x>x_begin;x--){
    hs_ff=hs[x-1][y][z_s];
    if(y) hs_fb=hs[x-1][y-1][z_s];
    if(z_sf>=0 && z_sf<nmesh_z){
      hs_uee=hs[x-1][y][z_sf];
      if(y) hs_ueb=hs[x-1][y-1][z_sf];
    }   
    alert1=t_1d(x-1,y,z,t[x][y][z],hs_fb,hs_ff,hs_ueb,hs_uee);
    alert0=t_2d(x-1,y,z,t[x][y][z0],t[x][y][z],hs_fb,hs_ff);
    if(alert1) longflags[x*ny-ny+y]=1;
    if(alert0) longflags[x*ny-ny+y]=0;
  }
 
  hs_bb=hs_ubb=hs_ueb=T3D_INF;
  for(x=x_start,y=y_start;y<y_end;y++){
    if(x) hs_bb=hs[x-1][y][z_s];
    hs_fb=hs[x][y][z_s];
    if(z_sf>=0 && z_sf<nmesh_z){
      if(x) hs_ubb=hs[x-1][y][z_sf];
      hs_ueb=hs[x][y][z_sf];
    }
    alert1=t_1d(x,y+1,z,t[x][y][z],hs_bb,hs_fb,hs_ubb,hs_ueb);
    alert0=t_2d(x,y+1,z,t[x][y][z0],t[x][y][z],hs_bb,hs_fb);
    if(alert1) longflags[x*ny+y+1]=1;
    if(alert0) longflags[x*ny+y+1]=0;
  }
 
  for(x=x_start;x>x_begin;x--){
    for(y=y_start;y<y_end;y++){

      hs_ff=hs[x-1][y+1][z_s];
      if(x>1) hs_bb=hs[x-2][y][z_s];
      else hs_bb=T3D_INF;
      hs_fb=hs[x-1][y][z_s];
      if(z_sf>=0 && z_sf<nmesh_z){
        hs_uee=hs[x-1][y+1][z_sf];
        if(x>1) hs_ubb=hs[x-2][y][z_sf];
        else hs_ubb=T3D_INF;
        hs_ueb=hs[x-1][y][z_sf];
      }
      else hs_ubb=hs_uee=hs_ueb=T3D_INF;

      /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
      alert0=edge_diff(x-1,y+1,z,t[x][y][z0],t[x][y][z],hs_fb)
        +t_3d(x-1,y+1,z,t[x][y][z0],
              t[x-1][y][z0],t[x][y][z],t[x-1][y][z],hs_fb)
        +t_3d(x-1,y+1,z,t[x][y][z0],
              t[x][y+1][z0],t[x][y][z],t[x][y+1][z],hs_fb);
      if(alert0){
        updated+=alert0;
        longflags[x*ny-ny+y+1]=0;
      }

      /* interface waves along x_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x-1,y+1,z,t[x][y+1][z],hs_ff,hs_fb,hs_uee,hs_ueb);
      alert0=t_2d(x-1,y+1,z,t[x][y+1][z0],t[x][y+1][z],hs_ff,hs_fb);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_bb++; /* this scan must be (re-)examined */
        if(x_start_bb<x) x_start_bb=x;
        if(y_start_bb<y+1) y_start_bb=y+1;
        if(alert1) longflags[x*ny-ny+y+1]=1;
        else       longflags[x*ny-ny+y+1]=0;
      }

      /* interface waves along y_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x-1,y+1,z,t[x-1][y][z],hs_bb,hs_fb,hs_ubb,hs_ueb);
      alert0=t_2d(x-1,y+1,z,t[x-1][y][z0],t[x-1][y][z],hs_bb,hs_fb);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_ff++; /* this scan must be (re-)examined */
        if(x_start_ff>x-1) x_start_ff=x-1;
        if(y_start_ff>y) y_start_ff=y;
        if(alert1) longflags[x*ny-ny+y+1]=1;
        else       longflags[x*ny-ny+y+1]=0;
      }

      /* interface waves along z_side : 2 2D transmission and 1 2D diffraction */
      alert1=diff_2d(x-1,y+1,z,t[x][y][z],hs_fb,hs_ueb)
        +t_2d(x-1,y+1,z,t[x][y][z],t[x-1][y][z],hs_fb,hs_ueb)
        +t_2d(x-1,y+1,z,t[x][y][z],t[x][y+1][z],hs_fb,hs_ueb);
      if(alert1){
        updated+=alert1;
        longflags[x*ny-ny+y+1]=1;
      }
    }
  }

  flag_bf=0;
  x_start_bf=x_begin;
  y_start_bf=y_end;

  return updated;
}

/*--------------------------------------Z_SIDE() : SCAN_Z_EB()--------------*/

static int scan_z_bb(int x_begin, int x_start, int y_begin, int y_start,
                     int z0, int z, int z_s)

{    int
    updated=0,
    alert0,alert1,
    x,y,
    z_sf;
  double
    hs_ff,hs_bf,hs_fb,hs_uee,hs_ube,hs_ueb;

  z_sf=z_s+z-z0;

  hs_ff=hs_uee=hs_ueb=T3D_INF; 
  for(x=x_start,y=y_start;x>x_begin;x--){ 
    if(y) hs_ff=hs[x-1][y-1][z_s]; 
    hs_fb=hs[x-1][y][z_s]; 
    if(z_sf>=0 && z_sf<nmesh_z){ 
      if(y) hs_uee=hs[x-1][y-1][z_sf]; 
      hs_ueb=hs[x-1][y][z_sf]; 
    } 
    alert1=t_1d(x-1,y,z,t[x][y][z],hs_fb,hs_ff,hs_ueb,hs_uee);
    alert0=t_2d(x-1,y,z,t[x][y][z0],t[x][y][z],hs_fb,hs_ff); 
    if(alert1) longflags[x*ny-ny+y]=1;
    if(alert0) longflags[x*ny-ny+y]=0;
  }      
  
  hs_ff=hs_uee=hs_ube=T3D_INF; 
  for(x=x_start,y=y_start;y>y_begin;y--){ 
    if(x) hs_ff=hs[x-1][y-1][z_s]; 
    hs_bf=hs[x][y-1][z_s]; 
    if(z_sf>=0 && z_sf<nmesh_z){ 
      if(x) hs_uee=hs[x-1][y-1][z_sf]; 
      hs_ube=hs[x][y-1][z_sf]; 
    }    
    alert1=t_1d(x,y-1,z,t[x][y][z],hs_bf,hs_ff,hs_ube,hs_uee); 
    alert0=t_2d(x,y-1,z,t[x][y][z0],t[x][y][z],hs_bf,hs_ff); 
    if(alert1) longflags[x*ny+y-1]=1;
    if(alert0) longflags[x*ny+y-1]=0;
  } 
 
  for(x=x_start;x>x_begin;x--){
    for(y=y_start;y>y_begin;y--){

      hs_ff=hs[x-1][y-1][z_s];
      if(x>1) hs_bf=hs[x-2][y-1][z_s];
      else hs_bf=T3D_INF;
      if(y>1) hs_fb=hs[x-1][y-2][z_s];
      else hs_fb=T3D_INF;
      if(z_sf>=0 && z_sf<nmesh_z){
        hs_uee=hs[x-1][y-1][z_sf];
        if(x>1) hs_ube=hs[x-2][y-1][z_sf];
        else hs_ube=T3D_INF;
        if(y>1) hs_ueb=hs[x-1][y-2][z_sf];
        else hs_ueb=T3D_INF;
      }
      else hs_uee=hs_ube=hs_ueb=T3D_INF;

      /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
      alert0=edge_diff(x-1,y-1,z,t[x][y][z0],t[x][y][z],hs_ff)
        +t_3d(x-1,y-1,z,t[x][y][z0],
              t[x-1][y][z0],t[x][y][z],t[x-1][y][z],hs_ff)
        +t_3d(x-1,y-1,z,t[x][y][z0],
              t[x][y-1][z0],t[x][y][z],t[x][y-1][z],hs_ff);
      if(alert0){
        updated+=alert0;
        longflags[x*ny-ny+y-1]=0;
      }

      /* interface waves along x_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x-1,y-1,z,t[x][y-1][z],hs_ff,hs_fb,hs_uee,hs_ueb);
      alert0=t_2d(x-1,y-1,z,t[x][y-1][z0],t[x][y-1][z],hs_ff,hs_fb);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_bf++; /* this scan must be (re-)examined */
        if(x_start_bf<x) x_start_bf=x;
        if(y_start_bf>y-1) y_start_bf=y-1;
        if(alert1) longflags[x*ny-ny+y-1]=1;
        else       longflags[x*ny-ny+y-1]=0;
      }

      /* interface waves along y_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x-1,y-1,z,t[x-1][y][z],hs_ff,hs_bf,hs_uee,hs_ube);
      alert0=t_2d(x-1,y-1,z,t[x-1][y][z0],t[x-1][y][z],hs_ff,hs_bf);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_fb++; /* this scan must be (re-)examined */
        if(x_start_fb>x-1) x_start_fb=x-1;
        if(y_start_fb<y) y_start_fb=y;
        if(alert1) longflags[x*ny-ny+y-1]=1;
        else       longflags[x*ny-ny+y-1]=0;
      }

      /* interface waves along z_side : 2 2D transmission and 1 2D diffraction */
      alert1=diff_2d(x-1,y-1,z,t[x][y][z],hs_ff,hs_uee)
        +t_2d(x-1,y-1,z,t[x][y][z],t[x-1][y][z],hs_ff,hs_uee)
        +t_2d(x-1,y-1,z,t[x][y][z],t[x][y-1][z],hs_ff,hs_uee);
      if(alert1){
        updated+=alert1;
        longflags[x*ny-ny+y-1]=1;
      }
    }
  }

  flag_bb=0;
  x_start_bb=x_begin;
  y_start_bb=y_begin;

  return updated;
}

/*--------------------------------------Z_SIDE() : SCAN_Z_EB()--------------*/

static int scan_z_fb(int x_start, int x_end, int y_begin, int y_start,
                     int z0, int z, int z_s)

{    int
    updated=0,
    alert0,alert1,
    x,y,
    z_sf;
  double
    hs_ff,hs_bb,hs_bf,hs_uee,hs_ubb,hs_ube;

  z_sf=z_s+z-z0;

  hs_bf=hs_ubb=hs_ube=T3D_INF; 
  for(x=x_start,y=y_start;x<x_end;x++){ 
    if(y) hs_bf=hs[x][y-1][z_s]; 
    hs_bb=hs[x][y][z_s]; 
    if(z_sf>=0 && z_sf<nmesh_z){ 
      if(y) hs_ube=hs[x][y-1][z_sf]; 
      hs_ubb=hs[x][y][z_sf]; 
    }    
    alert1=t_1d(x+1,y,z,t[x][y][z],hs_bf,hs_bb,hs_ube,hs_ubb);
    alert0=t_2d(x+1,y,z,t[x][y][z0],t[x][y][z],hs_bf,hs_bb); 
    if(alert1) longflags[x*ny+ny+y]=1;
    if(alert0) longflags[x*ny+ny+y]=0;
  }   
  
  hs_ff=hs_uee=hs_ube=T3D_INF; 
  for(x=x_start,y=y_start;y>y_begin;y--){ 
    if(x) hs_ff=hs[x-1][y-1][z_s]; 
    hs_bf=hs[x][y-1][z_s]; 
    if(z_sf>=0 && z_sf<nmesh_z){ 
      if(x) hs_uee=hs[x-1][y-1][z_sf]; 
      hs_ube=hs[x][y-1][z_sf]; 
    } 
    alert1=t_1d(x,y-1,z,t[x][y][z],hs_bf,hs_ff,hs_ube,hs_uee); 
    alert0=t_2d(x,y-1,z,t[x][y][z0],t[x][y][z],hs_bf,hs_ff); 
    if(alert1) longflags[x*ny+y-1]=1;
    if(alert0) longflags[x*ny+y-1]=0;
  }   
 
  for(x=x_start;x<x_end;x++){
    for(y=y_start;y>y_begin;y--){

      hs_bf=hs[x][y-1][z_s];
      if(y>1) hs_bb=hs[x][y-2][z_s];
      else hs_bb=T3D_INF;
      hs_ff=hs[x+1][y-1][z_s];
      if(z_sf>=0 && z_sf<nmesh_z){
        hs_ube=hs[x][y-1][z_sf];
        if(y>1) hs_ubb=hs[x][y-2][z_sf];
        else hs_ubb=T3D_INF;
        hs_uee=hs[x+1][y-1][z_sf];
      }
      else hs_ubb=hs_ube=hs_uee=T3D_INF;

      /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
      alert0=edge_diff(x+1,y-1,z,t[x][y][z0],t[x][y][z],hs_bf)
        +t_3d(x+1,y-1,z,t[x][y][z0],
              t[x+1][y][z0],t[x][y][z],t[x+1][y][z],hs_bf)
        +t_3d(x+1,y-1,z,t[x][y][z0],
              t[x][y-1][z0],t[x][y][z],t[x][y-1][z],hs_bf);
      if(alert0){
        updated+=alert0;
        longflags[x*ny+ny+y-1]=0;
      }

      /* interface waves along x_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x+1,y-1,z,t[x][y-1][z],hs_bb,hs_bf,hs_ubb,hs_ube);
      alert0=t_2d(x+1,y-1,z,t[x][y-1][z0],t[x][y-1][z],hs_bb,hs_bf);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_ff++; /* this scan must be (re-)examined */
        if(x_start_ff>x) x_start_ff=x;
        if(y_start_ff>y-1) y_start_ff=y-1;
        if(alert1) longflags[x*ny+ny+y-1]=1;
        else       longflags[x*ny+ny+y-1]=0;
      }

      /* interface waves along y_side: 1 1D transmission and 1 2D transmission */
      alert1=t_1d(x+1,y-1,z,t[x+1][y][z],hs_ff,hs_bf,hs_uee,hs_ube);
      alert0=t_2d(x+1,y-1,z,t[x+1][y][z0],t[x+1][y][z],hs_ff,hs_bf);
      if(alert0+alert1){
        updated+=alert0+alert1;
        flag_bb++; /* this scan must be (re-)examined */
        if(x_start_bb<x+1) x_start_bb=x+1;
        if(y_start_bb<y) y_start_bb=y;
        if(alert1) longflags[x*ny+ny+y-1]=1;
        else       longflags[x*ny+ny+y-1]=0;
      }

      /* interface waves along z_side : 2 2D transmission and 1 2D diffraction */
      alert1=diff_2d(x+1,y-1,z,t[x][y][z],hs_bf,hs_ube)
        +t_2d(x+1,y-1,z,t[x][y][z],t[x+1][y][z],hs_bf,hs_ube)
        +t_2d(x+1,y-1,z,t[x][y][z],t[x][y-1][z],hs_bf,hs_ube);
      if(alert1){
        updated+=alert1;
        longflags[x*ny+ny+y-1]=1;
      }
    }
  }

  flag_fb=0;
  x_start_fb=x_end;
  y_start_fb=y_begin;

  return updated;
}
/*--------------------------------------------END_Z_SIDE()--------------------*/
/*********************************************END_TIME_3D**********************/
