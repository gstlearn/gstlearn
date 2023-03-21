/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "geoslib_old_f.h"
#include "geoslib_enum.h"

#include "Enum/EJustify.hpp"

#include "Basic/Law.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Basic/OptDbg.hpp"

#include <math.h>

/*! \cond */
#define LTABLE 8
#define NTABLE 513

#define ARRAY(ix,iy)      (array[(iy) * NX + (ix)])
#define LIMIT(ix,iy)      (limit[(iy) * NX + (ix)])
#define VV(itr,iz)        (db_v->getArray(iatt_v,NTRACE * (iz) + itr) / VFACT)
#define COVTAB(ivar,jvar) (covtab[(jvar) * NVAR + (ivar)])
#define C00(ivar,jvar)    (c00   [(jvar) * NVAR + (ivar)])
#define LB(ivar,jvar)     (lb    [(jvar) * NVAR + (ivar)])
#define IND(iech,ivar)    ((iech) + (ivar) * nech)
#define LHS(i,iv,j,jv)    (lhs [IND(i,iv) + neqmax * IND(j,jv)])
#define RHS(i,iv,jv)      (rhs [IND(i,iv) + neqmax * (jv)])
/*! \endcond */

static double *covtab;
static double DX, DZ;
static int NX, NY, NZ, NVAR, NTRACE;
static double VFACT = 1000.;
static int IECH_OUT = -1;
typedef struct
{
  int nvois;
  int nactive;
  int n_v1;
  int n_v2;
  int *ix_ngh;
  int *iz_ngh;
  double *v1_ngh;
  double *v2_ngh;
} ST_Seismic_Neigh;

/****************************************************************************/
/*!
 **  Calculate the extrema of the velocity array vv[]
 **
 ** \return  Error return code
 **
 ** \param[in]  nech Number of samples in the array vv[]
 ** \param[in]  vv   Array of interest
 **
 ** \param[out]  v0   Initial velocity value
 ** \param[out]  v1   Last velocity value
 ** \param[out]  vmin Minimum value
 ** \param[out]  vmax Maximum value
 **
 *****************************************************************************/
static int st_velocity_minmax(int nech,
                              double *vv,
                              double *v0,
                              double *v1,
                              double *vmin,
                              double *vmax)
{
  int i, number;
  double vvdef, delta;

  (*v0) = 1.e30;
  (*v1) = 1.e30;
  (*vmin) = 1.e30;
  (*vmax) = -1.e30;

  /* Check the extreme (defined) values */

  for (i = number = 0; i < nech; i++)
  {
    if (FFFF(vv[i])) continue;
    if (vv[i] <= 0) continue;
    if (vv[i] < (*vmin)) (*vmin) = vv[i];
    if (vv[i] > (*vmax)) (*vmax) = vv[i];
    number++;
  }

  /* Final checks */

  if (number <= 0)
  {
    messerr("The velocity field is not defined: no active value");
    return (1);
  }

  /* All velocity values are correctly defined */

  if (number == nech) return (0);

  /* Calculate an arbitrary small velocity value */

  delta = (*vmax) - (*vmin);
  vvdef = (*vmin) - delta / 2.;
  if (delta <= 0 || vvdef <= 0.) vvdef = (*vmin) / 2.;

  /* Replace unknown or negative or zero values by an arbitrary small value */

  for (i = 0; i < nech; i++)
  {
    if (FFFF(vv[i]) || vv[i] <= 0.) vv[i] = vvdef;
  }
  (*v0) = vv[0];
  (*v1) = vv[nech - 1];

  return (0);
}

/****************************************************************************/
/*!
 **  Compute regularly-sampled monotonically increasing function x(y)
 **  from regularly-sampled monotonically increasing function y(x)
 **  by inverse linear interpolation
 **
 ** \param[in]  nx   number of samples of y(x)
 ** \param[in]  dx   x sampling interval; dx>0.0 is required
 ** \param[in]  x0   first x
 ** \param[in]  y    array[nx] of y(x) values; y[0] < ... < y[nx-1] required
 ** \param[in]  ny   number of samples of x(y)
 ** \param[in]  dy   y sampling interval; dy>0.0 is required
 ** \param[in]  y0   first y
 ** \param[in]  xylo x value assigned to x(y) when y is less than smallest y(x)
 ** \param[in]  xyhi x value assigned to x(y) when y is greater than largest y(x)
 **
 ** \param[out]  x   array[ny] of x(y) values
 **
 *****************************************************************************/
static void st_yxtoxy(int nx,
                      double dx,
                      double x0,
                      double *y,
                      int ny,
                      double dy,
                      double y0,
                      double xylo,
                      double xyhi,
                      double *x)
{
  int nxi, nyo, jxi1, jxi2, jyo;
  double dxi, fxi, dyo, fyo, fyi, yo, xi1, yi1, yi2;

  nxi = nx;
  dxi = dx;
  fxi = x0;
  nyo = ny;
  dyo = dy;
  fyo = y0;
  fyi = y[0];

  /* loop over output y less than smallest input y */
  for (jyo = 0, yo = fyo; jyo < nyo; jyo++, yo += dyo)
  {
    if (yo >= fyi) break;
    x[jyo] = xylo;
  }

  /* loop over output y between smallest and largest input y */
  if (jyo == nyo - 1 && yo == fyi)
  {
    x[jyo++] = fxi;
    yo += dyo;
  }
  jxi1 = 0;
  jxi2 = 1;
  xi1 = fxi;
  while (jxi2 < nxi && jyo < nyo)
  {
    yi1 = y[jxi1];
    yi2 = y[jxi2];
    if (yi1 <= yo && yo <= yi2)
    {
      x[jyo++] = xi1 + dxi * (yo - yi1) / (yi2 - yi1);
      yo += dyo;
    }
    else
    {
      jxi1++;
      jxi2++;
      xi1 += dxi;
    }
  }

  /* loop over output y greater than largest input y */
  while (jyo < nyo)
    x[jyo++] = xyhi;
}

/****************************************************************************/
/*!
 **  Function for tabulating dsinc()
 **
 ** \return  Value of the dsinc() function
 **
 ** \param[in]  x Argument of the function
 **
 *****************************************************************************/
static double dsinc(double x)

{
  double pix;

  if (x == 0.0) return 1.0;
  pix = GV_PI * x;
  return (sin(pix) / pix);
}

/****************************************************************************/
/*!
 **  Solve a symmetric Toeplitz linear system of equations Rf=g for f
 **
 ** \param[in]  n dimension of system
 ** \param[in]  r array[n] of top row of Toeplitz matrix
 ** \param[in]  g array[n] of right-hand-side column vector
 **
 ** \param[out]  f array[n] of solution (left-hand-side) column vector
 ** \param[out]  a array[n] of solution to Ra=v (Claerbout, FGDP, p. 57)
 **
 ** \remark This routine does NOT solve the case when the main diagonal is
 ** \remark zero, it just silently returns.
 ** \remark The left column of the Toeplitz matrix is assumed to be equal to
 ** \remark the top row (as specified in r); i.e., the Toeplitz matrix is
 ** \remark assumed symmetric.
 **
 *****************************************************************************/
static void stoepd(int n, double *r, double *g, double *f, double *a)
{
  int i, j;
  double v, e, c, w, bot;

  if (r[0] == 0.0) return;

  a[0] = 1.0;
  v = r[0];
  f[0] = g[0] / r[0];

  for (j = 1; j < n; j++)
  {

    /* solve Ra=v as in Claerbout, FGDP, p. 57 */
    a[j] = 0.0;
    f[j] = 0.0;
    for (i = 0, e = 0.0; i < j; i++)
      e += a[i] * r[j - i];
    c = e / v;
    v -= c * e;
    for (i = 0; i <= j / 2; i++)
    {
      bot = a[j - i] - c * a[i];
      a[i] -= c * a[j - i];
      a[j - i] = bot;
    }

    /* use a and v above to get f[i], i = 0,1,2,...,j */
    for (i = 0, w = 0.0; i < j; i++)
      w += f[i] * r[j - i];
    c = (w - g[j]) / v;
    for (i = 0; i <= j; i++)
      f[i] -= c * a[j - i];
  }
}

/****************************************************************************/
/*!
 **  Compute least-squares optimal sinc interpolation coefficients.
 **
 ** \param[in]  d     fractional distance to interpolation point; 0.0<=d<=1.0
 ** \param[in]  lsinc length of sinc approximation; lsinc%2==0 and lsinc<=20
 **
 ** \param[out]  sinc array[lsinc] containing interpolation coefficients
 **
 ** \remark The coefficients are a least-squares-best approximation to the
 ** \remark ideal sinc function for frequencies from zero up to a computed
 ** \remark maximum frequency.
 ** \remark For a given interpolator length, lsinc, mksinc computes the
 ** \remark maximum frequency, fmax (expressed as a fraction of the Nyquist
 ** \remark frequency), using the following empirically derived relation (from
 ** \remark a Western Geophysical Technical Memorandum by Ken Larner):
 ** \remark        fmax = min(0.066+0.265*log(lsinc),1.0)
 ** \remark Note that fmax increases as lsinc increases, up to a maximum of 1.
 ** \remark Use the coefficients to interpolate a uniformly-sampled function
 ** \remark y(i) as follows:
 ** \remark             lsinc-1
 ** \remark     y(i+d) =  sum  sinc[j]*y(i+j+1-lsinc/2)
 ** \remark               j=0
 ** \remark Interpolation error is greatest for d=0.5, but for frequencies less
 ** \remark than fmax, the error should be less than 1.0 percent.

 *****************************************************************************/
static void st_mksinc(double d, int lsinc, double *sinc)
{
  int j;
  double s[20], a[20], c[20], work[20], fmax;

  /* compute auto-correlation and cross-correlation arrays */
  fmax = 0.066 + 0.265 * log((double) lsinc);
  fmax = (fmax < 1.0) ? fmax :
                        1.0;
  for (j = 0; j < lsinc; j++)
  {
    a[j] = dsinc(fmax * j);
    c[j] = dsinc(fmax * (lsinc / 2 - j - 1 + d));
  }

  /* solve symmetric Toeplitz system for the sinc approximation */
  stoepd(lsinc, a, c, s, work);
  for (j = 0; j < lsinc; j++)
    sinc[j] = s[j];
}

/****************************************************************************/
/*!
 **  Interpolation of a uniformly-sampled complex function y(x)
 **  via a table of 8-coefficient interpolators
 **
 ** \param[in]  ntable  number of tabulated interpolation operators; ntable>=2
 ** \param[in]  table   array of tabulated 8-point interpolation operators
 ** \param[in]  nxin    number of x values at which y(x) is input
 ** \param[in]  dxin    x sampling interval for input y(x)
 ** \param[in]  fxin    x value of first sample input
 ** \param[in]  yin     array of input y(x) values:  yin[0] = y(fxin), etc.
 ** \param[in]  yinl    value used to extrapolate yin values to left of yin[0]
 ** \param[in]  yinr    value used to extrapolate yin values to right of
 **                     yin[nxin-1]
 ** \param[in]  nxout   number of x values a which y(x) is output
 ** \param[in]  xout    array of x values at which y(x) is output
 **
 ** \param[out]  yout   array of output y(x) values:  yout[0] = y(xout[0]), etc.
 **
 ** \remark ntable must not be less than 2.
 ** \remark The table of interpolation operators must be as follows:
 ** \remark Let d be the distance, expressed as a fraction of dxin, from a
 ** \remark particular xout value to the sampled location xin just to the left
 ** \remark of xout.  Then:
 ** \remark for d = 0., table[0][0:7] = 0., 0., 0., 1., 0., 0., 0., 0.
 ** \remark are the weights applied to the 8 input samples nearest xout.
 ** \remark for d = 1., table[ntable-1][0:7] = 0., 0., 0., 0., 1., 0., 0., 0.
 ** \remark are the weights applied to the 8 input samples nearest xout.
 ** \remark for d = (float)itable/(float)(ntable-1), table[itable][0:7] are the
 ** \remark weights applied to the 8 input samples nearest xout.
 ** \remark If the actual sample distance d does not exactly equal one of the
 ** \remark values for which interpolators are tabulated, then the interpolator
 ** \remark corresponding to the nearest value of d is used.
 ** \remark Because extrapolation of the input function y(x) is defined by the
 ** \remark left and right values yinl and yinr, the xout values are not
 ** \remark restricted to lie within the range of sample locations defined by
 ** \remark nxin, dxin, and fxin.
 *****************************************************************************/
static void st_intt8r(int ntable,
                      double table[][LTABLE],
                      int nxin,
                      double dxin,
                      double fxin,
                      double *yin,
                      double yinl,
                      double yinr,
                      int nxout,
                      double *xout,
                      double *yout)
{
  int ioutb, nxinm8, ixout, ixoutn, kyin, ktable, itable;
  double xoutb, xoutf, xoutn, frac, fntablem1, yini, sum;

  /* compute constants */
  ioutb = -3 - 8;
  xoutf = fxin;
  xoutb = 8.0 - xoutf / dxin;
  fntablem1 = (double) (ntable - 1);
  nxinm8 = nxin - 8;

  /* loop over output samples */
  for (ixout = 0; ixout < nxout; ixout++)
  {

    /* determine pointers into table and yin */
    sum = 0.;
    xoutn = xoutb + xout[ixout] / dxin;
    ixoutn = (int) xoutn;
    kyin = ioutb + ixoutn;
    frac = xoutn - (double) ixoutn;
    ktable = (int) ((frac >= 0.0) ? frac * fntablem1 + 0.5 :
                                    (frac + 1.0) * fntablem1 - 0.5);

    /* if totally within input array, use fast method */
    if (kyin >= 0 && kyin <= nxinm8)
    {
      sum = yin[kyin + 0] * table[ktable][0] + yin[kyin + 1] * table[ktable][1]
            + yin[kyin + 2] * table[ktable][2]
            + yin[kyin + 3] * table[ktable][3]
            + yin[kyin + 4] * table[ktable][4]
            + yin[kyin + 5] * table[ktable][5]
            + yin[kyin + 6] * table[ktable][6]
            + yin[kyin + 7] * table[ktable][7];

      /* else handle end effects with care */
    }
    else
    {

      /* sum over 8 tabulated coefficients */
      for (itable = 0; itable < 8; itable++, kyin++)
      {
        if (kyin < 0)
          yini = yinl;
        else if (kyin >= nxin)
          yini = yinr;
        else
          yini = yin[kyin];
        sum += yini * table[ktable][itable];
      }
    }
    yout[ixout] = sum;
  }
}

/****************************************************************************/
/*!
 **  Interpolation of a uniformly-sampled real function y(x) via a
 **  table of 8-coefficient sinc approximations
 **
 ** \param[out]  table array of weights
 **
 *****************************************************************************/
static void st_weights(double table[][LTABLE])

{
  int jtable;
  double frac;

  /* Tabulate sinc interpolation coefficients */

  for (jtable = 1; jtable < NTABLE - 1; jtable++)
  {
    frac = (double) jtable / (double) (NTABLE - 1);
    st_mksinc(frac, LTABLE, &table[jtable][0]);
  }
  for (jtable = 0; jtable < LTABLE; jtable++)
  {
    table[0][jtable] = 0.0;
    table[NTABLE - 1][jtable] = 0.0;
  }
  table[0][LTABLE / 2 - 1] = 1.0;
  table[NTABLE - 1][LTABLE / 2] = 1.0;

  return;
}

/****************************************************************************/
/*!
 **  Debugging printout for the seismic time-to-depth conversion
 **
 ** \param[in]  rankz   Rank for Depth (0 for Input; 1 for Output)
 ** \param[in]  nz      number of depth samples
 ** \param[in]  z0      first depth value
 ** \param[in]  dz      depth sampling interval
 ** \param[in]  rankt   Rank for Time (0 for Input; 1 for Output)
 ** \param[in]  nt      number of time samples
 ** \param[in]  t0      first time value
 ** \param[in]  dt      time sampling interval
 ** \param[in]  vmin    Minimum velocity value
 ** \param[in]  vmax    Maximum velocity value
 **
 ****************************************************************************/
static void st_seismic_debug(int rankz,
                             int nz,
                             double z0,
                             double dz,
                             int rankt,
                             int nt,
                             double t0,
                             double dt,
                             double vmin,
                             double vmax)
{
  int i;

  for (i = 0; i < 2; i++)
  {
    if (i == 0)
      message("Input:\n");
    else
      message("Output:\n");

    if (i == rankz)
    {
      message("\tNumber of depth samples = %d\n", nz);
      message("\tDepth sampling interval = %g (m)\n", dz);
      message("\tDepth of first sample   = %g (m)\n", z0);
      message("\tDepth of last sample    = %g (m)\n", z0 + (nz - 1) * dz);
    }

    if (i == rankt)
    {
      message("\tNumber of time samples = %d\n", nt);
      message("\tTime sampling interval = %g (ms)\n", dt);
      message("\tTime of first sample   = %g (ms)\n", t0);
      message("\tTime of last sample    = %g (ms)\n", t0 + (nt - 1) * dt);
    }
    if (i == 0)
    {
      message("Velocity:\n");
      message("\tMinimum value          = %g (m/s)\n", vmin);
      message("\tMaximum value          = %g (m/s)\n", vmax);
    }
  }
  message("\n");
}

/****************************************************************************/
/*!
 **  Define the Time Grid characteristics from the Depth Grid
 **
 ** \return  Error return code
 **
 ** \param[in]  verbose Verbose flag
 ** \param[in]  db_z    Depth Grid structure
 ** \param[in]  iatt_v  Attribute address of the Velocity (int Depth Grid)
 **
 ** \param[out] nx      Number of grid nodes along each direction
 ** \param[out] x0      Origin of the grid along each direction
 ** \param[out] dx      Mesh of the grid along each direction
 **
 *****************************************************************************/
int seismic_z2t_grid(int verbose,
                     DbGrid *db_z,
                     int iatt_v,
                     int *nx,
                     double *x0,
                     double *dx)
{
  double *vv, z0, t0, v0, v1, dz, dt, vmin, vmax;
  int ndim, nech, nt, nz, error, i;

  /* Initializations */

  if (! db_z->isGrid())
  {
    messerr("This procedure requires an input Grid Db");
    return (1);
  }
  error = 1;
  vv = nullptr;
  ndim = db_z->getNDim();
  nech = db_z->getSampleNumber();

  /* Core allocation */

  vv = (double*) mem_alloc(nech * sizeof(double), 0);
  if (vv == nullptr) goto label_end;
  for (i = 0; i < ndim; i++)
  {
    nx[i] = db_z->getNX(i);
    dx[i] = db_z->getDX(i);
    x0[i] = db_z->getX0(i);
  }

  /* Read the velocity variable */

  if (db_vector_get_att(db_z, iatt_v, vv)) goto label_end;
  if (st_velocity_minmax(nech, vv, &v0, &v1, &vmin, &vmax)) goto label_end;

  /* Update the vertical direction */

  nz = db_z->getNX(ndim - 1);
  z0 = db_z->getX0(ndim - 1);
  dz = db_z->getDX(ndim - 1);
  dt = 2. * dz / vmin;
  t0 = 2. * z0 / v0;
  nt = (int) (1 + (nz - 1) * (2. * dz) / (dt * vmax));
  dt *= VFACT;
  t0 *= VFACT;
  dx[ndim - 1] = dt;
  x0[ndim - 1] = t0;
  nx[ndim - 1] = nt;

  /* Optional printout */

  if (verbose) st_seismic_debug(0, nz, z0, dz, 1, nt, t0, dt, vmin, vmax);

  /* Set the error return code */

  error = 0;

  label_end: vv = (double*) mem_free((char* ) vv);
  return (error);
}

/****************************************************************************/
/*!
 **  Define the Depth Grid characteristics from the Time Grid
 **
 ** \return  Error return code
 **
 ** \param[in]  verbose Verbose flag
 ** \param[in]  db_t    Time Grid structure
 ** \param[in]  iatt_v  Attribute address of the Velocity (in Time grid)
 **
 ** \param[out] nx      Number of grid nodes along each direction
 ** \param[out] x0      Origin of the grid along each direction
 ** \param[out] dx      Mesh of the grid along each direction
 **
 *****************************************************************************/
int seismic_t2z_grid(int verbose,
                     DbGrid *db_t,
                     int iatt_v,
                     int *nx,
                     double *x0,
                     double *dx)
{
  double *vv, z0, t0, v0, v1, dz, dt, vmin, vmax;
  int ndim, nech, nt, nz, error, i;

  /* Initializations */

  if (! db_t->isGrid())
  {
    messerr("This procedure requires an input Grid Db");
    return (1);
  }
  error = 1;
  vv = nullptr;
  ndim = db_t->getNDim();
  nech = db_t->getSampleNumber();

  /* Core allocation */

  vv = (double*) mem_alloc(nech * sizeof(double), 0);
  if (vv == nullptr) goto label_end;
  for (i = 0; i < ndim; i++)
  {
    nx[i] = db_t->getNX(i);
    dx[i] = db_t->getDX(i);
    x0[i] = db_t->getX0(i);
  }

  /* Read the velocity variable */

  if (db_vector_get_att(db_t, iatt_v, vv)) goto label_end;
  if (st_velocity_minmax(nech, vv, &v0, &v1, &vmin, &vmax)) goto label_end;

  /* Update the vertical direction */

  nt = db_t->getNX(ndim - 1);
  t0 = db_t->getX0(ndim - 1);
  dt = db_t->getDX(ndim - 1);
  dz = vmin * dt / 2.;
  z0 = v0 * t0 / 2.;
  nz = (int) (1 + (nt - 1) * (dt * vmax) / (2. * dz));
  dz /= VFACT;
  z0 /= VFACT;
  dx[ndim - 1] = dz;
  x0[ndim - 1] = z0;
  nx[ndim - 1] = nz;

  /* Optional printout */

  if (verbose) st_seismic_debug(1, nz, z0, dz, 0, nt, t0, dt, vmin, vmax);

  /* Set the error return code */

  error = 0;

  label_end: vv = (double*) mem_free((char* ) vv);
  return (error);
}

/****************************************************************************/
/*!
 **  Copy a trace between the Db and a given array
 **
 ** \param[in]  mode  Type of operation
 ** \li               0 : Copy from Db into array
 ** \li               1 : Copy from array into Db
 ** \param[in]  db    Db structure
 ** \param[in]  iatt  Rank of the variable
 ** \param[in]  ival  Rank of the column
 ** \param[in]  tab   Array containing the column of values
 **
 *****************************************************************************/
static void st_copy(int mode, DbGrid *db, int iatt, int ival, double *tab)
{
  int ndim = db->getNDim();
  int nech = db->getSampleNumber();
  int nval = db->getNX(ndim - 1);
  int nby = nech / nval;

  /* Dispatch */

  switch (mode)
  {
    case 0: /* Copy from Db into array */
      for (int i = 0; i < nval; i++)
        tab[i] = db->getArray(nby * i + ival, iatt);
      break;

    case 1: /* Copy from array into Db */
      for (int i = 0; i < nval; i++)
        db->setArray(iatt + nby * i + ival, iatt, tab[i]);
      break;
  }
  return;
}

/****************************************************************************/
/*!
 **  Check that the Depth and Time grid coincide
 **
 ** \return  1 if the two files coincide; 0 otherwise
 **
 ** \param[in]  db_z    Depth Grid structure
 ** \param[in]  db_t    Time Grid structure (optional)
 **
 ** \remark  In case of error, a message is displayed
 **
 *****************************************************************************/
static int st_match(DbGrid *db_z, DbGrid *db_t)
{
  int idim, ndim, nech, nz, error;

  /* Initializations */

  NTRACE = 0;
  error = 1;
  nech = db_z->getSampleNumber();
  ndim = db_z->getNDim();
  nz = db_z->getNX(ndim - 1);

  /* Check the grid characteristics */

  if (db_t != nullptr)
  {
    if (!db_t->isSameGrid(db_z->getGrid())) goto label_end;
  }

  /* Set the error return code */

  error = 0;

  label_end: if (error)
  {
    messerr("Error for one of the following reasons:");
    messerr("- Different space dimensions: Depth(%d) - Time(%d)",
            db_z->getNDim(), db_t->getNDim());
    for (idim = 0; idim < ndim - 1; idim++)
    {
      messerr("- Number of mesh (axe #%d): Depth(%d) - Time(%d)",
              db_z->getNX(idim), db_t->getNX(idim));
      messerr("- Origin of the Grid (axe #%d): Depth(%d) - Time(%d)",
              db_z->getX0(idim), db_t->getX0(idim));
      messerr("- Mesh of the Grid (axe #%d): Depth(%d) - Time(%d)",
              db_z->getDX(idim), db_t->getDX(idim));
    }
  }
  else
    NTRACE = nech / nz;

  return (error);
}

/****************************************************************************/
/*!
**  Resample from depth to time
**
** \param[in]  db_z     Depth Db
** \param[in]  iatt_z   Address of the first variable of the Depth Db
** \param[in]  nz       Number of meshes of the Depth Db
** \param[in]  z0       First Depth
** \param[in]  dz       Mesh of the Depth Db
** \param[in]  db_t     Time Db
** \param[in]  iatt_t   Address of the first variable of the Time Db
** \param[in]  nt       Number of meshes of the Time Db
** \param[in]  t0       First Time
** \param[in]  t1       Last Time
** \param[in]  dt       Mesh of the Time Db
** \param[in]  db_v     Velocity Db
** \param[in]  iatt_v   Address of the velocity variable in Depth Db
** \param[in]  natt     Number of seismic attributes
**
** \param[out] tz       Working array (dimension: nt)
** \param[out] zt       Working array (dimension: nz)
** \param[out] at       Working array (dimension: nt)
** \param[out] az       Working array (dimension: nz)
**
** \remark Linear interpolation and constant extrapolation is used to
** \remark determine interval velocities at times not specified.
**
*****************************************************************************/
static void st_seismic_z2t_convert(DbGrid *db_z,
                                   int    iatt_z,
                                   int    nz,
                                   double z0,
                                   double /*z1*/,
                                   double dz,
                                   DbGrid *db_t,
                                   int iatt_t,
                                   int nt,
                                   double t0,
                                   double t1,
                                   double dt,
                                   DbGrid *db_v,
                                   int iatt_v,
                                   int natt,
                                   double *tz,
                                   double *zt,
                                   double *at,
                                   double *az)
{
  double t, vz0, vz1, table[NTABLE][LTABLE];
  int itrace, iz, it, iatt;

  /* Calculate the interpolation weights */

  st_weights(table);

  /* Loop on the vertical lines */

  for (itrace = 0; itrace < NTRACE; itrace++)
  {

    /* Transform v(z) into t(z) */

    tz[0] = 2. * z0 / VV(itrace, 0);
    for (iz = 1; iz < nz; ++iz)
      tz[iz] = tz[iz - 1] + 2. * dz / VV(itrace, iz - 1);

    /* Transform t(z) into z(t) */

    vz0 = VV(itrace, 0);
    vz1 = VV(itrace, nz - 1);
    st_yxtoxy(nz, dz, z0, tz, nt, dt, t0, 0., 0., zt);

    /* for t values before t0, use first velocity to calculate z(t) */
    for (it = 0, t = t0; t <= tz[0]; ++it, t += dt)
      zt[it] = vz0 * t / 2.;

    /* for t values from t1 down, calculate z(t) */
    for (it = nt - 1, t = t1; t >= tz[nz - 1]; it--, t -= dt)
      zt[it] = vz1 * t1 / 2. + vz1 * (t - tz[nz - 1]) / 2.;

    /* Convert attributes from Depth to Time */

    for (iatt = 0; iatt < natt; iatt++)
    {
      st_copy(0, db_z, iatt_z + iatt, itrace, az);
      st_intt8r(NTABLE, table, nz, dz, z0, az, 0., 0., nt, zt, at);
      st_copy(1, db_t, iatt_t + iatt, itrace, at);
    }
  }

  return;
}

/****************************************************************************/
/*!
 **  Resample from time to depth
 **
 ** \param[in]  db_t     Time Db
 ** \param[in]  iatt_t   Address of the first variable of the Time Db
 ** \param[in]  nt       Number of meshes of the Time Db
 ** \param[in]  t0       First Time
 ** \param[in]  t1       Last Time
 ** \param[in]  dt       Mesh of the Time Db
 ** \param[in]  db_z     Depth Db
 ** \param[in]  iatt_z   Address of the first variable of the Depth Db
 ** \param[in]  nz       Number of meshes of the Depth Db
 ** \param[in]  z0       First Depth
 ** \param[in]  z1       Last Depth
 ** \param[in]  dz       Mesh of the Depth Db
 ** \param[in]  db_v     Velocity Db
 ** \param[in]  iatt_v   Address of the velocity variable in Depth Db
 ** \param[in]  natt     Number of seismic attributes
 **
 ** \param[out] tz       Working array (dimension: nt)
 ** \param[out] zt       Working array (dimension: nz)
 ** \param[out] at       Working array (dimension: nt)
 ** \param[out] az       Working array (dimension: nz)
 **
 ** \remark Linear interpolation and constant extrapolation is used to
 ** \remark determine interval velocities at times not specified.
 **
 *****************************************************************************/
static void st_seismic_t2z_convert(DbGrid *db_t,
                                   int iatt_t,
                                   int nt,
                                   double t0,
                                   double t1,
                                   double dt,
                                   DbGrid *db_z,
                                   int iatt_z,
                                   int nz,
                                   double z0,
                                   double z1,
                                   double dz,
                                   DbGrid *db_v,
                                   int iatt_v,
                                   int natt,
                                   double *tz,
                                   double *zt,
                                   double *at,
                                   double *az)
{
  double z, vt0, vt1, table[NTABLE][LTABLE];
  int itrace, it, iz, iatt;

  /* Calculate the interpolation weights */

  st_weights(table);

  /* Loop on the vertical traces */

  for (itrace = 0; itrace < NTRACE; itrace++)
  {

    /* Transform v(t) into z(t) */

    zt[0] = t0 * VV(itrace, 0) / 2.;
    for (it = 1; it < nt; it++)
      zt[it] = zt[it - 1] + dt * VV(itrace, it - 1) / 2.;

    /* Transform z(t) into t(z) */

    vt0 = VV(itrace, 0);
    vt1 = VV(itrace, nt - 1);
    st_yxtoxy(nt, dt, t0, zt, nz, dz, z0, 0., 0., tz);

    /* for z values before z0, use first velocity to calculate t(z) */
    for (iz = 0, z = z0; z <= zt[0]; iz++, z += dz)
      tz[iz] = 2. * z / vt0;

    /* for z values from z1 down, calculate t(z) */
    for (iz = nz - 1, z = z1; z >= zt[nt - 1]; iz--, z -= dz)
      tz[iz] = t1 + 2. * (z - zt[nt - 1]) / vt1;

    /* Convert attributes from Time to Depth */

    for (iatt = 0; iatt < natt; iatt++)
    {
      st_copy(0, db_t, iatt_t + iatt, itrace, at);
      st_intt8r(NTABLE, table, nt, dt, t0, at, 0., 0., nz, tz, az);
      st_copy(1, db_z, iatt_z + iatt, itrace, az);
    }
  }

  return;
}

/****************************************************************************/
/*!
 **  Local function to read the Data Base
 **
 ** \return  Returned value
 **
 ** \param[in]  db        Db structure
 ** \param[in]  iatt_in   Address of the first input attribute
 ** \param[in]  iatt      Rank of the internal attribute
 ** \param[in]  itrace    Rank of the trace
 ** \param[in]  it        Rank of the sample on the trace
 **
 *****************************************************************************/
static double TR_IN(Db *db, int iatt_in, int iatt, int itrace, int it)
{
  return (db->getArray(iatt_in + iatt, NTRACE * (it) + (itrace)));
}

/****************************************************************************/
/*!
 **  Local function to write into the Data Base
 **
 ** \param[in]  db        Db structure
 ** \param[in]  iatt_out  Address of the first input attribute
 ** \param[in]  iatt      Rank of the internal attribute
 ** \param[in]  itr       Rank of the trace
 ** \param[in]  it        Rank of the sample on the trace
 ** \param[in]  value     Value to be stored
 **
 *****************************************************************************/
static void TR_OUT(Db *db,
                   int iatt_out,
                   int iatt,
                   int itr,
                   int it,
                   double value)
{
  db->setArray(iatt_out + iatt, NTRACE * (it) + (itr), value);
}

/****************************************************************************/
/*!
 **  Do unary arithmetic operation on traces
 **
 ** \return  Error return code
 **
 ** \param[in]  db        Db structure
 ** \param[in]  oper      Operator flag (::ENUM_SEISMICS)
 ** \param[in]  natt      Number of seismic attributes
 ** \param[in]  nt        Number of samples on input trace
 ** \param[in]  iatt_in   Address of the first input attribute
 ** \param[in]  iatt_out  Address of the first output attribute
 ** \param[in]  dt        time sampling interval
 **
 ** \remark Operations inv, slog and slog10 are "punctuated", meaning that if,
 ** \remark the input contains 0 values, 0 values are returned.
 **
 *****************************************************************************/
static int st_seismic_operate(Db *db,
                              int oper,
                              int natt,
                              int nt,
                              int iatt_in,
                              int iatt_out,
                              double dt)
{
  int it, iatt, itrace, count;
  double x, y, max, denom, sum, avg;

  /* Dispatch */

  for (itrace = 0; itrace < NTRACE; itrace++)
  {
    switch (oper)
    {
      case SEISMIC_NOP:
        break;

      case SEISMIC_FABS:
        for (iatt = 0; iatt < natt; iatt++)
          for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            TR_OUT(db, iatt_out, iatt, itrace, it,
                   (FFFF(x)) ? TEST :
                               ABS(TR_IN(db, iatt_in, iatt, itrace, it)));
          }
        break;

      case SEISMIC_SSQRT:
        for (iatt = 0; iatt < natt; iatt++)
          for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            TR_OUT(db, iatt_out, iatt, itrace, it,
                   FFFF(x) ? TEST :
                             SIGN(x, sqrt(ABS(x))));
          }
        break;

      case SEISMIC_SQR:
        for (iatt = 0; iatt < natt; iatt++)
          for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            TR_OUT(db, iatt_out, iatt, itrace, it, (FFFF(x)) ? TEST :
                                                               x * x);
          }
        break;

      case SEISMIC_SSQR:
        for (iatt = 0; iatt < natt; iatt++)
          for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            TR_OUT(db, iatt_out, iatt, itrace, it, (FFFF(x)) ? TEST :
                                                               SIGN(x, x * x));
          }
        break;

      case SEISMIC_SIGN:
        for (iatt = 0; iatt < natt; iatt++)
          for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            TR_OUT(db, iatt_out, iatt, itrace, it, (FFFF(x)) ? TEST :
                                                               SIGN(x, 1));
          }
        break;

      case SEISMIC_EXP:
        for (iatt = 0; iatt < natt; iatt++)
          for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            TR_OUT(db, iatt_out, iatt, itrace, it, (FFFF(x)) ? TEST :
                                                               exp(x));
          }
        break;

      case SEISMIC_SLOG:
        for (iatt = 0; iatt < natt; iatt++)
          for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            TR_OUT(db, iatt_out, iatt, itrace, it,
                   (FFFF(x) || ABS(x) <= 0) ? TEST :
                                              SIGN(x, log(ABS(x))));
          }
        break;

      case SEISMIC_SLOG10:
        for (iatt = 0; iatt < natt; iatt++)
          for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            TR_OUT(db, iatt_out, iatt, itrace, it,
                   (FFFF(x) || ABS(x) <= 0) ? TEST :
                                              SIGN(x, log10(ABS(x))));
          }
        break;

      case SEISMIC_COS:
        for (iatt = 0; iatt < natt; iatt++)
          for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            TR_OUT(db, iatt_out, iatt, itrace, it, (FFFF(x)) ? TEST :
                                                               cos(x));
          }
        break;

      case SEISMIC_SIN:
        for (iatt = 0; iatt < natt; iatt++)
          for (it = 0; it < nt; it++)
            TR_OUT(db, iatt_out, iatt, itrace, it,
                   sin(TR_IN(db, iatt_in, iatt, itrace, it)));
        break;

      case SEISMIC_TAN:
        for (iatt = 0; iatt < natt; iatt++)
          for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            TR_OUT(db, iatt_out, iatt, itrace, it, (FFFF(x)) ? TEST :
                                                               tan(x));
          }
        break;

      case SEISMIC_COSH:
        for (iatt = 0; iatt < natt; iatt++)
          for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            TR_OUT(db, iatt_out, iatt, itrace, it, (FFFF(x)) ? TEST :
                                                               cosh(x));
          }
        break;

      case SEISMIC_SINH:
        for (iatt = 0; iatt < natt; iatt++)
          for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            TR_OUT(db, iatt_out, iatt, itrace, it, (FFFF(x)) ? TEST :
                                                               sinh(x));
          }
        break;

      case SEISMIC_TANH:
        for (iatt = 0; iatt < natt; iatt++)
          for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            TR_OUT(db, iatt_out, iatt, itrace, it, (FFFF(x)) ? TEST :
                                                               tanh(x));
          }
        break;

      case SEISMIC_NORM:
        for (iatt = 0; iatt < natt; iatt++)
        {
          max = 0.0;
          for (it = 0; it < nt; it++)
          {
            x = ABS(TR_IN(db, iatt_in, iatt, itrace, it));
            if (!FFFF(x) && max < x) max = x;
          }
          if (max != 0.0) for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            TR_OUT(db, iatt_out, iatt, itrace, it, (FFFF(x)) ? TEST :
                                                               x / max);
          }
        }
        break;

      case SEISMIC_DB:
        for (iatt = 0; iatt < natt; iatt++)
          for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            TR_OUT(db, iatt_out, iatt, itrace, it,
                   (FFFF(TEST) || ABS(x) <= 0) ? TEST :
                                                 20.0 * SIGN(x, log10(ABS(x))));
          }
        break;

      case SEISMIC_NEG:
        for (iatt = 0; iatt < natt; iatt++)
          for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            TR_OUT(db, iatt_out, iatt, itrace, it, (FFFF(x)) ? TEST :
                                                               -x);
          }
        break;

      case SEISMIC_ONLY_POS:
        for (iatt = 0; iatt < natt; iatt++)
          for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            TR_OUT(db, iatt_out, iatt, itrace, it, (FFFF(x) || x <= 0) ? TEST :
                                                                         x);
          }
        break;

      case SEISMIC_ONLY_NEG:
        for (iatt = 0; iatt < natt; iatt++)
          for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            TR_OUT(db, iatt_out, iatt, itrace, it, (FFFF(x) || x >= 0) ? TEST :
                                                                         x);
          }
        break;

      case SEISMIC_SUM:
        for (iatt = 0; iatt < natt; iatt++)
        {
          for (it = 0; it < nt; it++)
            TR_OUT(db, iatt_out, iatt, itrace, it,
                   TR_IN(db, iatt_in, iatt, itrace, it));
          for (it = 1; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            y = TR_IN(db, iatt_in, iatt, itrace, it - 1);
            TR_OUT(db, iatt_out, iatt, itrace, it,
                   (FFFF(x) || FFFF(y)) ? TEST :
                                          (x + y) * (2. * dt));
          }
        }
        break;

      case SEISMIC_DIFF:
        for (iatt = 0; iatt < natt; iatt++)
        {

          /* do centered differences for the rest */
          for (it = 2; it < nt - 2; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it + 1);
            y = TR_IN(db, iatt_in, iatt, itrace, it - 1);
            TR_OUT(db, iatt_out, iatt, itrace, it,
                   (FFFF(x) || FFFF(y)) ? TEST :
                                          (x - y) / (2. * dt));
          }

          /* simple difference for tr.data[0] */
          x = TR_IN(db, iatt_in, iatt, itrace, 1);
          y = TR_IN(db, iatt_in, iatt, itrace, 0);
          TR_OUT(db, iatt_out, iatt, itrace, 0,
                 (FFFF(x) || FFFF(y)) ? TEST :
                                        (x - y) / dt);
          x = TR_IN(db, iatt_in, iatt, itrace, nt - 1);
          y = TR_IN(db, iatt_in, iatt, itrace, nt - 2);
          TR_OUT(db, iatt_out, iatt, itrace, nt - 1,
                 (FFFF(x) || FFFF(y)) ? TEST :
                                        (x - y) / dt);

          /* centered difference for tr.data[1] */
          x = TR_IN(db, iatt_in, iatt, itrace, 2);
          y = TR_IN(db, iatt_in, iatt, itrace, 0);
          TR_OUT(db, iatt_out, iatt, itrace, 1,
                 (FFFF(x) || FFFF(y)) ? TEST :
                                        (x - y) / (2. * dt));
          x = TR_IN(db, iatt_in, iatt, itrace, nt - 1);
          y = TR_IN(db, iatt_in, iatt, itrace, nt - 3);
          TR_OUT(db, iatt_out, iatt, itrace, nt - 2,
                 (FFFF(x) || FFFF(y)) ? TEST :
                                        (x - y) / (2. * dt));
        }
        break;

      case SEISMIC_REFL:
        for (iatt = 0; iatt < natt; iatt++)
        {
          for (it = nt - 1; it > 0; --it)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            y = TR_IN(db, iatt_in, iatt, itrace, it - 1);
            if (FFFF(x) || FFFF(y))
              TR_OUT(db, iatt_out, iatt, itrace, it, TEST);
            else
            {
              denom = (x + y);
              TR_OUT(db, iatt_out, iatt, itrace, it,
                     (denom == 0.) ? TEST :
                                     (x - y) / denom);
            }
          }
          TR_OUT(db, iatt_out, iatt, itrace, 0, TEST);
        }
        break;

      case SEISMIC_MOD_2PI:
        for (iatt = 0; iatt < natt; iatt++)
          for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            if (FFFF(x))
              TR_OUT(db, iatt_out, iatt, itrace, it, TEST);
            else
            {
              y = x;
              while (y < 0)
                y += 2. * GV_PI;
              while (y >= 2. * GV_PI)
                y -= 2. * GV_PI;
              TR_OUT(db, iatt_out, iatt, itrace, it, y);
            }
          }
        break;

      case SEISMIC_INV:
        for (iatt = 0; iatt < natt; iatt++)
          for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            TR_OUT(db, iatt_out, iatt, itrace, it,
                   (FFFF(x) || x == 0.) ? TEST :
                                          1. / x);
          }
        break;

      case SEISMIC_AVG:
        for (iatt = 0; iatt < natt; iatt++)
        {
          sum = avg = count = 0;
          for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            if (FFFF(x)) continue;
            count++;
            sum += x;
          }
          for (it = 0; it < nt; it++)
          {
            x = TR_IN(db, iatt_in, iatt, itrace, it);
            TR_OUT(db, iatt_out, iatt, itrace, it,
                   (FFFF(x) || count <= 0) ? TEST :
                                             x - sum / count);
          }
        }
        break;

      default:
        messerr("Non recognized operation");
        return (1);
    }
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Defining some wavelets for modeling of seimic data
 **
 ** \return  Array containing the discretized wavelet, allocated here
 ** \return  Dimension = 2 * ntw + 1
 **
 ** \param[in]  verbose   1 for a verbose output; 0 otherwise
 ** \param[in]  type      Type of the wavelet (::ENUM_WAVELETS)
 ** \param[in]  ntw       half-length of the wavelet excluding center (samples)
 ** \param[in]  tindex    time index to locate the spike (Spike)
 ** \param[in]  dt        time step
 ** \param[in]  fpeak     peak frequency of the Ricker wavelet
 ** \param[in]  period    wavelet period (s) (Ricker)
 ** \param[in]  amplitude wavelet amplitude (Ricker)
 ** \param[in]  distort   wavelet distortion factor (Ricker)
 **
 *****************************************************************************/
static double* st_seismic_wavelet(int verbose,
                                  int type,
                                  int ntw,
                                  int tindex,
                                  double dt,
                                  double fpeak,
                                  double period,
                                  double amplitude,
                                  double distort)
{
  int it, ntw2;
  double *wavelet, t, t1, t0, tnorm, value, wsym;

  /* Initializations */

  ntw2 = 2 * ntw + 1;
  t0 = ntw * dt;
  wavelet = (double*) mem_alloc(sizeof(double) * ntw2, 0);
  if (wavelet == nullptr) return (wavelet);
  for (it = 0; it < ntw2; it++)
    wavelet[it] = 0.;

  /* Dispatch */

  switch (type)
  {
    case WAVELET_NONE:
      wavelet[0] = 1.;
      break;

    case WAVELET_RICKER1:
      for (it = 0; it < ntw2; it++)
      {
        t1 = it * dt;
        value = GV_PI * fpeak * (t1 - t0);
        value = value * value;
        wavelet[it] = (1. - 2. * value) * exp(-value);
      }
      break;

    case WAVELET_RICKER2:
      t = 0.;
      for (it = 0; it <= ntw; it++, t += dt)
      {
        tnorm = t / period;
        value = (2.5 * tnorm) * (2.5 * tnorm);
        wsym = amplitude * (1. - 2. * value) * exp(-value);
        wavelet[ntw + it] = wsym * (1. - 2. * distort * tnorm);
        wavelet[ntw - it] = wsym * (1. + 2. * distort * tnorm);
      }
      break;

    case WAVELET_AKB:
      for (it = 0; it < ntw2; it++)
      {
        t1 = it * dt;
        value = fpeak * (t1 - t0);
        value = value * value;
        wavelet[it] = -amplitude * (t1 - t0) * exp(-2. * value);
      }
      break;

    case WAVELET_SPIKE:
      wavelet[tindex] = 1.0;
      break;

    case WAVELET_UNIT:
      for (it = 0; it < ntw2; it++)
        wavelet[it] = 1.0;
      break;

    default:
      messerr("This wavelet type (%d) is unknown", type);
      break;
  }

  /* Print the wavelet */

  if (verbose) print_vector("Wavelet", 0, ntw2, wavelet);

  return (wavelet);
}

/****************************************************************************/
/*!
 **  Compute the convolution of two input vector arrays
 **
 ** \param[in]  nx  length of x array
 ** \param[in]  ix0 sample index of first x
 ** \param[in]  x   array[nx] to be convolved with y
 ** \param[in]  ny  length of y array
 ** \param[in]  iy0 sample index of first y
 ** \param[in]  y   array[ny] with which x is to be convolved
 ** \param[in]  nz  length of z array
 ** \param[in]  iz0 sample index of first z
 **
 ** \param[out]  z  array[nz] containing x convolved with y
 **
 ** \remark The operation z = x convolved with y is defined to be
 ** \remark            ix0+nx-1
 ** \remark     z[i] =   sum    x[j]*y[i-j]  ;  i = iz0,...,iz0+nz-1
 ** \remark             j=ix0
 ** \remark
 ** \remark The x samples are contained in x[0], x[1], ..., x[nx-1];
 ** \remark likewise for the y and z samples.  The sample indices of the
 ** \remark first x, y, and z values determine the location of the origin
 ** \remark for each array.  For example, if z is to be a weighted average
 ** \remark of the nearest 5 samples of y, one might use
 ** \remark         x[0] = x[1] = x[2] = x[3] = x[4] = 1/5
 ** \remark         conv(5,-2,x,ny,0,y,nz,0,z)
 ** \remark In this example, the filter x is symmetric, with index of first
 ** \remark sample = -2.
 **
 *****************************************************************************/
static void st_seismic_convolve(int nx,
                                int ix0,
                                double *x,
                                int ny,
                                int iy0,
                                double *y,
                                int nz,
                                int iz0,
                                double *z)
{
  int ix1, iy1, iz1, i, j, j0, j1;
  double sum;

  ix1 = ix0 + nx - 1;
  iy1 = iy0 + ny - 1;
  iz1 = iz0 + nz - 1;
  x -= ix0;
  y -= iy0;
  z -= iz0;

  for (i = iz0; i <= iz1; ++i)
  {
    j0 = i - iy1;
    if (j0 < ix0) j0 = ix0;
    j1 = i - iy0;
    if (j1 > ix1) j1 = ix1;
    sum = 0.;
    for (j = j0; j <= j1; ++j)
      sum += x[j] * y[i - j];
    z[i] = sum;
  }
  return;
}

/****************************************************************************/
/*!
 **  Resample from depth to time
 **
 ** \return  Error return code
 **
 ** \param[in]  db_z     Depth Grid structure
 ** \param[in]  iatt_v   Address of the Velocity variable (in Depth grid)
 **
 ** \param[out]  db_t    Time Grid structure
 **
 ** \remark Linear interpolation and constant extrapolation is used to
 ** \remark determine interval velocities at times not specified.
 **
 *****************************************************************************/
int seismic_z2t_convert(DbGrid *db_z, int iatt_v, DbGrid *db_t)
{
  DbGrid *db_v;
  double *tz, *zt, *az, *at, z0, z1, t0, t1, dz, dt;
  int nz, nt, ndim, natt, iatt_t, iatt_z, error;

  /* Initializations */

  if (st_match(db_z, db_t)) return (1);
  error = 1;
  db_v = db_z;
  ndim = db_z->getNDim();
  natt = db_z->getLocNumber(ELoc::Z);
  nz = db_z->getNX(ndim - 1);
  nt = db_t->getNX(ndim - 1);
  z0 = db_z->getX0(ndim - 1);
  t0 = db_t->getX0(ndim - 1);
  dz = db_z->getDX(ndim - 1);
  dt = db_t->getDX(ndim - 1);
  t1 = t0 + (nt - 1) * dt;
  z1 = z0 + (nz - 1) * dz;
  zt = tz = az = at = nullptr;

  /* Create the output variables */

  iatt_t = db_t->addColumnsByConstant(natt, 0.);
  if (iatt_t < 0) goto label_end;
  iatt_z = db_attribute_identify(db_z, ELoc::Z, 0);

  /* Core allocation */

  zt = (double*) mem_alloc(nt * sizeof(double), 0);
  if (zt == nullptr) goto label_end;
  tz = (double*) mem_alloc(nz * sizeof(double), 0);
  if (tz == nullptr) goto label_end;
  az = (double*) mem_alloc(nz * sizeof(double), 0);
  if (az == nullptr) goto label_end;
  at = (double*) mem_alloc(nt * sizeof(double), 0);
  if (at == nullptr) goto label_end;

  /* Perform the conversion from Depth to Time */

  st_seismic_z2t_convert(db_z, iatt_z, nz, z0, z1, dz, db_t, iatt_t, nt, t0, t1,
                         dt, db_v, iatt_v, natt, tz, zt, at, az);

  /* Set the error return code */

  error = 0;

  label_end: zt = (double*) mem_free((char* ) zt);
  tz = (double*) mem_free((char* ) tz);
  az = (double*) mem_free((char* ) az);
  at = (double*) mem_free((char* ) at);
  return (error);
}

/****************************************************************************/
/*!
 **  Resample from time to depth
 **
 ** \return  Error return code
 **
 ** \param[in]  db_t     Time Grid structure
 ** \param[in]  iatt_v   Address of the Velocity variable (in Time grid)
 **
 ** \param[out]  db_z    Depth Grid structure
 **
 ** \remark Linear interpolation and constant extrapolation is used to
 ** \remark determine interval velocities at times not specified.
 **
 *****************************************************************************/
int seismic_t2z_convert(DbGrid *db_t, int iatt_v, DbGrid *db_z)
{
  DbGrid *db_v;
  double *zt, *tz, *az, *at, z0, z1, t0, t1, dz, dt;
  int nz, nt, ndim, natt, iatt_z, iatt_t, error;

  /* Initializations */

  if (st_match(db_z, db_t)) return (1);
  error = 1;
  db_v = db_t;
  ndim = db_t->getNDim();
  natt = db_t->getLocNumber(ELoc::Z);
  nz = db_z->getNX(ndim - 1);
  nt = db_t->getNX(ndim - 1);
  z0 = db_z->getX0(ndim - 1);
  t0 = db_t->getX0(ndim - 1);
  dz = db_z->getDX(ndim - 1);
  dt = db_t->getDX(ndim - 1);
  t1 = t0 + (nt - 1) * dt;
  z1 = z0 + (nz - 1) * dz;
  zt = tz = az = at = nullptr;

  /* Create the output variables */

  iatt_z = db_z->addColumnsByConstant(natt, 0.);
  if (iatt_z < 0) goto label_end;
  iatt_t = db_attribute_identify(db_t, ELoc::Z, 0);

  /* Core allocation */

  zt = (double*) mem_alloc(nt * sizeof(double), 0);
  if (zt == nullptr) goto label_end;
  tz = (double*) mem_alloc(nz * sizeof(double), 0);
  if (tz == nullptr) goto label_end;
  az = (double*) mem_alloc(nz * sizeof(double), 0);
  if (az == nullptr) goto label_end;
  at = (double*) mem_alloc(nt * sizeof(double), 0);
  if (at == nullptr) goto label_end;

  /* Perform the conversion from Time to Depth */

  st_seismic_t2z_convert(db_t, iatt_t, nt, t0, t1, dt, db_z, iatt_z, nz, z0, z1,
                         dz, db_v, iatt_v, natt, tz, zt, at, az);

  /* Set the error return code */

  error = 0;

  label_end: tz = (double*) mem_free((char* ) tz);
  zt = (double*) mem_free((char* ) zt);
  az = (double*) mem_free((char* ) az);
  at = (double*) mem_free((char* ) at);
  return (error);
}

/****************************************************************************/
/*!
 **  Do unary arithmetic operation on traces
 **
 ** \return  Error return code
 **
 ** \param[in]  db      Db structure
 ** \param[in]  oper    Operator flag (::ENUM_SEISMICS)
 **
 ** \remark Operations inv, slog and slog10 are "punctuated", meaning that if,
 ** \remark the input contains 0 values, 0 values are returned.
 **
 *****************************************************************************/
int seismic_operate(DbGrid *db, int oper)
{
  int ndim, natt, nt, iatt_in, iatt_out;
  double dt;

  /* Initializations */

  if (st_match(db, nullptr)) return (1);
  ndim = db->getNDim();
  natt = db->getLocNumber(ELoc::Z);
  nt = db->getNX(ndim - 1);
  dt = db->getDX(ndim - 1);

  /* Create the output variables */

  iatt_in = db_attribute_identify(db, ELoc::Z, 0);
  if (iatt_in < 0) return (1);
  iatt_out = db->addColumnsByConstant(natt, 0.);
  if (iatt_out < 0) return (1);

  if (st_seismic_operate(db, oper, natt, nt, iatt_in, iatt_out, dt)) return (1);

  return (0);
}

/****************************************************************************/
/*!
**  Loads the values of a Db column in the output array
**  Completes undefined values
**
** \param[in]  nz          Number of elements in the column
** \param[in]  shift       The shift of the trace before convolution
** \param[in]  val_before  Replacement value for undefined element
**                         before first defined sample
** \param[in]  val_middle  Replacement value for undefined element
**                         between defined samples
** \param[in]  val_after   Replacement value for undefined element
**                         after last defined sample
** \param[in]  tab0        Input array (Dimension: nz)
**
** \param[out] tab1        Output array (Dimension: nz)
**
*****************************************************************************/
static void st_seismic_affect(Db     * /*db*/,
                              int     nz,
                              int     shift,
                              double  val_before,
                              double  val_middle,
                              double  val_after,
                              double *tab0,
                              double *tab1)
{
  int iz, flag_already;
  double value;

  /* Set the complete column to val_before and val_after */

  for (iz = 0; iz < shift; iz++)
    tab1[iz] = val_before;
  for (iz = shift; iz < nz + 2 * shift; iz++)
    tab1[iz] = val_after;

  /* Loop on the elements of the column */

  flag_already = 0;
  for (iz = 0; iz < nz; iz++)
  {
    value = tab0[iz];
    if (FFFF(value))
    {
      value = (flag_already) ? val_middle :
                               val_before;
    }
    else
    {
      flag_already = 1;
    }
    tab1[iz + shift] = value;
  }

  return;
}

/****************************************************************************/
/*!
 **  Converts an impedance into a contrast
 **
 ** \param[in]  nz       Number of elements in the column
 **
 ** \param[out] tab      Output array (Dimension: nz)
 **
 *****************************************************************************/
static void st_seismic_contrast(int nz, double *tab)
{
  int iz;
  double denom, x, y;

  /* Loop on the elements of the column */

  for (iz = nz - 1; iz > 0; --iz)
  {
    x = tab[iz];
    y = tab[iz - 1];
    if (FFFF(x) || FFFF(y))
      tab[iz] = TEST;
    else
    {
      denom = (x + y);
      tab[iz] = (denom == 0.) ? TEST :
                                (x - y) / denom;
    }
  }
  tab[0] = TEST;
  return;
}

/****************************************************************************/
/*!
 **  Convolve with a given wavelet
 **
 ** \return  Array containing the discretized wavelet, allocated here
 **
 ** \param[in]  db            Db structure
 ** \param[in]  flag_operate  1 to perform the convolution; 0 otherwise
 ** \param[in]  flag_contrast 1 to perform contrast; 0 otherwise
 ** \param[in]  type        Type of the wavelet (::ENUM_WAVELETS)
 ** \param[in]  ntw         half-length of the wavelet excluding center (samples)
 ** \param[in]  option      option used to perform the convolution
 ** \li                     -1 : erode the edge (on ntw pixels)
 ** \li                      0 : truncate the wavelet on the edge
 ** \li                     +1 : extend the trace with padding before convolution
 ** \li                     +2 : extend the trace with the last informed values
 ** \param[in]  tindex      time index to locate the spike (Spike)
 ** \param[in]  fpeak       peak frequency of the Ricker wavelet
 ** \param[in]  period      wavelet period (s) (Ricker)
 ** \param[in]  amplitude   wavelet amplitude (Ricker)
 ** \param[in]  distort     wavelet distortion factor (Ricker)
 ** \param[in]  val_before  Replacement value for undefined element
 **                         before first defined sample
 ** \param[in]  val_middle  Replacement value for undefined element
 **                         between defined samples
 ** \param[in]  val_after   Replacement value for undefined element
 **                         after last defined sample
 ** \param[in]  wavelet     Wavelet defined as input (Dimension: 2*ntw+1)
 **
 *****************************************************************************/
int seismic_convolve(DbGrid *db,
                     int flag_operate,
                     int flag_contrast,
                     int type,
                     int ntw,
                     int option,
                     int tindex,
                     double fpeak,
                     double period,
                     double amplitude,
                     double distort,
                     double val_before,
                     double val_middle,
                     double val_after,
                     double *wavelet)
{
  int ndim, iatt, natt, itrace, iz, nz, iatt_in, iatt_out, error, size, shift;
  double *tab0, *tab1, *tab2, dz;

  /* Initializations */

  if (st_match(db, nullptr)) return (1);
  tab0 = tab1 = tab2 = nullptr;
  ndim = db->getNDim();
  natt = db->getLocNumber(ELoc::Z);
  nz = db->getNX(ndim - 1);
  dz = db->getDX(ndim - 1);
  error = size = 0;

  /* Generation of the wavelet */

  if (type >= 0)
  {
    wavelet = st_seismic_wavelet(!flag_operate, type, ntw, tindex, dz, fpeak,
                                 period, amplitude, distort);
  }
  if (wavelet == nullptr) goto label_end;
  if (!flag_operate) goto label_end;

  /* Create the output variables */

  error = 1;
  iatt_in = db_attribute_identify(db, ELoc::Z, 0);
  if (iatt_in < 0) return (1);
  iatt_out = db->addColumnsByConstant(natt, 0.);
  if (iatt_out < 0) return (1);

  /* Core allocation */

  shift = (option > 0) ? ntw :
                         0;
  size = nz + 2 * shift;
  tab0 = (double*) mem_alloc(size * sizeof(double), 0);
  if (tab0 == nullptr) goto label_end;
  tab1 = (double*) mem_alloc(size * sizeof(double), 0);
  if (tab1 == nullptr) goto label_end;
  tab2 = (double*) mem_alloc(size * sizeof(double), 0);
  if (tab2 == nullptr) goto label_end;

  /* Loop on the attributes */

  for (itrace = 0; itrace < NTRACE; itrace++)
  {
    for (iatt = 0; iatt < natt; iatt++)
    {

      for (iz = 0; iz < nz; iz++)
        tab0[iz] = TR_IN(db, iatt_in, iatt, itrace, iz);

      /* Perform the constrast (if required) */

      if (flag_contrast) st_seismic_contrast(nz, tab0);

      /* Set the padding values, if needed */

      if (option == 2)
      {
        val_before = tab0[0];
        val_after = tab0[nz - 1];
      }

      /* Load the input attribute */

      st_seismic_affect(db, nz, shift, val_before, val_middle, val_after, tab0,
                        tab1);

      /* Perform the convolution */

      st_seismic_convolve(2 * ntw + 1, -ntw, wavelet, size, 0, tab1, size, 0,
                          tab2);

      /* In case of edge erosion, set the edge to FFFF */

      if (option < 0) for (iz = 0; iz < ntw; iz++)
        tab2[iz] = tab2[nz - iz - 1] = TEST;

      /* Save the output attribute */

      for (iz = 0; iz < nz; iz++)
        TR_OUT(db, iatt_out, iatt, itrace, iz,
               (FFFF(TR_IN(db, iatt_in, iatt, itrace, iz))) ? TEST :
                                                              tab2[iz + shift]);
    }
  }

  /* Set the error return code */

  error = 0;

  label_end: tab0 = (double*) mem_free((char* ) tab0);
  tab1 = (double*) mem_free((char* ) tab1);
  tab2 = (double*) mem_free((char* ) tab2);
  if (type >= 0) wavelet = (double*) mem_free((char* ) wavelet);
  return (error);
}

/****************************************************************************/
/*!
 **  Returns the absolute index of the sample in the grid
 **
 ** \return  Returned absolute address
 **
 ** \param[in]  db        Grid Db structure
 ** \param[in]  ix        Rank of the target trace
 ** \param[in]  iz        Rank of the target sample within the target trace
 **
 *****************************************************************************/
static int st_absolute_index(DbGrid *db, int ix, int iz)
{
  int indg[3];

  indg[0] = ix;
  indg[1] = 0;
  indg[2] = iz;
  return (db_index_grid_to_sample(db, indg));
}

/****************************************************************************/
/*!
 **  Remove the sample which coincides with the target site
 **  from the neighborhood
 **
 ** \param[in,out]  ngh        ST_Seismic_Neigh structure
 **
 *****************************************************************************/
static void st_sample_remove_central(ST_Seismic_Neigh *ngh)

{
  int lec, ecr;

  ngh->n_v1 = ngh->n_v2 = 0;
  for (lec = ecr = 0; lec < ngh->nactive; lec++)
  {
    if (ngh->ix_ngh[lec] == 0 && ngh->iz_ngh[lec] == 0) continue;
    ngh->ix_ngh[ecr] = ngh->ix_ngh[lec];
    ngh->iz_ngh[ecr] = ngh->iz_ngh[lec];
    ngh->v1_ngh[ecr] = ngh->v1_ngh[lec];
    ngh->v2_ngh[ecr] = ngh->v2_ngh[lec];
    if (!FFFF(ngh->v1_ngh[lec])) ngh->n_v1++;
    if (!FFFF(ngh->v2_ngh[lec])) ngh->n_v2++;
    ecr++;
  }
  ngh->nactive = ecr;
  return;
}

/****************************************************************************/
/*!
 **  Add the sample to the neighborhood
 **
 ** \param[in]  db         Grid Db structure
 ** \param[in]  iatt_z1    Address of the first data attribute
 ** \param[in]  iatt_z2    Address of the second data attribute
 ** \param[in]  flag_test  Flag to check if the sample must be added
 ** \li                    0 : only if at least one variable is defined
 ** \li                    1 : only if the first variable is defined
 ** \li                    2 : only if the second variable is defined
 ** \param[in]  ix         Rank of the trace
 ** \param[in]  iz         Rank of the target sample within the target trace
 ** \param[in,out] ngh     ST_Seismic_Neigh structure
 **
 *****************************************************************************/
static void st_sample_add(DbGrid *db,
                          int iatt_z1,
                          int iatt_z2,
                          int flag_test,
                          int ix,
                          int iz,
                          ST_Seismic_Neigh *ngh)
{
  int iech, i, found;
  double v1, v2;

  if (ngh->nactive >= ngh->nvois) messageAbort("Overflow in st_sample_add");

  /* Check if the sample has already been allocated */

  for (i = 0, found = -1; i < ngh->nactive && found < 0; i++)
    if (ngh->ix_ngh[i] == ix && ngh->iz_ngh[i] == iz) found = i;
  if (found >= 0) return;

  /* Calculate the absolute sample address */

  iech = st_absolute_index(db, ix, iz);
  if (!db->isActive(iech)) return;

  /* Read the new values */

  v1 = db->getArray(iech, iatt_z1);
  if (flag_test == 1 && FFFF(v1)) return;
  v2 = db->getArray(iech, iatt_z2);
  if (flag_test == 2 && FFFF(v2)) return;
  if (flag_test == 0 && FFFF(v1) && FFFF(v2)) return;

  /* Store the new information in the neighborhood */

  ngh->ix_ngh[ngh->nactive] = ix;
  ngh->iz_ngh[ngh->nactive] = iz;
  ngh->v1_ngh[ngh->nactive] = v1;
  ngh->v2_ngh[ngh->nactive] = v2;
  if (!FFFF(v1)) ngh->n_v1++;
  if (!FFFF(v2)) ngh->n_v2++;
  ngh->nactive++;

  return;
}

/****************************************************************************/
/*!
 **  Check the presence of neighboring traces
 **
 ** \param[in]  db        Grid Db structure
 ** \param[in]  ivar      Rank of the variable of interest
 **
 ** \param[out] npres     Number of informed traces
 ** \param[out] presence  Array giving the number of valid samples per trace
 **
 *****************************************************************************/
static void st_estimate_check_presence(DbGrid *db,
                                       int ivar,
                                       int *npres,
                                       int *presence)
{
  int ix, iz, iech;

  /* Loop on the traces */

  *npres = 0;
  for (ix = 0; ix < NX; ix++)
  {

    /* Loop on the values along the trace */

    presence[ix] = 0;
    for (iz = 0; iz < NZ; iz++)
    {
      iech = st_absolute_index(db, ix, iz);
      if (!FFFF(db->getLocVariable(ELoc::Z,iech, ivar))) presence[ix]++;
    }
    if (presence[ix] > 0) (*npres)++;
  }

  return;
}

/****************************************************************************/
/*!
 **  Initialize the ST_Seismic_Neigh structure
 **
 ** \param[in]  ngh      ST_Seismic_Neigh structure to be freed (if mode<0)
 **
 *****************************************************************************/
static void st_estimate_neigh_init(ST_Seismic_Neigh *ngh)

{
  int i;

  ngh->nactive = 0;
  ngh->n_v1 = 0;
  ngh->n_v2 = 0;
  for (i = 0; i < ngh->nvois; i++)
  {
    ngh->ix_ngh[i] = ITEST;
    ngh->iz_ngh[i] = ITEST;
    ngh->v1_ngh[i] = TEST;
    ngh->v2_ngh[i] = TEST;
  }
}

/****************************************************************************/
/*!
 **  Manage the ST_Seismic_Neigh structure
 **
 ** \return  Address of the newly managed ST_Seismic_Neigh structure
 **
 ** \param[in]  mode     Type of management
 ** \li                  1 for allocation
 ** \li                 -1 for deallocation
 ** \param[in]  nvois    Maximum number of samples in the neighborhood
 ** \param[in]  ngh      ST_Seismic_Neigh structure to be freed (if mode<0)
 **
 *****************************************************************************/
static ST_Seismic_Neigh* st_estimate_neigh_management(int mode,
                                                      int nvois,
                                                      ST_Seismic_Neigh *ngh)
{
  /* Dispatch */

  if (mode > 0)
  {

    /* Allocation */

    ngh = (ST_Seismic_Neigh*) mem_alloc(sizeof(ST_Seismic_Neigh), 0);
    if (ngh == nullptr) return (ngh);
    ngh->nvois = nvois;
    ngh->nactive = 0;
    ngh->n_v1 = 0;
    ngh->n_v2 = 0;
    ngh->ix_ngh = nullptr;
    ngh->iz_ngh = nullptr;
    ngh->v1_ngh = nullptr;
    ngh->v2_ngh = nullptr;

    ngh->ix_ngh = (int*) mem_alloc(sizeof(int) * nvois, 0);
    if (ngh->ix_ngh == nullptr) goto label_free;
    ngh->iz_ngh = (int*) mem_alloc(sizeof(int) * nvois, 0);
    if (ngh->iz_ngh == nullptr) goto label_free;
    ngh->v1_ngh = (double*) mem_alloc(sizeof(double) * nvois, 0);
    if (ngh->v1_ngh == nullptr) goto label_free;
    ngh->v2_ngh = (double*) mem_alloc(sizeof(double) * nvois, 0);
    if (ngh->v2_ngh == nullptr) goto label_free;

    st_estimate_neigh_init(ngh);
  }
  else
  {
    label_free: if (ngh == nullptr) return (ngh);
    ngh->ix_ngh = (int*) mem_free((char* ) ngh->ix_ngh);
    ngh->iz_ngh = (int*) mem_free((char* ) ngh->iz_ngh);
    ngh->v1_ngh = (double*) mem_free((char* ) ngh->v1_ngh);
    ngh->v2_ngh = (double*) mem_free((char* ) ngh->v2_ngh);
    ngh = (ST_Seismic_Neigh*) mem_free((char* ) ngh);
  }
  return (ngh);
}

/****************************************************************************/
/*!
 **  Check if the neighborhood has changed
 **
 ** \return  1 if the Neighborhood is unchanged
 **
 ** \param[in]  ngh_old   Previous ST_Seismic_Neigh structure
 ** \param[in]  ngh_cur   Current ST_Seismic_Neigh structure
 **
 *****************************************************************************/
static int st_estimate_neigh_unchanged(ST_Seismic_Neigh *ngh_old,
                                       ST_Seismic_Neigh *ngh_cur)
{
  int i, flag_unchanged;

  flag_unchanged = 0;
  if (ngh_old == nullptr || ngh_cur == nullptr) goto label_end;
  if (ngh_old->nactive != ngh_cur->nactive) goto label_end;
  if (ngh_cur->nactive <= 0) goto label_end;

  /* Loop on the samples of the Old and the New Neighborhood */

  for (i = 0; i < ngh_old->nactive; i++)
  {
    if (ngh_cur->ix_ngh[i] != ngh_old->ix_ngh[i] || ngh_cur->iz_ngh[i]
        != ngh_old->iz_ngh[i]) goto label_end;
  }

  /* The neighborhood is unchanged */

  flag_unchanged = 1;

  label_end:

  /* Optional printout */

  if (OptDbg::query(EDbg::NBGH) && !OptDbg::force() && flag_unchanged)
    message("The neighborhood is unchanged\n");

  return (flag_unchanged);
}

/****************************************************************************/
/*!
 **  Copy the current neighborhood into the old one
 **
 ** \param[in]  ngh_cur   Current ST_Seismic_Neigh structure
 **
 ** \param[out]  ngh_old  Old ST_Seismic_Neigh structure
 **
 *****************************************************************************/
static void st_estimate_neigh_copy(ST_Seismic_Neigh *ngh_cur,
                                   ST_Seismic_Neigh *ngh_old)
{
  int i;

  /* Blank out the old structure */

  st_estimate_neigh_init(ngh_old);

  /* Copy the current neighborhood into the old one */

  for (i = 0; i < ngh_cur->nactive; i++)
  {
    ngh_old->ix_ngh[i] = ngh_cur->ix_ngh[i];
    ngh_old->iz_ngh[i] = ngh_cur->iz_ngh[i];
    ngh_old->v1_ngh[i] = ngh_cur->v1_ngh[i];
    ngh_old->v2_ngh[i] = ngh_cur->v2_ngh[i];
  }
  ngh_old->nactive = ngh_cur->nactive;
  return;
}

/****************************************************************************/
/*!
 **  Print the Kriging Information
 **
 ** \param[in]  ngh       Current ST_Seismic_Neigh structure
 ** \param[in]  ix0       Rank of the target trace
 ** \param[in]  iz0       Rank of the target sample within the target trace
 **
 *****************************************************************************/
static void st_estimate_neigh_print(ST_Seismic_Neigh *ngh, int ix0, int iz0)
{
  int i;

  /* Header */

  mestitle(0, "Neighborhood information");

  message("For (ix0=%d - iz0=%d) - Number = %d\n", ix0 + 1, iz0 + 1,
          ngh->nactive);

  if (ngh->nactive <= 0) return;
  tab_prints(NULL, "Sample");
  tab_prints(NULL, "Delta-X");
  tab_prints(NULL, "Delta-Z");
  tab_prints(NULL, "V1");
  tab_prints(NULL, "V2");
  message("\n");
  for (i = 0; i < ngh->nactive; i++)
  {
    tab_printi(NULL, i + 1);
    tab_printi(NULL, ngh->ix_ngh[i]);
    tab_printi(NULL, ngh->iz_ngh[i]);
    tab_printg(NULL, ngh->v1_ngh[i]);
    tab_printg(NULL, ngh->v2_ngh[i]);
    message("\n");
  }
}

/****************************************************************************/
/*!
 **  Establish the neighborhood
 **
 ** \return  Error return code
 **
 ** \param[in]  db        Grid Db structure
 ** \param[in]  flag_exc  1 if the target sample must be excluded
 ** \param[in]  iatt_z1   Address of the first data attribute
 ** \param[in]  iatt_z2   Address of the second data attribute
 ** \param[in]  ix0       Rank of the target trace
 ** \param[in]  iz0       Rank of the target sample within the target trace
 ** \param[in]  nbench    Vertical Radius of the neighborhood (center excluded)
 ** \param[in]  nv2max    Maximum number of traces of second variable on each
 **                       side of the target trace
 ** \param[in]  presence  Array giving the number of valid samples per trace
 **
 ** \param[out]  ngh      Current ST_Seismic_Neigh structure
 **
 *****************************************************************************/
static int st_estimate_neigh_create(DbGrid *db,
                                    int flag_exc,
                                    int iatt_z1,
                                    int iatt_z2,
                                    int ix0,
                                    int iz0,
                                    int nbench,
                                    int nv2max,
                                    int /*npres*/[2],
                                    int *presence[2],
                                    ST_Seismic_Neigh *ngh)
{
  int i, idx, ix, iz, jz, count, flag_valid;

  /* Blank out the ST_Seismic_Neigh structure */

  st_estimate_neigh_init(ngh);

  /* Load the primary variable first (and its collocated variable) */

  for (ix = 0; ix < NX; ix++)
  {
    if (presence[0][ix] <= 0) continue;

    for (iz = -nbench; iz <= nbench; iz++)
    {
      jz = iz + iz0;
      if (jz < 0 || jz >= NZ) continue;
      st_sample_add(db, iatt_z1, iatt_z2, 1, ix, jz, ngh);
    }
  }

  /* Load the information from the current trace (if any) */

  for (iz = -nbench; iz <= nbench; iz++)
  {
    jz = iz + iz0;
    if (jz < 0 || jz >= NZ) continue;
    st_sample_add(db, iatt_z1, iatt_z2, 0, ix0, jz, ngh);
  }

  /* Load the information from the closest traces (for V2) */
  /* located before the current trace */

  count = 0;
  for (idx = 1; idx < NX; idx++)
  {
    ix = ix0 - idx;
    if (ix < 0 || ix >= NX) continue;
    if (presence[1][ix] <= 0) continue;
    count++;
    if (count > nv2max) break;

    /* Add the samples of the trace of interest */

    for (iz = -nbench; iz <= nbench; iz++)
    {
      jz = iz + iz0;
      if (jz < 0 || jz >= NZ) continue;
      st_sample_add(db, iatt_z1, iatt_z2, 0, ix, jz, ngh);
    }
  }

  /* Load the information from the closest traces (for V2) */
  /* located before the current trace */

  count = 0;
  for (idx = 1; idx < NX; idx++)
  {
    ix = ix0 + idx;
    if (ix < 0 || ix >= NX) continue;
    if (presence[1][ix] <= 0) continue;
    count++;
    if (count > nv2max) break;

    /* Add the samples of the trace of interest */

    for (iz = -nbench; iz <= nbench; iz++)
    {
      jz = iz + iz0;
      if (jz < 0 || jz >= NZ) continue;
      st_sample_add(db, iatt_z1, iatt_z2, 0, ix, jz, ngh);
    }
  }

  /* Convert location indices into index distances from the target */

  for (i = 0; i < ngh->nactive; i++)
  {
    ngh->ix_ngh[i] -= ix0;
    ngh->iz_ngh[i] -= iz0;
  }

  /* Optional printout */

  if (OptDbg::query(EDbg::NBGH)) st_estimate_neigh_print(ngh, ix0, iz0);

  /* Check the validity of the neighborhood */

  flag_valid = (ngh->n_v1 >= 1 && ngh->n_v2 >= 1);

  /* Exclude the target site (simulation case) */

  if (flag_exc) st_sample_remove_central(ngh);

  /* In case of error, clean the neighborhood */

  if (!flag_valid) st_estimate_neigh_init(ngh);
  return (!flag_valid);
}

/****************************************************************************/
/*!
 **  Establish the Kriging flag
 **
 ** \param[in]  ngh       Current ST_Seismic_Neigh structure
 ** \param[in]  nfeq      Number of drift condition per variable
 **
 ** \param[out] flag      Array containing the Kriging flag
 ** \param[out] nred      Reduced number of equations
 **
 *****************************************************************************/
static void st_estimate_flag(ST_Seismic_Neigh *ngh,
                             int nfeq,
                             int *flag,
                             int *nred)
{
  double value;
  int i, ivar, nech, neqmax, nvalid;

  /* Initializations */

  nech = ngh->nactive;
  neqmax = NVAR * (nech + nfeq);

  /* Initialize the flag array */

  for (i = 0; i < neqmax; i++)
    flag[i] = 1;

  /* Cancel the flag array */

  for (ivar = 0; ivar < NVAR; ivar++)
    for (i = 0; i < nech; i++)
    {
      value = (ivar == 0) ? ngh->v1_ngh[i] :
                            ngh->v2_ngh[i];
      if (FFFF(value)) flag[ivar * nech + i] = 0;
    }

  /* Cancel the drift */

  for (ivar = 0; ivar < NVAR; ivar++)
  {
    /* Count the number of valid information per variable */

    nvalid = 0;
    for (i = 0; i < nech; i++)
      nvalid += flag[ivar * nech + i];
    if (nvalid <= 0) flag[NVAR * nech + ivar] = 0;
  }

  /* Count the number of valid equations */

  *nred = 0;
  for (i = 0; i < neqmax; i++)
    (*nred) += flag[i];

  return;
}

/****************************************************************************/
/*!
 **  Establish the constant terms for the variance calculations
 **
 ** \param[in]  model     Model structure
 **
 ** \param[out] var0      Array containing the C00 terms (Dimension: NVAR)
 **
 *****************************************************************************/
static void st_estimate_var0(Model *model, double *var0)
{
  CovCalcMode mode;

  VectorDouble d1(model->getDimensionNumber());
  mode.setMember(ECalcMember::VAR);
  model_calcul_cov(NULL,model, mode, 1, 1., d1, covtab);

  for (int ivar = 0; ivar < NVAR; ivar++)
    var0[ivar] = COVTAB(ivar, ivar);

  return;
}

/****************************************************************************/
/*!
 **  Establish the constant terms for the variance calculations
 **
 ** \param[in]  model    Model structure
 **
 ** \param[out] c00      Array containing the C00 terms (Dimension: NVAR * NVAR)
 **
 *****************************************************************************/
static void st_estimate_c00(Model *model, double *c00)
{
  CovCalcMode mode;
  VectorDouble d1(model->getDimensionNumber());
  mode.setMember(ECalcMember::VAR);
  model_calcul_cov(NULL,model, mode, 1, 1., d1, covtab);

  for (int ivar = 0; ivar < NVAR; ivar++)
    for (int jvar = 0; jvar < NVAR; jvar++)
      C00(ivar,jvar) = COVTAB(ivar, jvar);

  return;
}

/****************************************************************************/
/*!
 **  Establish the Kriging L.H.S.
 **
 ** \param[in]  ngh       Current ST_Seismic_Neigh structure
 ** \param[in]  model     Model structure
 ** \param[in]  nfeq      Number of drift condition per variable
 ** \param[in]  nred      Reduced number of equations
 ** \param[in]  flag      Array giving the flag
 **
 ** \param[out] lhs       Array containing the Kriging L.H.S.
 **
 *****************************************************************************/
static void st_estimate_lhs(ST_Seismic_Neigh *ngh,
                            Model *model,
                            int nfeq,
                            int nred,
                            int *flag,
                            double *lhs)
{
  int iech, jech, ivar, jvar, nech, neqmax, i, j, lec, ecr;
  VectorDouble d1;
  CovCalcMode mode;

  /* Initializations */

  d1.resize(3, 0.);
  nech = ngh->nactive;
  neqmax = NVAR * (nech + nfeq);

  /* Initialize the L.H.S. array */

  for (i = 0; i < neqmax * neqmax; i++)
    lhs[i] = 0.;

  /* Establish the covariance part */

  for (iech = 0; iech < nech; iech++)
    for (jech = 0; jech < nech; jech++)
    {
      d1[0] = DX * (ngh->ix_ngh[iech] - ngh->ix_ngh[jech]);
      d1[2] = DZ * (ngh->iz_ngh[iech] - ngh->iz_ngh[jech]);
      model_calcul_cov(NULL,model, mode, 1, 1., d1, covtab);

      for (ivar = 0; ivar < NVAR; ivar++)
        for (jvar = 0; jvar < NVAR; jvar++)
          LHS(iech,ivar,jech,jvar) = COVTAB(ivar, jvar);
    }

  /* Establish the drift part */

  if (nfeq)
  {
    for (iech = 0; iech < nech; iech++)
      for (ivar = 0; ivar < NVAR; ivar++)
        for (jvar = 0; jvar < NVAR; jvar++)
        {
          LHS(iech,ivar,jvar,NVAR) = (ivar == jvar);
          LHS(jvar,NVAR,iech,ivar) = (ivar == jvar);
        }
  }

  /* Compress from isotopic to heterotopic */

  lec = ecr = 0;
  for (i = 0; i < neqmax; i++)
    for (j = 0; j < neqmax; j++, lec++)
      if (flag[i] && flag[j]) lhs[ecr++] = lhs[lec];

  if (OptDbg::query(EDbg::KRIGING)) krige_lhs_print(nech, neqmax, nred, flag, lhs);

  return;
}

/****************************************************************************/
/*!
 **  Establish the Kriging R.H.S.
 **
 ** \param[in]  ngh       Current ST_Seismic_Neigh structure
 ** \param[in]  model     Model structure
 ** \param[in]  nfeq      Number of drift condition per variable
 ** \param[in]  nred      Reduced number of equations
 ** \param[in]  flag      Array giving the flag
 **
 ** \param[out] rhs       Array containing the Kriging R.H.S.
 **
 *****************************************************************************/
static void st_estimate_rhs(ST_Seismic_Neigh *ngh,
                            Model *model,
                            int nfeq,
                            int nred,
                            int *flag,
                            double *rhs)
{
  int iech, nech, ivar, jvar, neqmax, lec, ecr, i;
  VectorDouble d1;
  CovCalcMode mode;

  /* Initializations */

  d1.resize(3, 0.);
  nech = ngh->nactive;
  neqmax = NVAR * (nech + nfeq);
  mode.setMember(ECalcMember::RHS);

  /* Initialize the L.H.S. array */

  for (i = 0; i < neqmax * NVAR; i++)
    rhs[i] = 0.;

  /* Establish the covariance part */

  for (iech = 0; iech < nech; iech++)
  {
    d1[0] = DX * ngh->ix_ngh[iech];
    d1[2] = DZ * ngh->iz_ngh[iech];
    model_calcul_cov(NULL,model, mode, 1, 1., d1, covtab);

    for (ivar = 0; ivar < NVAR; ivar++)
      for (jvar = 0; jvar < NVAR; jvar++)
        RHS(iech,ivar,jvar) = COVTAB(ivar, jvar);
  }

  /* Establish the drift part */

  if (nfeq)
  {
    for (ivar = 0; ivar < NVAR; ivar++)
      for (jvar = 0; jvar < NVAR; jvar++)
        RHS(jvar,NVAR,ivar) = (ivar == jvar);
  }

  /* Compress from isotopic to heterotopic */

  ecr = lec = 0;
  for (ivar = 0; ivar < NVAR; ivar++)
    for (i = 0; i < neqmax; i++, lec++)
      if (flag[i]) rhs[ecr++] = rhs[lec];

  if (OptDbg::query(EDbg::KRIGING))
    krige_rhs_print(NVAR, nech, neqmax, nred, flag, rhs);

  return;
}

/****************************************************************************/
/*!
 **  Print the kriging weights
 **
 ** \param[in]  ngh     ST_Seismic_Neigh structure
 ** \param[in]  nvar    Number of variables
 ** \param[in]  nech    Number of active points
 ** \param[in]  nred    Reduced number of equations
 ** \param[in]  flag    Flag array
 ** \param[in]  wgt     Array of Kriging weights
 **
 *****************************************************************************/
static void st_wgt_print(ST_Seismic_Neigh *ngh,
                         int nvar,
                         int nech,
                         int nred,
                         int *flag,
                         double *wgt)
{
  double sum[2], value;
  int iwgt, ivar, jvar, iech, lec, cumflag;
  char string[10];

  /* Header */

  mestitle(0, "(Co-) Kriging weights");

  /* First line */

  tab_prints(NULL, "Rank");
  tab_prints(NULL, "Delta-X");
  tab_prints(NULL, "Delta-Z");
  tab_prints(NULL, "Data");
  for (ivar = 0; ivar < nvar; ivar++)
  {
    (void) gslSPrintf(string, "Z%d*", ivar + 1);
    tab_prints(NULL, string);
  }
  message("\n");

  /* Display the information and the weights */

  for (jvar = lec = cumflag = 0; jvar < nvar; jvar++)
  {
    if (nvar > 1) message("Using variable Z%-2d\n", jvar + 1);

    /* Loop on the samples */

    for (ivar = 0; ivar < nvar; ivar++)
      sum[ivar] = 0.;
    for (iech = 0; iech < nech; iech++, lec++)
    {
      tab_printi(NULL, iech + 1);
      tab_printi(NULL, ngh->ix_ngh[iech]);
      tab_printi(NULL, ngh->iz_ngh[iech]);
      if (jvar == 0)
        tab_printg(NULL, ngh->v1_ngh[iech]);
      else
        tab_printg(NULL, ngh->v2_ngh[iech]);

      for (ivar = 0; ivar < nvar; ivar++)
      {
        iwgt = nred * ivar + cumflag;
        value = (flag[lec]) ? wgt[iwgt] :
                              TEST;
        if (!FFFF(value)) sum[ivar] += value;
        tab_printg(NULL, value);
      }
      if (flag[lec]) cumflag++;
      message("\n");
    }

    tab_prints(NULL, "Sum of weights", 4, EJustify::LEFT);
    for (ivar = 0; ivar < nvar; ivar++)
      tab_printg(NULL, sum[ivar]);
    message("\n");
  }

  return;
}

/****************************************************************************/
/*!
**  Establish the Kriging Weights
**
** \return  Error return code
**
** \param[in]  ngh       Current ST_Seismic_Neigh structure
** \param[in]  nred      Reduced number of Kriging equations
** \param[in]  flag      Array giving the flag
** \param[in]  lhs       Array of Kriging L.H.S.
** \param[in]  rhs       Array of Kriging R.H.S.
**
** \param[out] wgt       Array containing the Kriging Weights
**
*****************************************************************************/
static int st_estimate_wgt(ST_Seismic_Neigh *ngh,
                           Model  * /*model*/,
                           int     nred,
                           int    *flag,
                           double *lhs,
                           double *rhs,
                           double *wgt)
{
  int nech;

  /* Initializations */

  nech = ngh->nactive;
  if (nech <= 0) return (0);

  /* Calculate the kriging weights */

  if (matrix_invert(lhs, nred, IECH_OUT)) return (1);
  matrix_product(nred, nred, 2, lhs, rhs, wgt);

  if (OptDbg::query(EDbg::KRIGING)) st_wgt_print(ngh, NVAR, nech, nred, flag, wgt);

  return (0);
}

/****************************************************************************/
/*!
**  Perform the Kriging
**
** \param[in]  db        Grid Db structure
** \param[in]  ngh       Current ST_Seismic_Neigh structure
** \param[in]  flag_std  1 for the calculation of the St. Dev.
** \param[in]  nred      Reduced number of equations
** \param[in]  flag      Array giving the flag
** \param[in]  wgt       Array containing the Kriging weights
** \param[in]  rhs       Array containing the R.H.S. member
** \param[in]  var0      Array containing the C00 terms
** \param[in]  iatt_est  Array of pointers to the kriged result
** \param[in]  iatt_std  Array of pointers to the variance of estimation
**
*****************************************************************************/
static void st_estimate_result(Db     *db,
                               ST_Seismic_Neigh *ngh,
                               int     flag_std,
                               int     /*nfeq*/,
                               int     nred,
                               int    *flag,
                               double *wgt,
                               double *rhs,
                               double *var0,
                               int *iatt_est,
                               int *iatt_std)
{
  int i, ivar, jvar, nech, lec;
  double result, value, stdev;

  /* Initializations */

  nech = ngh->nactive;
  if (OptDbg::query(EDbg::RESULTS)) mestitle(0, "(Co-) Kriging results");

  /* Loop on the variables */

  for (ivar = 0; ivar < NVAR; ivar++)
  {

    /* Estimation */

    result = stdev = 0.;
    lec = ivar * nred;
    for (jvar = 0; jvar < NVAR; jvar++)
      for (i = 0; i < nech; i++)
      {
        if (!flag[jvar * nech + i]) continue;
        value = (jvar == 0) ? ngh->v1_ngh[i] :
                              ngh->v2_ngh[i];
        result += wgt[lec++] * value;
      }
    db->setArray(IECH_OUT, iatt_est[ivar], result);

    /* Variance of estimation */

    if (flag_std)
    {
      stdev = var0[ivar];
      lec = ivar * nred;
      for (i = 0; i < nred; i++)
        stdev -= rhs[lec + i] * wgt[lec + i];
      stdev = (stdev > 0) ? sqrt(stdev) :
                            0.;
      db->setArray(IECH_OUT, iatt_std[ivar], stdev);
    }

    if (OptDbg::query(EDbg::RESULTS))
    {
      tab_printi(NULL, ivar + 1);
      tab_printg(" - Estimate  = ", result);
      if (flag_std) tab_printg(" - St. Dev.  = ", stdev);
      message("\n");
    }
  }
  return;
}

/****************************************************************************/
/*!
**  Perform the Simulation
**
** \param[in]  db        Grid Db structure
** \param[in]  ix0       Rank of the target trace
** \param[in]  iz0       Rank of the target sample within the target trace
** \param[in]  ngh       Current ST_Seismic_Neigh structure
** \param[in]  nbsimu    Number of simulations
** \param[in]  nred      Reduced number of equations
** \param[in]  flag      Array giving the flag
** \param[in]  wgt       Array containing the Kriging weights
** \param[in]  rhs       Array containing the R.H.S. member
** \param[in]  c00       Array containing the C00 terms
** \param[in]  iatt_sim  Array of pointers to the simulated results
**
*****************************************************************************/
static void st_simulate_result(DbGrid *db,
                               int     ix0,
                               int     iz0,
                               ST_Seismic_Neigh *ngh,
                               int     nbsimu,
                               int     /*nfeq*/,
                               int     nred,
                               int    *flag,
                               double *wgt,
                               double *rhs,
                               double *c00,
                               int *iatt_sim)
{
  int i, ivar, jvar, nech, lec, isimu, iech;
  double result[2], value, sigma0, sigma1, sigma2, z1, z2, l0, lb[4];

  /* Initializations */

  nech = ngh->nactive;
  if (OptDbg::query(EDbg::RESULTS)) mestitle(0, "(Co-) Simulation results");

  /* Pre-calculations */

  for (ivar = 0; ivar < NVAR; ivar++)
    for (jvar = 0; jvar < NVAR; jvar++)
    {
      value = 0.;
      for (i = 0; i < nred; i++)
        value += wgt[nred * ivar + i] * rhs[nred * jvar + i];
      LB(ivar,jvar) = value;
    }

  /* Loop on the simulations */

  for (isimu = 0; isimu < nbsimu; isimu++)
  {

    /* Estimation of the variables (processed individually) */

    for (ivar = 0; ivar < NVAR; ivar++)
    {

      /* Estimation */

      result[ivar] = 0.;
      lec = ivar * nred;
      for (jvar = 0; jvar < NVAR; jvar++)
        for (i = 0; i < nech; i++)
        {
          if (!flag[jvar * nech + i]) continue;
          iech = st_absolute_index(db, ix0 + ngh->ix_ngh[i],
                                   iz0 + ngh->iz_ngh[i]);
          value = db->getArray(iech, iatt_sim[jvar] + isimu);
          result[ivar] += wgt[lec++] * value;
        }
    }

    /* Link the two variables */

    sigma1 = C00(0,0) - LB(0, 0);
    sigma2 = C00(1,1) - LB(1, 1);

    z2 = db->getArray(IECH_OUT, iatt_sim[1] + isimu);
    if (FFFF(z2)) z2 = result[1] + sqrt(MAX(0, sigma2)) * law_gaussian();
    z1 = db->getArray(IECH_OUT, iatt_sim[0] + isimu);
    if (FFFF(z1) && sigma2 > 0)
    {
      l0 = (C00(0,1) - LB(0, 1)) / sigma2;
      sigma0 = sigma1 - l0 * (C00(0,1) - LB(1, 0));
      sigma0 = sqrt(MAX(0, sigma0));
      z1 = result[0] + (z2 - result[1]) * l0 + sigma0 * law_gaussian();
    }
    result[0] = z1;
    result[1] = z2;

    /* Store the results */

    for (ivar = 0; ivar < NVAR; ivar++)
    {
      db->setArray(IECH_OUT, iatt_sim[ivar] + isimu, result[ivar]);

      if (OptDbg::query(EDbg::RESULTS))
      {
        message("Simulation #%d of Z%-2d : ", isimu + 1, ivar + 1);
        tab_printg(" = ", result[ivar]);
        message("\n");
      }
    }
  }

  return;
}

/****************************************************************************/
/*!
 **  Order the traces to be processed
 **
 ** \return Error return code
 **
 ** \param[in] presence  Array giving the number of valid samples per trace
 **
 ** \param[out] rank     Array giving the order of the traces
 **
 *****************************************************************************/
static int st_estimate_sort(int *presence, int *rank)
{
  double *dist, distmin, distval;
  int ix, jx;

  /* Initializations */

  dist = nullptr;

  /* Core allocation */

  dist = (double*) mem_alloc(sizeof(double) * NX, 0);
  if (dist == nullptr) return (1);

  /* Sort the traces */
  for (ix = 0; ix < NX; ix++)
  {
    rank[ix] = ix;
    distmin = 1.e30;
    for (jx = 0; jx < NX; jx++)
    {
      if (presence[jx] <= 0) continue;
      distval = ABS(ix - jx);
      if (distval <= distmin) distmin = distval;
    }
    dist[ix] = distmin;
  }
  ut_sort_double(1, NX, rank, dist);

  /* Core deallocation */

  dist = (double*) mem_free((char* ) dist);
  return (0);
}

/****************************************************************************/
/*!
 **  Perform a bivariate estimation on a grid
 **
 ** \return  Error return code
 **
 ** \param[in,out]  db     Grid Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  nbench     Vertical Radius of the neighborhood (center excluded)
 ** \param[in]  nv2max     Maximum number of traces of second variable on each
 **                        side of the target trace
 ** \param[in]  flag_ks    1 for a Simple Kriging; otherwise Ordinary Kriging
 ** \param[in]  flag_std   1 for the calculation of the St. Dev.
 ** \param[in]  flag_sort  1 if the traces to be treated are sorted
 **                        by increasing distance to trace where first
 **                        variable is defined
 ** \param[in]  flag_stat  1 for producing final statistics
 **
 *****************************************************************************/
int seismic_estimate_XZ(DbGrid *db,
                        Model *model,
                        int nbench,
                        int nv2max,
                        int flag_ks,
                        int flag_std,
                        int flag_sort,
                        int flag_stat)
{
  int *flag, *rank, *presence[2], npres[2], iatt_est[2], iatt_std[2];
  int i, ix0, jx0, iz0, nvois, size, error, nred, nfeq, iatt_z1, iatt_z2;
  int nb_total, nb_process, nb_calcul;
  double *lhs, *rhs, *wgt, *var0;
  ST_Seismic_Neigh *ngh_cur, *ngh_old;

  /* Initializations */

  error = 1;
  nvois = size = nred = nb_total = nb_process = nb_calcul = 0;
  lhs = rhs = wgt = var0 = covtab = nullptr;
  flag = rank = nullptr;
  nfeq = (flag_ks) ? 0 :
                     1;
  ngh_cur = ngh_old = nullptr;
  for (i = 0; i < 2; i++)
  {
    iatt_est[i] = -1;
    iatt_std[i] = -1;
    presence[i] = nullptr;
    npres[i] = 0;
  }
  if (krige_koption_manage(1, 1, EKrigOpt::PONCTUAL, 1, VectorInt()))
    goto label_end;

  /* Check that the grid is XZ */

  if (! db->isGrid())
  {
    messerr("The Db structure must be a Grid Db");
    return (1);
  }
  NX = db->getNX(0);
  NY = db->getNX(1);
  NZ = db->getNX(2);
  DX = db->getDX(0);
  DZ = db->getDX(2);
  NVAR = db->getLocNumber(ELoc::Z);
  iatt_z1 = db_attribute_identify(db, ELoc::Z, 0);
  iatt_z2 = db_attribute_identify(db, ELoc::Z, 1);

  if (!(NX > 1 && NY == 1 && NZ > 1))
  {
    messerr("The Db grid is not XZ");
    goto label_end;
  }
  if (NVAR != 2)
  {
    messerr("This function is restricted to the use of two variable");
    goto label_end;
  }

  /* Core allocation (first) */

  for (i = 0; i < 2; i++)
  {
    presence[i] = (int*) mem_alloc(sizeof(int) * NX, 0);
    if (presence[i] == nullptr) goto label_end;
  }

  /* Look for columns where each variable is defined */

  for (i = 0; i < 2; i++)
    st_estimate_check_presence(db, i, &npres[i], presence[i]);

  /* Maximum dimension of the neighborhood */

  nvois = (2 * nbench + 1) * (npres[0] + 2 * nv2max + 1);
  size = NVAR * (nvois + nfeq);

  /* Core allocation (second) */

  ngh_cur = st_estimate_neigh_management(1, nvois, ngh_cur);
  if (ngh_cur == nullptr) goto label_end;
  ngh_old = st_estimate_neigh_management(1, nvois, ngh_old);
  if (ngh_old == nullptr) goto label_end;
  covtab = (double*) mem_alloc(sizeof(double) * NVAR * NVAR, 0);
  if (covtab == nullptr) goto label_end;
  lhs = (double*) mem_alloc(sizeof(double) * size * size, 0);
  if (lhs == nullptr) goto label_end;
  rhs = (double*) mem_alloc(sizeof(double) * size * NVAR, 0);
  if (rhs == nullptr) goto label_end;
  wgt = (double*) mem_alloc(sizeof(double) * size * NVAR, 0);
  if (wgt == nullptr) goto label_end;
  flag = (int*) mem_alloc(sizeof(int) * size, 0);
  if (flag == nullptr) goto label_end;
  rank = (int*) mem_alloc(sizeof(int) * NX, 0);
  if (rank == nullptr) goto label_end;
  if (flag_std)
  {
    var0 = (double*) mem_alloc(sizeof(double) * NVAR, 0);
    if (var0 == nullptr) goto label_end;
  }

  /* Add the resulting variables */

  for (i = 0; i < 2; i++)
  {
    iatt_est[i] = db->addColumnsByConstant(1, TEST);
    if (iatt_est[i] < 0) goto label_end;
  }
  if (flag_std) for (i = 0; i < 2; i++)
  {
    iatt_std[i] = db->addColumnsByConstant(1, TEST);
    if (iatt_std[i] < 0) goto label_end;
  }

  /* Calculate the constant terms for the variances */

  if (flag_std) st_estimate_var0(model, var0);

  /* Calculate the order of the columns to be treated */

  if (st_estimate_sort(presence[0], rank)) goto label_end;

  /* Loop on the grid nodes */

  for (jx0 = 0; jx0 < NX; jx0++)
  {
    ix0 = (flag_sort) ? rank[jx0] :
                        jx0;
    for (iz0 = 0; iz0 < NZ; iz0++)
    {
      nb_total++;
      IECH_OUT = st_absolute_index(db, ix0, iz0);
      OptDbg::setCurrentIndex(IECH_OUT + 1);
      if (!db->isActive(IECH_OUT)) continue;

      /* Look for the neighborhood */

      if (st_estimate_neigh_create(db, 0, iatt_z1, iatt_z2, ix0, iz0, nbench,
                                   nv2max, npres, presence, ngh_cur)) continue;

      /* Check if the neighborhood is unchanged */

      if (!st_estimate_neigh_unchanged(ngh_old, ngh_cur) || OptDbg::force())
      {
        nb_calcul++;

        /* Establish the flag */

        st_estimate_flag(ngh_cur, nfeq, flag, &nred);

        /* Establish the kriging L.H.S. */

        st_estimate_lhs(ngh_cur, model, nfeq, nred, flag, lhs);

        /* Establish the kriging R.H.S. */

        st_estimate_rhs(ngh_cur, model, nfeq, nred, flag, rhs);

        /* Derive the kriging weights */

        if (st_estimate_wgt(ngh_cur, model, nred, flag, lhs, rhs, wgt))
          continue;
      }

      /* Perform the estimation */

      st_estimate_result(db, ngh_cur, flag_std, nfeq, nred, flag, wgt, rhs,
                         var0, iatt_est, iatt_std);

      /* Save the neighborhood */

      st_estimate_neigh_copy(ngh_cur, ngh_old);
      nb_process++;
    }
  }

  /* Set the error return code */

  error = 0;

  label_end: OptDbg::setCurrentIndex(0);
  if (flag_stat)
  {
    message("Statistics on the number of nodes:\n");
    message("- Total number of grid nodes  = %d\n", nb_total);
    message("- Nodes sucessfully processed = %d\n", nb_process);
    message("- Kriging system calculated   = %d\n", nb_calcul);
  }
  for (i = 0; i < 2; i++)
  {
    presence[i] = (int*) mem_free((char* ) presence[i]);
    if (error && iatt_est[i] >= 0) db->deleteColumnByUID(iatt_est[i]);
    if (error && iatt_std[i] >= 0) db->deleteColumnByUID(iatt_std[i]);
  }
  (void) krige_koption_manage(-1, 1, EKrigOpt::PONCTUAL, 1, VectorInt());
  flag = (int*) mem_free((char* ) flag);
  lhs = (double*) mem_free((char* ) lhs);
  rhs = (double*) mem_free((char* ) rhs);
  wgt = (double*) mem_free((char* ) wgt);
  var0 = (double*) mem_free((char* ) var0);
  covtab = (double*) mem_free((char* ) covtab);
  ngh_cur = st_estimate_neigh_management(-1, nvois, ngh_cur);
  ngh_old = st_estimate_neigh_management(-1, nvois, ngh_old);

  return (error);
}

/****************************************************************************/
/*!
 **  Copy the main attributes in the simulations
 **
 ** \param[in]  db       Grid Db structure
 ** \param[in]  nbsimu   Number of simulations
 ** \param[in]  iatt     Address of the first item for each variable
 **
 *****************************************************************************/
static void st_copy_attribute(Db *db, int nbsimu, int *iatt)
{
  int ivar, isimu, iech, nech;
  double value;

  /* Initializations */

  nech = db->getSampleNumber();

  /* Loop on the variables */

  for (ivar = 0; ivar < NVAR; ivar++)
  {

    /* Loop on the samples */

    for (iech = 0; iech < nech; iech++)
    {
      value = db->getLocVariable(ELoc::Z,iech, ivar);

      /* Loop on the simulations */

      for (isimu = 0; isimu < nbsimu; isimu++)
        db->setArray(iech, iatt[ivar] + isimu, value);
    }
  }
  return;
}

/****************************************************************************/
/*!
 **  Perform a bivariate cosimulation on a grid
 **
 ** \return  Error return code
 **
 ** \param[in,out] db      Grid Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  nbench     Vertical Radius of the neighborhood (center excluded)
 ** \param[in]  nv2max     Maximum number of traces of second variable on each
 **                        side of the target trace
 ** \param[in]  nbsimu     Number of simulations
 ** \param[in]  seed       Seed for the random number generator
 ** \param[in]  flag_ks    1 for a Simple Kriging; otherwise Ordinary Kriging
 ** \param[in]  flag_sort  1 if the traces to be treated are sorted
 **                        by increasing distance to trace where first
 **                        variable is defined
 ** \param[in]  flag_stat  1 for producing final statistics
 **
 *****************************************************************************/
int seismic_simulate_XZ(DbGrid *db,
                        Model *model,
                        int nbench,
                        int nv2max,
                        int nbsimu,
                        int seed,
                        int flag_ks,
                        int flag_sort,
                        int flag_stat)
{
  int *flag, *rank, *presence[2], npres[2], iatt_sim[2];
  int i, isimu, ix0, iz0, nvois, size, error, nred, nfeq, jx0;
  int nb_total, nb_process, nb_calcul;
  double *lhs, *rhs, *wgt, *c00;
  ST_Seismic_Neigh *ngh_cur, *ngh_old;

  /* Initializations */

  error = 1;
  nvois = size = nred = nb_total = nb_process = nb_calcul = 0;
  lhs = rhs = wgt = c00 = covtab = nullptr;
  flag = rank = nullptr;
  nfeq = (flag_ks) ? 0 :
                     1;
  ngh_cur = ngh_old = nullptr;
  for (i = 0; i < 2; i++)
  {
    iatt_sim[i] = -1;
    presence[i] = nullptr;
    npres[i] = 0;
  }

  /* Check that the grid is XZ */

  if (! db->isGrid())
  {
    messerr("The Db structure must be a Grid Db");
    return (1);
  }
  NX = db->getNX(0);
  NY = db->getNX(1);
  NZ = db->getNX(2);
  DX = db->getDX(0);
  DZ = db->getDX(2);
  NVAR = db->getLocNumber(ELoc::Z);

  if (!(NX > 1 && NY == 1 && NZ > 1))
  {
    messerr("The Db grid is not XZ");
    goto label_end;
  }
  if (NVAR != 2)
  {
    messerr("This function is restricted to the use of two variable");
    goto label_end;
  }

  /* Core allocation (first) */

  for (i = 0; i < 2; i++)
  {
    presence[i] = (int*) mem_alloc(sizeof(int) * NX, 0);
    if (presence[i] == nullptr) goto label_end;
  }

  /* Look for columns where each variable is defined */

  law_set_random_seed(seed);
  for (i = 0; i < 2; i++)
    st_estimate_check_presence(db, i, &npres[i], presence[i]);

  /* Maximum dimension of the neighborhood */

  nvois = (2 * nbench + 1) * (npres[0] + 2 * nv2max + 1);
  size = NVAR * (nvois + nfeq);

  /* Core allocation (second) */

  ngh_cur = st_estimate_neigh_management(1, nvois, ngh_cur);
  if (ngh_cur == nullptr) goto label_end;
  ngh_old = st_estimate_neigh_management(1, nvois, ngh_old);
  if (ngh_old == nullptr) goto label_end;
  covtab = (double*) mem_alloc(sizeof(double) * NVAR * NVAR, 0);
  if (covtab == nullptr) goto label_end;
  lhs = (double*) mem_alloc(sizeof(double) * size * size, 0);
  if (lhs == nullptr) goto label_end;
  rhs = (double*) mem_alloc(sizeof(double) * size * NVAR, 0);
  if (rhs == nullptr) goto label_end;
  wgt = (double*) mem_alloc(sizeof(double) * size * NVAR, 0);
  if (wgt == nullptr) goto label_end;
  flag = (int*) mem_alloc(sizeof(int) * size, 0);
  if (flag == nullptr) goto label_end;
  c00 = (double*) mem_alloc(sizeof(double) * NVAR * NVAR, 0);
  if (c00 == nullptr) goto label_end;
  rank = (int*) mem_alloc(sizeof(int) * NX, 0);
  if (rank == nullptr) goto label_end;

  /* Add the resulting variables */

  for (i = 0; i < NVAR; i++)
  {
    iatt_sim[i] = db->addColumnsByConstant(nbsimu, TEST);
    if (iatt_sim[i] < 0) goto label_end;
  }

  /* Copy the known values in the simulation arrays */

  st_copy_attribute(db, nbsimu, iatt_sim);

  /* Calculate the constant terms for the variances */

  st_estimate_c00(model, c00);

  /* Calculate the order of the columns to be treated */

  if (st_estimate_sort(presence[0], rank)) goto label_end;

  /* Loop on the grid nodes */

  for (jx0 = 0; jx0 < NX; jx0++)
  {
    ix0 = (flag_sort) ? rank[jx0] :
                        jx0;
    for (iz0 = 0; iz0 < NZ; iz0++)
    {
      nb_total++;
      IECH_OUT = st_absolute_index(db, ix0, iz0);
      OptDbg::setCurrentIndex(IECH_OUT + 1);
      if (!db->isActive(IECH_OUT)) continue;

      /* Look for the neighborhood */

      if (st_estimate_neigh_create(db, 1, iatt_sim[0], iatt_sim[1], ix0, iz0,
                                   nbench, nv2max, npres, presence, ngh_cur))
        continue;

      /* Check if the neighborhood is unchanged */

      if (!st_estimate_neigh_unchanged(ngh_old, ngh_cur) || OptDbg::force())
      {
        nb_calcul++;

        /* Establish the flag */

        st_estimate_flag(ngh_cur, nfeq, flag, &nred);

        /* Establish the kriging L.H.S. */

        st_estimate_lhs(ngh_cur, model, nfeq, nred, flag, lhs);

        /* Establish the kriging R.H.S. */

        st_estimate_rhs(ngh_cur, model, nfeq, nred, flag, rhs);

        /* Derive the kriging weights */

        if (st_estimate_wgt(ngh_cur, model, nred, flag, lhs, rhs, wgt))
          continue;
      }

      /* Perform the estimation */

      st_simulate_result(db, ix0, iz0, ngh_cur, nbsimu, nfeq, nred, flag, wgt,
                         rhs, c00, iatt_sim);

      /* Save the neighborhood */

      st_estimate_neigh_copy(ngh_cur, ngh_old);
      nb_process++;
    }
    presence[1][ix0] = 1;
  }

  /* Set the error return code */

  error = 0;

  label_end: OptDbg::setCurrentIndex(0);
  if (flag_stat)
  {
    message("Statistics on the number of nodes:\n");
    message("- Total number of grid nodes  = %d\n", nb_total);
    message("- Nodes sucessfully processed = %d\n", nb_process);
    message("- Kriging system calculated   = %d\n", nb_calcul);
  }
  for (i = 0; i < 2; i++)
  {
    presence[i] = (int*) mem_free((char* ) presence[i]);
    if (error) for (isimu = 0; isimu < nbsimu; isimu++)
      db->deleteColumnByUID(iatt_sim[i] + isimu);
  }
  flag = (int*) mem_free((char* ) flag);
  rank = (int*) mem_free((char* ) rank);
  lhs = (double*) mem_free((char* ) lhs);
  rhs = (double*) mem_free((char* ) rhs);
  wgt = (double*) mem_free((char* ) wgt);
  c00 = (double*) mem_free((char* ) c00);
  covtab = (double*) mem_free((char* ) covtab);
  ngh_cur = st_estimate_neigh_management(-1, nvois, ngh_cur);
  ngh_old = st_estimate_neigh_management(-1, nvois, ngh_old);

  return (error);
}
