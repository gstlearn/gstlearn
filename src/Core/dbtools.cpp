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
#include "Basic/Grid.hpp"
#include "Basic/Law.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/PolyLine2D.hpp"
#include "Basic/String.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Core/Keypair.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbHelper.hpp"
#include "LinearOp/OptimCostColored.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "LinearOp/ShiftOpMatrix.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Model/CovInternal.hpp"
#include "Model/Model.hpp"
#include "Polygon/Polygons.hpp"
#include "Stats/Classical.hpp"
#include "Tree/Ball.hpp"
#include "geoslib_define.h"
#include "geoslib_old_f.h"
#include <cmath>
#include <cstring>

// https://stackoverflow.com/a/26359433/3952924
#ifdef _MSC_VER
#  define strncasecmp _strnicmp
#  define strcasecmp  _stricmp
#  define strncasecmp _strnicmp
#  define strcasecmp  _stricmp
#endif

/*! \cond */
#define TRACE(i, iseg)        (trace[(i) * nseg + (iseg)])
#define LINE(nbline, i)       (line[npline * (nbline) + (i)])
#define PROP1(iz, iprop)      (prop1[(iz) * nprop + (iprop)])
#define PROP2(iz, iprop)      (prop2[(iz) * nprop + (iprop)])
#define WTAB(iz, icode, ivar) (wtab[(ivar) + nvar * ((iz) + nz * (icode))])
#define WCOR(iz, icode, idim) (wcor[(idim) + ndim * ((iz) + nz * (icode))])
#define WCNT(iz, icode)       (wcnt[(iz) + nz * (icode)])
#define TRACE(i, iseg)        (trace[(i) * nseg + (iseg)])
#define LINE(nbline, i)       (line[npline * (nbline) + (i)])
#define PROP1(iz, iprop)      (prop1[(iz) * nprop + (iprop)])
#define PROP2(iz, iprop)      (prop2[(iz) * nprop + (iprop)])
#define WTAB(iz, icode, ivar) (wtab[(ivar) + nvar * ((iz) + nz * (icode))])
#define WCOR(iz, icode, idim) (wcor[(idim) + ndim * ((iz) + nz * (icode))])
#define WCNT(iz, icode)       (wcnt[(iz) + nz * (icode)])

#define R(i, j) (R[(i) * n + (j)])
#define R(i, j) (R[(i) * n + (j)])

namespace gstlrn
{
/*! \endcond */

/****************************************************************************/
/*!
 **  Calculate the (discretized) surface of influence
 **
 ** \return  Error returned code
 **
 ** \param[in]  db_point Db containing the data points
 ** \param[in]  db_grid  Db containing the discretization grid
 ** \param[in]  dlim     Maximum distance (TEST if not defined)
 **
 ** \param[out]  dtab    Array containing the surface of influence
 **                      (Dimension = Number of samples in db_point)
 ** \param[out]  gtab    Array containing the surface of influence of the
 **                      polygon to which it belongs (or TEST)
 **                      (Dimension = Number of samples in db_grid)
 **
 *****************************************************************************/
Id surface(Db* db_point,
           DbGrid* db_grid,
           Id /*icol*/,
           double dlim,
           double* dtab,
           double* gtab)
{
  bool flagTest;

  if (!db_grid->hasSameDimension(db_point)) return (1);
  Id ndim = db_point->getNDim();

  /* Preliminary calculations */

  double d2max  = (FFFF(dlim)) ? MAXIMUM_BIG : dlim * dlim;
  double maille = db_grid->getCellSize();
  for (Id iech = 0; iech < db_point->getNSample(); iech++)
    dtab[iech] = 0.;

  /* Loop on the target points */

  VectorDouble vgrid(ndim);
  for (Id igrid = 0; igrid < db_grid->getNSample(); igrid++)
  {
    gtab[igrid] = -1;
    if (!db_grid->isActive(igrid)) continue;
    flagTest = false;
    for (Id idim = 0; idim < ndim && !flagTest; idim++)
    {
      vgrid[idim] = db_grid->getCoordinate(igrid, idim);
      if (FFFF(vgrid[idim])) flagTest = true;
    }
    if (flagTest) continue;

    /* Loop on the data points */

    double d2min = d2max;
    for (Id iech = 0; iech < db_point->getNSample(); iech++)
    {
      if (!db_point->isActive(iech)) continue;

      /* Calculate the distance between node and data */

      double dist = 0.;
      flagTest    = false;
      for (Id idim = 0; idim < ndim && !flagTest; idim++)
      {
        double v2 = db_point->getCoordinate(iech, idim);
        if (FFFF(v2))
          flagTest = true;
        else
        {
          double delta = v2 - vgrid[idim];
          dist += delta * delta;
        }
      }
      if (flagTest) continue;
      if (dist > d2max) continue;

      /* Keep the closest sample */

      if (dist > d2min) continue;
      gtab[igrid] = iech;
      d2min       = dist;
    }
  }

  /* Calculate the influence of each datum */

  for (Id igrid = 0; igrid < db_grid->getNSample(); igrid++)
  {
    Id jech = static_cast<Id>(gtab[igrid]);
    if (jech >= 0) dtab[jech]++;
  }
  for (Id iech = 0; iech < db_point->getNSample(); iech++)
    dtab[iech] *= maille;

  /* Evaluate each grid node with the size of the influence polygon */
  /* to which it belongs                                            */

  for (Id igrid = 0; igrid < db_grid->getNSample(); igrid++)
  {
    Id jech = static_cast<Id>(gtab[igrid]);
    if (jech >= 0)
      gtab[igrid] = dtab[jech];
    else
      gtab[igrid] = TEST;
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Generates the discretized points along the trace
 **
 ** \param[in]  nseg   Number of vertices along the trace
 ** \param[in]  trace  Array defining the trace
 **                    (Dimension: 2 * nseg)
 ** \param[in]  disc   Discretization distance
 **
 ** \param[out] np_arg   Number of discretized points
 ** \param[out] xp       Array of first coordinates
 ** \param[out] yp       Array of second coordinates
 ** \param[out] dd       Array of distances between discretized points
 ** \param[out] del      Array of distances between vertices
 ** \param[out] dist_arg Total distance of the trace
 **
 *****************************************************************************/
void ut_trace_discretize(Id nseg,
                         const double* trace,
                         double disc,
                         Id* np_arg,
                         VectorDouble& xp,
                         VectorDouble& yp,
                         VectorDouble& dd,
                         VectorDouble& del,
                         double* dist_arg)
{
  double deltax, deltay, x0, y0, x1, y1, dist;
  Id iseg, np, ecr, nloc, ip;

  /* Initializations */

  (*np_arg) = np = 0;
  (*dist_arg) = x1 = y1 = 0.;
  del.resize(nseg);
  del[0] = 0.;

  /* Loop on the trace segments */

  for (iseg = ecr = 0; iseg < nseg - 1; iseg++)
  {

    /* Consider a segment trace */

    x0     = TRACE(0, iseg);
    y0     = TRACE(1, iseg);
    x1     = TRACE(0, iseg + 1);
    y1     = TRACE(1, iseg + 1);
    deltax = x1 - x0;
    deltay = y1 - y0;
    dist   = sqrt(deltax * deltax + deltay * deltay);
    (*dist_arg) += dist;
    del[iseg + 1] = (*dist_arg);

    /* Discretize the trace segment */

    nloc = static_cast<Id>(floor(dist / disc));
    if (ABS(nloc * disc - dist) < dist / 1000) nloc--;
    np += nloc;
    xp.resize(np);
    yp.resize(np);

    for (ip = 0; ip < nloc; ip++, ecr++)
    {
      xp[ecr] = x0 + deltax * ip / nloc;
      yp[ecr] = y0 + deltay * ip / nloc;
    }
  }

  /* Adding the last vertex */

  np++;
  xp.resize(np);
  yp.resize(np);
  xp[ecr] = x1;
  yp[ecr] = y1;
  ecr++;

  /* Elaborate the vector of distances */

  dd.resize(np);
  dd[0] = 0.;
  for (ip = 0; ip < np - 1; ip++)
  {
    deltax     = xp[ip + 1] - xp[ip];
    deltay     = yp[ip + 1] - yp[ip];
    dd[ip + 1] = dd[ip] + sqrt(deltax * deltax + deltay * deltay);
  }
}

/*****************************************************************************/
/*!
 **  Sample the point Db close to discretized points of the trace
 **
 ** \param[in]  db     Db to be sampled
 ** \param[in]  ptype  Type of locator
 ** \param[in]  np     Number of discretized points
 ** \param[in]  xp     Array of first coordinates
 ** \param[in]  yp     Array of second coordinates
 ** \param[in]  dd     Array of distances
 ** \param[in]  radius Neighborhood radius
 **
 ** \param[out] ns_arg     Number of sampled points
 ** \param[out] xs         Array of first coordinates of sampled points
 ** \param[out] ys         Array of second coordinates of sampled points
 ** \param[out] rks        Array of sample indices (starting from 1)
 ** \param[out] lys        Array of layer indices of sampled points
 ** \param[out] typ        Array of data type
 **                        1 for hard data in Z or TIME
 **                        2 for lower bound
 **                        3 for upper bound
 **
 *****************************************************************************/
void ut_trace_sample(Db* db,
                     const ELoc& ptype,
                     Id np,
                     const double* xp,
                     const double* yp,
                     const double* dd,
                     double radius,
                     Id* ns_arg,
                     VectorDouble& xs,
                     VectorDouble& ys,
                     VectorInt& rks,
                     VectorInt& lys,
                     VectorInt& typ)
{
  Id iech, ip, ns, ipmin, nvar;
  double cote, layer, bound[2];
  double radcarre, xx, yy, delx, dely, dist, ddmin;

  /* Initializations */

  radcarre = radius * radius;
  ns       = 0;
  nvar     = db->getNInterval();

  /* Loop on the samples */

  for (iech = 0; iech < db->getNSample(); iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Coordinates of the sample point */

    xx = db->getCoordinate(iech, 0);
    yy = db->getCoordinate(iech, 1);

    /* Loop on the discretized samples */

    ipmin = -1;
    ddmin = MAXIMUM_BIG;
    for (ip = 0; ip < np; ip++)
    {
      delx = xx - xp[ip];
      dely = yy - yp[ip];
      dist = delx * delx + dely * dely;
      if (dist > radcarre || dist > ddmin) continue;
      ddmin = dist;
      ipmin = ip;
    }
    if (ipmin < 0) continue;

    /* Keep sample defined by locator */

    cote = db->getFromLocator(ptype, iech);
    if (!FFFF(cote))
    {
      layer = db->getFromLocator(ELoc::LAYER, iech);
      xs.resize(ns + 1);
      ys.resize(ns + 1);
      lys.resize(ns + 1);
      typ.resize(ns + 1);
      rks.resize(ns + 1);
      xs[ns]  = dd[ipmin];
      ys[ns]  = cote;
      lys[ns] = (FFFF(layer)) ? 1 : static_cast<Id>(layer) + 1;
      typ[ns] = 1;
      rks[ns] = iech + 1;
      ns++;
    }

    /* Keep sample defined by locator UP or LOW */

    for (Id ivar = 0; ivar < nvar; ivar++)
    {
      bound[0] = db->getLocVariable(ELoc::L, iech, ivar);
      bound[1] = db->getLocVariable(ELoc::U, iech, ivar);
      for (Id ib = 0; ib < 2; ib++)
      {
        if (FFFF(bound[ib])) continue;
        xs.resize(ns + 1);
        ys.resize(ns + 1);
        lys.resize(ns + 1);
        typ.resize(ns + 1);
        rks.resize(ns + 1);
        xs[ns]  = dd[ipmin];
        ys[ns]  = bound[ib];
        lys[ns] = ivar + 1;
        typ[ns] = ib + 2;
        rks[ns] = iech + 1;
        ns++;
      }
    }
  }

  /* Returning arguments */

  *ns_arg = ns;
}

/*****************************************************************************/
/*!
 **  Create a set of samples according to a Poisson process
 **
 ** \return  Array of returned values (or NULL)
 **
 ** \param[in]  number      Number of samples to be generated
 ** \param[in]  coormin     Vector of minimum coordinates
 ** \param[in]  coormax     Vector of maximum coordinates
 ** \param[in]  flag_repulsion True if the repulsion process is active
 ** \param[in]  range       Repulsion range
 ** \param[in]  beta        Repulsion beta coefficient
 **
 *****************************************************************************/
static VectorDouble st_point_init_homogeneous(Id number,
                                              const VectorDouble& coormin,
                                              const VectorDouble& coormax,
                                              bool flag_repulsion,
                                              double range,
                                              double beta)
{
  VectorDouble tab;

  if (coormin.empty() || coormax.empty())
  {
    messerr("This method requires 'coormin' and 'coormax' defined");
    return tab;
  }
  VectorDouble extend = VH::subtract(coormin, coormax);
  Id ndim             = static_cast<Id>(coormin.size());
  VectorDouble coor(ndim);
  VectorDouble delta(ndim);

  // Generate the samples

  tab.resize(ndim * number);

  Id ecr = 0;
  for (Id ip = 0; ip < number; ip++)
  {
    for (Id idim = 0; idim < ndim; idim++)
      coor[idim] = coormin[idim] + law_uniform(0., extend[idim]);

    // Check if the point is acceptable

    bool flag_drop = false;
    if (flag_repulsion)
    {

      // Calculate the shortest distance with the previous samples

      double ddmin = MAXIMUM_BIG;
      for (Id jp = 0; jp < ip; jp++)
      {
        for (Id idim = 0; idim < ndim; idim++)
          delta[idim] = (tab[ndim * jp + idim] - coor[idim]) / range;
        double dd = VH::norm(delta);
        if (dd < ddmin) ddmin = dd;
      }

      /* Check the rejection criterion */

      double proba = exp(-pow(ddmin, beta));
      double alea  = law_uniform(0., 1.);
      flag_drop    = (alea < proba);
    }

    // Add the new point
    if (flag_drop) continue;
    for (Id idim = 0; idim < ndim; idim++)
      tab[ndim * ecr + idim] = coor[idim];
    ecr++;
  }
  tab.resize(ndim * ecr);
  return (tab);
}

/*****************************************************************************/
/*!
 **  Create a set of samples according to a Poisson Regionalized process
 **
 ** \return  Array of returned values (or NULL)
 **
 ** \param[in]  number      Number of samples to be generated
 ** \param[in]  dbgrid      Descriptor of the Db grid parameters
 ** \param[in]  flag_repulsion True if the repulsion process is active
 ** \param[in]  range       Repulsion range
 ** \param[in]  beta        Repulsion beta coefficient
 **
 ** \remarks Thinning can only be defined in 2-D.
 ** \remarks If the thinning is regionalized, its parameters are stored
 ** \remarks as NOSTAT variables: Range-1, Range-2 and Angle
 **
 *****************************************************************************/
static VectorDouble st_point_init_inhomogeneous(Id number,
                                                DbGrid* dbgrid,
                                                bool flag_repulsion,
                                                double range,
                                                double beta)
{
  VectorDouble tab;

  Id ndim = dbgrid->getNDim();
  if (dbgrid == nullptr)
  {
    messerr("This function requires a DbGrid data base");
    return tab;
  }
  if (!dbgrid->isGrid())
  {
    messerr("This function requires the Db organized as a grid");
    return tab;
  }
  bool flag_dens   = (dbgrid->getNLoc(ELoc::Z) == 1);
  bool flag_region = (ndim == 2 && dbgrid->getNLoc(ELoc::NOSTAT) == (ndim + 1));

  VectorDouble coor(ndim);
  VectorDouble coorbis(ndim);
  VectorDouble delta(ndim);
  VectorDouble radius(ndim);
  VectorDouble radip(ndim);
  double angip = 0.;

  /* Evaluate the density */

  Id ngrid = dbgrid->getNSample(true);
  VectorDouble dens;
  dens.resize(ngrid, 0.);
  double denstot = 0.;
  if (flag_dens)
  {
    Id ig = 0;
    for (Id jg = 0, ng = dbgrid->getNSampleActive(); jg < ng; jg++)
    {
      if (!dbgrid->isActiveAndDefined(jg, 0)) continue;
      double densloc = dbgrid->getZVariable(jg, 0);
      if (densloc >= 0) denstot += densloc;
      dens[ig++] = denstot;
    }
  }
  else
  {
    denstot = dbgrid->getNSample(true);
  }

  /* Point generation */

  Id ecr    = 0;
  Id indip  = 0;
  Id indjp  = 0;
  Id ntrial = 0;
  while (number - ecr > ntrial / 10)
  {
    // Draw a probability

    double proba = law_uniform(0., denstot);
    ntrial++;

    // Draw a grid cell at random

    if (flag_dens)
    {
      double denscum = 0.;
      indip          = -1;
      for (Id ig = 0; ig < ngrid && indip < 0; ig++)
      {
        if (!dbgrid->isActive(ig)) continue;
        double densloc = dbgrid->getZVariable(ig, 0);
        if (FFFF(densloc) || densloc < 0) continue;
        denscum += densloc;
        if (denscum > proba) indip = ig;
      }
      if (indip < 0) indip = ngrid - 1;
    }
    else
    {
      indip = static_cast<Id>(proba);
    }

    // Draw the point within the elected cell

    dbgrid->rankToCoordinatesInPlace(indip, coor);
    for (Id idim = 0; idim < ndim; idim++)
      coor[idim] += law_uniform(0., dbgrid->getDX(idim));
    if (flag_region)
    {
      for (Id idim = 0; idim < ndim; idim++)
        radip[idim] = dbgrid->getFromLocator(ELoc::NOSTAT, indip, idim);
      angip = dbgrid->getFromLocator(ELoc::NOSTAT, indip, ndim);
    }

    // Check if the point is acceptable

    bool flag_drop = false;
    if (flag_repulsion)
    {

      // Calculate the shortest distance with the previous samples

      flag_drop = false;
      for (Id jp = 0; jp < ecr && !flag_drop; jp++)
      {
        double dd = 0.;
        for (Id idim = 0; idim < ndim; idim++)
        {
          coorbis[idim] = tab[ndim * jp + idim];
          delta[idim]   = (coorbis[idim] - coor[idim]);
        }

        if (!flag_region)
        {
          dd = VH::norm(delta) / range;
        }
        else
        {
          indjp = dbgrid->coordinateToRank(coorbis);
          for (Id idim = 0; idim < ndim; idim++)
            radius[idim] = 2. / (radip[idim] + dbgrid->getFromLocator(ELoc::NOSTAT, indjp, idim));
          double angle = (angip + dbgrid->getFromLocator(ELoc::NOSTAT, indjp, ndim)) / 2.;
          Tensor tensor(ndim);
          tensor.setRotationAngle(0, angle);
          tensor.setRadiusVec(radius);
          dd = VH::norm(tensor.applyInverse(delta));
        }

        // Check if the point 'ip' must be dropped
        proba     = exp(-pow(dd, beta));
        flag_drop = (law_uniform(0., 1.) < proba);
      }
    }
    if (flag_drop) continue;

    // Add the new point
    for (Id idim = 0; idim < ndim; idim++)
      tab.push_back(coor[idim]);
    ecr++;
  }

  return tab;
}

/****************************************************************************/
/*!
 **  Create indicator residual variables
 **
 ** \return  Error returned code
 **
 ** \param[in]  db      Db structure
 ** \param[in]  ivar    Index of the target variable
 ** \param[in]  zcut    Array containing the cutoffs
 **
 ** \remarks The array 'zcut' must be provided in increasing order
 **
 *****************************************************************************/
Id db_resind(Db* db, Id ivar, const VectorDouble& zcut)
{
  Id nech = db->getNSample();
  Id ncut = static_cast<Id>(zcut.size());
  if (!VH::isSorted(zcut, true))
  {
    messerr("The cutoffs must be provided in increasing order");
    return 1;
  }

  /* Calculate the tonnages */

  Id ntot = 0;
  VectorDouble tonnage(ncut, 0);
  for (Id iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    double value = db->getArray(iech, ivar);
    if (FFFF(value)) continue;
    ntot++;

    for (Id icut = 0; icut < ncut; icut++)
      if (value >= zcut[icut]) tonnage[icut]++;
  }
  for (Id icut = 0; icut < ncut; icut++)
    tonnage[icut] /= static_cast<double>(ntot);

  /* Create the variables */

  Id iptr = db->addColumnsByConstant(ncut, TEST);
  if (iptr < 0) return 1;

  /* Loop on the samples */

  for (Id iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    double value = db->getArray(iech, ivar);
    if (FFFF(value)) continue;

    /* Loop on the cutoffs */

    for (Id icut = 0; icut < ncut; icut++)
    {
      double zval     = zcut[icut];
      Id ind_cut0     = (value > zval);
      zval            = (icut > 0) ? zcut[icut - 1] : 0.;
      Id ind_cut1     = (value > zval);
      double ton_cut0 = tonnage[icut];
      double ton_cut1 = (icut > 0) ? tonnage[icut - 1] : 1.;
      Id ir           = ind_cut0 / ton_cut0 - ind_cut1 / ton_cut1;
      db->setArray(iech, iptr + icut, ir);
    }
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Normalize the gradient components
 **
 ** \param[in]  dbgrid  Db structure (grid organized)
 **
 *****************************************************************************/
static void st_gradient_normalize(Db* dbgrid)

{
  double norme, grad;
  Id ndim;

  /* Initializations */

  ndim = dbgrid->getNDim();

  /* Loop on the samples */

  for (Id iech = 0; iech < dbgrid->getNSample(); iech++)
  {

    norme = 0.;
    for (Id idim = 0; idim < ndim; idim++)
    {
      grad = dbgrid->getLocVariable(ELoc::G, iech, idim);
      norme += grad * grad;
    }

    if (norme <= 0) continue;
    norme = sqrt(norme);

    for (Id idim = 0; idim < ndim; idim++)
    {
      grad = dbgrid->getLocVariable(ELoc::G, iech, idim);
      dbgrid->setLocVariable(ELoc::G, iech, idim, grad / norme);
    }
  }
}

/****************************************************************************/
/*!
 **  Calculate the gradient over a grid
 **
 ** \return  Rank of the newly created variables (or -1)
 **
 ** \param[in]  dbgrid  Db structure (grid organized)
 **
 *****************************************************************************/
Id db_gradient_components(DbGrid* dbgrid)

{
  Id iptrz, iptr, nx, ny, nz, nmax, error, ndim, j1, j2, number;
  double dinc, v1, v2, delta;
  VectorInt indg;

  /* Preliminary check */

  error = number = 1;
  iptr           = -1;
  ndim           = dbgrid->getNDim();
  if (!dbgrid->isGrid())
  {
    messerr("The Db should be organized as a Grid");
    goto label_end;
  }
  if (!dbgrid->isNVarComparedTo(1)) goto label_end;
  if (ndim > 3)
  {
    messerr("This function is limited to Space Dimension <= 3");
    goto label_end;
  }

  /* Initializations */

  nx = dbgrid->getNX(0);
  ny = dbgrid->getNX(1);
  nz = dbgrid->getNX(2);
  indg.resize(ndim, 0);

  /* Create the new variable */

  iptrz = dbgrid->getColIdxByLocator(ELoc::Z, 0);
  if (iptrz < 0) goto label_end;
  iptr = dbgrid->addColumnsByConstant(ndim, TEST, String(), ELoc::G);

  /* Calculate the Gradient components */

  for (Id ix = 0; ix < nx; ix++)
    for (Id iy = 0; iy < ny; iy++)
      for (Id iz = 0; iz < nz; iz++)
      {
        for (Id idim = 0; idim < ndim; idim++)
        {
          nmax = dbgrid->getNX(idim);
          dinc = dbgrid->getDX(idim);

          v1 = v2 = 0.;
          if (idim == 0)
          {
            j1 = (ix + 1 > nmax - 1) ? ix : ix + 1;
            v1 = get_grid_value(dbgrid, iptrz, indg, j1, iy, iz);
            if (FFFF(v1)) continue;
            j2 = (ix - 1 < 0) ? ix : ix - 1;
            v2 = get_grid_value(dbgrid, iptrz, indg, j2, iy, iz);
            if (FFFF(v2)) continue;
            number = j1 - j2;
          }

          if (idim == 1)
          {
            j1 = (iy + 1 > nmax - 1) ? iy : iy + 1;
            v1 = get_grid_value(dbgrid, iptrz, indg, ix, j1, iz);
            if (FFFF(v1)) continue;
            j2 = (iy - 1 < 0) ? iy : iy - 1;
            v2 = get_grid_value(dbgrid, iptrz, indg, ix, j2, iz);
            if (FFFF(v2)) continue;
            number = j1 - j2;
          }

          if (idim == 2)
          {
            j1 = (iz + 1 > nmax - 1) ? iz : iz + 1;
            v1 = get_grid_value(dbgrid, iptrz, indg, ix, iy, j1);
            if (FFFF(v1)) continue;
            j2 = (iz - 1 < 0) ? iz : iz - 1;
            v2 = get_grid_value(dbgrid, iptrz, indg, ix, iy, j2);
            if (FFFF(v2)) continue;
            number = j1 - j2;
          }

          delta = (v1 - v2) / (number * dinc);
          set_grid_value(dbgrid, iptr + idim, indg, ix, iy, iz, delta);
        }
      }

  /* Set the error return code */

  error = 0;

label_end:
  if (error)
  {
    dbgrid->deleteColumnsByUIDRange(iptr, ndim);
    iptr = -1;
  }
  return (iptr);
}

/****************************************************************************/
/*!
 **  Check if one (at least) of the gradient components is undefined
 **
 ** \return  1 If one component (at least) is undefined
 **
 ** \param[in]  dbgrid    Db structure (grid organized)
 ** \param[in]  iptr_grad Rank of the first gradient component
 ** \param[in]  iech      Sample rank
 **
 *****************************************************************************/
static Id st_is_undefined(Db* dbgrid, Id iptr_grad, Id iech)
{
  Id ndim;

  ndim = dbgrid->getNDim();
  for (Id idim = 0; idim < ndim; idim++)
  {
    if (FFFF(dbgrid->getArray(iech, iptr_grad + idim))) return (1);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Check if the gradient is null
 **
 ** \return  1 If gradient is null
 **
 ** \param[in]  dbgrid    Db structure (grid organized)
 ** \param[in]  iptr_grad Rank of the first gradient component
 ** \param[in]  iech      Sample rank
 **
 *****************************************************************************/
static Id st_is_zero(Db* dbgrid, Id iptr_grad, Id iech)
{
  double grad, delta;
  Id ndim;

  grad = 0.;
  ndim = dbgrid->getNDim();
  for (Id idim = 0; idim < ndim; idim++)
  {
    delta = dbgrid->getArray(iech, iptr_grad + idim);
    grad += delta * delta;
  }
  return (grad < 1.e-5);
}

/****************************************************************************/
/*!
 **  Get the next gradient-based data
 **
 ** \return  1 Error return code
 **
 ** \param[in]  dbgrid    Db structure (grid organized)
 ** \param[in]  iptr_grad Rank of the first gradient component
 ** \param[in]  coor      Array of coordinates
 **
 ** \param[out] knd       Absolute index
 ** \param[out] surf      Local value of the surface
 **
 *****************************************************************************/
static Id st_get_next(DbGrid* dbgrid,
                      Id iptr_grad,
                      VectorDouble& coor,
                      Id* knd,
                      double* surf)
{
  Id knd_loc;
  double surf_loc;

  knd_loc = dbgrid->coordinateToRank(coor);
  if (knd_loc < 0) return 1;
  if (!dbgrid->isActive(knd_loc)) return 1;
  surf_loc = dbgrid->getZVariable(knd_loc, 0);
  if (FFFF(surf_loc) || st_is_undefined(dbgrid, iptr_grad, knd_loc)) return (1);
  if (st_is_zero(dbgrid, iptr_grad, knd_loc)) return (1);
  *knd  = knd_loc;
  *surf = surf_loc;
  return 0;
}

/****************************************************************************/
/*!
 **  Calculate the streamlines
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid    Db structure (grid organized)
 ** \param[in]  dbpoint   Db structure for control points
 ** \param[in]  niter     Maximum number of iterations
 ** \param[in]  step      Progress step
 ** \param[in]  flag_norm 1 if the gradients must be normalized
 ** \param[in]  use_grad  1 if the existing gradients must be used
 **                       0 the gradients must be calculated here
 ** \param[in]  save_grad 1 save the gradients generated in this function
 **                       0 delete gradients when calculated here
 **
 ** \param[out] nbline_loc Number of streamline steps
 ** \param[out] npline_loc Number of values per line vertex
 ** \param[out] line       Array of streamline steps (Dimension: 5 * nbline)
 **
 ** \remarks The returned array 'line_loc' must be freed by the calling function
 ** \remarks Use get_keypone("Streamline_Skip",1) to define the skipping ratio
 **
 *****************************************************************************/
Id db_streamline(DbGrid* dbgrid,
                 Db* dbpoint,
                 Id niter,
                 double step,
                 Id flag_norm,
                 Id use_grad,
                 Id save_grad,
                 Id* nbline_loc,
                 Id* npline_loc,
                 VectorDouble& line)
{
  Id error, npline, idim, ecr;
  Id iptr_time, iptr_accu, iptr_grad, nbline, knd, nquant, nbyech, ndim;
  double surf, date;
  static Id quant = 1000;
  VectorDouble coor;
  VectorDouble coor0;

  /* Initializations */

  error  = 1;
  nbline = nquant = 0;
  iptr_grad       = -1;
  if (dbpoint == nullptr) dbpoint = dbgrid;
  nbyech = static_cast<Id>(get_keypone("Streamline_Skip", 1.));

  /* Preliminary checks */

  ndim = dbgrid->getNDim();
  if (ndim < 2 || ndim > 3)
  {
    messerr("This function is limited to 2-D or 3-D case");
    goto label_end;
  }
  npline = ndim + 3;

  /* Core allocation on the Grid Db */

  coor.resize(ndim);
  coor0.resize(ndim);
  iptr_time = dbgrid->addColumnsByConstant(1, TEST);
  if (iptr_time < 0) goto label_end;
  iptr_accu = dbgrid->addColumnsByConstant(1, 0.);
  if (iptr_accu < 0) goto label_end;

  /* Calculate the gradient */

  if (use_grad)
  {
    if (dbgrid->getNLoc(ELoc::G) != ndim)
    {
      messerr("When using the option 'use.grad'");
      messerr("the number of gradients should be %d (%d)", ndim,
              dbgrid->getNLoc(ELoc::G));
      goto label_end;
    }
    iptr_grad = dbgrid->getColIdxByLocator(ELoc::G, 0);
  }
  else
  {
    iptr_grad = db_gradient_components(dbgrid);
  }
  if (iptr_grad < 0) goto label_end;

  /* Normalize the gradient (optional) */

  if (flag_norm) st_gradient_normalize(dbgrid);

  /* Loop on the drop points */

  for (Id iech = 0; iech < dbpoint->getNSample(); iech++)
  {
    if (!dbpoint->isActive(iech)) continue;
    if (iech % nbyech != 0) continue;
    dbpoint->getCoordinatesInPlace(coor, iech);
    if (st_get_next(dbgrid, iptr_grad, coor, &knd, &surf)) break;

    /* Store the new point in the Streamline */

    if (nbline >= nquant * quant)
    {
      nquant++;
      line.resize(npline * nquant * quant);
    }
    for (idim = ecr = 0; idim < ndim; idim++)
      LINE(nbline, ecr++) = coor[idim];
    LINE(nbline, ecr++) = surf;
    LINE(nbline, ecr++) = knd + 1.;
    LINE(nbline, ecr++) = 0.;
    nbline++;

    for (Id i = 0; i < niter; i++)
    {
      for (idim = 0; idim < ndim; idim++)
        coor[idim] -= step * dbgrid->getArray(knd, iptr_grad + idim);
      if (st_get_next(dbgrid, iptr_grad, coor, &knd, &surf)) break;

      /* Store the new point in the Streamline */

      if (nbline >= nquant * quant)
      {
        nquant++;
        line.resize(npline * nquant * quant);
      }
      for (idim = ecr = 0; idim < ndim; idim++)
        LINE(nbline, ecr++) = coor[idim];
      LINE(nbline, ecr++) = surf;
      LINE(nbline, ecr++) = knd + 1.;
      LINE(nbline, ecr++) = i + 1.;
      nbline++;

      /* Update variables in the grid Db */

      date = MIN(dbgrid->getArray(knd, iptr_time), i + 1.);
      dbgrid->setArray(knd, iptr_time, date);
      dbgrid->updArray(knd, iptr_accu, EOperator::ADD, 1.);
    }

    /* Add the endpoint */

    if (nbline >= nquant * quant)
    {
      nquant++;
      line.resize(npline * nquant * quant);
    }
    for (idim = ecr = 0; idim < ndim; idim++)
      LINE(nbline, ecr++) = TEST;
    LINE(nbline, ecr++) = TEST;
    LINE(nbline, ecr++) = TEST;
    LINE(nbline, ecr++) = TEST;
    nbline++;
  }

  /* Final reallocation */

  line.resize(npline * nbline);

  /* Set the error return code */

  *nbline_loc = nbline;
  *npline_loc = npline;
  error       = 0;

label_end:
  if (!use_grad && !save_grad && iptr_grad >= 0)
    dbgrid->deleteColumnsByUIDRange(iptr_grad, ndim);
  return (error);
}

/*****************************************************************************/
/*!
 **  Smooth out the VPC
 **
 ** \return  Distance value
 **
 ** \param[in]  db          3-D Db structure containing the VPCs
 ** \param[in]  width       Width of the Filter
 ** \param[in]  range       Range of the Gaussian Weighting Function
 **
 ** \remarks Work is performed IN PLACE
 **
 *****************************************************************************/
Id db_smooth_vpc(DbGrid* db, Id width, double range)
{
  Id iz, nz, nprop, ecr, nkern, jz, error;
  double total, propval, dz, quant, quant0;
  VectorDouble prop1;
  VectorDouble prop2;
  VectorDouble kernel;

  /* Initializations */

  error = 1;
  nprop = db->getNLoc(ELoc::P);
  nz    = db->getNX(2);
  dz    = db->getDX(2);

  /* Core allocation */

  quant0 = law_invcdf_gaussian(0.975);
  if (FFFF(range))
    range = dz * width / quant0;
  else if (IFFFF(width))
    width = static_cast<Id>(range * quant0 / dz);
  else
  {
    messerr("You must define either 'width' or 'range'");
    goto label_end;
  }
  nkern = 2 * width + 1;
  prop1.resize(nz * nprop);
  prop2.resize(nz * nprop);
  kernel.resize(nkern);

  /* Establish the Kernel */

  for (Id i = -width; i <= width; i++)
  {
    quant             = (i * dz) / range;
    kernel[i + width] = law_df_gaussian(quant) / range;
  }

  /* Preliminary checks */

  if (!db->isGrid() || db->getNDim() != 3) goto label_end;

  /* Loop on the 2-D grid cells */

  for (Id ix = 0; ix < db->getNX(0); ix++)
    for (Id iy = 0; iy < db->getNX(1); iy++)
    {

      /* Load the proportions */

      if (db_prop_read(db, ix, iy, prop1.data())) goto label_end;

      /* Loop on the proportions */

      for (Id iprop = 0; iprop < nprop; iprop++)
      {

        /* Loop on the levels of the VPC */

        for (iz = ecr = 0; iz < nz; iz++)
        {

          /* Loop on the kernel items */

          total = 0.;
          for (Id i = -width; i <= width; i++)
          {
            jz      = Grid::generateMirrorIndex(nz, iz + i);
            propval = PROP1(jz, iprop);
            total += kernel[i + width] * propval;
          }
          PROP2(iz, iprop) = total;
        }
      }
      if (db_prop_write(db, ix, iy, prop2.data())) goto label_end;
    }

  /* Set the error return code */

  error = 0;

label_end:
  return (error);
}

/*****************************************************************************/
/*!
 **  Regularize variables along vertical wells
 **
 ** \return  Pointer to the newly created Db
 **
 ** \param[in]  db          Initial Db
 ** \param[in]  dbgrid      Reference Grid
 ** \param[in]  flag_center When TRUE, the sample is centered in the layer
 **                         to which it belongs
 **
 ** \remarks This function requires the input well ('db') and the grid to be
 ** \remarks defined in space >= 3D
 ** \remarks It requires a CODE variable to be defined in the input 'db'
 ** \remarks This function regularizes all the variables marked with a Z-locator
 ** \remarks This function takes a sample into account only if isotopic
 **
 *****************************************************************************/
Db* db_regularize(Db* db, DbGrid* dbgrid, Id flag_center)
{
  Db* dbnew = nullptr;
  if (db == nullptr) return dbnew;

  // Preliminary checks */

  if (!dbgrid->isGrid())
  {
    messerr("This function requires 'dbgrid' to correspond to a Grid");
    return dbnew;
  }

  if (db->getNDim() < 3)
  {
    messerr("This function requires the 'db' to be defined in 3D or more");
    return dbnew;
  }

  if (dbgrid->getNDim() < 3)
  {
    messerr("This function requires the 'dbgrid' to be defined in 3D or more");
    return dbnew;
  }

  if (!db->hasLocVariable(ELoc::C))
  {
    messerr("This function requires the definition of a CODE variable in 'db'");
    return dbnew;
  }

  if (!db->isNVarComparedTo(1, 1))
  {
    messerr("You should define some Z-variables in input 'db'");
    return dbnew;
  }

  // Core allocation

  Id iz   = 0;
  Id nz   = dbgrid->getNX(2);
  Id nvar = db->getNLoc(ELoc::Z);
  Id ndim = db->getNDim();
  Id size = ndim + nvar + 1;

  VectorDouble codes = db->getCodeList();
  Id ncode           = static_cast<Id>(codes.size());
  VectorDouble coor(ndim, 0);
  VectorDouble wcnt(ncode * nz, 0);
  VectorDouble wcor(ncode * nz * ndim, 0);
  VectorDouble wtab(ncode * nz * nvar, 0);

  // Loop on the different samples

  Id ntot = db->getNSample();

  for (Id iech = 0; iech < ntot; iech++)
  {
    if (!db->isActive(iech)) continue;
    mes_process("Regularize Wells", ntot, iech);
    Id code = db->getLocVariable(ELoc::C, iech, 0);

    // Identify the rank of the code

    Id icode = -1;
    for (Id i = 0; i < ncode && icode < 0; i++)
      if (isEqual(code, codes[i])) icode = i;
    if (icode < 0) continue;

    // Load the coordinates

    for (Id idim = 0; idim < ndim; idim++)
      coor[idim] = db->getCoordinate(iech, idim);

    Id err = point_to_bench(dbgrid, coor.data(), 0, &iz);
    if (err < 0) continue;
    if (iz < 0 || iz >= nz) continue;

    // Check if all variables are defined

    Id not_defined = 0;
    for (Id ivar = 0; ivar < nvar && not_defined == 0; ivar++)
      if (FFFF(db->getZVariable(iech, ivar))) not_defined = 1;
    if (not_defined) continue;

    // Cumulate this sample

    WCNT(iz, icode) += 1.;
    for (Id idim = 0; idim < ndim; idim++)
      WCOR(iz, icode, idim) += db->getCoordinate(iech, idim);
    for (Id ivar = 0; ivar < nvar; ivar++)
      WTAB(iz, icode, ivar) += db->getZVariable(iech, ivar);
  }

  // Normalization

  Id nech = 0;
  for (Id icode = 0; icode < ncode; icode++)
    for (iz = 0; iz < nz; iz++)
    {
      double ratio = WCNT(iz, icode);
      if (ratio <= 0) continue;
      for (Id idim = 0; idim < ndim; idim++)
        WCOR(iz, icode, idim) /= ratio;
      if (flag_center)
        WCOR(iz, icode, 2) = dbgrid->getX0(2) + (0.5 + iz) * dbgrid->getDX(2);
      for (Id ivar = 0; ivar < nvar; ivar++)
        WTAB(iz, icode, ivar) /= ratio;
      nech++;
    }

  // Load in storing array

  VectorDouble wecr(size * nech);

  Id ecr = 0;
  for (Id icode = 0; icode < ncode; icode++)
    for (iz = 0; iz < nz; iz++)
    {
      double ratio = WCNT(iz, icode);
      if (ratio <= 0) continue;
      for (Id idim = 0; idim < ndim; idim++)
        wecr[ecr++] = WCOR(iz, icode, idim);
      wecr[ecr++] = codes[icode];
      for (Id ivar = 0; ivar < nvar; ivar++)
        wecr[ecr++] = WTAB(iz, icode, ivar);
    }

  // Create the new db

  dbnew = Db::createFromSamples(nech, ELoadBy::SAMPLE, wecr, VectorString(), VectorString(), 0);
  if (dbnew == nullptr) goto label_end;

  ecr = 0;
  dbnew->setLocatorsByUID(ndim, ecr, ELoc::X, 0);
  ecr += ndim;
  dbnew->setLocatorByUID(ecr, ELoc::C, 0);
  ecr += 1;
  dbnew->setLocatorsByUID(nvar, ecr, ELoc::Z, 0);
  ecr += nvar;
  DECLARE_UNUSED(ecr);

label_end:
  return dbnew;
}

/*****************************************************************************/
/*!
 **  Sampling a fine grid in a coarser set of cells
 **
 ** \return  Error returned code
 **
 ** \param[in]  dbgrid      reference Grid
 ** \param[in]  nvar        Number of variables
 ** \param[in]  vars        Array of variable ranks
 ** \param[in]  npacks      Vector of packing factors
 ** \param[in]  npcell      Number of samples per cell
 ** \param[in]  nmini       Minimum number of nodes before drawing
 **
 ** \param[out] nech_ret    Number of selected samples
 ** \param[out] coor        Array of coordinates
 ** \param[out] data        Array of variables
 **
 ** \remarks The returned arrays 'coor' and 'data' must be freed by
 ** \remarks the calling function
 **
 *****************************************************************************/
Id db_grid2point_sampling(DbGrid* dbgrid,
                          Id nvar,
                          Id* vars,
                          const Id* npacks,
                          Id npcell,
                          Id nmini,
                          Id* nech_ret,
                          VectorDouble& coor,
                          VectorDouble& data)
{
  Id ndim, ntotal, nech, nret, nfine, iech, ecrc, ecrd, error;
  VectorInt ranks;
  VectorInt retain;
  VectorDouble rndval;

  // Initializations

  *nech_ret = 0;

  error = 1;
  ndim  = dbgrid->getNDim();
  nfine = dbgrid->getNSample();
  nmini = MAX(nmini, npcell);
  VectorInt indg(ndim, 0);
  if (ndim > 3)
  {
    messerr("This function is limited to 3D or less");
    goto label_end;
  }

  // Core allocation

  ntotal = 1;
  for (Id idim = 0; idim < ndim; idim++)
    ntotal *= npacks[idim];
  rndval.resize(ntotal);
  ranks.resize(ntotal);
  retain.resize(nfine);

  // Dispatch

  nret = 0;
  if (ndim == 1)
  {
    for (Id ixcell = 0; ixcell < dbgrid->getNX(0); ixcell += npacks[0])
    {

      // Collect eligible samples

      nech = 0;
      for (Id ix = 0; ix < npacks[0]; ix++)
      {
        indg[0] = ixcell + ix;
        if (indg[0] >= dbgrid->getNX(0)) break;
        iech = dbgrid->indiceToRank(indg);
        if (dbgrid->isActive(iech)) ranks[nech++] = iech;
      }
      if (nech < nmini) continue;

      // Draw sample(s) at random

      for (Id i = 0; i < nech; i++)
        rndval[i] = law_uniform(0., 1.);
      VH::arrangeInPlace(0, ranks, rndval, true, nech);
      for (Id i = 0; i < npcell; i++)
        retain[nret++] = ranks[i];
    }
  }
  else if (ndim == 2)
  {
    for (Id ixcell = 0; ixcell < dbgrid->getNX(0); ixcell += npacks[0])
      for (Id iycell = 0; iycell < dbgrid->getNX(1); iycell += npacks[1])
      {

        // Collect eligible samples

        nech = 0;
        for (Id ix = 0; ix < npacks[0]; ix++)
          for (Id iy = 0; iy < npacks[1]; iy++)
          {
            indg[0] = ixcell + ix;
            if (indg[0] >= dbgrid->getNX(0)) break;
            indg[1] = iycell + iy;
            if (indg[1] >= dbgrid->getNX(1)) break;
            iech = dbgrid->indiceToRank(indg);
            if (dbgrid->isActive(iech)) ranks[nech++] = iech;
          }
        if (nech < nmini) continue;

        // Draw sample(s) at random

        for (Id i = 0; i < nech; i++)
          rndval[i] = law_uniform(0., 1.);
        VH::arrangeInPlace(0, ranks, rndval, true, nech);
        for (Id i = 0; i < npcell; i++)
          retain[nret++] = ranks[i];
      }
  }
  else
  {
    for (Id ixcell = 0; ixcell < dbgrid->getNX(0); ixcell += npacks[0])
      for (Id iycell = 0; iycell < dbgrid->getNX(1); iycell += npacks[1])
        for (Id izcell = 0; izcell < dbgrid->getNX(2); izcell += npacks[2])
        {

          // Collect eligible samples

          nech = 0;
          for (Id ix = 0; ix < npacks[0]; ix++)
            for (Id iy = 0; iy < npacks[1]; iy++)
              for (Id iz = 0; iz < npacks[2]; iz++)
              {
                indg[0] = ixcell + ix;
                if (indg[0] >= dbgrid->getNX(0)) break;
                indg[1] = iycell + iy;
                if (indg[1] >= dbgrid->getNX(1)) break;
                indg[2] = izcell + iz;
                if (indg[2] >= dbgrid->getNX(2)) break;
                iech = dbgrid->indiceToRank(indg);
                if (dbgrid->isActive(iech)) ranks[nech++] = iech;
              }
          if (nech < nmini) continue;

          // Draw sample(s) at random

          for (Id i = 0; i < nech; i++)
            rndval[i] = law_uniform(0., 1.);
          VH::arrangeInPlace(0, ranks, rndval, true, nech);
          for (Id i = 0; i < npcell; i++)
            retain[nret++] = ranks[i];
        }
  }

  // Allocate the array for coordinates and data

  coor.resize(ndim * nret);
  data.resize(nvar * nret);

  // Load the returned arrays

  ecrc = ecrd = 0;
  for (Id i = 0; i < nret; i++)
  {
    iech = retain[i];
    for (Id idim = 0; idim < ndim; idim++)
      coor[ecrc++] = dbgrid->getCoordinate(iech, idim);
    for (Id ivar = 0; ivar < nvar; ivar++)
      data[ecrd++] = dbgrid->getArray(iech, vars[ivar]);
  }

  // Set the error return code

  *nech_ret = nret;
  error     = 0;

  // Core deallocation

label_end:
  return (error);
}

/*****************************************************************************/
/*!
 **  Create a new Data Base with points generated at random
 **
 ** \return  Pointer for the new Db structure
 **
 ** \param[in]  nech        Expected number of samples
 ** \param[in]  coormin     Vector of lower coordinates
 ** \param[in]  coormax     Vector of upper coordinates
 ** \param[in]  dbgrid      Descriptor of the Db grid parameters
 ** \param[in]  flag_exact  True if the number of samples is dran from Poisson
 ** \param[in]  flag_repulsion Use repulsion (need: 'range' and 'beta')
 ** \param[in]  range       Repulsion range
 ** \param[in]  beta        Bending coefficient
 ** \param[in]  extend      Extension of the bounding box (when positive)
 ** \param[in]  seed        Seed for the random number generator
 ** \param[in]  flagAddSampleRank true if the Rank must be generated in the output Db
 **
 ** \remarks Arguments 'extend' is only valid when 'dbgrid' is not defined
 **
 *****************************************************************************/
Db* db_point_init(Id nech,
                  const VectorDouble& coormin,
                  const VectorDouble& coormax,
                  DbGrid* dbgrid,
                  bool flag_exact,
                  bool flag_repulsion,
                  double range,
                  double beta,
                  double extend,
                  Id seed,
                  bool flagAddSampleRank)
{
  VectorDouble tab;
  Db* db  = nullptr;
  Id ndim = 0;
  if (dbgrid == nullptr)
    ndim = static_cast<Id>(coormin.size());
  else
    ndim = dbgrid->getNDim();
  if (ndim <= 0) return db;

  // Initiate the pseudo-random number generator

  law_set_random_seed(seed);

  // Process the bounding box extension (optional)

  VectorDouble locmin = coormin;
  VectorDouble locmax = coormax;
  if (extend > 0)
  {
    for (Id idim = 0; idim < ndim; idim++)
    {
      locmin[idim] -= extend;
      locmax[idim] += extend;
    }
  }

  // Draw the number of data to be generated in the Poisson process

  Id number = nech;
  if (!flag_exact) law_poisson(nech);

  // Dispatch

  if (number > 0)
  {
    if (dbgrid == nullptr)
    {
      tab = st_point_init_homogeneous(number, locmin, locmax,
                                      flag_repulsion, range, beta);
    }
    else
    {
      tab = st_point_init_inhomogeneous(number, dbgrid,
                                        flag_repulsion, range, beta);
    }
  }

  /* Allocate the main structure */

  number = static_cast<Id>(tab.size()) / ndim;
  db     = Db::createFromSamples(number, ELoadBy::SAMPLE, tab, VectorString(),
                                 VectorString(), flagAddSampleRank);

  /* Set the locators */

  VectorString names = generateMultipleNames("x", ndim);
  for (Id idim = 0; idim < ndim; idim++)
  {
    Id jdim = (flagAddSampleRank) ? idim + 1 : idim;
    db->setNameByUID(jdim, names[idim]);
    db->setLocatorByUID(jdim, ELoc::X, idim);
  }
  return (db);
}

/****************************************************************************/
/*!
 **  Standard Kriging
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  model       Model structure
 ** \param[in]  niter       Number of iterations
 ** \param[in]  verbose     Verbose flag
 ** \param[in]  namconv     Naming convention
 **
 ** \remarks The procedure uses the FIRST covariance of the Model
 ** \remarks to describe the spatial structure
 **
 *****************************************************************************/
Id db_proportion_estimate(Db* dbin,
                          DbGrid* dbout,
                          Model* model,
                          Id niter,
                          bool verbose,
                          const NamingConvention& namconv)
{
  VectorVectorInt splits;

  // Preliminary checks

  if (dbin == nullptr)
  {
    messerr("This method requires a 'dbin' argument");
    return 1;
  }
  if (dbout == nullptr)
  {
    messerr("This method requires a 'dbout' argument");
    return 1;
  }
  if (model == nullptr)
  {
    messerr("This method requires a 'model' argument");
    return 1;
  }
  if (dbin->getNLoc(ELoc::Z) != 1)
  {
    messerr("The argument 'dbin' should have a single variable");
    return 1;
  }

  // Define the environment

  MeshETurbo mesh(dbout);
  ShiftOpMatrix S(&mesh, model->getCovAniso(0), dbout);
  PrecisionOp Qprop(&S, model->getCovAniso(0));
  ProjMatrix AprojDat(dbin, &mesh);
  ProjMatrix AprojOut(dbout, &mesh);

  // Invoke the calculation

  VectorDouble propGlob = dbStatisticsFacies(dbin);
  Id ncat               = static_cast<Id>(propGlob.size());
  OptimCostColored Oc(ncat, &Qprop, &AprojDat);

  VectorDouble facies      = dbin->getColumnByLocator(ELoc::Z);
  VectorVectorDouble props = Oc.minimize(facies, splits, propGlob, verbose, niter);

  // Loading the resulting results in the output 'dbout'

  Id iptr0 = -1;
  VectorDouble propout(dbout->getNSample(true));
  for (Id i = 0; i < ncat; i++)
  {
    AprojOut.mesh2point(props[i], propout);
    Id iptr = dbout->addColumns(propout, String(), ELoc::UNKNOWN, 0, true);
    if (i == 0) iptr0 = iptr;
    namconv.setNamesAndLocators(nullptr, VectorString(), ELoc::UNKNOWN, -1, dbout, iptr,
                                concatenateStrings("-", toString(i + 1)));
  }
  namconv.setLocators(dbout, iptr0, 1, ncat);

  return 0;
}
} // namespace gstlrn
