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
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "Variogram/Vario.hpp"
#include "Variogram/VarioParam.hpp"
#include "Anamorphosis/Anam.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Polynomials/Hermite.hpp"
#include "Morpho/Morpho.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"
#include "Drifts/EDrift.hpp"
#include "Basic/EJustify.hpp"
#include "Basic/File.hpp"
#include "Basic/String.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Stats/PCA.hpp"

#include <string.h>
#include <math.h>

/*! \cond */
#define VARS(ivar,jvar)     (vario->vars[(ivar) * vario->getNVar() + (jvar)])
#define POISSON_MEANS(ivar) (VARIO->means[(ivar)])
#define DATES(idate,i)      (vario->dates[2 * (idate) + (i)])
#define IAD(ivar,jvar)      ((ivar) * nvar + (jvar))
#define ADD(ix,iy,iz,nx)    ((iz) + nx[2] * ((iy) + nx[1] * (ix)))
#define OPP(idim,i)         (dims[idim] - i - 1)
#define X_MATDRF(il,jl)     (MATDRF[(il) * nbfl + (jl)])
#define X_DRFXGX(il,jl)     (DRFXGX[(il) * nbfl + (jl)])
#define X_DRFTAB(il,iech)   (DRFTAB[(il) * nech + (iech)])
#define X_DRFXA(il,iech)    (DRFXA [(il) * nech + (iech)])
#define X_DRFGX(il,iech)    (DRFGX [(il) * nech + (iech)])

/*! \endcond */

static Vario *VARIO;
static Model *MODEL;
static Db *DBMAP;
static int IDIRLOC;
static double *BETA, *DRFLOC, *DRFTAB, *MATDRF, *DRFXA, *DRFGX, *DRFXGX, *DRFDIAG;
static int IECH1, IECH2, IPTV, IPTW;

static int NWGT[4] = { 2, 3, 4, 5 };
static int NORWGT[4] = { 2, 6, 20, 70 };
static int VARWGT[4][5] = { { 1, -1, 0, 0, 0 },
                            { 1, -2, 1, 0, 0 },
                            { 1, -3, 3, -1, 0 },
                            { 1, -4, 6, -4, 1 } };

/****************************************************************************/
/*!
 **  Fix plausible values for the Direction coefficients.
 **  They must be defined and with norm equal to 1
 **
 ** \param[in]  ndim      Space dimension
 ** \param[in,out]  codir Input/Output Direction coefficients
 **
 *****************************************************************************/
void vario_fix_codir(int ndim, VectorDouble &codir)
{
  double norme;

  if (codir.empty()) return;
  norme = matrix_norm(codir.data(), ndim);
  if (norme <= 0.)
  {
    for (int idim = 0; idim < ndim; idim++)
      codir[idim] = 0.;
    codir[0] = 1.;
  }
  else
  {
    norme = sqrt(norme);
    for (int idim = 0; idim < ndim; idim++)
      codir[idim] /= norme;
  }
}

/****************************************************************************/
/*!
 **  Get variable order
 **
 ** \return  Rank of the pair of variables (-1 if incorrect arguments)
 **
 ** \param[in]  nvar      Number of variables
 ** \param[in]  ivar0     Rank of the first variable
 ** \param[in]  jvar0     Rank of the second variable
 **
 *****************************************************************************/
static int st_get_variable_order(int nvar, int ivar0, int jvar0)
{
  int rank, ivar, jvar;

  for (ivar = rank = 0; ivar < nvar; ivar++)
    for (jvar = 0; jvar <= ivar; jvar++, rank++)
    {
      if (ivar == ivar0 && jvar == jvar0) return (rank);
      if (ivar == jvar0 && jvar == ivar0) return (rank);
    }
  return (-1);
}

/****************************************************************************/
/*!
 **  Return the Order of the generalized variogram
 **
 ** \return Order of the Generalized covariance
 **
 ** \param[in]  vario   Vario structure
 **
 *****************************************************************************/
static int st_get_generalized_variogram_order(const Vario *vario)
{
  int norder;

  norder = 0;
  if (vario->getCalculType() == ECalcVario::GENERAL1) norder = 1;
  if (vario->getCalculType() == ECalcVario::GENERAL2) norder = 2;
  if (vario->getCalculType() == ECalcVario::GENERAL3) norder = 3;
  return (norder);
}

/****************************************************************************/
/*!
 **  Manage the drift removal option
 **
 ** \param[in]  type    Type of action
 **                     0 (initialization)
 **                     +1 (allocation)
 **                     -1 (deallocation)
 ** \param[in]  db      Db structure
 ** \param[in]  model   Model structure
 **
 ** \remark If the drift removal is not called, set argument 'model' to NULL
 **
 *****************************************************************************/
static void st_manage_drift_removal(int type, Db *db, Model *model)
{
  int nbfl, nech, i;

  /* Dispatch */

  switch (type)
  {
    case 0:
      MODEL = nullptr;
      BETA = nullptr;
      DRFLOC = nullptr;
      DRFTAB = nullptr;
      MATDRF = nullptr;
      DRFXA = nullptr;
      DRFGX = nullptr;
      DRFXGX = nullptr;
      DRFDIAG = nullptr;
      break;

    case 1:
      if (model != nullptr)
      {
        nbfl = model->getDriftNumber();
        nech = db->getActiveAndDefinedNumber(0);
        MODEL = model;
        DRFLOC = (double*) mem_alloc(sizeof(double) * nbfl, 1);
        BETA = (double*) mem_alloc(sizeof(double) * nbfl, 1);
        for (i = 0; i < nbfl; i++)
          BETA[i] = 0.;
        MATDRF = (double*) mem_alloc(sizeof(double) * nbfl * nbfl, 1);
        for (i = 0; i < nbfl * nbfl; i++)
          MATDRF[i] = 0.;
        DRFTAB = (double*) mem_alloc(sizeof(double) * nech * nbfl, 1);
        for (i = 0; i < nech * nbfl; i++)
          DRFTAB[i] = 0.;
        DRFXA = (double*) mem_alloc(sizeof(double) * nech * nbfl, 1);
        for (i = 0; i < nech * nbfl; i++)
          DRFXA[i] = 0.;
        DRFGX = (double*) mem_alloc(sizeof(double) * nech * nbfl, 1);
        for (i = 0; i < nech * nbfl; i++)
          DRFGX[i] = 0.;
        DRFXGX = (double*) mem_alloc(sizeof(double) * nbfl * nbfl, 1);
        for (i = 0; i < nbfl * nbfl; i++)
          DRFXGX[i] = 0.;
        DRFDIAG = (double*) mem_alloc(sizeof(double) * nech, 1);
        for (i = 0; i < nech; i++)
          DRFDIAG[i] = 0.;
      }
      break;

    case -1:
      MODEL = nullptr;
      BETA = (double*) mem_free((char* ) BETA);
      DRFLOC = (double*) mem_free((char* ) DRFLOC);
      DRFTAB = (double*) mem_free((char* ) DRFTAB);
      MATDRF = (double*) mem_free((char* ) MATDRF);
      DRFXA = (double*) mem_free((char* ) DRFXA);
      DRFGX = (double*) mem_free((char* ) DRFGX);
      DRFXGX = (double*) mem_free((char* ) DRFXGX);
      DRFDIAG = (double*) mem_free((char* ) DRFDIAG);
      break;
  }
}

/****************************************************************************/
/*!
 **  Checks if the maximum variogram distance has been passed
 **
 ** \return    1 if the maximum distance has been passed and 0 otherwise      IDIRLOC = idir;
 **
 ** \param[in]  db      Db descriptor
 ** \param[in]  iech    Rank of the first sample
 ** \param[in]  jech    Rank of the second sample
 ** \param[in]  maxdist Maximum distance
 **
 *****************************************************************************/
int variogram_maximum_dist1D_reached(Db *db, int iech, int jech, double maxdist)
{
  double dist = db->getDistance1D(iech, jech, 0, true);
  return (dist > maxdist);
}

/****************************************************************************/
/*!
 **  Calculate the data value (possibly after removing the global trend)
 **
 ** \return    The data value (or the residual)
 **
 ** \param[in]  db      Db descriptor
 ** \param[in]  iech    Rank of the sample
 ** \param[in]  ivar    Rank of the variable
 **
 ** \remark  The trend removal only applies on the first variable
 ** \remark  Therefore, if applied on any variable rank other than 0,
 ** \remark  TEST is returned
 **
 *****************************************************************************/
static double st_get_IVAR(const Db *db, int iech, int ivar)
{
  double zz, drfval;

  zz = db->getVariable(iech, ivar);
  if (FFFF(zz)) return (TEST);
  if (MODEL == nullptr) return (zz);
  if (ivar != 0) return (TEST);
  drfval = model_drift_evaluate(0, MODEL, db, iech, 0, BETA, DRFLOC);
  if (FFFF(drfval)) return (TEST);
  return (zz - drfval);
}

/****************************************************************************/
/*!
 **  Calculates the statistics for the variogram calculations
 **
 ** \param[in]  db      Db descriptor
 ** \param[in]  vario   Vario structure
 **
 *****************************************************************************/
static void st_variogram_stats(Db *db, Vario *vario)
{
  double z1, z2, ww;

  /* Initializations */

  for (int ivar = 0; ivar < db->getVariableNumber(); ivar++)
  {
    vario->setMean(ivar, 0.);
    for (int jvar = 0; jvar < db->getVariableNumber(); jvar++)
      vario->setVar(ivar, jvar, 0.);
  }

  /* Loop on the variables */

  for (int ivar = 0; ivar < db->getVariableNumber(); ivar++)
  {
    double s1w = 0.;
    double s1z = 0.;
    for (int iech = 0; iech < db->getSampleNumber(); iech++)
    {
      if (!db->isActive(iech)) continue;
      ww = db->getWeight(iech);
      if (FFFF(ww) || ww < 0.) continue;
      z1 = st_get_IVAR(db, iech, ivar);
      s1w += ww;
      s1z += z1;
    }

    if (s1w <= 0.) continue;
    vario->setMean(ivar, s1z / s1w);
  }

  for (int ivar = 0; ivar < db->getVariableNumber(); ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
    {

      /* Loop on the samples */

      double s12w = 0.;
      double s12wz1 = 0.;
      double s12wz2 = 0.;
      double s12wzz = 0.;

      for (int iech = 0; iech < db->getSampleNumber(); iech++)
      {
        if (!db->isActive(iech)) continue;
        ww = db->getWeight(iech);
        if (FFFF(ww) || ww < 0.) continue;

        z1 = st_get_IVAR(db, iech, ivar);
        z2 = st_get_IVAR(db, iech, jvar);
        if (FFFF(z1) || FFFF(z2)) continue;

        s12w += ww;
        s12wz1 += ww * z1;
        s12wz2 += ww * z2;
        s12wzz += ww * z1 * z2;
      }
      if (s12w <= 0.) continue;

      if (vario->getCalculType() == ECalcVario::COVARIOGRAM)
      {
        vario->setVar(ivar, jvar, s12wzz);
        vario->setVar(jvar, ivar, s12wzz);
      }
      else
      {
        vario->setVar(ivar, jvar,
                       s12wzz / s12w - (s12wz1 / s12w) * (s12wz2 / s12w));
        vario->setVar(jvar, ivar,
                       s12wzz / s12w - (s12wz1 / s12w) * (s12wz2 / s12w));
      }
    }

  // Modification when the ultimate variogram is a transformed one

  if (vario->getCalculType() == ECalcVario::TRANS1)
  {
    for (int ivar = 0; ivar < db->getVariableNumber(); ivar++)
      for (int jvar = 0; jvar < ivar; jvar++)
      {
        double value = -vario->getVar(ivar, jvar) / vario->getVar(jvar, jvar);
        vario->setVar(ivar, jvar, value);
        vario->setVar(jvar, ivar, value);
      }
  }
  else if (vario->getCalculType() == ECalcVario::TRANS2)
  {
    for (int ivar = 0; ivar < db->getVariableNumber(); ivar++)
      for (int jvar = 0; jvar < ivar; jvar++)
      {
        double value = -vario->getVar(ivar, jvar) / vario->getVar(ivar, ivar);
        vario->setVar(ivar, jvar, value);
        vario->setVar(jvar, ivar, value);
      }
  }
  else if (vario->getCalculType() == ECalcVario::BINORMAL)
  {
    for (int ivar = 0; ivar < db->getVariableNumber(); ivar++)
      for (int jvar = 0; jvar < db->getVariableNumber(); jvar++)
        if (ivar != jvar)
          vario->setVar(
              ivar,
              jvar,
              vario->getVar(ivar, jvar) / sqrt(
                  vario->getVar(ivar, ivar) * vario->getVar(jvar, jvar)));
  }
  return;
}

/****************************************************************************/
/*!
 **  Calculates the statistics for the covariance calculations
 **  for directional variables
 **
 ** \param[in]  db     Db descriptor
 ** \param[in]  vario  Vario structure
 ** \param[in]  ncomp  Number of components
 **
 *****************************************************************************/
static void st_variovect_stats(Db *db, Vario *vario, int ncomp)
{
  double vi, vj, vij, s12ww, s12wzz, zi, zj, ww;
  int iech, ivar, jvar, nb_neg, icomp;

  /* Loop on the variables */

  nb_neg = 0;
  for (ivar = 0; ivar < vario->getVariableNumber(); ivar++)
    for (jvar = 0; jvar <= ivar; jvar++)
    {

      /* Loop on the samples */

      s12ww = s12wzz = 0.;

      for (iech = 0; iech < db->getSampleNumber(); iech++)
      {
        if (!db->isActive(iech)) continue;
        ww = db->getWeight(iech);
        if (FFFF(ww)) continue;
        if (ww < 0.)
        {
          nb_neg++;
          continue;
        }

        /* Loop on the components */

        vi = vj = vij = 0;
        for (icomp = 0; icomp < ncomp; icomp++)
        {
          zi = st_get_IVAR(db, iech, ivar * ncomp + icomp);
          zj = st_get_IVAR(db, iech, jvar * ncomp + icomp);
          if (!FFFF(zi) && !FFFF(zj))
          {
            vi += zi * zi;
            vj += zj * zj;
            vij += zi * zj;
          }
          else
          {
            vij = TEST;
            break;
          }
        }
        if (ABS(vi * vj) < 1.e-10 || FFFF(vij)) continue;
        vij = ABS(vij) / sqrt(vi * vj);

        s12ww += ww * ww;
        s12wzz += ww * ww * vij;
      }

      vario->setVar(ivar, jvar, (s12ww > 0) ? s12wzz / s12ww :
                                               0.);
      vario->setVar(jvar, ivar, (s12ww > 0) ? s12wzz / s12ww :
                                               0.);
    }

  if (nb_neg > 0)
    message("There were %d negative weights. They have been set to zero\n",
            nb_neg);
  return;
}

/****************************************************************************/
/*!
 **  Returns the cosine of the angular tolerance
 **
 ** \param[in]  tolang Angular tolerance
 **
 ** This method is not documented on purpose. It should remain private
 **
 *****************************************************************************/
double _variogram_convert_angular_tolerance(double tolang)

{
  double psval;

  if (FFFF(tolang))
    psval = 0.;
  else if (tolang == 00.)
    psval = 1.;
  else if (tolang == 90.)
    psval = 0.;
  else
    psval = ABS(cos(ut_deg2rad(tolang)));
  return (psval);
}

/****************************************************************************/
/*!
 **  Return the rank of the lag
 **
 ** \return  Rank of the lag or ITEST
 **
 ** \param[in]  vario        Vario structure
 ** \param[in]  idir         Dir rank
 ** \param[in]  ps           Cosinus of the angle
 ** \param[in]  psmin        Angular tolerance
 ** \param[in]  dist         Distance
 **
 *****************************************************************************/
int variogram_get_lag(Vario *vario,
                      int idir,
                      double ps,
                      double psmin,
                      double *dist)
{
  int k, ilag;
  const DirParam &dirparam = vario->getDirParam(idir);

  /* Determine the rank of the lag */

  ilag = -1;
  if (dirparam.getFlagRegular())
  {
    ilag = (int) floor((*dist) / dirparam.getDPas() + 0.5);
    if (ABS((*dist) - ilag * dirparam.getDPas()) > dirparam.getTolDist()
        * dirparam.getDPas()) return (ITEST);
  }
  else
  {
    for (k = 0, ilag = -1; k < dirparam.getLagNumber() && ilag < 0; k++)
      if ((*dist) > dirparam.getBreaks()[k] && (*dist)
          <= dirparam.getBreaks()[k + 1]) ilag = k;
  }
  if (ilag < 0 || ilag >= dirparam.getLagNumber()) return (ITEST);

  /* For asymmetric function, set the distance to negative value */
  /* for pairs of samples opposite to the calculation direction */
  /* Be sure to consider its absolute value */

  if (vario->getFlagAsym())
  {
    if (ps < psmin) (*dist) = -(*dist);
  }

  return (ilag);
}

/****************************************************************************/
/*!
 **  Printout function for Debug case
 **
 ** \param[in]  iech1 Rank of the first sample
 ** \param[in]  iech2 Rank of the second sample
 ** \param[in]  ivar  Rank of the first variable
 ** \param[in]  jvar  Rank of the second variable
 ** \param[in]  ilag  Rank of the Lag
 ** \param[in]  scale Weighting factor
 ** \param[in]  value Variogram value
 **
 *****************************************************************************/
static void st_print_debug(int iech1,
                           int iech2,
                           int ivar,
                           int jvar,
                           int ilag,
                           double scale,
                           double value)
{
  message(
      "Samples: %d/%d - Variables: %d/%d - Weight: %lf - Lag: %d - Variogram: %lf\n",
      iech1 + 1, iech2 + 1, ivar + 1, jvar + 1, scale, ilag, value);
  return;
}

/****************************************************************************/
/*!
 **  Internal function for setting a variogram value
 **
 ** \param[in]  calcul_type Type of calculation (ECalcVario)
 ** \param[in]  ipas        Rank of the variogram lag
 ** \param[in]  ivar        Index of the first variable
 ** \param[in]  jvar        Index of the second variable
 ** \param[in]  orient      Orientation
 ** \param[in]  ww          Weight
 ** \param[in]  dist        Distance
 ** \param[in]  value       Variogram value
 **
 *****************************************************************************/
static void st_variogram_set(const ECalcVario &calcul_type,
                             int /*nvar*/,
                             int ipas,
                             int ivar,
                             int jvar,
                             int orient,
                             double ww,
                             double dist,
                             double value)
{
  int i = VARIO->getDirAddress(IDIRLOC, ivar, jvar, ipas, false, orient);

  VARIO->updateGgByIndex(IDIRLOC, i, ww * value);
  if (calcul_type == ECalcVario::POISSON)
    VARIO->updateGgByIndex(IDIRLOC, i, -VARIO->getMean(ivar) / 2.);
  VARIO->updateHhByIndex(IDIRLOC, i, ww * dist);
  VARIO->updateSwByIndex(IDIRLOC, i, ww);
  if (debug_query("variogram"))
    st_print_debug(IECH1, IECH2, ivar, jvar, i, ww, value);
  return;
}

/****************************************************************************/
/*!
 **  Internal function for setting a VMAP value
 **
 ** \param[in]  nvar        Number of variables
 ** \param[in]  ipas        Rank of the variogram lag
 ** \param[in]  ivar        Index of the first variable
 ** \param[in]  jvar        Index of the second variable
 ** \param[in]  ww          Weight
 ** \param[in]  value       Variogram value
 **
 *****************************************************************************/
static void st_vmap_set(const ECalcVario& /*calcul_type*/,
                        int nvar,
                        int ipas,
                        int ivar,
                        int jvar,
                        int /*orient*/,
                        double ww,
                        double /*dist*/,
                        double value)
{
  int ijvar;

  ijvar = st_get_variable_order(nvar, ivar, jvar);

  DBMAP->updArray(ipas, IPTV + ijvar, 0, ww * value);
  DBMAP->updArray(ipas, IPTW + ijvar, 0, ww);

  if (debug_query("variogram"))
    st_print_debug(IECH1, IECH2, ivar, jvar, ipas, ww, value);
  return;
}

/****************************************************************************/
/*!
 **  Update the variogram values
 **
 ** \param[in]  db             Db descriptor
 ** \param[in]  calcul_type    Type of calculation (ECalcVario)
 ** \param[in]  nvar           Number of variables
 ** \param[in]  iech1          Rank of the first sample
 ** \param[in]  iech2          Rank of the second sample
 ** \param[in]  ipas           Rank of the lag
 ** \param[in]  dist           Distance value
 ** \param[in]  do_asym        When FALSE, do not perform the symmetry
 ** \param[in]  st_generic_set Generic function for setting the variogram
 **
 ** \remarks: The argument 'do_asym' allows performing the double assignment
 ** \remarks: for the asymmetric functions (such as covariance)
 ** \remarks: This is due to the fact that the double is called entirely
 ** \remarks: in the calling function
 **
 *****************************************************************************/
static void st_variogram_evaluate(Db *db,
                                  const ECalcVario &calcul_type,
                                  int nvar,
                                  int iech1,
                                  int iech2,
                                  int ipas,
                                  double dist,
                                  int do_asym,
                                  void (*st_generic_set)(const ECalcVario &calcul_type,
                                                         int nvar,
                                                         int iadlag,
                                                         int ivar,
                                                         int jvar,
                                                         int orient,
                                                         double ww,
                                                         double dist,
                                                         double value))
{
  double w1, w2, z11, z12, z21, z22, scale, value;
  int ivar, jvar, orient;

  w1 = db->getWeight(iech1);
  w2 = db->getWeight(iech2);
  if (FFFF(w1) || FFFF(w2)) return;
  orient = (dist > 0) ? 1 :
                        -1;
  dist = ABS(dist);

  switch (calcul_type.toEnum())
  {
    case ECalcVario::E_VARIOGRAM:
    case ECalcVario::E_TRANS1:
    case ECalcVario::E_TRANS2:
    case ECalcVario::E_BINORMAL:
      scale = w1 * w2;
      for (ivar = 0; ivar < nvar; ivar++)
        for (jvar = 0; jvar <= ivar; jvar++)
        {
          z11 = st_get_IVAR(db, iech1, ivar);
          z12 = st_get_IVAR(db, iech2, ivar);
          z21 = st_get_IVAR(db, iech1, jvar);
          z22 = st_get_IVAR(db, iech2, jvar);
          if (!FFFF(z11) && !FFFF(z21) && !FFFF(z12) && !FFFF(z22))
          {
            value = (z12 - z11) * (z22 - z21) / 2.;
            st_generic_set(calcul_type, nvar, ipas, ivar, jvar, 0, scale, dist,
                           value);
          }
        }
      break;

    case ECalcVario::E_MADOGRAM:
      scale = w1 * w2;
      for (ivar = 0; ivar < nvar; ivar++)
        for (jvar = 0; jvar <= ivar; jvar++)
        {
          z11 = st_get_IVAR(db, iech1, ivar);
          z12 = st_get_IVAR(db, iech2, ivar);
          z21 = st_get_IVAR(db, iech1, jvar);
          z22 = st_get_IVAR(db, iech2, jvar);
          if (!FFFF(z11) && !FFFF(z21) && !FFFF(z12) && !FFFF(z22))
          {
            value = sqrt(ABS((z12 - z11) * (z22 - z21))) / 2.;
            st_generic_set(calcul_type, nvar, ipas, ivar, jvar, 0, scale, dist,
                           value);
          }
        }
      break;

    case ECalcVario::E_RODOGRAM:
      scale = w1 * w2;
      for (ivar = 0; ivar < nvar; ivar++)
        for (jvar = 0; jvar <= ivar; jvar++)
        {
          ;
          z11 = st_get_IVAR(db, iech1, ivar);
          z12 = st_get_IVAR(db, iech2, ivar);
          z21 = st_get_IVAR(db, iech1, jvar);
          z22 = st_get_IVAR(db, iech2, jvar);
          if (!FFFF(z11) && !FFFF(z21) && !FFFF(z12) && !FFFF(z22))
          {
            value = pow(ABS((z12 - z11) * (z22 - z21)), 0.25) / 2.;
            st_generic_set(calcul_type, nvar, ipas, ivar, jvar, 0, scale, dist,
                           value);
          }
        }
      break;

    case ECalcVario::E_POISSON:
      scale = (w1 * w2) / (w1 + w2);
      for (ivar = 0; ivar < nvar; ivar++)
        for (jvar = 0; jvar <= ivar; jvar++)
        {
          z11 = st_get_IVAR(db, iech1, ivar);
          z12 = st_get_IVAR(db, iech2, ivar);
          z21 = st_get_IVAR(db, iech1, jvar);
          z22 = st_get_IVAR(db, iech2, jvar);
          if (!FFFF(z11) && !FFFF(z21) && !FFFF(z12) && !FFFF(z22)
              && (w1 > 0. && w2 > 0.))
          {
            value = (z12 - z11) * (z22 - z21) / 2.;
            st_generic_set(calcul_type, nvar, ipas, ivar, jvar, 0, scale, dist,
                           value);
          }
        }
      break;

    case ECalcVario::E_COVARIANCE:
    case ECalcVario::E_COVARIANCE_NC:
      scale = w1 * w2;
      for (ivar = 0; ivar < nvar; ivar++)
        for (jvar = 0; jvar <= ivar; jvar++)
        {
          z11 = st_get_IVAR(db, iech1, ivar);
          z12 = st_get_IVAR(db, iech2, ivar);
          z21 = st_get_IVAR(db, iech1, jvar);
          z22 = st_get_IVAR(db, iech2, jvar);
          if (!FFFF(z11) && !FFFF(z22))
          {
            value = z11 * z22;
            st_generic_set(calcul_type, nvar, ipas, ivar, jvar, orient, scale,
                           dist, value);
          }
          if (!FFFF(z12) && !FFFF(z21) && do_asym)
          {
            value = z12 * z21;
            st_generic_set(calcul_type, nvar, ipas, ivar, jvar, -orient, scale,
                           dist, value);
          }
        }
      break;

    case ECalcVario::E_COVARIOGRAM:
      scale = w2;
      for (ivar = 0; ivar < nvar; ivar++)
        for (jvar = 0; jvar <= ivar; jvar++)
        {
          z11 = st_get_IVAR(db, iech1, ivar);
          z12 = st_get_IVAR(db, iech2, ivar);
          z21 = st_get_IVAR(db, iech1, jvar);
          z22 = st_get_IVAR(db, iech2, jvar);
          if (!FFFF(z11) && !FFFF(z22))
          {
            value = z11 * z22;
            st_generic_set(calcul_type, nvar, ipas, ivar, jvar, orient, scale,
                           dist, value);
          }
          if (!FFFF(z12) && !FFFF(z21) && do_asym)
          {
            value = z12 * z21;
            st_generic_set(calcul_type, nvar, ipas, ivar, jvar, -orient, scale,
                           dist, value);
          }
        }
      break;

    case ECalcVario::E_ORDER4:
      scale = w1 * w2;
      for (ivar = 0; ivar < nvar; ivar++)
        for (jvar = 0; jvar <= ivar; jvar++)
        {
          z11 = st_get_IVAR(db, iech1, ivar);
          z12 = st_get_IVAR(db, iech2, ivar);
          z21 = st_get_IVAR(db, iech1, jvar);
          z22 = st_get_IVAR(db, iech2, jvar);
          if (!FFFF(z11) && !FFFF(z21) && !FFFF(z12) && !FFFF(z22))
          {
            value = (z12 - z11) * (z22 - z21);
            value = value * value / 2.;
            st_generic_set(calcul_type, nvar, ipas, ivar, jvar, 0, scale, dist,
                           value);
          }
        }
      break;

    default:
      messageAbort("st_variogram_evaluate");
      break;
  }

  return;
}

/****************************************************************************/
/*!
 **  Scale the variogram calculations
 **
 ** \param[in]  vario Vario structure
 ** \param[in]  idir  Rank of the Direction
 **
 *****************************************************************************/
void variogram_scale(Vario *vario, int idir)
{
  int nvar = vario->getVariableNumber();

  /* Scale the experimental variogram quantities */

  int ecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
    {
      for (int i = 0; i < vario->getLagTotalNumber(idir); i++, ecr++)
      {
        int j = vario->getDirAddress(idir, ivar, jvar, i, true, 0);
        if (vario->getSwByIndex(idir, j) <= 0)
        {
          vario->setHhByIndex(idir, j, TEST);
          vario->setGgByIndex(idir, j, TEST);
        }
        else
        {
          vario->setHhByIndex(idir, j, vario->getHhByIndex(idir, j) / vario->getSwByIndex(idir, j));
          if (vario->getFlagAsym() && i < vario->getLagNumber(idir))
            vario->setHhByIndex(idir, j, -ABS(vario->getHhByIndex(idir, j)));
          if (vario->getCalculType() != ECalcVario::COVARIOGRAM)
            vario->setGgByIndex(idir, j,
                         vario->getGgByIndex(idir, j) / vario->getSwByIndex(idir, j));
        }
      }
    }

  // Process the variogram transformations

  if (vario->getCalculType() == ECalcVario::TRANS1)
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < ivar; jvar++)
      {
        for (int i = 0; i < vario->getLagTotalNumber(idir); i++, ecr++)
        {
          int j = vario->getDirAddress(idir, ivar, jvar, i, true, 0);
          int j0 = vario->getDirAddress(idir, jvar, jvar, i, true, 0);
          vario->setGgByIndex(idir, j,
                       -vario->getGgByIndex(idir, j) / vario->getGgByIndex(idir, j0));
        }
      }
  }
  else if (vario->getCalculType() == ECalcVario::TRANS2)
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < ivar; jvar++)
      {
        for (int i = 0; i < vario->getLagTotalNumber(idir); i++, ecr++)
        {
          int j = vario->getDirAddress(idir, ivar, jvar, i, true, 0);
          int j0 = vario->getDirAddress(idir, ivar, ivar, i, true, 0);
          vario->setGgByIndex(idir, j,
                       -vario->getGgByIndex(idir, j) / vario->getGgByIndex(idir, j0));
        }
      }
  }
  else if (vario->getCalculType() == ECalcVario::BINORMAL)
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < ivar; jvar++)
      {
        for (int i = 0; i < vario->getLagTotalNumber(idir); i++, ecr++)
        {
          int j = vario->getDirAddress(idir, ivar, jvar, i, true, 0);
          int j1 = vario->getDirAddress(idir, ivar, ivar, i, true, 0);
          int j2 = vario->getDirAddress(idir, jvar, jvar, i, true, 0);
          vario->setGgByIndex(
              idir,
              j,
              vario->getGgByIndex(idir, j) / sqrt(
                  vario->getGgByIndex(idir, j1) * vario->getGgByIndex(idir, j2)));
        }
      }
  }
  return;
}

/****************************************************************************/
/*!
 **  Center the covariance calculations
 **
 ** \param[in]  db    Db descriptor
 ** \param[in]  vario Vario structure
 ** \param[in]  idir  Rank of the direction
 **
 *****************************************************************************/
static void st_covariance_center(Db *db, Vario *vario, int idir)
{
  int i, j, ivar, jvar, iech;
  double m1, m2, sumw, z1, z2, ww;

  if (!vario->getFlagAsym()) return;

  /* Scale the experimental variogram quantities */

  for (ivar = 0; ivar < vario->getVariableNumber(); ivar++)
    for (jvar = 0; jvar <= ivar; jvar++)
    {
      /* Calculate the mean for each variable */

      m1 = m2 = sumw = 0.;
      for (iech = 0; iech < db->getSampleNumber(); iech++)
      {
        if (!db->isActive(iech)) continue;
        ww = db->getWeight(iech);
        if (FFFF(ww) || ww < 0.) continue;
        z1 = st_get_IVAR(db, iech, ivar);
        z2 = st_get_IVAR(db, iech, jvar);
        if (FFFF(z1) || FFFF(z2)) continue;
        m1 += ww * z1;
        m2 += ww * z2;
        sumw += ww;
      }

      if (sumw > 0 && (vario->getCalculType() == ECalcVario::COVARIANCE
          || vario->getCalculType() == ECalcVario::COVARIANCE_NC))
      {
        m1 /= sumw;
        m2 /= sumw;
      }

      /* Perform the Centering */

      if (!(vario->getCalculType() == ECalcVario::COVARIOGRAM
          || vario->getCalculType() == ECalcVario::COVARIANCE_NC))
        for (i = 0; i < vario->getLagTotalNumber(idir); i++)
        {
          j = vario->getDirAddress(idir, ivar, jvar, i, true, 0);
          if (vario->getSwByIndex(idir, j) > 0)
            vario->setGgByIndex(idir, j, vario->getGgByIndex(idir, j) - m1 * m2);
        }
    }
  return;
}

/****************************************************************************/
/*!
 **  Patch the value of C(0) for covariances
 **
 ** \param[in]  db    Db descriptor
 ** \param[in]  vario Vario structure
 ** \param[in]  idir  Rank of the Direction
 **
 *****************************************************************************/
static void st_variogram_patch_c00(Db *db, Vario *vario, int idir)
{
  int i, ivar, jvar, iech;
  double z1, z2, s12w, s12wzz, ww, scale, value, m1, m2, sumw;

  /* Initializations */

  if (!vario->getFlagAsym()) return;

  /* Calculate the C00 term */

  for (ivar = 0; ivar < db->getVariableNumber(); ivar++)
    for (jvar = 0; jvar <= ivar; jvar++)
    {
      i = vario->getDirAddress(idir, ivar, jvar, 0, false, 0);
      vario->setHhByIndex(idir, i, 0.);

      scale = 1.;
      m1 = m2 = s12w = s12wzz = sumw = 0.;

      /* Calculate the statistics for each variable */

      for (iech = 0; iech < db->getSampleNumber(); iech++)
      {
        if (!db->isActive(iech)) continue;
        ww = db->getWeight(iech);
        if (FFFF(ww) || ww < 0.) continue;
        z1 = st_get_IVAR(db, iech, ivar);
        z2 = st_get_IVAR(db, iech, jvar);
        if (FFFF(z1) || FFFF(z2)) continue;
        m1 += ww * z1;
        m2 += ww * z2;
        sumw += ww;
        value = z1 * z2;
        if (vario->getCalculType() == ECalcVario::COVARIOGRAM)
        {
          scale = ww;
        }
        else
        {
          scale = ww * ww;
          s12w += scale;
        }
        s12wzz += scale * value;
        if (debug_query("variogram"))
          st_print_debug(iech, iech, ivar, jvar, i, scale, value);
      }

      if (sumw > 0 && (vario->getCalculType() == ECalcVario::COVARIANCE
          || vario->getCalculType() == ECalcVario::COVARIANCE_NC))
      {
        m1 /= sumw;
        m2 /= sumw;
      }

      /* Final centering and normation */

      vario->setSwByIndex(idir, i, sumw);
      if (vario->getCalculType() == ECalcVario::COVARIOGRAM)
        vario->setGgByIndex(idir, i, s12wzz);
      else if (vario->getCalculType() == ECalcVario::COVARIANCE_NC)
        vario->setGgByIndex(idir, i, s12wzz / s12w);
      else
        vario->setGgByIndex(idir, i, s12wzz / s12w - m1 * m2);
    }
  return;
}

/****************************************************************************/
/*!
 **  Check if a pair must be kept according to code criterion
 **
 ** \return  1 if the codes are not comparable
 **
 ** \param[in]  db1        First Db structure
 ** \param[in]  db2        Second Db structure
 ** \param[in]  iech       Rank of the first sample
 ** \param[in]  jech       Rank of the second sample
 ** \param[in]  opt_code   code selection option
 ** \li                     0 : no use of the code selection
 ** \li                     1 : codes must be close enough
 ** \li                     2 : codes must be different
 ** \param[in]  tolcode    Code tolerance
 **
 ** \remarks When used in variogram calculation, pairs are discarded then the
 ** \remarks resulting value is 1.
 **
 *****************************************************************************/
int code_comparable(const Db *db1,
                    const Db *db2,
                    int iech,
                    int jech,
                    int opt_code,
                    int tolcode)
{
  double code1, code2;

  /* Dispatch */

  switch (opt_code)
  {
    case 0:
      break;

    case 1: /* Code must be close */
      code1 = db1->getCode(iech);
      code2 = db2->getCode(jech);
      if (ABS(code1 - code2) > tolcode) return (1);
      break;

    case 2: /* Code must be different */
      code1 = db1->getCode(iech);
      code2 = db2->getCode(jech);
      if (code1 == code2) return (1);
      break;
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Check if dates are involved in the variogram calculation
 **
 ** \return  1 if dates are used; 0 otherwise
 **
 ** \param[in]  varioparam VarioParam structure
 ** \param[in]  db1        First Db structure
 ** \param[in]  db2        Second Db structure
 **
 *****************************************************************************/
static int st_date_is_used(const VarioParam *varioparam,
                           const Db *db1,
                           const Db *db2)
{
  if (varioparam->getDates().empty()) return (0);
  if (!db1->hasDate()) return (0);
  if (!db2->hasDate()) return (0);
  return (1);
}

/****************************************************************************/
/*!
 **  Check if a pair must be kept according to date criterion
 **
 ** \return  1 if the dates are not comparable
 **
 ** \param[in]  varioparam VarioParam structure
 ** \param[in]  db1        First Db structure
 ** \param[in]  db2        Second Db structure
 ** \param[in]  iech       Rank of the first sample
 ** \param[in]  jech       Rank of the second sample
 ** \param[in]  idate      Index of the Date interval
 **
 ** \remarks When pairs are discarded then the resulting value is 1.
 **
 *****************************************************************************/
static int st_date_comparable(const VarioParam *varioparam,
                              const Db *db1,
                              const Db *db2,
                              int iech,
                              int jech,
                              int idate)
{
  double date1, date2, delta;

  /* Dispatch */

  if (!varioparam->hasDate()) return (0);
  date1 = db1->getDate(iech);
  date2 = db2->getDate(jech);
  if (FFFF(date1) || FFFF(date2)) return (0);

  delta = date2 - date1;
  if (delta < varioparam->getDate(idate, 0)) return (1);
  if (delta >= varioparam->getDate(idate, 1)) return (1);
  return (0);
}

/****************************************************************************/
/*!
 **  Local function defined as the half-sum of the Variance of Measurement
 **  Error on the two points constituting a pair
 **
 ** \return  Returned value
 **
 ** \param[in]  db     Db description
 ** \param[in]  iech   Rank of the first sample
 ** \param[in]  jech   Rank of the second sample
 **
 *****************************************************************************/
static double st_s(Db *db, int iech, int jech)
{
  double value;

  value = 0.5 * (db->getVarianceError(iech, 0) + db->getVarianceError(jech, 0));
  return (value);
}

/****************************************************************************/
/*!
 **  Local function defined as the half-sum squared difference of values
 **  between the two points constituting a pair
 **
 ** \return  Returned value
 **
 ** \param[in]  db     Db description
 ** \param[in]  iech   Rank of the first sample
 ** \param[in]  jech   Rank of the second sample
 **
 *****************************************************************************/
static double st_g(Db *db, int iech, int jech)
{
  double value;

  value = st_get_IVAR(db, iech, 0) - st_get_IVAR(db, jech, 0);
  value = value * value / 2.;
  return (value);
}

/****************************************************************************/
/*!
 **  Update the Variogram when Variance of measurement error is available
 **
 ** \return  Error returned code
 **
 ** \param[in]  db        Db description
 ** \param[in]  vario     Vario structure
 ** \param[in]  idir      Rank of the direction
 ** \param[in]  vorder    Vario_Order structure
 ** \param[in]  verr_mode Mode of variogram correction (1, 2 or 3)
 **
 *****************************************************************************/
static int st_update_variogram_verr(Db *db,
                                    Vario *vario,
                                    int idir,
                                    Vario_Order *vorder,
                                    int verr_mode)
{
  int ipas, npas, ipair, ifirst, ilast, iech, jech, number, nfois;
  double dist, value, g_old, diff, sumt, sumb, wgt, sval, gval;
  static double tol = 1.e-05;
  static int maxiter = 100;

  /* Initializations */

  npas = vario->getLagNumber(idir);

  /* Loop on the lags */

  for (ipas = 0; ipas < npas; ipas++)
  {
    vario_order_get_bounds(vorder, idir, ipas, &ifirst, &ilast);
    if (ifirst > ilast) continue;

    /* Dispatch according to the method */

    switch (verr_mode)
    {

      case 1: /* Simple bias correction */
        number = 0;
        value = 0.;
        for (ipair = ifirst; ipair < ilast; ipair++)
        {
          vario_order_get_indices(vorder, ipair, &iech, &jech, &dist);
          value += st_s(db, iech, jech);
          number++;
        }
        value = (number > 0) ? value / number :
                               0.;
        vario->setGg(idir, 0, 0, ipas,
                     MAX(0, vario->getGg(idir, 0, 0, ipas) - value));
        break;

      case 2:
        nfois = 0;
        diff = 1.e30;
        while (diff > tol && nfois < maxiter)
        {
          nfois++;
          g_old = vario->getGg(idir, 0, 0, ipas);
          sumt = sumb = 0.;
          for (ipair = ifirst; ipair < ilast; ipair++)
          {
            vario_order_get_indices(vorder, ipair, &iech, &jech, &dist);
            sval = st_s(db, iech, jech);
            gval = st_g(db, iech, jech);
            value = sval + vario->getGg(idir, 0, 0, ipas);
            wgt = 1. / (value * value);
            sumt += wgt * (gval - sval);
            sumb += wgt;
          }
          vario->setGg(idir, 0, 0, ipas, sumt / sumb);
          diff = ABS(vario->getGg(idir, 0, 0, ipas) - g_old);
        }
        vario->setGg(idir, 0, 0, ipas, MAX(0, vario->getGg(idir, 0, 0, ipas)));
        if (nfois == maxiter && debug_query("converge"))
          message("Convergence not reached for lag %d\n", ipas + 1);
        break;

      case 3:
        nfois = 0;
        diff = 1.e30;
        while (diff > tol && nfois < maxiter)
        {
          g_old = vario->getGg(idir, 0, 0, ipas);
          sumt = sumb = 0.;
          for (ipair = ifirst; ipair < ilast; ipair++)
          {
            vario_order_get_indices(vorder, ipair, &iech, &jech, &dist);
            sval = st_s(db, iech, jech);
            gval = st_g(db, iech, jech);
            value = sval + vario->getGgByIndex(idir, ipas);
            wgt = vario->getGg(idir, 0, 0, ipas) / value;
            sumt += wgt * gval;
            sumb += 1.;
          }
          vario->setGg(idir, 0, 0, ipas, sumt / sumb);
          diff = ABS(vario->getGg(idir, 0, 0, ipas) - g_old);
        }
        if (nfois == maxiter && debug_query("converge"))
          message("Convergence not reached for lag %d\n", ipas + 1);
        break;

      default:
        messerr("The method (%d) for updating the Variogram", verr_mode);
        messerr("calculated with Variance of Measurement Error is unknown");
        return (1);
    }
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Return the relative sample rank corresponding to an abolsute number
 **
 ** \return The corresponding relative sample rank
 **
 ** \param[in]  db     Db structure
 ** \param[in]  iech0  Absolute rank of the first sample
 **
 *****************************************************************************/
static int st_get_relative_sample_rank(Db *db, int iech0)
{
  int iech, iiech;

  for (iech = iiech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActiveAndDefined(iech, 0)) continue;
    if (iech == iech0) return (iiech);
    iiech++;
  }
  return (-1);
}

/****************************************************************************/
/*!
 **  Calculate the a bias term between samples iech and jech
 **
 ** \param[in]  db    Db structure
 ** \param[in]  nbfl  Number of drift terms
 ** \param[in]  iiech Relative rank of the first sample
 ** \param[in]  jjech Relative rank of the second sample
 **
 *****************************************************************************/
static double st_get_bias_value(Db *db, int nbfl, int iiech, int jjech)
{
  double bias0, bias1, bias2;
  int il, jl, nech;

  /* Get the relative sample ranks */

  nech = db->getActiveAndDefinedNumber(0);
  bias0 = bias1 = bias2 = 0.;

  for (il = 0; il < nbfl; il++)
    for (jl = 0; jl < nbfl; jl++)
      bias0 += X_DRFXA(il,iiech) * X_DRFXGX(il, jl) * X_DRFXA(jl, jjech);

  for (il = 0; il < nbfl; il++)
    bias1 += X_DRFXA(il,iiech) * X_DRFGX(il, jjech);

  for (il = 0; il < nbfl; il++)
    bias2 += X_DRFGX(il,iiech) * X_DRFXA(il, jjech);

  return (bias0 - (bias1 + bias2));
}

/****************************************************************************/
/*!
 **  Calculate the local bias terms
 **
 ** \param[in]  db        Db description
 ** \param[in]  vorder    Vario_Order structure
 ** \param[in]  ifirst    Rank of the first lag
 ** \param[in]  ilast     Rank of the last lag
 **
 *****************************************************************************/
static double st_calculate_bias_local(Db *db,
                                      Vario_Order *vorder,
                                      int ifirst,
                                      int ilast)
{
  int ipair, iech, jech, nbfl, iiech, jjech;
  double dist, diff, tot0, tot1, tot2, totnum, result, v1, v2;

  /* Initializations */

  nbfl = MODEL->getDriftNumber();

  /* Calculate the first corrected term */

  tot0 = tot1 = tot2 = totnum = 0.;
  for (ipair = ifirst; ipair < ilast; ipair++)
  {
    vario_order_get_indices(vorder, ipair, &iech, &jech, &dist);
    v1 = st_get_IVAR(db, iech, 0);
    v2 = st_get_IVAR(db, jech, 0);
    if (FFFF(v1) || FFFF(v2)) continue;

    iiech = st_get_relative_sample_rank(db, iech);
    jjech = st_get_relative_sample_rank(db, jech);

    diff = v1 - v2;
    tot0 += diff * diff;
    tot1 += st_get_bias_value(db, nbfl, iiech, jjech);
    tot2 += (DRFDIAG[iiech] + DRFDIAG[jjech]) / 2.;
    totnum += 1.;
  }

  tot0 /= 2.;
  result = tot0 - tot1 + tot2;
  if (totnum > 0.) result /= totnum;

  return (result);
}

/****************************************************************************/
/*!
 **  Calculate the global bias terms
 **
 ** \param[in]  db        Db description
 ** \param[in]  d1        Working vector (Dimension: ndim)
 **
 *****************************************************************************/
static void st_calculate_bias_global(Db *db, VectorDouble d1)
{
  double c00, covtab, value;
  int idim, ndim, il, jl, nech, iech, iiech, jech, jjech, nbfl;
  CovCalcMode mode;

  /* Initializations */

  nbfl = MODEL->getDriftNumber();
  ndim = MODEL->getDimensionNumber();
  nech = db->getActiveAndDefinedNumber(0);

  /* Calculate the c00 term */

  for (idim = 0; idim < ndim; idim++)
    d1[idim] = 0.;
  model_calcul_cov(MODEL, mode, 1, 1., d1, &c00);

  /* Calculate the term: G %*% X */

  for (iech = iiech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActiveAndDefined(iech, 0)) continue;
    for (il = 0; il < nbfl; il++)
    {
      value = 0;
      for (jech = jjech = 0; jech < db->getSampleNumber(); jech++)
      {
        if (!db->isActiveAndDefined(jech, 0)) continue;
        for (idim = 0; idim < ndim; idim++)
          d1[idim] = db->getDistance1D(iech, jech, idim);
        model_calcul_cov(MODEL, mode, 1, 1., d1, &covtab);
        value += (c00 - covtab) * X_DRFTAB(il, jjech);
        jjech++;
      }
      X_DRFGX(il,iiech) = value;
    }
    iiech++;
  }

  /* Calculate the term: t(X) %*% G %*% X */

  for (il = 0; il < nbfl; il++)
    for (jl = 0; jl < nbfl; jl++)
    {
      value = 0;
      for (iech = 0; iech < nech; iech++)
        value += X_DRFGX(il,iech) * X_DRFTAB(jl, iech);
      X_DRFXGX(il,jl) = value;
    }

  /* Calculate the term: diag(bias) */

  for (iech = iiech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActiveAndDefined(iech, 0)) continue;
    DRFDIAG[iiech] = st_get_bias_value(db, nbfl, iiech, iiech);
    iiech++;
  }
}

/****************************************************************************/
/*!
 **  Update the Variogram of Residuals when Drift has been removed
 **
 ** \return  Error returned code
 **
 ** \param[in]  db         Db description
 ** \param[in]  vario      Vario structure
 ** \param[in]  vorder     Vario_Order structure
 ** \param[in]  verbose    Verbose flag
 **
 ** \remark The number of iterations used for the debiasing procedure
 ** \remark in presence of a drift can be defined using:
 ** \remark set_keypair("KU_Niter",newval)
 **
 *****************************************************************************/
static int st_update_variogram_ku(Db *db,
                                  Vario *vario,
                                  Vario_Order *vorder,
                                  int verbose)
{
  Option_VarioFit optvar;
  Option_AutoFit mauto;
  double newval;
  int error, iter, ndim, idir, ipas, ifirst, ilast, niter_ku, nbfl;
  VectorDouble d1;
  Constraints constraints;

  /* Initializations */

  error = 1;
  ndim = MODEL->getDimensionNumber();
  nbfl = MODEL->getDriftNumber();
  niter_ku = (int) get_keypone("KU_Niter", 0);

  /* Core allocation */

  d1.resize(ndim);

  /* Loop on the iterations */

  for (iter = 0; iter < niter_ku; iter++)
  {

    /* Perform the Automatic structure recognition */

    if (MODEL != nullptr)
    {
      if (model_auto_fit(vario, MODEL, verbose, mauto, constraints, optvar))
        goto label_end;
    }

    /* Calculate the global bias correction terms */

    st_calculate_bias_global(db, d1);

    /* Optional printout */

    if (verbose)
    {
      message("Drift removal at iteration #d/%d\n", iter + 1, niter_ku);
      print_matrix("Drift Coefficients Matrix", 0, 1, nbfl, nbfl, NULL, DRFXGX);
    }

    /* Loop on the directions */

    for (idir = 0; idir < vario->getDirectionNumber(); idir++)
    {

      /* Loop on the lags */

      for (ipas = 0; ipas < vario->getLagNumber(idir); ipas++)
      {
        vario_order_get_bounds(vorder, idir, ipas, &ifirst, &ilast);
        if (ifirst > ilast) continue;

        /* Calculate the local bias correction terms */

        newval = st_calculate_bias_local(db, vorder, ifirst, ilast);

        /* Patch the new value */

        vario->setGg(idir, 0, 0, ipas, newval);
      }
    }
  }
  error = 0;

  label_end: return (error);
}

/****************************************************************************/
/*!
 **  Perform the Variogram evaluation when the Geometry has been established
 **
 ** \param[in]  db     Db description
 ** \param[in]  vario  Vario structure
 ** \param[in]  idir   Rank of the direction
 ** \param[in]  vorder Vario_Order structure
 **
 *****************************************************************************/
static void st_variogram_calcul_internal(Db *db,
                                         Vario *vario,
                                         int idir,
                                         Vario_Order *vorder)
{
  int iech, jech, ipas, npas, ifirst, ilast, ipair;
  double dist;

  /* Initializations */

  npas = vario->getLagNumber(idir);

  /* Loop on the lags */

  for (ipas = 0; ipas < npas; ipas++)
  {
    vario_order_get_bounds(vorder, idir, ipas, &ifirst, &ilast);

    /* Loop on the pairs contributing to this lag */

    for (ipair = ifirst; ipair < ilast; ipair++)
    {
      vario_order_get_indices(vorder, ipair, &iech, &jech, &dist);

      /* Evaluate the variogram */

      VARIO = vario;
      IDIRLOC = idir;
      IECH1 = iech;
      IECH2 = jech;
      st_variogram_evaluate(db, vario->getCalculType(),
                            vario->getVariableNumber(), iech, jech, ipas, dist,
                            1, st_variogram_set);
    }
  }

  /* Scale the variogram calculations */

  variogram_scale(vario, idir);

  /* Center the covariance function */

  st_covariance_center(db, vario, idir);

  /* Patch the central value */

  st_variogram_patch_c00(db, vario, idir);

  return;
}

/****************************************************************************/
/*!
 **  Check if a pair must be rejected or not
 **
 ** \return  1 if the pair must be rejected
 **
 ** \param[in]  db     Db description
 ** \param[in]  iech   Rank of the first sample
 ** \param[in]  jech   Rank of the second sample
 ** \param[in]  dist   Distance between the two samples
 ** \param[in]  psmin  Direction cosine
 ** \param[in]  bench  Slicing bench
 ** \param[in]  cylrad Slicing radius
 ** \param[in]  codir  Direction
 **
 ** \param[out] ps     The cosine between the vector of data and the direction
 **
 *****************************************************************************/
int variogram_reject_pair(const Db *db,
                          int iech,
                          int jech,
                          double dist,
                          double psmin,
                          double bench,
                          double cylrad,
                          const VectorDouble &codir,
                          double *ps)
{
  /* When the distance is zero, the pair is always accepted */

  if (dist <= 0.) return (0);

  /* Calculate the cosine between the pairs and the direction */

  *ps = db->getCosineToDirection(iech, jech, codir);

  /* Angular tolerance */

  if (ABS(*ps) < psmin) return (1);

  /* Vertical slicing */

  if (!FFFF(bench) && bench > 0.)
  {
    if (bench_distance(db, iech, jech) > bench) return (1);
  }

  /* Cylinder rejection */

  if (!FFFF(cylrad) && cylrad > 0.)
  {
    if (cylinder_radius(db, iech, jech, codir) > cylrad) return (1);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Evaluate the variogram using the traditional method
 **
 ** \return  Error return code
 **
 ** \param[in]  db     Db description
 ** \param[in]  vario  Vario structure
 ** \param[in]  idir   Rank of the direction
 ** \param[in]  rindex Array of sorted samples
 ** \param[in]  vorder Vario_Order structure
 **
 *****************************************************************************/
static int st_variogram_calcul1(Db *db,
                                Vario *vario,
                                int idir,
                                int *rindex,
                                Vario_Order *vorder)
{
  int iiech, iech, jjech, jech, nech, ipas, npair, ideb;
  double psmin, ps, dist, maxdist;

  ps = 0.;
  psmin = _variogram_convert_angular_tolerance(
      vario->getDirParam(idir).getTolAngle());
  nech = db->getSampleNumber();
  maxdist = vario->getMaximumDistance(idir);
  const VarioParam &varioparam = vario->getVarioParam();

  /* Loop on the first point */

  for (iiech = 0; iiech < nech - 1; iiech++)
  {
    iech = rindex[iiech];
    if (!db->isActive(iech)) continue;
    if (FFFF(db->getWeight(iech))) continue;

    ideb = (st_date_is_used(&varioparam, db, db)) ? 0 :
                                                    iiech + 1;
    for (jjech = ideb; jjech < nech; jjech++)
    {
      jech = rindex[jjech];
      if (variogram_maximum_dist1D_reached(db, iech, jech, maxdist)) break;
      if (!db->isActive(jech)) continue;
      if (FFFF(db->getWeight(jech))) continue;

      /* Check if the pair must be kept (Code criterion) */

      if (code_comparable(db, db, iech, jech,
                          vario->getDirParam(idir).getOptionCode(),
                          (int) vario->getDirParam(idir).getTolCode()))
        continue;

      /* Check if the pair must be kept (Date criterion) */

      if (st_date_comparable(&varioparam, db, db, iech, jech,
                             vario->getIdate(idir))) continue;

      /* Check if the pair must be kept */

      dist = distance_intra(db, iech, jech, NULL);
      if (variogram_reject_pair(db, iech, jech, dist, psmin,
                                vario->getDirParam(idir).getBench(),
                                vario->getDirParam(idir).getCylRad(),
                                vario->getDirParam(idir).getCodir(), &ps))
        continue;

      /* Get the rank of the lag */

      ipas = variogram_get_lag(vario, idir, ps, psmin, &dist);
      if (IFFFF(ipas)) continue;

      /* Case of internal storage */

      if (vorder != (Vario_Order*) NULL)
        vario_order_add(vorder, iech, jech, NULL, NULL, ipas, idir, dist);
      else
      {

        /* Evaluate the variogram */

        VARIO = vario;
        IDIRLOC = idir;
        IECH1 = iech;
        IECH2 = jech;
        st_variogram_evaluate(db, vario->getCalculType(),
                              vario->getVariableNumber(), iech, jech, ipas,
                              dist, 1, st_variogram_set);
      }
    }
  }

  /* Internal storage */

  if (vorder != (Vario_Order*) NULL)
  {
    vorder = vario_order_final(vorder, &npair);
  }
  else
  {

    /* Scale the variogram calculations */

    variogram_scale(vario, idir);

    /* Center the covariance function */

    st_covariance_center(db, vario, idir);

    /* Patch the central value */

    st_variogram_patch_c00(db, vario, idir);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Evaluate the variogram by sample
 **
 ** \return  Error return code
 **
 ** \param[in]  db     Db descriptor
 ** \param[in]  vario  Vario structure
 ** \param[in]  idir   Rank of the direction
 ** \param[in]  rindex Array of sorted samples
 **
 *****************************************************************************/
static int st_variogram_calcul2(Db *db, Vario *vario, int idir, int *rindex)
{
  int iiech, iech, jjech, jech, i, ipas, nech, size, error, ideb;
  double *gg_sum, *hh_sum, *sw_sum, ps, psmin, dist, w1, maxdist;

  /* Initializations */

  error = 1;
  ps = 0.;
  gg_sum = hh_sum = sw_sum = nullptr;
  const VarioParam &varioparam = vario->getVarioParam();
  const DirParam &dirparam = vario->getDirParam(idir);
  psmin = _variogram_convert_angular_tolerance(dirparam.getTolAngle());
  nech = db->getSampleNumber();
  size = vario->getDirSize(idir);
  maxdist = vario->getMaximumDistance(idir);

  /* Core allocation */

  gg_sum = (double*) mem_alloc(size * sizeof(double), 0);
  if (gg_sum == nullptr) goto label_end;
  hh_sum = (double*) mem_alloc(size * sizeof(double), 0);
  if (hh_sum == nullptr) goto label_end;
  sw_sum = (double*) mem_alloc(size * sizeof(double), 0);
  if (sw_sum == nullptr) goto label_end;
  for (i = 0; i < size; i++)
  {
    gg_sum[i] = 0.;
    hh_sum[i] = 0.;
    sw_sum[i] = 0.;
  }

  /* Loop on the first sample */

  for (iiech = 0; iiech < nech; iiech++)
  {
    iech = rindex[iiech];
    if (!db->isActive(iech)) continue;
    w1 = db->getWeight(iech);
    if (FFFF(w1)) continue;

    /* Looking for the second sample */

    ideb = (st_date_is_used(&varioparam, db, db)) ? 0 :
                                                    iiech + 1;
    for (jjech = ideb; jjech < nech; jjech++)
    {
      jech = rindex[jjech];
      if (variogram_maximum_dist1D_reached(db, iech, jech, maxdist)) break;
      if (!db->isActive(jech)) continue;
      if (FFFF(db->getWeight(jech))) continue;

      /* Check if the pair must be kept (Code criterion) */

      if (code_comparable(db, db, iech, jech, dirparam.getOptionCode(),
                          (int) dirparam.getTolCode())) continue;

      /* Check if the pair must be kept (Date criterion) */

      if (st_date_comparable(&varioparam, db, db, iech, jech,
                             dirparam.getIdate())) continue;

      /* Check if the pair must be kept */

      dist = distance_intra(db, iech, jech, NULL);
      if (variogram_reject_pair(db, iech, jech, dist, psmin,
                                dirparam.getBench(), dirparam.getCylRad(),
                                dirparam.getCodir(), &ps)) continue;

      /* Get the rank of the lag */

      ipas = variogram_get_lag(vario, idir, ps, psmin, &dist);
      if (IFFFF(ipas)) continue;

      /* Evaluate the variogram */

      VARIO = vario;
      IECH1 = iech;
      IECH2 = jech;
      st_variogram_evaluate(db, vario->getCalculType(),
                            vario->getVariableNumber(), iech, jech, ipas, dist,
                            1, st_variogram_set);
    }

    /* Cumulate to the global variogram */

    for (i = 0; i < size; i++)
    {
      if (vario->getSwByIndex(idir, i) <= 0) continue;
      sw_sum[i] += w1;
      gg_sum[i] += w1 * vario->getGgByIndex(idir, i) / vario->getSwByIndex(idir, i);
      hh_sum[i] += w1 * vario->getHhByIndex(idir, i) / vario->getSwByIndex(idir, i);
    }
  }

  /* Copy the cumulated variogram into the Vario structure */

  for (i = 0; i < size; i++)
  {
    vario->setGgByIndex(idir, i, gg_sum[i]);
    vario->setHhByIndex(idir, i, hh_sum[i]);
    vario->setSwByIndex(idir, i, sw_sum[i]);
  }

  /* Scale the variogram calculations */

  variogram_scale(vario, idir);

  /* Center the covariance function */

  st_covariance_center(db, vario, idir);

  /* Patch the central value */

  st_variogram_patch_c00(db, vario, idir);

  /* Set the error return flag */

  error = 0;

  label_end: gg_sum = (double*) mem_free((char* ) gg_sum);
  hh_sum = (double*) mem_free((char* ) hh_sum);
  sw_sum = (double*) mem_free((char* ) sw_sum);
  return (error);
}

/****************************************************************************/
/*!
 **  Evaluate the covariance for directional variables
 **
 ** \return  Error return code
 **
 ** \param[in]  db     Db description
 ** \param[in]  vario  Vario structure
 ** \param[in]  idir   Rank of the Direction
 ** \param[in]  ncomp  Number of components
 ** \param[in]  rindex Array of sorted samples
 **
 *****************************************************************************/
static int st_variovect_calcul(Db *db,
                               Vario *vario,
                               int idir,
                               int ncomp,
                               int *rindex)
{
  int iiech, iech, jjech, jech, nech, ipas, ivar, jvar, i, icomp, nvar;
  double psmin, ps, dist, w1, w2, zi1, zi2, zj1, zj2, v12, v21;
  double di1, di2, dj1, dj2, maxdist;

  const DirParam &dirparam = vario->getDirParam(idir);
  ps = 0.;
  psmin = _variogram_convert_angular_tolerance(dirparam.getTolAngle());
  nech = db->getSampleNumber();
  nvar = vario->getVariableNumber();
  maxdist = vario->getMaximumDistance(idir);

  /* Loop on the first point */

  for (iiech = 0; iiech < nech - 1; iiech++)
  {
    iech = rindex[iiech];
    if (!db->isActive(iech)) continue;
    if (FFFF(db->getWeight(iech))) continue;

    for (jjech = iiech + 1; jjech < nech; jjech++)
    {
      jech = rindex[jjech];
      if (variogram_maximum_dist1D_reached(db, iech, jech, maxdist)) break;
      if (!db->isActive(jech)) continue;
      if (FFFF(db->getWeight(jech))) continue;

      /* Check if the pair must be kept (distance criterion) */

      dist = distance_intra(db, iech, jech, NULL);
      if (variogram_reject_pair(db, iech, jech, dist, psmin,
                                dirparam.getBench(), dirparam.getCylRad(),
                                dirparam.getCodir(), &ps)) continue;

      /* Get the rank of the lag */

      ipas = variogram_get_lag(vario, idir, ps, psmin, &dist);
      if (IFFFF(ipas)) continue;

      w1 = db->getWeight(iech);
      w2 = db->getWeight(jech);

      for (ivar = 0; ivar < nvar; ivar++)
        for (jvar = 0; jvar <= ivar; jvar++)
        {

          /* Evaluate the variogram */

          v12 = v21 = di1 = di2 = dj1 = dj2 = 0.;
          for (icomp = 0; icomp < ncomp; icomp++)
          {
            zi1 = st_get_IVAR(db, iech, ivar * ncomp + icomp);
            zi2 = st_get_IVAR(db, iech, jvar * ncomp + icomp);
            zj1 = st_get_IVAR(db, jech, ivar * ncomp + icomp);
            zj2 = st_get_IVAR(db, jech, jvar * ncomp + icomp);
            if (FFFF(zi1) || FFFF(zi2) || FFFF(zj1) || FFFF(zj2))
            {
              v12 = v21 = TEST;
              break;
            }
            else
            {
              v12 += zi1 * zj2;
              v21 += zi2 * zj1;
              di1 += zi1 * zi1;
              di2 += zi2 * zi2;
              dj1 += zj1 * zj1;
              dj2 += zj2 * zj2;
            }
          }
          if (FFFF(v12) || FFFF(v21)) continue;
          if (ABS(di1) < EPSILON8 || ABS(di2) < EPSILON8) continue;
          if (ABS(dj1) < EPSILON8 || ABS(dj2) < EPSILON8) continue;
          di1 = sqrt(di1);
          di2 = sqrt(di2);
          dj1 = sqrt(dj1);
          dj2 = sqrt(dj2);
          v12 = ABS(v12) / (di1 * dj2);
          v21 = ABS(v21) / (di2 * dj1);

          i = vario->getDirAddress(idir, ivar, jvar, ipas, false, 1);
          vario->setGgByIndex(idir, i, vario->getGgByIndex(idir, i) + w1 * w2 * v12);
          vario->setHhByIndex(idir, i, vario->getHhByIndex(idir, i) + w1 * w2 * dist);
          vario->setSwByIndex(idir, i, vario->getSwByIndex(idir, i) + w1 * w2);

          i = vario->getDirAddress(idir, ivar, jvar, ipas, false, -1);
          vario->setGgByIndex(idir, i, vario->getGgByIndex(idir, i) + w1 * w2 * v21);
          vario->setHhByIndex(idir, i, vario->getHhByIndex(idir, i) + w1 * w2 * dist);
          vario->setSwByIndex(idir, i, vario->getSwByIndex(idir, i) + w1 * w2);
        }
    }
  }

  /* Scale the variogram calculations */

  variogram_scale(vario, idir);

  /* Center the covariance function */

  st_covariance_center(db, vario, idir);

  /* Patch the central value */

  st_variogram_patch_c00(db, vario, idir);

  return (0);
}

/****************************************************************************/
/*!
 **  Evaluate the variogram on a grid
 **
 ** \return  Error return code
 **
 ** \param[in]  db    Db description
 ** \param[in]  vario Vario structure
 ** \param[in]  idir  Rank of the Direction
 **
 *****************************************************************************/
static int st_variogram_grid(Db *db, Vario *vario, int idir)
{
  int *indg1, *indg2, iech, jech, nech, ipas, idim, error, npas;
  double dist;

  /* Initializations */

  error = 1;
  nech = db->getSampleNumber();
  indg1 = indg2 = nullptr;
  npas = vario->getLagNumber(idir);
  const DirParam &dirparam = vario->getDirParam(idir);
  const VarioParam &varioparam = vario->getVarioParam();

  /* Core allocation */

  indg1 = db_indg_alloc(db);
  if (indg1 == nullptr) goto label_end;
  indg2 = db_indg_alloc(db);
  if (indg2 == nullptr) goto label_end;

  /* Loop on the first point */

  for (iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    if (FFFF(db->getWeight(iech))) continue;
    db_index_sample_to_grid(db, iech, indg1);

    for (ipas = 1; ipas < npas; ipas++)
    {
      for (idim = 0; idim < db->getNDim(); idim++)
        indg2[idim] = indg1[idim] + (int) (ipas * vario->getGrincr(idir, idim));
      jech = db_index_grid_to_sample(db, indg2);
      if (jech < 0) continue;

      if (!db->isActive(jech)) continue;
      if (FFFF(db->getWeight(jech))) continue;

      /* Check if the pair must be kept (Code criterion) */

      if (code_comparable(db, db, iech, jech, dirparam.getOptionCode(),
                          (int) dirparam.getTolCode())) continue;

      /* Check if the pair must be kept (Date criterion) */

      if (st_date_comparable(&varioparam, db, db, iech, jech,
                             dirparam.getIdate())) continue;

      /* Evaluate the variogram */

      dist = ipas * dirparam.getDPas();
      VARIO = vario;
      IDIRLOC = idir;
      IECH1 = iech;
      IECH2 = jech;
      st_variogram_evaluate(db, vario->getCalculType(),
                            vario->getVariableNumber(), iech, jech, ipas, dist,
                            1, st_variogram_set);
    }
  }

  /* Scale the variogram calculations */

  variogram_scale(vario, idir);

  /* Center the covariance function */

  st_covariance_center(db, vario, idir);

  /* Patch the central value */

  st_variogram_patch_c00(db, vario, idir);

  /* Set the error return status */

  error = 0;

  label_end:

  /* Core deallocation */

  indg1 = db_indg_free(indg1);
  indg2 = db_indg_free(indg2);
  return (error);
}

/****************************************************************************/
/*!
 **  Evaluate the generalized variogram along lines
 **
 ** \param[in]  db      Db description
 ** \param[in]  vario   Vario structure
 ** \param[in]  idir    Rank of the Direction
 ** \param[in]  norder  Order of the generalized variogram
 **
 *****************************************************************************/
static void st_variogen_line(Db *db, Vario *vario, int idir, int norder)
{
  int iech, jech, nech, ipas, npas, iwgt, keep, nvar;
  double dist, dist0, value, zz, psmin, ps;

  const DirParam &dirparam = vario->getDirParam(idir);
  psmin = _variogram_convert_angular_tolerance(dirparam.getTolAngle());
  nech = db->getSampleNumber();
  npas = vario->getLagNumber(idir);

  /* Loop on the first point */

  for (iech = 0; iech < nech - 1; iech++)
  {
    if (!db->isActive(iech)) continue;

    for (ipas = 1; ipas < npas; ipas++)
    {
      value = st_get_IVAR(db, iech, 0);
      if (FFFF(value)) break;
      dist0 = 0.;

      for (iwgt = keep = 1; iwgt < NWGT[norder] && keep; iwgt++)
      {
        keep = 0;
        jech = iech + iwgt * ipas;
        if (jech < 0 || jech > nech) break;
        if (!db->isActive(jech)) break;

        /* Calculate the distance */

        dist = distance_intra(db, iech, jech, NULL);
        if (variogram_reject_pair(db, iech, jech, dist, psmin,
                                  dirparam.getBench(), dirparam.getCylRad(),
                                  dirparam.getCodir(), &ps)) continue;
        if (iwgt == 1) dist0 = dist;

        /* Check if the next sample belongs to the same line (same code) */

        if (code_comparable(db, db, iech, jech, 1, 0)) break;

        /* Evaluate the variogram */

        zz = st_get_IVAR(db, jech, 0);
        if (FFFF(zz)) break;
        keep = 1;
        value += zz * VARWGT[norder][iwgt];
      }
      if (keep)
      {
        value = value * value / NORWGT[norder];
        VARIO = vario;
        nvar = vario->getVariableNumber();
        st_variogram_set(vario->getCalculType(), nvar, ipas, 0, 0, 0, 1., dist0,
                         value);
      }
    }
  }

  /* Scale the variogram calculations */

  variogram_scale(vario, idir);

  /* Center the covariance function */

  st_covariance_center(db, vario, idir);

  /* Patch the central value */

  st_variogram_patch_c00(db, vario, idir);

  return;
}

/****************************************************************************/
/*!
 **  Evaluate the generalized variogram on a grid
 **
 ** \return  Error return code
 **
 ** \param[in]  db      Db description
 ** \param[in]  vario   Vario structure
 ** \param[in]  idir    Rank of the direction
 ** \param[in]  norder  Order of the generalized variogram
 **
 *****************************************************************************/
static int st_variogen_grid(Db *db, Vario *vario, int idir, int norder)
{
  int *indg1, *indg2;
  int iwgt, iech, jech, nech, ipas, idim, error, npas, keep, nvar;
  double dist, zz, value;

  /* Initializations */

  error = 1;
  nech = db->getSampleNumber();
  indg1 = indg2 = nullptr;
  npas = vario->getLagNumber(idir);
  const DirParam &dirparam = vario->getDirParam(idir);
  const VarioParam &varioparam = vario->getVarioParam();

  /* Core allocation */

  indg1 = db_indg_alloc(db);
  if (indg1 == nullptr) goto label_end;
  indg2 = db_indg_alloc(db);
  if (indg2 == nullptr) goto label_end;

  /* Loop on the first point */

  for (iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    db_index_sample_to_grid(db, iech, indg1);

    for (ipas = 1; ipas < npas; ipas++)
    {
      value = st_get_IVAR(db, iech, 0);
      if (FFFF(value)) break;
      dist = ipas * vario->getDPas(idir);

      for (iwgt = keep = 1; iwgt < NWGT[norder] && keep; iwgt++)
      {
        keep = 0;
        for (idim = 0; idim < db->getNDim(); idim++)
          indg2[idim] = indg1[idim]
              + (int) (ipas * iwgt * dirparam.getGrincr(idim));

        jech = db_index_grid_to_sample(db, indg2);
        if (jech < 0) continue;
        if (!db->isActive(jech)) continue;

        /* Check if the pair must be kept (Code criterion) */

        if (code_comparable(db, db, iech, jech, dirparam.getOptionCode(),
                            (int) dirparam.getTolCode())) continue;

        /* Check if the pair must be kept (Date criterion) */

        if (st_date_comparable(&varioparam, db, db, iech, jech,
                               dirparam.getIdate())) continue;

        /* Evaluate the variogram */

        zz = st_get_IVAR(db, jech, 0);
        if (FFFF(zz)) break;
        keep = 1;
        value += zz * VARWGT[norder][iwgt];
      }
      if (keep)
      {
        value = value * value / NORWGT[norder];
        VARIO = vario;
        nvar = vario->getVariableNumber();
        st_variogram_set(vario->getCalculType(), nvar, ipas, 0, 0, 0, 1., dist,
                         value);
      }
    }
  }

  /* Scale the variogram calculations */

  variogram_scale(vario, idir);

  /* Center the covariance function */

  st_covariance_center(db, vario, idir);

  /* Patch the central value */

  st_variogram_patch_c00(db, vario, idir);

  /* Set the error return status */

  error = 0;

  label_end:

  /* Core deallocation */

  indg1 = db_indg_free(indg1);
  indg2 = db_indg_free(indg2);
  return (error);
}

/****************************************************************************/
/*!
 **  Evaluate the experimental generalized variogram on grid
 **
 ** \return  Error return code
 **
 ** \param[in]  db         Db descriptor
 ** \param[in]  vario      Vario structure
 **
 *****************************************************************************/
static int st_variogen_grid_calcul(Db *db, Vario *vario)
{
  int idir, error, norder;

  /* Initializations */

  error = 1;
  if (db == nullptr) return (1);
  if (vario == nullptr) return (1);
  norder = st_get_generalized_variogram_order(vario);
  st_manage_drift_removal(0, NULL, NULL);

  /* Preliminary checks */

  if (db->getNDim() != vario->getDimensionNumber() || db->getVariableNumber()
      != vario->getVariableNumber())
  {
    messerr("Inconsistent parameters:");
    messerr("Data Base: NDIM=%d NVAR=%d", db->getNDim(),
            db->getVariableNumber());
    messerr("Variogram: NDIM=%d NVAR=%d", vario->getDimensionNumber(),
            vario->getVariableNumber());
    return (1);
  }
  if (vario->getVariableNumber() != 1)
  {
    messerr("The generalized variogram requires a single variable");
    return (1);
  }
  if (st_get_generalized_variogram_order(vario) == 0)
  {
    messerr("This calculation requires a generalized variogram definition");
    return (1);
  }
  if (!is_grid(db))
  {
    messerr("This calculation facility is dedicated to grid architecture");
    return (1);
  }

  /* Update the global statistics */

  st_variogram_stats(db, vario);

  /* Loop on the directions to evaluate */

  for (idir = 0; idir < vario->getDirectionNumber(); idir++)
  {
    error = st_variogen_grid(db, vario, idir, norder);
    if (error) break;
  }

  return (error);
}

/****************************************************************************/
/*!
 **  Evaluate the experimental generalized variogram along lines
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  vario        Vario structure
 **
 *****************************************************************************/
static int st_variogen_line_calcul(Db *db, Vario *vario)
{
  int idir, error, norder;

  /* Initializations */

  error = 1;
  if (db == nullptr) return (1);
  if (vario == nullptr) return (1);
  norder = st_get_generalized_variogram_order(vario);
  st_manage_drift_removal(0, NULL, NULL);

  /* Preliminary checks */

  if (db->getNDim() != vario->getDimensionNumber() || db->getVariableNumber()
      != vario->getVariableNumber())
  {
    messerr("Inconsistent parameters:");
    messerr("Data Base: NDIM=%d NVAR=%d", db->getNDim(),
            db->getVariableNumber());
    messerr("Variogram: NDIM=%d NVAR=%d", vario->getDimensionNumber(),
            vario->getVariableNumber());
    return (1);
  }
  if (vario->getVariableNumber() != 1)
  {
    messerr("The generalized variogram requires a single variable");
    return (1);
  }
  if (st_get_generalized_variogram_order(vario) == 0)
  {
    messerr("This calculation requires a generalized variogram definition");
    return (1);
  }
  if (!is_grid(db))
  {
    messerr("This calculation facility is dedicated to line architecture");
    return (1);
  }
  if (!db->hasCode())
  {
    messerr("This calculation facility requires the definition of a CODE");
    return (1);
  }

  /* Update the global statistics */

  st_variogram_stats(db, vario);

  /* Loop on the directions to evaluate */

  for (idir = 0; idir < vario->getDirectionNumber(); idir++)
    st_variogen_line(db, vario, idir, norder);

  /* Set the error return code */

  error = 0;

  return (error);
}

/****************************************************************************/
/*!
 **  Print the parameters for calculating a variogram in a direction
 **
 *****************************************************************************/
static void st_vario_params_print(int ndim,
                                  const VectorDouble &codir,
                                  double tolang,
                                  double bench,
                                  double cylrad)
{
  VectorDouble angles(ndim);

  message("Direction coefficients      = (");
  for (int idim = 0; idim < ndim; idim++)
    tab_printg(NULL, 1, EJustify::LEFT, codir[idim]);
  message(")\n");
  if (ndim > 1)
  {
    (void) ut_angles_from_codir(ndim, 1, codir, angles);
    message("Direction angles (degrees)  = (");
    for (int idim = 0; idim < ndim; idim++)
      tab_printg(NULL, 1, EJustify::LEFT, angles[idim]);
    message(")\n");
  }
  if (!FFFF(tolang))
    message("Tolerance on direction      = %lf (deg)\n", tolang);
  if (!FFFF(bench) && bench > 0.)
    message("Slice bench                 = %lf\n", bench);
  if (!FFFF(cylrad) && cylrad > 0.)
    message("Slice radius                = %lf\n", cylrad);
}

/****************************************************************************/
/*!
 **  Print the experimental variograms in one direction
 **
 ** \param[in]  vario     Vario structure
 ** \param[in]  idir      Rank of the direction
 ** \param[in]  verbose   0 for brief output; 1 for a long output
 **
 *****************************************************************************/
void vardir_print(Vario *vario, int idir, int verbose)
{
  if (vario == nullptr) return;
  if (idir < 0 || idir >= vario->getDirectionNumber()) return;
  message(vario->getDirParam(idir).toString(verbose).c_str());
  return;
}

/****************************************************************************/
/*!
 **  Print the experimental variograms
 **
 ** \param[in]  vario     Vario structure
 ** \param[in]  verbose 0 for brief output; 1 for a long output
 **
 *****************************************************************************/
void variogram_print(const Vario *vario, int verbose)
{
  if (vario != nullptr) messageFlush(vario->toString(verbose));
}

/****************************************************************************/
/*!
 **  Estimate the coefficients of the global drift
 **
 ** \return Error return code
 **
 ** \param[in]  db      Db descriptor
 ** \param[in]  verbose Verbose flag
 **
 *****************************************************************************/
static int st_estimate_drift_coefficients(Db *db, int verbose)

{
  int nbfl, nech, i, il, jl, iech, iiech, error;
  double *b;

  /* Initializations */

  if (MODEL == nullptr) return (0);
  error = 1;
  nbfl = MODEL->getDriftNumber();
  nech = db->getActiveAndDefinedNumber(0);
  b = nullptr;

  /* Core allocation */

  b = (double*) mem_alloc(sizeof(double) * nbfl, 1);
  for (i = 0; i < nbfl; i++)
    b[i] = 0.;

  /* Calculate: t(X) %*% X */

  for (iech = iiech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActiveAndDefined(iech, 0)) continue;
    model_calcul_drift(MODEL, ECalcMember::LHS, db, iech, DRFLOC);

    for (il = 0; il < nbfl; il++)
    {
      if (FFFF(DRFLOC[il]))
      {
        messerr(
            "Drift cannot be calculated: term (%d) is undefined at sample (%d)",
            il + 1, iech + 1);
        goto label_end;
      }
      X_DRFTAB(il,iiech) = DRFLOC[il];
      b[il] += DRFLOC[il] * db->getVariable(iech, 0);
      for (jl = 0; jl < nbfl; jl++)
        X_MATDRF(il,jl) += DRFLOC[il] * DRFLOC[jl];
    }
    iiech++;
  }

  /* Calculate: A = (t(X) %*% X)-1 */

  if (matrix_invert(MATDRF, nbfl, 0)) goto label_end;

  /* Calculate: A %*% t(X) %*% Y */

  matrix_product(nbfl, nbfl, 1, MATDRF, b, BETA);

  /* Optional printout */

  if (verbose)
  {
    message("Drift removal initial step\n");
    print_matrix("Drift Coefficients Matrix", 0, 1, nbfl, nbfl, NULL, MATDRF);
  }

  /* Pre-process the vector X %*% A */

  matrix_product(nech, nbfl, nbfl, DRFTAB, MATDRF, DRFXA);

  /* Set the error return code */

  error = 0;

  label_end: b = (double*) mem_free((char* ) b);
  return (error);
}

/****************************************************************************/
/*!
 **  Evaluate the experimental variogram on irregular data
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  vario        Vario structure
 ** \param[in]  model        Model structure (triggers the KU option)
 ** \param[in]  flag_sample  calculate the variogram per sample
 ** \param[in]  verr_mode    Mode of variogram correction (1, 2 or 3)
 ** \param[in]  verbose      Verbose flag
 **
 ** \note The number of iterations in KU is given can be updated using the
 ** \note keypair technique with name "KU_Niter".
 **
 *****************************************************************************/
static int st_variogram_general(Db *db,
                                Vario *vario,
                                Model *model,
                                int flag_sample,
                                int verr_mode,
                                int verbose)
{
  int idir, error, flag_verr, flag_ku, nbfl;
  Vario_Order *vorder;
  VectorInt rindex;

  /* Initializations */

  error = 1;
  flag_verr = flag_ku = nbfl = 0;
  if (db == nullptr) return (1);
  if (vario == nullptr) return (1);
  vorder = (Vario_Order*) NULL;
  st_manage_drift_removal(0, NULL, NULL);

  /* Preliminary checks */

  if (db->getNDim() != vario->getDimensionNumber() || db->getVariableNumber()
      != vario->getVariableNumber())
  {
    messerr("Inconsistent parameters:");
    messerr("Data Base: NDIM=%d NVAR=%d", db->getNDim(),
            db->getVariableNumber());
    messerr("Variogram: NDIM=%d NVAR=%d", vario->getDimensionNumber(),
            vario->getVariableNumber());
    goto label_end;
  }
  if (st_get_generalized_variogram_order(vario) > 0)
  {
    messerr("This calculation does not allow generalized variogram definition");
    goto label_end;
  }

  /* Particular case of Transitive Covariogram */
  /* It is only coded in the by_sample case and uses the regression technique */

  if (vario->getCalculType() == ECalcVario::COVARIOGRAM) flag_sample = 1;

  /* Auxiliary check for Variance Measurement Error */

  if (db->getVarianceErrorNumber() > 0 && verr_mode > 0)
  {
    vorder = vario_order_manage(1, 1, 0, vorder);
    flag_verr = 1;
  }

  /* Auxiliary check for Drift removal */

  if (model != nullptr && (model->getDriftNumber() > 1
      || model->getDriftType(0) != EDrift::UC))
  {
    if (vorder == (Vario_Order*) NULL)
      vorder = vario_order_manage(1, 1, 0, vorder);
    flag_ku = 1;
    st_manage_drift_removal(1, db, model);
  }

  /* Complementary checks */

  if (flag_verr && flag_ku)
  {
    messerr("These two options are incompatible");
    messerr("- Correction for the Variance of Error Measurements");
    messerr("- Correction for bias when removing the Drift");
    goto label_end;
  }
  if (flag_verr || flag_ku)
  {
    if (flag_sample)
    {
      messerr("The special Variogram option is incompatible with flag.sample");
      goto label_end;
    }
    if (!db->isVariableNumberComparedTo(1)) goto label_end;
  }

  /* Evaluate the drift coefficients */

  if (flag_ku)
  {
    if (st_estimate_drift_coefficients(db, verbose)) goto label_end;
  }

  /* Update the global statistics */

  st_variogram_stats(db, vario);

  /* Loop on the directions to evaluate */

  rindex = db->getSortArray();
  for (idir = 0; idir < vario->getDirectionNumber(); idir++)
  {
    if (!flag_sample)
    {
      if (st_variogram_calcul1(db, vario, idir, rindex.data(), vorder))
        goto label_end;
    }
    else
    {
      if (st_variogram_calcul2(db, vario, idir, rindex.data())) goto label_end;
    }

    if (vorder != (Vario_Order*) NULL)
      st_variogram_calcul_internal(db, vario, idir, vorder);
  }

  /* Posterior calculations when presence of Variance of Measurement errors */

  if (flag_verr)
  {
    for (idir = 0; idir < vario->getDirectionNumber(); idir++)
    {
      if (st_update_variogram_verr(db, vario, idir, vorder, verr_mode))
        goto label_end;
    }
  }

  /* Posterior update when filtering the bias attached to drift removal */

  if (flag_ku)
  {
    if (st_update_variogram_ku(db, vario, vorder, verbose)) goto label_end;
  }

  /* Set the error return code */

  error = 0;

  label_end: vorder = vario_order_manage(-1, 1, 0, vorder);
  st_manage_drift_removal(-1, db, model);
  return (error);
}

/****************************************************************************/
/*!
 **  Evaluate the experimental covariance for directional variables
 **
 ** \return  Error return code
 **
 ** \param[in]  db     Db descriptor
 ** \param[in]  vario  Vario structure
 ** \param[in]  ncomp  Number of components
 **
 *****************************************************************************/
int variovect_compute(Db *db, Vario *vario, int ncomp)
{
  int idir, error;
  VectorInt rindex;

  /* Initializations */

  error = 0;
  if (db == nullptr) return (1);
  if (vario == nullptr) return (1);
  st_manage_drift_removal(0, NULL, NULL);

  /* Preliminary checks */

  if (db->getNDim() != vario->getDimensionNumber() || db->getVariableNumber()
      != vario->getVariableNumber() * ncomp)
  {
    messerr("Inconsistent parameters:");
    messerr("Data Base: NDIM=%d NVAR=%d", db->getNDim(),
            db->getVariableNumber());
    messerr("Variogram: NDIM=%d NVAR=%d", vario->getDimensionNumber(),
            vario->getVariableNumber());
    messerr("Number of components = %d", ncomp);
    return (1);
  }

  /* Update the global statistics */

  st_variovect_stats(db, vario, ncomp);

  /* Loop on the directions to evaluate */

  rindex = db->getSortArray();
  for (idir = 0; idir < vario->getDirectionNumber(); idir++)
  {
    error = st_variovect_calcul(db, vario, idir, ncomp, rindex.data());
    if (error) break;
  }
  return (error);
}

/****************************************************************************/
/*!
 **  Scale the variogram maps
 **
 ** \param[in]  dbmap     Db of Grid type containing the Variogram Maps
 ** \param[in]  nv2       nvar ( nvar + 1) /2
 **
 *****************************************************************************/
static void st_vmap_scale(Db *dbmap, int nv2)
{
  for (int iech = 0; iech < dbmap->getSampleNumber(); iech++)
  {
    for (int ijvar = 0; ijvar < nv2; ijvar++)
    {
      double value = dbmap->getArray(iech, IPTW + ijvar);
      if (value <= 0.)
        dbmap->setArray(iech, IPTV + ijvar, TEST);
      else
        dbmap->updArray(iech, IPTV + ijvar, 3, value);
    }
  }
}

/****************************************************************************/
/*!
 **  Get the absolute index of a grid node, shifted within a neighborhood
 **
 ** \return  Returned absolute index (<0 is not relevant)
 **
 ** \param[in]  dbmap     Db of Grid type containing the Variogram Maps
 ** \param[in]  indg0     Index decomposition for the central cell
 ** \param[in]  neigh     Array describing the neighborhooed
 **                       (Dimension: ndim * nbmax)
 ** \param[in]  rank      Rank of the cell within the neghborhood
 **
 ** \param[out] indg1     Working array for grid indices
 **
 *****************************************************************************/
static int st_find_neigh_cell(Db *dbmap,
                              int *indg0,
                              int *neigh,
                              int rank,
                              int *indg1)
{
  int ndim;

  // Initializations

  ndim = dbmap->getNDim();

  // Get the indices of the neighboring cell

  for (int idim = 0; idim < ndim; idim++)
    indg1[idim] = indg0[idim] + neigh[rank * ndim + idim];

  return (db_index_grid_to_sample(dbmap, indg1));
}

/****************************************************************************/
/*!
 **  Calculate the variogram map when data are isolated points
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db containing the data
 ** \param[in]  dbmap        Db of Grid type containing the Variogram Maps
 ** \param[in]  calcul_type  Type of calculation (ECalcVario)
 ** \param[in]  radius       Dilation radius (used to smooth the resulting maps)
 ** \param[in]  namconv      Naming convention
 **
 *****************************************************************************/
static int st_vmap_general(Db *db,
                           Db *dbmap,
                           const ECalcVario &calcul_type,
                           int radius,
                           const NamingConvention& namconv)
{
  int error, nvar, nv2, i, idim, flag_out, nbmax;
  int *indg0, *indg1, *ind1, iech0, iech1, iech2, jech1, jech2, nech, ndim;
  double *delta, *mid, *coor, x0;
  VectorInt neigh;

  /* Preliminary checks */

  error = 0;
  indg0 = indg1 = ind1 = nullptr;
  delta = coor = mid = nullptr;
  if (db == nullptr) return (1);
  if (dbmap == nullptr) return (1);

  if (!is_grid(dbmap))
  {
    messerr("This feature requires a Grid Data Base");
    messerr("to store the Variogram Maps");
    return (1);
  }
  if (db->getNDim() != 2 && db->getNDim() != 3)
  {
    messerr("The Variogram Map can only be calculated on a grid data set");
    messerr("with dimension equal to 2 or 3");
    return (1);
  }
  if (dbmap->getNDim() > db->getNDim())
  {
    messerr("The space dimension of the VMAP (%d)", dbmap->getNDim());
    messerr(
        "must not be larger than the space dimension of the input Grid (%d)",
        db->getNDim());
    return (1);
  }

  /* Initializations */

  ndim = dbmap->getNDim();
  nvar = db->getVariableNumber();
  nech = db->getSampleNumber();
  nv2 = nvar * (nvar + 1) / 2;

  /* Core allocation */

  indg0 = db_indg_alloc(dbmap);
  if (indg0 == nullptr) goto label_end;
  indg1 = db_indg_alloc(dbmap);
  if (indg1 == nullptr) goto label_end;
  ind1 = (int*) mem_alloc(sizeof(int) * nech, 0);
  if (ind1 == nullptr) goto label_end;
  delta = db_sample_alloc(db, ELoc::X);
  if (delta == nullptr) goto label_end;
  mid = db_sample_alloc(db, ELoc::X);
  if (mid == nullptr) goto label_end;
  coor = db_vector_alloc(db);
  if (coor == nullptr) goto label_end;

  /* Create the variables in the Variogram Map file */

  IPTV = dbmap->addFieldsByConstant(nv2, 0.);
  if (IPTV < 0) goto label_end;
  IPTW = dbmap->addFieldsByConstant(nv2, 0.);
  if (IPTW < 0) goto label_end;

  /* Calculate a neighborhood (if radius > 0) */

  neigh = gridcell_neigh(ndim, 1, radius, 0, 0, &nbmax);

  /* Calculate the VMAP half-extension */

  for (idim = 0; idim < ndim; idim++)
    mid[idim] = dbmap->getNX(idim) * dbmap->getDX(idim) / 2;

  /* Sorting the samples according to their first coordinate */

  if (db_coorvec_get(db, 0, coor)) goto label_end;
  for (i = 0; i < nech; i++)
    ind1[i] = i;
  ut_sort_double(1, nech, ind1, coor);

  /* Loop on the first data */

  for (jech1 = 0; jech1 < nech; jech1++)
  {
    iech1 = ind1[jech1];
    if (!db->isActive(iech1)) continue;
    x0 = db->getCoordinate(iech1, 0);

    /* Loop on the second data */

    for (jech2 = jech1; jech2 < nech; jech2++)
    {
      iech2 = ind1[jech2];
      if (!db->isActive(iech2)) continue;
      delta[0] = db->getCoordinate(iech2, 0) - x0;
      if (delta[0] > mid[0]) break;

      for (idim = 1, flag_out = 0; idim < ndim && flag_out == 0; idim++)
      {
        delta[idim] = db->getDistance1D(iech2, iech1, idim);
        if (delta[idim] > mid[idim]) flag_out = 1;
      }
      if (flag_out) continue;

      /* Evaluate the variogram map */

      DBMAP = dbmap;

      // Apply to the target cell
      if (point_to_grid(dbmap, delta, 0, indg0)) continue;
      for (int in = 0; in < nbmax; in++)
      {
        iech0 = st_find_neigh_cell(dbmap, indg0, neigh.data(), in, indg1);
        if (iech0 < 0) continue;
        st_variogram_evaluate(db, calcul_type, nvar, iech1, iech2, iech0, TEST,
                              0, st_vmap_set);
      }

      // Avoid symmetry if point is compared to itself
      if (iech1 == iech2) continue;

      // Apply to the opposite target cell
      for (idim = 0; idim < ndim; idim++)
        delta[idim] = -delta[idim];
      if (point_to_grid(dbmap, delta, 0, indg0)) continue;
      for (int in = 0; in < nbmax; in++)
      {
        iech0 = st_find_neigh_cell(dbmap, indg0, neigh.data(), in, indg1);
        if (iech0 < 0) continue;
        st_variogram_evaluate(db, calcul_type, nvar, iech1, iech2, iech0, TEST,
                              0, st_vmap_set);
      }
    }
  }

  /* Normalization */

  st_vmap_scale(dbmap, nv2);

  /* Set the error return code */

  namconv.setNamesAndLocators(db, ELoc::Z, -1, dbmap, IPTW, "Nb", 1, false);
  namconv.setNamesAndLocators(db, ELoc::Z, -1, dbmap, IPTV, "Var");
  error = 0;

  label_end: indg0 = db_indg_free(indg0);
  indg1 = db_indg_free(indg1);
  ind1 = (int*) mem_free((char* ) ind1);
  delta = db_sample_free(delta);
  mid = db_sample_free(mid);
  coor = db_vector_free(coor);
  return (error);
}

/****************************************************************************/
/*!
 **  Calculate the variogram map when data are defined on a Grid
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid       Db of Grid type containing the data
 ** \param[in]  dbmap        Db of Grid type containing the Variogram Maps
 ** \param[in]  calcul_type  Type of calculation (ECalcVario)
 ** \param[in]  namconv      Naming Convention
 **
 *****************************************************************************/
static int st_vmap_grid(Db *dbgrid,
                        Db *dbmap,
                        const ECalcVario &calcul_type,
                        const NamingConvention& namconv)
{
  int error, nvar, nv2, idim, delta;
  int *ind1, *ind2, *ind0, iech0, iech1, iech2, flag_out, ndim;

  /* Preliminary checks */

  error = 0;
  ind0 = ind1 = ind2 = nullptr;
  if (dbgrid == nullptr) return (1);
  if (dbmap == nullptr) return (1);

  if (!is_grid(dbgrid))
  {
    messerr("This Variogram Map is defined for Grid Data Base only");
    return (1);
  }
  if (!is_grid(dbmap))
  {
    messerr("This feature requires a Grid Data Base");
    messerr("to store the Variogram Maps");
    return (1);
  }
  if (dbgrid->getNDim() != 2 && dbgrid->getNDim() != 3)
  {
    messerr("The Variogram Map can only be calculated on a grid data set");
    messerr("with dimension equal to 2 or 3");
    return (1);
  }
  if (dbmap->getNDim() > dbgrid->getNDim())
  {
    messerr("The space dimension of the VMAP (%d)", dbmap->getNDim());
    messerr(
        "must not be larger than the space dimension of the input Grid (%d)",
        dbgrid->getNDim());
    return (1);
  }
  for (idim = 0; idim < dbmap->getNDim(); idim++)
  {
    if (ABS(dbmap->getDX(idim) - dbgrid->getDX(idim)) > 1.e-03)
    {
      messerr("The grid mesh in the direction %d (dx=%lf)", idim + 1,
              dbgrid->getDX(idim));
      messerr("must match the mesh of the Variogram Map grid (dx=%lf)",
              dbmap->getDX(idim));
      return (1);
    }
  }

  /* Initializations */

  ndim = dbmap->getNDim();
  nvar = dbgrid->getVariableNumber();
  nv2 = nvar * (nvar + 1) / 2;

  /* Core allocation */

  ind0 = db_indg_alloc(dbmap);
  if (ind0 == nullptr) goto label_end;
  ind1 = db_indg_alloc(dbgrid);
  if (ind1 == nullptr) goto label_end;
  ind2 = db_indg_alloc(dbgrid);
  if (ind2 == nullptr) goto label_end;

  /* Create the variables in the Variogram Map file */

  IPTV = dbmap->addFieldsByConstant(nv2, 0.);
  if (IPTV < 0) goto label_end;
  IPTW = dbmap->addFieldsByConstant(nv2, 0.);
  if (IPTW < 0) goto label_end;

  /* Loop on the first data */

  for (iech1 = 0; iech1 < dbgrid->getSampleNumber(); iech1++)
  {
    if (!dbgrid->isActive(iech1)) continue;
    db_index_sample_to_grid(dbgrid, iech1, ind1);

    /* Loop on the second data */

    for (iech2 = 0; iech2 < dbgrid->getSampleNumber(); iech2++)
    {
      if (!dbgrid->isActive(iech2)) continue;
      db_index_sample_to_grid(dbgrid, iech2, ind2);

      for (idim = flag_out = 0; idim < ndim && flag_out == 0; idim++)
      {
        delta = ind1[idim] - ind2[idim];
        int moitie = (dbmap->getNX(idim) - 1) / 2;
        if (delta < -moitie || delta > moitie) flag_out = 1;
        ind0[idim] = delta + moitie;
      }
      if (flag_out) continue;
      iech0 = db_index_grid_to_sample(dbmap, ind0);

      /* Evaluate the variogram map */

      DBMAP = dbmap;
      st_variogram_evaluate(dbgrid, calcul_type, nvar, iech1, iech2, iech0,
      TEST,
                            0, st_vmap_set);
    }
  }

  /* Normalization */

  st_vmap_scale(dbmap, nv2);

  /* Set the error return code */

  namconv.setNamesAndLocators(dbgrid, ELoc::Z, -1, dbmap, IPTW, "Nb", 1, false);
  namconv.setNamesAndLocators(dbgrid, ELoc::Z, -1, dbmap, IPTV, "Var");
  error = 0;

  label_end: ind0 = db_indg_free(ind0);
  ind1 = db_indg_free(ind1);
  ind2 = db_indg_free(ind2);
  return (error);
}

/****************************************************************************/
/*!
 **  Calculate the variogram extension for a pair of variables
 **
 ** \param[in]  vario     Vario structure
 ** \param[in]  ivar      Rank of the first variable
 ** \param[in]  jvar      Rank of the second variable
 ** \param[in]  idir0     Rank of the direction (-1 for all)
 ** \param[in]  flag_norm 1 if the variogram must be normalized by variance
 ** \param[in]  flag_vars 1 if the global statistics must be taken into account
 ** \param[in]  distmin   Minimum along the distance axis
 ** \param[in]  distmax   Maximum along the distance axis
 ** \param[in]  varmin    Minimum along the variogram (or covariance) axis
 ** \param[in]  varmax    Maximum along the variogram (or covariance) axis
 **
 ** \param[out]  flag_hneg 1 if the distance scale can be negative
 ** \param[out]  flag_gneg 1 if the variogram scale can be negative
 ** \param[out]  c0        Value of the variogram at the origin
 ** \param[out]  hmin      Minimum distance
 ** \param[out]  hmax      Maximum distance
 ** \param[out]  gmin      Minimum variogram value
 ** \param[out]  gmax      Maximum variogram value
 **
 *****************************************************************************/
void variogram_extension(const Vario *vario,
                         int ivar,
                         int jvar,
                         int idir0,
                         int flag_norm,
                         int flag_vars,
                         double distmin,
                         double distmax,
                         double varmin,
                         double varmax,
                         int *flag_hneg,
                         int *flag_gneg,
                         double *c0,
                         double *hmin,
                         double *hmax,
                         double *gmin,
                         double *gmax)
{
  double hh, gg;
  int i, j, idir, jdir, ndir;
  double tol = 0.1;

  /* Initializations */

  (*hmin) = 0.;
  (*gmin) = 0.;
  (*hmax) = -1.e30;
  (*gmax) = -1.e30;
  if (vario->getFlagAsym())
  {
    *c0 = vario->getGgByIndex(0, vario->getDirAddress(0, ivar, jvar, 0, false, 0));
  }
  else
  {
    *c0 = vario->getVar(ivar, jvar);
  }
  if (st_get_generalized_variogram_order(vario) > 0) (*c0) = TEST;
  if (FFFF(*c0) && flag_norm)
  {
    messerr("The Normalization option is discarded for this variogram");
    messerr("probably as it corresponds to a generalized variogram");
    flag_norm = 0;
  }
  (*flag_hneg) = (ivar != jvar && vario->getFlagAsym());
  (*flag_gneg) = (ivar != jvar || vario->getFlagAsym());

  /* Loop on the directions */

  ndir = vario->getDirectionNumber();
  if (idir0 >= 0) ndir = 1;
  for (jdir = 0; jdir < ndir; jdir++)
  {
    idir = (idir0 >= 0) ? idir0 :
                          jdir;
    for (i = 0; i < vario->getLagTotalNumber(idir); i++)
    {
      j = vario->getDirAddress(idir, ivar, jvar, i, true, 0);
      if (vario->getSwByIndex(idir, j) <= 0) continue;
      hh = vario->getHhByIndex(idir, j);
      gg = vario->getGgByIndex(idir, j);
      if (FFFF(hh) || FFFF(gg)) continue;
      if (flag_norm) gg /= (*c0);
      if (!FFFF(distmin) && hh < distmin) continue;
      if (!FFFF(distmax) && hh > distmax) continue;
      if (hh < (*hmin)) (*hmin) = hh;
      if (hh > (*hmax)) (*hmax) = hh;
      if (gg < (*gmin)) (*gmin) = gg;
      if (gg > (*gmax)) (*gmax) = gg;
    }
  }

  if (flag_norm) (*c0) = 1.;
  if (!FFFF(*c0) && flag_vars)
  {
    if ((*c0) < (*gmin)) (*gmin) = (*c0);
    if ((*c0) > (*gmax)) (*gmax) = (*c0);
  }

  /* Expand the vertical graphic scales */

  (*gmax) *= (1. + tol);
  if ((*gmin) < 0) (*gmin) *= (1. + tol);

  /* Correction due to variogram type */

  if (*flag_hneg)
  {
    (*hmax) = MAX(ABS(*hmin), ABS(*hmax));
    (*hmin) = -(*hmax);
  }

  if (*flag_gneg)
  {
    (*gmax) = MAX(ABS(*gmin), ABS(*gmax));
    if (ivar != jvar) (*gmin) = -(*gmax);
  }
  else
  {
    (*gmin) = 0.;
  }

  /* Truncation to limits provided by the user */

  if (!FFFF(distmax)) (*hmax) = distmax;
  if (!FFFF(distmin)) (*hmin) = distmin;
  if (!FFFF(varmax)) (*gmax) = varmax;
  if (!FFFF(varmin)) (*gmin) = varmin;

  return;
}

/****************************************************************************/
/*!
 **  Evaluate the experimental variogram on grid
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  vario        Vario structure
 **
 *****************************************************************************/
static int st_variogrid_calcul(Db *db, Vario *vario)
{
  int idir, error, iadd_new, iatt_old, iech;
  double maille;

  /* Initializations */

  error = 1;
  iadd_new = iatt_old = -1;
  if (db == nullptr) return (1);
  if (vario == nullptr) return (1);
  st_manage_drift_removal(0, NULL, NULL);

  /* Preliminary checks */

  if (db->getNDim() != vario->getDimensionNumber() || db->getVariableNumber()
      != vario->getVariableNumber())
  {
    messerr("Inconsistent parameters:");
    messerr("Data Base: NDIM=%d NVAR=%d", db->getNDim(),
            db->getVariableNumber());
    messerr("Variogram: NDIM=%d NVAR=%d", vario->getDimensionNumber(),
            vario->getVariableNumber());
    goto label_end;
  }
  if (!is_grid(db))
  {
    messerr("This calculation facility is dedicated to grid architecture");
    goto label_end;
  }
  if (st_get_generalized_variogram_order(vario) > 0)
  {
    messerr("This calculation does not allow generalized variogram definition");
    goto label_end;
  }

  /* In the case of Covariogram, add the weight set to the scale */

  if (vario->getCalculType() == ECalcVario::COVARIOGRAM)
  {
    iatt_old = db_attribute_identify(db, ELoc::W, 0);
    iadd_new = db->addFieldsByConstant(1, 0.);
    if (iadd_new < 0) goto label_end;
    db->setLocatorByAttribute(iadd_new, ELoc::W);
    maille = db_grid_maille(db);
    for (iech = 0; iech < db->getSampleNumber(); iech++)
      db->setWeight(iech, maille);
  }

  /* Update the global statistics */

  st_variogram_stats(db, vario);

  /* Loop on the directions to evaluate */

  for (idir = 0; idir < vario->getDirectionNumber(); idir++)
  {
    error = st_variogram_grid(db, vario, idir);
    if (error) break;
  }

  /* Delete the additional weight variable (optional) */

  if (vario->getCalculType() == ECalcVario::COVARIOGRAM)
  {
    if (iadd_new > 0) db->deleteFieldByAttribute(iadd_new);
    if (iatt_old > 0) db->setLocatorByAttribute(iatt_old, ELoc::W);
  }

  /* Set the error return code */

  error = 0;

  label_end: return (error);
}

/****************************************************************************/
/*!
 **  Initialize a new calculation direction
 **
 ** \return  Error return code
 **
 ** \param[in]  varioparam   VarioParam structure
 ** \param[in]  npas         number of lags
 ** \param[in]  opt_code     code selection option
 ** \li                      0 : no use of the code selection
 ** \li                      1 : codes must be close enough
 ** \li                      2 : codes must be different
 ** \param[in]  idate        Rank of the Date interval
 ** \param[in]  dpas         lag value
 ** \param[in]  toldis       tolerance on distance (proportion of the lag)
 ** \param[in]  tolang       angular tolerance (in degrees)
 ** \param[in]  bench        Slicing bench
 ** \param[in]  cylrad       Slicing radius
 ** \param[in]  tolcode      Tolerance on the code
 ** \param[in]  breaks       array for irregular lags
 ** \param[in]  codir        calculation direction (Dimension = ndim)
 ** \param[in]  grincr       direction increment only used for grid
 **                          (Dimension = ndim)
 **
 *****************************************************************************/
int variogram_direction_add(VarioParam *varioparam,
                            int npas,
                            int opt_code,
                            int idate,
                            double dpas,
                            double toldis,
                            double tolang,
                            double bench,
                            double cylrad,
                            double tolcode,
                            const VectorDouble &breaks,
                            const VectorDouble &codir,
                            const VectorInt &grincr)
{
  if (varioparam == (VarioParam*) NULL) return (1);
  DirParam dirparam = DirParam(varioparam->getDimensionNumber(), npas, dpas,
                               toldis, tolang, opt_code, idate, bench, cylrad,
                               tolcode, breaks, codir, grincr);
  varioparam->addDirs(dirparam);
  return (0);
}

/****************************************************************************/
/*!
 **  Free the Vario structure
 **
 ** \return  Pointer to the freed Vario structure
 **
 ** \param[in]  vario Vario structure to be freed
 **
 *****************************************************************************/
Vario* variogram_delete(Vario *vario)

{
  if (vario == nullptr) return (vario);
  delete vario;
  vario = nullptr;
  return (vario);
}

/****************************************************************************/
/*!
 **  Replace zero values by TEST values
 **
 ** \param[in]  db    Discretization Grid descriptor
 ** \param[in]  iptr  Pointer of the attribute to be modified
 **
 *****************************************************************************/
static void st_final_discretization_grid(Db *db, int iptr)
{
  int iech, nech;
  double value;

  nech = db->getSampleNumber();
  for (iech = 0; iech < nech; iech++)
  {
    value = db->getArray(iech, iptr);
    if (value != 0.) continue;
    db->setArray(iech, iptr, TEST);
  }
}

/****************************************************************************/
/*!
 **  Add one pick to the discretization grid
 **
 ** \return  Index of the grid cell (or -1)
 **
 ** \param[in]  db    Discretization Grid descriptor
 ** \param[in]  x     Coordinate along the first axis
 ** \param[in]  y     Coordinate along the first axis
 **
 *****************************************************************************/
static int st_update_discretization_grid(Db *db, double x, double y)
{
  int iech, ix, iy, indg[2];

  ix = (int) floor((x - db->getX0(0)) / db->getDX(0) + 0.5);
  iy = (int) floor((y - db->getX0(1)) / db->getDX(1) + 0.5);
  if (ix < 0 || ix >= db->getNX(0)) return (-1);
  if (iy < 0 || iy >= db->getNX(1)) return (-1);
  indg[0] = ix;
  indg[1] = iy;
  iech = db_index_grid_to_sample(db, indg);
  return (iech);
}

/****************************************************************************/
/*!
 **  Evaluate the correlation (according to argument flag_same):
 **  - the standard correlation (flag_same = 1)
 **     Correl(Z1(x) , Z2(x))
 **  - the shifted correlation calculated as follows:
 **     Correl(Z1(x) , Z2(x+h))
 **
 ** \return  Error return code
 **
 ** \param[in]  db1          Db descriptor (first variable)
 ** \param[in]  db2          Db descriptor (second variable for flag.same=T)
 ** \param[in]  dbgrid       Discretization Grid descriptor
 ** \param[in]  flag_same    1 if the two samples must coincide
 ** \param[in]  icol1        Rank of the first column
 ** \param[in]  icol2        Rank of the second column
 ** \param[in]  flag_verbose 1 for a verbose output
 ** \param[in]  dmin         Minimum distance
 ** \param[in]  dmax         Maximum distance
 ** \param[in]  tolang       Angular tolerance (in degrees)
 ** \param[in]  bench        Slicing bench
 ** \param[in]  cylrad       Slicing radius
 ** \param[in]  codir        Calculation direction
 ** \param[in]  opt_code     code selection option
 ** \li                       0 : no use of the code selection
 ** \li                       1 : codes must be close enough
 ** \li                       2 : codes must be different
 ** \param[in]  tolcode      Tolerance on the code
 **
 ** \param[out] nindice      Number of pairs
 ** \param[out] indices      Array of the indices of pairs of samples
 **                          (Dimension: 2 * nindice)
 ** \param[out] correl       Correlation coefficient
 **
 ** \remarks The two input Db must match exactly (same number of samples with
 ** \remarks same set of coordinates and same optional selection)
 **
 ** \remarks The returned array 'indices' contain the set of indices of the
 ** \remarks pairs of samples. Its contents is i1,j1,i2,j2,...
 ** \remarks The indices are numbered starting from 1
 ** \remarks The number of pairs is 'nindice'.
 ** \remarks This array must be freed by the calling function
 **
 *****************************************************************************/
int correlation_f(Db *db1,
                  Db *db2,
                  Db *dbgrid,
                  int flag_same,
                  int icol1,
                  int icol2,
                  int flag_verbose,
                  double dmin,
                  double dmax,
                  double tolang,
                  double bench,
                  double cylrad,
                  VectorDouble &codir,
                  int opt_code,
                  int tolcode,
                  int *nindice,
                  int **indices,
                  double *correl)
{
  int *ind, iech, jech, nech, iptr, igrid, nalloc, npair, error;
  double ps, psmin, dist, val1, val2, m1, m2, v1, v2, v12, nb, rho;

  /* Initializations */

  if (db1 == nullptr) return (1);
  if (db2 == nullptr) return (1);
  if (dbgrid == nullptr) return (1);
  nech = db1->getSampleNumber();
  error = 1;
  iptr = nalloc = 0;
  *correl = 0.;
  *nindice = 0;
  *indices = ind = nullptr;
  iptr = dbgrid->addFieldsByConstant(1, 0.);
  if (iptr < 0) return (1);
  m1 = m2 = v1 = v2 = v12 = nb = 0.;

  /* Initial core allocation */

  npair = 0;
  nalloc = nech;
  ind = (int*) mem_alloc(sizeof(int) * 2 * nalloc, 0);
  if (ind == nullptr) goto label_end;

  /* Dispatch */

  if (flag_same)
  {

    /* Regular correlation */

    for (iech = 0; iech < nech; iech++)
    {
      if (!db1->isActive(iech)) continue;
      val1 = db1->getArray(iech, icol1);
      if (FFFF(val1)) continue;
      val2 = db2->getArray(iech, icol2);
      if (FFFF(val2)) continue;

      /* Global statistics */

      nb += 1.;
      m1 += val1;
      m2 += val2;
      v1 += val1 * val1;
      v2 += val2 * val2;
      v12 += val1 * val2;

      /* Point update */

      if (npair >= nalloc)
      {
        nalloc += nech;
        ind = (int*) mem_realloc((char* ) ind, 2 * nalloc * sizeof(int), 0);
        if (ind == nullptr) goto label_end;
      }
      ind[2 * npair] = iech + 1;
      ind[2 * npair + 1] = iech + 1;
      npair++;

      /* Grid update */

      igrid = st_update_discretization_grid(dbgrid, val1, val2);
      if (igrid < 0) continue;
      dbgrid->updArray(igrid, iptr, 0, 1.);
    }
  }
  else
  {

    /* Shifted correlation */

    vario_fix_codir(db1->getNDim(), codir);
    psmin = _variogram_convert_angular_tolerance(tolang);

    for (iech = 0; iech < nech - 1; iech++)
    {
      if (!db1->isActive(iech)) continue;
      val1 = db1->getArray(iech, icol1);
      if (FFFF(val1)) continue;

      for (jech = iech + 1; jech < nech; jech++)
      {
        if (!db1->isActive(jech)) continue;
        val2 = db1->getArray(jech, icol2);
        if (FFFF(val2)) continue;

        /* Check if the pair must be kept (Code criterion) */

        if (code_comparable(db1, db1, iech, jech, opt_code, tolcode)) continue;

        /* Check if the pair must be kept */

        dist = distance_intra(db1, iech, jech, NULL);
        if (variogram_reject_pair(db1, iech, jech, dist, psmin, bench, cylrad,
                                  codir, &ps)) continue;

        /* Check the distance */

        if (dist < dmin || dist > dmax) continue;

        /* Global statistics */

        nb += 1.;
        m1 += val1;
        m2 += val2;
        v1 += val1 * val1;
        v2 += val2 * val2;
        v12 += val1 * val2;

        /* Point update */

        if (npair >= nalloc)
        {
          nalloc += nech;
          ind = (int*) mem_realloc((char* ) ind, 2 * nalloc * sizeof(int), 0);
          if (ind == nullptr) goto label_end;
        }
        ind[2 * npair] = iech + 1;
        ind[2 * npair + 1] = jech + 1;
        npair++;

        /* Grid update */

        igrid = st_update_discretization_grid(dbgrid, val1, val2);
        if (igrid < 0) continue;
        dbgrid->updArray(igrid, iptr, 0, 1.);
      }
    }
  }

  /* Calculate the correlation coefficient */

  *correl = TEST;
  if (nb > 0)
  {
    m1 /= nb;
    m2 /= nb;
    v1 = v1 / nb - m1 * m1;
    v2 = v2 / nb - m2 * m2;
    v12 = v12 / nb - m1 * m2;
    rho = (v1 * v2 > 0) ? v12 / sqrt(v1 * v2) :
                          1.;
    if (flag_verbose) message("Correlation coefficient = %lf\n", rho);
    *correl = rho;
  }
  else
  {
    messerr("No sample found where all variables are defined");
    error = 0;
    goto label_end;
  }

  /* Grid case: Convert zero values into TEST on the grid */

  st_final_discretization_grid(dbgrid, iptr);

  /* Point case: core reallocation */

  ind = (int*) mem_realloc((char* ) ind, 2 * npair * sizeof(int), 0);
  if (ind == nullptr) goto label_end;

  /* Messages */

  if (flag_verbose)
  {
    message("Total number of samples = %d\n", nech);
    if (flag_same)
      message("Number of samples defined = %d\n", (int) nb);
    else
      message("Number of pairs used for translated correlation = %d\n",
              (int) nb);
  }

  /* Set the error return code */

  error = 0;

  label_end: if (error)
  {
    ind = (int*) mem_free((char* ) ind);
    npair = 0;
  }
  *nindice = npair;
  *indices = ind;
  return (error);
}

/****************************************************************************/
/*!
 **  Identify samples from scatter plot when included within a polygon
 **
 ** \return  Error return code
 **
 ** \param[in]  db1          Db descriptor (first variable)
 ** \param[in]  db2          Db descriptor (second variable for flag.same=T)
 ** \param[in]  icol1        Rank of the first column
 ** \param[in]  icol2        Rank of the second column
 ** \param[in]  polygon      Polygons structure
 **
 ** \remarks The two input Db must match exactly (same number of samples with
 ** \remarks same set of coordinates and same optional selection)
 **
 *****************************************************************************/
int correlation_ident(Db *db1, Db *db2, int icol1, int icol2, Polygons *polygon)
{
  int iech, nech, number;
  double coor[2], val1, val2;

  /* Initializations */

  if (db1 == nullptr) return (1);
  if (db2 == nullptr) return (1);
  nech = db1->getSampleNumber();
  number = 0;

  /* Correlation */

  for (iech = 0; iech < nech; iech++)
  {
    if (!db1->isActive(iech)) continue;
    val1 = db1->getArray(iech, icol1);
    if (FFFF(val1)) continue;
    val2 = db2->getArray(iech, icol2);
    if (FFFF(val2)) continue;

    /* Check of the sample belongs to the polygon */

    coor[0] = val1;
    coor[1] = val2;
    if (!polygon_inside(coor[0], coor[1], TEST, 0, polygon)) continue;

    /* Print the reference of the sample */

    if (number == 0) mestitle(0, "Samples selected from scatter plot");
    message("Sample #%d - Variable #1=%lf - Variable #2=%lf\n", iech + 1, val1,
            val2);
    number++;
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Evaluate the variogram cloud
 **
 ** \param[in]  db      Db descriptor
 ** \param[in]  dbgrid  Discretization Grid descriptor
 ** \param[in]  iptr    Pointer for the variogram cloud (direction)
 ** \param[in]  varioparam VarioParam structure
 ** \param[in]  idir    Rank of the Direction
 **
 *****************************************************************************/
static void st_variogram_cloud(const Db *db,
                               Db *dbgrid,
                               int iptr,
                               const VarioParam *varioparam,
                               int idir)
{
  double ps, psmin, dist, value, w1, w2, z1, z2;
  int nech, iech, jech, igrid, ideb;

  /* Preliminary calculations */

  const DirParam &dirparam = varioparam->getDirParam(idir);
  psmin = _variogram_convert_angular_tolerance(dirparam.getTolAngle());
  nech = db->getSampleNumber();

  /* Loop on the first point */

  for (iech = 0; iech < nech - 1; iech++)
  {
    if (!db->isActive(iech)) continue;
    w1 = db->getWeight(iech);
    if (FFFF(w1)) continue;
    z1 = st_get_IVAR(db, iech, 0);
    if (FFFF(z1)) continue;

    ideb = (st_date_is_used(varioparam, db, db)) ? 0 :
                                                   iech + 1;
    for (jech = ideb; jech < nech; jech++)
    {
      if (!db->isActive(jech)) continue;
      w2 = db->getWeight(jech);
      if (FFFF(w2)) continue;
      z2 = st_get_IVAR(db, jech, 0);
      if (FFFF(z2)) continue;

      /* Check if the pair must be kept (Code criterion) */

      if (code_comparable(db, db, iech, jech, dirparam.getOptionCode(),
                          (int) dirparam.getTolCode())) continue;

      /* Check if the pair must be kept (Date criterion) */

      if (st_date_comparable(varioparam, db, db, iech, jech,
                             dirparam.getIdate())) continue;

      /* Check if the pair must be kept */

      dist = distance_intra(db, iech, jech, NULL);
      if (variogram_reject_pair(db, iech, jech, dist, psmin,
                                dirparam.getBench(), dirparam.getCylRad(),
                                dirparam.getCodir(), &ps)) continue;
      value = w1 * w2 * (z2 - z1) * (z2 - z1) / 2.;
      igrid = st_update_discretization_grid(dbgrid, dist, value);
      if (igrid < 0) continue;
      dbgrid->updArray(igrid, iptr, 0, 1.);
    }
  }
  return;
}

/****************************************************************************/
/*!
 **  Check the samples which are involved in the pairs which are located
 **  within the polygon
 **
 ** \param[in]  db      Db descriptor
 ** \param[in]  dbgrid  Discretization Grid descriptor
 ** \param[in]  vario   Vario structure
 ** \param[in]  polygon Polygons structure
 **
 *****************************************************************************/
void variogram_cloud_ident(Db *db, Db *dbgrid, Vario *vario, Polygons *polygon)
{
  double *ids, *coor, ps, psmin, dist, w1, w2, z1, z2, value, zcoor;
  int *indg, *rank, nech, iech, jech, igrid, idir, ideb;

  /* Initializations */

  indg = rank = nullptr;
  ids = coor = nullptr;
  const VarioParam &varioparam = vario->getVarioParam();

  /* Core allocation */

  nech = db->getSampleNumber();
  indg = db_indg_alloc(dbgrid);
  if (indg == nullptr) goto label_end;
  coor = db_sample_alloc(dbgrid, ELoc::X);
  if (coor == nullptr) goto label_end;
  rank = (int*) mem_alloc(sizeof(int) * nech, 0);
  if (rank == nullptr) goto label_end;
  ids = db_vector_alloc(db);
  if (ids == nullptr) goto label_end;
  for (iech = 0; iech < nech; iech++)
    ids[iech] = 0.;

  /* Loop on the first point */

  for (iech = 0; iech < nech - 1; iech++)
  {
    if (!db->isActive(iech)) continue;
    w1 = db->getWeight(iech);
    if (FFFF(w1)) continue;
    z1 = st_get_IVAR(db, iech, 0);
    if (FFFF(z1)) continue;

    ideb = (st_date_is_used(&varioparam, db, db)) ? 0 :
                                                    iech + 1;
    for (jech = ideb; jech < nech; jech++)
    {
      if (!db->isActive(jech)) continue;
      w2 = db->getWeight(jech);
      if (FFFF(w2)) continue;
      z2 = st_get_IVAR(db, jech, 0);
      if (FFFF(z2)) continue;

      /* Loop on the directions */

      for (idir = 0; idir < vario->getDirectionNumber(); idir++)
      {
        const DirParam &dirparam = vario->getDirParam(idir);

        /* Check if the pair must be kept (Code criterion) */

        if (code_comparable(db, db, iech, jech, dirparam.getOptionCode(),
                            (int) dirparam.getTolCode())) continue;

        /* Check if the pair must be kept (Date criterion) */

        if (st_date_comparable(&varioparam, db, db, iech, jech,
                               dirparam.getIdate())) continue;

        /* Check if the pair must be kept */

        psmin = _variogram_convert_angular_tolerance(dirparam.getTolAngle());
        dist = distance_intra(db, iech, jech, NULL);
        if (variogram_reject_pair(db, iech, jech, dist, psmin,
                                  dirparam.getBench(), dirparam.getCylRad(),
                                  dirparam.getCodir(), &ps)) continue;
        value = w1 * w2 * (z2 - z1) * (z2 - z1) / 2.;
        igrid = st_update_discretization_grid(dbgrid, dist, value);
        if (igrid < 0) continue;

        /* Check if the grid cell belongs to the polygon */

        db_index_sample_to_grid(dbgrid, igrid, indg);
        grid_to_point(dbgrid, indg, NULL, coor);
        zcoor = (dbgrid->getNDim() > 2) ? coor[2] :
                                          TEST;
        if (!polygon_inside(coor[0], coor[1], zcoor, 0, polygon)) continue;

        /* Add the references */

        ids[iech] += 1.;
        ids[jech] += 1.;
      }
    }
  }

  /* Printout the scores: they are ranked by decreasing number */

  mestitle(0, "Samples in variogram cloud (by decreasing order of occurence)");
  for (iech = 0; iech < nech; iech++)
    rank[iech] = iech;
  ut_sort_double(0, nech, rank, ids);

  for (iech = 0; iech < nech; iech++)
  {
    jech = nech - iech - 1;
    if (ids[jech] <= 0.) break;
    message("Sample #%3d: %d occurence(s)\n", rank[jech] + 1, (int) ids[jech]);
  }

  label_end: indg = db_indg_free(indg);
  coor = db_sample_free(coor);
  ids = db_vector_free(ids);
  rank = (int*) mem_free((char* ) rank);
  return;
}

/****************************************************************************/
/*!
 **  Evaluate the variogram cloud
 **
 ** \param[in]  db     Db descriptor
 ** \param[in]  varioparam VarioParam structure
 ** \param[in]  idir   Rank of the Direction
 **
 ** \param[out] vmax   Maximum variogram value
 **
 *****************************************************************************/
static void st_variogram_cloud_dim(Db *db,
                                   const VarioParam *varioparam,
                                   int idir,
                                   double *vmax)
{
  double ps, psmin, dist, value, w1, w2, z1, z2;
  int nech, iech, jech, ideb;

  /* Preliminary calculations */

  const DirParam &dirparam = varioparam->getDirParam(idir);
  psmin = _variogram_convert_angular_tolerance(dirparam.getTolAngle());
  nech = db->getSampleNumber();

  /* Loop on the first point */

  for (iech = 0; iech < nech - 1; iech++)
  {
    if (!db->isActive(iech)) continue;
    if (FFFF(db->getWeight(iech))) continue;

    ideb = (st_date_is_used(varioparam, db, db)) ? 0 :
                                                   iech + 1;
    for (jech = ideb; jech < nech; jech++)
    {
      if (!db->isActive(jech)) continue;
      if (FFFF(db->getWeight(jech))) continue;

      /* Check if the pair must be kept (Code criterion) */

      if (code_comparable(db, db, iech, jech, dirparam.getOptionCode(),
                          (int) dirparam.getTolCode())) continue;

      /* Check if the pair must be kept (Date criterion) */

      if (st_date_comparable(varioparam, db, db, iech, jech,
                             dirparam.getIdate())) continue;

      /* Check if the pair must be kept */

      dist = distance_intra(db, iech, jech, NULL);
      if (variogram_reject_pair(db, iech, jech, dist, psmin,
                                dirparam.getBench(), dirparam.getCylRad(),
                                dirparam.getCodir(), &ps)) continue;
      if (floor(dist / dirparam.getDPas() + 0.5) >= dirparam.getLagNumber())
        continue;

      w1 = db->getWeight(iech);
      w2 = db->getWeight(jech);
      if (FFFF(w1) || FFFF(w2)) continue;
      z1 = st_get_IVAR(db, iech, 0);
      z2 = st_get_IVAR(db, jech, 0);
      if (FFFF(z1) || FFFF(z2)) continue;
      value = w1 * w2 * (z2 - z1) * (z2 - z1) / 2.;
      if (value > (*vmax)) (*vmax) = value;
    }
  }
  return;
}

/****************************************************************************/
/*!
 **  Evaluate the experimental variogram cloud on irregular data
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  varioparam   VarioParam structure
 ** \param[in]  dbgrid       Output grid for storing the variogram cloud
 ** \param[in]  namconv      Naming convention
 **
 *****************************************************************************/
int variogram_cloud(const Db *db,
                    const VarioParam *varioparam,
                    Db *dbgrid,
                    const NamingConvention& namconv)
{
  int idir, iptr;

  /* Initializations */

  if (db == nullptr) return (1);
  if (dbgrid == nullptr) return (1);
  if (varioparam == (VarioParam*) NULL) return (1);
  st_manage_drift_removal(0, NULL, NULL);

  /* Preliminary checks */

  if (db->getNDim() != varioparam->getDimensionNumber())
  {
    messerr("Inconsistent parameters:");
    messerr("Data Base: NDIM=%d", db->getNDim());
    messerr("Variogram: NDIM=%d", varioparam->getDimensionNumber());
    return (1);
  }
  if (!db->isVariableNumberComparedTo(1)) return 1;
  if (dbgrid->getNDim() != 2)
  {
    messerr("The output Db for storing the variogram cloud must be 2-D");
    return (1);
  }

  /* Allocate new variables */

  int ndir = varioparam->getDirectionNumber();
  iptr = dbgrid->addFieldsByConstant(ndir, 0.);
  if (iptr < 0) return (1);

  /* Loop on the directions to evaluate */

  for (idir = 0; idir < ndir; idir++)
  {
    st_variogram_cloud(db, dbgrid, iptr + idir, varioparam, idir);

    /* Convert zero values into TEST */

    st_final_discretization_grid(dbgrid, iptr + idir);
  }

  // Naming of the newly created variables

  namconv.setNamesAndLocators(db, ELoc::Z, -1, dbgrid, iptr, String(), ndir,
                              false);

  return (0);
}

/****************************************************************************/
/*!
 **  Evaluate the bounds for the experimental variogram cloud
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  varioparam   VarioParam structure
 **
 ** \param[out] vmax         Maximum variogram value
 **
 *****************************************************************************/
int variogram_cloud_dim(Db *db, const VarioParam *varioparam, double *vmax)
{
  int idir;

  /* Initializations */

  if (db == nullptr) return (1);
  if (varioparam == (VarioParam*) NULL) return (1);

  /* Preliminary checks */

  if (db->getNDim() != varioparam->getDimensionNumber())
  {
    messerr("Inconsistent parameters:");
    messerr("Data Base: NDIM=%d", db->getNDim());
    messerr("Variogram: NDIM=%d", varioparam->getDimensionNumber());
    return (1);
  }
  if (!db->isVariableNumberComparedTo(1)) return 1;

  /* Loop on the directions to evaluate */

  *vmax = 0.;
  for (idir = 0; idir < varioparam->getDirectionNumber(); idir++)
    st_variogram_cloud_dim(db, varioparam, idir, vmax);

  return (0);
}

/****************************************************************************/
/*!
 **  Evaluate the regression
 **
 ** \return  Error return code
 **
 ** \param[in,out]  db1       Db descriptor (for target variable)
 ** \param[in]  db2           Db descriptor (for auxiliary variables)
 ** \param[in]  flag_mode     Type of calculation
 ** \li                       0 : standard multivariate case
 ** \li                       1 : using external drifts
 ** \li                       2 : using the Model definition
 ** \param[in]  icol          Rank of the target variable
 ** \param[in]  ncol          Number of explanatory variables
 ** \param[in]  icols         Array of ranks of the explanatory variables
 ** \param[in]  model         Model definition (only used if flag_mode == 2)
 ** \param[in]  flag_one      The constant is added as explanatory variable
 ** \param[in]  flag_verbose  0 no output
 ** \li                       1 for the output of the regression
 ** \li                       2 for the dump of the residuals
 **
 ** \param[out]  count     Number of active samples
 ** \param[out]  coeff     Coefficients of the regression
 **                        (Dimension: Number of explanatory variables + 1)
 ** \param[out]  variance  Variance of data
 ** \param[out]  varres    Variance of residuals
 ** \param[out]  correl    Correlation coefficient between response and
 **                        regression variable(s)
 **
 ** \remark  The flag_mode indicates the type of regression calculation:
 ** \remark  0 : V[icol] as a function of V[icols[i]]
 ** \remark  1 : Z1 as a function of the different Fi's
 ** \remark  2 : Z1 as a function of the drift functions
 **
 ** \remark  The Db structure has been modified: the last column of the Db1
 ** \remark  which has been added by this function, contains the value
 ** \remark  of the residuals at each datum (or TEST if the residual has not
 ** \remark  been calculated).
 **
 *****************************************************************************/
int regression_f(Db *db1,
                 Db *db2,
                 int flag_mode,
                 int icol,
                 int ncol,
                 int *icols,
                 Model *model,
                 int flag_one,
                 int flag_verbose,
                 int *count,
                 double *coeff,
                 double *variance,
                 double *varres,
                 double *correl)
{
  int nfex, nech, size, siztri, iech, i, j, error, iptr, ecr, pivot, number,
      nvar;
  int flag_test;
  double *a, *b, *x, value, prod, num, mean, drift, nb, m1, m2, v1, v2, v12,
      rho;
  char string[100];

  /* Initializations */

  error = 1;
  nvar = db1->getVariableNumber();
  nfex = db2->getExternalDriftNumber();
  nech = db1->getSampleNumber();
  size = 0;
  switch (flag_mode)
  {
    case 0:
      size = ncol;
      if (flag_one) size++;
      break;
    case 1:
      size = nfex;
      if (flag_one) size++;
      break;
    case 2:
      size = model->getDriftNumber();
      break;
  }

  siztri = size * (size + 1) / 2;
  a = b = x = nullptr;
  mean = *variance = *varres = 0.;
  *count = 0;
  *variance = *varres = *correl = TEST;
  for (i = 0; i < size; i++)
    coeff[i] = TEST;

  /* Preliminary checks */

  switch (flag_mode)
  {
    case 0:
      if (icol < 0 || icol >= db1->getFieldNumber())
      {
        messerr("The regression requires a valid target variable");
        goto label_end;
      }
      for (i = 0; i < ncol; i++)
      {
        if (icols[i] < 0 || icols[i] >= db2->getFieldNumber())
        {
          messerr("The regression requires a valid auxiliary variable (#%d)",
                  i + 1);
          goto label_end;
        }
      }
      break;

    case 1:
      if (nvar != 1 || nfex <= 0)
      {
        messerr("The multivariate regression is designated for one variable");
        messerr("as a function of several drift variables");
        messerr("The Db contains %d variables and %d drift variables", nvar,
                nfex);
        goto label_end;
      }
      break;

    case 2:
      model_setup(model);
      break;
  }

  /* Pointer allocation */

  iptr = db1->addFieldsByConstant(1, TEST);
  if (iptr < 0) goto label_end;

  /* Core allocation */

  x = (double*) mem_alloc(size * sizeof(double), 0);
  if (x == nullptr) goto label_end;
  b = (double*) mem_alloc(size * sizeof(double), 0);
  if (b == nullptr) goto label_end;
  a = (double*) mem_alloc(siztri * sizeof(double), 0);
  if (a == nullptr) goto label_end;
  for (i = 0; i < size; i++)
    b[i] = 0.;
  for (i = 0; i < siztri; i++)
    a[i] = 0.;
  prod = value = 0.;

  /* Loop on the samples */

  for (iech = number = 0; iech < nech; iech++)
  {
    if (!db1->isActive(iech)) continue;

    /* Get the information for the current sample */

    ecr = 0;
    switch (flag_mode)
    {
      case 0:
        value = db1->getArray(iech, icol);
        if (flag_one) x[ecr++] = 1.;
        for (i = 0; i < ncol; i++)
          x[ecr++] = db2->getArray(iech, icols[i]);
        break;
      case 1:
        value = st_get_IVAR(db1, iech, 0);
        if (flag_one) x[ecr++] = 1.;
        for (i = 0; i < nfex; i++)
          x[ecr++] = db2->getExternalDrift(iech, i);
        break;
      case 2:
        value = st_get_IVAR(db1, iech, 0);
        model_calcul_drift(model, ECalcMember::LHS, db2, iech, x);
        break;
    }

    for (i = flag_test = 0; i < size && !flag_test; i++)
      flag_test = FFFF(x[i]);
    if (FFFF(value) || flag_test) continue;
    prod += value * value;
    mean += value;
    number++;

    /* Update the matrices */

    for (i = ecr = 0; i < size; i++)
    {
      b[i] += value * x[i];
      for (j = 0; j <= i; j++, ecr++)
        a[ecr] += x[i] * x[j];
    }
  }

  if (number <= 0)
  {
    messerr("No sample found where variables are defined");
    error = 0;
    goto label_end;
  }
  mean /= number;
  *count = number;
  *variance = prod / number - mean * mean;

  /* Solve the regression system */

  if (matrix_solve(0, a, b, x, size, 1, &pivot))
  {
    messerr("Error during regression calculation: pivot %d is null", pivot);
    goto label_end;
  }

  /* Load the coefficients */

  for (i = 0; i < size; i++)
    coeff[i] = x[i];

  /* Calculate the residuals */

  for (i = ecr = 0; i < size; i++)
  {
    prod -= 2. * coeff[i] * b[i];
    for (j = 0; j <= i; j++, ecr++)
    {
      num = (i == j) ? 1. :
                       2.;
      prod += coeff[i] * coeff[j] * num * a[ecr];
    }
  }
  *varres = prod / number;

  /************/
  /* Printout */
  /************/

  if (flag_verbose > 0)
  {
    mestitle(1, "Linear regression:");
    for (i = 0; i < size; i++)
      message("Explanatory variable Aux.#%d - Coefficient = %lf\n", i + 1,
              coeff[i]);
    message("Variance of Residuals = %lf\n", (*varres));
  }

  /********************/
  /* Print the header */
  /********************/

  if (flag_verbose == 2)
  {
    message("\n");
    tab_prints(NULL, 1, EJustify::RIGHT, "Rank");
    tab_prints(NULL, 1, EJustify::RIGHT, "Target");
    for (i = 0; i < size; i++)
    {
      (void) gslSPrintf(string, "Aux.#%d", i + 1);
      tab_prints(NULL, 1, EJustify::RIGHT, string);
    }
    tab_prints(NULL, 1, EJustify::RIGHT, "Residuals");
    message("\n");
  }

  /* Store the regression error at sample points */

  nb = m1 = m2 = v1 = v2 = v12 = rho = 0.;
  for (iech = 0; iech < nech; iech++)
  {
    value = TEST;
    if (db1->isActive(iech))
    {

      /* Get the information for the current sample */

      ecr = 0;
      switch (flag_mode)
      {
        case 0:
          value = db1->getArray(iech, icol);
          if (flag_one) x[ecr++] = 1.;
          for (i = 0; i < ncol; i++)
            x[ecr++] = db2->getArray(iech, icols[i]);
          break;
        case 1:
          value = st_get_IVAR(db1, iech, 0);
          if (flag_one) x[ecr++] = 1.;
          for (i = 0; i < nfex; i++)
            x[ecr++] = db2->getExternalDrift(iech, i);
          break;
        case 2:
          value = st_get_IVAR(db1, iech, 0);
          model_calcul_drift(model, ECalcMember::LHS, db2, iech, x);
          break;
      }
      for (i = flag_test = 0; i < size && !flag_test; i++)
        flag_test = FFFF(x[i]);

      if (flag_verbose == 2)
      {
        tab_printi(NULL, 1, EJustify::RIGHT, iech + 1);
        tab_printg(NULL, 1, EJustify::RIGHT, value);
        for (i = 0; i < size; i++)
          tab_printg(NULL, 1, EJustify::RIGHT, x[i]);
      }

      if (FFFF(value) || flag_test)
      {
        value = TEST;
      }
      else
      {
        drift = 0.;
        for (i = 0; i < size; i++)
          drift += x[i] * coeff[i];
        nb += 1.;
        m1 += value;
        m2 += drift;
        v1 += value * value;
        v2 += drift * drift;
        v12 += value * drift;
        value -= drift;
      }

      if (flag_verbose == 2)
      {
        tab_printg(NULL, 1, EJustify::RIGHT, value);
        message("\n");
      }
    }
    db1->setArray(iech, iptr, value);
  }

  /* Calculate the correlation coefficient */

  *correl = 0.;
  if (nb > 0)
  {
    m1 /= nb;
    m2 /= nb;
    v1 = v1 / nb - m1 * m1;
    v2 = v2 / nb - m2 * m2;
    v12 = v12 / nb - m1 * m2;
    rho = v12 / sqrt(v1 * v2);
    *correl = rho;
  }

  /* Set the error returned code */

  error = 0;

  label_end: a = (double*) mem_free((char* ) a);
  b = (double*) mem_free((char* ) b);
  x = (double*) mem_free((char* ) x);
  return (error);
}

/****************************************************************************/
/*!
 **  Ask the characteristics of the Vario structure
 **
 ** \return  Error returned code
 **
 ** \param[in]  vario  Vario structure
 **
 ** \param[out]  calcul_type Type of calculation (ECalcVario)
 ** \param[out]  ndim        Space dimension
 ** \param[out]  nvar        Number of variables
 ** \param[out]  ndir        Number of calculation directions
 ** \param[out]  ndate       Number of Date Intervals
 ** \param[out]  scale       Scaling factor for the transitive covariogram
 ** \param[out]  dates       Array of bounds for Date Intervals
 **
 ** \remark The array 'dates' must be freed by calling function.
 ** \remark The following code shows how to extract the calculation results
 ** \remark from a variogram
 **
 *****************************************************************************/
int vario_extract(Vario *vario,
                  ECalcVario *calcul_type,
                  int *ndim,
                  int *nvar,
                  int *ndir,
                  int *ndate,
                  double *scale,
                  double **dates)
{
  double *date_loc;

  /* Returning arguments */

  *calcul_type = vario->getCalculType();
  *ndim = vario->getDimensionNumber();
  *nvar = vario->getVariableNumber();
  *ndir = vario->getDirectionNumber();
  *ndate = vario->getDateNumber();
  *scale = vario->getScale();
  date_loc = (double*) mem_alloc(sizeof(double) * vario->getDateNumber() * 2,
                                 1);
  int ecr = 0;
  for (int i = 0; i < vario->getDateNumber(); i++)
    for (int icas = 0; icas < 2; icas++)
      date_loc[ecr++] = vario->getDates(i, icas);
  *dates = date_loc;

  return (0);
}

/****************************************************************************/
/*!
 **  Ask for the rank of the 'vardir' structure, given direction and date
 **
 ** \return  Absolute rank (or -1 for error)
 **
 ** \param[in]  vario  Vario structure
 ** \param[in]  idir   Rank for the direction (starting from 0)
 ** \param[in]  idate  Rank for the Date (starting from 0)
 **
 ** \remark  An error occurs if 'idir' is negative or larger than 'ndir'
 ** \remark  or if 'idate' is negative or larger than 'ndate'
 **
 *****************************************************************************/
int vario_get_rank(Vario *vario, int idir, int idate)
{
  int ndir, ndate, rank;

  /* Initializations */

  rank = idir;
  ndir = vario->getDirectionNumber();
  ndate = vario->getDateNumber();
  if (idir < 0 || idir >= ndir) return (-1);
  if (ndate > 0)
  {
    if (idate < 0 || idate >= ndate) return (-1);
    rank = rank * idir + idate;
  }
  return (rank);
}

/****************************************************************************/
/*!
 **  Product of a vector by its conjugate
 **
 ** \param[in]  size      Dimension of the vectors
 ** \param[in]  coef      Multiplicative coefficient
 ** \param[in]  tab1      First complex array
 ** \param[in]  tab2      Second complex array
 **
 ** \param[out]  tab      Output complex array
 **
 *****************************************************************************/
static void st_product_conj(int size,
                            double coef,
                            double *tab1[2],
                            double *tab2[2],
                            double *tab[2])
{
  int i;

  for (i = 0; i < size; i++)
  {
    tab[0][i] += coef * (tab1[0][i] * tab2[0][i] + tab1[1][i] * tab2[1][i]);
    tab[1][i] += coef * (tab1[0][i] * tab2[1][i] - tab1[1][i] * tab2[0][i]);
  }
}

/****************************************************************************/
/*!
 **  Load the resulting array in the output VMAP grid
 **
 ** \param[in]  dbmap     Db containing the VMAP
 ** \param[in]  tab       Array to be stored (in FFT space)
 ** \param[in]  iptr      Pointer for storage
 **
 *****************************************************************************/
static void st_vmap_store(Db *dbmap, double *tab, int iptr)
{
  int ix, iy, iz, ecr, iech, indice[3];
  VectorDouble dims(3);
  int ndim = dbmap->getNDim();

  for (int idim = 0; idim < 3; idim++)
  {
    if (idim < ndim)
      dims[idim] = dbmap->getNX(idim);
    else
      dims[idim] = 1;
  }

  /* Loop on the sample (supposedly ordered) */

  ecr = 0;
  for (ix = 0; ix < dims[0]; ix++)
    for (iy = 0; iy < dims[1]; iy++)
      for (iz = 0; iz < dims[2]; iz++, ecr++)
      {
        indice[0] = ix;
        indice[1] = iy;
        indice[2] = iz;
        iech = db_index_grid_to_sample(dbmap, indice);
        dbmap->setArray(iech, iptr, tab[ecr]);
      }

  return;
}

/****************************************************************************/
/*!
 **  Blank the FFT arrays
 **
 ** \param[in] size    Dimension of the vectors
 ** \param[in] tab     Complex array to be blanked out
 **
 *****************************************************************************/
static void st_vmap_blank(int size, double *tab[2])
{
  int i, ic;

  for (ic = 0; ic < 2; ic++)
    for (i = 0; i < size; i++)
      tab[ic][i] = 0.;
  return;
}

/****************************************************************************/
/*!
 **  Scale the FFT arrays
 **
 ** \param[in] size     Dimension of the vectors
 ** \param[in] scale    Scaling factor
 ** \param[in,out] tab1 First complex array
 ** \param[in] tab2     Second complex array
 **
 *****************************************************************************/
static void st_vmap_rescale(int size, double scale, double *tab1, double *tab2)
{
  int i;
  double value;

  for (i = 0; i < size; i++)
  {
    value = tab2[i];
    if (value > EPSILON8) tab1[i] /= (scale * value);
  }
  return;
}

/****************************************************************************/
/*!
 **  Shift the product of means for the FFT arrays
 **
 ** \param[in] size    Dimension of the vectors
 ** \param[in,out] tab Input/Output complex array
 ** \param[in] tabm1   Complex array for mean of variable 1
 ** \param[in] tabm2   Complex array for mean of variable 2
 **
 *****************************************************************************/
static void st_vmap_shift(int size, double *tab, double *tabm1, double *tabm2)
{
  int i;

  for (i = 0; i < size; i++)
    tab[i] -= tabm1[i] * tabm2[i];

  return;
}

/****************************************************************************/
/*!
 **  Load the data into the FFT arrays
 **
 ** \return  Error returned code
 **
 ** \param[in] dbgrid    Db structure containing the input grid
 ** \param[in] ndim      Space dimension
 ** \param[in] sizetot   Dimension of the vectors
 ** \param[in] dims      Array of dimensions of the extended images
 ** \param[in] dinv      Array of dimensions of the extended images (inverted)
 ** \param[in] ivar      Rank of the first variable
 ** \param[in] jvar      Rank of the second variable
 **
 ** \param[out] i1i1    Pointer to the complex array 1(Z1)
 ** \param[out] z1i1    Pointer to the complex array Z1
 ** \param[out] i2i2    Pointer to the complex array 1(Z2)
 ** \param[out] z2i2    Pointer to the complex array Z2
 ** \param[out] i1i2    Pointer to the complex array 1(Z1) * 1(Z2)
 ** \param[out] z1i2    Pointer on the complex array Z1 * 1(Z2)
 ** \param[out] z2i1    Pointer on the complex array Z2 * 1(Z1)
 ** \param[out] z1z2    Pointer on the complex array Z1 * Z2
 **
 ** \remark The arrays are evaluated only if the input pointer is defined
 **
 *****************************************************************************/
static int st_vmap_load(Db *dbgrid,
                        int ndim,
                        int sizetot,
                        int *dims,
                        int *dinv,
                        int ivar,
                        int jvar,
                        double *i1i1[2],
                        double *z1i1[2],
                        double *i2i2[2],
                        double *z2i2[2],
                        double *i1i2[2],
                        double *z1i2[2],
                        double *z2i1[2],
                        double *z1z2[2])
{
  int ix, iy, iz, iech, ecr, ind1, ind2, indice[3];
  double val1, val2;

  /* Initialize the complex array */

  if (i1i1 != NULL) st_vmap_blank(sizetot, i1i1);
  if (z1i1 != NULL) st_vmap_blank(sizetot, z1i1);
  if (i2i2 != NULL) st_vmap_blank(sizetot, i2i2);
  if (z2i2 != NULL) st_vmap_blank(sizetot, z2i2);
  if (i1i2 != NULL) st_vmap_blank(sizetot, i1i2);
  if (z1i2 != NULL) st_vmap_blank(sizetot, z1i2);
  if (z2i1 != NULL) st_vmap_blank(sizetot, z2i1);
  if (z1z2 != NULL) st_vmap_blank(sizetot, z1z2);

  /* Loop on the grid cells */

  for (ix = ecr = 0; ix < dims[0]; ix++)
    for (iy = 0; iy < dims[1]; iy++)
      for (iz = 0; iz < dims[2]; iz++, ecr++)
      {
        if (ndim >= 1 && ix >= dbgrid->getNX(0)) continue;
        if (ndim >= 2 && iy >= dbgrid->getNX(1)) continue;
        if (ndim >= 3 && iz >= dbgrid->getNX(2)) continue;
        indice[0] = ix;
        indice[1] = iy;
        indice[2] = iz;
        iech = db_index_grid_to_sample(dbgrid, indice);
        if (!dbgrid->getSelection(iech)) continue;
        val1 = st_get_IVAR(dbgrid, iech, jvar);
        val2 = st_get_IVAR(dbgrid, iech, ivar);
        ind1 = (!FFFF(val1));
        ind2 = (!FFFF(val2));

        if (i1i1 != NULL) i1i1[0][ecr] = (ind1);
        if (z1i1 != NULL) z1i1[0][ecr] = (ind1) ? val1 :
                                                  0.;
        if (i2i2 != NULL) i2i2[0][ecr] = (ind2);
        if (z2i2 != NULL) z2i2[0][ecr] = (ind2) ? val2 :
                                                  0.;
        if (i1i2 != NULL) i1i2[0][ecr] = (ind1 && ind2);
        if (z1i2 != NULL) z1i2[0][ecr] = (ind1 && ind2) ? val1 :
                                                          0.;
        if (z2i1 != NULL) z2i1[0][ecr] = (ind1 && ind2) ? val2 :
                                                          0.;
        if (z1z2 != NULL) z1z2[0][ecr] = (ind1 && ind2) ? val1 * val2 :
                                                          0.;
      }

  /* Perform the Fourier transform of the complex arrays defined */

  if (i1i1 != NULL && fftn(ndim, dinv, i1i1[0], i1i1[1], 1, 1.)) return (1);
  if (z1i1 != NULL && fftn(ndim, dinv, z1i1[0], z1i1[1], 1, 1.)) return (1);
  if (i2i2 != NULL && fftn(ndim, dinv, i2i2[0], i2i2[1], 1, 1.)) return (1);
  if (z2i2 != NULL && fftn(ndim, dinv, z2i2[0], z2i2[1], 1, 1.)) return (1);
  if (i1i2 != NULL && fftn(ndim, dinv, i1i2[0], i1i2[1], 1, 1.)) return (1);
  if (z1i2 != NULL && fftn(ndim, dinv, z1i2[0], z1i2[1], 1, 1.)) return (1);
  if (z2i1 != NULL && fftn(ndim, dinv, z2i1[0], z2i1[1], 1, 1.)) return (1);
  if (z1z2 != NULL && fftn(ndim, dinv, z1z2[0], z1z2[1], 1, 1.)) return (1);
  return (0);
}

/*****************************************************************to be updated***********/
/*!
 **  Allocate an array of complex values
 **
 ** \return  Error returned code
 **
 ** \param[in]  size    Dimension of the complex array
 **
 ** \param[out] tab     Complex array to be allocated
 **
 *****************************************************************************/
static int st_complex_array_alloc(int size, double *tab[2])
{
  int ic;

  for (ic = 0; ic < 2; ic++)
    tab[ic] = nullptr;

  for (ic = 0; ic < 2; ic++)
  {
    tab[ic] = (double*) mem_alloc(sizeof(double) * size, 0);
    if (tab[ic] == nullptr) return (1);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Free an array of complex values
 **
 ** \param[out] tab     Complex array to be freed
 **
 *****************************************************************************/
static void st_complex_array_free(double *tab[2])
{
  int ic;

  if (tab == NULL) return;
  for (ic = 0; ic < 2; ic++)
  {
    if (tab[ic] != nullptr) tab[ic] = (double*) mem_free((char* ) tab[ic]);
  }
  return;
}

/****************************************************************************/
/*!
 **  Load the resulting array in the output VMAP grid
 **
 ** \param[in]  nxmap   Dimensions of the Vmap Grid
 ** \param[in]  nxgrid  Dimensions of the Db Grid
 ** \param[in]  dims    Array of dimensions of the extended images
 ** \param[in]  tabin   Array to be stored (in FFT space)
 **
 ** \param[out] tabout  Array containing the resulting VMAP
 **
 *****************************************************************************/
static void st_vmap_extract(int *nxmap,
                            int *nxgrid,
                            int *dims,
                            double *tabin,
                            double *tabout)
{
  int ix, iy, iz, jx, jy, jz, nxs2, nys2, nzs2, nxloc, nyloc, nzloc;

  /* Initializations */

  nxs2 = nxmap[0] / 2;
  nys2 = nxmap[1] / 2;
  nzs2 = nxmap[2] / 2;
  nxloc = MIN(1 + nxs2, nxgrid[0]);
  nyloc = MIN(1 + nys2, nxgrid[1]);
  nzloc = MIN(1 + nzs2, nxgrid[2]);

  /* Fill the array (IX+,IY+,IZ+) */

  for (ix = 0; ix < nxloc; ix++)
    for (iy = 0; iy < nyloc; iy++)
      for (iz = 0; iz < nzloc; iz++)
      {
        jx = nxs2 + ix;
        jy = nys2 + iy;
        jz = nzs2 + iz;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(ix, iy, iz, dims)];
      }

  /* Fill the array (IX-,IY+,IZ+) */

  for (ix = 0; ix < nxloc - 1; ix++)
    for (iy = 0; iy < nyloc; iy++)
      for (iz = 0; iz < nzloc; iz++)
      {
        jx = nxs2 - ix - 1;
        jy = nys2 + iy;
        jz = nzs2 + iz;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(OPP(0,ix), iy, iz, dims)];
      }

  /* Fill the array (IX+,IY-,IZ+) */

  for (ix = 0; ix < nxloc; ix++)
    for (iy = 0; iy < nyloc - 1; iy++)
      for (iz = 0; iz < nzloc; iz++)
      {
        jx = nxs2 + ix;
        jy = nys2 - iy - 1;
        jz = nzs2 + iz;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(ix, OPP(1,iy), iz, dims)];
      }

  /* Fill the array (IX-,IY-,IZ+) */

  for (ix = 0; ix < nxloc - 1; ix++)
    for (iy = 0; iy < nyloc - 1; iy++)
      for (iz = 0; iz < nzloc; iz++)
      {
        jx = nxs2 - ix - 1;
        jy = nys2 - iy - 1;
        jz = nzs2 + iz;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(OPP(0,ix), OPP(1,iy), iz,
                                                   dims)];
      }

  /* Fill the array (IX+,IY+,IZ-) */

  for (ix = 0; ix < nxloc; ix++)
    for (iy = 0; iy < nyloc; iy++)
      for (iz = 0; iz < nzloc - 1; iz++)
      {
        jx = nxs2 + ix;
        jy = nys2 + iy;
        jz = nzs2 - iz - 1;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(ix, iy, OPP(2,iz), dims)];
      }

  /* Fill the array (IX-,IY+,IZ-) */

  for (ix = 0; ix < nxloc - 1; ix++)
    for (iy = 0; iy < nyloc; iy++)
      for (iz = 0; iz < nzloc - 1; iz++)
      {
        jx = nxs2 - ix - 1;
        jy = nys2 + iy;
        jz = nzs2 - iz - 1;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(OPP(0,ix), iy, OPP(2,iz),
                                                   dims)];
      }

  /* Fill the array (IX+,IY-,IZ-) */

  for (ix = 0; ix < nxloc; ix++)
    for (iy = 0; iy < nyloc - 1; iy++)
      for (iz = 0; iz < nzloc - 1; iz++)
      {
        jx = nxs2 + ix;
        jy = nys2 - iy - 1;
        jz = nzs2 - iz - 1;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(ix, OPP(1,iy), OPP(2,iz),
                                                   dims)];
      }

  /* Fill the array (IX-,IY-,IZ-) */

  for (ix = 0; ix < nxloc - 1; ix++)
    for (iy = 0; iy < nyloc - 1; iy++)
      for (iz = 0; iz < nzloc - 1; iz++)
      {
        jx = nxs2 - ix - 1;
        jy = nys2 - iy - 1;
        jz = nzs2 - iz - 1;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(OPP(0,ix), OPP(1,iy),
                                                   OPP(2,iz), dims)];
      }

  return;
}

/*****************************************************************************/
/*!
 **  Calculates the 2-D variogram map on a grid using FFT
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid       Db of Grid type containing the data
 ** \param[in]  dbmap        Db of Grid type containing the Variogram Maps
 ** \param[in]  calcul_type  Type of calculation (ECalcVario)
 ** \param[in]  namconv      Naming convention
 **
 *****************************************************************************/
static int st_vmap_grid_fft(Db *dbgrid,
                            Db *dbmap,
                            const ECalcVario &calcul_type,
                            const NamingConvention& namconv)
{
  int dims[3], dinv[3], nxmap[3], nxgrid[3], sizetot, sizemap, sizegrid;
  int i, ndim, ic, error, idim, ivar, jvar, ijvar, nvar, nvs2;
  double *i1i1[2], *z1i1[2], *i2i2[2], *z2i2[2], *i1i2[2];
  double *z1i2[2], *z2i1[2], *z1z2[2], *ztab[2];
  double *res_nn, *res_gg, *res_m1, *res_m2;
  static int verbose = 0;

  /* Initializations */

  error = 1;
  sizetot = 0;
  res_nn = res_gg = res_m1 = res_m2 = nullptr;
  for (ic = 0; ic < 2; ic++)
  {
    i1i1[ic] = z1i1[ic] = i2i2[ic] = z2i2[ic] = nullptr;
    i1i2[ic] = z1i2[ic] = z2i1[ic] = z1z2[ic] = nullptr;
    ztab[ic] = nullptr;
  }

  /* Preliminary checks */

  if (dbgrid == nullptr) return (1);
  if (dbmap == nullptr) return (1);

  if (calcul_type != ECalcVario::VARIOGRAM && calcul_type
      != ECalcVario::COVARIOGRAM
      && calcul_type != ECalcVario::COVARIANCE
      && calcul_type != ECalcVario::COVARIANCE_NC)
  {
    messerr("This function is limited to the calculation of");
    messerr("Variogram, Covariance (centered or not) or Covariogram");
    return (1);
  }
  if (!is_grid(dbgrid))
  {
    messerr("This Variogram Map is defined for Grid Data Base only");
    return (1);
  }
  if (!is_grid(dbmap))
  {
    messerr("This feature requires a Grid Data Base");
    messerr("to store the Variogram Maps");
    return (1);
  }
  if (dbgrid->getNDim() != 2 && dbgrid->getNDim() != 3)
  {
    messerr("The Variogram Map can only be calculated on a grid data set");
    messerr("with dimension equal to 2 or 3");
    return (1);
  }
  if (dbmap->getNDim() > dbgrid->getNDim())
  {
    messerr("The space dimension of the VMAP (%d)", dbmap->getNDim());
    messerr(
        "must not be larger than the space dimension of the input Grid (%d)",
        dbgrid->getNDim());
    return (1);
  }
  for (idim = 0; idim < dbmap->getNDim(); idim++)
  {
    if (ABS(dbmap->getDX(idim) - dbgrid->getDX(idim)) > 1.e-03)
    {
      messerr("The grid mesh in the direction %d (dx=%lf)", idim + 1,
              dbgrid->getDX(idim));
      messerr("must match the mesh of the Variogram Map grid (dx=%lf)",
              dbmap->getDX(idim));
      return (1);
    }
  }

  for (idim = 0; idim < 3; idim++)
    nxgrid[idim] = nxmap[idim] = 1;
  for (idim = 0; idim < dbgrid->getNDim(); idim++)
    nxgrid[idim] = dbgrid->getNX(idim);
  for (idim = 0; idim < dbmap->getNDim(); idim++)
    nxmap[idim] = dbmap->getNX(idim);

  /* Preliminary calculations */

  nvar = dbgrid->getVariableNumber();
  nvs2 = nvar * (nvar + 1) / 2;
  ndim = 0;
  sizetot = sizemap = sizegrid = 1;
  for (i = 0; i < 3; i++)
  {
    dinv[i] = 1;
    if (nxgrid[i] <= 1)
    {
      dims[i] = 1;
    }
    else
    {
      dims[i] = (int) ceil((double) (nxgrid[i] + nxmap[i] - 1) / 8.) * 8;
      sizegrid *= nxgrid[i];
      sizemap *= nxmap[i];
      sizetot *= dims[i];
      ndim++;
    }
  }
  for (i = 0; i < ndim; i++)
    dinv[i] = dims[ndim - i - 1];
  if (verbose)
  {
    mestitle(0, "Simulation of a grid using FFT");
    message("Grid dimension:");
    for (idim = 0; idim < ndim; idim++)
      message(" %4d", nxgrid[idim]);
    message("\n");
    message("Variogram Map :");
    for (idim = 0; idim < ndim; idim++)
      message(" %4d", nxmap[idim]);
    message("\n");
    message("Working array :");
    for (idim = 0; idim < ndim; idim++)
      message(" %4d", dinv[idim]);
    message("\n");
  }

  /* Create the variables in the Variogram Map file */

  IPTV = dbmap->addFieldsByConstant(nvs2, TEST);
  if (IPTV < 0) goto label_end;
  IPTW = dbmap->addFieldsByConstant(nvs2, TEST);
  if (IPTW < 0) goto label_end;

  /* Core allocation */

  if (st_complex_array_alloc(sizetot, ztab)) goto label_end;
  if (calcul_type == ECalcVario::VARIOGRAM)
  {
    if (st_complex_array_alloc(sizetot, i1i2)) goto label_end;
    if (st_complex_array_alloc(sizetot, z1i2)) goto label_end;
    if (st_complex_array_alloc(sizetot, z2i1)) goto label_end;
    if (st_complex_array_alloc(sizetot, z1z2)) goto label_end;
  }
  else
  {
    if (st_complex_array_alloc(sizetot, i1i1)) goto label_end;
    if (st_complex_array_alloc(sizetot, z1i1)) goto label_end;
    if (st_complex_array_alloc(sizetot, i2i2)) goto label_end;
    if (st_complex_array_alloc(sizetot, z2i2)) goto label_end;
  }

  res_nn = (double*) mem_alloc(sizeof(double) * sizemap, 0);
  if (res_nn == nullptr) goto label_end;
  res_gg = (double*) mem_alloc(sizeof(double) * sizemap, 0);
  if (res_gg == nullptr) goto label_end;
  for (i = 0; i < sizemap; i++)
    res_nn[i] = res_gg[i] = 0.;
  if (calcul_type == ECalcVario::COVARIANCE || calcul_type
      == ECalcVario::COVARIOGRAM)
  {
    res_m1 = (double*) mem_alloc(sizeof(double) * sizemap, 0);
    if (res_m1 == nullptr) goto label_end;
    res_m2 = (double*) mem_alloc(sizeof(double) * sizemap, 0);
    if (res_m2 == nullptr) goto label_end;
    for (i = 0; i < sizemap; i++)
      res_m1[i] = res_m2[i] = 0.;
  }

  /* Loop on the variables */

  for (ivar = ijvar = 0; ivar < nvar; ivar++)
    for (jvar = 0; jvar <= ivar; jvar++, ijvar++)
    {

      /* Calculate the structural function */

      if (calcul_type == ECalcVario::VARIOGRAM)
      {
        if (st_vmap_load(dbgrid, ndim, sizetot, dims, dinv, ivar, jvar,
        NULL,
                         NULL, NULL, NULL, i1i2, z1i2, z2i1, z1z2)) continue;

        /* Calculate the number of pairs */
        st_vmap_blank(sizetot, ztab);
        st_product_conj(sizetot, 1., i1i2, i1i2, ztab);
        if (fftn(ndim, dinv, ztab[0], ztab[1], -1, -1.)) continue;
        st_vmap_extract(nxmap, nxgrid, dims, ztab[0], res_nn);

        /* Structural function */
        st_vmap_blank(sizetot, ztab);
        st_product_conj(sizetot, 1., z1z2, i1i2, ztab);
        st_product_conj(sizetot, 1., i1i2, z1z2, ztab);
        st_product_conj(sizetot, -1., z1i2, z2i1, ztab);
        st_product_conj(sizetot, -1., z2i1, z1i2, ztab);
        if (fftn(ndim, dinv, ztab[0], ztab[1], -1, -1.)) continue;
        st_vmap_extract(nxmap, nxgrid, dims, ztab[0], res_gg);
        st_vmap_rescale(sizemap, 2., res_gg, res_nn);
      }
      else
      {
        if (st_vmap_load(dbgrid, ndim, sizetot, dims, dinv, ivar, jvar, i1i1,
                         z1i1, i2i2, z2i2, NULL, NULL, NULL, NULL)) continue;

        /* Calculate the number of pairs */
        st_vmap_blank(sizetot, ztab);
        st_product_conj(sizetot, 1., i1i1, i2i2, ztab);
        if (fftn(ndim, dinv, ztab[0], ztab[1], -1, -1.)) continue;
        st_vmap_extract(nxmap, nxgrid, dims, ztab[0], res_nn);

        if (calcul_type == ECalcVario::COVARIANCE || calcul_type
            == ECalcVario::COVARIOGRAM)
        {
          /* Calculate the means */
          st_vmap_blank(sizetot, ztab);
          st_product_conj(sizetot, 1., z1i1, i2i2, ztab);
          if (fftn(ndim, dinv, ztab[0], ztab[1], -1, -1.)) continue;
          st_vmap_extract(nxmap, nxgrid, dims, ztab[0], res_m1);
          st_vmap_rescale(sizemap, 1., res_m1, res_nn);

          st_vmap_blank(sizetot, ztab);
          st_product_conj(sizetot, 1., i1i1, z2i2, ztab);
          if (fftn(ndim, dinv, ztab[0], ztab[1], -1, -1.)) continue;
          st_vmap_extract(nxmap, nxgrid, dims, ztab[0], res_m2);
          st_vmap_rescale(sizemap, 1., res_m2, res_nn);
        }

        /* Structural function */
        st_vmap_blank(sizetot, ztab);
        st_product_conj(sizetot, 1., z1i1, z2i2, ztab);
        if (fftn(ndim, dinv, ztab[0], ztab[1], -1, -1.)) continue;
        st_vmap_extract(nxmap, nxgrid, dims, ztab[0], res_gg);
        st_vmap_rescale(sizemap, 1., res_gg, res_nn);

        if (calcul_type == ECalcVario::COVARIANCE || calcul_type
            == ECalcVario::COVARIOGRAM)
          st_vmap_shift(sizemap, res_gg, res_m1, res_m2);
      }

      /* Store the results */
      st_vmap_store(dbmap, res_nn, IPTW + ijvar);
      st_vmap_store(dbmap, res_gg, IPTV + ijvar);
    }

  /* Set the error return code */

  namconv.setNamesAndLocators(dbgrid, ELoc::Z, -1, dbmap, IPTW, "Nb", 1, false);
  namconv.setNamesAndLocators(dbgrid, ELoc::Z, -1, dbmap, IPTV, "Var", 1);
  error = 0;

  label_end: st_complex_array_free(ztab);
  if (calcul_type == ECalcVario::VARIOGRAM)
  {
    st_complex_array_free(i1i2);
    st_complex_array_free(z1i2);
    st_complex_array_free(z2i1);
    st_complex_array_free(z1z2);
  }
  else
  {
    st_complex_array_free(i1i1);
    st_complex_array_free(z1i1);
    st_complex_array_free(i2i2);
    st_complex_array_free(z2i2);
  }
  res_nn = (double*) mem_free((char* ) res_nn);
  res_gg = (double*) mem_free((char* ) res_gg);
  res_m1 = (double*) mem_free((char* ) res_m1);
  res_m1 = (double*) mem_free((char* ) res_m2);
  return (error);
}

/****************************************************************************/
/*!
 **  Copy a direction from one Vario to another Vario
 **
 ** \param[in]  vario_in     Input Vario structure
 ** \param[in]  idir_in      Rank of the Input Direction
 ** \param[in]  vario_out    Output Vario structure
 ** \param[in]  idir_out     Rank of the Output Direction
 **
 *****************************************************************************/
void vardir_copy(VarioParam *vario_in,
                 int idir_in,
                 VarioParam *vario_out,
                 int idir_out)
{
  if (vario_in == (VarioParam*) NULL) return;
  if (idir_in < 0 || idir_in >= vario_in->getDirectionNumber()) return;
  if (vario_out == (VarioParam*) NULL) return;
  if (idir_out < 0 || idir_out >= vario_in->getDirectionNumber()) return;

  DirParam dir_in = vario_in->getDirParam(idir_in);
  DirParam dir_out = vario_out->getDirParam(idir_out);
  dir_out = dir_in;
}

/****************************************************************************/
/*!
 **  Deallocate a PCA
 **
 ** \return  Pointer to the newly freed PCA structure
 **
 ** \param[in]  pca          PCA structure to be freed
 **
 *****************************************************************************/
PCA* pca_free(PCA *pca)

{
  if (pca == nullptr) return (pca);
  free(pca);
  pca = nullptr;
  return (pca);
}

/****************************************************************************/
/*!
 **  Allocate a PCA
 **
 ** \return  Pointer to the newly allocated PCA structure or NULL (if problem)
 **
 ** \param[in]  nvar         Number of variables
 **
 *****************************************************************************/
PCA* pca_alloc(int nvar)
{
  PCA *pca;

  /* Initializations */

  pca = new (PCA);
  pca->init(nvar);

  return (pca);
}

/****************************************************************************/
/*!
 **  Normalize the isotropic array of values
 **
 ** \param[in]  nvar         Number of variables
 ** \param[in,out] data      Array of information
 ** \param[in]  mean         Array containing the mean
 ** \param[in]  sigma        Array containing the standard deviation
 **
 *****************************************************************************/
static void st_center(int nvar,
                      double *data,
                      const VectorDouble &mean,
                      const VectorDouble &sigma)
{
  int ivar;

  for (ivar = 0; ivar < nvar; ivar++)
  {
    if (sigma[ivar] <= 0.) continue;
    data[ivar] = (data[ivar] - mean[ivar]) / sigma[ivar];
  }
}

/****************************************************************************/
/*!
 **  Fill the mean and variance arrays
 **
 ** \param[in] flag_normalize 1 if the normalization must be calculated
 ** \param[in] verbose     Verbose flag
 ** \param[in] db          Db descriptor
 ** \param[in] data        Array containing variables for one sample
 **
 ** \param[out] mean       Array of means
 ** \param[out] sigma      Array of standard deviations
 **
 *****************************************************************************/
static void st_calculate_normalization(int flag_normalize,
                                       int verbose,
                                       Db *db,
                                       double *data,
                                       VectorDouble &mean,
                                       VectorDouble &sigma)
{
  int niso, nvar, nech;

  // Initializations

  niso = 0;
  nvar = db->getVariableNumber();
  nech = db->getSampleNumber();

  if (flag_normalize)
  {
    /* Initializations */

    for (int ivar = 0; ivar < nvar; ivar++)
      mean[ivar] = sigma[ivar] = 0.;

    /* Calculate the statistics */

    for (int iech = 0; iech < nech; iech++)
    {
      if (!db_is_isotropic(db, iech, data)) continue;
      niso++;
      for (int ivar = 0; ivar < nvar; ivar++)
      {
        mean[ivar] += data[ivar];
        sigma[ivar] += data[ivar] * data[ivar];
      }
    }

    /* Normalization */

    if (niso > 0)
    {
      for (int ivar = 0; ivar < nvar; ivar++)
      {
        mean[ivar] /= niso;
        sigma[ivar] = (sigma[ivar] / niso - mean[ivar] * mean[ivar]);
        sigma[ivar] = (sigma[ivar] > 0) ? sqrt(sigma[ivar]) :
                                          0.;
      }
    }

    if (verbose)
    {
      mestitle(1, "Normalization");
      message("Number of variables         = %d\n", nvar);
      message("Number of samples in the Db = %d\n", nech);
      message("Number of isotropic samples = %d\n", niso);
      message("\n");
    }
  }
  else
  {

    for (int ivar = 0; ivar < nvar; ivar++)
    {
      mean[ivar] = 0.;
      sigma[ivar] = 1.;
    }
    if (verbose)
    {
      mestitle(1, "Normalization");
      message("Normalization has been bypassed\n");
      message("\n");
    }
  }
}

/****************************************************************************/
/*!
 **  Internal PCA covariance calculation
 **
 ** \return The array of covariance (or NULL)
 **
 ** \param[in]  verbose     Verbose flag
 ** \param[in]  db          Db descriptor
 ** \param[in]  mean        Array containing the mean
 ** \param[in]  sigma       Array containing the standard deviation
 **
 ** \param[out] data1       Array containing variables for one sample
 **
 *****************************************************************************/
static VectorDouble st_pca_covariance0(int verbose,
                                       Db *db,
                                       VectorDouble &mean,
                                       VectorDouble &sigma,
                                       double *data1)
{
  int niso, nvar, nech;
  VectorDouble c0;

  /* Cleaning the previous contents of PCA structure */

  nvar = db->getVariableNumber();
  nech = db->getSampleNumber();
  niso = 0;

  /* Core allocation */

  c0.resize(nvar * nvar, 0);

  /* Calculate the variance-covariance matrix at distance 0 */

  niso = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    if (!db_is_isotropic(db, iech, data1)) continue;
    niso++;
    st_center(nvar, data1, mean, sigma);
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
        c0[ivar * nvar + jvar] += data1[ivar] * data1[jvar];
  }

  /* Normalization */

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      c0[ivar * nvar + jvar] /= niso;

  /* Printout of the covariance matrix (optional) */

  if (verbose)
  {
    mestitle(1, "PCA calculations");
    message("Number of variables         = %d\n", nvar);
    message("Number of samples in the Db = %d\n", nech);
    message("Number of isotropic samples = %d\n", niso);
    message("\n");
    print_matrix("Variance-Covariance matrix for distance 0", 0, 1, nvar, nvar,
    NULL,
                 c0.data());
  }
  return (c0);
}

/****************************************************************************/
/*!
 **  Internal PCA translated covariance calculation
 **
 ** \return The array of covariance (or NULL)
 **
 ** \param[in]  verbose     Verbose flag
 ** \param[in]  db          Db descriptor
 ** \param[in]  opt_code    code selection option
 ** \li                      0 : no use of the code selection
 ** \li                      1 : codes must be close enough
 ** \li                      2 : codes must be different
 ** \param[in]  tolcode     Code tolerance
 ** \param[in]  codir       Direction
 ** \param[in]  tolang      Angular tolerance
 ** \param[in]  bench       Slicing bench
 ** \param[in]  cylrad      Slicing radius
 ** \param[in]  h0          Reference distance
 ** \param[in]  dh          Tolerance on distance
 **
 ** \param[out] data1       Array containing variables for one sample
 ** \param[out] data2       Array containing variables for one sample
 **
 *****************************************************************************/
static VectorDouble st_pca_covarianceh(int verbose,
                                       Db *db,
                                       int opt_code,
                                       double tolcode,
                                       VectorDouble &codir,
                                       double tolang,
                                       double bench,
                                       double cylrad,
                                       double h0,
                                       double dh,
                                       double *data1,
                                       double *data2)
{
  int npairs, nech, nvar;
  double dist, psmin, ps, di, dj;
  VectorDouble ch;

  /* Initializations */

  nech = db->getSampleNumber();
  nvar = db->getVariableNumber();
  npairs = 0;
  psmin = _variogram_convert_angular_tolerance(tolang);

  /* Core allocation */

  ch.resize(nvar * nvar, 0);

  for (int iech = 0; iech < nech; iech++)
  {
    if (!db_is_isotropic(db, iech, data1)) continue;

    /* Loop on the second sample */

    for (int jech = 0; jech < iech; jech++)
    {
      if (!db_is_isotropic(db, jech, data2)) continue;

      /* Should the pair be retained */

      dist = distance_intra(db, iech, jech, NULL);
      if (dist < h0 - dh || dist > h0 + dh) continue;
      if (code_comparable(db, db, iech, jech, opt_code, (int) tolcode))
        continue;
      if (variogram_reject_pair(db, iech, jech, dist, psmin, bench, cylrad,
                                codir, &ps)) continue;

      /* Update the variance-covariance matrix at distance h */

      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar < nvar; jvar++)
        {
          di = data1[ivar] - data2[ivar];
          dj = data1[jvar] - data2[jvar];
          ch[ivar * nvar + jvar] += di * dj;
        }
      npairs++;
    }
  }

  /* Normation */

  if (npairs <= 1)
  {
    messerr("Number of pairs of isotropic samples is smaller than 2");
    return (ch);
  }
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      ch[ivar * nvar + jvar] /= npairs;

  /* Verbose printout */

  if (verbose)
  {
    mestitle(1, "MAF calculations");
    st_vario_params_print(db->getNDim(), codir, tolang, bench, cylrad);
    message("Reference Distance          = %lf\n", h0);
    message("Tolerance on Distance       = %lf\n", dh);
    message("Number of variables         = %d\n", nvar);
    message("Number of samples in the Db = %d\n", nech);
    message("Number of isotropic pairs   = %d\n", npairs);
    message("\n");
    print_matrix("Variance-Covariance matrix for distance h", 0, 1, nvar, nvar,
    NULL,
                 ch.data());
  }
  return (ch);
}

/****************************************************************************/
/*!
 **  Procedure for transforming the factors into variables using PCA
 **
 ** \param[in]  flag_norm_out 1 if the output variable must be denormalized
 ** \param[in]  iptr         Pointer to the storage
 ** \param[in]  db           Db descriptor
 ** \param[in]  pca          PCA descriptor
 ** \param[in]  data1        Input array
 **
 *****************************************************************************/
static void st_pca_f2z(int flag_norm_out,
                       int iptr,
                       Db *db,
                       PCA *pca,
                       double *data1)
{
  int nvar, nech, flag_anisotropy;
  double value;

  /* Initializations */

  nvar = db->getVariableNumber();
  nech = db->getSampleNumber();

  /* Loop on the samples */

  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    flag_anisotropy = !db_is_isotropic(db, iech, data1);

    /* Loop on the factors */

    for (int ivar = 0; ivar < nvar; ivar++)
    {
      if (flag_anisotropy)
        value = TEST;
      else
      {
        value = 0.;
        for (int ifac = 0; ifac < nvar; ifac++)
          value += pca->getF2Z(ivar, ifac) * data1[ifac];
        if (flag_norm_out)
          value = pca->getMean(ivar) + pca->getSigma(ivar) * value;
      }
      db->setArray(iech, ivar + iptr, value);
    }
  }
}

/****************************************************************************/
/*!
 **  Procedure for transforming the variables into factors using PCA
 **
 ** \param[in]  flag_norm_out 1 if the factors must be normalized
 ** \param[in]  iptr         Pointer for storing the result in db
 ** \param[in]  db           Db descriptor
 ** \param[in]  pca          PCA descriptor
 ** \param[in]  data1        Input array
 ** \param[in]  mean         Array containing the mean
 ** \param[in]  sigma        Array containing the standard deviation
 **
 ** \remarks The standardization statistics are passed as arguments
 ** \remarks The ones stored in PCA structure are not used.
 **
 *****************************************************************************/
static void st_pca_z2f(int flag_norm_out,
                       int iptr,
                       Db *db,
                       PCA *pca,
                       double *data1,
                       const VectorDouble &mean,
                       const VectorDouble &sigma)
{
  int nvar, nech, flag_anisotropy;
  double value;

  /* Initializations */

  nvar = db->getVariableNumber();
  nech = db->getSampleNumber();

  /* Loop on the samples */

  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    flag_anisotropy = !db_is_isotropic(db, iech, data1);
    st_center(nvar, data1, mean, sigma);

    /* Loop on the factors */

    for (int ifac = 0; ifac < nvar; ifac++)
    {
      if (flag_anisotropy)
        value = TEST;
      else
      {
        value = 0.;
        for (int ivar = 0; ivar < nvar; ivar++)
          value += pca->getZ2F(ifac, ivar) * data1[ivar];
        if (flag_norm_out) value /= sqrt(pca->getEigen(ifac));
      }
      db->setArray(iech, ifac + iptr, value);
    }
  }
}

/****************************************************************************/
/*!
 **  Internal function to calculate MAF
 **
 ** \return  Error return code
 **
 ** \param[in]  flag_norm  1 if the normalization must be calculated
 ** \param[in]  db         Db descriptor
 ** \param[in]  data1      Data vector
 ** \param[in]  verbose    Verbose flag
 **
 ** \param[out] pca        Output PCA structure
 **
 ** \remarks Note that the ELoc::Z is redefined in this function
 **
 *****************************************************************************/
static int st_pca_calculate(int flag_norm,
                            Db *db,
                            double *data1,
                            PCA *pca,
                            int verbose)
{
  VectorDouble c0;
  int error;
  VectorDouble mean, sigma;

  // Initializations

  error = 1;
  int nvar = pca->getNVar();

  // Preparing

  pca->clean();
  mean.resize(nvar, 0);
  sigma.resize(nvar, 0);
  st_calculate_normalization(flag_norm, verbose, db, data1, mean, sigma);
  c0 = st_pca_covariance0(verbose, db, mean, sigma, data1);
  pca->setMean(mean);
  pca->setSigma(sigma);

  // Processing

  if (pca->calculateEigen(nvar, c0)) goto label_end;
  if (verbose) pca->display(1, 1);

  // Set the error return code

  error = 0;

  label_end: return (error);
}

/****************************************************************************/
/*!
 **  Evaluate the MAF on irregular data
 **
 ** \return  Error return code
 **
 ** \param[in]  db         Db descriptor
 ** \param[in]  opt_code   code selection option
 ** \li                     0 : no use of the code selection
 ** \li                     1 : codes must be close enough
 ** \li                     2 : codes must be different
 ** \param[in]  tolcode    Code tolerance
 ** \param[in]  codir      Direction
 ** \param[in]  tolang     Angular tolerance
 ** \param[in]  bench      Slicing bench
 ** \param[in]  cylrad     Slicing radius
 ** \param[in]  h0         Reference distance
 ** \param[in]  dh         Tolerance on distance
 ** \param[in]  verbose    Verbose flag
 **
 ** \param[out] pca        Output PCA structure
 **
 ** \remarks Note that the ELoc::Z is redefined in this function
 **
 *****************************************************************************/
int maf_compute(Db *db,
                int opt_code,
                double tolcode,
                VectorDouble &codir,
                double tolang,
                double bench,
                double cylrad,
                double h0,
                double dh,
                int verbose,
                PCA *pca)
{
  int iptr, error, nvar, kvar;
  double *data1, *data2, *identity;
  PCA *pca2;
  VectorDouble mean, sigma, pcaz2f, pca2z2f, pcaf2z, ch;
  static int flag_norm = 1;

  /* Initializations */

  error = 0;
  iptr = -1;
  if (db == nullptr) return (1);
  data1 = data2 = identity = nullptr;
  pca2 = nullptr;

  /* Preliminary checks */

  nvar = pca->getNVar();
  if (nvar != db->getVariableNumber())
  {
    messerr(
        "Number of variables in the PCA (%d) and in the Db (%d) are inconsistent",
        nvar, db->getVariableNumber());
    goto label_end;
  }

  /* Allocate new variables */

  iptr = db->addFieldsByConstant(nvar, TEST);
  if (iptr < 0) goto label_end;

  /* Core allocation */

  data1 = (double*) mem_alloc(sizeof(double) * nvar, 0);
  if (data1 == nullptr) goto label_end;
  data2 = (double*) mem_alloc(sizeof(double) * nvar, 0);
  if (data2 == nullptr) goto label_end;
  identity = (double*) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (identity == nullptr) goto label_end;
  for (int i = 0; i < nvar * nvar; i++)
    identity[i] = 0.;

  /* Calculate the first PCA (centered and normalized) */

  if (st_pca_calculate(flag_norm, db, data1, pca, verbose)) goto label_end;

  /* Rotate the initial data in the PCA system */

  st_pca_z2f(1, iptr, db, pca, data1, pca->getMean(), pca->getSigma());
  db->setLocatorsByAttribute(nvar, iptr, ELoc::Z);

  /* Calculate the variance-covariance matrix at distance [h0-dh,h0+dh] */

  pca2 = pca_alloc(nvar);
  if (pca2 == nullptr) goto label_end;
  ch = st_pca_covarianceh(verbose, db, opt_code, tolcode, codir, tolang, bench,
                          cylrad, h0, dh, data1, data2);
  if (pca2->calculateEigen(nvar, ch)) goto label_end;
  if (verbose) pca2->display(0, 1);

  /* Rotate the initial data in the second PCA system */

  for (int ivar = 0; ivar < nvar; ivar++)
    identity[ivar * nvar + ivar] = 1. / sqrt(pca->getEigen(ivar));
  pcaz2f = pca->getZ2F();
  pca2z2f = pca2->getZ2F();
  pcaf2z = pca->getF2Z();
  matrix_product(nvar, nvar, nvar, pcaz2f.data(), identity, pcaz2f.data());
  matrix_product(nvar, nvar, nvar, pcaz2f.data(), pca2z2f.data(),
                 pca2z2f.data());
  pca->setPcaZ2F(pcaz2f);
  pca2->setPcaZ2F(pca2z2f);
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    kvar = nvar - ivar - 1;
    pca->setEigen(kvar, pca2->getEigen(ivar));
    for (int jvar = 0; jvar < nvar; jvar++)
      pca->setPcaZ2F(kvar, jvar, pca2->getZ2F(ivar, jvar));
  }
  if (matrix_invert_copy(pca->getZ2F().data(), nvar, pcaf2z.data()))
    goto label_end;
  pca->setPcaF2Z(pcaf2z);
  if (verbose) pca->display(1, 1);

  /* Set the error return code */

  error = 0;

  label_end: data1 = (double*) mem_free((char* ) data1);
  data2 = (double*) mem_free((char* ) data2);
  identity = (double*) mem_free((char* ) identity);
  pca2 = pca_free(pca2);
  if (iptr >= 0) (void) db_attribute_del_mult(db, iptr, nvar);
  return (error);
}

/****************************************************************************/
/*!
 **  Evaluate the PCA on irregular data
 **
 ** \return  Error return code
 **
 ** \param[in]  db          Db descriptor
 ** \param[in]  verbose     Verbose flag
 **
 ** \param[out] pca        Output PCA structure
 **
 *****************************************************************************/
int pca_compute(Db *db, int verbose, PCA *pca)
{
  int error, nvar;
  double *c0, *data1;
  static int flag_normalize = 1;

  /* Initializations */

  error = 0;
  if (db == nullptr) return (1);
  c0 = data1 = nullptr;

  /* Preliminary checks */

  nvar = pca->getNVar();
  if (nvar != db->getVariableNumber())
  {
    messerr(
        "Number of variables in the PCA (%d) and in the Db (%d) are inconsistent",
        nvar, db->getVariableNumber());
    goto label_end;
  }

  /* Core allocation */

  c0 = (double*) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (c0 == nullptr) goto label_end;
  data1 = (double*) mem_alloc(sizeof(double) * nvar, 0);
  if (data1 == nullptr) goto label_end;

  /* Calculate the PCA */

  /* Calculate the first PCA (centered and normalized) */

  if (st_pca_calculate(flag_normalize, db, data1, pca, verbose)) goto label_end;

  /* Set the error return code */

  error = 0;

  label_end: c0 = (double*) mem_free((char* ) c0);
  data1 = (double*) mem_free((char* ) data1);
  return (error);
}

/****************************************************************************/
/*!
 **  Transform the variables into factors using PCA
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  pca          PCA descriptor
 ** \param[in]  flag_norm_in 1 if the variable must be normalized
 ** \param[in]  verbose      1 for a verbose output
 **
 *****************************************************************************/
int pca_z2f(Db *db, PCA *pca, int flag_norm_in, int verbose)
{
  int iptr, nvar, error, ivar;
  double *data;
  VectorDouble mean, sigma;
  VectorInt cols;

  /* Preliminary checks */

  error = 1;
  nvar = db->getVariableNumber();
  data = nullptr;
  if (nvar <= 0)
  {
    messerr("The Transformation requires Z located variables");
    goto label_end;
  }
  if (nvar != pca->getNVar())
  {
    messerr(
        "The number of variables in PCA (%d) does not match the number of Z-variables in the Db (%d)",
        pca->getNVar(), nvar);
    goto label_end;
  }
  mean.resize(nvar);
  sigma.resize(nvar * nvar);

  /* Allocate new variables */

  iptr = db->addFieldsByConstant(nvar, TEST);
  if (iptr < 0) goto label_end;

  /* Core allocation */

  data = (double*) mem_alloc(sizeof(double) * nvar, 0);
  if (data == nullptr) goto label_end;

  /* Normalization (optional) */

  st_calculate_normalization(flag_norm_in, verbose, db, data, mean, sigma);

  /* Perform the normalization */

  st_pca_z2f(0, iptr, db, pca, data, mean, sigma);

  /* Optional printout */

  if (verbose)
  {
    cols.resize(nvar);
    for (ivar = 0; ivar < nvar; ivar++)
      cols[ivar] = iptr + ivar;
    db_stats_print(db, cols, VectorString(), 1, 1, "Statistics on Factors",
                   "Factor");
  }

  /* Set the error return code */

  error = 0;

  label_end: data = (double*) mem_free((char* ) data);
  return (error);
}

/****************************************************************************/
/*!
 **  Transform the factors into variables using PCA
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  pca          PCA descriptor
 ** \param[in]  flag_norm_out 1 if the output variable must be denormalized
 ** \param[in]  verbose      1 for a verbose output
 **
 *****************************************************************************/
int pca_f2z(Db *db, PCA *pca, int flag_norm_out, int verbose)
{
  int iptr, nvar, error;
  double *data;
  VectorInt cols;

  /* Preliminary checks */

  error = 1;
  nvar = db->getVariableNumber();
  data = nullptr;
  if (nvar <= 0)
  {
    messerr("The Transformation requires Z located variables");
    goto label_end;
  }
  if (nvar != pca->getNVar())
  {
    messerr(
        "The number of variables in PCA (%d) does not match the number of Z-variables in the Db (%d)",
        pca->getNVar(), nvar);
    goto label_end;
  }

  /* Allocate new variables */

  iptr = db->addFieldsByConstant(nvar, TEST);
  if (iptr < 0) goto label_end;

  /* Core allocation */

  data = (double*) mem_alloc(sizeof(double) * nvar, 0);
  if (data == nullptr) goto label_end;

  /* Rotate the factors into data in the PCA system */

  st_pca_f2z(flag_norm_out, iptr, db, pca, data);

  /* Optional printout */

  if (verbose)
  {
    cols.resize(nvar);
    for (int ivar = 0; ivar < nvar; ivar++)
      cols[ivar] = iptr + ivar;
    db_stats_print(db, cols, VectorString(), 1, 1, "Statistics on Variables",
                   "Variable");
  }

  /* Set the error return code */

  error = 0;

  label_end: data = (double*) mem_free((char* ) data);
  return (error);
}

/****************************************************************************/
/*!
 **  Calculate the geometry for a given direction
 **
 ** \return  Error return code
 **
 ** \param[in]  db     Db description
 ** \param[in]  vario  Vario structure
 **
 ** \param[out] vorder Vario_Order structure
 ** \param[out] npair  Number of pairs
 **
 *****************************************************************************/
int geometry_compute(Db *db, Vario *vario, Vario_Order *vorder, int *npair)
{
  double psmin, ps, dist, maxdist;
  int iiech, iech, jjech, jech, nech, ipas, error, idir, ideb;
  VectorInt rindex;

  /* Initializations */

  if (db == nullptr) return (1);
  if (vario == nullptr) return (1);
  const VarioParam &varioparam = vario->getVarioParam();
  error = 1;

  /* Preliminary checks */

  if (db->getNDim() != vario->getDimensionNumber() || db->getVariableNumber()
      != vario->getVariableNumber())
  {
    messerr("Inconsistent parameters:");
    messerr("Data Base: NDIM=%d NVAR=%d", db->getNDim(),
            db->getVariableNumber());
    messerr("Variogram: NDIM=%d NVAR=%d", vario->getDimensionNumber(),
            vario->getVariableNumber());
    goto label_end;
  }
  if (st_get_generalized_variogram_order(vario) > 0)
  {
    messerr("This calculation does not allow generalized variogram definition");
    goto label_end;
  }

  /* Sort the data */

  rindex = db->getSortArray();

  /* Loop on the directions */

  for (idir = 0; idir < vario->getDirectionNumber(); idir++)
  {
    const DirParam &dirparam = vario->getDirParam(idir);
    ps = 0.;
    psmin = _variogram_convert_angular_tolerance(dirparam.getTolAngle());
    nech = db->getSampleNumber();
    maxdist = vario->getMaximumDistance(idir);

    /* Loop on the first point */

    for (iiech = 0; iiech < nech - 1; iiech++)
    {
      iech = rindex[iiech];
      if (!db->isActive(iech)) continue;
      if (FFFF(db->getWeight(iech))) continue;

      ideb = (st_date_is_used(&varioparam, db, db)) ? 0 :
                                                      iiech + 1;
      for (jjech = ideb; jjech < nech; jjech++)
      {
        jech = rindex[jjech];
        if (variogram_maximum_dist1D_reached(db, iech, jech, maxdist)) break;
        if (!db->isActive(jech)) continue;
        if (FFFF(db->getWeight(jech))) continue;

        /* Check if the pair must be kept (Code criterion) */

        if (code_comparable(db, db, iech, jech, dirparam.getOptionCode(),
                            (int) dirparam.getTolCode())) continue;

        /* Check if the pair must be kept (Date criterion) */

        if (st_date_comparable(&varioparam, db, db, iech, jech,
                               dirparam.getIdate())) continue;

        /* Check if the pair must be kept */

        dist = distance_intra(db, iech, jech, NULL);
        if (variogram_reject_pair(db, iech, jech, dist, psmin,
                                  dirparam.getBench(), dirparam.getCylRad(),
                                  dirparam.getCodir(), &ps)) continue;

        /* Get the rank of the lag */

        ipas = variogram_get_lag(vario, idir, ps, psmin, &dist);
        if (IFFFF(ipas)) continue;

        /* Case of internal storage */

        vario_order_add(vorder, iech, jech, NULL, NULL, ipas, idir, dist);
      }
    }
  }

  /* Sort the geometry */

  vorder = vario_order_final(vorder, npair);

  /* Set the error return code */

  error = 0;

  label_end: return (error);
}

/****************************************************************************/
/*!
 **  Linear interpolation
 **
 ** \return  Interpolated value
 **
 ** \param[in]  n      Number of discretization steps
 ** \param[in]  x      Discretized X (sorted increasingly)
 ** \param[in]  y      Discretized Y
 ** \param[in]  x0     Origin
 **
 *****************************************************************************/
static double st_linear_interpolate(int n, double *x, double *y, double x0)
{
  int i;

  if (x0 < x[0]) return (y[0]);
  if (x0 > x[n - 1]) return (y[n - 1]);
  for (i = 1; i < n; i++)
  {
    if (x0 < x[i - 1]) continue;
    if (x0 > x[i]) continue;
    return (y[i - 1] + (y[i] - y[i - 1]) * (x0 - x[i - 1]) / (x[i] - x[i - 1]));
  }
  return (TEST);
}

/****************************************************************************/
/*!
 **  Calculate the experimental variogram of the completed variable starting
 **  from the experimental variogram of the truncated variable
 **
 ** \param[in,out] vario  Vario structure
 ** \param[in]  nh     Number of Hermite polynomials
 ** \param[in]  ycut   Truncation (lowest) value
 **
 *****************************************************************************/
void variogram_trans_cut(Vario *vario, int nh, double ycut)
{
  double variance, sum, cyp, cyy;
  int ih, idisc, ndisc, idir, ipas;
  static double disc = 0.01;

  /* Initializations */

  ndisc = (int) (2. / disc + 1.);

  /* Core allocation */

  VectorDouble fncoeff(nh);
  VectorDouble ro(ndisc);
  VectorDouble covyp(ndisc);
  for (idisc = 0; idisc < ndisc; idisc++)
    ro[idisc] = disc * idisc - 1.;

  /* Calculate the first normalized hermite polynomials for ycut */

  VectorDouble hermite = hermiteFunction(ycut, nh);

  /* Variance */

  variance = 0.;
  for (ih = 1; ih < nh; ih++)
    variance += fncoeff[ih] * fncoeff[ih];
  for (idisc = 0; idisc < ndisc; idisc++)
  {
    sum = 0.;
    for (ih = 1; ih < nh; ih++)
      sum += fncoeff[ih] * fncoeff[ih] * pow(ro[idisc], ih);
    covyp[idisc] = sum;
  }

  /* Loop on the directions */

  for (idir = 0; idir < vario->getDirectionNumber(); idir++)
  {

    /* Loop on the lags */

    for (ipas = 0; ipas < vario->getLagNumber(idir); ipas++)
    {
      cyp = variance - vario->getGg(idir, 0, 0, ipas);
      cyy = st_linear_interpolate(ndisc, covyp.data(), ro.data(), cyp);
      vario->setGg(idir, 0, 0, ipas, MAX(0, 1. - cyy));
    }
  }

  /* Set the variance */

  vario->setVar(0, 0, 1.);
}

/****************************************************************************/
/*!
 **  Determine the samples used for a variogram in multilayers framework
 **
 ** \return  Error return code
 **
 ** \param[in]  db     Db description
 ** \param[in]  seltab Number of sample definition (0, 1 or 2)
 ** \param[in]  vario  Vario structure
 **
 ** \param[out]  vorder Vario_Order struct
 ure
 **
 *****************************************************************************/
int variogram_mlayers(Db *db, int *seltab, Vario *vario, Vario_Order *vorder)
{
  int iiech, iech, jjech, jech, nech, ipas, npair, idir;
  double psmin, ps, dist;

  /* Initializations */

  if (db == nullptr) return (1);
  if (vario == nullptr) return (1);

  /* Loop on the directions */

  for (idir = 0; idir < vario->getDirectionNumber(); idir++)
  {
    const DirParam &dirparam = vario->getDirParam(idir);
    ps = 0.;
    psmin = _variogram_convert_angular_tolerance(dirparam.getTolAngle());
    nech = db->getSampleNumber();

    /* Loop on the first point */

    for (iech = iiech = 0; iech < nech; iech++)
    {
      if (!db->isActive(iech)) continue;
      if (seltab[iech] == 0) continue;
      for (int ifois = 0; ifois < seltab[iech]; ifois++, iiech++)
      {

        for (jech = jjech = 0; jech < nech; jech++)
        {
          if (!db->isActive(jech)) continue;
          if (seltab[jech] == 0) continue;
          for (int jfois = 0; jfois < seltab[jech]; jfois++, jjech++)
          {

            /* Check if the pair must be kept (Code criterion) */

            if (code_comparable(db, db, iech, jech, dirparam.getOptionCode(),
                                (int) dirparam.getTolCode())) continue;

            /* Check if the pair must be kept */

            dist = distance_intra(db, iech, jech, NULL);
            if (variogram_reject_pair(db, iech, jech, ABS(dist), psmin,
                                      dirparam.getBench(), dirparam.getCylRad(),
                                      dirparam.getCodir(), &ps)) continue;

            /* Get the rank of the lag */

            ipas = variogram_get_lag(vario, idir, ps, psmin, &dist);
            if (IFFFF(ipas)) continue;

            /* Internal storage */

            vario_order_add(vorder, iiech, jjech, &iech, &jech, ipas, idir,
                            ABS(dist));
          }
        }
      }
    }
  }
  vorder = vario_order_final(vorder, &npair);
  return (0);
}

/*****************************************************************************/
/*!
 **  Calculate the experimental variogram of the Raw starting from the Model
 **  of the Gaussian variable
 **
 ** \return  Error return code
 **
 ** \param[in,out] vario    Experimental variogram
 ** \param[in]  anam        Point anamorphosis
 ** \param[in]  model       Model of the Punctual Gaussian
 **
 ** \remark  At entrance, the input variogram only serves in providing
 ** \remark  the calculation parameters
 **
 *****************************************************************************/
int variogram_y2z(Vario *vario, Anam *anam, Model *model)
{
  int error, idir, ndim;
  double chh, varz, cov_value;
  VectorDouble d1;
  CovCalcMode mode;

  /* Preliminary checks */

  error = 1;
  if (vario == nullptr) return (error);
  if (anam == (Anam*) NULL) return (error);
  if (model == nullptr) return (error);
  if (anam->getType() != EAnam::HERMITIAN)
  {
    messerr("This function is restricted to Gaussian Anamorphosis");
    return (error);
  }
  AnamHermite *anam_hermite = dynamic_cast<AnamHermite*>(anam);
  if (anam_hermite->getRCoef() != 1.)
  {
    messerr("This function is restricted to Punctual Anamoprhosis");
    return (error);
  }
  if (vario == nullptr) return (error);
  if (vario->getVariableNumber() != 1)
  {
    messerr("This function is restricted to Monovariate Variogram");
    return (error);
  }
  if (model->getVariableNumber() != 1)
  {
    messerr("This function requires a Monovariate Model");
    return (error);
  }
  if (model->getDimensionNumber() != vario->getDimensionNumber())
  {
    messerr("Variogram and Model should share the same Space Dimension");
    return (error);
  }

  /* Initializations */

  ndim = vario->getDimensionNumber();

  /* Core allocation */

  d1.resize(ndim);

  /* Calculate the theoretical variance of Z */

  varz = anam_hermite->calculateVarianceFromPsi(1.);

  /* Loop on the directions of the variogram */

  for (idir = 0; idir < vario->getDirectionNumber(); idir++)
  {
    /* Loop on the lags */

    for (int ipas = 0; ipas < vario->getLagNumber(idir); ipas++)
    {
      for (int idim = 0; idim < ndim; idim++)
        d1[idim] = (ipas + 1) * vario->getDPas(idir)
                   * vario->getCodir(idir, idim);

      model_calcul_cov(model, mode, 1, 1., d1, &chh);
      if (chh < 0.)
      {
        messerr("Gaussian covariance is negative in direction %d for lag %d",
                idir + 1, ipas + 1);
        messerr("Calculation is impossible");
        return (error);
      }

      cov_value = anam_hermite->calculateVarianceFromPsi(chh);
      vario->setGg(idir, 0, 0, ipas, varz - cov_value);
      vario->setHh(idir, 0, 0, ipas, (ipas + 1) * vario->getDPas(idir));
      vario->setSw(idir, 0, 0, ipas, 1.);
    }
  }

  error = 0;
  return (error);
}

/****************************************************************************/
/*!
 **  Evaluate the experimental conditional expectation
 **
 ** \param[in]  db1           Db descriptor (for target variable)
 ** \param[in]  db2           Db descriptor (for auxiliary variables)
 ** \param[in]  icol1         Rank of the target variable
 ** \param[in]  icol2         Rank of the explanatory variable
 ** \param[in]  mini          Minimum value for the explanaroty variable
 ** \param[in]  maxi          Maximum value for the explanaroty variable
 ** \param[in]  nclass        Number of classes
 ** \param[in]  verbose       Verbose flag
 **
 ** \param[out] ncond         Array of number of samples per class
 ** \param[out] xcond         Array of conditional expectation along X
 ** \param[out] ycond         Array of conditional expectation along Y
 **
 *****************************************************************************/
void condexp(Db *db1,
             Db *db2,
             int icol1,
             int icol2,
             double mini,
             double maxi,
             int nclass,
             int verbose,
             int *ncond,
             double *xcond,
             double *ycond)
{
  int rank;
  double val1, val2;

  /* Initializations */

  for (int iclass = 0; iclass < nclass; iclass++)
  {
    ncond[iclass] = 0;
    xcond[iclass] = 0.;
    ycond[iclass] = 0.;
  }

  /* Loop on the samples */

  for (int iech = 0; iech < db1->getSampleNumber(); iech++)
  {
    if (!db1->isActive(iech)) continue;
    val1 = db1->getArray(iech, icol1);
    if (FFFF(val1)) continue;
    val2 = db2->getArray(iech, icol2);
    if (FFFF(val2)) continue;
    if (val2 < mini || val2 > maxi) continue;

    rank = int((nclass - 1.) * (val2 - mini) / (maxi - mini));

    ncond[rank]++;
    xcond[rank] += val2;
    ycond[rank] += val1;
  }

  /* Normation */

  for (int iclass = 0; iclass < nclass; iclass++)
  {
    if (ncond[iclass] <= 0)
    {
      xcond[iclass] = TEST;
      ycond[iclass] = TEST;
    }
    else
    {
      xcond[iclass] /= ncond[iclass];
      ycond[iclass] /= ncond[iclass];
    }
  }

  /* Optional printout */

  if (verbose)
  {
    message("Experimental Conditional Expectation\n");
    for (int iclass = 0; iclass < nclass; iclass++)
    {
      if (ncond[iclass] > 0)
        message("Class %2d : V1=%lf V2=%lf\n", iclass + 1, xcond[iclass],
                ycond[iclass]);
    }
  }
}

/****************************************************************************/
/*!
 **  Evaluate the experimental variogram
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  vario        Vario structure
 ** \param[in]  flag_grid    1 for calculation on a grid
 ** \param[in]  flag_gen     1 for calculation of generalized variogram
 ** \param[in]  flag_sample  calculate the variogram per sample
 ** \param[in]  verr_mode    Mode of variogram correction (1, 2 or 3)
 ** \param[in]  model        Model structure (triggers the KU option)
 ** \param[in]  verbose      Verbose flag
 **
 ** \note The number of iterations in KU is given can be updated using the
 ** \note keypair technique with name "KU_Niter".
 **
 *****************************************************************************/
int _variogram_compute(Db *db,
                       Vario *vario,
                       int flag_grid,
                       int flag_gen,
                       int flag_sample,
                       int verr_mode,
                       Model *model,
                       int verbose)
{
  int error;

  if (flag_grid)
  {
    if (flag_gen)
      error = st_variogen_grid_calcul(db, vario);
    else
      error = st_variogrid_calcul(db, vario);
  }
  else
  {
    if (flag_gen)
      error = st_variogen_line_calcul(db, vario);
    else
      error = st_variogram_general(db, vario, model, flag_sample, verr_mode,
                                   verbose);
  }
  return (error);
}

/****************************************************************************
 **
 ** FUNCTION: st_identify_calcul_type
 **
 ** PURPOSE:  Identify the type of variogram calculation
 **
 ** IN_ARGS:  calcul_name  : Type of the variogram
 **
 ** REMARKS:  In case the calculation type is not identified,
 ** REMARKS:  the routine returns ECalcVario::UNDEFINED
 ** REMARKS:  The error message is produced internally
 **
 *****************************************************************************/
ECalcVario vario_identify_calcul_type(const String &calcul_name)

{
  ECalcVario calcul_type;

  if (!strcmp(calcul_name.c_str(), "vg"))
    calcul_type = ECalcVario::VARIOGRAM;
  else if (!strcmp(calcul_name.c_str(), "cov"))
    calcul_type = ECalcVario::COVARIANCE;
  else if (!strcmp(calcul_name.c_str(), "covnc"))
    calcul_type = ECalcVario::COVARIANCE_NC;
  else if (!strcmp(calcul_name.c_str(), "covg"))
    calcul_type = ECalcVario::COVARIOGRAM;
  else if (!strcmp(calcul_name.c_str(), "mado"))
    calcul_type = ECalcVario::MADOGRAM;
  else if (!strcmp(calcul_name.c_str(), "rodo"))
    calcul_type = ECalcVario::RODOGRAM;
  else if (!strcmp(calcul_name.c_str(), "poisson"))
    calcul_type = ECalcVario::POISSON;
  else if (!strcmp(calcul_name.c_str(), "general1"))
    calcul_type = ECalcVario::GENERAL1;
  else if (!strcmp(calcul_name.c_str(), "general2"))
    calcul_type = ECalcVario::GENERAL2;
  else if (!strcmp(calcul_name.c_str(), "general3"))
    calcul_type = ECalcVario::GENERAL3;
  else if (!strcmp(calcul_name.c_str(), "order4"))
    calcul_type = ECalcVario::ORDER4;
  else if (!strcmp(calcul_name.c_str(), "trans1"))
    calcul_type = ECalcVario::TRANS1;
  else if (!strcmp(calcul_name.c_str(), "trans2"))
    calcul_type = ECalcVario::TRANS2;
  else if (!strcmp(calcul_name.c_str(), "binormal"))
    calcul_type = ECalcVario::BINORMAL;
  else
  {
    messerr("Invalid variogram calculation name : %s", calcul_name.c_str());
    messerr("The only valid names are:");
    messerr("vg       : Variogram");
    messerr("cov      : Covariance");
    messerr("covnc    : Non-centered ergodic covariance");
    messerr("covg     : Covariogram");
    messerr("mado     : Madogram");
    messerr("rodo     : Rodogram");
    messerr("poisson  : Poisson");
    messerr("general1 : Generalized variogram of order 1");
    messerr("general2 : Generalized variogram of order 2");
    messerr("general3 : Generalized variogram of order 3");
    messerr("order4   : Variogram of order 4");
    messerr("trans1   : Cross-to-Simple Variogram G12/G1");
    messerr("trans2   : Cross-to-Simple Variogram G12/G1");
    messerr("binormal : Cross-to-Simple Variogram G12/sqrt(G1*G2)");

    calcul_type = ECalcVario::UNDEFINED;
  }
  return (calcul_type);
}

/****************************************************************************/
/*!
 **  Evaluate the experimental variogram cloud
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  varioparam   VarioParam structure
 ** \param[in]  lagmax       Maximum distance
 ** \param[in]  varmax       Maximum Variance value
 ** \param[in]  lagnb        Number of discretization steps along distance axis
 ** \param[in]  varnb        Number of discretization steps along variance axis
 ** \param[in]  namconv      Naming convention
 **
 *****************************************************************************/
Db* db_variogram_cloud(Db *db,
                       const VarioParam *varioparam,
                       double lagmax,
                       double varmax,
                       int lagnb,
                       int varnb,
                       const NamingConvention& namconv)
{
  if (FFFF(lagmax)) lagmax = db->getFieldSize();
  if (FFFF(varmax)) (void) variogram_cloud_dim(db, varioparam, &varmax);

  // Create a grid as a support for the variogram cloud calculations

  VectorInt nx(2);
  nx[0] = lagnb;
  nx[1] = varnb;
  VectorDouble dx(2);
  dx[0] = lagmax / (double) lagnb;
  dx[1] = varmax / (double) varnb;
  VectorDouble x0(2);
  x0[0] = 0.;
  x0[1] = 0.;
  Db *dbgrid = new Db(nx, dx, x0);

  // Calling the variogram cloud calculation function

  int error = variogram_cloud(db, varioparam, dbgrid, namconv);

  // In case of error, free the newly created structure

  if (error) dbgrid = db_delete(dbgrid);

  return dbgrid;
}

/****************************************************************************/
/*!
 **  Calculate the variogram map
 **
 ** \return  Error return code
 **
 ** \param[in]  db          Db containing the data
 ** \param[in]  dbmap       VMAP grid structure
 ** \param[in]  calcul_type Type of calculation (ECalcVario)
 ** \param[in]  radius      Dilation radius (mooth resulting maps) only on points
 ** \param[in]  flag_FFT    Use FFT method (only valid on grid)
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int vmap_compute(Db *db,
                 Db *dbmap,
                 const ECalcVario &calcul_type,
                 int radius,
                 bool flag_FFT,
                 const NamingConvention& namconv)
{
  int error = 0;

  // Calculating the variogram map in different ways

  if (db->isGrid())
  {
    // Case where Data are on a regular grid

    if (flag_FFT)
      error = st_vmap_grid_fft(db, dbmap, calcul_type, namconv);
    else
      error = st_vmap_grid(db, dbmap, calcul_type, namconv);
  }
  else
  {
    // case where Data are on a set of points

    error = st_vmap_general(db, dbmap, calcul_type, radius, namconv);
  }

  return error;
}

/****************************************************************************/
/*!
 **  Calculate the variogram map
 **
 ** \return  Error return code
 **
 ** \param[in]  db          Db containing the data
 ** \param[in]  calcul_type Type of calculation (ECalcVario)
 ** \param[in]  nxx         (Half-) Number of nodes for the VMAP grid along X
 ** \param[in]  nyy         (Half-) Number of nodes for the VMAP grid along Y
 ** \param[in]  dxx         Mesh of the VMAP grid along X
 ** \param[in]  dyy         Mesh of the VMAP grid along Y
 ** \param[in]  radius      Dilation radius (mooth resulting maps) only on points
 ** \param[in]  flag_FFT    Use FFT method (only valid on grid)
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
Db* db_vmap_compute(Db *db,
                    const ECalcVario &calcul_type,
                    int nxx,
                    int nyy,
                    double dxx,
                    double dyy,
                    int radius,
                    bool flag_FFT,
                    const NamingConvention& namconv)
{
  int error = 0;

  // Creating the output Variogram Map grid

  VectorInt nx(2);
  nx[0] = 2 * nxx + 1;
  nx[1] = 2 * nyy + 1;
  VectorDouble dx(2);
  if (db->isGrid())
  {
    dx[0] = db->getDX(0);
    dx[1] = db->getDX(1);
  }
  else
  {
    dx[0] = (!FFFF(dxx)) ? dxx :
                           db->getExtension(0) / (double) nxx;
    dx[1] = (!FFFF(dyy)) ? dyy :
                           db->getExtension(1) / (double) nyy;
  }
  VectorDouble x0(2);
  x0[0] = -nxx * dx[0];
  x0[1] = -nyy * dx[1];

  Db *dbmap = new Db(nx, dx, x0);

  // Calculating the variogram map in different ways

  error = vmap_compute(db, dbmap, calcul_type, radius, flag_FFT, namconv);

  // In case of error, free the newly created VMAP structure

  if (error) dbmap = db_delete(dbmap);

  return dbmap;
}
