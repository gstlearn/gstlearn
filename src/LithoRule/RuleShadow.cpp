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
#include "geoslib_old_f.h"
#include "geoslib_enum.h"

#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"
#include "LithoRule/RuleShadow.hpp"
#include "LithoRule/Rule.hpp"
#include "LithoRule/Node.hpp"
#include "LithoRule/PropDef.hpp"
#include "Model/Model.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"

#include <sstream>
#include <math.h>

RuleShadow::RuleShadow()
    : Rule(),
      _shDsup(0.),
      _shDown(0.),
      _slope(0.),
      _shift(),
      _dMax(TEST),
      _tgte(TEST),
      _incr(TEST),
      _xyz(),
      _ind1(),
      _ind2()
{
  setModeRule(ERule::SHADOW);
}

/**
 * Defining the Lithotype Rule in the case of Shadow
 * @param slope   Slope definition
 * @param sh_dsup Maximum threshold
 * @param sh_down Minimum threshold
 * @param shift   Shift orientation
 */
RuleShadow::RuleShadow(double slope,
                       double sh_dsup,
                       double sh_down,
                       const VectorDouble &shift)
    :
    Rule(),
    _shDsup(sh_dsup),
    _shDown(sh_down),
    _slope(slope),
    _shift(shift),
    _dMax(TEST),
    _tgte(TEST),
    _incr(TEST),
    _xyz(),
    _ind1(),
    _ind2()
{
  setModeRule(ERule::SHADOW);
  VectorString nodnames = { "S", "T", "F1", "F2", "F3" };
  setMainNodeFromNodNames(nodnames);
  _normalizeShift();
}

RuleShadow::RuleShadow(const RuleShadow &m)
    :
    Rule(m),
    _shDsup(m._shDsup),
    _shDown(m._shDown),
    _slope(m._slope),
    _shift(m._shift),
    _dMax(m._dMax),
    _tgte(m._tgte),
    _incr(m._incr)
{
}

RuleShadow& RuleShadow::operator=(const RuleShadow &m)
{
  if (this != &m)
  {
    Rule::operator =(m);
    _shDsup = m._shDsup;
    _shDown = m._shDown;
    _slope = m._slope;
    _shift = m._shift;
    _dMax = m._dMax;
    _tgte = m._tgte;
    _incr = m._incr;
  }
  return *this;
}

RuleShadow::~RuleShadow()
{
}

bool RuleShadow::_deserialize(std::istream& is, bool /*verbose*/)
{
  bool ret = true;
  _shift.resize(3);

  ret = ret && Rule::_deserialize(is);

  ret = ret && _recordRead<double>(is, "Slope for Shadow Rule", _slope);
  ret = ret && _recordRead<double>(is, "Lower Threshold for Shadow Rule", _shDown);
  ret = ret && _recordRead<double>(is, "Upper Threshold for Shadow Rule", _shDsup);
  ret = ret && _recordRead<double>(is, "Shift along first direction", _shift[0]);
  ret = ret && _recordRead<double>(is, "Shift along second direction", _shift[1]);
  ret = ret && _recordRead<double>(is, "Shift along third direction", _shift[2]);
  return ret;
}

bool RuleShadow::_serialize(std::ostream& os, bool /*verbose*/) const
{
  double slope = (FFFF(_slope)) ? 0. : _slope;
  double shdown = (FFFF(_shDown)) ? 0. : _shDown;
  double shdsup = (FFFF(_shDsup)) ? 0. : _shDsup;
  VectorDouble shiftloc = _shift;
  shiftloc.resize(3);

  bool ret = true;

  ret = ret && Rule::_serialize(os);

  ret = ret && _recordWrite<double>(os, "", slope);
  ret = ret && _recordWrite<double>(os, "", shdown);
  ret = ret && _recordWrite<double>(os, "Parameters for Shadow option", shdsup);
  ret = ret && _recordWrite<double>(os, "", shiftloc[0]);
  ret = ret && _recordWrite<double>(os, "", shiftloc[1]);
  ret = ret && _recordWrite<double>(os, "Parameters for Shift option", shiftloc[2]);
  return ret;
}

String RuleShadow::displaySpecific() const
{
  std::stringstream sstr;
  sstr << toTitle(2, "Shadow Option");
  sstr << toVector("Normalized Translation Vector = ", _shift);
  sstr << "Slope for shadow                  = " << _slope << "(degrees)"
       << std::endl;
  sstr << "Upwards shift for the threshold   = " << _shDsup << std::endl;
  sstr << "Downwards shift for the threshold = " << _shDown << std::endl;
  sstr << std::endl;
  sstr << "Note for non-stationary case:" << std::endl;
  sstr << "- P1 gives the proportion of Island" << std::endl;
  sstr << "- P2 gives the value of Upwards shift" << std::endl;
  sstr << "- P3 gives the value of Downwards shift" << std::endl;
  sstr << "(With the 'Shadow' option, only the first GRF is used)" << std::endl;

  return sstr.str();
}

/****************************************************************************/
/*!
 **  Define the particularities of the PGS model (for Shadow)
 **
 ** \return  Error return code
 **
 ** \param[in]  db              Db structure
 ** \param[in]  dbprop          Db structure used for proportions
 ** \param[in]  model           Model structure (only used for shift option)
 ** \param[in]  flag_grid_check 1 if grid is compulsory; 0 otherwise
 **                             (only for SHIFT)
 ** \param[in]  flag_stat       1 for stationary; 0 otherwise
 **
 *****************************************************************************/
int RuleShadow::particularities(Db *db,
                                const Db *dbprop,
                                Model *model,
                                int /*flag_grid_check*/,
                                int flag_stat) const
{
  double sh_dsup_max, sh_down_max;
  int ndim = (model != nullptr) ? model->getDimensionNumber() : 0;

  _incr = 1.e30;
  for (int idim = 0; idim < ndim; idim++)
    if (_shift[idim] != 0) _incr = MIN(_incr, db->getUnit(idim));

  /* Calculate the maximum distance */

  _tgte = tan(ut_deg2rad(_slope));
  _st_shadow_max(dbprop, flag_stat, &sh_dsup_max, &sh_down_max);
  _dMax = (_tgte > 0) ? (sh_dsup_max + sh_down_max) / _tgte : 0.;

  return (0);
}

void RuleShadow::_st_shadow_max(const Db *dbprop,
                                int flag_stat,
                                double *sh_dsup_max,
                                double *sh_down_max) const
{
  int iech;
  double val2, val3;

  if (flag_stat || dbprop == nullptr)
  {
    /* Stationary case */

    *sh_dsup_max = _shDsup;
    *sh_down_max = _shDown;
  }
  else
  {
    *sh_dsup_max = *sh_down_max = 0.;
    for (iech = 0; iech < dbprop->getSampleNumber(); iech++)
    {
      val2 = dbprop->getLocVariable(ELoc::P,iech, 1);
      if (val2 > (*sh_dsup_max)) (*sh_dsup_max) = val2;
      val3 = dbprop->getLocVariable(ELoc::P,iech, 2);
      if (val3 > (*sh_down_max)) (*sh_down_max) = val3;
    }
  }
}

/****************************************************************************/
/*!
 **  Evaluation of the value from a grid by inverse square distance
 **  interpolation from the 4 surrounding grid nodes
 **
 ** \return  Interpolated value or FFFF if out of grid
 **
 ** \param[in]  dbgrid Db structure
 ** \param[in]  isimu  Rank of the simulation
 ** \param[in]  icase  Rank of the Simulation storage
 ** \param[in]  nbsimu Number of simulations
 **
 ** \param[out]  xyz0  Working array
 **
 *****************************************************************************/
double RuleShadow::_st_grid_eval(DbGrid *dbgrid,
                                 int isimu,
                                 int icase,
                                 int nbsimu,
                                 VectorDouble &xyz0) const
{
  double top = 0.;
  double bot = 0.;
  int ndim = dbgrid->getNDim();

  /* First point */
  int iech = dbgrid->indiceToRank(_ind2);
  double z = dbgrid->getSimvar(ELoc::SIMU, iech, isimu, 0, icase, nbsimu, 1);
  if (!FFFF(z))
  {
    double d2 = 0.;
    dbgrid->indicesToCoordinateInPlace(_ind2, xyz0);
    for (int idim = 0; idim < ndim; idim++)
    {
      double delta = _xyz[idim] - xyz0[idim];
      d2 += delta * delta;
    }
    if (d2 <= 0.) return (z);
    top += z / d2;
    bot += 1. / d2;
  }

  /* Second point */
  _ind2[0] += 1;
  iech = dbgrid->indiceToRank(_ind2);
  z = dbgrid->getSimvar(ELoc::SIMU, iech, isimu, 0, icase, nbsimu, 1);
  if (!FFFF(z))
  {
    double d2 = 0.;
    dbgrid->indicesToCoordinateInPlace(_ind2, xyz0);
    for (int idim = 0; idim < ndim; idim++)
    {
      double delta = _xyz[idim] - xyz0[idim];
      d2 += delta * delta;
    }
    if (d2 <= 0.) return (z);
    top += z / d2;
    bot += 1. / d2;
  }

  /* Third point */
  _ind2[1] += 1;
  iech = dbgrid->indiceToRank(_ind2);
  z = dbgrid->getSimvar(ELoc::SIMU, iech, isimu, 0, icase, nbsimu, 1);
  if (!FFFF(z))
  {
    double d2 = 0.;
    dbgrid->indicesToCoordinateInPlace(_ind2, xyz0);
    for (int idim = 0; idim < ndim; idim++)
    {
      double delta = _xyz[idim] - xyz0[idim];
      d2 += delta * delta;
    }
    if (d2 <= 0.) return (z);
    top += z / d2;
    bot += 1. / d2;
  }

  /* Fourth point */
  _ind2[0] -= 1;
  iech = dbgrid->indiceToRank(_ind2);
  z = dbgrid->getSimvar(ELoc::SIMU, iech, isimu, 0, icase, nbsimu, 1);
  if (!FFFF(z))
  {
    double d2 = 0.;
    dbgrid->indicesToCoordinateInPlace(_ind2, xyz0);
    for (int idim = 0; idim < ndim; idim++)
    {
      double delta = _xyz[idim] - xyz0[idim];
      d2 += delta * delta;
    }
    if (d2 <= 0.) return (z);
    top += z / d2;
    bot += 1. / d2;
  }

  /* Final interpolation */
  _ind2[1] -= 1;
  z = (bot != 0) ? top / bot :
                   TEST;
  return (z);
}

/****************************************************************************/
/*!
 **  Combine the underlying GRF into a facies value at data points (Shadow)
 **
 ** \return  Error return code
 **
 ** \param[in]  propdef    Props structure
 ** \param[in]  dbin       Db input structure
 ** \param[in]  dbout      Db output structure
 ** \param[in]  flag_used  1 if the gaussian is used; 0 otherwise
 ** \param[in]  ipgs       Rank of the PGS
 ** \param[in]  isimu      Rank of the simulation
 ** \param[in]  nbsimu     Number of simulations
 **
 ** \remark Attributes ELoc::GAUSFAC are mandatory
 ** \remark Attributes ELoc::FACIES are mandatory
 **
 *****************************************************************************/
int RuleShadow::gaus2facData(PropDef *propdef,
                             Db *dbin,
                             Db* /*dbout*/,
                             int *flag_used,
                             int ipgs,
                             int isimu,
                             int nbsimu)
{
  double y[2], facies, t1min, t1max, t2min, t2max, sh_dsup, sh_down;

  /* Initializations */

  check_mandatory_attribute("rule_gaus2fac_data_shadow", dbin, ELoc::GAUSFAC);

  /* Processing the translation */

  for (int iech = 0; iech < dbin->getSampleNumber(); iech++)
  {
    if (!dbin->isActive(iech)) continue;

    /* Initializations */

    facies = TEST;
    for (int igrf = 0; igrf < 2; igrf++)
      y[igrf] = TEST;

    if (rule_thresh_define_shadow(propdef, dbin, this, ITEST, iech, isimu,
                                  nbsimu, &t1min, &t1max, &t2min, &t2max,
                                  &sh_dsup, &sh_down)) return 1;

    for (int igrf = 0; igrf < 2; igrf++)
    {
      int icase = get_rank_from_propdef(propdef, ipgs, igrf);
      y[igrf] =
          (flag_used[igrf]) ? dbin->getSimvar(ELoc::GAUSFAC, iech, isimu, 0,
                                              icase, nbsimu, 1) :
                              0.;
    }
    facies = getFaciesFromGaussian(y[0], y[1]);

    /* Combine the underlying GRFs to derive Facies*/

    dbin->setSimvar(ELoc::FACIES, iech, isimu, 0, ipgs, nbsimu, 1, facies);
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Combine the underlying GRF into a facies value (Shadow case)
 **
 ** \return  Error return code
 **
 ** \param[in]  propdef    Props structure
 ** \param[in]  dbout      Db output structure
 ** \param[in]  flag_used  1 if the gaussian is used; 0 otherwise
 ** \param[in]  ipgs       Rank of the PGS
 ** \param[in]  isimu      Rank of the simulation
 ** \param[in]  nbsimu     Number of simulations
 **
 ** \remark Attributes ELoc::FACIES and ELoc::SIMU are mandatory
 **
 *****************************************************************************/
int RuleShadow::gaus2facResult(PropDef *propdef,
                               Db *dbout,
                               int* /*flag_used*/,
                               int ipgs,
                               int isimu,
                               int nbsimu) const
{
  int ndim, iech, jech, error, idim, nstep, istep, flag, flag_shadow, igrf, icase;
  double *del, y[2], facies, dinc, dy, ys, yc_dsup, yc_down;
  double t1min, t1max, t2min, t2max, s1min, s1max, s2min, s2max, sh_dsup, sh_down, seuil;

  /* Initializations */

  check_mandatory_attribute("rule_gaus2fac_result_shadow", dbout, ELoc::FACIES);
  check_mandatory_attribute("rule_gaus2fac_result_shadow", dbout, ELoc::SIMU);
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(dbout);
  if (dbgrid == nullptr) return 1;

  error = 1;
  del = nullptr;
  nstep = 0;
  ndim = dbgrid->getNDim();
  icase = get_rank_from_propdef(propdef, ipgs, 0);
  _xyz.resize(ndim);
  _ind1.resize(ndim);
  _ind2.resize(ndim);

  /* Initializations */

  del = db_vector_alloc(dbgrid);
  if (del == nullptr) goto label_end;
  dinc = getIncr();
  nstep = (int) floor(getDMax() / dinc);
  dy = dinc * getTgte();
  for (idim = 0; idim < ndim; idim++)
    del[idim] = dinc * _shift[idim];

  /* Processing the translation */

  for (iech = 0; iech < dbgrid->getSampleNumber(); iech++)
  {
    if (!dbgrid->isActive(iech)) continue;

    /* Initializations */

    facies = TEST;
    for (igrf = 0; igrf < 2; igrf++)
      y[igrf] = TEST;

    y[0] = dbgrid->getSimvar(ELoc::SIMU, iech, isimu, 0, icase, nbsimu, 1);
    if (FFFF(y[0])) break;
    if (rule_thresh_define_shadow(propdef, dbgrid, this, SHADOW_WATER, iech,
                                  isimu, nbsimu, &t1min, &t1max, &t2min, &t2max,
                                  &yc_dsup, &yc_down)) goto label_end;
    dbgrid->rankToIndice(iech, _ind2);
    dbgrid->indicesToCoordinateInPlace(_ind2, _xyz);

    if (y[0] >= t1max)
      facies = SHADOW_ISLAND;
    else
    {
      flag_shadow = 0;
      dbgrid->rankToIndice(iech, _ind2);
      dbgrid->indicesToCoordinateInPlace(_ind2, _xyz);
      for (istep = 1; istep <= nstep && flag_shadow == 0; istep++)
      {
        for (idim = 0; idim < ndim; idim++)
          _xyz[idim] -= del[idim];
        flag = point_to_grid(dbgrid, _xyz.data(), 0, _ind2.data());
        if (flag > 0) break;
        if (flag < 0) continue;
        ys = _st_grid_eval(dbgrid, isimu, icase, nbsimu, _xyz);
        if (FFFF(ys)) continue;
        jech = dbgrid->indiceToRank(_ind2);
        if (rule_thresh_define_shadow(propdef, dbgrid, this, SHADOW_WATER, jech,
                                      isimu, nbsimu, &s1min, &s1max, &s2min,
                                      &s2max, &sh_dsup, &sh_down)) return (1);
        if (ys < s1max) continue; /* Upstream point not in island */
        seuil = t1max - yc_down + dy * istep;
        flag_shadow = (MIN(ys,s1max + sh_dsup) > seuil);
      }
      facies = (flag_shadow) ? SHADOW_SHADOW :
                               SHADOW_WATER;
    }

    /* Combine the underlying GRFs to derive Facies */

    dbgrid->setSimvar(ELoc::FACIES, iech, isimu, 0, ipgs, nbsimu, 1, facies);
  }

  /* Set the error return code */

  error = 0;

  label_end:
  db_vector_free(del);
  return (error);
}

/****************************************************************************/
/*!
 **  Set the bounds and possibly add replicates (Shadow)
 **
 ** \return  Error return code
 **
 ** \param[in]  propdef    PropDef structure
 ** \param[in]  dbin       Db structure
 ** \param[in]  dbout      Db grid structure
 ** \param[in]  isimu      Rank of the simulation (if EProcessOper::CONDITIONAL)
 ** \param[in]  igrf       Rank of the GRF
 ** \param[in]  ipgs       Rank of the GS
 ** \param[in]  nbsimu     Number of simulations (if EProcessOper::CONDITIONAL)
 **
 *****************************************************************************/
int RuleShadow::evaluateBounds(PropDef *propdef,
                               Db *dbin,
                               Db *dbout,
                               int isimu,
                               int igrf,
                               int ipgs,
                               int nbsimu) const
{
  int iech, jech, nadd, nech, idim, facies, nstep, istep, valid;
  double dist, t1min, t1max, t2min, t2max, s1min, s1max, s2min, s2max;
  double dinc, seuil, alea, sh_dsup, sh_down, yc_down, dval;
  const DbGrid* dbgrid;

  /* Initializations */

  if (dbin == nullptr) return (0);
  nadd = jech = 0;
  nech = dbin->getSampleNumber();
  dist = 0.;
  dinc = getIncr();
  nstep = (int) floor(getDMax() / dinc);
  seuil = s1min = s1max = s2min = s2max = TEST;
  dbgrid = dynamic_cast<const DbGrid*>(dbout);

  /* Case of the shadow */

  if (igrf == 1) return (0);

  /* Loop on the data */
  for (iech = 0; iech < nech; iech++)
  {
    /* Convert the proportions into thresholds for data point */
    if (!dbin->isActive(iech)) continue;
    if (!point_inside_grid(dbin, iech, dbgrid)) continue;
    facies = (int) dbin->getZVariable(iech, 0);
    if (rule_thresh_define_shadow(propdef, dbin, this, facies, iech, isimu,
                                  nbsimu, &t1min, &t1max, &t2min, &t2max,
                                  &sh_dsup, &sh_down)) return (1);
    yc_down = sh_down;
    dbin->setLocVariable(ELoc::L,iech, get_rank_from_propdef(propdef, ipgs, igrf),
                        t1min);
    dbin->setLocVariable(ELoc::U,iech, get_rank_from_propdef(propdef, ipgs, igrf),
                        t1max);

    /* The data belongs to the island, no replicate */

    if (facies == SHADOW_ISLAND) continue;

    /* In the case of data belonging to the SHADOW */
    /* Generate one replicate in ISLAND at uniform distance upstream */

    if (facies == SHADOW_SHADOW)
    {
      /* Add one replicate */
      jech = dbin->addSamples(1, 0.);
      if (jech < 0) return (1);

      /* Set the coordinates of the replicate */
      /* - at a point where proportions are known */
      /* - after truncation, the point can create shadow at target */

      alea = 1;
      valid = 0;
      while (!valid)
      {
        dist = 0.;
        alea = law_uniform(0., 1.);
        for (idim = 0; idim < dbin->getNDim(); idim++)
        {
          dval = alea * getDMax() * _shift[idim];
          dbin->setCoordinate(jech, idim,
                              dbin->getCoordinate(iech, idim) - dval);
          dist += dval * dval;
        }
        dist = sqrt(dist);

        /* Can the replicate be added */

        if (replicateInvalid(dbin, dbout, jech))
        {
          (void) dbin->deleteSample(jech);
          continue;
        }

        /* Get proportion at the tentative replicate */
        if (rule_thresh_define_shadow(propdef, dbin, this, facies, jech, isimu,
                                      nbsimu, &s1min, &s1max, &s2min, &s2max,
                                      &sh_dsup, &sh_down))
        {
          (void) dbin->deleteSample(jech);
          return (1);
        }
        seuil = t1max - yc_down + dist * getTgte();
        if (seuil > s1max + sh_dsup) continue;
        valid = 1;
      }

      /* Set the attributes of the replicate */
      dbin->setLocVariable(ELoc::Z,jech, 0, SHADOW_ISLAND);
      dbin->setLocVariable(ELoc::L,jech, get_rank_from_propdef(propdef, ipgs, igrf),
                          MAX(seuil, s1max));
      dbin->setLocVariable(ELoc::U,jech, get_rank_from_propdef(propdef, ipgs, igrf),
      THRESH_SUP);
      nadd++;
    }

    /* In the case of data belonging to the WATER            */
    /* Generate series of replicates (maximum number: nstep) */
    /* whose "elevation" which will not create any shadow    */

    if (facies == SHADOW_WATER)
    {
      /* Loop on the replicates */
      for (istep = 1; istep <= nstep; istep++)
      {
        jech = dbin->addSamples(1, 0.);
        if (jech < 0) return (1);

        /* Set the coordinates of the replicate */
        dist = 0.;
        for (idim = 0; idim < dbin->getNDim(); idim++)
        {
          dval = dinc * _shift[idim] * istep;
          dbin->setCoordinate(jech, idim,
                              dbin->getCoordinate(iech, idim) - dval);
          dist += dval * dval;
        }
        dist = sqrt(dist);

        /* Can the replicate be added */
        if (replicateInvalid(dbin, dbout, jech))
        {
          (void) dbin->deleteSample(jech);
          continue;
        }

        /* Get proportion at the tentative replicate */
        if (rule_thresh_define_shadow(propdef, dbin, this, facies, jech, isimu,
                                      nbsimu, &s1min, &s1max, &s2min, &s2max,
                                      &sh_dsup, &sh_down))
        {
          (void) dbin->deleteSample(jech);
          return (1);
        }

        /* Set the attributes of the replicate */
        seuil = t1max - yc_down + dist * getTgte();
        if (seuil > s1max + sh_dsup)
        {
          /* The replicate is not necessary */
          (void) dbin->deleteSample(jech);
          continue;
        }

        dbin->setLocVariable(ELoc::Z,jech, 0, SHADOW_IDLE);
        dbin->setLocVariable(ELoc::L,jech, get_rank_from_propdef(propdef, ipgs, igrf),
                             THRESH_INF);
        dbin->setLocVariable(ELoc::U,jech, get_rank_from_propdef(propdef, ipgs, igrf),
                            MAX(seuil, s1max));
        nadd++;
      }
    }
    jech++;
  }

  if (igrf == 0 && nadd > 0)
  {
    message("Initial count of data = %d\n", nech);
    message("Number of replicates  = %d\n", nadd);
  }
  return (0);
}

void RuleShadow::_normalizeShift()
{
  if (!_shift.empty()) VH::normalize(_shift);
}
