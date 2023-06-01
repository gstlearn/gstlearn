/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_old_f.h"
#include "geoslib_enum.h"

#include "Enum/ERule.hpp"

#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Basic/AException.hpp"
#include "LithoRule/RuleShift.hpp"
#include "LithoRule/Rule.hpp"
#include "LithoRule/Node.hpp"
#include "Model/Model.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"

#include <sstream>
#include <math.h>

RuleShift::RuleShift()
    : Rule(),
      _shDsup(0.),
      _shDown(0.),
      _slope(0.),
      _shift(),
      _incr(TEST),
      _xyz(),
      _ind1(),
      _ind2()
{
  setModeRule(ERule::SHIFT);
}

RuleShift::RuleShift(const RuleShift& m)
    :  Rule(m),
      _shDsup(m._shDsup),
      _shDown(m._shDown),
      _slope(m._slope),
      _shift(m._shift),
      _incr(m._incr)
{
}

RuleShift& RuleShift::operator=(const RuleShift& m)
{
  if (this != &m)
  {
    Rule::operator =(m);
    _shDsup = m._shDsup;
    _shDown = m._shDown;
    _slope = m._slope;
    _shift = m._shift;
    _incr = m._incr;
  }
  return *this;
}

RuleShift::~RuleShift()
{
}

/**
 * Definition of the Lithotype RuleShift
 * @param nodes List of "integer" nodes (should only include "S", no "T")
 * @param shift Vector defining the Shift
 */
int RuleShift::resetFromNodes(const VectorInt& nodes, const VectorDouble& shift)
{
  _shift = shift;
  setModeRule(ERule::SHIFT);
  setMainNodeFromNodNames(nodes);
  return 0;
}

int RuleShift::resetFromNames(const VectorString& nodnames, const VectorDouble& shift)
{
  _shift = shift;
  setModeRule(ERule::SHIFT);
  setMainNodeFromNodNames(nodnames);
  return 0;
}

/**
 * Definition of the Lithotype RuleShift
 * @param nfacies Number of facies
 * @param shift Vector defining the shift
 */
int RuleShift::resetFromFaciesCount(int nfacies, const VectorDouble& shift)
{
  _shift = shift;
  setModeRule(ERule::SHIFT);
  VectorString nodnames = buildNodNames(nfacies);
  setMainNodeFromNodNames(nodnames);
  return 0;
}

int RuleShift::resetFromNumericalCoding(const VectorInt& n_type,
                                        const VectorInt& n_facs,
                                        const VectorDouble& shift)
{
  _shift = shift;
  setModeRule(ERule::SHIFT);
  setMainNodeFromNodNames(n_type, n_facs);
  return 0;
}

bool RuleShift::_deserialize(std::istream& is, bool /*verbose*/)
{
  _shift.resize(3);
  bool ret = true;

  ret = ret && Rule::_deserialize(is);

  ret = ret && _recordRead<double>(is, "Slope for Shadow Rule", _slope);
  ret = ret && _recordRead<double>(is, "Lower Threshold for Shadow Rule", _shDown);
  ret = ret && _recordRead<double>(is, "Upper Threshold for Shadow Rule", _shDsup);
  ret = ret && _recordRead<double>(is, "Shift along first direction", _shift[0]);
  ret = ret && _recordRead<double>(is, "Shift along second direction", _shift[1]);
  ret = ret && _recordRead<double>(is, "Shift along third direction", _shift[2]);
  return ret;
}

bool RuleShift::_serialize(std::ostream& os, bool /*verbose*/) const
{
  double slope  = (FFFF(_slope)) ? 0. : _slope;
  double shdown = (FFFF(_shDown)) ? 0. : _shDown;
  double shdsup = (FFFF(_shDsup)) ? 0. : _shDsup;
  VectorDouble shiftloc = _shift;
  shiftloc.resize(3);

  bool ret = true;

  ret = ret && Rule::_serialize(os);

  ret = ret && _recordWrite<double>(os, "Slope for Shadow Rule", slope);
  ret = ret && _recordWrite<double>(os, "Lower Threshold for Shadow Rule", shdown);
  ret = ret && _recordWrite<double>(os, "Upper Threshold for Shadow Rule", shdsup);
  ret = ret && _recordWrite<double>(os, "Shift along first direction", shiftloc[0]);
  ret = ret && _recordWrite<double>(os, "Shift along second direction", shiftloc[1]);
  ret = ret && _recordWrite<double>(os, "Shift along third direction", shiftloc[2]);
  return ret;
}

String RuleShift::displaySpecific() const
{
  std::stringstream sstr;
  sstr << toTitle(2,"Shift Option");
  sstr << toVector("Translation Vector",_shift);
  sstr << "(With the 'Shift' option, only the first GRF is used)" << std::endl;
  return sstr.str();
}

/****************************************************************************/
/*!
**  Define the particularities of the PGS model
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
int RuleShift::particularities(Db* db,
                               const Db* /*dbprop*/,
                               Model* model,
                               int flag_grid_check,
                               int /*flag_stat*/) const
{
  int ndim = (model != nullptr) ? model->getDimensionNumber() : 0;
  VectorDouble wxyz(ndim);
  CovCalcMode mode;
  double rhoval;

  /* Dispatch */

  _xyz.resize(ndim);
  double hval = 0.;
  for (int idim = 0; idim < ndim; idim++)
  {
    _xyz[idim] = _shift[idim];
    hval += _xyz[idim] * _xyz[idim];
  }
  hval = sqrt(hval);

  /* Calculate the covariance between the two GRF */

  for (int idim = 0; idim < ndim; idim++)
    wxyz[idim] = _xyz[idim];
  model_evaluate(model, 0, 0, mode, 1, wxyz, &hval, &rhoval);
  setRho(rhoval);

  /* Translate the shift into grid increments */

  if (_st_shift_on_grid(db, ndim, flag_grid_check)) return (1);
  return (0);
}

int RuleShift::_st_shift_on_grid(Db *db, int ndim, int flag_grid_check) const
{
  _xyz.resize(ndim);
  _ind1.resize(ndim);

  DbGrid* dbgrid = dynamic_cast<DbGrid*>(db);
  if (dbgrid == nullptr)
  {
    if (! flag_grid_check) return(0);
    messerr("The shift Rule requires a Grid Db");
    return(1);
  }

  for (int idim=0; idim<ndim; idim++)
    _xyz[idim] = _shift[idim] + dbgrid->getX0(idim);

  (void) point_to_grid(dbgrid,_xyz.data(),-1,_ind1.data());

  /* Check that the translation is significant */

  int ntot = 0;
  for (int idim=0; idim<ndim; idim++)
    ntot += ABS(_ind1[idim]);
  if (ntot <= 0)
  {
    messerr("The shift of the Lithotype Rule cannot be rendered");
    messerr("using the Output Grid characteristics");
    return(1);
  }
  return(0);
}

bool RuleShift::checkModel(const Model* model, int nvar) const
{
  if (model == nullptr)
  {
    messerr("No Model is provided");
    return false;
  }
  if (nvar > 0 && model->getVariableNumber() != nvar)
  {
    messerr("The number of variables in the Model (%d) does not match",
            model->getVariableNumber());
    messerr(" the number of variables in the Db (%d)", nvar);
    return false;
  }
  return true;
}

/****************************************************************************/
/*!
**  Combine the underlying GRF into a facies value
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
int RuleShift::gaus2facResult(PropDef* propdef,
                              Db* dbout,
                              int* /*flag_used*/,
                              int ipgs,
                              int isimu,
                              int nbsimu) const
{
  int    ndim,iech,jech,idim,igrf,icase;
  double t1min,t1max,t2min,t2max,facies,y[2];

  /* Initializations */

  check_mandatory_attribute("rule_gaus2fac_result",dbout,ELoc::FACIES);
  check_mandatory_attribute("rule_gaus2fac_result",dbout,ELoc::SIMU);
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(dbout);
  if (dbgrid == nullptr) return 1;
  ndim   = dbgrid->getNDim();
  _xyz.resize(ndim);
  _ind1.resize(ndim);
  _ind2.resize(ndim);

  /* Processing the translation */

  for (iech=0; iech<dbgrid->getSampleNumber(); iech++)
  {
    if (! dbgrid->isActive(iech)) continue;

    facies = TEST;
    for (igrf=0; igrf<2; igrf++) y[igrf] = TEST;
    icase = get_rank_from_propdef(propdef, ipgs, 0);
    y[0] = dbgrid->getSimvar(ELoc::SIMU, iech, isimu, 0, icase, nbsimu, 1);
    if (FFFF(y[0])) break;

    if (rule_thresh_define(propdef, dbgrid, this, ITEST, iech, isimu, nbsimu, 1,
                           &t1min, &t1max, &t2min, &t2max)) return 1;
    db_index_sample_to_grid(dbgrid, iech, _ind2.data());
    for (idim = 0; idim < ndim; idim++)
      _ind2[idim] -= _ind1[idim];
    jech = db_index_grid_to_sample(dbgrid, _ind2.data());
    if (jech >= 0)
      y[1] = dbgrid->getSimvar(ELoc::SIMU, jech, isimu, 0, icase, nbsimu, 1);
    else
      y[1] = TEST;
    facies = getFaciesFromGaussian(y[0], y[1]);

    /* Combine the underlying GRFs to derive Facies */

    dbgrid->setSimvar(ELoc::FACIES,iech,isimu,0,ipgs,nbsimu,1,facies);
  }
  return 0;
}

/****************************************************************************/
/*!
**  Set the bounds and possibly add replicates
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
int RuleShift::evaluateBounds(PropDef *propdef,
                              Db *dbin,
                              Db *dbout,
                              int isimu,
                              int igrf,
                              int ipgs,
                              int nbsimu) const
{
  int    iech,jech,nadd,nech,idim,facies,nstep;
  double t1min,t1max,t2min,t2max,s1min,s1max,s2min,s2max;

  /* Initializations */

  if (dbin == nullptr) return(0);
  nadd = nstep = 0;
  nech = dbin->getSampleNumber();

  /* Dispatch */

  if (igrf == 1) return (0);

  /* Loop on the data */
  for (iech = 0; iech < nech; iech++)
  {
    /* Convert the proportions into thresholds for data point */
    if (!dbin->isActive(iech)) continue;
    facies = (int) dbin->getLocVariable(ELoc::Z,iech, 0);
    if (rule_thresh_define(propdef, dbin, this, facies, iech, isimu, nbsimu, 1,
                           &t1min, &t1max, &t2min, &t2max)) return (1);
    dbin->setLocVariable(ELoc::L,iech, get_rank_from_propdef(propdef, ipgs, igrf),
                        t1min);
    dbin->setLocVariable(ELoc::U,iech, get_rank_from_propdef(propdef, ipgs, igrf),
                        t1max);
    if (facies == SHADOW_ISLAND) continue;

    /* Add one replicate */
    jech = dbin->addSamples(1, 0.);
    if (jech < 0) return (1);

    /* Set the coordinates of the replicate */
    for (idim = 0; idim < dbin->getNDim(); idim++)
      dbin->setCoordinate(jech, idim,
                          dbin->getCoordinate(iech, idim) - _shift[idim]);

    /* Can the replicate be added */
    if (replicateInvalid(dbin, dbout, jech))
    {
      (void) dbin->deleteSample(jech);
      return (1);
    }

    /* Convert the proportions into thresholds for replicate */
    if (rule_thresh_define(propdef, dbin, this, facies, jech, isimu, nbsimu, 1,
                           &s1min, &s1max, &s2min, &s2max))
    {
      (void) dbin->deleteSample(jech);
      return (1);
    }

    /* Set the attributes of the replicate */
    if (facies == SHADOW_WATER) dbin->setLocVariable(ELoc::Z,jech, 0, SHADOW_WATER);
    if (facies == SHADOW_SHADOW) dbin->setLocVariable(ELoc::Z,jech, 0, SHADOW_ISLAND);
    dbin->setLocVariable(ELoc::L,jech, get_rank_from_propdef(propdef, ipgs, igrf),
                        s2min);
    dbin->setLocVariable(ELoc::U,jech, get_rank_from_propdef(propdef, ipgs, igrf),
                        s2max);
    nadd++;
  }

  if (igrf == 0 && nadd > 0)
  {
    message("Initial count of data = %d\n", nech);
    message("Number of replicates  = %d\n", nadd);
  }
  return (0);
}

RuleShift* RuleShift::createFromNodes(const VectorInt& nodes,
                                      const VectorDouble& shift)
{
  RuleShift* ruleshift = new RuleShift();
  if (ruleshift->resetFromNodes(nodes, shift))
  {
    messerr("Problem when creating RuleShift from Nodes");
    delete ruleshift;
    return nullptr;
  }
  return ruleshift;
}
RuleShift* RuleShift::createFromNames(const VectorString& nodnames,
                                      const VectorDouble& shift)
{
  RuleShift* ruleshift = new RuleShift();
  if (ruleshift->resetFromNames(nodnames, shift))
  {
    messerr("Problem when creating RuleShift from Node Names");
    delete ruleshift;
    return nullptr;
  }
  return ruleshift;
}
RuleShift* RuleShift::createFromFaciesCount(int nfacies,
                                            const VectorDouble& shift)
{
  RuleShift* ruleshift = new RuleShift();
  if (ruleshift->resetFromFaciesCount(nfacies, shift))
  {
    messerr("Problem when creating RuleShift from Count of Facies");
    delete ruleshift;
    return nullptr;
  }
  return ruleshift;
}
RuleShift* RuleShift::createFromNumericalCoding(const VectorInt& n_type,
                                                const VectorInt& n_facs,
                                                const VectorDouble& shift)
{
  RuleShift* ruleshift = new RuleShift();
  if (ruleshift->resetFromNumericalCoding(n_type, n_facs, shift))
  {
    messerr("Problem when creating RuleShift from Numerical Coding");
    delete ruleshift;
    return nullptr;
  }
  return ruleshift;
}
