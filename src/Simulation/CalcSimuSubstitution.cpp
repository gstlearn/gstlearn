/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Boolean/AShape.hpp"
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Simulation/CalcSimuSubstitution.hpp"
#include "Simulation/SimuSubstitutionParam.hpp"
#include "Simulation/ACalcSimulation.hpp"
#include "Basic/Law.hpp"

#include <math.h>

#define TRANS(i,j)     (trans[(j) + nfacies * (i)])

CalcSimuSubstitution::CalcSimuSubstitution(int nbsimu, int seed, bool verbose)
    : ACalcSimulation(nbsimu, seed),
      _verbose(verbose),
      _iattOut(-1),
      _subparam(),
      _planes()
{
}

CalcSimuSubstitution::~CalcSimuSubstitution()
{
}

/*****************************************************************************
 **
 ** Generate a simulation on a regular 3D grid using substitution method
 **
 ** \returns Error return code
 **
 ** \param[in]  dbgrid      Db structure (should be a grid)
 ** \param[in]  subparam    SimuSubstitutionParam structure
 ** \param[in]  iptr        Vector for storing the simulation
 ** \param[in]  verbose     Verbose option
 **
 *****************************************************************************/
bool CalcSimuSubstitution::_simulate()
{
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  int np = 0;

  /***********************/
  /* Information process */
  /***********************/

  if (_subparam.isFlagDirect())
  {

    /* Calculate the number of planes */

    double diagonal = dbgrid->getExtensionDiagonal();
    np = law_poisson(diagonal * _subparam.getIntensity() * GV_PI);
    if (np <= 0) return false;

    /* Generate the Poisson planes */

    _planes = Plane::poissonPlanesGenerate(dbgrid,np);

    /* Assigning a value to the half-space that contains the center */

    for (int ip = 0; ip < np; ip++)
    {
      if (! _subparam.isFlagOrient())
      {
        int ival = (_planes[ip].getRndval() > 0.5) ? -1 : 1;
        _planes[ip].setValue((double) ival);
      }
      else if (! _subparam.isLocal())
      {
        _calculValue(ip, _subparam.getFactor(), _subparam.getVector());
      }
    }

    /* Simulating the directing function */

    for (int iech = 0; iech < dbgrid->getSampleNumber(); iech++)
    {
      VectorDouble cen = dbgrid->getSampleCoordinates(iech);

      /* Loop on the planes */

      double valtot = 0.;
      for (int ip = 0; ip < np; ip++)
      {
        double prod = 0.;
        for (int i = 0; i < (int) cen.size(); i++)
          prod += _planes[ip].getCoor(i) * cen[i];

        if (_subparam.isLocal())
        {
          double factor = 0.;
          VectorDouble vector(3);
          if (_subparam.getColfac() >= 0)
          {
            factor = dbgrid->getArray(iech, _subparam.getColfac());
            _subparam.isValidFactor(&factor);
          }
          if (_subparam.isAngleLocal())
          {
            for (int i = 0; i < 3; i++)
              if (_subparam.getColang(i) >= 0)
                vector[i] = dbgrid->getArray(iech, _subparam.getColang(i));
            _subparam.isValidOrientation(vector);
          }
          _calculValue(ip, factor, vector);
        }

        double valloc = _planes[ip].getValue() / 2.;
        valtot += (prod + _planes[ip].getIntercept() > 0) ? valloc : -valloc;
      }
      dbgrid->setArray(iech, _iattOut, valtot);
    }

    /* Printout statistics on the information process */

    if (_verbose)
    {
      message("\n");
      message("Information Process:\n");
      message("Number of planes generated = %d\n", np);
    }
  }

  /***************************************/
  /* Determination of the extreme values */
  /***************************************/

  double vmin =  1.e30;
  double vmax = -1.e30;
  for (int iech = 0; iech < dbgrid->getSampleNumber(); iech++)
  {
    if (!dbgrid->isActive(iech)) continue;
    double value = (_subparam.isFlagDirect()) ?
        dbgrid->getArray(iech, _iattOut) : dbgrid->getLocVariable(ELoc::Z,iech, 0);
    if (value < vmin) vmin = value;
    if (value > vmax) vmax = value;
  }
  if (vmin > vmax)
  {
    messerr("No Direction Function has been coded");
    messerr("before the Coding Process takes place");
    return false;
  }
  np = (int) (vmax - vmin + 0.5);

  /******************/
  /* Coding process */
  /******************/

  int nstates = _subparam.getNstates();
  if (_subparam.isFlagCoding())
  {
    if (_subparam.isFlagAuto()) nstates = np;
    if (_subparam.isFlagDirect() && nstates != np)
    {
      message("You have used the internal information process\n");
      message("The number of states should be equal to %d\n", np);
      message("Nevertheless, your choice prevails\n");
    }
    VectorDouble props = _transToProp(_subparam, _verbose);
    VectorInt status(nstates);

    /* Simulation of the initial state */

    double u = law_uniform(0., 1.);
    double w0 = 0.;
    int ie = 0;
    while (w0 < u)
      w0 += props[ie++];
    status[0] = ie - 1;

    /* Simulation of the current state */

    VectorDouble trans = _subparam.getTrans();
    int nfacies = _subparam.getNfacies();
    for (int ip = 1; ip < nstates; ip++)
    {
      u = law_uniform(0., 1.);
      double p0 = 0.;
      int je = 0;
      ie = status[ip - 1];
      while (p0 < u)
      {
        p0 += TRANS(ie, je);
        je++;
      }
      status[ip] = je - 1;
    }

    /* Simulating the directing function */

    for (int iech = 0; iech < dbgrid->getSampleNumber(); iech++)
    {
      if (!dbgrid->isActive(iech)) continue;
      double value = (_subparam.isFlagDirect()) ? dbgrid->getArray(iech, _iattOut) :
          dbgrid->getLocVariable(ELoc::Z,iech, 0);
      int ival = (int) ((value - vmin) / (vmax - vmin) * nstates);
      if (ival < 0) ival = 0;
      if (ival >= nstates) ival = nstates - 1;
      dbgrid->setArray(iech, _iattOut, 1 + status[ival]);
    }

    /* Printout statistics */

    if (_verbose)
    {
      message("\nCoding process: \n");
      message("Number of coded states     = %d \n", nstates);
      message("Minimum information value  = %lf\n", vmin);
      message("Maximum information value  = %lf\n", vmax);
    }
  }
  return true;
}

/*****************************************************************************
 **
 ** Calculate the projected value
 **
 ** \param[in,out]  ip      Rank of the Plane
 **
 *****************************************************************************/
void CalcSimuSubstitution::_calculValue(int ip,
                                    double factor,
                                    const VectorDouble& vector)
{
  int ival = ((2. * _planes[ip].getRndval()) > (1. + factor)) ? -1 : 1;
  double cossin = 0.;
  for (int i = 0; i < (int) vector.size(); i++)
    cossin += _planes[ip].getCoor(i) * vector[i];
  if (cossin < 0) ival = -ival;
  _planes[ip].setValue((double) ival);
}

/****************************************************************************/
/*!
 **  Derive proportions from the transition matrix
 **
 ** \return The proportion matrix
 **
 ** \param[in]  subparam SimuSubstitutionParam structure
 ** \param[in]  verbose  Verbose option
 ** \param[in]  eps      Tolerance
 **
 *****************************************************************************/
VectorDouble CalcSimuSubstitution::_transToProp(const SimuSubstitutionParam &subparam,
                                                bool verbose,
                                                double eps)
{
  VectorDouble props;
  VectorDouble propold;
  VectorDouble trans = subparam.getTrans();

  /* Initializations */

  int nfacies = subparam.getNfacies();
  if (nfacies <= 0 || subparam.getTrans().empty()) return props;
  props.resize(nfacies);
  propold.resize(nfacies);

  /* Checks the transition matrix */

  bool flag_error = false;
  for (int i = 0; i < nfacies && !flag_error; i++)
  {
    double total = 0.;
    for (int j = 0; j < nfacies; j++)
      total += ABS(TRANS(i,j));

    if (total <= 0.)
      flag_error = 1;
    else
      for (int j = 0; j < nfacies; j++)
        TRANS(i,j) = ABS(TRANS(i,j)) / total;
  }

  /* Wrong transition matrix: initialize it */

  if (flag_error)
  {
    for (int i = 0; i < nfacies; i++)
      for (int j = 0; j < nfacies; j++)
        TRANS(i,j) = 1. / nfacies;
  }

  /* Initialize the proportion matrix */

  for (int i = 0; i < nfacies; i++)
    props[i] = 1. / nfacies;

  /* Loop to reach the stationarity of the proportions */

  double diff = 2. * eps;
  while (diff > eps)
  {
    double w0 = 0.;
    diff = 0.;
    for (int i = 0; i < nfacies; i++)
      propold[i] = props[i];
    for (int i = 0; i < nfacies; i++)
    {
      double val = 0.;
      for (int j = 0; j < nfacies; j++)
        val += TRANS(j,i) * propold[j];
      props[i] = val;
      w0 += val;
    }
    if (w0 == 0.) w0 = 1.;
    for (int i = 0; i < nfacies; i++)
    {
      props[i] /= w0;
      diff += ABS(propold[i] - props[i]);
    }
  }

  /* Printout the proportions */

  if (verbose)
    print_matrix("Proportions", 0, 1, 1, nfacies, NULL, props.data());

  return props;
}

bool CalcSimuSubstitution::_check()
{
  if (! ACalcSimulation::_check()) return false;

  if (! hasDbout()) return false;
  int ndim = _getNDim();
  if (ndim > 3)
  {
    messerr("The Substitution Method is not a relevant simulation model");
    messerr("for this Space Dimension (%d)", ndim);
    return false;
  }
  if (! getDbout()->isGrid())
  {
    messerr("The argument 'dbout'  should be a grid");
    return false;
  }
  if (! _subparam.isValid(_verbose)) return false;
  return true;
}

bool CalcSimuSubstitution::_preprocess()
{
  if (!ACalcInterpolator::_check()) return false;

  _iattOut = _addVariableDb(2, 1, ELoc::SIMU, 0, 1);
  if (_iattOut < 0) return false;
  return true;
}

bool CalcSimuSubstitution::_run()
{
  law_set_random_seed(getSeed());

  return (_simulate());
}

bool CalcSimuSubstitution::_postprocess()
{
  /* Free the temporary variables */
  _cleanVariableDb(2);

  _renameVariable(2, 1, _iattOut, String(), getNbSimu());
  return true;
}

void CalcSimuSubstitution::_rollback()
{
  _cleanVariableDb(1);
}

/*****************************************************************************
 **
 ** Generate a simulation on a regular 3D grid using substitution method
 **
 ** \returns Error return code
 **
 ** \param[in]  dbgrid      Db structure (should be a grid)
 ** \param[in]  subparam    SimuSubstitutionParam structure
 ** \param[in]  seed        Seed
 ** \param[in]  verbose     Verbose option
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int substitution(DbGrid *dbgrid,
                 SimuSubstitutionParam& subparam,
                 int seed,
                 int verbose,
                 const NamingConvention& namconv)
{
  CalcSimuSubstitution simsub(1, seed, verbose);
  simsub.setDbout(dbgrid);
  simsub.setNamingConvention(namconv);
  simsub.setSubparam(subparam);

  // Run the calculator
  int error = (simsub.run()) ? 0 : 1;
  return error;
}
