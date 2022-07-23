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
#include "Boolean/AShape.hpp"
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Simulation/SimuSubstitution.hpp"
#include "Simulation/SimuSubstitutionParam.hpp"
#include "Simulation/ACalcSimulation.hpp"
#include "Basic/Law.hpp"

#include <math.h>

#define TRANS(i,j)     (trans[(j) + nfacies * (i)])

SimuSubstitution::SimuSubstitution(int nbsimu, int seed)
    : ACalcSimulation(nbsimu, seed)
{
}

SimuSubstitution::~SimuSubstitution()
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
int SimuSubstitution::simulate(DbGrid *dbgrid,
                               const SimuSubstitutionParam& subparam,
                               int iptr,
                               bool verbose)
{
  law_set_random_seed(getSeed());
  int np = 0;

  /***********************/
  /* Information process */
  /***********************/

  if (subparam.isFlagDirect())
  {

    /* Calculate the number of planes */

    double diagonal = dbgrid->getExtensionDiagonal();
    np = law_poisson(diagonal * subparam.getIntensity() * GV_PI);
    if (np <= 0) return 1;

    /* Generate the Poisson planes */

    _planes = Plane::poissonPlanesGenerate(dbgrid,np);

    /* Assigning a value to the half-space that contains the center */

    for (int ip = 0; ip < np; ip++)
    {
      if (! subparam.isFlagOrient())
      {
        int ival = (_planes[ip].getRndval() > 0.5) ? -1 : 1;
        _planes[ip].setValue((double) ival);
      }
      else if (! subparam.isLocal())
      {
        _calculValue(ip, subparam.getFactor(), subparam.getVector());
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
        for (int i = 0; i < 3; i++)
          prod += _planes[ip].getCoor(i) * cen[i];

        if (subparam.isLocal())
        {
          double factor = 0.;
          VectorDouble vector(3);
          if (subparam.getColfac() >= 0)
          {
            factor = dbgrid->getArray(iech, subparam.getColfac());
            subparam.isValidFactor(&factor);
          }
          if (subparam.isAngleLocal())
          {
            for (int i = 0; i < 3; i++)
              if (subparam.getColang(i) >= 0)
                vector[i] = dbgrid->getArray(iech, subparam.getColang(i));
            subparam.isValidOrientation(vector);
          }
          _calculValue(ip, factor, vector);
        }

        double valloc = _planes[ip].getValue() / 2.;
        valtot += (prod + _planes[ip].getIntercept() > 0) ? valloc : -valloc;
      }
      dbgrid->setArray(iech, iptr, valtot);
    }

    /* Printout statistics on the information process */

    if (verbose)
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
    double value = (subparam.isFlagDirect()) ?
        dbgrid->getArray(iech, iptr) : dbgrid->getVariable(iech, 0);
    if (value < vmin) vmin = value;
    if (value > vmax) vmax = value;
  }
  if (vmin > vmax)
  {
    messerr("No Direction Function has been coded");
    messerr("before the Coding Process takes place");
    return 1;
  }
  np = (int) (vmax - vmin + 0.5);

  /******************/
  /* Coding process */
  /******************/

  int nstates = subparam.getNstates();
  if (subparam.isFlagCoding())
  {
    if (subparam.isFlagAuto()) nstates = np;
    if (subparam.isFlagDirect() && nstates != np)
    {
      message("You have used the internal information process\n");
      message("The number of states should be equal to %d\n", np);
      message("Nevertheless, your choice prevails\n");
    }
    VectorDouble props = _transToProp(subparam, verbose);
    VectorInt status(nstates);

    /* Simulation of the initial state */

    double u = law_uniform(0., 1.);
    double w0 = 0.;
    double ie = 0;
    while (w0 < u)
      w0 += props[ie++];
    status[0] = ie - 1;

    /* Simulation of the current state */

    VectorDouble trans = subparam.getTrans();
    int nfacies = subparam.getNfacies();
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
      double value = (subparam.isFlagDirect()) ? dbgrid->getArray(iech, iptr) :
          dbgrid->getVariable(iech, 0);
      int ival = (int) ((value - vmin) / (vmax - vmin) * nstates);
      if (ival < 0) ival = 0;
      if (ival >= nstates) ival = nstates - 1;
      dbgrid->setArray(iech, iptr, 1 + status[ival]);
    }

    /* Printout statistics */

    if (verbose)
    {
      message("\nCoding process: \n");
      message("Number of coded states     = %d \n", nstates);
      message("Minimum information value  = %lf\n", vmin);
      message("Maximum information value  = %lf\n", vmax);
    }
  }

  return 0;
}

/*****************************************************************************
 **
 ** Calculate the projected value
 **
 ** \param[in,out]  ip      Rank of the Plane
 **
 *****************************************************************************/
void SimuSubstitution::_calculValue(int ip,
                                    double factor,
                                    const VectorDouble& vector)
{
  int ival = ((2. * _planes[ip].getRndval()) > (1. + factor)) ? -1 : 1;
  double cossin = 0.;
  for (int i = 0; i < 3; i++)
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
VectorDouble SimuSubstitution::_transToProp(const SimuSubstitutionParam& subparam,
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

bool SimuSubstitution::_run()
{
  return true;
}
