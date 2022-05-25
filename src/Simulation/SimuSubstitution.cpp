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
#include "Simulation/ASimulation.hpp"
#include "Simulation/SimuSubstitution.hpp"
#include "Simulation/SimuSubstitutionParam.hpp"
#include "Basic/Law.hpp"

#include <math.h>

SimuSubstitution::SimuSubstitution(int nbsimu, int seed)
    : ASimulation(nbsimu, seed),
      AStringable()
{
}

SimuSubstitution::SimuSubstitution(const SimuSubstitution &r)
    : ASimulation(r),
      AStringable(r)
{
}

SimuSubstitution& SimuSubstitution::operator=(const SimuSubstitution &r)
{
  if (this != &r)
  {
    ASimulation::operator=(r);
    AStringable::operator =(r);
  }
  return *this;
}

SimuSubstitution::~SimuSubstitution()
{
}

String SimuSubstitution::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  return sstr.str();
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
                               int verbose)
{
  law_set_random_seed(getSeed());
  int error = 1;
  VectorDouble props;
  bool flag_local = false;
  bool flag_angloc = false;
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

    _planes.resize(np);
    _poissonPlanesGenerate(dbgrid);

    /* Assigning a value to the half-space that contains the center */

    for (int ip = 0; ip < np; ip++)
    {
      if (! subparam.isFlagOrient())
      {
        int ival = (_planes[ip].getRndval() > 0.5) ? -1 : 1;
        _planes[ip].setValue((double) ival);
      }
      else if (! flag_local)
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

        if (flag_local)
        {
          double factor;
          VectorDouble vector(3);
          if (subparam.getColfac() >= 0)
          {
            factor = dbgrid->getArray(iech, subparam.getColfac());
            subparam.isValidFactor(&factor);
          }
          if (flag_angloc)
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

//  if (subparam.isFlagCoding())
//  {
//    if (subparam.isFlagAuto()) nstates = np;
//    if (subparam.isFlagDirect() && nstates != np)
//    {
//      message("You have used the internal information process\n");
//      message("The number of states should be equal to %d\n", np);
//      message("Nevertheless, your choice prevails\n");
//    }
//    props = trans_to_props(nfacies, verbose, trans);
//    status = (int*) mem_alloc(sizeof(int) * nstates, 1);
//
//    /* Simulation of the initial state */
//
//    u = law_uniform(0., 1.);
//    w0 = ie = 0;
//    while (w0 < u)
//      w0 += props[ie++];
//    status[0] = ie - 1;
//
//    /* Simulation of the current state */
//
//    for (ip = 1; ip < nstates; ip++)
//    {
//      u = law_uniform(0., 1.);
//      p0 = je = 0;
//      ie = status[ip - 1];
//      while (p0 < u)
//      {
//        p0 += TRANS(ie, je);
//        je++;
//      }
//      status[ip] = je - 1;
//    }
//
//    /* Simulating the directing function */
//
//    for (iech = 0; iech < dbgrid->getSampleNumber(); iech++)
//    {
//      if (!dbgrid->isActive(iech)) continue;
//      value = (flag_direct) ? dbgrid->getArray(iech, iptr) :
//                              dbgrid->getVariable(iech, 0);
//      ival = (int) ((value - vmin) / (vmax - vmin) * nstates);
//      if (ival < 0) ival = 0;
//      if (ival >= nstates) ival = nstates - 1;
//      dbgrid->setArray(iech, iptr, 1 + status[ival]);
//    }
//
//    /* Core deallocation */
//
//    status = (int*) mem_free((char* ) status);
//    props = (double*) mem_free((char* ) props);
//
//    /* Printout statistics */
//
//    if (verbose)
//    {
//      message("\nCoding process: \n");
//      message("Number of coded states     = %d \n", nstates);
//      message("Minimum information value  = %lf\n", vmin);
//      message("Maximum information value  = %lf\n", vmax);
//    }
//  }
//
//  /* Set the error return code */
//
//  error = 0;
//
//  /* Core deallocation */
//
//  label_end:
//  indg = db_indg_free(indg);
//  status = (int*) mem_free((char* ) status);
//  props = (double*) mem_free((char* ) props);
  return (error);
}

/****************************************************************************/
/*!
 **  Generate the Poisson planes that cover the grid
 **
 ** \param[in]  dbgrid   Db corresponding to the target grid
 **
 ** \param[out] splanes  SubPlanes structure
 **
 ** \remarks  The array 'planes' contains successively a,b,c,d such that
 ** \remarks  ax + by + cz + d = 0
 ** \remarks  The valuation of each line is assigned a uniform value [0,1]
 **
 *****************************************************************************/
void SimuSubstitution::_poissonPlanesGenerate(DbGrid *dbgrid)
{
  double ap[3];

  VectorDouble center = dbgrid->getCenter();
  double diagonal = dbgrid->getExtensionDiagonal();

  /* Loop on the planes to be generated */

  for (int ip = 0; ip < (int) _planes.size(); ip++)
  {
    double d0 = diagonal * law_uniform(-1., 1.) / 2.;
    double u = 0.;
    for (int idim = 0; idim < 3; idim++)
    {
      ap[idim] = law_gaussian();
      u += ap[idim] * ap[idim];
    }
    u = sqrt(u);
    for (int idim = 0; idim < 3; idim++)
    {
      ap[idim] /= u;
      d0 -= ap[idim] * center[idim];
    }
    if (d0 < 0)
    {
      for (int idim = 0; idim < 3; idim++)
        ap[idim] = -ap[idim];
      d0 = -d0;
    }

    /* Storing the plane */

    for (int idim = 0; idim < 3; idim++)
      _planes[ip].setCoor(idim, ap[idim]);
    _planes[ip].setIntercept(d0);
    _planes[ip].setRndval(law_uniform(0., 1.));
  }
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

