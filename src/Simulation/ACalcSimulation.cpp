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
#include "Simulation/ACalcSimulation.hpp"
#include "Calculators/ACalcInterpolator.hpp"

namespace gstlrn
{

ACalcSimulation::ACalcSimulation(Id nbsimu, Id seed)
  : ACalcInterpolator()
  , _nbsimu(nbsimu)
  , _seed(seed)
{
}

ACalcSimulation::~ACalcSimulation()
{
}

bool ACalcSimulation::_check()
{
  if (!ACalcInterpolator::_check()) return false;

  if (getNbSimu() <= 0)
  {
    messerr("You must define 'nbsimu' and 'nbtuba'");
    return false;
  }
  return true;
}

bool ACalcSimulation::_preprocess()
{
  return ACalcInterpolator::_preprocess();
}

Id ACalcSimulation::_compute(Db* db, Id icase, Id shift)
{
  DECLARE_UNUSED(db);
  DECLARE_UNUSED(icase);
  DECLARE_UNUSED(shift);
  messerr("The method 'compute' should be implemented in derived classes");
  return 1;
}

/*****************************************************************************/
/*!
 **  Perform non-conditional simulations on a set of gradient points using
 **  Turning Bands method.
 **
 ** \param[in]  dbgrd      Gradient Db structure
 ** \param[in]  delta      Value of the increment
 **
 ** \remarks The simulated gradients are stored as follows:
 ** \remarks idim * nbsimu + isimu (for simulation at first point)
 ** \remarks idim * nbsimu + isimu + ndim * nbsimu (for simulation at 2nd point)
 ** \remarks At the end, the simulated gradient is stored at first point
 **
 *****************************************************************************/
void ACalcSimulation::_computeGradient(Db* dbgrd, double delta)
{
  Id jsimu;
  Id icase               = 0;
  Id ndim                = dbgrd->getNDim();
  auto nbsimu            = getNbSimu();
  VectorBool activeArray = dbgrd->getActiveArray();

  for (Id idim = 0; idim < ndim; idim++)
  {

    /* Simulation at the initial location */

    for (Id isimu = 0; isimu < nbsimu; isimu++)
    {
      jsimu = isimu + idim * nbsimu;
      _compute(dbgrd, icase, jsimu);
    }

    /* Shift the information */

    for (Id iech = 0; iech < dbgrd->getNSample(); iech++)
      if (activeArray[iech])
        dbgrd->setCoordinate(iech, idim, dbgrd->getCoordinate(iech, idim) + delta);

    /* Simulation at the shift location */

    for (Id isimu = 0; isimu < nbsimu; isimu++)
    {
      jsimu = isimu + idim * nbsimu + ndim * nbsimu;
      _compute(dbgrd, icase, jsimu);
    }

    /* Un-Shift the information */

    for (Id iech = 0; iech < dbgrd->getNSample(); iech++)
      if (activeArray[iech])
        dbgrd->setCoordinate(iech, idim, dbgrd->getCoordinate(iech, idim) - delta);

    /* Scaling */

    for (Id isimu = 0; isimu < nbsimu; isimu++)
      for (Id iech = 0; iech < dbgrd->getNSample(); iech++)
      {
        if (!activeArray[iech]) continue;
        jsimu         = isimu + idim * nbsimu + ndim * nbsimu;
        double value2 = dbgrd->getSimvar(ELoc::SIMU, iech, jsimu, 0, icase,
                                         2 * ndim * nbsimu, 1);
        jsimu         = isimu + idim * nbsimu;
        double value1 = dbgrd->getSimvar(ELoc::SIMU, iech, jsimu, 0, icase,
                                         2 * ndim * nbsimu, 1);
        dbgrd->setSimvar(ELoc::SIMU, iech, jsimu, 0, icase,
                         2 * ndim * nbsimu,
                         1, (value2 - value1) / delta);
      }
  }
}

/*****************************************************************************/
/*!
 **  Perform non-conditional simulations on a set of tangent points using
 **  Turning Bands method.
 **
 ** \param[in]  dbtgt      Tangent Db structure
 ** \param[in]  delta      Value of the increment
 **
 ** \remarks Warning: To perform the simulation of the tangent, we must
 ** \remarks simulated the gradients first. So we need to dimension the
 ** \remarks simulation outcome variables as for the gradients
 **
 *****************************************************************************/
void ACalcSimulation::_computeTangent(Db* dbtgt, double delta)
{
  Id icase               = 0;
  auto nvar              = _getNVar();
  auto nbsimu            = getNbSimu();
  VectorBool activeArray = dbtgt->getActiveArray();

  /* Perform the simulation of the gradients at tangent points */

  _computeGradient(dbtgt, delta);

  /* Calculate the simulated tangent */

  for (Id isimu = 0; isimu < nbsimu; isimu++)
    for (Id iech = 0; iech < dbtgt->getNSample(); iech++)
    {
      if (!activeArray[iech]) continue;

      double value = 0.;
      for (Id idim = 0; idim < dbtgt->getNDim(); idim++)
        value += dbtgt->getLocVariable(ELoc::TGTE, iech, idim) *
                 dbtgt->getSimvar(ELoc::SIMU, iech, isimu, 0, icase, nbsimu, nvar);
      dbtgt->setSimvar(ELoc::SIMU, iech, isimu, 0, icase, nbsimu, nvar, value);
    }
}

} // namespace gstlrn