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
#include "Estimation/CalcKrigingGradient.hpp"
#include "Anamorphosis/AAnam.hpp"
#include "Covariances/CovAnisoList.hpp"
#include "Covariances/CovGradientAnalytic.hpp"
#include "Covariances/CovGradientGeneric.hpp"
#include "Db/Db.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Estimation/KrigOpt.hpp"
#include "Estimation/KrigingSystem.hpp"
#include "Model/Model.hpp"

namespace gstlrn
{
CalcKrigingGradient::CalcKrigingGradient(bool flag_est,
                                         bool flag_std,
                                         double ball_radius)
  : ACalcInterpolator()
  , _dbGradient(nullptr)
  , _modelGradient(nullptr)
  , _ballRadius(ball_radius)
  , _flagEst(flag_est)
  , _flagStd(flag_std)
  , _flagForceNumeric(false)
{
}

CalcKrigingGradient::~CalcKrigingGradient()
{
  delete _dbGradient;
  delete _modelGradient;
}

bool CalcKrigingGradient::_check()
{
  if (!ACalcInterpolator::_check()) return false;

  if (!hasDbin()) return false;
  if (!hasDbout()) return false;
  if (!hasModel()) return false;
  if (!hasNeigh()) return false;

  // Must be Monovariate
  if (_getNVar() != 1)
  {
    messerr("The Depth and Gradient method is limited to the Monovariate case");
    return false;
  }
  // Must have as many gradient components as the Space dimension
  Id ndim  = _getNDim();
  Id ngrad = getDbin()->getNLoc(ELoc::G);
  if (ndim != ngrad)
  {
    messerr("Number of Gradient components (%d) must be equal to the Space dimension (%d)",
            ngrad, ndim);
    return false;
  }
  return true;
}

/****************************************************************************/
/*!
 **  Update the Input Data Base for Kriging with Gradient components
 **
 ** \remark  In the case of Kriging with Gradient, the gradient
 ** \remark  components are transformed into additional variables
 **
 *****************************************************************************/
void CalcKrigingGradient::_updateDbin()

{
  // Duplicate the Input Db
  _dbGradient = getDbin()->clone();

  // Turn the locators from Gradients into Variables
  _dbGradient->switchLocator(ELoc::G, ELoc::Z);
}

/****************************************************************************/
/*!
 **  Duplicates a Model from another Model for Gradients
 **
 *****************************************************************************/
Id CalcKrigingGradient::_updateModel()

{
  Id ndim = _getNDim();

  // Check if the Model contains a ACov
  const ACov* acov = dynamic_cast<const ACov*>(getModel()->getCov());
  if (acov == nullptr)
  {
    messerr("The input Model should contain an 'ACov'");
    return 0;
  }

  // Create the new Model
  CovContext ctxt(*getModel()->getContext());
  Id new_nvar = 1 + ndim;
  ctxt.setNVar(new_nvar);
  _modelGradient = new ModelGeneric(ctxt);
  if (_modelGradient == nullptr) return 1;

  // Create the Covariance Part
  const auto* covlist = dynamic_cast<const CovAnisoList*>(acov);
  if (covlist == nullptr || _flagForceNumeric)
  {

    // Case where the ACov is not a CovAnisoList: numeric processing
    auto* covnew = new CovGradientGeneric(*acov, _ballRadius);
    _modelGradient->setCov(covnew);
    delete covnew;
  }
  else
  {

    // Case where the Acov is a CovAnisoList: analytical processing
    Id ncov = covlist->getNCov();
    if (ncov != 1)
    {
      messerr("The Gradient Model requires a single Covariance to be defined");
      return 1;
    }

    auto* covnew = new CovGradientAnalytic(*covlist->getCovAniso(0));
    _modelGradient->setCov(covnew);
    delete covnew;
  }

  // *********************************
  // Create the basic drift structures
  // *********************************

  DriftList* drifts = DriftFactory::createDriftListForGradients(getModel()->getDriftList(), ctxt);
  _modelGradient->setDriftList(drifts);
  delete drifts;

  // Set the context
  _modelGradient->setContext(ctxt);
  return 0;
}

bool CalcKrigingGradient::_preprocess()
{
  if (!ACalcInterpolator::_preprocess()) return false;

  // Duplicate the Input Data Base
  // Upgrade the gradient components as new variables
  _updateDbin();

  // Create a New Model for Depth and Gradient components
  return !_updateModel();
}

bool CalcKrigingGradient::_postprocess()
{
  /* Free the temporary variables */
  _cleanVariableDb(2);

  return true;
}

void CalcKrigingGradient::_rollback()
{
  _cleanVariableDb(1);
}

/****************************************************************************/
/*!
 **  Standard Kriging
 **
 ** \return  Error return code
 **
 *****************************************************************************/
bool CalcKrigingGradient::_run()
{
  (void)kriging(_dbGradient, getDbout(), _modelGradient, getNeigh(),
                _flagEst, _flagStd, 0,
                KrigOpt(), getNamingConvention());
  return true;
}

Id krigingGradient(Db* dbin,
                   Db* dbout,
                   ModelGeneric* model,
                   ANeigh* neigh,
                   bool flag_est,
                   bool flag_std,
                   double ball_radius,
                   bool flagForceNumeric,
                   const NamingConvention& namconv)
{
  CalcKrigingGradient krigeGradient(flag_est, flag_std, ball_radius);
  krigeGradient.setDbin(dbin);
  krigeGradient.setDbout(dbout);
  krigeGradient.setModel(model);
  krigeGradient.setNeigh(neigh);
  krigeGradient.setFlagForceNumeric(flagForceNumeric);
  krigeGradient.setNamingConvention(namconv);

  // Run the calculator
  Id error = (krigeGradient.run()) ? 0 : 1;
  return error;
}

} // namespace gstlrn
