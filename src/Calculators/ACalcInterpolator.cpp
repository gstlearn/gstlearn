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
#include "Calculators/ACalcInterpolator.hpp"
#include "Calculators/ACalculator.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Neigh/ANeighParam.hpp"

ACalcInterpolator::ACalcInterpolator()
    : ACalcDbToDb(),
      _model(nullptr),
      _neighparam(nullptr)
{
}

ACalcInterpolator::~ACalcInterpolator()
{
}

int ACalcInterpolator::_getNDim() const
{
  int ndim = ACalcDbToDb::_getNDim();
  if (_model != nullptr)
  {
    if (ndim > 0)
    {
      if (ndim != _model->getDimensionNumber()) return -1;
    }
    else
    {
      ndim = _model->getDimensionNumber();
    }
  }
  return ndim;
}

int ACalcInterpolator::_getNVar() const
{
  int nvar = ACalcDbToDb::_getNVar();
  if (_model != nullptr)
  {
    if (nvar > 0)
    {
      if (nvar != _model->getVariableNumber()) return -1;
    }
    else
    {
      nvar = _model->getVariableNumber();
    }
  }
  return nvar;
}

int ACalcInterpolator::_getNCova() const
{
  if (_model == nullptr) return -1;
  return _model->getCovaNumber();
}

bool ACalcInterpolator::_check()
{
  if (! ACalcDbToDb::_check()) return false;

  /**************************************************/
  /* Cross-checking the Space Dimension consistency */
  /**************************************************/

  int ndim = ACalcDbToDb::_getNDim();
  if (_model != nullptr)
  {
    if (ndim > 0)
    {
      if (ndim != _model->getDimensionNumber())
      {
        messerr("Inconsistent Space dimension:");
        messerr("- Current dimension = %d",ndim);
        messerr("- Space Dimension of 'model' = %d",_model->getDimensionNumber());
        return false;
      }
    }
    else
    {
      ndim = _model->getDimensionNumber();
    }
  }

  if (_neighparam != nullptr)
  {
    if (ndim > 0)
    {
      if (ndim != _neighparam->getNDim())
      {
        messerr("Inconsistent Space dimension:");
        messerr("- Current dimension = %d",ndim);
        messerr("- Space Dimension of '_neighparam' = %d",_neighparam->getNDim());
        return false;
      }
    }
    else
    {
      ndim = _neighparam->getNDim();
    }
  }

  /**************************************************/
  /* Cross-Checking the Variable Number consistency */
  /**************************************************/

  int nvar = ACalcDbToDb::_getNVar();
  if (_model != nullptr)
  {
    if (nvar > 0)
    {
      if (nvar != _model->getVariableNumber())
      {
        messerr("Inconsistent the Variable Number:");
        messerr("- Current number = %d",nvar);
        messerr("- Number of variables in 'model' = %d",_model->getVariableNumber());
        return false;
      }
    }
    else
    {
      nvar = _model->getVariableNumber();
    }
  }

  /************************************************************/
  /* Cross-Checking the number of External Drifts consistency */
  /************************************************************/

  int nfex = 0;
  if (_model != nullptr)
  {
    nfex = _model->getExternalDriftNumber();
    if (nfex > 0)
    {
      if (hasDbin(false))
      {
        if (getDbin()->getExternalDriftNumber() != nfex)
        {
          messerr("The _model requires %d external drift(s)", nfex);
          messerr("but the input Db refers to %d external drift variables",
                  getDbin()->getExternalDriftNumber());
          return false;
        }
      }

      if (hasDbout(false))
      {
        if (getDbout()->getExternalDriftNumber() != nfex)
        {
          messerr("The _model requires %d external drift(s)", nfex);
          messerr("but the output Db refers to %d external drift variables",
                  getDbout()->getExternalDriftNumber());
          return false;
        }
      }
    }
  }

  /**************************************/
  /* Checking the Validity of the _model */
  /**************************************/

  if (_model != nullptr)
  {
    if (_model->getCovaNumber() <= 0)
    {
      messerr("The number of covariance must be positive");
      return false;
    }
  }

  /*********************************/
  /* Calculate the field extension */
  /*********************************/

  if (_model != nullptr)
  {
    VectorDouble db_mini(ndim, TEST);
    VectorDouble db_maxi(ndim, TEST);
    if (hasDbin(false))  getDbin()->getExtensionInPlace(db_mini, db_maxi);
    if (hasDbout(false)) getDbout()->getExtensionInPlace(db_mini, db_maxi);
    _model->setField(ut_vector_extension_diagonal(db_mini, db_maxi));
  }
  return true;
}

bool ACalcInterpolator::hasModel(bool verbose) const
{
  if (_model == nullptr)
  {
    if (verbose)
      messerr("The argument 'model' must be defined");
    return false;
  }
  return true;
}
bool ACalcInterpolator::hasNeighParam(bool verbose) const
{
  if (_neighparam == nullptr)
  {
    if (verbose)
      messerr("The argument 'neighparam' must be defined");
    return false;
  }
  return true;
}

