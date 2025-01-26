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
#include "Calculators/ACalcInterpolator.hpp"
#include "Calculators/ACalcDbToDb.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Db/Db.hpp"
#include "Db/DbHelper.hpp"
#include "Model/Model.hpp"
#include "Neigh/ANeigh.hpp"

ACalcInterpolator::ACalcInterpolator()
  : ACalcDbToDb()
  , _model(nullptr)
  , _neigh(nullptr)
  , _ncova(0)
{
}

ACalcInterpolator::~ACalcInterpolator()
{
}

bool ACalcInterpolator::_setNCov(int ncova)
{
  if (ncova <= 0) return true;
  if (_ncova <= 0)
    _ncova = ncova;
  else
  {
    if (_ncova != ncova)
    {
      messerr("Inconsistent Covariance function Number:");
      messerr("- Number already defined = %d", _ncova);
      messerr("- Number of variables newly declared = %d", ncova);
      return false;
    }
  }
  return true;
}

bool ACalcInterpolator::_check()
{
  if (! ACalcDbToDb::_check()) return false;

  /**************************************************/
  /* Cross-checking the Space Dimension consistency */
  /**************************************************/

  int ndim = _getNDim();
  if (_model != nullptr)
  {
    if (ndim > 0)
    {
      if (ndim != _model->getNDim())
      {
        messerr("Inconsistent Space dimension:");
        messerr("- Current dimension = %d",ndim);
        messerr("- Space Dimension of 'model' = %d",_model->getNDim());
        return false;
      }
    }
    else
    {
      ndim = _model->getNDim();
    }
  }

  if (_neigh != nullptr)
  {
    if (ndim > 0)
    {
      if (ndim != (int) _neigh->getNDim())
      {
        messerr("Inconsistent Space dimension:");
        messerr("- Current dimension = %d",ndim);
        messerr("- Space Dimension of '_neigh' = %d",_neigh->getNDim());
        return false;
      }
    }
    else
    {
      ndim = (int) _neigh->getNDim();
    }
    // Attach the input and output files
    _neigh->attach(getDbin(), getDbout());
  }

  /**************************************************/
  /* Cross-Checking the Variable Number consistency */
  /**************************************************/

  int nvar = _getNVar();
  if (_model != nullptr)
  {
    if (nvar > 0)
    {
      if (nvar != _model->getNVar())
      {
        messerr("Inconsistent Variable Number:");
        messerr("- Current number = %d",nvar);
        messerr("- Number of variables in 'model' = %d",_model->getNVar());
        return false;
      }
    }
    else
    {
//      nvar = _model->getNVar(); // Never reached
    }
  }

  /************************************************************/
  /* Cross-Checking the number of External Drifts consistency */
  /************************************************************/

  int nfex = 0;
  if (_model != nullptr)
  {
    nfex = _model->getNExtDrift();
    if (nfex > 0)
    {
      // No check needs to be performed on the Input file as
      // the possibly missing variables will be expanded from the Output File
      // during the preprocessing step

      if (hasDbout(false))
      {
        if (getDbout()->getNLoc(ELoc::F) != nfex)
        {
          messerr("The model requires %d external drift(s)", nfex);
          messerr("but the output Db refers to %d external drift variables",
                  getDbout()->getNLoc(ELoc::F));
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
    if (_model->getNCov() <= 0)
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
    if (hasDbin(false))  getDbin()->getExtensionInPlace(db_mini, db_maxi, true);
    if (hasDbout(false)) getDbout()->getExtensionInPlace(db_mini, db_maxi, true);
    _model->setField(VH::extensionDiagonal(db_mini, db_maxi));
  }
  return true;
}

bool ACalcInterpolator::_preprocess()
{
  if (!ACalcDbToDb::_preprocess()) return false;

  // Space dimension

  if (_model != nullptr)
  {
    if (!_setNdim(_model->getNDim())) return false;
  }
  if (_neigh != nullptr)
  {
    if (!_setNdim((int)_neigh->getNDim())) return false;
  }

  // Number of variables

  if (_model != nullptr)
  {
    if (!_setNvar(_model->getNVar())) return false;
  }

  // Number of covariance functions

  if (_model != nullptr)
  {
     if (!_setNCov(_model->getNCov())) return false;
  }

  // Expand information amongst Db if necessary

  if (_model != nullptr && _model->getNExtDrift() > 0)
  {
    if (_expandInformation(1, ELoc::F)) return false;
  }
  if (_expandInformation(1, ELoc::NOSTAT)) return false;

  return true;
}

bool ACalcInterpolator::hasModel(bool verbose) const
{
  if (_model == nullptr)
  {
    if (verbose) messerr("The argument 'model' must be defined");
    return false;
  }
  return true;
}
bool ACalcInterpolator::hasNeigh(bool verbose) const
{
  if (_neigh == nullptr)
  {
    if (verbose) messerr("The argument 'neigh' must be defined");
    return false;
  }
  return true;
}

int ACalcInterpolator::_centerDataToGrid(DbGrid* dbgrid)
{
  int iuid_out = _addVariableDb(1, 2, ELoc::UNKNOWN, 0, _getNDim(), TEST);
  for (int idim = 0; idim < _getNDim(); idim++)
  {
    int iuid_in = getDbin()->getUIDByLocator(ELoc::X, idim);
    getDbin()->duplicateColumnByUID(iuid_in, iuid_out + idim);
    getDbin()->setLocatorByUID(iuid_out + idim, ELoc::X, idim);
  }
  return DbH::centerPointToGrid(getDbin(), dbgrid, 0.);
}
