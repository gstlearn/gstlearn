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
#include "Basic/VectorHelper.hpp"
#include "Calculators/ACalcDbToDb.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Db/Db.hpp"
#include "Db/DbHelper.hpp"
#include "Model/Model.hpp"
#include "Neigh/ANeigh.hpp"

namespace gstlrn
{
  ACalcInterpolator::ACalcInterpolator()
    : ACalcDbToDb()
    , _modelGeneric(nullptr)
    , _neigh(nullptr)
    , _krigopt()
    , _ncova(0)
  {
  }

  ACalcInterpolator::~ACalcInterpolator() {}

  bool ACalcInterpolator::_check()
  {
    if (!ACalcDbToDb::_check()) return false;

    /**************************************************/
    /* Cross-checking the Space Dimension consistency */
    /**************************************************/

    size_t ndim = _getNDim();
    if (_modelGeneric != nullptr)
    {
      if (ndim > 0)
      {
        if (ndim != _modelGeneric->getNDim())
        {
          messerr("Inconsistent Space dimension:");
          messerr("- Current dimension = %d", ndim);
          messerr("- Space Dimension of 'model' = %d",
                  _modelGeneric->getNDim());
          return false;
        }
      }
      else
      {
        ndim = _modelGeneric->getNDim();
      }
    }

    if (_neigh != nullptr)
    {
      if (ndim > 0)
      {
        if (ndim != _neigh->getNDim())
        {
          messerr("Inconsistent Space dimension:");
          messerr("- Current dimension = %d", ndim);
          messerr("- Space Dimension of '_neigh' = %d", _neigh->getNDim());
          return false;
        }
      }
      else
      {
        ndim = static_cast<Id>(_neigh->getNDim());
      }
      // Attach the input and output files
      _neigh->attach(getDbin(), getDbout());
    }

    /**************************************************/
    /* Cross-Checking the Variable Number consistency */
    /**************************************************/

    auto nvar = _getNVar();
    if (_modelGeneric != nullptr)
    {
      Id nvarModel = _modelGeneric->getNVar();
      if (nvar > 0)
      {
        if (nvar != nvarModel)
        {
          messerr("Inconsistent Variable Number:");
          messerr("- Current number = %d", nvar);
          messerr("- Number of variables in 'model' = %d", nvarModel);
          return false;
        }
      }

      // Check consistency with the number of Z-variables present in _dbin (if defined)
      // already stored in _getNVar() and the number of variables in the Model: they should match
      if (getDbin() != nullptr)
      {
        if (nvar != nvarModel)
        {
          messerr("The number of variables defined in Dbin (%d)", nvar);
          messerr(
            "should match the number of variables defined in the Model (%d)",
            nvarModel);
          return false;
        }
      }
      _setNvar(nvarModel);
    }

    /************************************************************/
    /* Cross-Checking the number of External Drifts consistency */
    /************************************************************/

    Id nfex = 0;
    if (_modelGeneric != nullptr)
    {
      nfex = _modelGeneric->getNExtDrift();
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

    /***************************************/
    /* Checking the Validity of the _model */
    /***************************************/

    if (_modelGeneric != nullptr)
    {
      if (_modelGeneric->getCov() == nullptr)
      {
        messerr("The Model must contain at least one Covariance");
        return false;
      }
    }

    /******************************/
    /* Cross-checking the Krigopt */
    /******************************/

    if (!_krigopt.isCorrect(getDbout(), _neigh, _modelGeneric)) return false;

    /*********************************/
    /* Calculate the field extension */
    /*********************************/

    if (_modelGeneric != nullptr)
    {
      VectorDouble db_mini(ndim, TEST);
      VectorDouble db_maxi(ndim, TEST);
      if (hasDbin(false))
        getDbin()->getExtensionInPlace(db_mini, db_maxi, true);
      if (hasDbout(false))
        getDbout()->getExtensionInPlace(db_mini, db_maxi, true);
      _modelGeneric->setField(VH::extensionDiagonal(db_mini, db_maxi));
    }
    return true;
  }

  bool ACalcInterpolator::_preprocess()
  {
    if (!ACalcDbToDb::_preprocess()) return false;

    // Space dimension

    if (_modelGeneric != nullptr)
    {
      if (!_setNdim(static_cast<Id>(_modelGeneric->getNDim()))) return false;
    }
    if (_neigh != nullptr)
    {
      if (!_setNdim(static_cast<Id>(_neigh->getNDim()))) return false;
    }

    // Number of variables

    if (_modelGeneric != nullptr)
    {
      if (!_setNvar(_modelGeneric->getNVar())) return false;
    }

    // Number of covariance functions

    if (_modelGeneric != nullptr) _ncova = _calculateNCova();

    // Expand information amongst Db if necessary

    if (_modelGeneric != nullptr && _modelGeneric->getNExtDrift() > 0)
    {
      if (_expandInformation(1, ELoc::F)) return false;
    }
    if (_expandInformation(1, ELoc::NOSTAT)) return false;

    return true;
  }

  Id ACalcInterpolator::_calculateNCova()
  {
    const auto* modelcovlist = dynamic_cast<const ModelCovList*>(_modelGeneric);
    if (modelcovlist == nullptr) return 0;
    const auto* covanisolist =
      dynamic_cast<const CovAnisoList*>(modelcovlist->getCovList());
    if (covanisolist == nullptr) return 0;
    _ncova = modelcovlist->getNCov();
    return _ncova;
  }

  void ACalcInterpolator::setModelGeneric(ModelGeneric* modelGeneric)
  {
    _modelGeneric = modelGeneric;
  }

  bool ACalcInterpolator::hasModelGeneric(bool verbose) const
  {
    if (_modelGeneric == nullptr)
    {
      if (verbose) messerr("The argument 'modelGeneric' must be defined");
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

  Id ACalcInterpolator::_centerDataToGrid(DbGrid* dbgrid)
  {
    Id iuid_out = _addVariableDb(1, 2, ELoc::UNDEFINED, 0, _getNDim(), TEST);
    for (Id idim = 0; idim < _getNDim(); idim++)
    {
      Id iuid_in = getDbin()->getUIDByLocator(ELoc::X, idim);
      getDbin()->duplicateColumnByUID(iuid_in, iuid_out + idim);
      getDbin()->setLocatorByUID(iuid_out + idim, ELoc::X, idim);
    }
    return DbH::centerPointToGrid(getDbin(), dbgrid, 0.);
  }
} // namespace gstlrn
