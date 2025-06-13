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
#include "Estimation/AModelOptimFactory.hpp"

#include "Covariances/CorAniso.hpp"
#include "Enum/EConsElem.hpp"
#include "Enum/EConsType.hpp"
#include "Estimation/Vecchia.hpp"
#include "Estimation/Likelihood.hpp"
#include "Estimation/AModelOptim.hpp"
#include "Model/Constraints.hpp"
#include "Model/ConsItem.hpp"
#include "Model/Model.hpp"
#include "Model/ModelOptimVMap.hpp"
#include "Model/ModelOptimVario.hpp"
#include "Model/ModelCovList.hpp"
#include "Covariances/CovLMCAnamorphosis.hpp"
#include "Covariances/CovAniso.hpp"
#include "Variogram/Vario.hpp"
#include "Db/DbGrid.hpp"
#include "geoslib_define.h"

static void _modifyMopForAnam(ModelGeneric* model,
                              ModelOptimParam& mop)
{
  const CovLMCAnamorphosis* covanam = dynamic_cast<const CovLMCAnamorphosis*>(model->_getCovModify());
  if (covanam != nullptr)
  {
    EAnam anamtype = covanam->getAnamType();
    if (anamtype != EAnam::HERMITIAN && mop.getFlagGoulard())
      mop.setFlagGoulard(false);
  }
}

/****************************************************************************/
/*!
 **  Define the options of the structure ModelOptimParam
 **
 ** \param[in]  dbmap       Db Grid structure containing the Vmap
 ** \param[in]  model       ModelGeneric structure
 ** \param[in]  constraints Constraints structure
 **
 ** \param[out]  mop        ModelOptimParam structure
 **
 *****************************************************************************/
static int _modifyMopForVMap(const DbGrid* dbmap,
                             ModelGeneric* model,
                             Constraints* constraints,
                             ModelOptimParam& mop)
{
  // Clever setting of options 
  mop.setAuthAniso(1);
  mop.setAuthRotation(1);
  mop.setLockNo3d(dbmap->getNDim() <= 2);

  // Case when properties are defined: Goulard is switch off
  _modifyMopForAnam(model, mop);

  // Case when constraints involve sill(s)
  if (constraints != nullptr && constraints->isDefinedForSill() && mop.getFlagGoulard())
  {
    mop.setFlagGoulard(false);
    if (modify_constraints_on_sill(*constraints)) return (1);
  }

  // Return an error if Goulard is not used in multivariate case 
  if (model->getNVar() > 1 && !mop.getFlagGoulard())
  {
    messerr("In Multivariate case, Goulard option is mandatory");
    messerr("It seems that it has been switched OFF. This is an error");
    return 1;
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Define the options of the structure Opt_Vario
 **
 ** \return Error return code
 **
 ** \param[in]  vario       Vario structure containing the exp. variogram
 ** \param[in]  model       ModelGeneric structure
 ** \param[in]  constraints Constraints structure
 **
 ** \param[out]  mop        ModelOptimParam structure
 **
 *****************************************************************************/
static int _modifyMopForVario(const Vario* vario,
                              ModelGeneric* model,
                              Constraints* constraints,
                              ModelOptimParam& mop)
{
  int ndim = model->getNDim();
  int ndir = vario->getNDir();
  int n_2d = 0;
  int n_3d = 0;

  /* 2-D case */
  if (ndim == 2)
  {
    n_2d = ndir;
    n_3d = 0;
  }

  /* 3-D case */
  if (ndim == 3)
  {
    for (int idir = 0; idir < ndir; idir++)
    {
      if (isZero(vario->getCodir(idir, 2)))
        n_2d++;
      else
        n_3d++;
    }
    mop.setLockNo3d(n_3d <= 0);
    mop.setLockIso2d(n_2d <= 0);
  }

  /* Clever setting of options */

  if (ndir <= ndim) mop.setAuthRotation(0);
  if (ndir <= 1 || ndim <= 1) mop.setAuthAniso(0);
  if (ndir <= 1 || ndim <= 1) mop.setAuthRotation(0);

  if (ndim == 3 && n_3d <= 0) mop.setLockNo3d(1);
  if (n_2d <= 1) mop.setLockIso2d(1);
  if (mop.getLockIso2d()) mop.setAuthRotation(0);
  if (mop.getLockNo3d()) mop.setLockRot2d(1);

  /* Consequences of no anisotropy */

  if (!mop.getAuthAniso())
  {
    mop.setAuthRotation(0);
    mop.setLockSamerot(0);
    mop.setLockRot2d(0);
    mop.setLockNo3d(0);
    mop.setLockIso2d(0);
  }

  /* Case when properties are defined: Goulard is switch off */

  _modifyMopForAnam(model, mop);

  /* Case when constraints involve sill(s) */

  if (constraints != nullptr && constraints->isDefinedForSill() && mop.getFlagGoulard())
  {
    if (modify_constraints_on_sill(*constraints)) return (1);
    mop.setFlagGoulard(false);
  }

  /* Return an error if Goulard is not used in multivariate case */

  if (model->getNVar() > 1 && !mop.getFlagGoulard())
  {
    messerr("In Multivariate case, Goulard option is mandatory");
    messerr("It seems that it has been switched OFF. This is an error");
    return 1;
  }
  return 0;
}

static void _modifyOneParam(const EConsType& cas,
                            ParamInfo* param,
                            double value)
{
  if (cas == EConsType::EQUAL)
  {
    param->setValue(value);
    param->setFixed(true);
  }
  else if (cas == EConsType::LOWER)
  {
    param->increaseMin(value);
  }
  else
  {
    param->decreaseMax(value);
  }
}

static int _modifyModelForConstraints(Constraints* constraints,
                                      ModelGeneric* model)
{
  // Check the constraints
  if (constraints == nullptr) return 0;
  int ncons = constraints->getNConsItem();
  if (ncons <= 0) return 0;

  // Check the ModelCovList
  ModelCovList* mcv = dynamic_cast<ModelCovList*>(model);
  if (mcv == nullptr) return 1;

  int ndim = model->getNDim();
  int nvar = model->getNVar();
  for (int i = 0; i < ncons; i++)
  {
    const ConsItem* consitem = constraints->getConsItems(i);
    const EConsElem type     = consitem->getType();
    const EConsType cas      = consitem->getIcase();
    double value             = consitem->getValue();
    int igrf                 = consitem->getIGrf();
    int icov                 = consitem->getICov();
    int iv1                  = consitem->getIV1();
    int iv2                  = consitem->getIV2();

    CovBase* covbase = mcv->getCovBase(icov);
    CovAniso* covaniso = dynamic_cast<CovAniso*>(covbase);
    ParamInfo* param = nullptr;

    if (igrf != 0)
    {
      messerr("Setting constraint for IGRF(%d) is not possible", igrf);
      return 1;
    }
    if (type == EConsElem::ANGLE)
    {
      if (iv1 < 0 || iv1 >= ndim)
      {
        messerr("Setting Angle(%d) not possible as ndim=%d", iv1, ndim);
        return 1;
      }
      param = &covaniso->getCorAniso()->getParamInfoAngle(iv1);
    }
    else if (type == EConsElem::RANGE)
    {
      if (iv1 < 0 || iv1 >= ndim)
      {
        messerr("Setting Range(%d) not possible as ndim=%d", iv1, ndim);
        return 1;
      }
      // Convert range into scale (using the current value for 'param')
      double scadef = covaniso->getCorAniso()->getScadef();
      value *= scadef;
      param = &covaniso->getCorAniso()->getParamInfoScale(iv1);
    }
    else if (type == EConsElem::SCALE)
    {
      if (iv1 < 0 || iv1 >= ndim)
      {
        messerr("Setting Scale(%d) not possible as ndim=%d", iv1, ndim);
        return 1;
      }
      param = &covaniso->getCorAniso()->getParamInfoScale(iv1);
    }
    else if (type == EConsElem::SILL)
    {
      if (nvar > 1)
      {
        messerr("Setting Sill is impossible as nvar = %d", nvar);
        return 1;
      }
      if (iv1 != 0 && iv2 != 0)
      {
        messerr("Setting Sill(%d,%d) not possible as nvar = %d", iv1, iv2, nvar); 
        return 1;
      }
      param = &covbase->getParamInfoCholSills(iv1, iv2);
    }
    else
    {
      messerr("Unknown Parameter");
    }
    _modifyOneParam(cas, param, value);
  }
  return 0;
}

static int _modifyModelForMop(const ModelOptimParam& mop,
                              ModelGeneric* model)
{
  ModelCovList* mcv = dynamic_cast<ModelCovList*>(model);
  if (mcv == nullptr) return 0;

  // Note that the creation of the dedicated area in case of Goulard will
  // be performed by the calling function as this area depends on the
  // calling function (Vario or Vmap)

  // Loop on the structures
  int nvar = model->getNVar();
  for (int icov = 0, ncov = mcv->getNCov(); icov < ncov; icov++)
  {
    CovBase* covbase = mcv->getCovBase(icov);
    if (covbase == nullptr) continue;


    // Set the Goulard constraints
    if (mop.getFlagGoulard())
    {
      // Fix the sills
      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar < nvar; jvar++)
        {
          ParamInfo paraminfo = covbase->getParamInfoCholSills(ivar, jvar);
          paraminfo.setFixed(true);
        }
    }

    CovAniso* covaniso = dynamic_cast<CovAniso*>(covbase);
    if (covaniso == nullptr) continue;
    CorAniso* coraniso = covaniso->getCorAniso();
    if (coraniso == nullptr) continue;

    // Anisotropy
    if (mop.getAuthAniso())
    {

      // Anisotropy is authorized
      if (mop.getLockIso2d())
        coraniso->setOptimLockIso2d(true);

      if (mop.getLockNo3d())
      {
        std::vector<ParamInfo>& params = coraniso->getParamInfoScales();
        for (int ipar = 2, npar = (int) params.size(); ipar < npar; ipar++)
        {
          params[ipar].setValue(0.);
          params[ipar].setFixed(true);
        }
      }

      // Anisotropy rotation
      if (mop.getAuthRotation())
      {
        // Anisotropy rotation is authorized

        if (mop.getLockRot2d())
        {
          std::vector<ParamInfo>& params = coraniso->getParamInfoAngles();
          for (int ipar = 1, npar = (int)params.size(); ipar < npar; ipar++)
          {
            params[ipar].setValue(0.);
            params[ipar].setFixed(true);
          }
        }
      }
      else
      {
        std::vector<ParamInfo>& params = coraniso->getParamInfoAngles();
        for (int ipar = 0, npar = (int) params.size(); ipar < npar; ipar++)
        {
          params[ipar].setFixed(true);
        }
      }
    }
    else 
    {

      // Anisotropy forbidden
      if (coraniso == nullptr) continue;
      coraniso->setOptimNoAniso(true);
    }
  }

  return 0;
}

AModelOptim* AModelOptimFactory::create(ModelGeneric* model,
                                        const Db* db,
                                        Vario* vario,
                                        const DbGrid* dbmap,
                                        Constraints* constraints,
                                        const ModelOptimParam& mop,
                                        int nb_neighVecchia)
{
  ModelOptimParam mopLocal = mop;

  // Fitting from LogLikelihood
  if (db != nullptr)
  {
    if ((int)model->getNDim() != db->getNDim()) return nullptr;
    if (nb_neighVecchia != ITEST) return Vecchia::createForOptim(model, db, nb_neighVecchia);
    return Likelihood::createForOptim(model, db);
  }

  // Fitting from a Variogram Map
  if (dbmap != nullptr)
  {
    if ((int)model->getNDim() != dbmap->getNDim()) return nullptr;
    if (_modifyModelForConstraints(constraints, model)) return nullptr;
    if (_modifyMopForVMap(dbmap, model, constraints, mopLocal)) return nullptr;
    if (_modifyModelForMop(mop, model)) return nullptr;
    return ModelOptimVMap::createForOptim(model, dbmap, constraints, mopLocal);
  }

  // Fitting from an experimental Variogram
  if (vario != nullptr)
  {
    if ((int)model->getNDim() != vario->getNDim()) return nullptr;
    if (_modifyMopForVario(vario, model, constraints, mopLocal)) return nullptr;
    if (_modifyModelForConstraints(constraints, model)) return nullptr;
    if (_modifyModelForMop(mop, model)) return nullptr;
    return ModelOptimVario::createForOptim(model, vario, constraints, mopLocal);
  }
  return nullptr;
}
