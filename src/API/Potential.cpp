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
#include "API/Potential.hpp"
#include "Basic/Law.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/Utilities.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovGradientAnalytic.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbHelper.hpp"
#include "Drifts/DriftList.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Model/Model.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "geoslib_old_f.h"

#include <cmath>
#include <cstring>

#define DRF(il)           (TAB_DRF[il])
#define MAT(i, j)         (mat[(i) * _nequa + (j)])
#define B(isol, i)        (b[(isol) * _number + (i)])
#define POTVAL(isimu, il) (potval[(isimu) * _nlayers + (il)])
#define POTSIM(isimu, il) (potsim[(isimu) * _nlayers + (il)])

namespace gstlrn
{
static Id TAB_DRF[9];

Potential::Potential(Db* dbiso,
                     Db* dbgrd,
                     Db* dbtgt,
                     ModelGeneric* model,
                     double nugget_grd,
                     double nugget_tgt)
  : _dbiso(dbiso)
  , _dbgrd(dbgrd)
  , _dbtgt(dbtgt)
  , _dbExt(nullptr)
  , _model(model)
  , _modelExt(nullptr)
  , _ndim(0)
  , _niso(0)
  , _nlayers(0)
  , _ngrd(0)
  , _ntgt(0)
  , _next(0)
  , _nequa(0)
  , _order(0)
  , _nring(0)
  , _nfull(0)
  , _sizeIso(0)
  , _sizeGrd(0)
  , _sizeTgt(0)
  , _sizeDrf(0)
  , _sizeExt(0)
  , _startIso(0)
  , _startGrd(0)
  , _startTgt(0)
  , _startDrf(0)
  , _startExt(0)
  , _nugget_grd(nugget_grd)
  , _nugget_tgt(nugget_tgt)
  , _rangeExt(0)
  , _nbPerLayer()
  , _ptrPerLayer()
  , _rankIso()
  , _rankGrd()
  , _rankTgt()
  , _optionPart(0)
  , _flagPot(false)
  , _flagGrad(false)
  , _flagTrans(false)
  , _verbose(false)
  , _indg()
  , _indg0()
  , _dataExt()
  , _wgtExt()
  , _p1()
  , _p2()
{
}

Potential::~Potential()
{
  delete _dbExt;
  delete _modelExt;
}

bool Potential::_isEnvironmentValid(DbGrid* dbout, Id nring)
{
  if (_ndim > 3)
  {
    messerr("The input Db must be defined in Space with dimension < 3");
    return false;
  }
  if (_dbiso == nullptr)
  {
    messerr("The IsoPotential Data Db is compulsory");
    return false;
  }
  if (_dbgrd != nullptr && _dbgrd->getNDim() != _ndim)
  {
    messerr("The Gradient and Data Db must share the same space dimension");
    return false;
  }
  if (_dbtgt != nullptr && _dbtgt->getNDim() != _ndim)
  {
    messerr("The Tangent and Data Db must share the same space dimension");
    return false;
  }
  if (static_cast<Id>(_model->getNDim()) != _ndim)
  {
    messerr("The Model and Data Db must have the same space dimension");
    return false;
  }
  if (dbout != NULL && dbout->getNDim() != _ndim)
  {
    messerr("The Db files 'dbin' and 'dbout' should have the same dimension");
    return false;
  }

  if (!_dbiso->hasLocator(ELoc::LAYER))
  {
    messerr("The input Db must contain a LAYER locator");
    return false;
  }
  const auto* covpot = dynamic_cast<const CovGradientAnalytic*>(_model->getCov());
  if (covpot == nullptr)
  {
    messerr("The Model is invalid for Potential calculations");
    messerr("It may only contain CovGradientAnalytic(s)");
    return false;
  }
  if (_model->getNVar() != 1)
  {
    messerr("The Model must be monovariate");
    return false;
  }
  Id next = _model->getNExtDrift();
  if (dbout != NULL && next != dbout->getNLoc(ELoc::F))
  {
    messerr("Inconsistency for External Drift between Model and Dbout");
    return false;
  }
  if (dbout == NULL && next > 0)
  {
    messerr("Usage of External drift is forbidden without Output Grid");
  }
  if (next > 0)
  {
    if (next > 1)
    {
      messerr("This application cannot deal with more than 1 External Drift");
      messerr("Check your output file");
      return false;
    }
    if (!dbout->isGrid())
    {
      messerr("The External Drift requires an Output Grid File");
      return false;
    }
    _rangeExt = 3. * MAX(dbout->getDX(0), dbout->getDX(1));
    _nring    = nring;

    if (_extdriftEstablish(dbout)) return (1);
  }

  return true;
}

/****************************************************************************/
/*!
 **  Manage the class members
 **
 ** \param[in]      flag_pot   True if the potential must be calculated
 ** \param[in]      flag_grad  True if the gradients must be calculated
 ** \param[in]      flag_trans True if the estimation result must be translated
 **                            into layer number
 ** \param[in]      option_part   Option to exhibit only a part of estimation:
 ** \li                        0 : the whole estimation
 ** \li                        1 : the gradient contribution only
 ** \li                        2 : the tangent contribution only
 ** \li                        3 : the isovalues contribution only
 ** \li                        4 : the drift contribution only
 ** \li                        5 : the external drift contribution only
 **
 *****************************************************************************/
void Potential::_environmentManage(bool flag_pot,
                                   bool flag_grad,
                                   bool flag_trans,
                                   Id option_part)
{
  if (option_part) flag_trans = false;

  _niso        = 0;
  _nlayers     = 0;
  _ngrd        = 0;
  _ntgt        = 0;
  _next        = 0;
  _nequa       = 0;
  _order       = 0;
  _sizeIso     = 0;
  _sizeGrd     = 0;
  _sizeTgt     = 0;
  _sizeDrf     = 0;
  _sizeExt     = 0;
  _startIso    = 0;
  _startGrd    = 0;
  _startTgt    = 0;
  _startDrf    = 0;
  _startExt    = 0;
  _nbPerLayer  = VectorInt();
  _ptrPerLayer = VectorInt();
  _rankIso     = VectorInt();
  _rankGrd     = VectorInt();
  _rankTgt     = VectorInt();
  _flagPot     = flag_pot;
  _flagGrad    = flag_grad;
  _flagTrans   = flag_trans;
  _optionPart  = option_part;
  _ndim        = _dbiso->getNDim();

  _p1 = SpacePoint(_ndim);
  _p2 = SpacePoint(_ndim);
}

Id Potential::_updateIsopot()
{
  if (_dbiso == nullptr) return (0);
  Id nech    = _dbiso->getNSample();
  Id nlayers = 0;
  Id niso    = 0;
  VectorInt laycnt;
  VectorInt layval;

  // Count the number of different iso-potential values

  for (Id iech = 0; iech < nech; iech++)
  {
    if (!_dbiso->isActive(iech)) continue;
    double value = _dbiso->getFromLocator(ELoc::LAYER, iech);
    if (FFFF(value)) continue;
    Id ival = static_cast<Id>(value);

    // Look for an already registered layer value

    Id found = -1;
    for (Id i = 0; i < nlayers && found < 0; i++)
      if (ival == layval[i]) found = i;

    if (found < 0)
    {
      layval.push_back(ival);
      laycnt.push_back(1);
      nlayers++;
    }
    else
    {
      laycnt[found]++;
    }
  }

  // Eliminate layers with not enough samples

  niso = 0;
  Id j = 0;
  for (Id i = 0; i < nlayers; i++)
  {
    if (laycnt[i] < 2) continue;
    layval[j] = layval[i];
    laycnt[j] = laycnt[i];
    niso += laycnt[i];
    j++;
  }
  _nlayers = nlayers = j;
  _niso              = niso;
  _sizeIso           = niso - nlayers;

  // Core allocation

  _nbPerLayer.resize(nlayers);
  _ptrPerLayer.resize(nlayers);
  _rankIso.resize(niso);
  if (_rankIso.empty()) return 1;

  // Set the final length and pointers

  Id ipos = 0;
  for (Id i = 0; i < nlayers; i++)
  {
    _nbPerLayer[i]  = laycnt[i];
    _ptrPerLayer[i] = ipos;
    ipos += _nbPerLayer[i];
  }

  // Sort the samples per iso-potential value

  Id ecr = 0;
  for (Id i = 0; i < nlayers; i++)
  {
    for (Id iech = 0; iech < nech; iech++)
    {
      if (!_dbiso->isActive(iech)) continue;
      double value = _dbiso->getFromLocator(ELoc::LAYER, iech);
      if (FFFF(value)) continue;
      Id ival = static_cast<Id>(value);
      if (ival != layval[i]) continue;
      _rankIso[ecr++] = iech;
    }
  }

  // Reading failure

  if (niso < 1 || nlayers < 1)
  {
    messerr("The number of iso-potential informations cannot be null");
    return 1;
  }

  return 0;
}

Id Potential::_updateGradient()
{
  if (_dbgrd == nullptr) return (0);
  Id nech = _dbgrd->getNSample();
  Id ngrd = 0;
  _rankGrd.resize(nech);

  // Loop on the gradients

  for (Id iech = 0; iech < nech; iech++)
  {
    if (!_dbgrd->isActive(iech)) continue;
    Id found = 0;
    for (Id idim = 0; idim < _ndim && found == 0; idim++)
      if (FFFF(_dbgrd->getLocVariable(ELoc::G, iech, idim))) found = 1;
    if (found) continue;
    set_IAD_GRD(ngrd++, iech);
  }

  // Core reallocation

  _rankGrd.resize(ngrd);
  _ngrd    = ngrd;
  _sizeGrd = ngrd * _ndim;

  if (ngrd < 1)
  {
    messerr("The number of gradient informations cannot be null");
    return (1);
  }

  return (0);
}

Id Potential::_updateTangent()
{
  if (_dbtgt == nullptr) return (0);
  Id nech = _dbtgt->getNSample();
  Id ntgt = 0;
  _rankTgt.resize(nech);

  // Loop on the tangents

  for (Id iech = 0; iech < nech; iech++)
  {
    if (!_dbtgt->isActive(iech)) continue;
    Id found = 0;
    for (Id idim = 0; idim < _ndim && found == 0; idim++)
      if (FFFF(_dbtgt->getLocVariable(ELoc::TGTE, iech, idim))) found = 1;
    if (found) continue;
    set_IAD_TGT(ntgt++, iech);
  }

  // Core reallocation

  _rankTgt.resize(ntgt);
  _ntgt    = ntgt;
  _sizeTgt = ntgt;

  return (0);
}

Id Potential::_updateModel()
{
  Id nbfl = _model->getNDrift();
  if (_model->isDriftDefined(VectorInt(), 0)) nbfl--;
  _order   = _model->getDriftMaxIRFOrder();
  _sizeDrf = nbfl;
  _next = _sizeExt = _model->getNExtDrift();

  return (0);
}

Id Potential::_updateFinal()
{
  // Compute the starting addresses

  Id pos    = 0;
  _startGrd = pos;
  pos += _sizeGrd;
  _startTgt = pos;
  pos += _sizeTgt;
  _startIso = pos;
  pos += _sizeIso;
  _startDrf = pos;
  pos += _sizeDrf;
  _startExt = pos;

  // Compute the number of equations in the CoKriging System

  _nequa = (_sizeIso + _sizeGrd + _sizeTgt + _sizeDrf + _sizeExt);

  // Define the addresses for the drift functions

  pos = _startDrf;
  for (Id i = 0; i < 9; i++)
    TAB_DRF[i] = -1;

  if (_model->isDriftDefined(VectorInt {1})) TAB_DRF[0] = pos++;
  if (_model->isDriftDefined(VectorInt {0, 1})) TAB_DRF[1] = pos++;
  if (_model->isDriftDefined(VectorInt {0, 0, 1})) TAB_DRF[2] = pos++;
  if (_model->isDriftDefined(VectorInt {2})) TAB_DRF[3] = pos++;
  if (_model->isDriftDefined(VectorInt {0, 2})) TAB_DRF[4] = pos++;
  if (_model->isDriftDefined(VectorInt {0, 0, 2})) TAB_DRF[5] = pos++;
  if (_model->isDriftDefined(VectorInt {1, 1})) TAB_DRF[6] = pos++;
  if (_model->isDriftDefined(VectorInt {1, 0, 1})) TAB_DRF[7] = pos++;
  if (_model->isDriftDefined(VectorInt {0, 1, 1})) TAB_DRF[8] = pos++;

  /* Optional output */

  if (_verbose)
  {
    mestitle(0, "Environment summary");
    message("Space dimension         = %d\n", _ndim);
    message("Number of Iso-Potential = %d\n", _nlayers);
    message("Number of Gradients     = %d\n", _ngrd);
    message("Number of Tangents      = %d\n", _ntgt);
    message("Number of Isovalues     = %d\n", _niso);
    message("Order of the drift      = %d\n", _order);
    message("Number of Drifts        = %d\n", _sizeDrf);
    message("Number of Ext. Drifts   = %d\n", _sizeExt);
    message("Number of Equations     = %d\n", _nequa);
  }

  return (0);
}

void Potential::_saveResultData(Db* db,
                                Id nvar,
                                double value,
                                const ELoc& loctype_pot,
                                const ELoc& loctype_grad,
                                VectorInt& uid_pot,
                                VectorInt& uid_grad) const
{
  uid_pot.clear();
  uid_grad.clear();
  if (db == nullptr) return;

  if (_flagPot)
  {
    Id uid = db->addColumnsByConstant(nvar, value, "Potential", loctype_pot);
    uid_pot.push_back(uid);
  }
  if (_flagGrad)
  {
    Id uid = db->addColumnsByConstant(nvar * _ndim, value, "Gradients", loctype_grad);
    for (Id idim = 0; idim < _ndim; idim++)
      uid_grad.push_back(uid + idim);
  }
}

/****************************************************************************/
/*!
 **  Calculate the covariance and the derivatives
 **
 ** \param[in] model     Model to be used for the calculations
 ** \param[in] flag_grad True if the Gradients must be calculated
 ** \param[in] x1,y1,z1  First point
 ** \param[in] x2,y2,z2  Second point
 **
 ** \param[out] covar    Covariance of potential
 ** \param[out] covGp    Covariance between potential and gradient
 ** \param[out] covGG    Covariance between gradient and gradient
 **
 *****************************************************************************/
void Potential::_calculateCovs(ModelGeneric* model,
                               bool flag_grad,
                               double x1,
                               double y1,
                               double z1,
                               double x2,
                               double y2,
                               double z2,
                               double& covar,
                               VectorDouble& covGp,
                               VectorDouble& covGG) const
{
  auto* covpot = dynamic_cast<CovGradientAnalytic*>(model->_getCovModify());
  if (covpot == nullptr) return;

  // Define a new pair of points
  if (_ndim >= 1) _p1.setCoord(0, x1);
  if (_ndim >= 2) _p1.setCoord(1, y1);
  if (_ndim >= 3) _p1.setCoord(2, z1);
  if (_ndim >= 1) _p2.setCoord(0, x2);
  if (_ndim >= 2) _p2.setCoord(1, y2);
  if (_ndim >= 3) _p2.setCoord(2, z2);

  covpot->setFlagCalculateGG(flag_grad);

  covar = covpot->evalCov(_p1, _p2, 0, 0);
  for (Id idim = 0, ecr = 0; idim < 3; idim++)
  {
    covGp[idim] = covpot->evalCov(_p1, _p2, 1 + idim, 0);
    if (flag_grad)
    {
      for (Id jdim = 0; jdim < 3; jdim++, ecr++)
      {
        covGG[ecr] = covpot->evalCov(_p1, _p2, 1 + idim, 1 + jdim);
      }
    }
  }
}

/****************************************************************************/
/*!
 **  Set and element in the Kriging L.H.S. matrix
 **
 ** \param[in] lhs      Matrix to be filled
 ** \param[in] i        Row number
 ** \param[in] j        Column number
 ** \param[in] value    Value to be assigned to this cell
 **
 *****************************************************************************/
void Potential::_setLHS(MatrixSymmetric& lhs, Id i, Id j, double value)
{
  if (i < 0 || j < 0) return;
  lhs.setValue(i, j, value);
}

/****************************************************************************/
/*!
 **  Get one element from the Kriging L.H.S. matrix
 **
 ** \return The returned value
 **
 ** \param[in] lhs      Matrix to be filled
 ** \param[in] i        Row number
 ** \param[in] j        Column number
 **
 *****************************************************************************/
double Potential::_getLHS(MatrixSymmetric& lhs, Id i, Id j)
{
  if (i < 0 || j < 0) return (0.);
  return lhs.getValue(i, j);
}

/****************************************************************************/
/*!
 **  Set and element in the Kriging R.H.S. vector
 **
 ** \param[in] rhs      Vector to be filled
 ** \param[in] i        Row number
 ** \param[in] isol     Column number
 ** \param[in] value    Value to be assigned to this cell
 **
 *****************************************************************************/
void Potential::_setRHS(MatrixDense& rhs, Id i, Id isol, double value)
{
  if (i < 0 || isol < 0) return;
  rhs.setValue(i, isol, value);
}

/****************************************************************************/
/*!
 **  Calculate the inner product of two vectors
 **
 **  \param[in]   ux     : First coordinate of the first vector
 **  \param[in]   uy     : Second coordinate of the first vector
 **  \param[in]   uz     : Third coordinate of the first vector
 **  \param[in]   vx     : First coordinate of the second vector
 **  \param[in]   vy     : Second coordinate of the second vector
 **  \param[in]   vz     : Third coordinate of the second vector
 **
 *****************************************************************************/
double Potential::_setMatUV(double ux,
                            double uy,
                            double uz,
                            double vx,
                            double vy,
                            double vz) const
{
  double prod = 0.;
  if (_ndim >= 1 && !FFFF(ux) && !FFFF(vx)) prod += ux * vx;
  if (_ndim >= 2 && !FFFF(uy) && !FFFF(vy)) prod += uy * vy;
  if (_ndim >= 3 && !FFFF(uz) && !FFFF(vz)) prod += uz * vz;
  return prod;
}

/****************************************************************************/
/*!
 **  Calculate the norm product of two vectors by a matrix
 **
 **  \param[in]   a      : Matrix
 **  \param[in]   ux     : First coordinate of the first vector
 **  \param[in]   uy     : Second coordinate of the first vector
 **  \param[in]   uz     : Third coordinate of the first vector
 **  \param[in]   vx     : First coordinate of the second vector
 **  \param[in]   vy     : Second coordinate of the second vector
 **  \param[in]   vz     : Third coordinate of the second vector
 **
 *****************************************************************************/
double Potential::_setMatUAV(const VectorDouble& a,
                             double ux,
                             double uy,
                             double uz,
                             double vx,
                             double vy,
                             double vz) const
{
  double prod = 0.;
  if (_ndim >= 1 && !FFFF(ux) && !FFFF(vx)) prod += ux * vx * a[0];
  if (_ndim >= 2 && !FFFF(ux) && !FFFF(vy)) prod += ux * vy * a[1];
  if (_ndim >= 3 && !FFFF(ux) && !FFFF(vz)) prod += ux * vz * a[2];
  if (_ndim >= 2 && !FFFF(uy) && !FFFF(vx)) prod += uy * vx * a[3];
  if (_ndim >= 2 && !FFFF(uy) && !FFFF(vy)) prod += uy * vy * a[4];
  if (_ndim >= 3 && !FFFF(uy) && !FFFF(vz)) prod += uy * vz * a[5];
  if (_ndim >= 3 && !FFFF(uz) && !FFFF(vx)) prod += uz * vx * a[6];
  if (_ndim >= 3 && !FFFF(uz) && !FFFF(vy)) prod += uz * vy * a[7];
  if (_ndim >= 3 && !FFFF(uz) && !FFFF(vz)) prod += uz * vz * a[8];
  return (prod);
}

/****************************************************************************/
/*!
 **  Establish the cokriging system
 **
 ** \return  LHS Matrix
 **
 ** \param[in]  dbout         Output Db (for external drift)
 ** \param[out] lhs           Cokriging LHS matrix
 **
 ** \remark   Organization of the cokriging system
 ** \remark
 ** \remark   |   A11 =        |  A12  =        | A13  =              | F1 =  |
 ** \remark   | <Gu(i),Gv(i')> | <Gu(i) ,(T,G)> | <Gu(i),Pl(j)-Pl(0)> | G(Fl) |
 ** \remark   -----------------------------------------------------------------
 ** \remark   |   A21  =       |  A22  =        | A23  =              | F2 =  |
 ** \remark   | <(T,G),Gv(i')  | <(T,G) ,(T,G)> | <(T,G),Pl(j)-Pl(0)> | (T,Fl)|
 ** \remark   -----------------------------------------------------------------
 ** \remark   |   A31  =       |  A32  =        | A33 =               | F3 =  |
 ** \remark   | <Pl(i)-Pl(0),G>| <Pli-Pl0,(T,g)>| <Pli-Pl0,Pl'i-Pl'0> |Fli-Fl0|
 ** \remark   -----------------------------------------------------------------
 ** \remark   | F1t            | F2t            | F3t                 |  0    |
 ** \remark
 ** \remark   The matrix A11 is subdivided as follows:
 ** \remark
 ** \remark         |  <Gx,Gx> |  <Gx,Gy>  |  <Gx,Gz> |
 ** \remark         |          |           |          |
 ** \remark   A11 = |  <Gy,Gx> |  <Gy,Gy>  |  ....    |
 ** \remark         |          |           |          |
 ** \remark         |  <Gz,Gx> |  ...      |  ..      |
 ** \remark
 ** \remark   each one of the 9 blocks has dimension = the number of gradients
 **
 *****************************************************************************/
Id Potential::_buildLHS(Db* dbout, MatrixSymmetric& lhs)
{
  lhs.resize(_nequa, _nequa);
  lhs.fill(0.);

  VectorDouble covGp(3, 0.);
  VectorDouble covGG(9, 0.);
  VectorDouble cov2Gp(3, 0.);
  VectorDouble cov2GG(9, 0.);
  VectorDouble center(3, 0.);
  VectorDouble extgrd(3, 0.);
  double covar   = 0.;
  double covar1  = 0.;
  double covar2  = 0.;
  double covar3  = 0.;
  double covar4  = 0.;
  double extval  = 0.;
  double extval1 = 0.;
  double extval2 = 0.;

  // Blank out the cokriging matrix

  /******************************/
  /* PART RELATIVE TO GRADIENTS */
  /******************************/

  for (Id ig = 0; ig < _ngrd; ig++)
  {
    for (Id jg = 0; jg < ig; jg++)
    {
      _calculateCovs(_model, true,
                     GRD_COO(ig, 0), GRD_COO(ig, 1), GRD_COO(ig, 2),
                     GRD_COO(jg, 0), GRD_COO(jg, 1), GRD_COO(jg, 2),
                     covar, covGp, covGG);

      _setLHS(lhs, GRX(ig), GRX(jg), covGG[0]);
      _setLHS(lhs, GRX(ig), GRY(jg), covGG[1]);
      _setLHS(lhs, GRX(ig), GRZ(jg), covGG[2]);
      _setLHS(lhs, GRY(ig), GRX(jg), covGG[3]);
      _setLHS(lhs, GRY(ig), GRY(jg), covGG[4]);
      _setLHS(lhs, GRY(ig), GRZ(jg), covGG[5]);
      _setLHS(lhs, GRZ(ig), GRX(jg), covGG[6]);
      _setLHS(lhs, GRZ(ig), GRY(jg), covGG[7]);
      _setLHS(lhs, GRZ(ig), GRZ(jg), covGG[8]);
    }
    _calculateCovs(_model, true,
                   GRD_COO(ig, 0), GRD_COO(ig, 1), GRD_COO(ig, 2),
                   GRD_COO(ig, 0), GRD_COO(ig, 1), GRD_COO(ig, 2),
                   covar, covGp, covGG);
    _setLHS(lhs, GRX(ig), GRX(ig), covGG[0] + _nugget_grd);
    _setLHS(lhs, GRX(ig), GRY(ig), covGG[1]);
    _setLHS(lhs, GRX(ig), GRZ(ig), covGG[2]);
    _setLHS(lhs, GRY(ig), GRX(ig), covGG[3]);
    _setLHS(lhs, GRY(ig), GRY(ig), covGG[4] + _nugget_grd);
    _setLHS(lhs, GRY(ig), GRZ(ig), covGG[5]);
    _setLHS(lhs, GRZ(ig), GRX(ig), covGG[6]);
    _setLHS(lhs, GRZ(ig), GRY(ig), covGG[7]);
    _setLHS(lhs, GRZ(ig), GRZ(ig), covGG[8] + _nugget_grd);
  }

  /*****************************/
  /* PART RELATIVE TO TANGENTS */
  /*****************************/

  for (Id it = 0; it < _ntgt; it++)
  {

    /* block tangents-gradients */

    for (Id ig = 0; ig < _ngrd; ig++)
    {
      _calculateCovs(_model, true,
                     TGT_COO(it, 0), TGT_COO(it, 1), TGT_COO(it, 2),
                     GRD_COO(ig, 0), GRD_COO(ig, 1), GRD_COO(ig, 2),
                     covar, covGp, covGG);

      _setLHS(lhs, TGT(it), GRX(ig),
              _setMatUV(TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                        covGG[0], covGG[1], covGG[2]));
      _setLHS(lhs, TGT(it), GRY(ig),
              _setMatUV(TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                        covGG[3], covGG[4], covGG[5]));
      _setLHS(lhs, TGT(it), GRZ(ig),
              _setMatUV(TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                        covGG[6], covGG[7], covGG[8]));
    }

    /* block diagonal tangents */

    for (Id jt = 0; jt < it; jt++)
    {
      _calculateCovs(_model, true,
                     TGT_COO(it, 0), TGT_COO(it, 1), TGT_COO(it, 2),
                     TGT_COO(jt, 0), TGT_COO(jt, 1), TGT_COO(jt, 2),
                     covar, covGp, covGG);

      _setLHS(lhs, TGT(it), TGT(jt),
              _setMatUAV(covGG,
                         TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                         TGT_VAL(jt, 0), TGT_VAL(jt, 1), TGT_VAL(jt, 2)));
    }
    _calculateCovs(_model, true,
                   TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                   TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                   covar, covGp, covGG);
    _setLHS(lhs, TGT(it), TGT(it),
            _setMatUAV(covGG,
                       TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                       TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2)) +
              _nugget_tgt);
  }

  /***********************************/
  /* PART RELATIVE TO ISO-POTENTIALS */
  /***********************************/

  for (Id ic1 = 0; ic1 < _nlayers; ic1++)
  {
    for (Id j1 = 1; j1 < _nbPerLayer[ic1]; j1++)
    {

      /* Interactions isopotentials - gradients */

      for (Id ig = 0; ig < _ngrd; ig++)
      {
        _calculateCovs(_model, true,
                       GRD_COO(ig, 0), GRD_COO(ig, 1), GRD_COO(ig, 2),
                       ISO_COO(ic1, 0, 0), ISO_COO(ic1, 0, 1), ISO_COO(ic1, 0, 2),
                       covar, covGp, covGG);
        _calculateCovs(_model, true,
                       GRD_COO(ig, 0), GRD_COO(ig, 1), GRD_COO(ig, 2),
                       ISO_COO(ic1, j1, 0), ISO_COO(ic1, j1, 1), ISO_COO(ic1, j1, 2),
                       covar, cov2Gp, cov2GG);
        _setLHS(lhs, ISC(ic1, j1), GRX(ig), cov2Gp[0] - covGp[0]);
        _setLHS(lhs, ISC(ic1, j1), GRY(ig), cov2Gp[1] - covGp[1]);
        _setLHS(lhs, ISC(ic1, j1), GRZ(ig), cov2Gp[2] - covGp[2]);
      }

      /* Interactions increments-tangentes */

      for (Id it = 0; it < _ntgt; it++)
      {
        _calculateCovs(_model, true,
                       TGT_COO(it, 0), TGT_COO(it, 1), TGT_COO(it, 2),
                       ISO_COO(ic1, 0, 0), ISO_COO(ic1, 0, 1), ISO_COO(ic1, 0, 2),
                       covar, covGp, covGG);
        _calculateCovs(_model, true,
                       TGT_COO(it, 0), TGT_COO(it, 1), TGT_COO(it, 2),
                       ISO_COO(ic1, j1, 0), ISO_COO(ic1, j1, 1), ISO_COO(ic1, j1, 2),
                       covar, cov2Gp, cov2GG);
        _setLHS(lhs, ISC(ic1, j1), TGT(it),
                _setMatUV(TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                          cov2Gp[0] - covGp[0], cov2Gp[1] - covGp[1], cov2Gp[2] - covGp[2]));
      }

      /* Block diagonal for iso-potentials */

      for (Id ic2 = 0; ic2 <= ic1; ic2++)
      {
        for (Id j2 = 1; j2 < _nbPerLayer[ic2]; j2++)
        {
          _calculateCovs(_model, false,
                         ISO_COO(ic2, j2, 0), ISO_COO(ic2, j2, 1), ISO_COO(ic2, j2, 2),
                         ISO_COO(ic1, j1, 0), ISO_COO(ic1, j1, 1), ISO_COO(ic1, j1, 2),
                         covar1, covGp, covGG);
          _calculateCovs(_model, false,
                         ISO_COO(ic2, j2, 0), ISO_COO(ic2, j2, 1), ISO_COO(ic2, j2, 2),
                         ISO_COO(ic1, 00, 0), ISO_COO(ic1, 00, 1), ISO_COO(ic1, 00, 2),
                         covar2, covGp, covGG);
          _calculateCovs(_model, false,
                         ISO_COO(ic2, 00, 0), ISO_COO(ic2, 00, 1), ISO_COO(ic2, 00, 2),
                         ISO_COO(ic1, j1, 0), ISO_COO(ic1, j1, 1), ISO_COO(ic1, j1, 2),
                         covar3, covGp, covGG);
          _calculateCovs(_model, false,
                         ISO_COO(ic2, 00, 0), ISO_COO(ic2, 00, 1), ISO_COO(ic2, 00, 2),
                         ISO_COO(ic1, 00, 0), ISO_COO(ic1, 00, 1), ISO_COO(ic1, 00, 2),
                         covar4, covGp, covGG);
          _setLHS(lhs, ISC(ic1, j1), ISC(ic2, j2), covar1 - covar2 - covar3 + covar4);
        }
      }
    }
  }

  /****************************/
  /* Part linked to the drift */
  /****************************/

  /* Part relative to gradients */

  for (Id ig = 0; ig < _ngrd; ig++)
  {
    if (_order >= 1)
    {

      /* Function x, y and z */

      _setLHS(lhs, DRF(0), GRX(ig), 1.);
      _setLHS(lhs, DRF(0), GRY(ig), 0.);
      _setLHS(lhs, DRF(0), GRZ(ig), 0.);

      _setLHS(lhs, DRF(1), GRX(ig), 0.);
      _setLHS(lhs, DRF(1), GRY(ig), 1.);
      _setLHS(lhs, DRF(1), GRZ(ig), 0.);

      _setLHS(lhs, DRF(2), GRX(ig), 0.);
      _setLHS(lhs, DRF(2), GRY(ig), 0.);
      _setLHS(lhs, DRF(2), GRZ(ig), 1.);
    }

    if (_order >= 2)
    {

      /* Functions x^2, y^2 et z^2 */

      _setLHS(lhs, DRF(3), GRX(ig), 2. * GRD_COO(ig, 0));
      _setLHS(lhs, DRF(3), GRY(ig), 0.);
      _setLHS(lhs, DRF(3), GRZ(ig), 0.);

      _setLHS(lhs, DRF(4), GRX(ig), 0.);
      _setLHS(lhs, DRF(4), GRY(ig), 2. * GRD_COO(ig, 1));
      _setLHS(lhs, DRF(4), GRZ(ig), 0.);

      _setLHS(lhs, DRF(5), GRX(ig), 0.);
      _setLHS(lhs, DRF(5), GRY(ig), 0.);
      _setLHS(lhs, DRF(5), GRZ(ig), 2. * GRD_COO(ig, 2));

      /* Functions xy, xz, et yz */

      _setLHS(lhs, DRF(6), GRX(ig), GRD_COO(ig, 1));
      _setLHS(lhs, DRF(6), GRY(ig), GRD_COO(ig, 0));
      _setLHS(lhs, DRF(6), GRZ(ig), 0.);

      _setLHS(lhs, DRF(7), GRX(ig), GRD_COO(ig, 2));
      _setLHS(lhs, DRF(7), GRY(ig), 0.);
      _setLHS(lhs, DRF(7), GRZ(ig), GRD_COO(ig, 0));

      _setLHS(lhs, DRF(8), GRX(ig), 0.);
      _setLHS(lhs, DRF(8), GRY(ig), GRD_COO(ig, 2));
      _setLHS(lhs, DRF(8), GRZ(ig), GRD_COO(ig, 1));
    }

    /* External drift(s) */

    for (Id iext = 0; iext < _next; iext++)
    {
      if (_extdriftEval(GRD_COO(ig, 0), GRD_COO(ig, 1), GRD_COO(ig, 2),
                        dbout, &extval, extgrd)) return (1);
      _setLHS(lhs, EXT(iext), GRX(ig), extgrd[0]);
      _setLHS(lhs, EXT(iext), GRY(ig), extgrd[1]);
      _setLHS(lhs, EXT(iext), GRZ(ig), extgrd[2]);
    }
  }

  /* Part relative to tangents : Tx*f'x +Ty*f'y +Tz*f'z  */

  for (Id it = 0; it < _ntgt; it++)
  {
    if (_order >= 1)
    {

      /* Derivates f = x, y, et z */

      _setLHS(lhs, DRF(0), TGT(it), TGT_VAL(it, 0));
      _setLHS(lhs, DRF(1), TGT(it), TGT_VAL(it, 1));
      _setLHS(lhs, DRF(2), TGT(it), TGT_VAL(it, 2));
    }

    if (_order >= 2)
    {

      /* Derivates f = x^2, y^2, et z^2 */

      _setLHS(lhs, DRF(3), TGT(it),
              2. * TGT_COO(it, 0) * TGT_VAL(it, 0));
      _setLHS(lhs, DRF(4), TGT(it),
              2. * TGT_COO(it, 1) * TGT_VAL(it, 1));
      _setLHS(lhs, DRF(5), TGT(it),
              2. * TGT_COO(it, 2) * TGT_VAL(it, 2));

      /* Derivates f = xy, xz, et yz */

      _setLHS(lhs, DRF(6), TGT(it), (TGT_COO(it, 1) * TGT_VAL(it, 0) + TGT_COO(it, 0) * TGT_VAL(it, 1)));
      _setLHS(lhs, DRF(7), TGT(it), (TGT_COO(it, 2) * TGT_VAL(it, 0) + TGT_COO(it, 0) * TGT_VAL(it, 2)));
      _setLHS(lhs, DRF(8), TGT(it), (TGT_COO(it, 2) * TGT_VAL(it, 1) + TGT_COO(it, 1) * TGT_VAL(it, 2)));
    }

    /* External drift(s) */

    for (Id iext = 0; iext < _next; iext++)
    {
      if (_extdriftEval(TGT_COO(it, 0), TGT_COO(it, 1),
                        TGT_COO(it, 2), dbout, &extval, extgrd))
        return (1);
      _setLHS(lhs, EXT(iext), TGT(it),
              _setMatUV(TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                        extgrd[0], extgrd[1], extgrd[2]));
    }
  }

  /* Part relative to the iso-potentials */

  for (Id ic1 = 0; ic1 < _nlayers; ic1++)
  {
    for (Id j = 1; j < _nbPerLayer[ic1]; j++)
    {
      if (_order >= 1)
      {

        /* Functions x, y, z */

        _setLHS(lhs, DRF(0), ISC(ic1, j),
                ISO_COO(ic1, j, 0) - ISO_COO(ic1, 0, 0));
        _setLHS(lhs, DRF(1), ISC(ic1, j),
                ISO_COO(ic1, j, 1) - ISO_COO(ic1, 0, 1));
        _setLHS(lhs, DRF(2), ISC(ic1, j),
                ISO_COO(ic1, j, 2) - ISO_COO(ic1, 0, 2));
      }

      /* Functions x^2, y^2, z^2 */

      if (_order >= 2)
      {
        _setLHS(lhs, DRF(3), ISC(ic1, j),
                ISO_COO(ic1, j, 0) * ISO_COO(ic1, j, 0) -
                  ISO_COO(ic1, 0, 0) * ISO_COO(ic1, 0, 0));
        _setLHS(lhs, DRF(4), ISC(ic1, j),
                ISO_COO(ic1, j, 1) * ISO_COO(ic1, j, 1) -
                  ISO_COO(ic1, 0, 1) * ISO_COO(ic1, 0, 1));
        _setLHS(lhs, DRF(5), ISC(ic1, j),
                ISO_COO(ic1, j, 2) * ISO_COO(ic1, j, 2) -
                  ISO_COO(ic1, 0, 2) * ISO_COO(ic1, 0, 2));

        /* Functions xy,xz, yz */

        _setLHS(lhs, DRF(6), ISC(ic1, j),
                ISO_COO(ic1, j, 0) * ISO_COO(ic1, j, 1) -
                  ISO_COO(ic1, 0, 0) * ISO_COO(ic1, 0, 1));
        _setLHS(lhs, DRF(7), ISC(ic1, j),
                ISO_COO(ic1, j, 0) * ISO_COO(ic1, j, 2) -
                  ISO_COO(ic1, 0, 0) * ISO_COO(ic1, 0, 2));
        _setLHS(lhs, DRF(8), ISC(ic1, j),
                ISO_COO(ic1, j, 1) * ISO_COO(ic1, j, 2) -
                  ISO_COO(ic1, 0, 1) * ISO_COO(ic1, 0, 2));
      }

      /* External drift(s) */

      for (Id iext = 0; iext < _next; iext++)
      {
        if (_extdriftEval(ISO_COO(ic1, j, 0),
                          ISO_COO(ic1, j, 1), ISO_COO(ic1, j, 2), dbout,
                          &extval2, extgrd)) return (1);
        if (_extdriftEval(ISO_COO(ic1, 0, 0),
                          ISO_COO(ic1, 0, 1), ISO_COO(ic1, 0, 2), dbout,
                          &extval1, extgrd)) return (1);
        _setLHS(lhs, EXT(iext), ISC(ic1, j), extval2 - extval1);
      }
    }
  }

  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    print_matrix("LHS", 0, 1, _nequa, _nequa, NULL, lhs.getValues().data());

  return (0);
}

/****************************************************************************/
/*!
 **  Calculate the estimation and gradient components at one target location
 **
 ** \param[in]  flag_grad     True if the gradients must also be calculated
 ** \param[in]  dbgrid        Output Grid Db structure (for External Drift)
 ** \param[in]  coor          Coordinates of the target
 **
 ** \param[out] rhs           R.H.S. vector (Dimension: nequa * 4)
 **
 *****************************************************************************/
void Potential::_buildRHS(bool flag_grad,
                          DbGrid* dbgrid,
                          VectorDouble& coor,
                          MatrixDense& rhs)
{
  Id nsol = (flag_grad) ? 1 + _ndim : 1;
  rhs.fill(0.);

  double extval = 0.;
  double covar  = 0.;
  double covar1 = 0.;
  VectorDouble covGp(3, 0);
  VectorDouble cov1Gp(3, 0.);
  VectorDouble center(3, 0.);
  VectorDouble ccor(3, 0.);
  VectorDouble extgrd(3, 0.);
  VectorDouble covGG(9, 0.);
  VectorDouble cov1GG(9, 0.);

  /*******************/
  /* Covariance part */
  /*******************/

  /* Part relative to gradients */

  for (Id ig = 0; ig < _ngrd; ig++)
  {
    _calculateCovs(_model, flag_grad,
                   GRD_COO(ig, 0), GRD_COO(ig, 1), GRD_COO(ig, 2),
                   coor[0], coor[1], coor[2],
                   covar, covGp, covGG);
    _setRHS(rhs, GRX(ig), 0, covGp[0]);
    _setRHS(rhs, GRY(ig), 0, covGp[1]);
    _setRHS(rhs, GRZ(ig), 0, covGp[2]);
    if (flag_grad)
    {
      _setRHS(rhs, GRX(ig), 1, covGG[0]);
      _setRHS(rhs, GRY(ig), 1, covGG[1]);
      _setRHS(rhs, GRZ(ig), 1, covGG[2]);
      _setRHS(rhs, GRX(ig), 2, covGG[3]);
      _setRHS(rhs, GRY(ig), 2, covGG[4]);
      _setRHS(rhs, GRZ(ig), 2, covGG[5]);
      _setRHS(rhs, GRX(ig), 3, covGG[6]);
      _setRHS(rhs, GRY(ig), 3, covGG[7]);
      _setRHS(rhs, GRZ(ig), 3, covGG[8]);
    }
  }

  /* Part relative to tangents */

  for (Id it = 0; it < _ntgt; it++)
  {
    _calculateCovs(_model, flag_grad,
                   TGT_COO(it, 0), TGT_COO(it, 1), TGT_COO(it, 2),
                   coor[0], coor[1], coor[2],
                   covar, covGp, covGG);
    _setRHS(rhs, TGT(it), 0,
            _setMatUV(TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                      covGp[0], covGp[1], covGp[2]));
    if (flag_grad)
    {
      _setRHS(rhs, TGT(it), 1,
              _setMatUV(TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                        covGG[0], covGG[1], covGG[2]));
      _setRHS(rhs, TGT(it), 2,
              _setMatUV(TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                        covGG[3], covGG[4], covGG[5]));
      _setRHS(rhs, TGT(it), 3,
              _setMatUV(TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                        covGG[6], covGG[7], covGG[8]));
    }
  }

  /* Part relative to iso-potentials */

  for (Id ic = 0; ic < _nlayers; ic++)
  {
    for (Id j = 1; j < _nbPerLayer[ic]; j++)
    {
      _calculateCovs(_model, flag_grad,
                     ISO_COO(ic, j, 0), ISO_COO(ic, j, 1), ISO_COO(ic, j, 2),
                     coor[0], coor[1], coor[2],
                     covar1, cov1Gp, cov1GG);
      _calculateCovs(_model, flag_grad,
                     ISO_COO(ic, 0, 0), ISO_COO(ic, 0, 1), ISO_COO(ic, 0, 2),
                     coor[0], coor[1], coor[2],
                     covar, covGp, covGG);
      _setRHS(rhs, ISC(ic, j), 0, covar1 - covar);
      if (flag_grad)
      {
        _setRHS(rhs, ISC(ic, j), 1, -(cov1Gp[0] - covGp[0]));
        _setRHS(rhs, ISC(ic, j), 2, -(cov1Gp[1] - covGp[1]));
        _setRHS(rhs, ISC(ic, j), 3, -(cov1Gp[2] - covGp[2]));
      }
    }
  }

  /****************************/
  /* Part linked to the drift */
  /****************************/

  for (Id i = 0; i < 3; i++)
    ccor[i] = coor[i] - center[i];

  if (_order >= 1)
  {
    _setRHS(rhs, DRF(0), 0, ccor[0]);
    _setRHS(rhs, DRF(1), 0, ccor[1]);
    _setRHS(rhs, DRF(2), 0, ccor[2]);
    if (flag_grad)
    {
      _setRHS(rhs, DRF(0), 1, 1.);
      _setRHS(rhs, DRF(1), 2, 1.);
      _setRHS(rhs, DRF(2), 3, 1.);
    }
  }

  if (_order >= 2)
  {
    _setRHS(rhs, DRF(3), 0, ccor[0] * ccor[0]);
    _setRHS(rhs, DRF(4), 0, ccor[1] * ccor[1]);
    _setRHS(rhs, DRF(5), 0, ccor[2] * ccor[2]);
    _setRHS(rhs, DRF(6), 0, ccor[0] * ccor[1]);
    _setRHS(rhs, DRF(7), 0, ccor[0] * ccor[2]);
    _setRHS(rhs, DRF(8), 0, ccor[1] * ccor[2]);
    if (flag_grad)
    {
      _setRHS(rhs, DRF(3), 1, ccor[0] * 2.);
      _setRHS(rhs, DRF(4), 2, ccor[1] * 2.);
      _setRHS(rhs, DRF(5), 3, ccor[2] * 2.);
      _setRHS(rhs, DRF(6), 1, ccor[1]);
      _setRHS(rhs, DRF(6), 2, ccor[0]);
      _setRHS(rhs, DRF(7), 1, ccor[2]);
      _setRHS(rhs, DRF(7), 3, ccor[0]);
      _setRHS(rhs, DRF(8), 2, ccor[2]);
      _setRHS(rhs, DRF(8), 3, ccor[1]);
    }
  }

  for (Id iext = 0; iext < _next; iext++)
  {
    if (_extdriftEval(coor[0], coor[1], coor[2], dbgrid,
                      &extval, extgrd)) return;
    _setRHS(rhs, EXT(iext), 0, extval);
    if (flag_grad)
    {
      _setRHS(rhs, EXT(iext), 1, extgrd[0]);
      _setRHS(rhs, EXT(iext), 2, extgrd[1]);
      _setRHS(rhs, EXT(iext), 3, extgrd[2]);
    }
  }

  // Blank out the R.H.S. according to masking option

  _blankPartRHS(rhs);

  // Printout (optional)

  if (OptDbg::query(EDbg::KRIGING))
    print_matrix("RHS", 0, 1, _nequa, nsol, NULL, rhs.getValues().data());
}

/****************************************************************************/
/*!
 **  Blank out part the R.H.S. according to 'flag.part'
 **
 ** \param[in,out] rhs        Array for the R.H.S.
 **
 *****************************************************************************/
void Potential::_blankPartRHS(MatrixDense& rhs) const
{
  Id ideb = 0;
  Id ifin = _nequa;
  if (_optionPart == 0) return;

  /* Dispatch */

  switch (_optionPart)
  {
    case 1: /* Reveal Gradient */
      ideb = _startGrd;
      ifin = ideb + _sizeGrd;
      break;

    case 2: /* Reveal Tangent */
      ideb = _startTgt;
      ifin = ideb + _sizeTgt;
      break;

    case 3: /* Reveal Isovalues */
      ideb = _startIso;
      ifin = ideb + _sizeIso;
      break;

    case 4: /* Reveal internal drift */
      ideb = _startDrf;
      ifin = ideb + _sizeDrf;
      break;

    case 5: /* Reveal external drift */
      ideb = _startExt;
      ifin = ideb + _sizeExt;
      break;

    default:
      break;
  }

  /* Blank out the R.H.S. */

  for (Id i = 0; i < _nequa; i++)
  {
    if (i >= ideb && i < ifin) continue;
    _setRHS(rhs, i, 0, 0.);
    if (_flagGrad)
      for (Id igrad = 1; igrad < 4; igrad++)
        _setRHS(rhs, i, igrad, 0.);
  }
}

Id Potential::_extdriftCreateDb(DbGrid* dbout)
{
  VectorInt nx(_ndim);
  VectorDouble x0(_ndim);

  /* Creating the attributes from the output grid */

  Id nech = 1;
  for (Id idim = 0; idim < _ndim; idim++)
  {
    nx[idim] = 2 * _nring + 1;
    x0[idim] = -dbout->getDX(idim) * _nring;
    nech *= nx[idim];
  }

  /* Creating the data grid */

  _dbExt = DbGrid::create(nx, dbout->getDXs(), x0, dbout->getAngles(),
                          ELoadBy::COLUMN, VectorDouble(),
                          VectorString(), VectorString(), 1);
  if (_dbExt == nullptr) return 1;
  _nfull = nech;

  /* Add the selection */

  _dbExt->addColumnsByConstant(1, 0., String(), ELoc::SEL);

  /* Complementary core allocation */

  _dataExt.resize(nech);
  _wgtExt.resize(nech, 4);
  _indg0.resize(3);
  _indg.resize(3);
  return 0;
}

Id Potential::_extdriftCreateModel()
{
  CovContext ctxt(1, _ndim, 1.);
  _modelExt = new ModelGeneric(ctxt);

  // Covariance part
  const CovAniso* cova = CovAniso::createFromParam(ECov::CUBIC, _rangeExt, 1.);
  auto covp            = CovGradientAnalytic(*cova);
  _modelExt->setCov(&covp);

  // Drift part
  DriftList drifts(ctxt);
  drifts.setFlagLinked(true);
  _modelExt->setDriftList(&drifts);

  return 0;
}

MatrixDense Potential::_extdriftBuildRHS()
{
  double covar = 0.;
  VectorDouble covGp(3, 0.);
  VectorDouble covGG(9, 0.);

  /* Establish the kriging matrix */

  MatrixDense b(_nfull, 4);

  /* Establish the Right-Hand side */

  Id ecr = 0;
  for (Id iech = 0; iech < _nfull; iech++)
  {
    if (!_dbExt->isActive(iech)) continue;
    // TODO: check the validity of the bi-point for covariance calculation
    _calculateCovs(_modelExt, true,
                   _dbExt->getCoordinate(iech, 0),
                   _dbExt->getCoordinate(iech, 1),
                   _dbExt->getCoordinate(iech, 2),
                   0., 0., 0.,
                   covar, covGp, covGG);
    b.setValue(ecr, 0, covar);
    b.setValue(ecr, 1, -covGp[0]);
    b.setValue(ecr, 2, -covGp[1]);
    b.setValue(ecr, 3, -covGp[2]);
    ecr++;
  }
  return b;
}

Id Potential::_extdriftEstablish(DbGrid* dbout)
{
  /* Creating the Db for neighborhood */

  if (_extdriftCreateDb(dbout)) return 1;

  /* Creating the model */

  if (_extdriftCreateModel()) return 1;

  /* Solve the kriging system */

  MatrixSymmetric a = _modelExt->evalCovMatSym(_dbExt);
  if (a.invert()) return 1;

  MatrixDense b = _extdriftBuildRHS();

  a.prodMatMatInPlace(&b, &_wgtExt);

  return 0;
}

/****************************************************************************/
/*!
 **  Establish the local neighborhood
 **
 ** \return  Error return code (Neighborhood not complete)
 **
 ** \param[in]  dbgrid     Output Db Grid structure
 **
 *****************************************************************************/
Id Potential::_extdriftNeigh(DbGrid* dbgrid)
{
  /* Loop on the neighboring samples defined in the neighboring grid */

  Id ecr = 0;
  for (Id iz = 0; iz < _dbExt->getNX(2); iz++)
    for (Id iy = 0; iy < _dbExt->getNX(1); iy++)
      for (Id ix = 0; ix < _dbExt->getNX(0); ix++)
      {

        /* Calculate the index of the sample within the Ext Drift grid */

        _indg[0] = _indg0[0] + ix - _nring;
        if (_indg[0] < 0 || _indg[0] > dbgrid->getNX(0)) return (1);
        _indg[1] = _indg0[1] + iy - _nring;
        if (_indg[1] < 0 || _indg[1] > dbgrid->getNX(1)) return (1);
        _indg[2] = _indg0[2] + iz - _nring;
        if (_indg[2] < 0 || _indg[2] > dbgrid->getNX(2)) return (1);
        Id iech = dbgrid->indiceToRank(_indg);

        /* Check that the external drift value is defined */

        double drift = dbgrid->getLocVariable(ELoc::F, iech, 0);
        if (FFFF(drift)) return (1);
        _dataExt[ecr] = drift;
        ecr++;
      }
  return (0);
}

/****************************************************************************/
/*!
 **  Calculate the external drift contribution
 **
 ** \return  Error return code (target not within the grid or target on the
 ** \return  edge of the grid of external drift definition)
 **
 ** \param[in] x0       Coordinate along X
 ** \param[in] y0       Coordinate along Y
 ** \param[in] z0       Coordinate along Z
 ** \param[in] db       Output Db structure (must be a grid)
 **
 ** \param[out] extval  Value of the external drift
 ** \param[out] extgrd  Gradient components of the external drift
 **
 *****************************************************************************/
Id Potential::_extdriftEval(double x0,
                            double y0,
                            double z0,
                            Db* db,
                            double* extval,
                            VectorDouble& extgrd)
{
  if (db == nullptr) return 1;
  auto* dbgrid = dynamic_cast<DbGrid*>(db);
  if (dbgrid == nullptr) return 1;
  VectorDouble coor(3);
  coor[0] = x0;
  coor[1] = y0;
  coor[2] = z0;

  /* Find the location of the target within the external drift grid */

  if (point_to_grid(dbgrid, coor.data(), 0, _indg0.data()) < 0) return 1;

  /* Find the neighborhood around the target grid node */

  if (_extdriftNeigh(dbgrid)) return 1;

  /* Perform the estimation */

  VectorDouble result(4);
  _wgtExt.prodVecMatInPlace(_dataExt, result);

  /* Retrieve the results */

  *extval = result[0];
  for (Id idim = 0; idim < _ndim; idim++)
    extgrd[idim] = result[1 + idim];

  return 0;
}

void Potential::_fillDual(VectorDouble& zval)
{
  zval.resize(_nequa);
  zval.fill(0.);

  for (Id ig = 0; ig < _ngrd; ig++)
  {
    if (GRX(ig) >= 0) zval[GRX(ig)] = GRD_VAL(ig, 0);
    if (GRY(ig) >= 0) zval[GRY(ig)] = GRD_VAL(ig, 1);
    if (GRZ(ig) >= 0) zval[GRZ(ig)] = GRD_VAL(ig, 2);
  }
}

/****************************************************************************/
/*!
 **  Establish the simulation errors
 **
 ** \param[in]  nbsimu        Number of simulations
 **
 ** \param[out] zvals         Simulated errors
 **
 *****************************************************************************/
void Potential::_fillDualSimulation(Id nbsimu, MatrixDense& zvals)
{
  zvals.fill(0.);

  // Loop on the simulations */

  for (Id isimu = 0; isimu < nbsimu; isimu++)
  {

    // Load the gradient simulation errors

    for (Id ig = 0; ig < _ngrd; ig++)
    {
      if (_ndim >= 1)
        zvals.setValue(GRX(ig), isimu, _dbgrd->getSimvar(ELoc::SIMU, IAD_GRD(ig), isimu + 0 * nbsimu, 0, 0, _ndim * nbsimu, 1) - GRD_VAL(ig, 0));
      if (_ndim >= 2)
        zvals.setValue(GRY(ig), isimu, _dbgrd->getSimvar(ELoc::SIMU, IAD_GRD(ig), isimu + 1 * nbsimu, 0, 0, _ndim * nbsimu, 1) - GRD_VAL(ig, 1));
      if (_ndim >= 3)
        zvals.setValue(GRZ(ig), isimu, _dbgrd->getSimvar(ELoc::SIMU, IAD_GRD(ig), isimu + 2 * nbsimu, 0, 0, _ndim * nbsimu, 1) - GRD_VAL(ig, 2));
    }

    // Load the tangent simulation errors

    for (Id it = 0; it < _ntgt; it++)
    {
      zvals.setValue(TGT(it), isimu, _dbtgt->getSimvar(ELoc::SIMU, IAD_TGT(it), isimu, 0, 0, nbsimu, 1));
    }

    // Load the iso-potential simulation errors

    for (Id ic = 0; ic < _nlayers; ic++)
      for (Id j = 1; j < _nbPerLayer[ic]; j++)
      {
        zvals.setValue(ISC(ic, j), isimu, _dbiso->getSimvar(ELoc::SIMU, IAD_ISO(ic, j), isimu, 0, 0, nbsimu, 1) - _dbiso->getSimvar(ELoc::SIMU, IAD_ISO(ic, 0), isimu, 0, 0, nbsimu, 1));
      }
  }

  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    print_matrix("\n[Simu-Err]", 0, 1, nbsimu, _nequa, NULL, zvals.getValues().data());
}

void Potential::_calculatePoint(bool flag_grad,
                                DbGrid* dbgrid,
                                VectorDouble& zdual,
                                MatrixDense& rhs,
                                Db* db_target,
                                Id iech0,
                                VectorDouble& result)
{
  VectorDouble coor(3, 0.);

  /* Initializations */

  Id nsol = (flag_grad) ? 1 + _ndim : 1;

  /* Load the coordinates */

  for (Id idim = 0; idim < _ndim; idim++)
    coor[idim] = db_target->getCoordinate(iech0, idim);

  /* Optional printout */

  if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH))
  {
    mestitle(1, "Target location");
    db_sample_print(db_target, iech0, 1, 0, 0, 0);
  }

  /* Establish the R.H.S */

  _buildRHS(flag_grad, dbgrid, coor, rhs);

  /* Perform the estimation */

  result.fill(TEST);
  rhs.prodVecMatInPlace(zdual, result);

  // Printout (optional)

  if (OptDbg::query(EDbg::KRIGING))
  {
    print_matrix("Results", 0, 1, 1, nsol, NULL, result.data());
    message("\n");
  }
}

/****************************************************************************/
/*!
 **  Translate potential into layer value or center it
 **
 ** \param[in]  isimu         Rank of the simulation
 ** \param[in]  potval        Array of potential values at different layers
 ** \param[in]  result        Resulting value (in potential scale)
 **                           On output, Resulting value in layer scale
 **
 ** \remarks The potential values at iso-potential samples are assumed
 ** \remarks to be ordered
 ** \remarks It is assumed that the potential has already been centered
 ** \remarks Therefore the 'potval' values must also be centered (locally)
 **
 *****************************************************************************/
void Potential::_convertPotentialToLayer(Id isimu,
                                         const double* potval,
                                         VectorDouble& result) const
{
  double minval = MINIMUM_BIG;
  double potref = POTVAL(isimu, 0);

  Id ilayer = -1;
  for (Id i = 0; i < _nlayers && ilayer < 0; i++)
  {
    if (result[0] > minval && result[0] <= (POTVAL(isimu, i) - potref))
      ilayer = i;
    minval = (POTVAL(isimu, i) - potref);
  }
  result[0] = ilayer + 1;
}

/****************************************************************************/
/*!
 **  Calculate the estimation and/or gradient components at target samples
 **
 ** \param[in]  flag_grad     True if the gradients must also be calculated
 ** \param[in]  dbout         Output Db structure
 ** \param[in]  refpot        Potential at Reference point
 ** \param[in]  zdual         Dual vector (Dimension: nequa)
 ** \param[in]  rhs           R.H.S. vector (Dimension: nequa * 4)
 ** \param[in]  potval        Potential values at iso-potential samples
 **
 *****************************************************************************/
void Potential::_estimateResult(bool flag_grad,
                                DbGrid* dbout,
                                double refpot,
                                VectorDouble& zdual,
                                MatrixDense& rhs,
                                double* potval)
{
  VectorDouble result(4);

  for (Id iech = 0; iech < dbout->getNSample(); iech++)
  {
    mes_process("Potential Estimation on Grid", dbout->getNSample(), iech);
    OptDbg::setCurrentIndex(iech);
    if (!dbout->isActive(iech)) continue;

    // Perform the estimation

    _calculatePoint(flag_grad, dbout, zdual, rhs, dbout, iech, result);

    // Center to the reference potential

    result[0] -= refpot;

    // Printout (optional)

    if (OptDbg::query(EDbg::KRIGING))
    {
      message("Estimation (centered) = %lf\n", result[0]);
      for (int idim = 0; idim < _ndim; idim++)
        message("Gradient %d = %lf\n", idim + 1, result[idim + 1]);
    }

    // Translate from potential into layer

    if (_flagTrans)
      _convertPotentialToLayer(0, potval, result);

    // Store the results

    dbout->setLocVariable(ELoc::Z, iech, 0, result[0]);
    if (flag_grad)
      for (Id idim = 0; idim < _ndim; idim++)
        dbout->setLocVariable(ELoc::G, iech, idim, result[idim + 1]);
  }
  OptDbg::setCurrentIndex(-1);
}

void Potential::_estimateData(DbGrid* dbout,
                              double refpot,
                              VectorDouble& zdual,
                              MatrixDense& rhs,
                              Db* db_target,
                              VectorInt& uid_pot,
                              VectorInt& uid_grad)
{
  if (db_target == nullptr) return;
  VectorDouble result(4);

  for (Id iech = 0; iech < db_target->getNSample(); iech++)
  {
    if (!db_target->isActive(iech)) continue;

    // Perform the estimation

    _calculatePoint(true, dbout, zdual, rhs, db_target, iech, result);

    // Center to the reference potential

    result[0] -= refpot;

    // Store the results

    if (!uid_pot.empty())
    {
      db_target->setArray(iech, uid_pot[0], result[0]);
      db_target->setLocatorsByUID(uid_pot, ELoc::Z, 0);
    }
    if (!uid_grad.empty())
    {
      for (Id idim = 0; idim < _ndim; idim++)
        db_target->setArray(iech, uid_grad[idim], result[idim + 1]);
      db_target->setLocatorsByUID(uid_grad, ELoc::G, 0);
    }
  }
}

/****************************************************************************/
/*!
 **  Calculate the cross-validation at the iso-potential samples
 **
 ** \param[in]  ic0           Rank of the isoline
 ** \param[in]  j0            Rank of the sample within this isoline
 ** \param[in]  zval          Data vector
 ** \param[in]  lhs_orig_arg  Copy of the initial LHS (non inverted)
 ** \param[in]  rhs_arg       Right-hand side
 **
 ** \param[out] dist_euc      Error converted into Euclidean distance
 ** \param[out] dist_geo      Error converted into along surface distance
 **
 ** \remark We assume that the new data set (with one sample OFF) is still
 ** \remark contained in 'zval' as the Gradient information (coming first in
 ** \remark this vector) is never excluded.
 **
 *****************************************************************************/
void Potential::_convertDistance(Id ic0,
                                 Id j0,
                                 VectorDouble& zval,
                                 MatrixSymmetric& lhs_orig_arg,
                                 MatrixDense& rhs_arg,
                                 double* dist_euc,
                                 double* dist_geo)
{
  VectorDouble result(4);
  static Id niter_max = 50;
  static double eps   = 1.e-3;

  VectorDouble coor(3, 0.);
  VectorDouble coor0(3, 0.);
  Id neqm1 = _nequa - 1;
  Id icol0 = ISC(ic0, j0);
  VectorDouble deuc(_ndim, 0.);
  VectorDouble dgeo(_ndim, 0.);

  VectorDouble lhs_orig = lhs_orig_arg.getValues();
  VectorDouble rhs;
  MatrixSymmetric* lhs_aux = nullptr;
  MatrixDense* rhs_red     = nullptr;

  /* Update the L.H.S. by dropping the current data point */

  lhs_aux = dynamic_cast<MatrixSymmetric*>(MatrixFactory::createReduceOne(&lhs_orig_arg, icol0, icol0, false, false));

  /* Invert the new LHS */

  if (lhs_aux->invert()) return;

  /* Calculate the dual system */

  VectorDouble zdual_red(neqm1);
  VectorDouble zval_red = VH::reduceOne(zval, icol0);
  lhs_aux->prodMatVecInPlace(zval_red, zdual_red);
  delete lhs_aux;

  /* Evaluate the reference point */

  for (Id idim = 0; idim < _ndim; idim++)
    coor0[idim] = ISO_COO(ic0, 0, idim);
  _buildRHS(false, nullptr, coor0, rhs_arg);
  rhs_red = dynamic_cast<MatrixDense*>(MatrixFactory::createReduceOne(&rhs_arg, icol0, -1, false, false));
  rhs_red->prodVecMatInPlace(zdual_red, result);
  //  double potval = result[0]; // TODO: check why is potval not used
  delete rhs_red;

  /* Evaluate the target point */

  for (Id idim = 0; idim < _ndim; idim++)
    coor0[idim] = coor[idim] = ISO_COO(ic0, j0, idim);
  _buildRHS(true, nullptr, coor0, rhs_arg);
  rhs_red = dynamic_cast<MatrixDense*>(MatrixFactory::createReduceOne(&rhs_arg, icol0, -1, false, false));
  rhs_red->prodVecMatInPlace(zdual_red, result);
  delete rhs_red;

  if (OptDbg::query(EDbg::CONVERGE))
  {
    message("Sample:%2d/%2d Iter:%2d Potential:%lf", j0 + 1, ic0 + 1, 0,
            result[0]);
    for (Id idim = 0; idim < _ndim; idim++)
      message(" %lf", coor0[idim]);
    message("\n");
  }

  /* Move the target and estimate again */

  for (Id iter = 0; iter < niter_max; iter++)
  {
    if (ABS(result[0]) < eps) break;
    for (Id idim = 0; idim < _ndim; idim++)
    {
      if (ABS(result[1 + idim]) < eps) continue;
      double delta = 0.1 * result[0] / result[1 + idim];
      coor[idim] -= delta;
      dgeo[idim] += delta * delta;
    }
    _buildRHS(true, nullptr, coor, rhs_arg);
    rhs_red = dynamic_cast<MatrixDense*>(MatrixFactory::createReduceOne(&rhs_arg, icol0, -1, false, false));
    rhs_red->prodVecMatInPlace(zdual_red, result);
    delete rhs_red;

    if (OptDbg::query(EDbg::CONVERGE))
    {
      message("Sample:%2d/%2d Iter:%2d Potential:%lf", j0 + 1, ic0 + 1, iter,
              result[0]);
      for (Id idim = 0; idim < _ndim; idim++)
        message(" %lf", coor[idim]);
      message("\n");
    }
  }

  /* Determine the euclidean distance */

  for (Id idim = 0; idim < _ndim; idim++)
  {
    double delta = coor[idim] - coor0[idim];
    deuc[idim]   = delta * delta;
  }

  /* Find both distances */

  (*dist_euc) = (*dist_geo) = 0.;
  for (Id idim = 0; idim < _ndim; idim++)
  {
    (*dist_euc) += deuc[idim];
    (*dist_geo) += dgeo[idim];
  }
  (*dist_euc) = sqrt(*dist_euc);
  (*dist_geo) = sqrt(*dist_geo);
}

/****************************************************************************/
/*!
 **  Calculate the cross-validation at the iso-potential samples
 **
 ** \param[in]  lhs            Inverted Kriging system
 ** \param[in]  flag_dist_conv Flag for converting into distance
 ** \param[in]  zval           Data vector
 ** \param[in]  lhs_orig       Copy of the Initial LHS
 ** \param[in]  rhs            Right-hand side
 **
 ** \remarks Arguments from 'zval' are only used to convert into distance
 **
 *****************************************************************************/
void Potential::_xvalidCalculate(MatrixSymmetric& lhs,
                                 bool flag_dist_conv,
                                 VectorDouble& zval,
                                 MatrixSymmetric& lhs_orig,
                                 MatrixDense& rhs)
{
  double stats[4][2];
  VectorDouble result(4);

  /* Initializations */

  Id nitem = (flag_dist_conv) ? 4 : 2;
  for (Id i = 0; i < nitem; i++)
    for (Id j = 0; j < 2; j++)
      stats[i][j] = 0.;

  /* Loop on the Iso-potential samples */

  for (Id ic = 0; ic < _nlayers; ic++)
  {
    Id number = 0;
    for (Id j = 1; j < _nbPerLayer[ic]; j++)
    {
      Id iech0 = IAD_ISO(ic, j);
      mes_process("Potential Estimation on Iso-Potential %d of %d", j + 1, ic + 1);
      OptDbg::setCurrentIndex(iech0);

      // Get the variance and the weights from the inverted L.H.S.

      Id icol0        = ISC(ic, j);
      double variance = 1. / _getLHS(lhs, icol0, icol0);
      double stdev    = sqrt(variance);
      double dist_geo = 0.;
      double dist_euc = 0.;

      // Perform the estimation

      double value = 0.;
      for (Id ig = 0; ig < _ngrd; ig++)
      {
        value += _getLHS(lhs, icol0, GRX(ig)) * GRD_VAL(ig, 0);
        value += _getLHS(lhs, icol0, GRY(ig)) * GRD_VAL(ig, 1);
        value += _getLHS(lhs, icol0, GRZ(ig)) * GRD_VAL(ig, 2);
      }
      result[0] = -value * variance;

      // Finding the closest distance to the Isoline

      if (flag_dist_conv)
        _convertDistance(ic, j, zval, lhs_orig, rhs, &dist_geo, &dist_euc);

      // Debugging option

      if (OptDbg::query(EDbg::RESULTS))
      {
        message("Sample %d/%d (%d): Error=%lf - Variance=%lf", j + 1, ic + 1,
                iech0 + 1, result[0], variance);
        if (flag_dist_conv)
          message(" - D-Geo=%lf - D-Surf=%lf", dist_euc, dist_geo);
        message("\n");
      }

      // Storing the results

      _dbiso->setLocVariable(ELoc::Z, iech0, 0, result[0]);
      _dbiso->setLocVariable(ELoc::Z, iech0, 1, variance);
      if (flag_dist_conv)
      {
        _dbiso->setLocVariable(ELoc::Z, iech0, 2, dist_euc);
        _dbiso->setLocVariable(ELoc::Z, iech0, 3, dist_geo);
      }

      // Update statistics */

      value = result[0];
      stats[0][0] += value;
      stats[0][1] += value * value;
      value = result[0] / stdev;
      stats[1][0] += value;
      stats[1][1] += value * value;
      if (flag_dist_conv)
      {
        value = dist_geo;
        stats[2][0] += value;
        stats[2][1] += value * value;
        value = dist_euc;
        stats[3][0] += value;
        stats[3][1] += value * value;
      }
      number++;
    }

    // Print the global statistics (optinal)

    if (_verbose && number > 0)
    {
      for (Id i = 0; i < nitem; i++)
      {
        for (Id j = 0; j < 2; j++)
          stats[i][j] /= number;
        stats[i][1] -= stats[i][0] * stats[i][0];
        stats[i][1] = (stats[i][1] > 0) ? sqrt(stats[i][1]) : 0.;
      }
      message("\nIso-Potential #%d\n", ic + 1);
      message("Cross-validation Error: Mean=%lf St. Dev.=%lf\n", stats[0][0],
              stats[0][1]);
      message("Standardized Error    : Mean=%lf St. Dev.=%lf\n", stats[1][0],
              stats[1][1]);
      if (flag_dist_conv)
      {
        message("Euclidean Distance    : Mean=%lf St. Dev.=%lf\n", stats[2][0],
                stats[2][1]);
        message("Geodetic Distance     : Mean=%lf St. Dev.=%lf\n", stats[3][0],
                stats[3][1]);
      }
    }
  }
  OptDbg::setCurrentIndex(-1);
}

/****************************************************************************/
/*!
 **  Amortize the conditional simulations
 **
 ** \param[in]  dbout         Output Db structure
 ** \param[in]  iech          Rank of the sample
 ** \param[in]  dist_tempere  Distance for tempering simulations (or TEST)
 ** \param[in]  reskrige      Kriging result
 ** \param[in]  result        Conditional Simulation result
 **                           On output, Conditional Simulation tempered result
 **
 ** \remarks This function does nothing if 'dist_tempere' is undefined
 **
 *****************************************************************************/
void Potential::_tempere(DbGrid* dbout,
                         Id iech,
                         double dist_tempere,
                         double reskrige,
                         VectorDouble& result)
{
  static Id test = 0;

  double simerr = result[0] - reskrige;
  double kdist  = dbout->getZVariable(iech, 0);

  switch (test)
  {
    case 0: /* Simulation amortie */
    {
      double amortval = MIN(1., exp(-kdist / dist_tempere));
      result[0]       = reskrige + simerr * amortval;
      break;
    }

    case 1: /* Distance normee */
    {
      result[0] = kdist / dist_tempere;
      break;
    }

    case 2: /* Simulation Conditionnelle */
    {
      break;
    }

    case 3: /* Simulation non-conditionnelle */
    {
      result[0] = simerr;
      break;
    }

    case 4: /* Krigeage */
    {
      result[0] = reskrige;
      break;
    }

    default:
      break;
  }
}

/****************************************************************************/
/*!
 **  Calculate the conditional simulation at target samples
 **
 ** \param[in]  dist_tempere  Distance for tempering simulations (or TEST)
 ** \param[in]  flag_trans    True if the estimation result must be translated
 **                           into layer number
 ** \param[in]  nbsimu        Number of simulation
 ** \param[in]  dbout         Output Db structure
 ** \param[in]  refpot        Potential at Reference point
 ** \param[in]  potsim        Potential simulated values at iso-potential samples
 ** \param[in]  zdual         Dual estimated vector (Dimension: nequa)
 ** \param[in]  zduals        Dual simulated vector (Dimension: nequa * nbsimu)
 ** \param[in]  rhs           R.H.S. vector (Dimension: nequa * 4)
 **
 *****************************************************************************/
void Potential::_simcond(double dist_tempere,
                         bool flag_trans,
                         Id nbsimu,
                         DbGrid* dbout,
                         double refpot,
                         double* potsim,
                         VectorDouble& zdual,
                         MatrixDense& zduals,
                         MatrixDense& rhs)
{
  VectorDouble resest(4), result(4);

  Id ndim = _dbgrd->getNDim();
  for (Id iech = 0; iech < dbout->getNSample(); iech++)
  {
    mes_process("Potential Simulation on Grid", dbout->getNSample(), iech);
    OptDbg::setCurrentIndex(iech);
    if (!dbout->isActive(iech)) continue;

    if (!FFFF(dist_tempere))
    {

      // Perform the estimation

      _calculatePoint(true, dbout, zdual, rhs, dbout, iech, resest);

      // Center to the reference potential

      resest[0] -= refpot;
    }

    for (Id isimu = 0; isimu < nbsimu; isimu++)
    {

      // Perform the estimation of the simulated error

      VectorDouble zdual_loc = zduals.getColumn(isimu);
      _calculatePoint(false, dbout, zdual_loc, rhs, dbout, iech, result);

      // Convert into simulation error

      result[0] = (dbout->getSimvar(ELoc::SIMU, iech, isimu, 0, 0, nbsimu, 1) - result[0]);
      for (Id idim = 0; idim < ndim; idim++)
        result[1 + idim] = (_dbgrd->getSimvar(ELoc::SIMU, iech,
                                              isimu + idim * nbsimu, 0, 0,
                                              ndim * nbsimu, 1) -
                            result[1 + idim]);

      // Center to the reference potential

      result[0] -= refpot;

      // Amortize the variance for conditional simulation

      if (!FFFF(dist_tempere))
        _tempere(dbout, iech, dist_tempere, resest[0], result);

      // Translate from potential into layer

      if (flag_trans)
        _convertPotentialToLayer(isimu, potsim, result);

      // Store the results

      dbout->setSimvar(ELoc::SIMU, iech, isimu, 0, 0, nbsimu, 1, result[0]);
    }
  }
  OptDbg::setCurrentIndex(-1);
}

/****************************************************************************/
/*!
 **  Print the estimation at a target sample
 **
 ** \param[in]  isimu      Rank of the simulation (or -1)
 ** \param[in]  result     Array of results
 ** \param[in]  tgtval     Value of the tangent (or TEST)
 **
 *****************************************************************************/
void Potential::_printResult(Id isimu, double* result, double tgtval) const
{
  if (isimu >= 0) message("Simulation %2d - ", isimu + 1);

  message(" - Pot* =%10.5lf", roundZero(result[0]));

  message(" - Grad* =");
  for (Id idim = 0; idim < _ndim; idim++)
    message(" %10.5lf", roundZero(result[1 + idim]));

  if (!FFFF(tgtval)) message(" - Tang* = %10.5lf", roundZero(tgtval));

  message("\n");
}

/****************************************************************************/
/*!
 **  Calculate the estimation and gradient components at data information
 **
 ** \param[in]  dbgrid        Output Db structure (for External drift)
 ** \param[in]  isimu         Rank of the simulation (or -1)
 ** \param[in]  nbsimu        Number of simulations
 ** \param[in]  refpot        Potential at Reference point
 ** \param[in]  zdual         Dual vector (Dimension: nequa)
 ** \param[in]  rhs           R.H.S. vector (Dimension: nequa * 4)
 **
 *****************************************************************************/
void Potential::_checkData(DbGrid* dbgrid,
                           Id isimu,
                           Id nbsimu,
                           double refpot,
                           VectorDouble& zdual,
                           MatrixDense& rhs)
{
  VectorDouble result(4);

  /* Preliminary check */

  if (_verbose) mestitle(0, "Information completed at Data Points");

  /* For the Iso-Potential file */

  if (_dbiso != nullptr)
  {
    if (_verbose) mestitle(1, "Iso-Potential Information");

    Id rank = 0;
    for (Id ic = 0; ic < _nlayers; ic++)
    {
      for (Id j = 0; j < _nbPerLayer[ic]; j++, rank++)
      {
        Id iech = _dbiso->getRankRelativeToAbsolute(rank);
        _calculatePoint(true, dbgrid, zdual, rhs, _dbiso, iech, result);
        result[0] -= refpot;

        // Printout (conditional)

        if (_verbose)
        {
          // Convert into simulation error

          if (nbsimu > 0)
            result[0] = _dbiso->getSimvar(ELoc::SIMU, iech, isimu, 0, 0, nbsimu, 1) - result[0];

          // Print the results

          message(" %2d - %2d - Coor =", ic + 1, j + 1);
          for (Id idim = 0; idim < _ndim; idim++)
            message(" %7.4lf", ISO_COO(ic, j, idim));
          _printResult(isimu, result.data(), TEST);
        }
      }
    }
  }

  /* For the Gradient file */

  if (_dbgrd != nullptr)
  {
    if (_verbose) mestitle(1, "Gradient Information");

    for (Id ig = 0; ig < _ngrd; ig++)
    {
      Id iech = _dbgrd->getRankRelativeToAbsolute(ig);
      _calculatePoint(true, dbgrid, zdual, rhs, _dbgrd, iech, result);
      result[0] -= refpot;

      // Printout (optional)

      if (_verbose)
      {

        // Convert into simulation error

        if (nbsimu > 0)
        {
          for (Id idim = 0; idim < _ndim; idim++)
            result[1 + idim] = (_dbgrd->getSimvar(ELoc::SIMU, iech,
                                                  isimu + idim * nbsimu, 0, 0,
                                                  _ndim * nbsimu, 1) -
                                result[1 + idim]);
        }

        // Print the results

        message(" %2d - Coor =", ig + 1);
        for (Id idim = 0; idim < _ndim; idim++)
          message(" %7.4lf", GRD_COO(ig, idim));
        _printResult(isimu, result.data(), TEST);
      }
    }
  }

  /* For the Tangent file */

  if (_dbtgt != nullptr)
  {
    if (_verbose) mestitle(1, "Tangent Information");

    for (Id it = 0; it < _ntgt; it++)
    {
      Id iech = _dbtgt->getRankRelativeToAbsolute(it);
      if (!_dbtgt->isActive(iech)) continue;
      _calculatePoint(true, dbgrid, zdual, rhs, _dbtgt, iech, result);
      result[0] -= refpot;

      // Printout (conditional)

      if (_verbose)
      {
        double tgte = 0.;
        for (Id idim = 0; idim < _ndim; idim++)
          tgte += result[1 + idim] * _dbtgt->getLocVariable(ELoc::TGTE, iech, idim);

        // Print the results

        message(" %2d - Coor =", it + 1);
        for (Id idim = 0; idim < _ndim; idim++)
          message(" %7.4lf", TGT_COO(it, idim));
        _printResult(isimu, result.data(), tgte);
      }
    }
  }
}

/****************************************************************************/
/*!
 **  Calculate the estimation at the potential at first point of first potential
 **
 ** \return The Potential value at first point of first iso-potential
 **
 ** \param[in]  dbgrid        Ouput Db structure (for external drift)
 ** \param[in]  zdual         Dual vector (Dimension: nequa)
 ** \param[in]  rhs           R.H.S. vector (Dimension: nequa * 4)
 **
 *****************************************************************************/
double Potential::_evaluateRefPot(DbGrid* dbgrid,
                                  VectorDouble& zdual,
                                  MatrixDense& rhs)
{
  if (_dbiso == nullptr) return (TEST);
  VectorDouble result(4);

  // Calculate the reference values for iso-potentials

  Id ic  = 0;
  Id ip1 = IAD_ISO(ic, 0);
  _calculatePoint(false, dbgrid, zdual, rhs, _dbiso, ip1, result);
  return (result[0]);
}

/****************************************************************************/
/*!
 **  Calculate the estimation at the iso-potential samples
 **
 ** \param[in]  dbgrid        Ouput Db structure (for external drift)
 ** \param[in]  refpot        Potential at Reference point
 ** \param[in]  isimu         Rank of the simulation (or -1)
 ** \param[in]  nbsimu        Number of simulations (or 0)
 ** \param[in]  zdual         Dual vector (Dimension: nequa)
 ** \param[in]  rhs           R.H.S. vector (Dimension: nequa * 4)
 **
 ** \param[out] potval        Array of Potential values
 **
 *****************************************************************************/
void Potential::_evaluatePotential(DbGrid* dbgrid,
                                   double refpot,
                                   Id isimu,
                                   Id nbsimu,
                                   VectorDouble& zdual,
                                   MatrixDense& rhs,
                                   double* potval)
{
  if (_dbiso == nullptr) return;
  VectorDouble result(4);

  // Calculate the reference values for isopotentials

  for (Id ic = 0; ic < _nlayers; ic++)
  {
    Id ip1 = IAD_ISO(ic, 0);
    _calculatePoint(false, dbgrid, zdual, rhs, _dbiso, ip1, result);

    // Convert into simulation error

    if (nbsimu > 0)
      result[0] = (_dbiso->getSimvar(ELoc::SIMU, ip1, isimu, 0, 0, nbsimu, 1) - result[0]);

    // Center to the reference potential

    result[0] -= refpot;

    // Store in the 'potval' array

    potval[ic] = result[0];
  }

  // Sort them by ascending order

  ut_sort_double(0, _nlayers, NULL, potval);
}

/****************************************************************************/
/*!
 **  Transform the Estimation variable into a distance to the isoline
 **
 ** \return  Error return code
 **
 ** \param[in]  dbout        Output Db structure
 **
 *****************************************************************************/
Id Potential::_distanceToIsoline(DbGrid* dbout)

{
  Id radius  = 1;
  Id seed    = 3432521;
  auto memo  = law_get_random_seed();
  double eps = 1.e-3;

  // Highlight the isoline of interest
  for (Id iech = 0; iech < dbout->getNSample(); iech++)
  {
    double value = dbout->getZVariable(iech, 0);
    if (!FFFF(value) && ABS(value) > eps) dbout->setLocVariable(ELoc::Z, iech, 0, TEST);
  }

  // Calculate the distance
  if (DbHelper::dbgrid_filling(dbout, 3, seed, radius)) return 1;

  law_set_random_seed(memo);
  return (0);
}

/****************************************************************************/
/*!
 **  Potential estimation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbout      Output Db structure
 ** \param[in]  flag_pot   True if the Potential must be estimated
 ** \param[in]  flag_grad  True if the gradient must also be estimated
 ** \param[in]  flag_trans True if the estimation result must be translated
 **                        into layer number
 ** \param[in]  flag_save_data True if the Potential / Gradient must be
 **                        saved on any Information file
 ** \param[in]  option_part   Option to exhibit only a part of estimation:
 ** \li                       0 : the whole estimation
 ** \li                       1 : the gradient contribution only
 ** \li                       2 : the tangent contribution only
 ** \li                       3 : the isovalues contribution only
 ** \li                       4 : the drift contribution only
 ** \li                       5 : the external drift contribution only
 **
 ** \remark The results will be stored in the dbout file
 ** \remark - the estimation in the variable ELoc::Z
 ** \remark - the gradient components in the variables ELoc::GRD
 **
 *****************************************************************************/
Id Potential::kriging(DbGrid* dbout,
                      bool flag_pot,
                      bool flag_grad,
                      bool flag_trans,
                      bool flag_save_data,
                      Id option_part)
{
  VectorInt uid_iso_pot, uid_iso_grad;
  VectorInt uid_grd_pot, uid_grd_grad;
  VectorInt uid_tgt_pot, uid_tgt_grad;
  VectorInt uid_out_pot, uid_out_grad;

  // Initialization
  _environmentManage(flag_pot, flag_grad, flag_trans, option_part);
  if (!_isEnvironmentValid(dbout)) return 1;

  // Count the gradients and the tangents
  if (_updateIsopot()) return 1;
  if (_updateGradient()) return 1;
  if (_updateTangent()) return 1;
  if (_updateModel()) return 1;
  if (_updateFinal()) return 1;

  // Allocating the output variables

  _saveResultData(dbout, 1, TEST, ELoc::Z, ELoc::G, uid_out_pot, uid_out_grad);
  if (flag_save_data)
  {
    _saveResultData(_dbiso, 1, TEST, ELoc::UNKNOWN, ELoc::UNKNOWN, uid_iso_pot, uid_iso_grad);
    _saveResultData(_dbgrd, 1, TEST, ELoc::UNKNOWN, ELoc::UNKNOWN, uid_grd_pot, uid_grd_grad);
    _saveResultData(_dbtgt, 1, TEST, ELoc::UNKNOWN, ELoc::UNKNOWN, uid_tgt_pot, uid_tgt_grad);
  }

  // Core allocation
  MatrixDense rhs(_nequa, 4);
  VectorDouble potval(_nlayers);

  // Establish the cokriging system
  MatrixSymmetric lhs;
  if (_buildLHS(dbout, lhs)) return 1;

  // Invert the matrix
  if (lhs.invert()) return 1;
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    print_matrix("[LHS]-1", 0, 1, _nequa, _nequa, NULL, lhs.getValues().data());

  // Establish the data vector and get the dual form
  VectorDouble zval;
  _fillDual(zval);
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    print_matrix("\n[Z]", 0, 1, 1, _nequa, NULL, zval.data());
  VectorDouble zdual(_nequa);
  lhs.prodMatVecInPlace(zval, zdual);
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    print_matrix("\n[Z] * [LHS]-1", 0, 1, 1, _nequa, NULL, zdual.data());

  // Evaluate Potential at Reference point
  double refpot = _evaluateRefPot(dbout, zdual, rhs);

  // Check that the information is fulfilled correctly
  if (_verbose)
    _checkData(dbout, -1, 0, refpot, zdual, rhs);

  // Get the Potential value at the iso-potential samples
  _evaluatePotential(dbout, refpot, -1, 0, zdual, rhs, potval.data());

  // Perform the estimation
  _estimateResult(flag_grad, dbout, refpot, zdual, rhs, potval.data());
  if (flag_save_data)
  {
    _estimateData(dbout, refpot, zdual, rhs, _dbiso, uid_iso_pot, uid_iso_grad);
    _estimateData(dbout, refpot, zdual, rhs, _dbgrd, uid_grd_pot, uid_grd_grad);
    _estimateData(dbout, refpot, zdual, rhs, _dbtgt, uid_tgt_pot, uid_tgt_grad);
  }

  return 0;
}

/****************************************************************************/
/*!
 **  Potential simulations
 **
 ** \return  Error return code
 **
 ** \param[in]  dbout        Output Db structure
 ** \param[in]  dist_tempere Distance for tempering simulations (or TEST)
 ** \param[in]  flag_trans   True if the estimation result must be translated
 **                          into layer number
 ** \param[in]  seed         Seed for the random number generator
 ** \param[in]  nbsimu       Number of simulations
 ** \param[in]  nbtuba       Number of turning bands
 **
 ** \remark The simulations will be stored in the dbout file (ELoc::SIMU)
 **
 *****************************************************************************/
Id Potential::simulate(DbGrid* dbout,
                       double dist_tempere,
                       bool flag_trans,
                       Id seed,
                       Id nbsimu,
                       Id nbtuba)
{
  VectorInt uid_out_pot, uid_out_grad;
  VectorInt uid_iso_pot, uid_iso_grad;
  VectorInt uid_grd_pot, uid_grd_grad;
  VectorInt uid_tgt_pot, uid_tgt_grad;

  // Initialization
  _environmentManage(true, false, flag_trans, 0);
  if (!_isEnvironmentValid(dbout)) return 1;

  law_set_random_seed(seed);
  bool flag_tempere = !FFFF(dist_tempere);
  double refpot     = 0.;
  double delta      = _dbiso->getExtensionDiagonal() / 1000;

  // Count the gradients and the tangents
  if (_updateIsopot()) return 1;
  if (_updateGradient()) return 1;
  if (_updateTangent()) return 1;
  if (_updateModel()) return 1;
  if (_updateFinal()) return 1;

  /* Add the attributes for storing the results */
  _saveResultData(dbout, nbsimu, 0., ELoc::SIMU, ELoc::UNKNOWN,
                  uid_out_pot, uid_out_grad);
  _saveResultData(_dbiso, 2 * nbsimu, 0., ELoc::SIMU, ELoc::UNKNOWN,
                  uid_iso_pot, uid_iso_grad);
  _saveResultData(_dbgrd, 2 * nbsimu, 0., ELoc::SIMU, ELoc::UNKNOWN,
                  uid_grd_pot, uid_grd_grad);
  _saveResultData(_dbtgt, 2 * nbsimu, 0., ELoc::SIMU, ELoc::UNKNOWN,
                  uid_tgt_pot, uid_tgt_grad);
  if (flag_tempere)
    (void)dbout->addColumnsByConstant(1, TEST, String(), ELoc::Z);

  /* Processing the non-conditional simulation over the iso-values */
  {
    gstlrn::CalcSimuTurningBands situba_new(nbsimu, nbtuba, seed);
    if (situba_new.simulatePotential(_dbiso, _dbgrd, _dbtgt, dbout, _model, delta))
      return 1;
  }

  // Core allocation
  VectorDouble zval(_nequa);
  VectorDouble potsim(_nlayers * nbsimu);
  VectorDouble potval(_nlayers);
  MatrixDense zvals(_nequa, nbsimu);
  MatrixDense zduals(_nequa, nbsimu);
  MatrixDense rhs(_nequa, 4);
  VectorDouble zdual;
  if (flag_tempere) zdual.resize(_nequa);

  // Establish the cokriging system
  MatrixSymmetric lhs;
  if (_buildLHS(dbout, lhs)) return 1;

  // Invert the matrix
  if (lhs.invert()) return 1;
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    print_matrix("Inverted LHS", 0, 1, _nequa, _nequa, NULL, lhs.getValues().data());

  if (flag_tempere)
  {

    // Establish the data vector and get the dual form
    _fillDual(zval);
    if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
      print_matrix("\n[Z]", 0, 1, 1, _nequa, NULL, zval.data());
    lhs.prodMatVecInPlace(zval, zdual);
    if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
      print_matrix("\n[Z] *%* [A]-1", 0, 1, 1, _nequa, NULL, zdual.data());

    // Evaluate Potential at Reference point
    refpot = _evaluateRefPot(dbout, zdual, rhs);

    // Get the Estimated Potential value at the iso-potential samples
    _evaluatePotential(dbout, refpot, -1, 0, zdual, rhs, potval.data());

    // Perform the estimation
    _estimateResult(false, dbout, refpot, zdual, rhs, potval.data());

    // Transform the Estimation variable into a distance ,
    if (_distanceToIsoline(dbout)) return 1;
  }

  // Establish the simulated error vector and get the dual form
  _fillDualSimulation(nbsimu, zvals);
  lhs.prodMatMatInPlace(&zvals, &zduals);
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    print_matrix("\n[Simu-Err] *%* [A]-1", 0, 1, nbsimu, _nequa, NULL, zduals.getValues().data());

  // Get the Simulated Potential value at the iso-potential samples
  for (Id isimu = 0; isimu < nbsimu; isimu++)
  {
    VectorDouble zdual_loc = zduals.getColumn(isimu);

    // Calculate the simulated value at the reference point
    refpot = _evaluateRefPot(dbout, zdual_loc, rhs);

    // Check that the information is fulfilled correctly
    _checkData(dbout, isimu, nbsimu, refpot, zdual_loc, rhs);

    // Calculate the simulated iso-value
    _evaluatePotential(dbout, refpot, isimu, nbsimu, zdual_loc, rhs,
                       &potsim[isimu * _nlayers]);
  }

  // Perform the conditional simulations on the grid
  _simcond(dist_tempere, flag_trans, nbsimu, dbout, refpot, potsim.data(), zdual, zduals, rhs);

  if (flag_tempere) dbout->deleteColumnsByLocator(ELoc::Z);
  return 0;
}

/****************************************************************************/
/*!
 **  Potential cross-validation
 **
 ** \return  Error return code
 **
 ** \param[in]  flag_dist_conv Flag for converting into distance
 **
 *****************************************************************************/
Id Potential::xvalid(bool flag_dist_conv)
{
  // Initialization
  _environmentManage(true, false, false, 0);
  if (!_isEnvironmentValid(nullptr)) return 1;

  // Count the gradients and the tangents
  if (_updateIsopot()) return 1;
  if (_updateGradient()) return 1;
  if (_updateTangent()) return 1;
  if (_updateModel()) return 1;
  if (_updateFinal()) return 1;

  // Allocating the output variables
  int nvar = (flag_dist_conv) ? 4 : 2;
  (void)_dbiso->addColumnsByConstant(nvar, TEST, String(), ELoc::Z);

  // Core allocation
  VectorDouble zval(_nequa);
  VectorDouble zdual(_nequa);
  MatrixDense rhs(_nequa, 4);
  MatrixSymmetric lhs_orig;
  MatrixSymmetric lhs_aux;
  if (flag_dist_conv)
  {
    lhs_orig.resize(_nequa, _nequa);
    lhs_aux.resize(_nequa, _nequa);
  }

  // Establish the cokriging system
  MatrixSymmetric lhs;
  if (_buildLHS(nullptr, lhs)) return 1;

  // Save the matrix (used for converting into distance)
  if (flag_dist_conv) lhs_orig = lhs;

  // Invert the matrix
  if (lhs.invert()) return 1;
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    print_matrix("Inverted LHS", 0, 1, _nequa, _nequa, NULL, lhs.getValues().data());

  // Establish the data vector and get the dual form
  _fillDual(zval);
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    print_matrix("\n[Z]", 0, 1, 1, _nequa, NULL, zval.data());
  lhs.prodMatVecInPlace(zval, zdual);
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    print_matrix("\n[Z] *%* [A]-1", 0, 1, 1, _nequa, NULL, zdual.data());

  /* Process the estimate at masked-off isovalues */
  _xvalidCalculate(lhs, flag_dist_conv, zval, lhs_orig, rhs);

  return 0;
}

/****************************************************************************/
/*!
 **  Potential estimation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbiso      Iso-potential Db structure
 ** \param[in]  dbgrd      Gradient Db structure
 ** \param[in]  dbtgt      Tangent Db structure (optional)
 ** \param[in]  dbout      Output Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  nugget_grd Nugget effect for Gradients
 ** \param[in]  nugget_tgt Nugget effect for Tangents
 ** \param[in]  flag_pot   True if the Potential must be estimated
 ** \param[in]  flag_grad  True if the gradient must also be estimated
 ** \param[in]  flag_trans True if the estimation result must be translated
 **                        into layer number
 ** \param[in]  flag_save_data True if the Potential / Gradient must be
 **                        saved on any Information file
 ** \param[in]  option_part   Option to exhibit only a part of estimation:
 ** \li                       0 : the whole estimation
 ** \li                       1 : the gradient contribution only
 ** \li                       2 : the tangent contribution only
 ** \li                       3 : the isovalues contribution only
 ** \li                       4 : the drift contribution only
 ** \li                       5 : the external drift contribution only
 ** \param[in]  verbose    Verbose option
 **
 ** \remark The results will be stored in the dbout file
 ** \remark - the estimation in the variable ELoc::Z
 ** \remark - the gradient components in the variables ELoc::GRD
 **
 *****************************************************************************/
Id krigingPotential(Db* dbiso,
                    Db* dbgrd,
                    Db* dbtgt,
                    DbGrid* dbout,
                    ModelGeneric* model,
                    double nugget_grd,
                    double nugget_tgt,
                    bool flag_pot,
                    bool flag_grad,
                    bool flag_trans,
                    bool flag_save_data,
                    Id option_part,
                    bool verbose)
{
  Potential pot(dbiso, dbgrd, dbtgt, model, nugget_grd, nugget_tgt);
  pot.setVerbose(verbose);
  Id error = pot.kriging(dbout, flag_pot, flag_grad, flag_trans, flag_save_data, option_part);

  return error;
}

/****************************************************************************/
/*!
 **  Potential simulations
 **
 ** \return  Error return code
 **
 ** \param[in]  dbiso        Iso-potential Db structure
 ** \param[in]  dbgrd        Gradient Db structure
 ** \param[in]  dbtgt        Tangent Db structure (optional)
 ** \param[in]  dbout        Output Db structure
 ** \param[in]  model        Model structure
 ** \param[in]  nugget_grd   Nugget effect for Gradients
 ** \param[in]  nugget_tgt   Nugget effect for Tangents
 ** \param[in]  dist_tempere Distance for tempering simulations (or TEST)
 ** \param[in]  flag_trans   True if the estimation result must be translated
 **                          into layer number
 ** \param[in]  seed         Seed for the random number generator
 ** \param[in]  nbsimu       Number of simulations
 ** \param[in]  nbtuba       Number of turning bands
 ** \param[in]  verbose      Verbose option
 **
 ** \remark The simulations will be stored in the dbout file (ELoc::SIMU)
 **
 *****************************************************************************/
Id simulatePotential(Db* dbiso,
                     Db* dbgrd,
                     Db* dbtgt,
                     DbGrid* dbout,
                     ModelGeneric* model,
                     double nugget_grd,
                     double nugget_tgt,
                     double dist_tempere,
                     bool flag_trans,
                     Id seed,
                     Id nbsimu,
                     Id nbtuba,
                     bool verbose)
{
  Potential pot(dbiso, dbgrd, dbtgt, model, nugget_grd, nugget_tgt);
  pot.setVerbose(verbose);
  Id error = pot.simulate(dbout, dist_tempere, flag_trans, seed, nbsimu, nbtuba);

  return error;
}

/****************************************************************************/
/*!
 **  Potential cross-validation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbiso          Iso-potential Db structure
 ** \param[in]  dbgrd          Gradient Db structure
 ** \param[in]  dbtgt          Tangent Db structure (optional)
 ** \param[in]  model          Model structure
 ** \param[in]  nugget_grd     Nugget effect for Gradients
 ** \param[in]  nugget_tgt     Nugget effect for Tangents
 ** \param[in]  flag_dist_conv Flag for converting into distance
 ** \param[in]  verbose        Verbose option
 **
 *****************************************************************************/
Id xvalidPotential(Db* dbiso,
                   Db* dbgrd,
                   Db* dbtgt,
                   ModelGeneric* model,
                   double nugget_grd,
                   double nugget_tgt,
                   bool flag_dist_conv,
                   bool verbose)
{
  Potential pot(dbiso, dbgrd, dbtgt, model, nugget_grd, nugget_tgt);
  pot.setVerbose(verbose);
  Id error = pot.xvalid(flag_dist_conv);

  return error;
}

} // namespace gstlrn