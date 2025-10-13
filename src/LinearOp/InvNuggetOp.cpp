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

#include "LinearOp/InvNuggetOp.hpp"
#include "Basic/AStringable.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Model/Model.hpp"
#include "geoslib_define.h"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace gstlrn
{

InvNuggetOp::InvNuggetOp(const Db* dbin, Model* model, const SPDEParam& params, bool flagEigVals)
  : ASimulable()
  , MatrixSparse()
  , _logDeterminant(TEST)
  , _rangeEigenVal(INF, 0.)
  , _flagEigVals(flagEigVals)
{
  _buildInvNugget(dbin, model, params);
}

Id InvNuggetOp::getSize() const
{
  return getNRows();
}

InvNuggetOp::~InvNuggetOp()
{
}

Id InvNuggetOp::_addSimulateToDest(const constvect whitenoise, vect outv) const
{
  return _cholNuggetMatrix->addToDest(whitenoise, outv);
}

static void _addVerrConstant(MatrixSymmetric& sills, const VectorDouble& verrDef)
{
  Id nverr = static_cast<Id>(verrDef.size());
  if (nverr > 0)
  {
    for (Id iverr = 0; iverr < nverr; iverr++)
      sills.updValue(iverr, iverr, EOperator::ADD, verrDef[iverr]);
  }
}

static void _checkMinNugget(MatrixSymmetric& sills, const VectorDouble& minNug)
{
  Id nvar = static_cast<Id>(minNug.size());

  // Check that the diagonal of the Sill matrix is large enough
  for (Id ivar = 0; ivar < nvar; ivar++)
    sills.setValue(ivar, ivar, MAX(sills.getValue(ivar, ivar), minNug[ivar]));
}

static MatrixSymmetric _buildSillPartialMatrix(const MatrixSymmetric& sillsRef,
                                               Id nvar,
                                               Id ndef,
                                               const VectorInt& identity)
{
  MatrixSymmetric sills;
  if (ndef == nvar)
    sills = sillsRef;
  else
  {
    sills = MatrixSymmetric(ndef);
    for (Id idef = 0; idef < ndef; idef++)
      for (Id jdef = 0; jdef <= idef; jdef++)
        sills.setValue(idef, jdef, sillsRef.getValue(identity[idef], identity[jdef]));
  }
  return sills;
}

static Id _loadPositions(Id iech,
                         const VectorVectorInt& index1,
                         const VectorInt& cumul,
                         VectorInt& positions,
                         VectorInt& identity,
                         Id* rank_arg)
{
  Id nvar = static_cast<Id>(cumul.size());
  Id ndef = 0;
  Id rank = 0;
  for (Id ivar = 0; ivar < nvar; ivar++)
  {
    rank    = 2 * rank;
    Id ipos = VH::whereElement(index1[ivar], iech);
    if (ipos < 0)
      positions[ivar] = -1;
    else
    {
      positions[ivar] = ipos + cumul[ivar];
      identity[ndef]  = ivar;
      ndef++;
      rank += 1;
    }
  }
  *rank_arg = rank;
  return ndef;
}

/**
 * Compute the log det and update the range of eigenvalues
 *
 * @param sillsinv The inverse of the Nugget Effect matrix
 * @return The log determinant of the inverse Nugget Effect matrix
 */

double InvNuggetOp::_updateQuantities(MatrixSymmetric& sillsinv)
{
  sillsinv.computeEigen();
  auto eigenvals        = sillsinv.getEigenValues();
  auto rangevals        = std::minmax_element(eigenvals.begin(), eigenvals.end());
  _rangeEigenVal.first  = MIN(_rangeEigenVal.first, *rangevals.first);
  _rangeEigenVal.second = MAX(_rangeEigenVal.second, *rangevals.second);
  std::transform(eigenvals.begin(), eigenvals.end(), eigenvals.begin(), [](double x)
                 { return std::log(x); });
  return std::accumulate(eigenvals.begin(), eigenvals.end(), 0.);
}

/**
 * Build the inverse of the Nugget Effect matrix
 * It is established for:
 * - the number of variables defined in 'dbin' (and in 'Model')
 * - the active samples of 'dbin'
 * - the samples where Z-variable (and possibly V-variable) is defined
 *
 * @param db Input Db structure
 * @param model Input Model structure
 * @param params A structure for ruling the parameters of SPDE
 */
void InvNuggetOp::_buildInvNugget(const Db* db, Model* model, const SPDEParam& params)
{
  CholeskyDense chol;
  if (db == nullptr) return;
  Id nech = db->getNSample();
  if (model == nullptr) return;
  Id nvar = db->getNLoc(ELoc::Z);
  if (nvar != model->getNVar())
  {
    messerr("'db' and 'model' should have the same number of variables");
    return;
  }
  bool hasnugget = false;
  CovAniso* cova = nullptr;

  for (Id icov = 0; icov < model->getNCov(); icov++)
  {
    if (model->getCovAniso(icov)->getType() == ECov::NUGGET)
    {
      cova      = model->getCovAniso(icov);
      hasnugget = true;
      break;
    }
  }
  // If no Nugget Effect is defined, create a nugget to
  // integrate potential measurement error or minimum nugget
  if (!hasnugget)
  {
    MatrixSymmetric sills(model->getNVar());
    cova = CovAniso::createIsotropicMulti(*model->getContext(), ECov::NUGGET, 0, sills);
  }
  VectorInt ivars = VH::sequence(nvar);

  // Get the minimum value for diagonal terms
  double eps = params.getEpsNugget();
  VectorDouble minNug(nvar);
  for (Id ivar = 0; ivar < nvar; ivar++)
    minNug[ivar] = eps * model->getTotalSill(ivar, ivar);

  // Play the non-stationarity (if needed)
  bool flag_nostat_sill = cova->isNoStatForVariance();
  if (flag_nostat_sill)
    cova->informDbInForSills(db);

  // Create sets of Vector of valid sample indices per variable (not masked and defined)
  VectorVectorInt index1 = db->getSampleRanks(ivars);
  // 'cumul' counts the number of valid positions for all variables before 'ivar'
  VectorInt cumul = VH::cumulIncrement(index1);

  Id sizetot = VH::count(index1);
  // Convert from triplet to sparse matrix
  _allocate(sizetot, sizetot, nvar);
  _cholNuggetMatrix = std::make_shared<MatrixSparse>(sizetot, sizetot, nvar);

  // Check the various possibilities
  // - flag_verr: True if Variance of Measurement Error variable is defined
  // - flag_isotropic: True in Isotopic case
  // - flag_uniqueVerr: True if the Variance of Measurement Error is constant per variable
  // - flag_nostat: True is some non-stationarity is defined
  Id nverr           = db->getNLoc(ELoc::V);
  bool flag_verr     = (nverr > 0);
  bool flag_isotopic = true;
  for (Id ivar = 1; ivar < nvar && flag_isotopic; ivar++)
    if (!VH::isEqual(index1[ivar], index1[0])) flag_isotopic = false;
  bool flag_uniqueVerr = true;
  VectorDouble verrDef(nverr, 0.);
  if (flag_verr)
  {
    for (Id iverr = 0; iverr < nverr && flag_uniqueVerr; iverr++)
    {
      VectorDouble verr = db->getColumnByLocator(ELoc::V, iverr);
      if (static_cast<Id>(VH::unique(verr).size()) > 1) flag_uniqueVerr = false;
      verrDef[iverr] = verr[0];
    }
  }
  bool flag_constant = (!flag_nostat_sill && (!flag_verr || flag_uniqueVerr));

  // Elaborate the Sill matrix for the Nugget Effect component
  MatrixSymmetric sillsRef = cova->getSill();
  Id count                 = static_cast<Id>(pow(2, nvar));
  std::vector<MatrixSymmetric> sillsInv(count);

  // Pre-calculate the inverse of the sill matrix (if constant)

  std::vector<double> cacheLogDet;
  std::vector<CholeskyDense> cholcache;

  if (flag_constant)
  {
    cacheLogDet.resize(count, 0.);
    cholcache.resize(count);
    // In case of (Unique) Variance of measurement error, patch sill matrix
    if (flag_verr) _addVerrConstant(sillsRef, verrDef);

    // Check that the diagonal of the Sill matrix is large enough
    _checkMinNugget(sillsRef, minNug);
  }

  // Loop on the samples
  Id rank;
  Id ndef = nvar;
  VectorInt position(nvar);
  VectorInt identity(nvar);
  for (Id iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;

    // Count the number of variables for which current sample is valid
    ndef = _loadPositions(iech, index1, cumul, position, identity, &rank);
    if (ndef <= 0) continue;

    // If all samples are defined, in the stationary case, use the inverted sill matrix
    if (flag_constant)
    {
      if (sillsInv[rank].empty())
      {
        sillsInv[rank] = _buildSillPartialMatrix(sillsRef, nvar, ndef, identity);
        cholcache[rank].setMatrix(sillsInv[rank]);

        if (sillsInv[rank].invert() != 0) return;
        if (_flagEigVals)
        {
          if (FFFF(_logDeterminant))
            _logDeterminant = 0.;

          cacheLogDet[rank] = _updateQuantities(sillsInv[rank]);
        }
      }
      if (_flagEigVals)
      {
        _logDeterminant += cacheLogDet[rank];
      }
      _updateMatrix(sillsInv[rank], cholcache[rank], ndef, position, identity);
    }
    else
    {
      // Update due to non-stationarity (optional)
      if (flag_nostat_sill)
      {
        cova->updateCovByPoints(1, iech, 1, iech);
        sillsRef = cova->getSill();
      }

      // Establish a local matrix
      MatrixSymmetric local(ndef);
      for (Id idef = 0; idef < ndef; idef++)
        for (Id jdef = 0; jdef <= idef; jdef++)
        {
          // Load the sill value of the Nugget Effect component
          double value = sillsRef.getValue(identity[idef], identity[jdef]);

          // Patch the diagonal term of the local matrix
          if (idef == jdef)
          {
            // Add the Variance of measurement error (optional)
            if (flag_verr && idef < nverr)
              value += db->getFromLocator(ELoc::V, iech, identity[idef]);

            // Check the minimum values over the diagonal
            value = MAX(value, MAX(local.getValue(idef, idef), minNug[idef]));
          }

          local.setValue(idef, jdef, value);
        }
      CholeskyDense chollocal(local);
      if (local.invert() != 0) return;
      if (_flagEigVals)
      {
        if (FFFF(_logDeterminant))
          _logDeterminant = 0.;
        _logDeterminant += _updateQuantities(local);
      }

      _updateMatrix(local, chollocal, ndef, position, identity);
    }
  }
  // Free the non-stationary specific allocation
  if (!hasnugget)
    delete cova;
}
void InvNuggetOp::_updateMatrix(MatrixSymmetric& invsill, CholeskyDense& cholsill, Id ndef, const VectorInt& position, const VectorInt& identity)
{
  {
    for (Id idef = 0; idef < ndef; idef++)
      for (Id jdef = 0; jdef < ndef; jdef++)
      {
        Id pidef = position[identity[idef]];
        Id pjdef = position[identity[jdef]];
        setValue(pidef, pjdef, invsill.getValue(idef, jdef));
        if (idef >= jdef)
        {
          _cholNuggetMatrix->setValue(pidef, pjdef, cholsill.getLowerTriangle(idef, jdef));
        }
      }
  }
}

std::shared_ptr<MatrixSparse> buildInvNugget(Db* dbin, Model* model, const SPDEParam& params)
{
  InvNuggetOp invnugg(dbin, model, params);
  return std::make_shared<MatrixSparse>(*invnugg.getInvNuggetMatrix());
}

double InvNuggetOp::computeLogDet(Id nMC) const
{
  DECLARE_UNUSED(nMC)

  if (FFFF(_logDeterminant))
  {
    CholeskyDense cholesky(*this);
    return cholesky.computeLogDeterminant();
  }
  return _logDeterminant;
}

const MatrixSparse* InvNuggetOp::cloneInvNuggetMatrix() const
{
  return this->clone();
}

const MatrixSparse* InvNuggetOp::cloneCholNuggetMatrix() const
{
  return _cholNuggetMatrix->clone();
}
} // namespace gstlrn