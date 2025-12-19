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
#include "Variogram/AVario.hpp"

#include "Basic/String.hpp"
#include "Basic/VectorNumT.hpp"
#include "Db/Db.hpp"
#include "Enum/ECalcVario.hpp"

namespace gstlrn
{
AVario::AVario()
  : AStringable()
  , _calcul(ECalcVario::UNDEFINED)
  , _flagStorage(false)
  , _tabStorage()
{
}

AVario::AVario(const AVario& r)
  : AStringable(r)
  , _calcul(r._calcul)
  , _flagStorage(r._flagStorage)
  , _tabStorage(r._tabStorage)
{
}

AVario& AVario::operator=(const AVario& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _calcul      = r._calcul;
    _flagStorage = r._flagStorage;
    _tabStorage  = r._tabStorage;
  }
  return *this;
}

AVario::~AVario() {}

void AVario::_evaluateVariogram(
  Db* db,
  Id nvar,
  Id iech1,
  Id iech2,
  Id ilag,
  double dist,
  bool do_asym)
{
  DECLARE_UNUSED(do_asym);
  double w1 = db->getWeight(iech1);
  double w2 = db->getWeight(iech2);
  if (FFFF(w1) || FFFF(w2)) return;
  dist         = ABS(dist);
  double scale = w1 * w2;
  for (Id ivar = 0; ivar < nvar; ivar++)
  {
    double z11 = _getIVAR(db, iech1, ivar);
    double z12 = _getIVAR(db, iech2, ivar);
    if (FFFF(z11) || FFFF(z12)) continue;
    for (Id jvar = 0; jvar <= ivar; jvar++)
    {
      double z21 = _getIVAR(db, iech1, jvar);
      double z22 = _getIVAR(db, iech2, jvar);
      if (FFFF(z21) || FFFF(z22)) continue;
      double value = (z12 - z11) * (z22 - z21) / 2.;
      _setResult(iech1, iech2, nvar, ilag, ivar, jvar, 0, scale, dist, value);
    }
  }
}

void AVario::_evaluateMadogram(
  Db* db,
  Id nvar,
  Id iech1,
  Id iech2,
  Id ilag,
  double dist,
  bool do_asym)
{
  DECLARE_UNUSED(do_asym);
  double w1 = db->getWeight(iech1);
  double w2 = db->getWeight(iech2);
  if (FFFF(w1) || FFFF(w2)) return;
  dist         = ABS(dist);
  double scale = w1 * w2;
  for (Id ivar = 0; ivar < nvar; ivar++)
  {
    double z11 = _getIVAR(db, iech1, ivar);
    double z12 = _getIVAR(db, iech2, ivar);
    if (FFFF(z11) || FFFF(z12)) continue;
    for (Id jvar = 0; jvar <= ivar; jvar++)
    {
      double z21 = _getIVAR(db, iech1, jvar);
      double z22 = _getIVAR(db, iech2, jvar);
      if (FFFF(z21) || FFFF(z22)) continue;
      double value = sqrt(ABS((z12 - z11) * (z22 - z21))) / 2.;
      _setResult(iech1, iech2, nvar, ilag, ivar, jvar, 0, scale, dist, value);
    }
  }
}

void AVario::_evaluateRodogram(
  Db* db,
  Id nvar,
  Id iech1,
  Id iech2,
  Id ilag,
  double dist,
  bool do_asym)
{
  DECLARE_UNUSED(do_asym);
  double w1 = db->getWeight(iech1);
  double w2 = db->getWeight(iech2);
  if (FFFF(w1) || FFFF(w2)) return;
  dist         = ABS(dist);
  double scale = w1 * w2;
  for (Id ivar = 0; ivar < nvar; ivar++)
  {
    double z11 = _getIVAR(db, iech1, ivar);
    double z12 = _getIVAR(db, iech2, ivar);
    if (FFFF(z11) || FFFF(z12)) continue;
    for (Id jvar = 0; jvar <= ivar; jvar++)
    {
      double z21 = _getIVAR(db, iech1, jvar);
      double z22 = _getIVAR(db, iech2, jvar);
      if (FFFF(z21) || FFFF(z22)) continue;
      double value = pow(ABS((z12 - z11) * (z22 - z21)), 0.25) / 2.;
      _setResult(iech1, iech2, nvar, ilag, ivar, jvar, 0, scale, dist, value);
    }
  }
}

void AVario::_evaluatePoisson(
  Db* db,
  Id nvar,
  Id iech1,
  Id iech2,
  Id ilag,
  double dist,
  bool do_asym)
{
  DECLARE_UNUSED(do_asym);
  double w1 = db->getWeight(iech1);
  double w2 = db->getWeight(iech2);
  if (FFFF(w1) || FFFF(w2)) return;
  dist         = ABS(dist);
  double scale = (w1 * w2) / (w1 + w2);
  for (Id ivar = 0; ivar < nvar; ivar++)
  {
    double z11 = _getIVAR(db, iech1, ivar);
    double z12 = _getIVAR(db, iech2, ivar);
    if (FFFF(z11) || FFFF(z12)) continue;
    for (Id jvar = 0; jvar <= ivar; jvar++)
    {
      double z21 = _getIVAR(db, iech1, jvar);
      double z22 = _getIVAR(db, iech2, jvar);
      if (FFFF(z21) || FFFF(z22)) continue;
      double value = (z12 - z11) * (z22 - z21) / 2.;
      _setResult(iech1, iech2, nvar, ilag, ivar, jvar, 0, scale, dist, value);
    }
  }
}

void AVario::_evaluateCovariance(
  Db* db,
  Id nvar,
  Id iech1,
  Id iech2,
  Id ilag,
  double dist,
  bool do_asym)
{
  double w1 = db->getWeight(iech1);
  double w2 = db->getWeight(iech2);
  if (FFFF(w1) || FFFF(w2)) return;
  Id orient    = (dist > 0) ? 1 : -1;
  dist         = ABS(dist);
  double scale = w1 * w2;
  for (Id ivar = 0; ivar < nvar; ivar++)
  {
    double z11 = _getIVAR(db, iech1, ivar);
    double z12 = _getIVAR(db, iech2, ivar);
    if (FFFF(z11) || FFFF(z12)) continue;
    for (Id jvar = 0; jvar <= ivar; jvar++)
    {
      double z21 = _getIVAR(db, iech1, jvar);
      double z22 = _getIVAR(db, iech2, jvar);
      if (!FFFF(z22))
      {
        double value = z11 * z22;
        _setResult(iech1, iech2, nvar, ilag, ivar, jvar, orient, scale, dist, value);
      }
      if (!FFFF(z21) && do_asym)
      {
        double value = z12 * z21;
        _setResult(iech1, iech2, nvar, ilag, ivar, jvar, -orient, scale, dist, value);
      }
    }
  }
}

void AVario::_evaluateCovariogram(
  Db* db,
  Id nvar,
  Id iech1,
  Id iech2,
  Id ilag,
  double dist,
  bool do_asym)
{
  double w1 = db->getWeight(iech1);
  double w2 = db->getWeight(iech2);
  if (FFFF(w1) || FFFF(w2)) return;
  Id orient    = (dist > 0) ? 1 : -1;
  dist         = ABS(dist);
  double scale = w2;
  for (Id ivar = 0; ivar < nvar; ivar++)
  {
    double z11 = _getIVAR(db, iech1, ivar);
    double z12 = _getIVAR(db, iech2, ivar);
    if (FFFF(z11) || FFFF(z12)) continue;
    for (Id jvar = 0; jvar <= ivar; jvar++)
    {
      double z21 = _getIVAR(db, iech1, jvar);
      double z22 = _getIVAR(db, iech2, jvar);
      if (!FFFF(z22))
      {
        double value = z11 * z22;
        _setResult(iech1, iech2, nvar, ilag, ivar, jvar, orient, scale, dist, value);
      }
      if (!FFFF(z21) && do_asym)
      {
        double value = z12 * z21;
        _setResult(iech1, iech2, nvar, ilag, ivar, jvar, -orient, scale, dist, value);
      }
    }
  }
}

void AVario::_evaluateOrder4(
  Db* db,
  Id nvar,
  Id iech1,
  Id iech2,
  Id ilag,
  double dist,
  bool do_asym)
{
  DECLARE_UNUSED(do_asym);
  double w1 = db->getWeight(iech1);
  double w2 = db->getWeight(iech2);
  if (FFFF(w1) || FFFF(w2)) return;
  dist         = ABS(dist);
  double scale = w1 * w2;
  for (Id ivar = 0; ivar < nvar; ivar++)
  {
    double z11 = _getIVAR(db, iech1, ivar);
    double z12 = _getIVAR(db, iech2, ivar);
    if (FFFF(z11) || FFFF(z12)) continue;
    for (Id jvar = 0; jvar <= ivar; jvar++)
    {
      double z21 = _getIVAR(db, iech1, jvar);
      double z22 = _getIVAR(db, iech2, jvar);
      if (FFFF(z21) || FFFF(z22)) continue;
      double value = (z12 - z11) * (z22 - z21);
      value        = value * value / 2.;
      _setResult(iech1, iech2, nvar, ilag, ivar, jvar, 0, scale, dist, value);
    }
  }
}

String AVario::_elemString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);
  std::stringstream sstr;
  if (getCalcul() == ECalcVario::UNDEFINED)
    sstr << toStrTitle(0, "Undefined");
  else
  {
    String total_string = std::string(getCalcul().getDescr()) + " characteristics";
    sstr << toStrTitle(0, total_string.c_str());
  }
  return sstr.str();
}

void AVario::setCalcul(const ECalcVario& calcul)
{
  _calcul = calcul;

  // Define the internal evaluation function
  switch (_calcul.toEnum())
  {
    case ECalcVario::E_VARIOGRAM:
    case ECalcVario::E_TRANS1:
    case ECalcVario::E_TRANS2:
    case ECalcVario::E_BINORMAL:
    {
      _evaluate = &AVario::_evaluateVariogram;
      break;
    }

    case ECalcVario::E_MADOGRAM:
    {
      _evaluate = &AVario::_evaluateMadogram;
      break;
    }

    case ECalcVario::E_RODOGRAM:
    {
      _evaluate = &AVario::_evaluateRodogram;
      break;
    }

    case ECalcVario::E_POISSON:
    {
      _evaluate = &AVario::_evaluatePoisson;
      break;
    }

    case ECalcVario::E_COVARIANCE:
    case ECalcVario::E_COVARIANCE_NC:
    {
      _evaluate = &AVario::_evaluateCovariance;
      break;
    }

    case ECalcVario::E_COVARIOGRAM:
    {
      _evaluate = &AVario::_evaluateCovariogram;
      break;
    }

    case ECalcVario::E_ORDER4:
    {
      _evaluate = &AVario::_evaluateOrder4;
      break;
    }

    default:
    {
      messageAbort("AVario::evaluate() ignores current calculation type");
      break;
    }
  }
}

void AVario::setStorage(bool flag)
{
  if (flag == _flagStorage) return;

  if (flag)
  {
    // Activate storage
    _tabStorage.clear();
    _flagStorage = true;
  }
  else
  {
    // Deactivate storage
    _tabStorage.clear();
    _flagStorage = false;
  }
}

void AVario::_storage(Id iech1, Id iech2, double dist, double value)
{
  if (!_flagStorage) return;

  VectorDouble record(4);
  record[0] = static_cast<double>(iech1);
  record[1] = static_cast<double>(iech2);
  record[2] = dist;
  record[3] = value;
  _tabStorage.push_back(record);
}

} // namespace gstlrn
