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

#include "Enum/ECalcVario.hpp"
#include "Db/Db.hpp"
#include "Basic/Utilities.hpp"

AVario::AVario()
  : AStringable()
  , _calcul(ECalcVario::UNDEFINED)
{
}

AVario::AVario(const AVario& r)
  : AStringable(r)
  , _calcul(r._calcul)
{
}

AVario& AVario::operator=(const AVario& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _calcul = r._calcul;
  }
  return *this;
}

AVario::~AVario() {}

void AVario::_evaluateVariogram(
  Db* db, int nvar, int iech1, int iech2, int ipas, double dist, bool /*do_asym*/)
{
  double w1 = db->getWeight(iech1);
  double w2 = db->getWeight(iech2);
  if (FFFF(w1) || FFFF(w2)) return;
  dist         = ABS(dist);
  double scale = w1 * w2;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    double z11 = _getIVAR(db, iech1, ivar);
    double z12 = _getIVAR(db, iech2, ivar);
    if (FFFF(z11) || FFFF(z12)) continue;
    for (int jvar = 0; jvar <= ivar; jvar++)
    {
      double z21 = _getIVAR(db, iech1, jvar);
      double z22 = _getIVAR(db, iech2, jvar);
      if (FFFF(z21) || FFFF(z22)) continue;
      double value = (z12 - z11) * (z22 - z21) / 2.;
      _setResult(iech1, iech2, nvar, ipas, ivar, jvar, 0, scale, dist, value);
    }
  }
}

void AVario::_evaluateMadogram(
  Db* db, int nvar, int iech1, int iech2, int ipas, double dist, bool /*do_asym*/)
{
  double w1 = db->getWeight(iech1);
  double w2 = db->getWeight(iech2);
  if (FFFF(w1) || FFFF(w2)) return;
  dist         = ABS(dist);
  double scale = w1 * w2;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    double z11 = _getIVAR(db, iech1, ivar);
    double z12 = _getIVAR(db, iech2, ivar);
    if (FFFF(z11) || FFFF(z12)) continue;
    for (int jvar = 0; jvar <= ivar; jvar++)
    {
      double z21 = _getIVAR(db, iech1, jvar);
      double z22 = _getIVAR(db, iech2, jvar);
      if (FFFF(z21) || FFFF(z22)) continue;
      double value = sqrt(ABS((z12 - z11) * (z22 - z21))) / 2.;
      _setResult(iech1, iech2, nvar, ipas, ivar, jvar, 0, scale, dist, value);
    }
  }
}

void AVario::_evaluateRodogram(
  Db* db, int nvar, int iech1, int iech2, int ipas, double dist, bool /*do_asym*/)
{
  double w1 = db->getWeight(iech1);
  double w2 = db->getWeight(iech2);
  if (FFFF(w1) || FFFF(w2)) return;
  dist         = ABS(dist);
  double scale = w1 * w2;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    double z11 = _getIVAR(db, iech1, ivar);
    double z12 = _getIVAR(db, iech2, ivar);
    if (FFFF(z11) || FFFF(z12)) continue;
    for (int jvar = 0; jvar <= ivar; jvar++)
    {
      double z21 = _getIVAR(db, iech1, jvar);
      double z22 = _getIVAR(db, iech2, jvar);
      if (FFFF(z21) || FFFF(z22)) continue;
      double value = pow(ABS((z12 - z11) * (z22 - z21)), 0.25) / 2.;
      _setResult(iech1, iech2, nvar, ipas, ivar, jvar, 0, scale, dist, value);
    }
  }
}

void AVario::_evaluatePoisson(
  Db* db, int nvar, int iech1, int iech2, int ipas, double dist, bool /*do_asym*/)
{
  double w1 = db->getWeight(iech1);
  double w2 = db->getWeight(iech2);
  if (FFFF(w1) || FFFF(w2)) return;
  dist         = ABS(dist);
  double scale = (w1 * w2) / (w1 + w2);
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    double z11 = _getIVAR(db, iech1, ivar);
    double z12 = _getIVAR(db, iech2, ivar);
    if (FFFF(z11) || FFFF(z12)) continue;
    for (int jvar = 0; jvar <= ivar; jvar++)
    {
      double z21 = _getIVAR(db, iech1, jvar);
      double z22 = _getIVAR(db, iech2, jvar);
      if (FFFF(z21) || FFFF(z22)) continue;
      double value = (z12 - z11) * (z22 - z21) / 2.;
      _setResult(iech1, iech2, nvar, ipas, ivar, jvar, 0, scale, dist, value);
    }
  }
}

void AVario::_evaluateCovariance(
  Db* db, int nvar, int iech1, int iech2, int ipas, double dist, bool do_asym)
{
  double w1 = db->getWeight(iech1);
  double w2 = db->getWeight(iech2);
  if (FFFF(w1) || FFFF(w2)) return;
  int orient = (dist > 0) ? 1 : -1;
  dist       = ABS(dist);
  double scale = w1 * w2;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    double z11 = _getIVAR(db, iech1, ivar);
    double z12 = _getIVAR(db, iech2, ivar);
    if (FFFF(z11) || FFFF(z12)) continue;
    for (int jvar = 0; jvar <= ivar; jvar++)
    {
      double z21 = _getIVAR(db, iech1, jvar);
      double z22 = _getIVAR(db, iech2, jvar);
      if (!FFFF(z22))
      {
        double value = z11 * z22;
        _setResult(iech1, iech2, nvar, ipas, ivar, jvar, orient, scale, dist, value);
      }
      if (!FFFF(z21) && do_asym)
      {
        double value = z12 * z21;
        _setResult(iech1, iech2, nvar, ipas, ivar, jvar, -orient, scale, dist, value);
      }
    }
  }
}

void AVario::_evaluateCovariogram(
  Db* db, int nvar, int iech1, int iech2, int ipas, double dist, bool do_asym)
{
  double w1 = db->getWeight(iech1);
  double w2 = db->getWeight(iech2);
  if (FFFF(w1) || FFFF(w2)) return;
  int orient = (dist > 0) ? 1 : -1;
  dist       = ABS(dist);
  double scale = w2;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    double z11 = _getIVAR(db, iech1, ivar);
    double z12 = _getIVAR(db, iech2, ivar);
    if (FFFF(z11) || FFFF(z12)) continue;
    for (int jvar = 0; jvar <= ivar; jvar++)
    {
      double z21 = _getIVAR(db, iech1, jvar);
      double z22 = _getIVAR(db, iech2, jvar);
      if (!FFFF(z22))
      {
        double value = z11 * z22;
        _setResult(iech1, iech2, nvar, ipas, ivar, jvar, orient, scale, dist, value);
      }
      if (!FFFF(z21) && do_asym)
      {
        double value = z12 * z21;
        _setResult(iech1, iech2, nvar, ipas, ivar, jvar, -orient, scale, dist, value);
      }
    }
  }
}

void AVario::_evaluateOrder4(
  Db* db, int nvar, int iech1, int iech2, int ipas, double dist, bool /*do_asym*/)
{
  double w1 = db->getWeight(iech1);
  double w2 = db->getWeight(iech2);
  if (FFFF(w1) || FFFF(w2)) return;
  dist         = ABS(dist);
  double scale = w1 * w2;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    double z11 = _getIVAR(db, iech1, ivar);
    double z12 = _getIVAR(db, iech2, ivar);
    if (FFFF(z11) || FFFF(z12)) continue;
    for (int jvar = 0; jvar <= ivar; jvar++)
    {
      double z21 = _getIVAR(db, iech1, jvar);
      double z22 = _getIVAR(db, iech2, jvar);
      if (FFFF(z21) || FFFF(z22)) continue;
      double value = (z12 - z11) * (z22 - z21);
      value        = value * value / 2.;
      _setResult(iech1, iech2, nvar, ipas, ivar, jvar, 0, scale, dist, value);
    }
  }
}

String AVario::_elemString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);
  std::stringstream sstr;

  // Print the calculation type

  switch (getCalcul().toEnum())
  {
    case ECalcVario::E_UNDEFINED: sstr << toTitle(0, "Undefined"); break;

    case ECalcVario::E_VARIOGRAM:
      sstr << toTitle(0, "Variogram characteristics");
      break;

    case ECalcVario::E_MADOGRAM:
      sstr << toTitle(0, "Madogram characteristics");
      break;

    case ECalcVario::E_RODOGRAM:
      sstr << toTitle(0, "Rodogram characteristics");
      break;

    case ECalcVario::E_POISSON:
      sstr << toTitle(0, "Poisson variogram characteristics");
      break;

    case ECalcVario::E_COVARIANCE:
      sstr << toTitle(0, "Covariance characteristics");
      break;

    case ECalcVario::E_COVARIANCE_NC:
      sstr << toTitle(0, "Non-centered Covariance characteristics");
      break;

    case ECalcVario::E_COVARIOGRAM:
      sstr << toTitle(0, "Transitive Covariogram characteristics");
      break;

    case ECalcVario::E_GENERAL1:
      sstr << toTitle(0, "Generalized Variogram of order 1 characteristics");
      break;

    case ECalcVario::E_GENERAL2:
      sstr << toTitle(0, "Generalized Variogram of order 2 characteristics");
      break;

    case ECalcVario::E_GENERAL3:
      sstr << toTitle(0, "Generalized Variogram of order 3 characteristics");
      break;

    case ECalcVario::E_ORDER4: sstr << toTitle(0, "Order-4 Variogram"); break;

    case ECalcVario::E_TRANS1:
      sstr << toTitle(0, "Cross-to_simple Variogram ratio G12/G1");
      break;

    case ECalcVario::E_TRANS2:
      sstr << toTitle(0, "Cross-to_simple Variogram ratio G12/G2");
      break;

    case ECalcVario::E_BINORMAL:
      sstr << toTitle(0, "Cross-to_simple Variogram ratio G12/sqrt(G1*G2)");
      break;

    default: break;
  }
  return sstr.str();
}

/**
 * Convert the Calculation Name into a Calculation Type (ECalcVario)
 *
 * @return The corresponding ECalcVario enum
 */
ECalcVario AVario::getCalculType(const String& calcul_name)
{
  ECalcVario calcul_type;

  if (calcul_name == "undefined")
    calcul_type = ECalcVario::UNDEFINED;
  else if (calcul_name == "vg")
    calcul_type = ECalcVario::VARIOGRAM;
  else if (calcul_name == "cov")
    calcul_type = ECalcVario::COVARIANCE;
  else if (calcul_name == "covnc")
    calcul_type = ECalcVario::COVARIANCE_NC;
  else if (calcul_name == "covg")
    calcul_type = ECalcVario::COVARIOGRAM;
  else if (calcul_name == "mado")
    calcul_type = ECalcVario::MADOGRAM;
  else if (calcul_name == "rodo")
    calcul_type = ECalcVario::RODOGRAM;
  else if (calcul_name == "poisson")
    calcul_type = ECalcVario::POISSON;
  else if (calcul_name == "general1")
    calcul_type = ECalcVario::GENERAL1;
  else if (calcul_name == "general2")
    calcul_type = ECalcVario::GENERAL2;
  else if (calcul_name == "general3")
    calcul_type = ECalcVario::GENERAL3;
  else if (calcul_name == "order4")
    calcul_type = ECalcVario::ORDER4;
  else if (calcul_name == "trans1")
    calcul_type = ECalcVario::TRANS1;
  else if (calcul_name == "trans2")
    calcul_type = ECalcVario::TRANS2;
  else if (calcul_name == "binormal")
    calcul_type = ECalcVario::BINORMAL;
  else
  {
    messerr("Invalid variogram calculation name : %s", calcul_name.c_str());
    messerr("The only valid names are:");
    messerr("vg       : Variogram");
    messerr("cov      : Covariance");
    messerr("covnc    : Non-centered ergodic covariance");
    messerr("covg     : Covariogram");
    messerr("mado     : Madogram");
    messerr("rodo     : Rodogram");
    messerr("poisson  : Poisson");
    messerr("general1 : Generalized variogram of order 1");
    messerr("general2 : Generalized variogram of order 2");
    messerr("general3 : Generalized variogram of order 3");
    messerr("order4   : Variogram of order 4");
    messerr("trans1   : Cross-to-Simple Variogram G12/G1");
    messerr("trans2   : Cross-to-Simple Variogram G12/G1");
    messerr("binormal : Cross-to-Simple Variogram G12/sqrt(G1*G2)");

    calcul_type = ECalcVario::UNDEFINED;
  }
  return calcul_type;
}

void AVario::setCalculByName(const String& calcul_name)
{
  const ECalcVario& calcul = getCalculType(calcul_name);
  setCalcul(calcul);
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

