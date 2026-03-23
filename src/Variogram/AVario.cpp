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
  static bool _flagAVarioCheck = false;

  AVario::AVario()
    : AStringable()
    , _calculType(ECalcVario::VARIOGRAM)
    , _flagNeedStats(false)
    , _flagErgodic(true)
    , _flagAsym(false)
    , _flagCentered(false)
    , _flagScaled(false)
    , _flagStorage(false)
    , _tabStorage()
  {
  }

  AVario::AVario(const AVario& r)
    : AStringable(r)
    , _calculType(r._calculType)
    , _flagNeedStats(r._flagNeedStats)
    , _flagErgodic(r._flagErgodic)
    , _flagAsym(r._flagAsym)
    , _flagCentered(r._flagCentered)
    , _flagScaled(r._flagScaled)
    , _flagStorage(r._flagStorage)
    , _tabStorage(r._tabStorage)
  {
  }

  AVario& AVario::operator=(const AVario& r)
  {
    if (this != &r)
    {
      AStringable::operator=(r);
      _calculType = r._calculType;
      _flagNeedStats = r._flagNeedStats;
      _flagErgodic = r._flagErgodic;
      _flagAsym = r._flagAsym;
      _flagCentered = r._flagCentered;
      _flagScaled = r._flagScaled;
      _flagStorage = r._flagStorage;
      _tabStorage = r._tabStorage;
    }
    return *this;
  }

  AVario::~AVario() {}

  void setFlagAVarioCheck(bool flag)
  {
    _flagAVarioCheck = flag;
  }

  bool getFlagAVarioCheck()
  {
    return _flagAVarioCheck;
  }

  double AVario::_evaluateScaleVariogram(double w1, double w2)
  {
    return w1 * w2;
  }

  double AVario::_evaluateScalePoisson(double w1, double w2)
  {
    return (w1 * w2) / (w1 + w2);
  }

  double AVario::_evaluateScaleCovariogram(double w1, double w2)
  {
    DECLARE_UNUSED(w1);
    return w2;
  }

  double AVario::_evaluateResultVariogram(Id icase,
                                          double z11,
                                          double z21,
                                          double z12,
                                          double z22)
  {
    DECLARE_UNUSED(icase);
    return (z12 - z11) * (z22 - z21) / 2.;
  }

  double AVario::_evaluateResultMadogram(Id icase,
                                         double z11,
                                         double z21,
                                         double z12,
                                         double z22)
  {
    DECLARE_UNUSED(icase);
    return sqrt(ABS((z12 - z11) * (z22 - z21))) / 2.;
  }

  double AVario::_evaluateResultRodogram(Id icase,
                                         double z11,
                                         double z21,
                                         double z12,
                                         double z22)
  {
    DECLARE_UNUSED(icase);
    return pow(ABS((z12 - z11) * (z22 - z21)), 0.25) / 2.;
  }

  double AVario::_evaluateResultCovariance(Id icase,
                                           double z11,
                                           double z21,
                                           double z12,
                                           double z22)
  {
    if (icase == 1) return z11 * z22;
    return z12 * z21;
  }

  double AVario::_evaluateResultOrder4(Id icase,
                                       double z11,
                                       double z21,
                                       double z12,
                                       double z22)
  {
    DECLARE_UNUSED(icase);
    double value = (z12 - z11) * (z22 - z21);
    return value * value / 2.;
  }

  void AVario::_evaluateSymmetric(Db* db,
                                  Id nvar,
                                  Id iech1,
                                  Id iech2,
                                  Id idir,
                                  Id ilag,
                                  double dist,
                                  bool do_asym)
  {
    DECLARE_UNUSED(do_asym);
    double w1 = db->getWeight(iech1);
    double w2 = db->getWeight(iech2);
    if (FFFF(w1) || FFFF(w2)) return;
    dist = ABS(dist);
    double scale = (this->*_evaluateScale)(w1, w2);
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
        double value = (this->*_evaluateResult)(0, z11, z21, z12, z22);
        _setAVarioResult(iech1,
                         iech2,
                         nvar,
                         idir,
                         ilag,
                         ivar,
                         jvar,
                         0,
                         scale,
                         w1,
                         w2,
                         TEST,
                         TEST,
                         dist,
                         value);
      }
    }
  }

  void AVario::_evaluateASymmetric(Db* db,
                                   Id nvar,
                                   Id iech1,
                                   Id iech2,
                                   Id idir,
                                   Id ilag,
                                   double dist,
                                   bool do_asym)
  {
    double w1 = db->getWeight(iech1);
    double w2 = db->getWeight(iech2);
    if (FFFF(w1) || FFFF(w2)) return;
    Id orient = (dist > 0) ? 1 : -1;
    dist = ABS(dist);
    double scale = (this->*_evaluateScale)(w1, w2);
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
          double value = (this->*_evaluateResult)(1, z11, z21, z12, z22);
          _setAVarioResult(iech1,
                           iech2,
                           nvar,
                           idir,
                           ilag,
                           ivar,
                           jvar,
                           orient,
                           scale,
                           w1,
                           w2,
                           z11,
                           z22,
                           dist,
                           value);
        }
        if (!FFFF(z21) && do_asym)
        {
          double value = (this->*_evaluateResult)(-1, z11, z21, z12, z22);
          _setAVarioResult(iech1,
                           iech2,
                           nvar,
                           idir,
                           ilag,
                           ivar,
                           jvar,
                           -orient,
                           scale,
                           w1,
                           w2,
                           z21,
                           z12,
                           dist,
                           value);
        }
      }
    }
  }

  String AVario::_elemString(const AStringFormat* strfmt) const
  {
    DECLARE_UNUSED(strfmt);
    std::stringstream sstr;
    String total_string =
      std::string(getCalcul().getDescr()) + " characteristics";
    sstr << toStrTitle(0, total_string.c_str());
    return sstr.str();
  }

  void AVario::setCalcul(const ECalcVario& calculType)
  {
    // Set the calculation type
    _calculType = calculType;

    // Define the internal evaluation function
    switch (_calculType.toEnum())
    {
      case ECalcVario::E_VARIOGRAM:
      case ECalcVario::E_TRANS1:
      case ECalcVario::E_TRANS2:
      case ECalcVario::E_BINORMAL:
      {
        _evaluate = &AVario::_evaluateSymmetric;
        _evaluateScale = &AVario::_evaluateScaleVariogram;
        _evaluateResult = &AVario::_evaluateResultVariogram;
        break;
      }

      case ECalcVario::E_MADOGRAM:
      {
        _evaluate = &AVario::_evaluateSymmetric;
        _evaluateScale = &AVario::_evaluateScaleVariogram;
        _evaluateResult = &AVario::_evaluateResultMadogram;
        break;
      }

      case ECalcVario::E_RODOGRAM:
      {
        _evaluate = &AVario::_evaluateSymmetric;
        _evaluateScale = &AVario::_evaluateScaleVariogram;
        _evaluateResult = &AVario::_evaluateResultRodogram;
        break;
      }

      case ECalcVario::E_POISSON:
      {
        _evaluate = &AVario::_evaluateSymmetric;
        _evaluateScale = &AVario::_evaluateScalePoisson;
        _evaluateResult = &AVario::_evaluateResultVariogram;
        break;
      }

      case ECalcVario::E_COVARIANCE:
      case ECalcVario::E_COVARIANCE_NC:
      case ECalcVario::E_CORRELOGRAM:
      {
        _evaluate = &AVario::_evaluateASymmetric;
        _evaluateScale = &AVario::_evaluateScaleVariogram;
        _evaluateResult = &AVario::_evaluateResultCovariance;
        break;
      }

      case ECalcVario::E_COVARIOGRAM:
      {
        _evaluate = &AVario::_evaluateASymmetric;
        _evaluateScale = &AVario::_evaluateScaleCovariogram;
        _evaluateResult = &AVario::_evaluateResultCovariance;
        break;
      }

      case ECalcVario::E_ORDER4:
      {
        _evaluate = &AVario::_evaluateSymmetric;
        _evaluateScale = &AVario::_evaluateScaleVariogram;
        _evaluateResult = &AVario::_evaluateResultOrder4;
        break;
      }

      default:
      {
        _evaluate = nullptr;
        messageAbort("AVario::evaluate() ignores current calculation type");
        break;
      }
    }

    // Define the different flags
    String cle = std::string(_calculType.getKey());
    _flagAsym = !ECalcVarioAttr.at(cle).isSymmetric;
    _flagCentered = ECalcVarioAttr.at(cle).isCentered;
    _flagScaled = ECalcVarioAttr.at(cle).isScaled;
    _flagNeedStats = _flagCentered || _flagScaled;
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
