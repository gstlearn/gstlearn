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
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/ECalcVario.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"

namespace gstlrn
{
class Db;
class ECalcVario;

class GSTLEARN_EXPORT AVario: public AStringable, public ICloneable
{
public:
  AVario();
  AVario(const AVario& r);
  AVario& operator=(const AVario& r);
  virtual ~AVario();

  const ECalcVario& getCalcul() const { return _calculType; }
  void setCalcul(const ECalcVario& calculType);
  void setStorage(bool flag);
  void setErgodic(bool flag) { _flagErgodic = flag; }
  bool getNeedStats() const { return _flagNeedStats; }
  bool getErgodic() const { return _flagErgodic; }
  bool getFlagAsym() const { return _flagAsym; }
  bool getCentered() const { return _flagCentered; }
  bool getScaled() const { return _flagScaled; }

protected:
  virtual double _getIVAR(const Db* db, Id iech, Id ivar) const = 0;

  virtual void _setAVarioResult(Id iech1,
                                Id iech2,
                                Id nvar,
                                Id idir,
                                Id ilag,
                                Id ivar,
                                Id jvar,
                                Id orient,
                                double ww,
                                double w1,
                                double w2,
                                double z1,
                                double z2,
                                double dist,
                                double value) = 0;

  String _elemString(const AStringFormat* strfmt) const;

  double _evaluateScaleVariogram(double w1, double w2);
  double _evaluateScalePoisson(double w1, double w2);
  double _evaluateScaleCovariogram(double w1, double w2);

  double _evaluateResultVariogram(Id icase, double z11, double z21, double z12, double z22);
  double _evaluateResultMadogram(Id icase, double z11, double z21, double z12, double z22);
  double _evaluateResultRodogram(Id icase, double z11, double z21, double z12, double z22);
  double _evaluateResultCovariance(Id icase, double z11, double z21, double z12, double z22);
  double _evaluateResultOrder4(Id icase, double z11, double z21, double z12, double z22);

  void _evaluateSymmetric(Db* db,
                          Id nvar,
                          Id iech1,
                          Id iech2,
                          Id idir,
                          Id ilag,
                          double dist,
                          bool do_asym);
  void _evaluateASymmetric(Db* db,
                           Id nvar,
                           Id iech1,
                           Id iech2,
                           Id idir,
                           Id ilag,
                           double dist,
                           bool do_asym);

  void (AVario::*_evaluate)(Db* db,
                            Id nvar,
                            Id iech1,
                            Id iech2,
                            Id idir,
                            Id ilag,
                            double dist,
                            bool do_asym);
  double (AVario::*_evaluateScale)(double w1, double w2);
  double (AVario::*_evaluateResult)(Id icase, double z11, double z21, double z12, double z22);

  void _storage(Id iech1, Id iech2, double dist, double value);
  bool _isStorage() const { return _flagStorage; }
  const VectorDouble& _getStorage(Id ipair) const { return _tabStorage[ipair]; }
  Id _getStorageSize() const { return static_cast<Id>(_tabStorage.size()); }

private:
  ECalcVario _calculType;
  mutable bool _flagNeedStats;
  mutable bool _flagErgodic;
  mutable bool _flagAsym;
  mutable bool _flagCentered;
  mutable bool _flagScaled;
  bool _flagStorage;              // Store information on each pair
  VectorVectorDouble _tabStorage; // Storage of the information on each pair; Dimension [npairs][4]
};

GSTLEARN_EXPORT bool getFlagAVarioCheck();
GSTLEARN_EXPORT void setFlagAVarioCheck(bool flag);

} // namespace gstlrn