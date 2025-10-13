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

class GSTLEARN_EXPORT AVario:  public AStringable, public ICloneable
{
public:
  AVario();
  AVario(const AVario& r);
  AVario& operator=(const AVario& r);
  virtual ~AVario();

  static ECalcVario getCalculType(const String& calcul_name);
  const ECalcVario& getCalcul() const { return _calcul; }
  void setCalcul(const ECalcVario& calcul);
  void setCalculByName(const String& calcul_name);

protected:
  virtual double _getIVAR(const Db* db, Id iech, Id ivar) const = 0;
  virtual void _setResult(Id iech1,
                          Id iech2,
                          Id nvar,
                          Id ilag,
                          Id ivar,
                          Id jvar,
                          Id orient,
                          double ww,
                          double dist,
                          double value) = 0;

  String _elemString(const AStringFormat* strfmt) const;

  void _evaluateVariogram(Db* db,
                          Id nvar,
                          Id iech1,
                          Id iech2,
                          Id ilag,
                          double dist,
                          bool do_asym);
  void _evaluateMadogram(Db* db,
                         Id nvar,
                         Id iech1,
                         Id iech2,
                         Id ilag,
                         double dist,
                         bool do_asym);
  void _evaluateRodogram(Db* db,
                         Id nvar,
                         Id iech1,
                         Id iech2,
                         Id ilag,
                         double dist,
                         bool do_asym);
  void _evaluatePoisson(Db* db,
                        Id nvar,
                        Id iech1,
                        Id iech2,
                        Id ilag,
                        double dist,
                        bool do_asym);
  void _evaluateCovariance(Db* db,
                           Id nvar,
                           Id iech1,
                           Id iech2,
                           Id ilag,
                           double dist,
                           bool do_asym);
  void _evaluateCovariogram(Db* db,
                            Id nvar,
                            Id iech1,
                            Id iech2,
                            Id ilag,
                            double dist,
                            bool do_asym);
  void _evaluateOrder4(Db* db,
                       Id nvar,
                       Id iech1,
                       Id iech2,
                       Id ilag,
                       double dist,
                       bool do_asym);

  void (AVario::*_evaluate)(Db* db,
                            Id nvar,
                            Id iech1,
                            Id iech2,
                            Id ilag
                              ,
                    double dist,
                    bool do_asym);

private:
  ECalcVario _calcul;
};
}