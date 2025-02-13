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
  virtual double _getIVAR(const Db* db, int iech, int ivar) const = 0;
  virtual void _setResult(int iech1,
                          int iech2,
                          int nvar,
                          int ilag,
                          int ivar,
                          int jvar,
                          int orient,
                          double ww,
                          double dist,
                          double value) = 0;

  String _elemString(const AStringFormat* strfmt) const;

  void _evaluateVariogram(Db* db,
                          int nvar,
                          int iech1,
                          int iech2,
                          int ilag,
                          double dist,
                          bool do_asym);
  void _evaluateMadogram(Db* db,
                         int nvar,
                         int iech1,
                         int iech2,
                         int ilag,
                         double dist,
                         bool do_asym);
  void _evaluateRodogram(Db* db,
                         int nvar,
                         int iech1,
                         int iech2,
                         int ilag,
                         double dist,
                         bool do_asym);
  void _evaluatePoisson(Db* db,
                        int nvar,
                        int iech1,
                        int iech2,
                        int ilag,
                        double dist,
                        bool do_asym);
  void _evaluateCovariance(Db* db,
                           int nvar,
                           int iech1,
                           int iech2,
                           int ilag,
                           double dist,
                           bool do_asym);
  void _evaluateCovariogram(Db* db,
                            int nvar,
                            int iech1,
                            int iech2,
                            int ilag,
                            double dist,
                            bool do_asym);
  void _evaluateOrder4(Db* db,
                       int nvar,
                       int iech1,
                       int iech2,
                       int ilag,
                       double dist,
                       bool do_asym);

  void (AVario::*_evaluate)(Db* db,
                            int nvar,
                            int iech1,
                            int iech2,
                            int ilag
                              ,
                    double dist,
                    bool do_asym);

private:
  ECalcVario _calcul;
};
