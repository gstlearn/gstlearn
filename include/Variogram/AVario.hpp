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

  const ECalcVario& getCalcul() const { return _calcul; }
  void setCalcul(const ECalcVario& calcul);
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
                                double z11,
                                double z21,
                                double z12,
                                double z22,
                                double dist,
                                double value) = 0;

  String _elemString(const AStringFormat* strfmt) const;

  void _evaluateVariogram(Db* db,
                          Id nvar,
                          Id iech1,
                          Id iech2,
                          Id idir,
                          Id ilag,
                          double dist,
                          bool do_asym);
  void _evaluateMadogram(Db* db,
                         Id nvar,
                         Id iech1,
                         Id iech2,
                         Id idir,
                         Id ilag,
                         double dist,
                         bool do_asym);
  void _evaluateRodogram(Db* db,
                         Id nvar,
                         Id iech1,
                         Id iech2,
                         Id idir,
                         Id ilag,
                         double dist,
                         bool do_asym);
  void _evaluatePoisson(Db* db,
                        Id nvar,
                        Id iech1,
                        Id iech2,
                        Id idir,
                        Id ilag,
                        double dist,
                        bool do_asym);
  void _evaluateCovariance(Db* db,
                           Id nvar,
                           Id iech1,
                           Id iech2,
                           Id idir,
                           Id ilag,
                           double dist,
                           bool do_asym);
  void _evaluateCovariogram(Db* db,
                            Id nvar,
                            Id iech1,
                            Id iech2,
                            Id idir,
                            Id ilag,
                            double dist,
                            bool do_asym);
  void _evaluateOrder4(Db* db,
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

  void _storage(Id iech1, Id iech2, double dist, double value);
  bool _isStorage() const { return _flagStorage; }
  const VectorDouble& _getStorage(Id ipair) const { return _tabStorage[ipair]; }
  Id _getStorageSize() const { return static_cast<Id>(_tabStorage.size()); }

private:
  ECalcVario _calcul;
  mutable bool _flagNeedStats;
  mutable bool _flagErgodic;
  mutable bool _flagAsym;
  mutable bool _flagCentered;
  mutable bool _flagScaled;
  bool _flagStorage;              // Store information on each pair
  VectorVectorDouble _tabStorage; // Storage of the information on each pair; Dimension [npairs][4]
};
} // namespace gstlrn