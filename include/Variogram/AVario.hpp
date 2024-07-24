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

  void evaluate(Db *db,
                int nvar,
                int iech1,
                int iech2,
                int ipas,
                double dist,
                int do_asym = 0);
  const ECalcVario& getCalcul() const { return _calcul; }
  void setCalcul(const ECalcVario &calcul) { _calcul = calcul; }

protected:
  virtual double _getIVAR(const Db *db, int iech, int ivar) const = 0;
  virtual void _setResult(int iech1,
                          int iech2,
                          int nvar,
                          int ipas,
                          int ivar,
                          int jvar,
                          int orient,
                          double ww,
                          double dist,
                          double value) = 0;

  String _elemString(const AStringFormat* strfmt) const;

protected:
  ECalcVario _calcul;
};
