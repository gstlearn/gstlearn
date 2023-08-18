/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Simulation/BooleanObject.hpp"
#include "Basic/AStringable.hpp"

class AShape;

class GSTLEARN_EXPORT ModelBoolean: public AStringable
{
public:
  ModelBoolean(double thetaCst = 1., bool flagStat = true);
  ModelBoolean(const ModelBoolean &r);
  ModelBoolean& operator=(const ModelBoolean &r);
  virtual ~ModelBoolean();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int getNbTokens() const { return (int) _shapes.size(); }
  void addToken(const AShape& token);
  void normalizeProportions();
  BooleanObject* generateObject(int ndim) const;
  const AShape* getToken(int itok) const { return _shapes[itok]; }

  bool   isFlagStat() const { return _flagStat; }
  double getThetaCst() const { return _thetaCst; }

  void setFlagStat(bool flagStat) { _flagStat = flagStat; }
  void setThetaCst(double thetaCst) { _thetaCst = thetaCst; }

private:
  bool   _flagStat;
  double _thetaCst;
  std::vector<AShape*> _shapes; // List of the Token
};
