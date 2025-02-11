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

#include "Enum/EKrigOpt.hpp"
#include "Covariances/CovCalcMode.hpp"

class Db;
class DbGrid;
class MatrixRectangular;

class GSTLEARN_EXPORT KrigOpt
{
public:
  KrigOpt(const EKrigOpt& calcul = EKrigOpt::POINT);
  KrigOpt(const KrigOpt &m);
  KrigOpt& operator=(const KrigOpt &m);
  virtual ~KrigOpt();

  int setMatLC(const MatrixRectangular* matLC, int nvar);
  int setKrigingOption(const EKrigOpt& calcul  = EKrigOpt::POINT,
                       DbGrid* dbgrid          = nullptr,
                       const VectorInt& ndiscs = VectorInt(),
                       bool flag_per_cell      = false);
  int setKrigingDGM(bool flag_dgm);
  void setMode(const CovCalcMode* mode);

  const CovCalcMode& getMode() const { return _mode; }
  const EKrigOpt& getCalcul() const { return _calcul; }
  int getNDisc() const { return _ndiscNumber; }
  bool isFlagCell() const { return _flagPerCell; }
  void blockDiscretize(int iechout, bool flagRandom = false, int seed = 1234546) const;
  VectorDouble getDisc1VD(int idisc) const;
  VectorDouble getDisc2VD(int idisc) const;
  VectorVectorDouble getDisc1VVD() const;
  VectorVectorDouble getDisc2VVD() const;
  const MatrixRectangular* getMatLC() const { return _matLC; }
  double getMatCLValue(int ivarcl, int ivar) const;
  bool isMatLC() const { return _matLC != nullptr; }
  int getNvarCL() const;

private:
  double _getDisc1(int idisc, int idim) const;
  double _getDisc2(int idisc, int idim) const;

private:
  EKrigOpt _calcul;
  CovCalcMode _mode;

  bool _flagPerCell;
  int _ndim;
  int _ndiscNumber;
  VectorInt _ndiscs;
  mutable VectorVectorDouble _disc1; // Dimension: ndiscNumber, ndim
  mutable VectorVectorDouble _disc2; // Dimension: ndiscNumber, ndim

  bool _flagDGM;

  const MatrixRectangular* _matLC; // Pointer not to be deleted
  DbGrid* _dbgrid; // Pointer to the DbGrid (not to be deleted)
};
