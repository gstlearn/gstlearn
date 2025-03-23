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
class ANeigh;
class ModelGeneric;

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
  int setRankColCok(const VectorInt& rank_colcok);
  void setMode(const CovCalcMode* mode);

  const CovCalcMode& getMode() const { return _mode; }
  const EKrigOpt& getCalcul() const { return _calcul; }

  int getNDisc() const { return _nDiscNumber; }
  const VectorInt& getDiscs() const { return _ndiscs; }
  int getDisc(int idim) const { return _ndiscs[idim]; }
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
  const VectorInt& getRankColcok() const { return _rankColcok; }
  int getRankColcok(int i) const { return _rankColcok[i]; }
  bool hasColcok() const { return _flagColcok; }

  bool isValid(const Db* dbout, const ANeigh* neigh, const ModelGeneric* model) const;
  void dumpOptions() const;

private:
  double _getDisc1(int idisc, int idim) const;
  double _getDisc2(int idisc, int idim) const;

private:
  // General information
  EKrigOpt _calcul;
  CovCalcMode _mode;

  // Discretization of cell for block calculation
  mutable bool _flagPerCell;
  int _nDiscDim;
  int _nDiscNumber;
  VectorInt _ndiscs;
  mutable VectorVectorDouble _disc1; // Dimension: ndiscNumber, ndim
  mutable VectorVectorDouble _disc2; // Dimension: ndiscNumber, ndim

  // Discrete Gaussian Model
  bool _flagDGM;

  // Colocated Kriging option
  bool _flagColcok;
  VectorInt _rankColcok;
  mutable VectorDouble _valuesColcok;

  // Matrix used for variable combination
  const MatrixRectangular* _matLC; // Pointer not to be deleted

  DbGrid* _dbgrid; // Pointer to the DbGrid (not to be deleted)
};
