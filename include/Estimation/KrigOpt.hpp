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
#include "Matrix/MatrixRectangular.hpp"

class Db;
class DbGrid;
class ANeigh;
class ModelGeneric;

class GSTLEARN_EXPORT KrigOpt
{
public:
  KrigOpt(const EKrigOpt& calcul = EKrigOpt::POINT);
  KrigOpt(const KrigOpt &m);
  KrigOpt& operator=(const KrigOpt &m);
  virtual ~KrigOpt();

  int setOptionCalcul(const EKrigOpt& calcul  = EKrigOpt::POINT,
                      const VectorInt& ndiscs = VectorInt(),
                      bool flag_per_cell      = false);
  int setColCok(const VectorInt& rank_colcok);
  int setMatLC(const MatrixRectangular* matLC);
  void setMode(const CovCalcMode* mode);
  void setOptionDGM(bool flag_dgm);

  const CovCalcMode& getMode() const { return _mode; }
  const EKrigOpt& getCalcul() const { return _calcul; }

  int getNDisc() const { return _nDiscNumber; }
  bool hasDiscs() const { return ! _ndiscs.empty(); }
  const VectorInt& getDiscs() const { return _ndiscs; }
  int getDisc(int idim) const { return _ndiscs[idim]; }
  bool hasFlagPerCell() const { return _flagPerCell; }
  void blockDiscretize(int iechout, bool flagRandom = false, int seed = 1234546) const;
  VectorDouble getDisc1VD(int idisc) const;
  VectorDouble getDisc2VD(int idisc) const;
  VectorVectorDouble getDisc1VVD() const;
  VectorVectorDouble getDisc2VVD() const;
  const MatrixRectangular* getMatLC() const { return _matLC; }
  double getMatLCValue(int ivarcl, int ivar) const;
  bool hasMatLC() const { return _matLC != nullptr; }
  int getMatLCNRows() const { return _matLC->getNRows(); }
  int getNvarLC() const;
  const VectorInt& getRankColcok() const { return _rankColcok; }
  int getRankColcok(int i) const { return _rankColcok[i]; }
  bool hasColcok() const { return _flagColcok; }
  bool hasFlagDGM() const { return _flagDGM; }

  bool isCorrect(const Db* dbout, const ANeigh* neigh, const ModelGeneric* model) const;
  void dumpOptions() const;

private:
  double _getDisc1(int idisc, int idim) const;
  double _getDisc2(int idisc, int idim) const;
  bool _isValidCalcul(const Db* dbout, const ANeigh* neigh) const;
  bool _isValidColcok(const Db* dbout, const ModelGeneric* model) const;
  bool _isValidMatLC(const ModelGeneric* model) const;
  bool _isValidDGM(const Db* dbout, const ModelGeneric* model) const;

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

  // Matrix used for variable combination
  const MatrixRectangular* _matLC; // Pointer not to be deleted

  mutable const DbGrid* _dbgrid; // Pointer to the DbGrid (not to be deleted)
};
