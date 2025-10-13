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
#include "Matrix/MatrixDense.hpp"

namespace gstlrn
{
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

  Id setOptionCalcul(const EKrigOpt& calcul  = EKrigOpt::POINT,
                      const VectorInt& ndiscs = VectorInt(),
                      bool flag_per_cell      = false);
  Id setColCok(const VectorInt& rank_colcok);
  Id setMatLC(const MatrixDense* matLC);
  void setMode(const CovCalcMode* mode);
  void setOptionDGM(bool flag_dgm);

  const CovCalcMode& getMode() const { return _mode; }
  const EKrigOpt& getCalcul() const { return _calcul; }

  Id getNDisc() const { return _nDiscNumber; }
  bool hasDiscs() const { return ! _ndiscs.empty(); }
  const VectorInt& getDiscs() const { return _ndiscs; }
  Id getDisc(Id idim) const { return _ndiscs[idim]; }
  bool hasFlagPerCell() const { return _flagPerCell; }
  void blockDiscretize(Id iechout, bool flagRandom = false, Id seed = 1234546) const;
  VectorDouble getDisc1VD(Id idisc) const;
  VectorDouble getDisc2VD(Id idisc) const;
  VectorVectorDouble getDisc1VVD() const;
  VectorVectorDouble getDisc2VVD() const;
  const MatrixDense* getMatLC() const { return _matLC; }
  double getMatLCValue(Id ivarcl, Id ivar) const;
  bool hasMatLC() const { return _matLC != nullptr; }
  Id getMatLCNRows() const { return _matLC->getNRows(); }
  Id getNvarLC() const;
  const VectorInt& getRankColcok() const { return _rankColcok; }
  Id getRankColcok(Id i) const { return _rankColcok[i]; }
  bool hasColcok() const { return _flagColcok; }
  bool hasFlagDGM() const { return _flagDGM; }

  bool isCorrect(const Db* dbout, const ANeigh* neigh, const ModelGeneric* model) const;
  void dumpOptions() const;

private:
  double _getDisc1(Id idisc, Id idim) const;
  double _getDisc2(Id idisc, Id idim) const;
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
  Id _nDiscDim;
  Id _nDiscNumber;
  VectorInt _ndiscs;
  mutable VectorVectorDouble _disc1; // Dimension: ndiscNumber, ndim
  mutable VectorVectorDouble _disc2; // Dimension: ndiscNumber, ndim

  // Discrete Gaussian Model
  bool _flagDGM;

  // Colocated Kriging option
  bool _flagColcok;
  VectorInt _rankColcok;

  // Matrix used for variable combination
  const MatrixDense* _matLC; // Pointer not to be deleted

  mutable const DbGrid* _dbgrid; // Pointer to the DbGrid (not to be deleted)
};
}