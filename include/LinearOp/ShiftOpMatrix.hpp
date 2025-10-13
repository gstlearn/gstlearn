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

#include "Basic/VectorNumT.hpp"
#include "Basic/VectorT.hpp"
#include "Enum/EPowerPT.hpp"
#include "LinearOp/AShiftOp.hpp"
#include "Mesh/AMesh.hpp"

#include "Matrix/MatrixSparse.hpp"

#include <map>
#include <memory>

#ifndef SWIG

#  include <Eigen/Core>
#  include <Eigen/Dense>
#endif

namespace gstlrn
{
class EConsElem;

class CovAniso;
class AMatrix;
class MatrixSquare;
class MatrixDense;
class MatrixSymmetric;
class MatrixSparse;

class GSTLEARN_EXPORT ShiftOpMatrix: public AShiftOp
{
public:
  ShiftOpMatrix();
  ShiftOpMatrix(const AMesh* amesh, const CovAniso* cova, const Db* dbout = nullptr, bool verbose = false);
  ShiftOpMatrix(const MatrixSparse* S, const VectorDouble& TildeC, const VectorDouble& Lambda, const CovAniso* cova, bool verbose = false);
  ShiftOpMatrix(const ShiftOpMatrix& shift);
  ShiftOpMatrix& operator=(const ShiftOpMatrix& shift);
  virtual ~ShiftOpMatrix();
  /// ICloneable interface
  IMPLEMENT_CLONING(ShiftOpMatrix)
  void normalizeLambdaBySills(const AMesh* mesh) override;
#ifndef SWIG
  Id _addToDest(const constvect inv, vect outv) const override;
#endif

  static ShiftOpMatrix* create(const AMesh* amesh, const CovAniso* cova, const Db* dbout = nullptr, bool verbose = false);
  static ShiftOpMatrix* createFromSparse(const MatrixSparse* S,
                                         const VectorDouble& TildeC,
                                         const VectorDouble& Lambda,
                                         const CovAniso* cova,
                                         bool verbose = false);
  Id initFromMesh(const AMesh* amesh, const CovAniso* cova, const Db* dbout = nullptr, bool flagAdvection = false, bool verbose = false);
  Id initGradFromMesh(const AMesh* amesh, const CovAniso* cova, bool verbose = false, double tol = EPSILON10);
  Id initFromCS(const MatrixSparse* S, const VectorDouble& TildeC, const VectorDouble& Lambda, const CovAniso* cova, bool verbose = false);

  Id getNDim() const
  {
    return _ndim;
  }
  Id getNCovAnisoGradParam() const
  {
    return _nCovAnisoGradParam;
  }
  void prodTildeC(const VectorDouble& x, VectorDouble& y, const EPowerPT& power) const;

  void prodLambdaOnSqrtTildeC(const VectorDouble& inv, VectorDouble& outv, double puis = 2) const;
  MatrixSparse* getS() const { return _S; }
  MatrixSparse* getTildeCGrad(Id iapex, Id igparam) const;
  MatrixSparse* getSGrad(Id iapex, Id igparam) const;

  const VectorDouble& getTildeC() const { return _TildeC; }
  const VectorDouble& getLambdaGrads(Id idim) const { return _LambdaGrad[idim]; }
  double getLambdaGrad(Id idim, Id iapex) const { return _LambdaGrad[idim][iapex]; }
  Id getSGradAddress(Id iapex, Id igparam) const;

  Id getLambdaGradSize() const;
  // void multiplyByValueAndAddDiagonal(double v1 = 1.,double v2 = 0.) override;

private:
  double _getMaxEigenValue() const override;
  Id _buildS(const AMesh* amesh, double tol = EPSILON10);
  Id _buildSGrad(const AMesh* amesh, double tol = EPSILON10);
  void _buildLambda(const AMesh* amesh);
  bool _buildLambdaGrad(const AMesh* amesh);

  static void _loadAux(VectorDouble& tab, const EConsElem& type, Id imesh = 0);
  void _loadHH(const AMesh* amesh, MatrixSymmetric& hh, Id imesh = 0);
  void _loadHHRegular(MatrixSymmetric& hh, Id imesh);
  void _loadHHVariety(MatrixSymmetric& hh, Id imesh);
  void _loadHHGrad(const AMesh* amesh, MatrixSymmetric& hh, Id igparam, Id ipref);
  double _computeGradLogDetHH(const AMesh* amesh, Id igparam, Id ipref, const MatrixSymmetric& HH, MatrixSymmetric& work, MatrixSymmetric& work2);

  void _reset();
  Id _resetGrad();
  void _reallocate(const ShiftOpMatrix& shift);
  static void _projectMesh(const AMesh* amesh,
                           const VectorDouble& srot,
                           Id imesh,
                           double coeff[3][2]);
  Id _preparMatrices(const AMesh* amesh, Id imesh, MatrixSquare& matu, MatrixDense& matw) const;
  Id _prepareMatricesSVariety(const AMesh* amesh,
                               Id imesh,
                               VectorVectorDouble& coords,
                               MatrixDense& matM,
                               MatrixSymmetric& matMtM,
                               AMatrix& matP,
                               double* deter) const;
  Id _prepareMatricesSphere(const AMesh* amesh,
                             Id imesh,
                             VectorVectorDouble& coords,
                             MatrixSquare& matMs,
                             double* deter) const;
  static void _updateCova(std::shared_ptr<CovAniso>& cova, Id imesh);
  VectorT<std::map<Id, double>> _mapCreate() const;
  VectorT<VectorT<std::map<Id, double>>> _mapVectorCreate() const;
  VectorT<std::map<Id, double>> _mapTildeCCreate() const;
  static void _mapTildeCUpdate(std::map<Id, double>& tab,
                               Id ip0,
                               double value,
                               double tol = EPSILON10);

  static void _mapGradUpdate(std::map<std::pair<Id, Id>, double>& tab,
                             Id ip0,
                             Id ip1,
                             double value,
                             double tol = EPSILON10);
  MatrixSparse* _BuildTildeCGradfromMap(std::map<Id, double>& tab) const;
  MatrixSparse* _BuildSGradfromMap(std::map<std::pair<Id, Id>, double>&
                                     tab) const;

  static bool _cond(Id indref, Id igparam, Id ipref);
  void _determineFlagNoStatByHH();
  void _updateHH(MatrixSymmetric& hh, Id imesh);
  static MatrixSparse* _prepareSparse(const AMesh* amesh);

private:
  VectorDouble _TildeC;
  MatrixSparse* _S;

  Id _nCovAnisoGradParam;
  VectorT<MatrixSparse*> _SGrad;
  VectorT<MatrixSparse*> _TildeCGrad;
  VectorVectorDouble _LambdaGrad;
  bool _flagNoStatByHH;
  std::vector<double> _detHH;
  mutable VectorDouble _diag;

  Id _ndim;
};

} // namespace gstlrn