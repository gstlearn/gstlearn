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

#include "Enum/EPowerPT.hpp"

#include "LinearOp/ALinearOp.hpp"
#include "Mesh/AMesh.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorT.hpp"
#include "Model/ANoStat.hpp"

#include "Matrix/MatrixSparse.hpp"

#include <map>

class Model;
class CovAniso;
class NoStatArray;
class EConsElem;
class AMatrix;
class AMatrixSquare;
class MatrixSquareGeneral;
class MatrixRectangular;
class MatrixSquareSymmetric;

/**
 * \brief Shift Operator for performing the basic tasks of SPDE
 */
class GSTLEARN_EXPORT ShiftOpCs: public ALinearOp
{

public:
  ShiftOpCs(const CGParam params = CGParam());
  ShiftOpCs(const AMesh* amesh,
            Model* model,
            const Db* dbout = nullptr,
            int igrf = 0,
            int icov = 0,
            const CGParam params = CGParam(),
            bool verbose = false);
#ifndef SWIG
  ShiftOpCs(const MatrixSparse* S,
            const VectorDouble& TildeC,
            const VectorDouble& Lambda,
            Model* model,
            const CGParam params = CGParam(),
            bool verbose = false);
#endif
  ShiftOpCs(const ShiftOpCs &shift);
  ShiftOpCs& operator=(const ShiftOpCs &shift);
  virtual ~ShiftOpCs();

  void _evalDirect(const VectorDouble& inv, VectorDouble& outv) const override;

  static ShiftOpCs* create(const AMesh *amesh,
                           Model *model,
                           const Db *dbout = nullptr,
                           int igrf = 0,
                           int icov = 0,
                           const CGParam params = CGParam(),
                           bool verbose = false);
#ifndef SWIG
  static ShiftOpCs* createFromSparse(const MatrixSparse *S,
                                     const VectorDouble &TildeC,
                                     const VectorDouble &Lambda,
                                     Model *model,
                                     const CGParam params = CGParam(),
                                     bool verbose = false);
#endif
  int initFromMesh(const AMesh* amesh,
                   Model* model,
                   const Db* dbout = nullptr,
                   int igrf = 0,
                   int icov = 0,
                   bool flagAdvection = false,
                   bool verbose = false);
  int initGradFromMesh(const AMesh* amesh,
                       Model* model,
                       int igrf = 0,
                       int icov = 0,
                       bool verbose = false,
                       double tol = EPSILON10);
  int initFromCS(const MatrixSparse* S,
                 const VectorDouble& TildeC,
                 const VectorDouble& Lambda,
                 Model* model,
                 bool verbose = false);
  int getSize() const override { return _napices; }
  int getNDim() const { return _ndim; }
  int getNModelGradParam() const { return _nModelGradParam; }
  void prodTildeC(const VectorDouble& x,
                  VectorDouble& y,
                  const EPowerPT& power) const;
  void prodLambda(const VectorDouble& x,
                  VectorDouble& y,
                  const EPowerPT& power) const;
  void prodLambdaOnSqrtTildeC(const VectorDouble& inv,
                              VectorDouble& outv,
                              double puis = 2) const;
  double getMaxEigenValue() const;
  MatrixSparse* getS() const { return _S; }
  MatrixSparse* getTildeCGrad(int iapex, int igparam) const;
  MatrixSparse* getSGrad(int iapex, int igparam) const;
  NF_Triplet getSToTriplet(bool flag_from_1 = false) const;
  NF_Triplet getTildeCGradToTriplet(int iapex, int igparam, bool flag_from_1 = false) const;
  NF_Triplet getSGradToTriplet(int iapex, int igparam, bool flag_from_1 = false) const;

  const VectorDouble& getTildeC() const { return _TildeC; }
  const VectorDouble& getLambdas() const { return _Lambda; }
  double getLambda(int iapex) const { return _Lambda[iapex]; }
  const VectorDouble& getLambdaGrads(int idim) const { return _LambdaGrad[idim]; }
  double getLambdaGrad(int idim,int iapex) const { return _LambdaGrad[idim][iapex]; }
  int getSGradAddress(int iapex, int igparam) const;

  int  getLambdaGradSize() const;

private:
  int  _getIcov() const { return _icov; }
  void _setIcov(int icov) { _icov = icov; }
  int  _getIgrf() const { return _igrf; }
  void _setIgrf(int igrf) { _igrf = igrf; }
  const Model* _getModel() const { return _model; }
  void _setModel(const Model* model) { _model = model; }
  bool _isNoStat();
  bool _isGlobalHH(int igrf, int icov);
  const CovAniso* _getCova();

  int  _buildS(const AMesh *amesh, double tol = EPSILON10);
  int  _buildSGrad(const AMesh *amesh, double tol = EPSILON10);
  void _buildLambda(const AMesh *amesh);
  bool _buildLambdaGrad(const AMesh *amesh);

  void _loadAux(VectorDouble &tab, const EConsElem &type, int imesh = 0);
  void _loadHH(const AMesh *amesh, MatrixSquareSymmetric &hh, int imesh = 0);
  void _loadHHRegular(MatrixSquareSymmetric &hh, int imesh);
  void _loadHHVariety(MatrixSquareSymmetric& hh, int imesh);
  void _loadHHGrad(const AMesh *amesh,
                   MatrixSquareSymmetric &hh,
                   int igparam,
                   int ipref);
  double _computeGradLogDetHH(const AMesh* amesh, int igparam,int ipref,
                              const MatrixSquareSymmetric& HH,
                              MatrixSquareSymmetric& work,
                              MatrixSquareSymmetric& work2);

  void _reset();
  void _resetGrad();
  void _reallocate(const ShiftOpCs& shift);
  void _projectMesh(const AMesh *amesh,
                    const VectorDouble& srot,
                    int imesh,
                    double coeff[3][2]);
  int _preparMatrices(const AMesh *amesh,
                      int imesh,
                      MatrixSquareGeneral& matu,
                      MatrixRectangular& matw) const;
  int _prepareMatricesSVariety(const AMesh *amesh,
                               int imesh,
                               VectorVectorDouble &coords,
                               MatrixRectangular& matM,
                               MatrixSquareSymmetric& matMtM,
                               AMatrix &matres,
                               double *deter);
  int _prepareMatricesSphere(const AMesh *amesh,
                             int imesh,
                             VectorVectorDouble &coords,
                             AMatrixSquare &matres,
                             double *deter);
  void _updateCova(CovAniso* cova, int imesh);
  VectorT<std::map<int, double>> _mapCreate() const;
  VectorT<VectorT<std::map<int, double>>> _mapVectorCreate() const;
  VectorT<std::map<int,double>> _mapTildeCCreate()const;
  void _mapTildeCUpdate(std::map<int, double>& tab,
                        int ip1,
                        double value,
                        double tol = EPSILON10) const;

  void _mapGradUpdate(std::map<std::pair<int, int>, double> &tab,
                      int ip0,
                      int ip1,
                      double value,
                      double tol = EPSILON10);
  MatrixSparse* _BuildTildeCGradfromMap(std::map< int, double> &tab) const;
  MatrixSparse* _BuildSGradfromMap(std::map<std::pair<int, int>, double> &tab);

  bool _cond(int indref, int igparam, int ipref);
  void _determineFlagNoStatByHH();
  void _updateHH(MatrixSquareSymmetric& hh, int imesh);
  MatrixSparse* _prepareSparse(const AMesh *amesh) const;

private:
  VectorDouble            _TildeC;
  VectorDouble            _Lambda;
  MatrixSparse*           _S;

  int                     _nModelGradParam;
  VectorT<MatrixSparse *> _SGrad;
  VectorT<MatrixSparse *> _TildeCGrad;
  VectorVectorDouble      _LambdaGrad;
  bool                    _flagNoStatByHH;

  // Following list of members are there to ease the manipulation and reduce argument list
  const Model* _model;
  int _igrf;
  int _icov;
  int _ndim;
  int _napices;
};
