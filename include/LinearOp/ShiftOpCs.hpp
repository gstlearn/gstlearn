/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/EPowerPT.hpp"

#include "LinearOp/ALinearOp.hpp"
#include "Mesh/AMesh.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorT.hpp"
#include "Model/ANoStat.hpp"

#include <map>

class Model;
class CovAniso;
class NoStatArray;
class EConsElem;

/**
 * \brief Shift Operator for performing the basic tasks of SPDE
 */
class GSTLEARN_EXPORT ShiftOpCs: public ALinearOp
{

public:
  ShiftOpCs();
  ShiftOpCs(const AMesh* amesh,
            Model* model,
            const Db* dbout = nullptr,
            int igrf = 0,
            int icov = 0,
            bool verbose = false);
  ShiftOpCs(const cs* S,
            const VectorDouble& TildeC,
            const VectorDouble& Lambda,
            Model* model,
            bool verbose = false);
  ShiftOpCs(const ShiftOpCs &shift);
  ShiftOpCs& operator=(const ShiftOpCs &shift);
  virtual ~ShiftOpCs();

  void _evalDirect(const VectorDouble& inv, VectorDouble& outv) const override;

  static ShiftOpCs* create(const AMesh *amesh,
                           Model *model,
                           const Db *dbout = nullptr,
                           int igrf = 0,
                           int icov = 0,
                           bool verbose = false);
  static ShiftOpCs* createFromSparse(const cs *S,
                                     const VectorDouble &TildeC,
                                     const VectorDouble &Lambda,
                                     Model *model,
                                     bool verbose = false);

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
  int initFromCS(const cs* S,
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
  int getVariety()const {return _variety;}
  cs* getS() const { return _S; }
  cs* getTildeCGrad(int iapex, int igparam) const;
  cs* getSGrad(int iapex, int igparam) const;
  const VectorDouble& getTildeC() const { return _TildeC; }
  const VectorDouble& getLambdas() const { return _Lambda; }
  double getLambda(int iapex) const { return _Lambda[iapex]; }
  const VectorDouble& getLambdaGrads(int idim) const { return _LambdaGrad[idim]; }
  double getLambdaGrad(int idim,int iapex) const { return _LambdaGrad[idim][iapex]; }
  int getSGradAddress(int iapex, int igparam) const;

  bool getFlagNoStatByHH() const { return _flagNoStatByHH; }
  void setFlagNoStatByHH(bool flagGradByHH) { _flagNoStatByHH = flagGradByHH; }
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
  bool _isVelocity();
  const CovAniso* _getCova();

  int  _buildSVel(const AMesh *amesh, double tol = EPSILON10);
  int  _buildSVariety(const AMesh *amesh, double tol = EPSILON10);
  int  _buildSSphere(const AMesh *amesh, double tol = EPSILON10);
  int  _buildSGrad(const AMesh *amesh, double tol = EPSILON10);
  int  _buildTildeC(const AMesh *amesh, const VectorDouble& units);
  void _buildLambda(const AMesh *amesh);

  void _loadAux(VectorDouble& tab,
                const EConsElem& type,
                int ip);
  void _loadAuxPerMesh(const AMesh* amesh,
                       VectorDouble& tab,
                       const EConsElem& type,
                       int imesh = 0);
  void _loadHHPerMesh(const AMesh* amesh,
                      MatrixSquareSymmetric& hh,
                      int imesh = 0);
  void _loadHHRegularByApex(MatrixSquareSymmetric& hh, int ip);
  void _loadHHVarietyByApex(MatrixSquareSymmetric& hh, int ip);
  void _loadHHByApex(const AMesh* amesh, MatrixSquareSymmetric& hh, int ip);
  void _loadHHGradByApex(MatrixSquareSymmetric& hh,
                         int igparam,
                         int ipref);
  double _computeGradLogDetHH(const AMesh* amesh, int igparam,int ipref,
                              const MatrixSquareSymmetric& HH,
                              MatrixSquareSymmetric& work,
                              MatrixSquareSymmetric& work2);
  void _loadHHGradPerMesh(MatrixSquareSymmetric& hh,
                          const AMesh* amesh,
                          int ipref,
                          int igparam);
  bool _buildLambdaGrad(const AMesh *amesh);

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
                               MatrixSquareSymmetric& matMtM,
                               AMatrix &matres,
                               double *deter);
  int _prepareMatricesSphere(const AMesh *amesh,
                             int imesh,
                             VectorVectorDouble &coords,
                             AMatrixSquare &matres,
                             double *deter);
  cs* _BuildSfromMap(VectorT<std::map<int, double>>& tab);
  cs* _BuildVecSfromMap(std::map<std::pair<int, int>, double>& tab);
  void _updateCova(CovAniso* cova, int ip);
  void _updateHH(MatrixSquareSymmetric& hh, int ip);
  VectorT<std::map<int, double>> _mapCreate() const;
  VectorT<VectorT<std::map<int, double>>> _mapVectorCreate() const;
  VectorT<std::map<int,double>> _mapTildeCCreate()const;
  void _determineFlagNoStatByHH();
  void _mapUpdate(std::map<int, double>& tab,
                  int ip1,
                  double value,
                  double tol = EPSILON10) const;
  void _mapTildeCUpdate(std::map<int, double>& tab,
                        int ip1,
                        double value,
                        double tol = EPSILON10) const;

  void _mapGradUpdate(std::map<std::pair<int, int>, double> &tab,
                      int ip0,
                      int ip1,
                      double value,
                      double tol = EPSILON10);
  cs* _BuildTildeCGradfromMap(std::map< int, double> &tab) const;
  cs* _BuildSGradfromMap(std::map<std::pair<int, int>, double> &tab);

  bool _cond(int indref, int igparam, int ipref);

private:
  VectorDouble       _TildeC;
  VectorDouble       _Lambda;
  cs*                _S;
  int                _nModelGradParam;
  VectorT<cs *>      _SGrad;
  VectorT<cs *>      _TildeCGrad;
  VectorVectorDouble _LambdaGrad;
  bool               _flagNoStatByHH;
  int                _variety;

  // Following list of members are there to ease the manipulation and reduce argument list
  const Model* _model;
  int _igrf;
  int _icov;
  int _ndim;
  int _napices;
};
