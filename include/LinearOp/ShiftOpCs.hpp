/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* Created on: 9 avr. 2019 by N. Desassis                                     */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/EPowerPT.hpp"

#include "LinearOp/ALinearOp.hpp"
#include "Mesh/AMesh.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Basic/Vector.hpp"
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

  void _evalDirect(const VectorDouble& in, VectorDouble& out) const override;

  int initFromOldMesh(SPDE_Mesh* s_mesh,
                      Model* model,
                      Db* dbout = nullptr,
                      bool flagAdvection = false,
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
                       Db* dbout,
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
  void prodTildeC(const VectorDouble& in,
                  VectorDouble& out,
                  const EPowerPT& power) const;
  void prodLambda(const VectorDouble& in,
                  VectorDouble& out,
                  const EPowerPT& power) const;
  void prodLambdaOnSqrtTildeC(const VectorDouble& out,
                              VectorDouble& in,
                              double puis = 2) const;
  double getMaxEigenValue() const;
  int getVariety()const {return _variety;}
  cs* getS() const { return _S; }
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
                         int ip);

  void _loadHHGradPerMesh(MatrixSquareSymmetric& hh,
                          const AMesh* amesh,
                          int igp0,
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
  cs* _BuildSfromMap(VectorT<std::map<int, double>>& tab, int nmax = -1);
  cs* _BuildVecSfromMap(std::map<std::pair<int, int>, double>& tab,
                        int nmax = -1);
  void _updateCova(CovAniso* cova, int ip);
  void _updateHH(MatrixSquareSymmetric& hh, int ip);
  VectorT<std::map<int, double>> _mapCreate() const;
  VectorT<std::map<std::pair<int, int>, double>> _mapVectorCreate() const;
  void _mapUpdate(std::map<int, double>& tab,
                  int ip1,
                  double value,
                  double tol = EPSILON10) const;
  void _mapVecUpdate(std::map<std::pair<int, int>, double>& tab,
                     int ip2,
                     int ip1,
                     double value,
                     double tol = EPSILON10) const;

  void _determineFlagNoStatByHH();
private:
  VectorDouble _TildeC;
  VectorDouble _Lambda;
  cs* _S;
  int _nModelGradParam;
  VectorT<cs *> _SGrad;
  VectorVectorDouble _LambdaGrad;
  bool _flagNoStatByHH;
  int _variety;
  // Following list of members are there to ease the manipulation and reduce argument list
  const Model* _model;
  int _igrf;
  int _icov;
  int _ndim;
  int _napices;
};
