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

#include "LinearOp/ALinearOp.hpp"
#include "Mesh/AMesh.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Basic/Vector.hpp"
#include "Model/ANoStat.hpp"
#include "LinearOp/EPowerPT.hpp"
#include <map>

class Model;
class CovAniso;
class NoStatArray;

/**
 * \brief Shift Operator for performing the basic tasks of SPDE
 */
class ShiftOpCs: public ALinearOp
{

public:
  ShiftOpCs();
  ShiftOpCs(AMesh* amesh,
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
  int initFromMesh(AMesh* amesh,
                   Model* model,
                   const Db* dbout = nullptr,
                   int igrf = 0,
                   int icov = 0,
                   bool flagAdvection = false,
                   bool verbose = false);
  int initGradFromMesh(AMesh* amesh,
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
  int getSize() const override { return _S->n; }
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

  cs* getS() const { return _S; }
  cs *getSGrad(int iapex, int igparam) const;
  const VectorDouble& getTildeC() const { return _TildeC; }
  const VectorDouble& getLambda() const { return _Lambda; }
  double getLambda(int iapex) const { return _Lambda[iapex]; }
  const VectorDouble& getLambdaGrad(int idim) const { return _LambdaGrad[idim]; }
  double getLambdaGrad(int idim,int iapex) const { return _LambdaGrad[idim][iapex]; }
  int getSGradAddress(int iapex, int igparam) const;

  bool getFlagGradByHh() const { return _flagGradByHH; }
  void setFlagGradByHh(bool flagGradByHh) { _flagGradByHH = flagGradByHh; }

private:
  int _getIcov() const { return _icov; }
  void _setIcov(int icov) { _icov = icov; }
  int _getIgrf() const { return _igrf; }
  void _setIgrf(int igrf) { _igrf = igrf; }
  Model* _getModel() const { return _model; }
  void _setModel(Model* model) { _model = model; }
  bool _isNoStat();
  bool _isNoStatByHH();
  bool _isVelocity();
  const CovAniso* _getCova();

  int  _buildS(AMesh *amesh, double tol = EPSILON10);
  int  _buildSVel(AMesh *amesh, double tol = EPSILON10);
  int  _buildSSphere(AMesh *amesh, double tol = EPSILON10);
  int  _buildSGrad(AMesh *amesh, double tol = EPSILON10);
  int  _buildTildeC(AMesh *amesh, const VectorDouble& units);
  void _buildLambda(AMesh *amesh);
  bool _buildLambdaGrad(AMesh *amesh);
  void _loadHHByApex(MatrixSquareSymmetric& hh, int ip);
  void _loadHHGradByApex(MatrixSquareSymmetric& hh,
                         int igparam,
                         int ip);
  void _loadAux(VectorDouble& tab,
                const EConsElem& type,
                int ip);
  void _loadHHPerMesh(MatrixSquareSymmetric& hh,
                      AMesh* amesh,
                      int imesh = 0);
  void _loadHHGradPerMesh(MatrixSquareSymmetric& hh,
                          AMesh* amesh,
                          int igp0,
                          int igparam,
                          int imesh = 0);
  void _loadAuxPerMesh(VectorDouble& tab,
                       AMesh* amesh,
                       const EConsElem& type,
                       int imesh = 0);
  void _reset();
  void _resetGrad();
  void _reallocate(const ShiftOpCs& shift);
  void _projectMesh(AMesh *amesh,
                    const VectorDouble& srot,
                    int imesh,
                    double coeff[3][2]);
  int _preparMatrices(AMesh *amesh,
                      int imesh,
                      MatrixSquareGeneral& matu,
                      MatrixRectangular& matw) const;
  cs* _BuildSfromMap(std::map<std::pair<int, int>, double> &tab);
  void _updateCova(CovAniso* cova, int ip);
  void _updateHH(MatrixSquareSymmetric& hh, int ip);
  void _mapUpdate(std::map<std::pair<int, int>, double>& tab, int ip1, int ip2, double vald, double tol=EPSILON10);

private:
  VectorDouble _TildeC;
  VectorDouble _Lambda;
  cs* _S;
  int _nModelGradParam;
  std::vector<cs *> _SGrad;
  std::vector<VectorDouble> _LambdaGrad;
  bool _flagGradByHH;

  // Following list of members are there to ease the manipulation and reduce argument list
  Model* _model;
  int _igrf;
  int _icov;
  int _ndim;
};
