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
#include "MatrixC/MatrixCSGeneral.hpp"
#include "MatrixC/MatrixCRectangular.hpp"
#include "Basic/Vector.hpp"
#include "Model/ANoStat.hpp"
#include <map>

class Model;
class CovAniso;

/**
 * \brief Shift Operator for performing the basic tasks of SPDE
 */
class ShiftOpCs: public ALinearOp
{

public:
  ShiftOpCs();
  ShiftOpCs(AMesh* amesh,
            Model* model,
            Db* dbout = nullptr,
            ANoStat* nostat = nullptr,
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
                      ANoStat* nostat = nullptr,
                      bool flagAdvection = false,
                      bool verbose = false);
  int initFromMesh(AMesh* amesh,
                   Model* model,
                   Db* dbout = nullptr,
                   ANoStat* nostat = nullptr,
                   int igrf = 0,
                   int icov = 0,
                   bool flagAdvection = false,
                   bool verbose = false);
  int initGradFromMesh(AMesh* amesh,
                       Model* model,
                       Db* dbout,
                       ANoStat* nostat = nullptr,
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

  void prodTildeC(const VectorDouble& in, VectorDouble& out, ENUM_POPTS power) const;
  void prodLambda(const VectorDouble& in, VectorDouble& out, ENUM_POPTS power) const;
  void prodLambdaOnSqrtTildeC(const VectorDouble& out,
                              VectorDouble& in,
                              double puis = 2) const;
//  void prodSGrad(int iapex, int iparam, const VectorDouble& in, VectorDouble& out) const {};
  double getMaxEigenValue() const;

  cs* getS() const { return _S; }
  cs *getSGrad(int iapex, int igparam) const;
  const VectorDouble& getTildeC() const { return _TildeC; }
  const VectorDouble& getLambda() const { return _Lambda; }
  const VectorDouble& getLambdaGrad(int idim) const { return _LambdaGrad[idim]; }

private:
  int _buildS(AMesh *amesh,
              const CovAniso& cova,
              int igrf,
              int icov,
              ANoStat* nostat,
              double tol = EPSILON10);
  int _buildSVel(AMesh *amesh,
                 const CovAniso& cova,
                 int igrf,
                 int icov,
                 ANoStat* nostat,
                 double tol = EPSILON10);
  int _buildSSphere(AMesh *amesh,
                    const CovAniso& cova,
                    int igrf,
                    int icov,
                    ANoStat* nostat,
                    double tol = EPSILON10);
  int _buildSGrad(AMesh *amesh,
                  const CovAniso& cova,
                  int igrf,
                  int icov,
                  ANoStat* nostat,
                  double tol = EPSILON10);
  int  _buildTildeC(AMesh *amesh, const VectorDouble& units);
  void _buildLambda(AMesh *amesh,
                    const CovAniso& cova,
                    int igrf,
                    int icov,
                    ANoStat* nostat);
  bool _buildLambdaGrad(AMesh *amesh,
                      CovAniso& cova,
                      int igrf,
                      int icov,
                      ANoStat* nostat);
  void _loadHHByApex(MatrixCSGeneral& hh,
                     const CovAniso& cova,
                     int igrf = 0,
                     int icov = 0,
                     int ip = 0,
                     ANoStat* nostat = nullptr);
  void _loadHHGradByApex(MatrixCSGeneral& hh,
                         const CovAniso& covini,
                         int igparam,
                         int igrf,
                         int icov,
                         int ip,
                         ANoStat* nostat,
                         bool flagFormal = true);
  void _loadAux(VectorDouble& tab,
                int igrf,
                int icov,
                ENUM_CONS type,
                int ip = 0,
                ANoStat *nostat = nullptr);
  void _loadHHPerMesh(MatrixCSGeneral& hh,
                      AMesh* amesh,
                      const CovAniso& cova,
                      int igrf = 0,
                      int icov = 0,
                      int imesh = 0,
                      ANoStat* nostat = nullptr);
  void _loadHHGradPerMesh(MatrixCSGeneral& hh,
                          AMesh* amesh,
                          const CovAniso& cova,
                          int igp0,
                          int igparam,
                          int igrf = 0,
                          int icov = 0,
                          int imesh = 0,
                          ANoStat* nostat = nullptr,
                          bool flagFormal = true);
  void _loadAuxPerMesh(VectorDouble& tab,
                       AMesh* amesh,
                       int igrf,
                       int icov,
                       ENUM_CONS type,
                       int imesh = 0,
                       ANoStat* nostat = nullptr);
  void _reset();
  void _resetGrad();
  void _reallocate(const ShiftOpCs& shift);
  void _projectMesh(AMesh *amesh,
                    const VectorDouble& srot,
                    int imesh,
                    double coeff[3][2]);
  int _preparMatrices(AMesh *amesh,
                      int imesh,
                      MatrixCSGeneral& matu,
                      MatrixCRectangular& matw) const;
  cs* _BuildSfromMap(std::map<std::pair<int, int>, double> &tab);
  int _getSGradAddress(int iapex, int igparam) const;
  void _updateCova(CovAniso& cova, int igrf, int icov, int ip, int ndim, ANoStat* nostat);
  void _mapUpdate(std::map<std::pair<int, int>, double>& tab, int ip1, int ip2, double vald, double tol=EPSILON10);

private:
  VectorDouble _TildeC;
  VectorDouble _Lambda;
  cs* _S;
  int _nModelGradParam;
  std::vector<cs *> _SGrad;
  std::vector<VectorDouble> _LambdaGrad;
};
