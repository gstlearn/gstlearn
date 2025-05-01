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

#include "Enum/EPowerPT.hpp"
#include "LinearOp/AShiftOp.hpp"
#include "Mesh/AMesh.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorT.hpp"

#include "Matrix/MatrixSparse.hpp"

#include <map>
#include <memory>

#ifndef SWIG

#include <Eigen/Core>
#include <Eigen/Dense>
#endif

class CovAniso;
class EConsElem;
class AMatrix;
class MatrixSquare;
class MatrixDense;
class MatrixSymmetric;


class GSTLEARN_EXPORT ShiftOpMatrix: public AShiftOp
{
  public:
    ShiftOpMatrix();
    ShiftOpMatrix(const AMesh* amesh, const CovAniso* cova, const Db* dbout = nullptr, bool verbose = false);
    ShiftOpMatrix(const MatrixSparse* S, const VectorDouble& TildeC,
              const VectorDouble& Lambda, const CovAniso* cova, bool verbose = false);
    ShiftOpMatrix(const ShiftOpMatrix& shift);
    ShiftOpMatrix& operator=(const ShiftOpMatrix& shift);
    virtual ~ShiftOpMatrix();
    /// ICloneable interface
    IMPLEMENT_CLONING(ShiftOpMatrix)
    void normalizeLambdaBySills(const AMesh* mesh) override;
#ifndef SWIG
    int _addToDest(const constvect inv, vect outv) const override;
#endif

    static ShiftOpMatrix* create(const AMesh* amesh, const CovAniso* cova,
                             const Db* dbout = nullptr, 
                             bool verbose = false);
    static ShiftOpMatrix* createFromSparse(
      const MatrixSparse* S, const VectorDouble& TildeC,
      const VectorDouble& Lambda, const CovAniso* cova, bool verbose = false);
    int initFromMesh(const AMesh* amesh, const CovAniso* cova,
                     const Db* dbout = nullptr, 
                     bool flagAdvection = false, bool verbose = false);
    int initGradFromMesh(const AMesh* amesh, const CovAniso* cova,
                         bool verbose = false,
                         double tol = EPSILON10);
    int initFromCS(const MatrixSparse* S, const VectorDouble& TildeC,
                   const VectorDouble& Lambda, const CovAniso* cova,
                   bool verbose = false);

    int getNDim() const
    {
      return _ndim;
    }
    int getNCovAnisoGradParam() const
    {
      return _nCovAnisoGradParam;
    }
    void prodTildeC(const VectorDouble& x, VectorDouble& y,
                    const EPowerPT& power) const;
  
    void prodLambdaOnSqrtTildeC(const VectorDouble& inv, VectorDouble& outv,
                                double puis = 2) const;
    double getMaxEigenValue() const override;
    MatrixSparse* getS() const { return _S; }
    MatrixSparse* getTildeCGrad(int iapex, int igparam) const;
    MatrixSparse* getSGrad(int iapex, int igparam) const;

    const VectorDouble& getTildeC() const
    {
      return _TildeC;
    }

    const VectorDouble& getLambdaGrads(int idim) const
    {
      return _LambdaGrad[idim];
    }
    double getLambdaGrad(int idim, int iapex) const
    {
      return _LambdaGrad[idim][iapex];
    }
    int getSGradAddress(int iapex, int igparam) const;

    int getLambdaGradSize() const;
    //void multiplyByValueAndAddDiagonal(double v1 = 1.,double v2 = 0.) override;
  private:

  
    int _buildS(const AMesh* amesh, double tol = EPSILON10);
    int _buildSGrad(const AMesh* amesh, double tol = EPSILON10);
    void _buildLambda(const AMesh* amesh);
    bool _buildLambdaGrad(const AMesh* amesh);

    static void _loadAux(VectorDouble & tab, const EConsElem& type, int imesh = 0);
    void _loadHH(const AMesh* amesh, MatrixSymmetric& hh, int imesh = 0);
    void _loadHHRegular(MatrixSymmetric & hh, int imesh);
    void _loadHHVariety(MatrixSymmetric & hh, int imesh);
    void _loadHHGrad(const AMesh* amesh, MatrixSymmetric& hh, int igparam,
                     int ipref);
    double _computeGradLogDetHH(const AMesh* amesh, int igparam, int ipref,
                                const MatrixSymmetric& HH,
                                MatrixSymmetric& work,
                                MatrixSymmetric& work2);

    void _reset();
    int _resetGrad();
    void _reallocate(const ShiftOpMatrix& shift);
    static void _projectMesh(const AMesh* amesh,
                             const VectorDouble& srot,
                             int imesh,
                             double coeff[3][2]);
    int _preparMatrices(const AMesh* amesh, int imesh,
                        MatrixSquare& matu, MatrixDense& matw)
      const;
    int _prepareMatricesSVariety(const AMesh* amesh,
                                 int imesh,
                                 VectorVectorDouble& coords,
                                 MatrixDense& matM,
                                 MatrixSymmetric& matMtM,
                                 AMatrix& matP,
                                 double* deter) const;
    int _prepareMatricesSphere(const AMesh* amesh,
                               int imesh,
                               VectorVectorDouble& coords,
                               MatrixSquare& matMs,
                               double* deter) const;
    static void _updateCova(std::shared_ptr<CovAniso> &cova, int imesh);
    VectorT<std::map<int, double>> _mapCreate() const;
    VectorT<VectorT<std::map<int, double>>> _mapVectorCreate() const;
    VectorT<std::map<int, double>> _mapTildeCCreate() const;
    static void _mapTildeCUpdate(std::map<int, double>& tab,
                                 int ip0,
                                 double value,
                                 double tol = EPSILON10);

    static void _mapGradUpdate(std::map<std::pair<int, int>, double>& tab,
                               int ip0,
                               int ip1,
                               double value,
                               double tol = EPSILON10);
    MatrixSparse* _BuildTildeCGradfromMap(std::map<int, double> & tab) const;
    MatrixSparse* _BuildSGradfromMap(std::map<std::pair<int, int>, double> &
                                     tab) const;

    static bool _cond(int indref, int igparam, int ipref);
    void _determineFlagNoStatByHH();
    void _updateHH(MatrixSymmetric & hh, int imesh);
    static MatrixSparse* _prepareSparse(const AMesh* amesh);

  private:
    VectorDouble _TildeC;
    MatrixSparse* _S;

    int _nCovAnisoGradParam;
    VectorT<MatrixSparse*> _SGrad;
    VectorT<MatrixSparse*> _TildeCGrad;
    VectorVectorDouble _LambdaGrad;
    bool _flagNoStatByHH;
    std::vector<double> _detHH;


    int _ndim;
  };

