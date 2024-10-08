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
class AMatrixSquare;
class MatrixSquareGeneral;
class MatrixRectangular;
class MatrixSquareSymmetric;

/**
 * \brief Shift Operator for performing the basic tasks of SPDE
 */

#ifndef SWIG
#  include "LinearOp/ALinearOpEigenCG.hpp"
DECLARE_EIGEN_TRAITS(ShiftOpCs)
#else
#  include "LinearOp/ALinearOp.hpp"
#endif

class GSTLEARN_EXPORT ShiftOpCs:
#ifndef SWIG
  public ALinearOpEigenCG<ShiftOpCs>
#else
  public ALinearOp
#endif
{
  public:
    ShiftOpCs();
    ShiftOpCs(const AMesh* amesh, const CovAniso* cova, const Db* dbout = nullptr, bool verbose = false);
    ShiftOpCs(const MatrixSparse* S, const VectorDouble& TildeC,
              const VectorDouble& Lambda, const CovAniso* cova, bool verbose = false);
    ShiftOpCs(const ShiftOpCs& shift);
    ShiftOpCs& operator=(const ShiftOpCs& shift);
    virtual ~ShiftOpCs();
    void normalizeLambdaBySills(const AMesh*);
#ifndef SWIG
    int _addToDest(const constvect inv, vect outv) const override;
#endif

    static ShiftOpCs* create(const AMesh* amesh, const CovAniso* cova,
                             const Db* dbout = nullptr, 
                             bool verbose = false);
    static ShiftOpCs* createFromSparse(
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
    int getSize() const override
    {
      return _napices;
    }

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
    void prodLambda(const VectorDouble& x, VectorDouble& y,
                    const EPowerPT& power) const;
  #ifndef SWIG
    void prodLambda(const constvect x, vect y, const EPowerPT& power) const;
    void prodLambda(const VectorDouble& x, vect y, const EPowerPT& power) const;
    void
    prodLambda(const constvect x, VectorDouble& y, const EPowerPT& power) const;
#endif
    void prodLambdaOnSqrtTildeC(const VectorDouble& inv, VectorDouble& outv,
                                double puis = 2) const;
    double getMaxEigenValue() const;
    MatrixSparse* getS() const
    {
      return _S;
    }
    MatrixSparse* getTildeCGrad(int iapex, int igparam) const;
    MatrixSparse* getSGrad(int iapex, int igparam) const;

    const VectorDouble& getTildeC() const
    {
      return _TildeC;
    }
    const VectorDouble& getLambdas() const
    {
      return _Lambda;
    }
    double getLambda(int iapex) const
    {
      return _Lambda[iapex];
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

  private:
    void _setCovAniso(const CovAniso* cova);
    bool _isNoStat();
    bool _isGlobalHH();
  
    std::shared_ptr<CovAniso>& _getCovAniso();
    int _buildS(const AMesh* amesh, double tol = EPSILON10);
    int _buildSGrad(const AMesh* amesh, double tol = EPSILON10);
    void _buildLambda(const AMesh* amesh);
    bool _buildLambdaGrad(const AMesh* amesh);

    void _loadAux(VectorDouble & tab, const EConsElem& type, int imesh = 0);
    void _loadHH(const AMesh* amesh, MatrixSquareSymmetric& hh, int imesh = 0);
    void _loadHHRegular(MatrixSquareSymmetric & hh, int imesh);
    void _loadHHVariety(MatrixSquareSymmetric & hh, int imesh);
    void _loadHHGrad(const AMesh* amesh, MatrixSquareSymmetric& hh, int igparam,
                     int ipref);
    double _computeGradLogDetHH(const AMesh* amesh, int igparam, int ipref,
                                const MatrixSquareSymmetric& HH,
                                MatrixSquareSymmetric& work,
                                MatrixSquareSymmetric& work2);

    void _reset();
    int _resetGrad();
    void _reallocate(const ShiftOpCs& shift);
    static void _projectMesh(const AMesh* amesh,
                             const VectorDouble& srot,
                             int imesh,
                             double coeff[3][2]);
    int _preparMatrices(const AMesh* amesh, int imesh,
                        MatrixSquareGeneral& matu, MatrixRectangular& matw)
      const;
    int _prepareMatricesSVariety(const AMesh* amesh,
                                 int imesh,
                                 VectorVectorDouble& coords,
                                 MatrixRectangular& matM,
                                 MatrixSquareSymmetric& matMtM,
                                 AMatrix& matP,
                                 double* deter) const;
    int _prepareMatricesSphere(const AMesh* amesh,
                               int imesh,
                               VectorVectorDouble& coords,
                               AMatrixSquare& matMs,
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
    void _updateHH(MatrixSquareSymmetric & hh, int imesh);
    static MatrixSparse* _prepareSparse(const AMesh* amesh);
    static std::shared_ptr<CovAniso> cloneAndCast(const CovAniso* cova);
    static std::shared_ptr<CovAniso> cloneAndCast(const std::shared_ptr<CovAniso> &cova);

  private:
    VectorDouble _TildeC;
    VectorDouble _Lambda;
    MatrixSparse* _S;

    int _nCovAnisoGradParam;
    VectorT<MatrixSparse*> _SGrad;
    VectorT<MatrixSparse*> _TildeCGrad;
    VectorVectorDouble _LambdaGrad;
    bool _flagNoStatByHH;

    // Following list of members are there to ease the manipulation and reduce
    // argument list
    std::shared_ptr<CovAniso> _cova;

    int _ndim;
    int _napices;
  };

#ifndef SWIG
  DECLARE_EIGEN_PRODUCT(ShiftOpCs)
#endif
