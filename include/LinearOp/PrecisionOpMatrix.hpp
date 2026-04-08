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

#include "LinearOp/ASimulableMatrix.hpp"
#include "LinearOp/CholeskySparse.hpp"
#include "LinearOp/PrecisionOp.hpp"

namespace gstlrn
{
  class AMesh;
  class CovAniso;
  class Model;
  class ShiftOpMatrix;

  /** This class is just a specialization of PrecisionOp when the shift
   * Operator is built with sparse matrices and therefore algebra can be performed with Cholesky.
   * It allows to return the precision matrix as a Sparse Matrix. */
  class GSTLEARN_EXPORT PrecisionOpMatrix: virtual public PrecisionOp,
                                           virtual public ASimulableMatrix
  {
  public:
    PrecisionOpMatrix(
      ShiftOpMatrix* shiftop = nullptr,
      const CovAniso* cova = nullptr,
      bool verbose = false);
    PrecisionOpMatrix(const AMesh* mesh, CovAniso* cova, bool verbose = false);
    PrecisionOpMatrix(const PrecisionOpMatrix& pmat) = default;
    PrecisionOpMatrix(PrecisionOpMatrix&& pmat) noexcept;
    PrecisionOpMatrix& operator=(const PrecisionOpMatrix& pmat) = default;
    PrecisionOpMatrix& operator=(PrecisionOpMatrix&& pmat) noexcept;
    virtual ~PrecisionOpMatrix();

    Id getSize() const override;

    const MatrixSparse& getQMat() const override { return *getQ(); }

    // Interface for PrecisionOp class
#ifndef SWIG
    void evalInverse(const constvect vecin, VectorDouble& vecout) override;
    Id _addSimulateToDest(const constvect whitenoise, vect outv) const override;
    Id _addToDest(const constvect inv, vect outv) const override;
#endif

    double computeLogDet(Id nMC = 1) const override
    {
      return ASimulableMatrix::computeLogDet(nMC);
    }

    VectorDouble extractDiag() const override;

    // void evalDerivPoly(const VectorDouble& inv, VectorDouble& outv,Id iapex,Id igparam) override;
#ifndef SWIG
    void evalDeriv(
      const constvect inv,
      vect outv,
      Id iapex,
      Id igparam,
      const EPowerPT& power) override;
    void evalDerivOptim(vect outv, Id iapex, Id igparam, const EPowerPT& power)
      override;
    void gradYQX(
      const constvect X,
      const constvect Y,
      vect result,
      const EPowerPT& power) override;
    void gradYQXOptim(
      const constvect X,
      const constvect Y,
      vect result,
      const EPowerPT& power) override;
#endif
    const MatrixSparse* getQ() const { return _Q.get(); }

    const MatrixSparse* getS() const;
    const VectorDouble& getTildeC() const;

  private:
    void _buildQ();
    MatrixSparse* _build_Q();

  private:
    std::shared_ptr<MatrixSparse> _Q;
    mutable CholeskySparse* _chol;
  };

} // namespace gstlrn
