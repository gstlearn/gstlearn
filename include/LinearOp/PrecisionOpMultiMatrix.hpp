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
#include "LinearOp/PrecisionOpMulti.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "geoslib_define.h"

namespace gstlrn
{
  class Model;

  /**
   * Class for the precision matrix of the latent field in SPDE (matricial form)
   */
  class GSTLEARN_EXPORT PrecisionOpMultiMatrix: virtual public PrecisionOpMulti,
                                                virtual public ASimulableMatrix
  {
  public:
    PrecisionOpMultiMatrix(
      Model* model = nullptr,
      const VectorMeshes& meshes = VectorMeshes());
    PrecisionOpMultiMatrix(const PrecisionOpMulti& m) = delete;
    PrecisionOpMultiMatrix& operator=(const PrecisionOpMulti& m) = delete;
    virtual ~PrecisionOpMultiMatrix();

    const MatrixSparse& getQMat() const override;
    const MatrixSparse* getQ() const;
    double computeLogDet(Id nMC = 1) const override;
    Id getSize() const override;

  protected:
    Id _addToDest(const constvect whitenoise, vect outv) const override;
    Id _addSimulateToDest(const constvect whitenoise, vect outv) const override;

  private:
    MatrixSparse _prepareMatrixNoStat(Id icov, const MatrixSparse* Q) const;
    MatrixSparse _prepareMatrixStationary(Id icov, const MatrixSparse* Q) const;
    void _prepareMatrix();
    void _buildQop(bool stencil = false) override;

    bool _isSingle() const { return _getNVar() == 1 && getNCov() == 1; }

  private:
    MatrixSparse _Q;
  };
} // namespace gstlrn
