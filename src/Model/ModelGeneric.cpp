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

#include "Model/ModelGeneric.hpp"


ModelGeneric::ModelGeneric(const CovContext &ctxt)
    : _cova(nullptr),
      _driftList(nullptr),
      _ctxt(ctxt)
{
 
}


ModelGeneric::~ModelGeneric()
{
}

MatrixRectangular ModelGeneric::evalDriftMatrix(const Db *db,
                                    int ivar0,
                                    const VectorInt& nbgh,
                                    const ECalcMember &member) const
{
    if (_driftList == nullptr) return MatrixRectangular();
    return _driftList->evalDriftMatrix(db, ivar0, nbgh, member);
}

MatrixRectangular ModelGeneric::evalCovMatrix(Db* db1,
                                  Db* db2 ,
                                  int ivar0,
                                  int jvar0,
                                  const VectorInt& nbgh1,
                                  const VectorInt& nbgh2,
                                  const CovCalcMode* mode)
  {
    if (_cova == nullptr) return MatrixRectangular();
    return _cova->evalCovMatrix(db1, db2, ivar0, jvar0, nbgh1, nbgh2, mode);
  }
  
  MatrixSquareSymmetric ModelGeneric::evalCovMatrixSymmetric(const Db *db1,
                                               int ivar0,
                                               const VectorInt &nbgh1,
                                               const CovCalcMode *mode)
  {
    if (_cova == nullptr) return MatrixSquareSymmetric();
    return _cova->evalCovMatrixSymmetric(db1, ivar0, nbgh1, mode);
  }

MatrixRectangular ModelGeneric::evalCovMatrixOptim(Db* db1,
                                  Db* db2 ,
                                  int ivar0,
                                  int jvar0,
                                  const VectorInt& nbgh1,
                                  const VectorInt& nbgh2,
                                  const CovCalcMode* mode)
  {
    if (_cova == nullptr) return MatrixRectangular();
    return _cova->evalCovMatrixOptim(db1, db2, ivar0, jvar0, nbgh1, nbgh2, mode);
  }
  
  MatrixSquareSymmetric ModelGeneric::evalCovMatrixSymmetricOptim(const Db *db1,
                                               int ivar0,
                                               const VectorInt &nbgh1,
                                               const CovCalcMode *mode)
  {
    if (_cova == nullptr) return MatrixSquareSymmetric();
    return _cova->evalCovMatrixSymmetricOptim(db1, ivar0, nbgh1, mode);
  }

  MatrixSparse* ModelGeneric::evalCovMatrixSparse(Db *db1,
                                    Db *db2,
                                    int ivar0,
                                    int jvar0,
                                    const VectorInt &nbgh1,
                                    const VectorInt &nbgh2,
                                    const CovCalcMode *mode,
                                    double eps)
  {
    if (_cova == nullptr) return nullptr;
    return _cova->evalCovMatrixSparse(db1, db2, ivar0, jvar0, nbgh1, nbgh2, mode, eps);
  }
  