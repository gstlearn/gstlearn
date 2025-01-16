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

#include "Enum/ECov.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/CovContext.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Drifts/DriftList.hpp"
#include "Basic/ICloneable.hpp"

class Model;
class Db;
class CovInternal;
class CovCalcMode;
class CovAnisoList;
class Vario;
class ADrift;
class AnamContinuous;
class AnamHermite;

typedef std::vector<ECov> VectorECov;

/**
 * \brief
 * Class containing the Model Information describing the formal Spatial (or Temporal) Characteristics
 * of the (set of) random variable(s) under study.
 *
 * The Model is essentially a container with two main contents:
 * - the **covariance** part: see ACov.hpp for more information
 * - the **drift** part: see DriftList.hpp for more information
 *
 * The additional member **CovContext** only serves in carrying the following information:
 * - the number of variables: if more than 1, the Model becomes multivariate
 * - the field extension: this information is needed to get a *stationary* version to any covariance
 * - the experimental mean vector and the variance-covariance matrix (used to calibrate the Model)
 */
class GSTLEARN_EXPORT ModelGeneric : public ICloneable
{
public:
  ModelGeneric(const CovContext& ctxt = CovContext());
  ModelGeneric(const ModelGeneric &m) = delete;
  ModelGeneric& operator= (const ModelGeneric &m) = delete;
  virtual ~ModelGeneric();

  MatrixRectangular evalDriftMatrix(const Db* db,
                                    int ivar0             = -1,
                                    const VectorInt& nbgh = VectorInt(),
                                    const ECalcMember& member = ECalcMember::fromKey("LHS")) const;
  MatrixRectangular evalCovMatrix(Db* db1,
                                  Db* db2                 = nullptr,
                                  int ivar0               = -1,
                                  int jvar0               = -1,
                                  const VectorInt& nbgh1  = VectorInt(),
                                  const VectorInt& nbgh2  = VectorInt(),
                                  const CovCalcMode* mode = nullptr);
  MatrixRectangular evalCovMatrixOptim(Db* db1,
                                       Db* db2                 = nullptr,
                                       int ivar0               = -1,
                                       int jvar0               = -1,
                                       const VectorInt& nbgh1  = VectorInt(),
                                       const VectorInt& nbgh2  = VectorInt(),
                                       const CovCalcMode* mode = nullptr,
                                       bool cleanOptim         = true);
  MatrixSquareSymmetric evalCovMatrixSymmetric(const Db *db1,
                                               int ivar0 = -1,
                                               const VectorInt &nbgh1 = VectorInt(),
                                               const CovCalcMode *mode = nullptr);
  MatrixSquareSymmetric
  evalCovMatrixSymmetricOptim(const Db* db1,
                              int ivar0               = -1,
                              const VectorInt& nbgh1  = VectorInt(),
                              const CovCalcMode* mode = nullptr,
                              bool cleanOptim         = true);
  
  MatrixSquareGeneral eval0Mat(const CovCalcMode* mode = nullptr) const;

  MatrixSparse* evalCovMatrixSparse(Db* db1,
                                    Db* db2                 = nullptr,
                                    int ivar0               = -1,
                                    int jvar0               = -1,
                                    const VectorInt& nbgh1  = VectorInt(),
                                    const VectorInt& nbgh2  = VectorInt(),
                                    const CovCalcMode* mode = nullptr,
                                    double eps              = EPSILON3);

  void setField(double field);
  bool isValid() const;

  // Pipes for the private members
  // Case of _ctxt
  int getVariableNumber() const;
  int getDimensionNumber() const;
  const VectorDouble& getMeans() const;

  // Case of _driftList
  int getDriftNumber() const;
  int getDriftEquationNumber() const;
  int getExternalDriftNumber() const;
  int getDriftMaxIRFOrder(void) const;
  void delAllDrifts();

    // Case of _cova
  const CovAnisoList* getCovAnisoList() const;
  CovAnisoList* getCovAnisoListModify() const;
  int getCovaMinIRFOrder() const;
  int getCovaNumber(bool skipNugget = false) const;
  void setActiveFactor(int iclass);
  int  getActiveFactor() const;

protected:               // TODO : pass into private to finish clean
  ACov* _cova;           /* Generic Covariance structure */
  DriftList* _driftList; /* Series of Drift functions */
  CovContext _ctxt;      /* Context */
};
