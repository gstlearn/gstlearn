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

#include "geoslib_define.h"
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

  //getters for member pointers
  const ACov*       getCov()             const { return  _cova;     }
  const CovContext* getContext()         const { return &_ctxt;     }
  const DriftList*  getDriftList()       const { return  _driftList;}
  ACov*             getCovModify()             { return  _cova;     }
  CovContext*       getContextModify()         { return &_ctxt;     }
  DriftList*        getDriftListModify()       { return  _driftList;}
  
  MatrixRectangular evalDriftMat(const Db* db,
                                    int ivar0             = -1,
                                    const VectorInt& nbgh = VectorInt(),
                                    const ECalcMember& member = ECalcMember::fromKey("LHS")) const;
  MatrixRectangular evalDriftMatByRanks(const Db* db,
                    const VectorVectorInt& sampleRanks,
                    int ivar0                 = -1,
                    const ECalcMember& member = ECalcMember::fromKey("LHS")) const;

  MatrixRectangular evalDriftMatByTarget(const Db* db,
                                          int ivar0 = -1,
                                          int iech2 = 0,
                                          const ECalcMember& member = ECalcMember::fromKey("LHS")) const;
  MatrixRectangular evalCovMat(Db* db1,
                               Db* db2                 = nullptr,
                               int ivar0               = -1,
                               int jvar0               = -1,
                               const VectorInt& nbgh1  = VectorInt(),
                               const VectorInt& nbgh2  = VectorInt(),
                               const CovCalcMode* mode = nullptr);
  MatrixRectangular evalCovMatOptim(Db* db1,
                                    Db* db2                 = nullptr,
                                    int ivar0               = -1,
                                    int jvar0               = -1,
                                    const VectorInt& nbgh1  = VectorInt(),
                                    const VectorInt& nbgh2  = VectorInt(),
                                    const CovCalcMode* mode = nullptr,
                                    bool cleanOptim         = true);
  MatrixRectangular evalCovMatOptimByRanks(const Db* db1,
                                           const Db* db2,
                                           const VectorVectorInt& sampleRanks1,
                                           int ivar0               = -1,
                                           int jvar0               = -1,
                                           const int iech2         = 0,
                                           const CovCalcMode* mode = nullptr,
                                           bool cleanOptim         = true) const;
  FORWARD_METHOD(getCov, evalCovMatSym,)
  FORWARD_METHOD(getCov, evalCovMatSymOptim,)

  MatrixSquareSymmetric evalCovMatSymOptimByRanks(const Db* db1,
                                                  const VectorVectorInt& sampleRanks1,
                                                  int ivar0               = -1,
                                                  const CovCalcMode* mode = nullptr,
                                                  bool cleanOptim         = true);
  MatrixSquareGeneral eval0Mat(const CovCalcMode* mode = nullptr) const;
  MatrixSparse* evalCovMatSparse(Db* db1,
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
  int getNVar() const;
  int getNDim() const;
  const VectorDouble& getMeans() const;

  // Case of _driftList
  int getNDrift() const;
  int getNDriftEquation() const;
  int getNExtDrift() const;
  int getDriftMaxIRFOrder(void) const;
  void delAllDrifts();

    // Case of _cova
  const CovAnisoList* getCovAnisoList() const;
  CovAnisoList* getCovAnisoListModify() const;
  int getCovaMinIRFOrder() const;
  int getNCov(bool skipNugget = false) const;
  void setActiveFactor(int iclass);
  int  getActiveFactor() const;

protected:               // TODO : pass into private to finish clean
  ACov* _cova;           /* Generic Covariance structure */
  DriftList* _driftList; /* Series of Drift functions */
  CovContext _ctxt;      /* Context */
};
