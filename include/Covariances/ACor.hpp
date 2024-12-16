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

#include "Basic/AStringable.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Space/ASpaceObject.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Space/SpacePoint.hpp"

#include <vector>

class Db;
class DbGrid;
class MatrixSquareGeneral;
class MatrixSparse;

/**
 * \brief
 * Class containing the Covariance part of the Model.
 *
 * It is the uppermost class of the Covariance Tree and is conceived as simple as possible on purpose
 * (in order to let the user defined its own version if necessary): it must simply be able to return its value
 * between two end-point (see eval method).
 *
 * It is mainly implemented in CovAniso.hpp or ACorAnisoList.hpp
 */
class GSTLEARN_EXPORT ACor : public ASpaceObject
{
public:
  ACor(const ASpace* space = nullptr, int nvar = 1);
  ACor(const ACor &r);
  ACor& operator=(const ACor &r);
  virtual ~ACor();

  /// ACor Interface
  virtual int getNVariables() const {return _nvar;};

 
  /// Calculate the covariance between two variables for 0-distance (stationary case)
  virtual double eval0(int ivar = 0,
                       int jvar = 0,
                       const CovCalcMode* mode = nullptr) const;
  /// Calculate the matrix of covariances for 0-distance (stationary case)
  
  virtual void eval0CovMatBiPointInPlace(MatrixSquareGeneral &mat,
                                 const CovCalcMode *mode = nullptr) const;
  virtual void addEval0CovMatBiPointInPlace(MatrixSquareGeneral &mat,
                                            const CovCalcMode *mode = nullptr) const;
  /// Calculate the covariance between two variables and two points (general case)
  virtual double eval(const SpacePoint& p1,
                      const SpacePoint& p2,
                      int ivar = 0,
                      int jvar = 0,
                      const CovCalcMode* mode = nullptr) const = 0;
  /// Calculate the matrix of covariances between two points (general case)
  virtual void evalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
                                        const SpacePoint &p1,
                                        const SpacePoint &p2,
                                        const CovCalcMode *mode = nullptr) const; 
                                        
  virtual void addEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
                               const SpacePoint& pwork1, 
                               const SpacePoint& pwork2,
                               const CovCalcMode *mode) const;
                               
  void evalCovKriging(MatrixSquareGeneral &mat,
                      SpacePoint &pwork1,
                      SpacePoint& pout, 
                      const CovCalcMode *mode = nullptr) const;
  virtual double evalCovOnSphere(double alpha,
                                 int degree = 50,
                                 bool flagScaleDistance = false,
                                 const CovCalcMode* mode = nullptr) const
  {
    DECLARE_UNUSED(alpha);
    DECLARE_UNUSED(degree);
    DECLARE_UNUSED(flagScaleDistance);
    DECLARE_UNUSED(mode);
    return TEST;
  }
  virtual VectorDouble evalSpectrumOnSphere(int n,
                                            bool flagNormDistance = false,
                                            bool flagCumul = false) const
  {
    DECLARE_UNUSED(n);
    DECLARE_UNUSED(flagNormDistance);
    DECLARE_UNUSED(flagCumul);
    return VectorDouble();
  }
  virtual double evalSpectrum(const VectorDouble &freq,
                              int ivar,
                              int jvar) const
  {
    DECLARE_UNUSED(freq);
    DECLARE_UNUSED(ivar);
    DECLARE_UNUSED(jvar);
    return TEST;
  }

  virtual void updateCovByPoints(int icas1, int iech1, int icas2, int iech2)
  {
    DECLARE_UNUSED(icas1);
    DECLARE_UNUSED(iech1);
    DECLARE_UNUSED(icas2);
    DECLARE_UNUSED(iech2);
  }

  /////////////////////////////////////////////////////////////////////////////////
  ///

  void optimizationSetTarget(const SpacePoint &pt) const;
  virtual void optimizationSetTargetByIndex(int iech) const {DECLARE_UNUSED(iech)};
  void optimizationPreProcess(const Db* db) const;
  void optimizationPreProcess(const std::vector<SpacePoint>& p) const;

  void optimizationPostProcess() const;
  virtual bool isOptimEnabled() const {return _isOptimEnabled();}

  VectorDouble eval(const std::vector<SpacePoint>& vec_p1,
                    const std::vector<SpacePoint>& vec_p2,
                    int ivar = 0,
                    int jvar = 0,
                    const CovCalcMode* mode = nullptr) const;
  MatrixSquareGeneral eval0Mat(const CovCalcMode* mode = nullptr) const;
  MatrixSquareGeneral evalMat(const SpacePoint& p1,
                              const SpacePoint& p2,
                              const CovCalcMode* mode = nullptr) const;

  

  MatrixRectangular evalCovMatrix(const Db* db1_arg,
                                  const Db* db2_arg = nullptr,
                                  int ivar0 = -1,
                                  int jvar0 = -1,
                                  const VectorInt& nbgh1 = VectorInt(),
                                  const VectorInt& nbgh2 = VectorInt(),
                                  const CovCalcMode* mode = nullptr);
  MatrixSquareSymmetric evalCovMatrixSymmetric(const Db *db1,
                                               int ivar0,
                                               const VectorInt &nbgh1,
                                               const CovCalcMode *mode);
  MatrixSparse* evalCovMatrixSparse(const Db *db1_arg,
                                    const Db *db2_arg = nullptr,
                                    int ivar0 = -1,
                                    int jvar0 = -1,
                                    const VectorInt &nbgh1 = VectorInt(),
                                    const VectorInt &nbgh2 = VectorInt(),
                                    const CovCalcMode *mode = nullptr,
                                    double eps = EPSILON3);
  


  void manage(const Db* db1,const Db* db2) const
  {
      _manage(db1, db2);
  }

  void load(const SpacePoint& p,bool case1) const;

  void loadAndAddEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,const SpacePoint& p1,const SpacePoint&p2,
                                              const CovCalcMode *mode = nullptr) const;

  double loadAndEval(const SpacePoint& p1,
                          const SpacePoint&p2,
                          int ivar,
                          int jvar,
                          const CovCalcMode *mode) const;
protected:
  virtual void _loadAndAddEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,const SpacePoint& p1,const SpacePoint&p2,
                                              const CovCalcMode *mode = nullptr) const;
  virtual void _optimizationSetTarget(const SpacePoint &pt) const;

  void _setOptimEnabled(bool enabled){ _optimEnabled = enabled;}
  VectorInt _getActiveVariables(int ivar0) const;
  static void _updateCovMatrixSymmetricVerr(const Db* db1,
                                            AMatrix* mat,
                                            const VectorInt& ivars,
                                            const VectorVectorInt& index1);

  virtual void _optimizationPreProcess(const std::vector<SpacePoint>& p) const;
  virtual void _addEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
                                            const SpacePoint& pwork1, 
                                            const SpacePoint& pwork2,
                                            const CovCalcMode *mode) const;
double _loadAndEval(const SpacePoint& p1,
                    const SpacePoint&p2,
                    int ivar = 0,
                    int jvar = 0,
                    const CovCalcMode *mode = nullptr) const;
private:
  virtual void _optimizationPostProcess() const; 
  virtual bool _isOptimEnabled() const {return _optimEnabled;}

  virtual void _manage(const Db* db1,const Db* db2) const 
  {
    DECLARE_UNUSED(db1)
    DECLARE_UNUSED(db2)
  }

  

protected:
  int _nvar;
  bool _optimEnabled;
  mutable bool _isOptimPreProcessed;
  mutable std::vector<SpacePoint> _p1As;
  mutable SpacePoint _p2A;
  const mutable SpacePoint* _pw1;
  const mutable SpacePoint* _pw2;
};
