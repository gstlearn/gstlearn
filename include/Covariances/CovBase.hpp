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
#include "Basic/AFunctional.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/TabNoStat.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Covariances/CovContext.hpp"
#include "Model/CovInternal.hpp"
#include "geoslib_define.h"
#include "Matrix/MatrixT.hpp"
#include "Basic/ParamInfo.hpp"
class AFunctional;
class CovInternal;

class GSTLEARN_EXPORT CovBase: public ACov
{
public:
    
  CovBase(ACov* cor = nullptr,const MatrixSquareSymmetric &sills = MatrixSquareSymmetric());
  CovBase(const CovBase &r) = delete;
  CovBase& operator=(const CovBase &r) = delete;
  virtual ~CovBase();

  ParamInfo createParamInfoForCholSill(int ivar = 0, int jvar = 0);

  virtual bool isConsistent(const ASpace* space) const override;
  virtual int getNVar() const override { return _ctxt.getNVar(); }
  bool isOptimizationInitialized(const Db* db = nullptr) const;
  
  void optimizationSetTargetByIndex(int iech) const override;
  void loadInfoValues() override;
  void setCholSill(int ivar, int jvar, double val) const;
  void setSill(double sill) const; /// Only valid when there is only one variable (in the context)
  void setSill(const MatrixSquareSymmetric& sill) const;
  void setSill(const VectorDouble& sill) const;
  void setSill(int ivar, int jvar, double sill) const;
  void initSill(double value = 0.);

  const MatrixSquareSymmetric& getSill() const { return _sillCur; }
  virtual void setCor(ACov* cor);
  const ACov* getCor() const { return _cor; }
  
  double getSill(int ivar, int jvar) const;
  void   attachNoStatDb(const Db* db);
  
  void   makeSillNoStatDb( const String &namecol, int ivar = 0, int jvar = 0,const Db* db = nullptr);
  void   makeSillStationary(int ivar = 0, int jvar = 0);
  void   makeSillsStationary(bool silent = false);
  void   makeSillNoStatFunctional(  const AFunctional *func, int ivar = 0, int jvar = 0);

  void   makeStationary() override;

  int getNSills()  const {return _tabNoStat.getNSills();}


  bool isNoStatForVariance()   const { return _tabNoStat.isDefinedForVariance();}

  void informMeshByMesh(const AMesh* amesh) const;
  void informMeshByApex(const AMesh* amesh) const;
  VectorDouble informCoords(const VectorVectorDouble& coords, 
                            const EConsElem& econs,
                            int iv1 = 0, int iv2 = 0) const;
  void informDbIn(const Db* dbin) const;
  void informDbOut(const Db* dbout) const;
  void informMeshByMeshForSills(const AMesh* amesh) const;
  void informMeshByApexForSills(const AMesh* amesh) const;
  void informDbInForSills(const Db* dbin) const;
  void informDbOutForSills(const Db* dbout) const;

  /// Tell if the use of Optimization is enabled or not

  void updateCovByPoints(int icas1, int iech1, int icas2, int iech2) const override;
  void updateCovByMesh(int imesh,bool aniso = true) const override;

  double getValue(const EConsElem& econs, int iv1, int iv2) const override;
  void nostatUpdate(CovInternal *covint) const;


protected:
    void _makeElemNoStat(const EConsElem &econs, int iv1, int iv2,
                      const AFunctional* func = nullptr, 
                      const Db* db = nullptr,const String& namecol = String());

  void _manage(const Db* db1,const Db* db2) const override;

  bool _checkSill(int ivar = 0, int jvar = 0) const;
  bool _checkDims(int idim, int jdim) const;

  void _setNoStatDbIfNecessary(const Db*& db);
  bool _checkAndManageNoStatDb(const Db*& db, const String& namecol);
  bool   _isVariableValid(int ivar) const;

  /// Update internal parameters consistency with the context
  virtual void _addEvalCovMatBiPointInPlace(MatrixSquareGeneral& mat,
                                            const SpacePoint& p1,
                                            const SpacePoint& p2,
                                            const CovCalcMode* mode = nullptr) const override;
  virtual void _updateFromContext() override;
  virtual void _initFromContext() override;
  void _copyCovContext(const CovContext& ctxt) override;

private:
  void _setContext(const CovContext& ctxt) override;

  void _optimizationPreProcess(const std::vector<SpacePoint>& p) const override;
  void _optimizationPostProcess() const override;

  void _evalOptim(SpacePoint* p1A,
                  SpacePoint* p2A,
                  MatrixSquareGeneral& mat,
                  const CovCalcMode* mode) const;

  void _load(const SpacePoint& p, bool case1) const override;
  void _optimizationSetTarget(const SpacePoint& pt) const override;

protected:
    MatrixT<ParamInfo> _cholSillsInfo;
    mutable MatrixSquareGeneral _cholSills;
    TabNoStat _tabNoStat;
    mutable MatrixSquareSymmetric _sillCur;
    mutable MatrixSquareGeneral _workMat;
private :
    ACov* _cor;
};
