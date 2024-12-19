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



class MatrixSquareGeneral;
class MatrixRectangular;
class CovInternal;
class ACov;
class ACor;
/**
 * \brief
 * This class describes an **elementary covariance**.
 *
 * This covariance is described through the following list of parameters:
 * - the covariance **type**: the list of these types is provided in ECov.hpp
 * - the largest set of parameters for any covariance: **range(s)**, **anisotropy angle(s)**, **third parameter**. Some of these parameters
 * do not make sense, depending on the covariance type: e.g. the range for nugget effect, the third parameter for a spherical
 * structure, ...
 * All these parameters are processed and stored as a **tensor** in order to avoid repetitive calculations.
 * - the **sill**. This comes as a square symmetric matrix whose dimension is equal to the number of variables.
 */
class GSTLEARN_EXPORT CovBase: public ACov, public ICloneable
{
public:
  CovBase(const ACor* = nullptr, MatrixSquareSymmetric sills);
  CovBase(const CovBase& r);
  CovBase& operator=(const CovBase& r);
  virtual ~CovBase();

  /// ICloneable Interface
  IMPLEMENT_CLONING(CovBase)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ACov Interface
  virtual double eval0(int ivar = 0,
                       int jvar = 0,
                       const CovCalcMode* mode = nullptr) const override;
  virtual double eval(const SpacePoint& p1,
                      const SpacePoint& p2,
                      int ivar = 0,
                      int jvar = 0,
                      const CovCalcMode* mode = nullptr) const override;
 
  virtual void addEval0CovMatBiPointInPlace(MatrixSquareGeneral &mat,
                               const CovCalcMode *mode = nullptr) const override;
 
  virtual double evalCovOnSphere(double alpha,
                                 int degree = 50,
                                 bool flagScaleDistance = true,
                                 const CovCalcMode* mode = nullptr) const override;
  virtual VectorDouble evalSpectrumOnSphere(int n,
                                            bool flagNormDistance = false,
                                            bool flagCumul = false) const override;
  virtual double evalSpectrum(const VectorDouble &freq,
                              int ivar = 0,
                              int jvar = 0) const override;



  void setSill(double sill); /// Only valid when there is only one variable (in the context)
  void setSill(const MatrixSquareSymmetric& sill);
  void setSill(const VectorDouble& sill);
  void setSill(int ivar, int jvar, double sill);
  void initSill(double value = 0.);

  const MatrixSquareSymmetric& getSill() const { return _sill; }
  double getSill(int ivar, int jvar) const;
  
  //////////////////////// New NoStat methods //////////////////////////
  void   attachNoStatDb(const Db* db);
 
  void   makeSillNoStatDb(  const String &namecol, int ivar = 0, int jvar = 0,const Db* db = nullptr);
  void   makeSillNoStatFunctional(  const AFunctional *func, int ivar = 0, int jvar = 0);
  void   makeSillStationary( int ivar = 0, int jvar = 0);
  void   makeStationary();



  int getNSills()  const {return _tabNoStat.getNSills();}
  bool isNoStat() const  override
  { 
    return _tabNoStat.isNoStat() || _cor.isNoStat(); 
  };
 
  bool isNoStatForVariance()   const { return _tabNoStat.isDefinedForVariance();}
 
  
  VectorDouble evalCovOnSphereVec(const VectorDouble &alpha,
                                  int degree = 50,
                                  bool flagScaleDistance = false,
                                  const CovCalcMode* mode = nullptr) const;
  Array evalCovFFT(const VectorDouble& hmax, int N = 128, int ivar = 0, int jvar = 0) const;
  
  //CovBase* createReduce(const VectorInt &validVars) const;

 
  void informDbInForSills(const Db* dbin) const;
  void informDbOutForSills(const Db* dbout) const;



protected:
  /// Update internal parameters consistency with the context
  virtual void _addEvalCovMatBiPointInPlace(MatrixSquareGeneral& mat,
                                            const SpacePoint& p1,
                                            const SpacePoint& p2,
                                            const CovCalcMode* mode = nullptr) const override;
  virtual void _updateFromContext();
  virtual void _initFromContext();
  void _optimizationSetTarget(const SpacePoint& pt) const override;

private:
void _optimizationPostProcess() const override; 

bool _isOptimEnabled() const override 
{ 
  return _optimEnabled && !isNoStatForAnisotropy(); 
}
void  _evalOptim(SpacePoint* p1A, SpacePoint* p2A,
                 MatrixSquareGeneral &mat,
                 const CovCalcMode *mode) const;
 void _makeElemNoStat(const EConsElem &econs, int iv1, int iv2,
                      const AFunctional* func = nullptr, 
                      const Db* db = nullptr,const String& namecol = String());

  void _manage(const Db* db1,const Db* db2) const override;

  bool _checkSill(int ivar = 0, int jvar = 0) const;
  bool _checkDims(int idim, int jdim) const;

  void _setNoStatDbIfNecessary(const Db*& db);
  bool _checkAndManageNoStatDb(const Db*& db, const String& namecol);
  bool   _isVariableValid(int ivar) const;
  void   _computeCorrec();
  double _getDetTensor() const;
  void   _optimizationTransformSP(const SpacePoint& ptin, SpacePoint& ptout) const;

private:

  ACor _cor;
  mutable MatrixSquareSymmetric _sill; /// Sill matrix (nvar x nvar)
  TabNoStat _tabNoStat;
  // These temporary information is used to speed up processing (optimization functions)
  // They are in a protected section as they may be modified by class hierarchy
};



  