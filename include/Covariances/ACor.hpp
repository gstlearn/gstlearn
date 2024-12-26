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
#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"
#include "Covariances/CovContext.hpp"
#include "Enum/EConsElem.hpp"
#include "Covariances/TabNoStat.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Mesh/AMesh.hpp"
#include "geoslib_define.h"
#include "Covariances/ACov.hpp"
#include "Space/ASpaceObject.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Space/SpacePoint.hpp"

#include <vector>

class CovContext;
class Db;
class DbGrid;
class MatrixSquareGeneral;
class MatrixSparse;
class AFunctional;
class EConsElem;
class AMesh;

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
  ACor(const CovContext& ctxt);
  ACor(const ACor &r);
  ACor& operator=(const ACor &r);
  virtual ~ACor();

  /// ACor Interface
  virtual int getNVariables() const {return _nvar;};

  void setContext(const CovContext& ctxt) { _ctxt = ctxt; }

  /// Calculate the covariance between two variables for 0-distance (stationary case)
  virtual double eval0(int ivar = 0,
                       int jvar = 0,
                       const CovCalcMode* mode = nullptr) const;
 
  virtual double eval(const SpacePoint& p1,
                      const SpacePoint& p2,
                      int ivar = 0,
                      int jvar = 0,
                      const CovCalcMode* mode = nullptr) const = 0;
virtual void copyCovContext(const CovContext& ctxt) ;
virtual void updateFromContext() {}

virtual void optimizationPostProcess() const {}
virtual int makeElemNoStat(const EConsElem &econs, int iv1, int iv2,
                     const AFunctional* func = nullptr, 
                     const Db* db = nullptr,const String& namecol = String());

  const CovContext& getContext() const { return _ctxt; }

virtual void initFromContext()
{
}

virtual void optimizationPreProcess(const std::vector<SpacePoint>& p,
                               std::vector<SpacePoint> &p1As) const
{
  DECLARE_UNUSED(p,p1As)
}

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
  bool _checkDims(int idim, int jdim) const;
  
  virtual VectorDouble evalSpectrumOnSphere(int n,
                                            bool flagNormDistance = false,
                                            bool flagCumul = false) const
  {
    DECLARE_UNUSED(n);
    DECLARE_UNUSED(flagNormDistance);
    DECLARE_UNUSED(flagCumul);
    return VectorDouble();
  }
  void   attachNoStatDb(const Db* db);
  virtual double evalSpectrum(const VectorDouble &freq,
                              int ivar,
                              int jvar) const
  {
    DECLARE_UNUSED(freq);
    DECLARE_UNUSED(ivar);
    DECLARE_UNUSED(jvar);
    return TEST;
  }

  bool checkAndManageNoStatDb(const Db*& db, const String& namecol);
  virtual void updateCovByMesh(int imesh,bool aniso = true) 
  {
    DECLARE_UNUSED(imesh,aniso)
  }
  virtual double getValue(const EConsElem &econs,int iv1,int iv2) const
  {
    DECLARE_UNUSED(econs,iv1,iv2)
    return TEST;
  }
  virtual void makeStationary();

  virtual void createNoStatTab();
  void informMeshByMesh(const AMesh* amesh) const;
  void informMeshByApex(const AMesh* amesh) const;
  VectorDouble informCoords(const VectorVectorDouble& coords, 
                            const EConsElem& econs,
                            int iv1 = 0, int iv2 = 0) const;
  void informDbIn(const Db* dbin) const;
  void informDbOut(const Db* dbout) const;

  virtual void updateCovByPoints(int icas1, int iech1, int icas2, int iech2)
  {
    DECLARE_UNUSED(icas1);
    DECLARE_UNUSED(iech1);
    DECLARE_UNUSED(icas2);
    DECLARE_UNUSED(iech2);
  }

  
  void manage(const Db* db1,const Db* db2) const
  {
      _manage(db1, db2);
  }
  int getDimensionNumber() const { return _ctxt.getNDim(); }

private:

  virtual void _manage(const Db* db1,const Db* db2) const 
  {
    DECLARE_UNUSED(db1)
    DECLARE_UNUSED(db2)
  }


  void setNoStatDbIfNecessary(const Db*& db);
private : 
 virtual void _copyCovContext(const CovContext &ctxt)
 {
  DECLARE_UNUSED(ctxt)
 }
protected:
  int _nvar;
  TabNoStat* _tabNoStat;

private:
  CovContext _ctxt;
};
