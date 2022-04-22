/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

// Enums
#include "Covariances/ECov.hpp"

#include "Basic/IClonable.hpp"

#include "Covariances/ACov.hpp"
#include "Covariances/CovCalcMode.hpp"

class ASpace;
class SpacePoint;
class MatrixSquareSymmetric;
class CovAniso;
class CovContext;
class AStringFormat;

class GSTLEARN_EXPORT ACovAnisoList : public ACov, public IClonable
{
public:
  ACovAnisoList(const ASpace* space = nullptr);
  ACovAnisoList(const ACovAnisoList &r);
  ACovAnisoList& operator= (const ACovAnisoList &r);
  virtual ~ACovAnisoList();

  /// Interface for IClonable
  virtual IClonable* clone() const override = 0;

  /// Interface for ASpaceObject
  virtual bool isConsistent(const ASpace* space) const override;

  /// Interface for ACov
  virtual int    getNVariables() const override;
  virtual double eval0(int ivar,
                       int jvar,
                       const CovCalcMode& mode = CovCalcMode()) const override;
  virtual double eval(int ivar,
                      int jvar,
                      const SpacePoint& p1,
                      const SpacePoint& p2,
                      const CovCalcMode& mode = CovCalcMode()) const override;

  /// Interface for AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ACovAnisoList Interface
  virtual void addCov(const CovAniso* cov);

  void addCovList(const ACovAnisoList* covs);
  // Remove an elementary covariance structure
  void delCov(unsigned int i);
  // Remove all elementary covariance structures
  void delAllCov();
  // Filter a covariance
  void setFiltered(unsigned int i, bool filtered);

  int             getCovNumber() const;
  bool            isFiltered(unsigned int i) const;
  bool            hasRange() const;
  bool            isStationary() const;
  double          getMaximumDistance() const;
  double          getTotalSill(int ivar, int jvar) const;
  MatrixSquareGeneral getTotalSill() const;

  /// TODO : to be removed (encapsulation)
  ////////////////////////////////////////////////
  const CovAniso*    getCova(int icov) const;
  CovAniso*          getCova(int icov); // TODO : beurk :(
  const ECov&        getType(int icov) const;
  String             getCovName(int icov) const;
  double             getParam(unsigned int icov) const;
  const MatrixSquareSymmetric& getSill(unsigned int icov) const;
  double             getSill(unsigned int icov, int ivar, int jvar) const;
  int                getGradParamNumber(unsigned int icov) const;
  void               setSill(unsigned int icov, int ivar, int jvar, double value);
  void               setType(unsigned int icov, const ECov& type);
  CovAniso           extractCova(int icov) const;
  ////////////////////////////////////////////////

  void copyCovContext(const CovContext& ctxt);

protected:
  bool   _isCovarianceIndexValid(unsigned int i) const;
  double _getNormalizationFactor(int ivar,
                                 int jvar,
                                 const CovCalcMode& mode) const;

#ifndef SWIG
protected:
  std::vector<CovAniso*> _covs;     /// Vector of elementary covariances
  VectorBool             _filtered; /// Vector of filtered flags (size is nb. cova)
#endif
};

