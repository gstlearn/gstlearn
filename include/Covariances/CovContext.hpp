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
#include "Basic/Vector.hpp"
#include "Space/ASpaceObject.hpp"

class ASpace;
class Vario;
class Db;

class GSTLEARN_EXPORT CovContext : public ASpaceObject
{
public:
  CovContext(int nvar = 1,
             const ASpace* space = nullptr,
             int irfMaxDegree = 1000,
             double field = 1);
  CovContext(int nvar,
             int ndim,
             int irfMaxDegree = 1000,
             double field = 1,
             double ballRadius = 0.,
             const VectorDouble& mean = VectorDouble(),
             const VectorDouble& covar0 = VectorDouble());
  CovContext(const Db *db,
             int irfMaxDegree = 1000,
             const ASpace* space = nullptr);
  CovContext(const Vario* vario,
             int irfMaxDegree = 1000,
             const ASpace* space = nullptr);
  CovContext(const CovContext &r);
  CovContext& operator= (const CovContext &r);
  virtual ~CovContext();

  /// AStringable interface
  virtual String toString(int level = 0) const override;

  /// Indicate if I am consistent with the provided space
  virtual bool isConsistent(const ASpace* space) const override;

  /// CovCotext equality
  bool isEqual(const CovContext &r) const;

  int                 getNVar()         const { return _nVar; }
  int                 getIrfMaxDegree() const { return _irfMaxDegree; }
  double              getField()        const { return _field; }
  double              getBallRadius()   const { return _ballRadius; }
  const VectorDouble& getMean()         const { return _mean; }
  const VectorDouble& getCovar0()       const { return _covar0; }

  double getMean(int ivar) const;
  double getCovar0(int ivar, int jvar) const;

  void setNVar(int nvar)                 { _nVar = nvar; _update(); }
  void setIrfMaxDegree(int irfMaxDegree) { _irfMaxDegree = irfMaxDegree; }
  void setField(double field)            { _field = field; }
  void setBallRadius(double ballRadius)  { _ballRadius = ballRadius; }

  void setMean(const VectorDouble& mean);
  void setMean(int ivar, const double mean);

  void setCovar0(const VectorDouble& covar0);
  void setCovar0(int ivar, int jvar, double covar0);

private:
  int           _nVar;         /*! Number of variables */
  int           _irfMaxDegree; /*! Current maximum admissible IRF degree */
  double        _field;        /*! Field maximum size */
  double        _ballRadius;   /*! Radius of the Ball for Numerical Gradient calculation */
  VectorDouble  _mean;         /*! Array of Variable Mean */
  VectorDouble  _covar0;       /*! Variance-Covariance matrix (used for covariances) */

private:
  int _getIndex(int ivar, int jvar) const;
  void _update();
};
