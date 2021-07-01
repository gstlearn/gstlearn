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

#include "Db/Db.hpp"
#include "Basic/Vector.hpp"
#include "Space/ASpace.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"

class CovContext : public AStringable
{
private:
  int           _nVar;         /*! Number of variables */
  int           _irfMaxDegree; /*! Current maximum admissible IRF degree */
  double        _field;        /*! Field maximum size */
  double        _ballRadius;   /*! Radius of the Ball for Numerical Gradient calculation */
  VectorDouble  _mean;         /*! Array of Variable Mean */
  VectorDouble  _covar0;       /*! Variance-Covariance matrix (used for covariances) */

  /// TODO : use shared pointer for ASpace* ?
  const ASpace* _space;        /*! Space context (Number of dimension, getDistance, etc...) */

public:
  /// TODO : default context (1 variable, big max IRF degree, and field size of 1) ok ?
  CovContext(int nvar = 1,
             int irfMaxDegree = 1000,
             double field = 1,
             const ASpace* space = nullptr);
  CovContext(const Db *db, int irfMaxDegree = 1000, const ASpace* space = nullptr);
  CovContext(const CovContext &r);
  CovContext& operator= (const CovContext &r);
  virtual ~CovContext();

  virtual std::string toString(int level = 0) const override;

  bool isEqual(const CovContext& r) const;

  int           getNVar()         const { return _nVar; }
  int           getIrfMaxDegree() const { return _irfMaxDegree; }
  double        getField()        const { return _field; }
  const ASpace* getSpace()        const { return _space; }
  unsigned int  getNDim()         const { return _space->getNDim(); }

  void setField(double field)
  {
    _field = field;
  }

  void setIrfMaxDegree(int irfMaxDegree)
  {
    _irfMaxDegree = irfMaxDegree;
  }

  void setNVar(int nvar)
  {
    _nVar = nvar;
    _update();
  }

  void setSpace(const ASpace* space)
  {
    _space = space;
  }

  const VectorDouble& getMean() const
  {
    return _mean;
  }
  const double getMean(int ivar) const
  {
    if (ivar < 0 || ivar >= (int) _mean.size())
      throw("Invalid argument in _getMean");
    return _mean[ivar];
  }

  void setMean(const VectorDouble& mean)
  {
    if (_mean.size() == mean.size())
      _mean = mean;
  }

  void setMean(int ivar, const double mean)
  {
    if (ivar < 0 || ivar >= (int) _mean.size())
      throw("Invalid argument in _setMean");
    _mean[ivar] = mean;
  }

  void setCovar0(int ivar, int jvar, double covar0)
  {
    int rank = _getIndex(ivar, jvar);
    if (rank < 0 || rank >= (int) _covar0.size())
      throw("Invalid argument in _setCovar0");
    _covar0[rank] = covar0;
  }
  void setCovar0(const VectorDouble& covar0)
  {
    if (_covar0.size() == covar0.size())
      _covar0 = covar0;
  }

  const VectorDouble& getCovar0() const
  {
    return _covar0;
  }
  const double getCovar0(int ivar, int jvar) const
  {
    int rank = _getIndex(ivar, jvar);
    if (rank < 0 || rank >= (int) _covar0.size())
      throw("Invalid argument in _setCovar0");
    return _covar0[rank];
  }

  double getBallRadius() const
  {
    return _ballRadius;
  }

  void setBallRadius(double ballRadius)
  {
    _ballRadius = ballRadius;
  }

private:
  int _getIndex(int ivar, int jvar) const
  {
    return ivar * getNVar() + jvar;
  }
  void _update();
};
