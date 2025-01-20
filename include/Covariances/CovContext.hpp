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

#include "Space/ASpaceObject.hpp"
#include "Space/ASpace.hpp"

class Vario;
class Db;

class GSTLEARN_EXPORT CovContext : public ASpaceObject
{
public:
  CovContext(int nvar = 1 , ASpaceSharedPtr space = ASpaceSharedPtr());
  CovContext(int nvar,
             int ndim,
             const VectorDouble& mean = VectorDouble(),
             const VectorDouble& covar0 = VectorDouble());
  CovContext(const Db *db, const ASpaceSharedPtr& space = ASpaceSharedPtr());
  CovContext(const Vario *vario, const ASpaceSharedPtr& space = ASpaceSharedPtr());
  CovContext(const CovContext &r);
  CovContext& operator= (const CovContext &r);
  virtual ~CovContext();

  /// AStringable interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Indicate if I am consistent with the provided space
  virtual bool isConsistent(const ASpace* space) const override;

  static CovContext* create(int nvar, int ndim);

  bool isEqual(const CovContext &r) const;

  int                             getNVar()      const { return _nVar; }
  double                          getField()     const { return _field; }
  const VectorDouble&             getMean()      const { return _mean; }
  const VectorDouble&             getCovar0()    const { return _covar0; }
  double getMean(int ivar) const;
  double getCovar0(int ivar, int jvar) const;

  void setNVar(int nvar)                 { _nVar = nvar; _update(); }
  void setField(double field)            { _field = field; }
  void setMean(const VectorDouble& mean);
  void setMean(const double mean,int ivar=0);
  void setCovar0(const VectorDouble& covar0);
  void setCovar0(int ivar, int jvar, double covar0);

  void copyCovContext(const CovContext& ctxt, bool severe = false);

  const CovContext* createReduce(const VectorInt &validVars) const;

private:
  int           _nVar;         /*! Number of variables */
  double        _field;        /*! Field maximum size */
  VectorDouble  _mean;         /*! Array of Variable Mean */
  VectorDouble  _covar0;       /*! Variance-Covariance matrix (used for covariances) */

private:
  int  _getIndex(int ivar, int jvar) const;
  void _update();
};
