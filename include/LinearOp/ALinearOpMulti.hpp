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
#include "Basic/VectorNumT.hpp"
#include "LinearOp/LogStats.hpp"
#ifndef SWIG
  #include <Eigen/Core>
  #include <Eigen/Dense>
#endif

#include <Eigen/src/Core/Matrix.h>

class ALinearOpMulti;

class GSTLEARN_EXPORT ALinearOpMulti {

public:
  ALinearOpMulti(int nitermax = 1000, double eps = EPSILON8);
  ALinearOpMulti(const ALinearOpMulti &m);
  ALinearOpMulti& operator=(const ALinearOpMulti &m);
  virtual ~ALinearOpMulti();

  void initLk(const std::vector<Eigen::VectorXd> &inv, std::vector<Eigen::VectorXd> &outv) const;
  virtual int sizes() const = 0;
  virtual int size(int) const = 0;

  void setNIterMax(int nitermax) { _nIterMax = nitermax; }
  void setNIterRestart(int niterrestart) { _nIterRestart = niterrestart; }
  void setEps(double eps) { _eps = eps; }
  void setPrecond(const ALinearOpMulti* precond, int status);

  const LogStats& getLogStats() const { return _logStats; }

  void prepare() const;

  void setUserInitialValue(bool b) { _userInitialValue = b; }

  #ifndef SWIG
  protected:

    virtual void _evalDirect(const std::vector<Eigen::VectorXd> &inv,
                             std::vector<Eigen::VectorXd> &outv) const = 0;
  public: 
  void evalDirect(const std::vector<Eigen::VectorXd> &inv,
                  std::vector<Eigen::VectorXd> &outv) const;
  virtual void evalInverse(const std::vector<Eigen::VectorXd> &vecin,
                           std::vector<Eigen::VectorXd> &vecout) const;
  #endif

  protected:
  void _updated() const;

private:
  int                       _nIterMax;
  int                       _nIterRestart;
  double                    _eps;
  bool                      _precondStatus;
  bool                      _userInitialValue;
  const ALinearOpMulti*     _precond;

  // Work arrays
  mutable bool                         _initialized;
  mutable std::vector<Eigen::VectorXd> _r;

public:
  mutable std::vector<Eigen::VectorXd> _temp;
  mutable std::vector<Eigen::VectorXd> _p;
  mutable std::vector<Eigen::VectorXd> _z;
  mutable double _nb;

protected:
  LogStats                   _logStats;
};
