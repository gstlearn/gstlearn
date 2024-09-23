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

#include "LinearOp/CGParam.hpp"
#include "LinearOp/LogStats.hpp"
#include "Basic/VectorNumT.hpp"

class GSTLEARN_EXPORT ALinearOpOld
{

public:
  ALinearOpOld(const CGParam& params = CGParam());
  ALinearOpOld(const ALinearOpOld &m);
  ALinearOpOld& operator=(const ALinearOpOld &m);
  virtual ~ALinearOpOld();

  virtual void evalInverse(const VectorDouble& inv, VectorDouble& outv) const;
  virtual int getSize() const = 0;

  void evalDirect(const VectorDouble& inv, VectorDouble& outv) const;

  void setX0(const VectorDouble& x0) { _params.setX0(x0); }
  void mustShowStats(bool status) { _logStats.mustShowStats(status); }

  const LogStats& getLogStats() const { return _logStats; }

protected:
  virtual void _evalDirect(const VectorDouble& inv, VectorDouble& outv) const = 0;

private:
  double _prod(const VectorDouble& x, const VectorDouble& y) const;

private:
  CGParam _params;

protected:
  LogStats _logStats;
};

