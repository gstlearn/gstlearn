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

#include "LinearOp/ALinearOp.hpp"

class GSTLEARN_EXPORT Identity: public ALinearOp
{

public:
  Identity(int n);
  virtual ~Identity();

  void evalInverse(const VectorDouble& inv, VectorDouble& outv) const override;
  int getSize() const override { return _n; }

protected:
  void _evalDirect(const VectorDouble& inv, VectorDouble& outv) const override;

private:
  int _n;
};
