/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
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
