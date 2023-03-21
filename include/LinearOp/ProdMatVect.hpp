/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "LinearOp/ALinearOp.hpp"

class GSTLEARN_EXPORT ProdMatVect : public ALinearOp {

public:
  ProdMatVect(int nx, int ny, double *A);
  virtual ~ProdMatVect();

  int getSize() const override
  {
    return _nx;
  }

protected:
  void _evalDirect(const VectorDouble& inv, VectorDouble& outv) const override;

private :
  int     _nx;
	int     _ny;
	double* _A;
};
