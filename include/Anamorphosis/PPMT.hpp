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

#include "Anamorphosis/AnamHermite.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/NamingConvention.hpp"

class Db;
class MatrixSquareSymmetric;
class MatrixRectangular;
class AMatrix;

class GSTLEARN_EXPORT PPMT : public ICloneable, public AStringable
{
public:
  PPMT(int nvar = 0);
  PPMT(const PPMT &m);
  PPMT& operator= (const PPMT &m);
  virtual ~PPMT();

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ICloneable Interface
  IMPLEMENT_CLONING(PPMT)

  int getNiter() const { return (int) _anams.size(); }

public:
  MatrixRectangular fillLegendre(const VectorDouble& r, int n) const;
  MatrixRectangular sphering(const MatrixRectangular& X);

  VectorDouble generateDirection(double angle) const;
  double getIndex(const MatrixRectangular &X,
                  const VectorDouble &direction,
                  int j = 2) const;
  VectorDouble optimize(const MatrixRectangular &X,
                        int j=2,
                        int N=10) const;
  MatrixRectangular rotate(const MatrixRectangular &X,
                           double alpha,
                           bool direct = true) const;
  void fit(const MatrixRectangular &X,
           int niter,
           int j = 2,
           int N = 10,
           int nbpoly = 20);
  MatrixRectangular eval(const MatrixRectangular& X);

private:
  int _nvar;
  MatrixSquareGeneral _S;
  std::vector<AnamHermite> _anams;
  VectorVectorDouble _directions;
};
