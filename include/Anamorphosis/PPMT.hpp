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
  PPMT(int nbpoly=30, int ndir=10, int legendre_order=2);
  PPMT(const PPMT &m);
  PPMT& operator= (const PPMT &m);
  virtual ~PPMT();

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ICloneable Interface
  IMPLEMENT_CLONING(PPMT)

  int getNiter() const { return (int) _anams.size(); }

public:
  MatrixRectangular fillLegendre(const VectorDouble& r) const;
  AMatrix* sphering(const AMatrix* X);

  VectorDouble generateDirection(double angle) const;
  double getIndex(const AMatrix *X, const VectorDouble &direction) const;
  VectorDouble optimize(const AMatrix *X) const;
  AMatrix* rotate(const AMatrix *X, double alpha, bool direct = true) const;

  int fit(const AMatrix* X, int niter);
  AMatrix* RawToTransform(const AMatrix* X);

private:
  int _nbpoly;
  int _ndir;
  int _legendreOrder;
  MatrixSquareGeneral _S;
  std::vector<AnamHermite> _anams;
  VectorVectorDouble _directions;

  int _nvar;
};
