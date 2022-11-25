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
#include "Basic/VectorNumT.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/NamingConvention.hpp"

class Db;
class MatrixSquareSymmetric;
class MatrixRectangular;

class GSTLEARN_EXPORT PPMT : public ICloneable, public AStringable
{
public:
  PPMT();
  PPMT(const PPMT &m);
  PPMT& operator= (const PPMT &m);
  virtual ~PPMT();

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ICloneable Interface
  IMPLEMENT_CLONING(PPMT)

  void addIteration(const AnamHermite& anam, const VectorDouble& dir);

public:
  MatrixRectangular fillLegendre(const VectorDouble& r, int n);
  MatrixSquareSymmetric sphering(const MatrixRectangular& X);

private:
  std::vector<AnamHermite> _anams;
  VectorVectorDouble _directions;
};
