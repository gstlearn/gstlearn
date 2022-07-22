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
#include "Matrix/MatrixRectangular.hpp"
#include "Mesh/AMesh.hpp"
#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/AFunctional.hpp"
#include "Model/ANoStat.hpp"

/**
 * This class concerns the non-stationarity defined as a function (hence its name).
 * It can be considered as an example of a 2-D implementation of a spirale with
 * a single non-stationary parameter (the angle)
 */
class GSTLEARN_EXPORT NoStatFunctional : public ANoStat
{
public:
	NoStatFunctional();
	NoStatFunctional(const AFunctional* func);
  NoStatFunctional(const NoStatFunctional &m);
  NoStatFunctional& operator=(const NoStatFunctional &m);
  virtual ~NoStatFunctional();

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  virtual ICloneable* clone() const override { return new NoStatFunctional(*this); };

  int  attachToMesh(const AMesh* mesh, bool verbose = false) const override;
  int  attachToDb(Db* db, int icas, bool verbose = false) const override;

  double getValue(int igrf, int icov, const EConsElem& type, int iv1, int iv2,
                  int icas, int rank) const override;
  double getValueByParam(int ipar, int icas, int rank) const override;

private:
  const AFunctional* _func;
};
