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
  PPMT(int ndir = 50,
       int niter = 10,
       double alpha = 2.,
       const String &method = "vdc");
  PPMT(const PPMT &m);
  PPMT& operator= (const PPMT &m);
  virtual ~PPMT();

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ICloneable Interface
  IMPLEMENT_CLONING(PPMT)

  static PPMT* create(int ndir = 50,
                      int niter = 10,
                      double alpha = 2.,
                      const String &method = "vdc");

  int rawToGaussian(Db *db,
                    const VectorString &names,
                    bool useSel = false,
                    bool verbose = false,
                    const NamingConvention &namconv = NamingConvention("Y"));

  int getNiter()    const { return _niter; }
  double getAlpha() const { return _alpha; }
  int getNdir()     const { return _ndir;  }
  int getNdim()     const { return _ndim;  }
  const String& getMethod() const { return _method; }

  VectorDouble getSerieAngle() const { return _serieAngle; }
  VectorDouble getSerieScore() const { return _serieScore; };

  void fit(AMatrix *X, bool verbose = false);

private:
  MatrixRectangular _fillLegendre(const VectorDouble& r, int legendreOrder) const;
  AMatrix* _sphering(const AMatrix* X);
  void _iteration(AMatrix *Y, const AMatrix *dir, int iter = 0);

private:
  int _niter;
  int _ndir;
  double _alpha;
  String _method;

  mutable int _ndim;
  mutable VectorDouble _serieAngle;
  mutable VectorDouble _serieScore;
  mutable VectorVectorDouble _directions;
};

