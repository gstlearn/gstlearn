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

#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"

class GSTLEARN_EXPORT SimuBooleanParam: public AStringable
{
public:
  SimuBooleanParam(int maxiter = 100000,
                   double tmax = 100.,
                   double background = TEST,
                   double facies = 1.,
                   const VectorDouble& dilate = VectorDouble());
  SimuBooleanParam(const SimuBooleanParam &r);
  SimuBooleanParam& operator=(const SimuBooleanParam &r);
  virtual ~SimuBooleanParam();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  double getBackground() const { return _background; }
  void setBackground(double background) { _background = background; }
  double getFacies() const { return _facies; }
  void setFacies(double facies) { _facies = facies; }
  int getMaxiter() const { return _maxiter; }
  void setMaxiter(int maxiter) { _maxiter = maxiter; }
  double getTmax() const { return _tmax; }
  void setTmax(double tmax) { _tmax = tmax; }
  const VectorDouble& getDilate() const { return _dilate; }
  void setDilate(const VectorDouble& dilate) { _dilate = dilate; }
  double getDilate(int idim) const;
  bool isDilate() const { return ! _dilate.empty(); }

private:
  int    _maxiter;
  double _tmax;
  double _background;
  double _facies;
  VectorDouble _dilate;
};
