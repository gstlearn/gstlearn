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

class GSTLEARN_EXPORT PolySet: public AStringable
{
public:
  PolySet();
  PolySet(const VectorDouble& x,
          const VectorDouble& y,
          double zmin = TEST,
          double zmax = TEST);
  PolySet(const PolySet& r);
  PolySet& operator=(const PolySet& r);
  virtual ~PolySet();

  virtual String toString(int level = 0) const override;

  int getNVertices() const { return static_cast<int>(_x.size()); }
  const VectorDouble& getX() const { return _x; }
  const VectorDouble& getY() const { return _y; }
  double getX(int i) const { return _x[i]; }
  double getY(int i) const { return _y[i]; }
  double getZmax() const { return _zmax; }
  double getZmin() const { return _zmin; }

  void init(const VectorDouble& x,
            const VectorDouble& y,
            double zmin = TEST,
            double zmax = TEST);
  void getExtension(double *xmin,
                    double *xmax,
                    double *ymin,
                    double *ymax) const;
  double getSurface() const;

  void setX(const VectorDouble& x) { _x = x; }
  void setY(const VectorDouble& y) { _y = y; }

private:
  VectorDouble _x;
  VectorDouble _y;
  double _zmin;
  double _zmax;
};
