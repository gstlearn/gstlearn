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
#include "Basic/ASerializable.hpp"
#include "Basic/Vector.hpp"

class GSTLEARN_EXPORT PolyLine2D : public AStringable, public ASerializable
{
public:
  PolyLine2D(const VectorDouble& x = VectorDouble(),
             const VectorDouble& y = VectorDouble());
  PolyLine2D(const PolyLine2D &m);
  PolyLine2D& operator=(const PolyLine2D &m);
  virtual ~PolyLine2D();

  /// Interface of AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static PolyLine2D* createFromNF(const String& neutralFilename,
                              bool verbose = false);
  static PolyLine2D* create(const VectorDouble& x = VectorDouble(),
                        const VectorDouble& y = VectorDouble());

  int getNPoints() const { return (int) _x.size(); }
  void init(const VectorDouble& x, const VectorDouble& y);
  const VectorDouble& getX() const { return _x; }
  const VectorDouble& getY() const { return _y; }
  double getX(int i) const { return _x[i]; }
  double getY(int i) const { return _y[i]; }
  double getXmin() const { return ut_vector_min(_x); }
  double getYmin() const { return ut_vector_min(_y); }
  double getXmax() const { return ut_vector_max(_x); }
  double getYmax() const { return ut_vector_max(_y); }
  void addPoint(double x, double y);

  void setX(const VectorDouble& x) { _x = x; }
  void setY(const VectorDouble& y) { _y = y; }

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "PolyLine2D"; }

private:
  VectorDouble _x;
  VectorDouble _y;
};
