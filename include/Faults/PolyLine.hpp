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

class GSTLEARN_EXPORT PolyLine: public AStringable, public ASerializable
{
public:
  PolyLine();
  PolyLine(const VectorDouble& x, const VectorDouble& y);
  PolyLine(const PolyLine& r);
  PolyLine& operator=(const PolyLine& r);
  virtual ~PolyLine();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int dumpToNF(const String& neutralFilename, bool verbose = false) const;
  static PolyLine* create(const VectorDouble& x = VectorDouble(),
                          const VectorDouble& y = VectorDouble());
  static PolyLine* createFromNF(const String& neutralFilename, bool verbose = false);

  void init(const VectorDouble& x, const VectorDouble& y);
  int getNVertices() const { return _x.size(); }

  const VectorDouble& getX() const { return _x; }
  const VectorDouble& getY() const { return _y; }
  double getX(int i) const { return _x[i]; }
  double getY(int i) const { return _y[i]; }

protected:
  virtual int _deserialize(std::istream& is, bool verbose = false) override;
  virtual int _serialize(std::ostream& os, bool verbose = false) const override;

private:
  VectorDouble _x;
  VectorDouble _y;
};
