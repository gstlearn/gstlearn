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

class GSTLEARN_EXPORT Line2D : public AStringable, public ASerializable
{
public:
  Line2D(const VectorDouble& x = VectorDouble(),
         const VectorDouble& y = VectorDouble());
  Line2D(const Line2D &m);
  Line2D& operator=(const Line2D &m);
  virtual ~Line2D();

  /// Interface of AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int dumpToNF(const String& neutralFilename, bool verbose = false) const;
  static Line2D* createFromNF(const String& neutralFilename,
                              bool verbose = false);
  static Line2D* create(const VectorDouble& x = VectorDouble(),
                        const VectorDouble& y = VectorDouble());

  int getNPoints() const { return (int) _x.size(); }
  void init(const VectorDouble& x, const VectorDouble& y);
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
