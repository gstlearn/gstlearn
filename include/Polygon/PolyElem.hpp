/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/PolyLine2D.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/AStringable.hpp"

class GSTLEARN_EXPORT PolyElem: public PolyLine2D
{
public:
  PolyElem(const VectorDouble &x = VectorDouble(),
           const VectorDouble &y = VectorDouble(),
           double zmin = TEST,
           double zmax = TEST);
  PolyElem(const PolyElem& r);
  PolyElem& operator=(const PolyElem& r);
  virtual ~PolyElem();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static PolyElem* create();
  static PolyElem* createFromNF(const String& neutralFilename, bool verbose = true);

  const VectorDouble& getX() const { return PolyLine2D::getX(); }
  const VectorDouble& getY() const { return PolyLine2D::getY(); }
  double getX(int i) const { return PolyLine2D::getX(i); }
  double getY(int i) const { return PolyLine2D::getY(i); }
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
  void closePolyElem();
  bool inside(const VectorDouble& coor);
  bool inside3D(double zz);

  PolyElem reduceComplexity(double distmin) const;

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "PolyElem"; }

private:
  bool _isClosed() const;

private:
  double _zmin;
  double _zmax;

  friend class Polygons; // TODO: to be improved (make serialize public)
};
