/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/PolyLine2D.hpp"
#include "Basic/VectorNumT.hpp"

namespace gstlrn
{
  class GSTLEARN_EXPORT PolyElem: public PolyLine2D
  {
  public:
    PolyElem(
      const VectorDouble& x = VectorDouble(),
      const VectorDouble& y = VectorDouble(),
      double zmin = TEST,
      double zmax = TEST);
    PolyElem(const PolyElem& r);
    PolyElem& operator=(const PolyElem& r);
    virtual ~PolyElem();

    /// ASerializable Interface
    String getNFName() const override { return "PolyElem"; }
#ifdef HDF5
    bool deserializeH5(H5::Group& grp) override;
    bool serializeH5(H5::Group& grp) const override;
#endif

    /// AStringable Interface
    String toString(const AStringFormat* strfmt = nullptr) const override;

    static PolyElem* create();
    static PolyElem*
      createFromNF(const String& NFFilename, bool verbose = true);

    const VectorDouble& getX() const { return PolyLine2D::getX(); }

    const VectorDouble& getY() const { return PolyLine2D::getY(); }

    double getX(Id i) const { return PolyLine2D::getX(i); }

    double getY(Id i) const { return PolyLine2D::getY(i); }

    double getZmax() const { return _zmax; }

    double getZmin() const { return _zmin; }

    void init(
      const VectorDouble& x,
      const VectorDouble& y,
      double zmin = TEST,
      double zmax = TEST);
#ifndef SWIG
    void getExtension(double& xmin, double& xmax, double& ymin, double& ymax)
      const;
#endif
    double getSurface() const;
    void closePolyElem();
    bool inside(const VectorDouble& coor);
    bool inside3D(double zz) const;
    VectorDouble getCentroid() const;

    PolyElem reduceComplexity(double distmin) const;

  protected:
    bool _deserializeAscii(std::istream& is) override;
    bool _serializeAscii(std::ostream& os) const override;

  private:
    bool _isClosed() const;

  private:
    double _zmin;
    double _zmax;

    friend class Polygons;
  };
} // namespace gstlrn
