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

#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/CSVformat.hpp"
#include "Basic/VectorNumT.hpp"
#include "Polygon/PolyElem.hpp"

namespace gstlrn
{
class Db;

class GSTLEARN_EXPORT Polygons: public AStringable, public ASerializable
{
public:
  Polygons();
  Polygons(const Polygons& r);
  Polygons& operator=(const Polygons& r);
  virtual ~Polygons();

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  Id resetFromDb(const Db* db, double dilate = 0., bool verbose = false);
  Id resetFromCSV(const String& filename,
                   const CSVformat& csv,
                   bool verbose = false,
                   Id ncol_max = -1,
                   Id nrow_max = -1);
  Id resetFromWKT(const String& filename,
                   const CSVformat& csv,
                   bool verbose = false,
                   Id ncol_max = -1,
                   Id nrow_max = -1);

  static Polygons* create();
  static Polygons* createFromNF(const String& NFFilename, bool verbose = false);
  static Polygons* createFromCSV(const String& filename,
                                 const CSVformat& csv = CSVformat(),
                                 Id verbose          = false,
                                 Id ncol_max         = -1,
                                 Id nrow_max         = -1);
  static Polygons* createFromWKT(const String& filename,
                                 const CSVformat& csv = CSVformat(),
                                 Id verbose          = false,
                                 Id ncol_max         = -1,
                                 Id nrow_max         = -1);
  static Polygons* createFromDb(const Db* db, double dilate = 0., bool verbose = false);

  Id getNPolyElem() const { return static_cast<Id>(_polyelems.size()); }
  void addPolyElem(const PolyElem& polyelem);

  const std::vector<PolyElem>& getPolyElems() const { return _polyelems; }
  const PolyElem& getPolyElem(Id ipol) const;
  PolyElem getClosedPolyElem(Id ipol) const;
  const VectorDouble& getX(Id ipol) const;
  const VectorDouble& getY(Id ipol) const;
  void setX(Id ipol, const VectorDouble& x);
  void setY(Id ipol, const VectorDouble& y);

  void getExtension(double* xmin,
                    double* xmax,
                    double* ymin,
                    double* ymax) const;
  double getSurface() const;
  bool inside(const VectorDouble& coor, bool flag_nested = false) const;

  Polygons reduceComplexity(double distmin) const;

protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "Polygon"; }

private:
  static PolyElem _extractFromTab(Id ideb, Id ifin, Id ncol, const VectorDouble& tab);
  static PolyElem _extractFromWKT(const CSVformat& csv, String& polye);
  bool _isValidPolyElemIndex(Id ipol) const;
  static VectorInt _getHullIndices(const VectorDouble& x, const VectorDouble& y);
  static void _getExtend(double ext, VectorDouble& x, VectorDouble& y, Id nsect = 16);
  Id _buildHull(const Db* db, double dilate, bool verbose);
  static void _polygonHullPrint(const VectorInt& index,
                                const VectorDouble& x,
                                const VectorDouble& y);

private:
  std::vector<PolyElem> _polyelems;

  VectorDouble _emptyVec; // dummy
  PolyElem _emptyElem;    // dummy
};

GSTLEARN_EXPORT void db_polygon(Db* db,
                                const Polygons* polygon,
                                bool flag_sel                   = false,
                                bool flag_period                = false,
                                bool flag_nested                = false,
                                const NamingConvention& namconv = NamingConvention("Polygon", true, true, true, ELoc::fromKey("SEL")));
GSTLEARN_EXPORT Id dbPolygonDistance(Db* db,
                                      Polygons* polygon,
                                      double dmax,
                                      Id scale,
                                      Id polin,
                                      const NamingConvention& namconv = NamingConvention("Distance"));
GSTLEARN_EXPORT Id db_selhull(Db* db1,
                               Db* db2,
                               double dilate                   = 0.,
                               bool verbose                    = false,
                               const NamingConvention& namconv = NamingConvention("Hull", true, true, true, ELoc::fromKey("SEL")));
} // namespace gstlrn