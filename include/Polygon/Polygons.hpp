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

#include "Polygon/PolyElem.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/CSVformat.hpp"

class Db;

class GSTLEARN_EXPORT Polygons: public AStringable, public ASerializable
{
public:
  Polygons();
  Polygons(const Polygons& r);
  Polygons& operator=(const Polygons& r);
  virtual ~Polygons();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int resetFromDb(const Db* db, double dilate=0., bool verbose=false);
  int resetFromCSV(const String& filename,
                   const CSVformat& csv,
                   int verbose = false,
                   int ncol_max = -1,
                   int nrow_max = -1);
  int resetFromWKT(const String& filename,
                   const CSVformat& csv,
                   int verbose = false,
                   int ncol_max = -1,
                   int nrow_max = -1);

  static Polygons* create();
  static Polygons* createFromNF(const String& neutralFilename, bool verbose = false);
  static Polygons* createFromCSV(const String& filename,
                                 const CSVformat& csv = CSVformat(),
                                 int verbose = false,
                                 int ncol_max = -1,
                                 int nrow_max = -1);
  static Polygons* createFromWKT(const String& filename,
                                 const CSVformat& csv = CSVformat(),
                                 int verbose = false,
                                 int ncol_max = -1,
                                 int nrow_max = -1);
  static Polygons* createFromDb(const Db* db, double dilate=0., bool verbose=false);

  int getPolyElemNumber() const { return static_cast<int>(_polyelems.size()); }
  void addPolyElem(const PolyElem& polyelem);

  const std::vector<PolyElem>& getPolyElems() const { return _polyelems; }
  const PolyElem& getPolyElem(int ipol) const;
  PolyElem getClosedPolyElem(int ipol) const;
  const VectorDouble& getX(int ipol) const;
  const VectorDouble& getY(int ipol) const;
  void setX(int ipol, const VectorDouble& x);
  void setY(int ipol, const VectorDouble& y);

  void getExtension(double *xmin,
                    double *xmax,
                    double *ymin,
                    double *ymax) const;
  double getSurface() const;
  bool inside(const VectorDouble& coor, bool flag_nested = false);

  Polygons reduceComplexity(double distmin) const;

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "Polygon"; }

private:
  PolyElem _extractFromTab(int ideb,
                          int ifin,
                          int ncol,
                          const VectorDouble &tab);
  PolyElem _extractFromWKT(const CSVformat& csv, String& polye);
  bool _isValidPolyElemIndex(int ipol) const;
  VectorInt _getHullIndices(const VectorDouble& x, const VectorDouble& y) const;
  void _getExtend(double ext, VectorDouble &x, VectorDouble &y, int nsect = 16);
  int  _buildHull(const Db *db, double dilate, bool verbose);
  void _polygonHullPrint(const VectorInt &index,
                         const VectorDouble &x,
                         const VectorDouble &y) const;

private:
  std::vector<PolyElem> _polyelems;

  VectorDouble _emptyVec; // dummy
  PolyElem     _emptyElem; // dummy
};

GSTLEARN_EXPORT void db_polygon(Db *db,
                                Polygons *polygon,
                                bool flag_sel = false,
                                bool flag_period = false,
                                bool flag_nested = false,
                                const NamingConvention& namconv = NamingConvention("Polygon", true, true, true,
                                                                                   ELoc::fromKey("SEL")));
GSTLEARN_EXPORT int dbPolygonDistance(Db *db,
                                      Polygons *polygon,
                                      double dmax,
                                      int scale,
                                      int polin,
                                      const NamingConvention &namconv = NamingConvention("Distance"));
GSTLEARN_EXPORT int db_selhull(Db *db1,
                               Db *db2,
                               double dilate = 0.,
                               bool verbose = false,
                               const NamingConvention& namconv = NamingConvention("Hull", true, true, true,
                                                                                  ELoc::fromKey("SEL")));
