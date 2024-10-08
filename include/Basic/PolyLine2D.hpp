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
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/NamingConvention.hpp"

class Db;
class DbGrid;

typedef struct
{
  int rank;
  double dist;
  VectorDouble coor;
} PolyPoint2D;

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

  static PolyLine2D* createFromNF(const String &neutralFilename,
                                  bool verbose = true);
  static PolyLine2D* create(const VectorDouble& x = VectorDouble(),
                            const VectorDouble& y = VectorDouble());

  int getNPoints() const { return (int) _x.size(); }
  void init(const VectorDouble& x, const VectorDouble& y);
  const VectorDouble& getX() const { return _x; }
  const VectorDouble& getY() const { return _y; }
  double getX(int i) const { return _x[i]; }
  double getY(int i) const { return _y[i]; }
  VectorDouble getPoint(int i) const;
  double getXmin() const { return VH::minimum(_x); }
  double getYmin() const { return VH::minimum(_y); }
  double getXmax() const { return VH::maximum(_x); }
  double getYmax() const { return VH::maximum(_y); }
  void addPoint(double x, double y);

  void setX(const VectorDouble& x) { _x = x; }
  void setY(const VectorDouble& y) { _y = y; }

  PolyPoint2D getPLIndex(const VectorDouble &xy0) const;
  double distanceBetweenPoints(double ap,
                               double al,
                               const VectorDouble &xy1,
                               const VectorDouble &xy2) const;
  double distanceAlongPolyline(const PolyPoint2D &pldist1,
                               const PolyPoint2D &pldist2) const;
  double angleAtPolyline(const PolyPoint2D &pldist, int nb_neigh = 0) const;
  double distanceAtPolyline(const PolyPoint2D &pldist,
                            const VectorDouble &target,
                            int nb_neigh) const;

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "PolyLine2D"; }

private:
  static void _shiftPoint(const VectorDouble& xy1,
                          const VectorDouble& xy2,
                          double ratio,
                          VectorDouble& xy0);
  void _getInterval(const PolyPoint2D &pldist,
                    int nb_neigh,
                    int *rfrom,
                    int *rto) const;

private:
  VectorDouble _x;
  VectorDouble _y;
};

double distanceBetweenPolylines(const PolyLine2D& poly1,
                                const PolyLine2D& poly2,
                                const PolyPoint2D& pldist1,
                                const PolyPoint2D& pldist2);
GSTLEARN_EXPORT int dbUnfoldPolyline(Db *db,
                                     const PolyLine2D &polyline,
                                     const NamingConvention &namconv = NamingConvention(
                                         "Unfold"));
GSTLEARN_EXPORT int dbFoldPolyline(DbGrid *dbin,
                                   Db *dbout,
                                   const VectorInt &cols,
                                   const PolyLine2D &polyline,
                                   const NamingConvention &namconv = NamingConvention(
                                       "Fold"));
GSTLEARN_EXPORT int dbFromPolylines(Db* db,
                                    const PolyLine2D &top,
                                    const PolyLine2D &bot,
                                    int nb_neigh = 0,
                                    bool flagMask = true,
                                    const NamingConvention &namconv = NamingConvention(
                                        "Lines"));
