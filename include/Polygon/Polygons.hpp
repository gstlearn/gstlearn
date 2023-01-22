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

#include "Polygon/PolySet.hpp"
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

  static Polygons* create();
  static Polygons* createFromNF(const String& neutralFilename, bool verbose = false);
  static Polygons* createFromCSV(const String& filename,
                                 const CSVformat& csv = CSVformat(),
                                 int verbose = false,
                                 int ncol_max = -1,
                                 int nrow_max = -1);
  static Polygons* createFromDb(const Db* db);

  int getPolySetNumber() const { return static_cast<int>(_polysets.size()); }
  void addPolySet(const PolySet& polyset);

  const std::vector<PolySet>& getPolySets() const { return _polysets; }
  const PolySet getPolySet(int ipol) const;
  PolySet getClosedPolySet(int ipol) const;
  const VectorDouble getX(int ipol) const;
  const VectorDouble getY(int ipol) const;
  void setX(int ipol, const VectorDouble& x);
  void setY(int ipol, const VectorDouble& y);

  void getExtension(double *xmin,
                    double *xmax,
                    double *ymin,
                    double *ymax) const;
  double getSurface() const;

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "Polygon"; }

private:
  PolySet _extractFromTab(int ideb,
                          int ifin,
                          int ncol,
                          const VectorDouble& tab);
  bool _isValidPolySetIndex(int ipol) const;

private:
  std::vector<PolySet> _polysets;
};
