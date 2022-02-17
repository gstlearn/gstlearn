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
#include "Basic/Vector.hpp"
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

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int resetFromDb(const Db* db);
  int resetFromCSV(const String& filename,
                   const CSVformat& csv,
                   int verbose = false,
                   int ncol_max = -1,
                   int nrow_max = -1);

  int dumpToNF(const String& neutralFilename, bool verbose = false) const;
  int dumpToNF2(const String& neutralFilename, bool verbose = false) const;
  static Polygons* create();
  static Polygons* createFromNF(const String& neutralFilename, bool verbose = false);
  static Polygons* createFromNF2(const String& neutralFilename, bool verbose = false);
  static Polygons* createFromCSV(const String& filename,
                                 const CSVformat& csv,
                                 int verbose = false,
                                 int ncol_max = -1,
                                 int nrow_max = -1);
  static Polygons* createFromDb(const Db* db);

  int getPolySetNumber() const { return static_cast<int>(_polysets.size()); }
  void addPolySet(const PolySet& polyset);

  const std::vector<PolySet>& getPolySets() const { return _polysets; }
  const PolySet& getPolySet(int ipol) const { return _polysets[ipol]; }
  const VectorDouble& getX(int ipol) const { return _polysets[ipol].getX(); }
  const VectorDouble& getY(int ipol) const { return _polysets[ipol].getY(); }
  void setX(int ipol, const VectorDouble& x) { _polysets[ipol].setX(x); }
  void setY(int ipol, const VectorDouble& y) { _polysets[ipol].setY(y); }

  void getExtension(double *xmin,
                    double *xmax,
                    double *ymin,
                    double *ymax) const;
  double getSurface() const;

protected:
  virtual int _deserialize(FILE* file, bool verbose = false);
  virtual int _serialize(FILE* file, bool verbose = false) const;

  virtual int _deserialize2(std::istream& is, bool verbose = false) override;
  virtual int _serialize2(std::ostream& os, bool verbose = false) const override;

private:
  PolySet _extractFromTab(int ideb,
                          int ifin,
                          int ncol,
                          const VectorDouble& tab);

private:
  std::vector<PolySet> _polysets;
};
