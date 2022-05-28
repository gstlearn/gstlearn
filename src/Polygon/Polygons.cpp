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
#include "Polygon/Polygons.hpp"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

#include "Db/Db.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/CSVformat.hpp"
#include "Basic/AException.hpp"

Polygons::Polygons()
  : _polysets()
{
}

Polygons::Polygons(const Polygons& r)
    : AStringable(r),
      ASerializable(r),
      _polysets(r._polysets)
{
}

Polygons& Polygons::operator=(const Polygons& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    ASerializable::operator=(r);
    _polysets = r._polysets;
  }
  return *this;
}

Polygons::~Polygons()
{
}

/**
 * Calculate the Polygon as the convex hull of the active samples of a Db
 * @param db
 */
int Polygons::resetFromDb(const Db* db)
{
  if (db == nullptr) return 1;

  // Clear previous contents
  _polysets.clear();

  // Calculate the hull
  VectorDouble x;
  VectorDouble y;
  if (polygon_hull(db, x, y)) return 1;

  PolySet polyset = PolySet(x, y, TEST, TEST);
  addPolySet(polyset);

  return 0;
}

/**
 * Reset the Polygon from a CSV file
 * @param filename Filename
 * @param csv      CSV characteristics
 * @param verbose  Verbose flag
 * @param ncol_max Maximum number of columns
 * @param nrow_max Maximum number of rows
 * @return
 */
int Polygons::resetFromCSV(const String& filename,
                           const CSVformat& csv,
                           int verbose,
                           int ncol_max,
                           int nrow_max)
{
  VectorString names;
  VectorDouble tab;
  int ncol, nrow;

  // Free the previous contents

  _polysets.clear();

  /* Reading the CSV file: the coordinates are supposed to be in the first two columns */

  if (csv_table_read(filename, csv, verbose, ncol_max, nrow_max, &ncol, &nrow, names, tab))
  {
    messerr("Problem when reading CSV file");
    return 1;
  }
  if (ncol < 2)
  {
    messerr("The CSV file must contain at least 2 columns");
    return 1;
  }

  // Loop on the contents of the first column to look for Polysets
  int ideb = 0;
  int ifin = nrow;
  for (int i = 0; i < nrow; i++)
  {
    if (FFFF(tab[ncol * i]))
    {
      PolySet polyset = _extractFromTab(ideb, i, ncol, tab);
      addPolySet(polyset);
      ideb = i + 1;
    }
  }
  if (ideb < ifin)
  {
    PolySet polyset = _extractFromTab(ideb, nrow, ncol, tab);
    addPolySet(polyset);
  }
  return 0;
 }

void Polygons::addPolySet(const PolySet& polyset)
{
  PolySet local = polyset;
  local.closePolySet();
  _polysets.push_back(local);
}

String Polygons::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  int npol = static_cast<int> (_polysets.size());

  sstr << toTitle(1, "Polygons");
  sstr << "Number of Polygon Sets = " << npol << std::endl;
  AStringFormat sf;
  if (strfmt != nullptr) sf = *strfmt;

  if (sf.getLevel() > 1)
  {
    for (int i=0; i<npol; i++)
    {
      sstr << toTitle(2, "Polyset #%d", i+1);
      sstr << _polysets[i].toString(strfmt);
    }
  }
  return sstr.str();
}

void Polygons::getExtension(double *xmin,
                            double *xmax,
                            double *ymin,
                            double *ymax) const
{
  double xmin_loc, xmax_loc, ymin_loc, ymax_loc;

  for (int ipol = 0; ipol < getPolySetNumber(); ipol++)
  {
    _polysets[ipol].getExtension(&xmin_loc, &xmax_loc, &ymin_loc, &ymax_loc);
    if (xmin_loc < *xmin) (*xmin) = xmin_loc;
    if (ymin_loc < *ymin) (*ymin) = ymin_loc;
    if (xmax_loc > *xmax) (*xmax) = xmax_loc;
    if (ymax_loc > *ymax) (*ymax) = ymax_loc;
  }
}

double Polygons::getSurface() const
{
  double surface = 0.;
  for (int ipol = 0; ipol < getPolySetNumber(); ipol++)
  {
    surface += _polysets[ipol].getSurface();
  }
  return (surface);
}

PolySet Polygons::_extractFromTab(int ideb,
                                  int ifin,
                                  int ncol,
                                  const VectorDouble& tab)
{
  int nval = ifin - ideb;
  VectorDouble x(nval);
  VectorDouble y(nval);
  for (int j = ideb; j < ifin; j++)
  {
    int i = j - ideb;
    x[i] = tab[ncol * j + 0];
    y[i] = tab[ncol * j + 1];
  }
  PolySet polyset = PolySet(x,y);
  return polyset;
}

int Polygons::_deserialize(std::istream& is, bool verbose)
{
  int npol;

  // Clear previous contents

  _polysets.clear();

  /* Create the Model structure */

  bool ret = _recordRead<int>(is, "Number of Polygons", npol);
  if (! ret) return 1;

  /* Loop on the PolySets */

  for (int ipol = 0; ipol < npol; ipol++)
  {
    PolySet polyset;
    if (polyset._deserialize(is, verbose)) return 1;
    addPolySet(polyset);
  }
  return 0;
}

int Polygons::_serialize(std::ostream& os, bool verbose) const
{

  /* Create the Model structure */

  bool ret = _recordWrite<int>(os, "Number of Polygons", getPolySetNumber());

  /* Writing the covariance part */

  for (int ipol = 0; ipol < getPolySetNumber(); ipol++)
  {
    const PolySet& polyset = getPolySet(ipol);
    if (polyset._serialize(os, verbose)) return 1;
  }
  return ret ? 0 : 1;
}

Polygons* Polygons::create()
{
  return new Polygons();
}

/**
 * Create a Polygon by loading the contents of a Neutral File
 * @param neutralFilename Name of the Neutral File
 * @param verbose         Verbose flag
 * @return
 */
Polygons* Polygons::createFromNF(const String& neutralFilename, bool verbose)
{
  Polygons* polygons = nullptr;
  std::ifstream is;
  if (_fileOpenRead(neutralFilename, "Polygon", is, verbose))
  {
    polygons = new Polygons();
    if (polygons->_deserialize(is, verbose))
    {
      if (verbose) messerr("Problem reading the Neutral File.");
      delete polygons;
      polygons = nullptr;
    }
    is.close();
  }
  return polygons;
}

int Polygons::dumpToNF(const String& neutralFilename, bool verbose) const
{
  std::ofstream os;
  int ret = 1;
  if (_fileOpenWrite(neutralFilename, "Polygon", os, verbose))
  {
    ret = _serialize(os, verbose);
    if (ret && verbose) messerr("Problem writing in the Neutral File.");
    os.close();
  }
  return ret;
}

Polygons* Polygons::createFromCSV(const String& filename,
                                  const CSVformat& csv,
                                  int verbose,
                                  int ncol_max,
                                  int nrow_max)
{
  Polygons* polygons = new Polygons();
  if (polygons->resetFromCSV(filename, csv, verbose, ncol_max, nrow_max))
  {
    if (verbose) messerr("Problem reading the CSV File.");
    delete polygons;
    return nullptr;
  }
  return polygons;
}
Polygons* Polygons::createFromDb(const Db* db)
{
  Polygons* polygons = new Polygons();
  if (polygons->resetFromDb(db))
  {
    messerr("Problem building Polygons from DB.");
    delete polygons;
    return nullptr;
  }
  return polygons;
}
