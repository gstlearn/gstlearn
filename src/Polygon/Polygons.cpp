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
#include "Db/Db.hpp"
#include "Basic/AStringable.hpp"
#include "Polygon/Polygons.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/CSVformat.hpp"
#include "Basic/AException.hpp"
#include "geoslib_f.h"

Polygons::Polygons()
  : _polysets()
{
}

Polygons::Polygons(const String& filename,
                   const CSVformat& csv,
                   int verbose,
                   int ncol_max,
                   int nrow_max,
                   int flag_add_rank)
    : AStringable(),
      ASerializable(),
      _polysets()
{
  VectorString names;
  VectorDouble tab;
  int ncol, nrow;

  /* Reading the CSV file: the coordinates are supposed to be in the first two columns */

  if (csv_table_read(filename.c_str(), verbose, csv.getFlagHeader(), csv.getNSkip(),
                     csv.getCharSep().c_str(), csv.getCharDec().c_str(), csv.getNaString().c_str(),
                     ncol_max, nrow_max, &ncol, &nrow, names, tab))
  {
    messerr("Problem when reading CSV file");
    return;
  }

  if (ncol < 2)
  {
    messerr("The CSV file must contain at least 2 columns");
    return;
  }

  // Loop on the contents of the first column to look for Polysets
  int ideb;
  int ifin;
  ideb = ifin = 0;
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
 }

/**
 * Constructor as the convex Hull of the points contained in the Db
 * @param db
 */
Polygons::Polygons(const Db* db)
    : AStringable(),
      ASerializable(),
      _polysets()
{
  VectorDouble x;
  VectorDouble y;
  (void) polygon_hull(db, x, y);

  PolySet polyset = PolySet(x, y, TEST, TEST);
  addPolySet(polyset);
}

Polygons::Polygons(const String& neutralFileName, bool verbose)
    : AStringable(),
      ASerializable(),
      _polysets()
{
  if (deSerialize(neutralFileName, verbose))
    my_throw("Problem reading the Neutral File");
}

Polygons::Polygons(const Polygons& r)
    : _polysets(r._polysets)
{
}

Polygons& Polygons::operator=(const Polygons& r)
{
  if (this != &r)
  {
    _polysets = r._polysets;
  }
  return *this;
}

Polygons::~Polygons()
{
}

void Polygons::addPolySet(const PolySet& polyset)
{
  _polysets.push_back(polyset);
}

std::string Polygons::toString(int level) const
{
  std::stringstream sstr;

  int npol = static_cast<int> (_polysets.size());

  sstr << toTitle(1, "Polygons");
  sstr << "Number of Polygon Sets = " << npol << std::endl;

  for (int i=0; i<npol; i++)
  {
    sstr << toTitle(2, "Polyset #%d", i+1);
    sstr << _polysets[i].toString(level);
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

int Polygons::deSerialize(const String& filename, bool verbose)
{
  int npol, nvert;
  double zmin = TEST;
  double zmax = TEST;

  /* Opening the Data file */

  if (_fileOpen(filename, "Polygon", "r", verbose)) return 1;

  /* Create the Model structure */

  if (_recordRead("Number of Polygons", "%d", &npol)) return 1;

  /* Loop on the PolySets */

  for (int ipol = 0; ipol < npol; ipol++)
  {
    if (_recordRead("Number of Vertices", "%d", &nvert)) return 1;
    VectorDouble x(nvert);
    VectorDouble y(nvert);

    /* Loop on the Vertices */

    for (int i = 0; i < nvert; i++)
    {
      if (_recordRead("X-Coordinate of a Polyset", "%lf", &x[i])) return 1;
      if (_recordRead("Y-Coordinate of a Polyset", "%lf", &y[i])) return 1;
    }

    /* Add the polyset */

    PolySet polyset = PolySet();
    polyset.init(x,y,zmin,zmax);
    addPolySet(polyset);
  }

  _fileClose(verbose);
  return 0;
}

int Polygons::serialize(const String& filename, bool verbose)
{
  if (_fileOpen(filename, "Polygon", "w", verbose)) return (1);

  /* Create the Model structure */

  _recordWrite("%d", getPolySetNumber());
  _recordWrite("#", "Number of Polygons");

  /* Writing the covariance part */

  for (int ipol = 0; ipol < getPolySetNumber(); ipol++)
  {
    const PolySet& polyset = getPolySet(ipol);
    _recordWrite("%d", polyset.getNVertices());
    _recordWrite("#", "Number of Vertices");

    for (int i = 0; i < polyset.getNVertices(); i++)
    {
      _recordWrite("%lf", polyset.getX(i));
      _recordWrite("%lf", polyset.getY(i));
      _recordWrite("\n");
    }
  }

  // Close the Neutral File
  _fileClose(verbose);

  return 0;
}
