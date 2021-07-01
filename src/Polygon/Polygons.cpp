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
#include "geoslib_f.h"

Polygons::Polygons()
  : _polysets()
{
}

Polygons::Polygons(const String& filename,
                   int flag_header,
                   int nskip,
                   const String& char_sep,
                   const String& char_dec,
                   const String& na_string,
                   int verbose,
                   int ncol_max,
                   int nrow_max,
                   int flag_add_rank)
    : _polysets()
{
  VectorString names;
  VectorDouble tab;
  int ncol, nrow;

  /* Reading the CSV file: the coordinates are supposed to be in the first two columns */

  if (csv_table_read(filename.c_str(), verbose, flag_header, nskip,
                     char_sep.c_str(), char_dec.c_str(), na_string.c_str(),
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
    : _polysets()
{
  VectorDouble x;
  VectorDouble y;
  (void) polygon_hull(db, x, y);

  PolySet polyset = PolySet(x, y, TEST, TEST);
  addPolySet(polyset);
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

  int npol = _polysets.size();

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
