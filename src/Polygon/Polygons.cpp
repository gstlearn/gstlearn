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
#include "geoslib_old_f.h"

#include "Core/CSV.hpp"
#include "Db/Db.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/File.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/CSVformat.hpp"
#include "Polygon/Polygons.hpp"

Polygons::Polygons()
  : _polyelems(),
    _emptyVec(),
    _emptyElem()
{
}

Polygons::Polygons(const Polygons& r)
    : AStringable(r),
      ASerializable(r),
      _polyelems(r._polyelems),
      _emptyVec(r._emptyVec),
      _emptyElem(r._emptyElem)
{
}

Polygons& Polygons::operator=(const Polygons& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    ASerializable::operator=(r);
    _polyelems = r._polyelems;
    _emptyVec = r._emptyVec;
    _emptyElem = r._emptyElem;
  }
  return *this;
}

Polygons::~Polygons()
{
}

int Polygons::resetFromDb(const Db* db, double dilate, bool verbose)
{
  if (db == nullptr) return 1;

  Polygons *polygons = new Polygons();
  if (polygons->_buildHull(db, dilate, verbose))
  {
    delete polygons;
    return 1;
  }

  *this = *polygons;
  delete polygons;

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

  _polyelems.clear();

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

  // Loop on the contents of the first column to look for polyelems
  int ideb = 0;
  int ifin = nrow;
  for (int i = 0; i < nrow; i++)
  {
    if (FFFF(tab[ncol * i]))
    {
      PolyElem polyelem = _extractFromTab(ideb, i, ncol, tab);
      addPolyElem(polyelem); // TODO : Prevent copy (optimization)
      ideb = i + 1;
    }
  }
  if (ideb < ifin)
  {
    PolyElem polyelem = _extractFromTab(ideb, nrow, ncol, tab);
    addPolyElem(polyelem); // TODO : Prevent copy (optimization)
  }
  return 0;
 }

/**
 * Reset the Polygon from a CSV file using WKT format (first column) exported by QGIS
 * @param filename Filename
 * @param csv      CSV characteristics
 * @param verbose  Verbose flag
 * @param ncol_max Maximum number of columns
 * @param nrow_max Maximum number of rows
 * @return
 */
int Polygons::resetFromWKT(const String& filename,
                           const CSVformat& csv,
                           int verbose,
                           int ncol_max,
                           int nrow_max)
{
  DECLARE_UNUSED(nrow_max);
  DECLARE_UNUSED(ncol_max);
  DECLARE_UNUSED(verbose);

  // Free the previous contents

  _polyelems.clear();

  /* Reading the file line by line, the first column is supposed to be WKT */
  // Open the journal
  std::fstream file(filename, std::ios_base::in);
  if (!file.is_open()) {
    return 1;
  }

  std::string line;
  std::string poly;
  std::string polyelem;
  // Has header ?
  if (csv.getFlagHeader())
  {
    int i = 0;
    // Skip first lines
    //while(std::getline(file, line) && i < csv.getNSkip())
    while(gslSafeGetline(file, line) && i < csv.getNSkip())
      i++;
    // Header too short
    if (i < csv.getNSkip())
      return 1;
  }
  // Parse next lines
  while(std::getline(file, line)) {
    // Cleanup line
    line = trim(line);
    // Ignore empty line
    if (line.length() == 0) continue;
    // Ignore comments
    if (line.find_first_of('#') == 0) continue;
    // Cleanup column header
    size_t found = line.find("(((");
    // If there is no "((("
    if(found == std::string::npos)
      return 1;
    poly = line.substr(found+3);
    // Cleanup trailing stuff
    found = poly.find(")))");
    // If there is no ")))"
    if(found == std::string::npos)
      return 1;
    poly = poly.substr(0,found);
    // Read each polygon elements from the line
    do
    {
      found = poly.find(")),((");
      if (found != std::string::npos)
        polyelem = poly.substr(0, found);
      else
        polyelem = poly;
      // Extract coordinates
      PolyElem pe = _extractFromWKT(csv, polyelem);
      addPolyElem(pe); // TODO : Prevent copy (optimization)
      // Look for next polygon element
      if (found != std::string::npos)
        poly = poly.substr(found+5);
    }
    while(found != std::string::npos);
  }

  return 0;
 }

/**
 * Add the PolyElem to the list of polygons.
 * This is performed only if the current PolyElem contains at least 3 vertices.
 * @param polyelem
 */
void Polygons::addPolyElem(const PolyElem& polyelem)
{
  if (polyelem.getNPoints() >= 3) _polyelems.push_back(polyelem);
}

String Polygons::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  int npol = static_cast<int> (_polyelems.size());

  sstr << toTitle(1, "Polygons");
  sstr << "Number of Polygon Sets = " << npol << std::endl;
  AStringFormat sf;
  if (strfmt != nullptr) sf = *strfmt;

  if (sf.getLevel() > 1)
    for (int i=0; i<npol; i++)
    {
      sstr << toTitle(2, "PolyElem #%d", i+1);
      sstr << _polyelems[i].toString(strfmt);
    }
  return sstr.str();
}

void Polygons::getExtension(double *xmin,
                            double *xmax,
                            double *ymin,
                            double *ymax) const
{
  double xmin_loc, xmax_loc, ymin_loc, ymax_loc;

  for (int ipol = 0; ipol < getPolyElemNumber(); ipol++)
  {
    _polyelems[ipol].getExtension(&xmin_loc, &xmax_loc, &ymin_loc, &ymax_loc);
    if (xmin_loc < *xmin) (*xmin) = xmin_loc;
    if (ymin_loc < *ymin) (*ymin) = ymin_loc;
    if (xmax_loc > *xmax) (*xmax) = xmax_loc;
    if (ymax_loc > *ymax) (*ymax) = ymax_loc;
  }
}

double Polygons::getSurface() const
{
  double surface = 0.;
  for (int ipol = 0; ipol < getPolyElemNumber(); ipol++)
  {
    surface += _polyelems[ipol].getSurface();
  }
  return (surface);
}

PolyElem Polygons::_extractFromTab(int ideb,
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
  PolyElem polyelem = PolyElem(x,y);
  return polyelem;
}

PolyElem Polygons::_extractFromWKT(const CSVformat& csv, String& polye)
{
  VectorDouble x;
  VectorDouble y;
  size_t found, found2;
  std::string coords;
  do
  {
    found = polye.find_first_of(csv.getCharSep());
    if (found != std::string::npos)
      coords = polye.substr(0, found);
    else
      coords = polye;
    found2 = coords.find_first_of(' ');
    if (found2 == std::string::npos)
    {
      x.clear();
      y.clear();
      break;
    }
    x.push_back(toDouble(coords.substr(0, found2), csv.getCharDec()));
    coords = coords.substr(found2+1);
    found2 = coords.find_first_of(' ');
    if (found2 != std::string::npos)
      coords = coords.substr(0, found2);
    y.push_back(toDouble(coords, csv.getCharDec()));
    if (found != std::string::npos)
      polye = polye.substr(found+1);
  }
  while (found != std::string::npos);

  PolyElem polyelem = PolyElem(x,y);
  return polyelem;
}


bool Polygons::_deserialize(std::istream& is, bool verbose)
{
  int npol = 0;

  // Clear previous contents

  _polyelems.clear();

  /* Create the Model structure */

  bool ret = true;
  ret = ret && _recordRead<int>(is, "Number of Polygons", npol);

  /* Loop on the PolyElems */

  for (int ipol = 0; ret && ipol < npol; ipol++)
  {
    PolyElem polyelem;
    ret = ret && polyelem._deserialize(is, verbose);
    if (ret)
    {
      addPolyElem(polyelem); // TODO : Prevent copy (optimization)
      if (verbose) message("PolyElem #%d - Number of vertices = %d\n",ipol+1, polyelem.getNPoints());
    }
    else
      messerr("Error when reading PolyElem #%d",ipol+1);
  }
  return ret;
}

bool Polygons::_serialize(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret = ret && _recordWrite<int>(os, "Number of Polygons", getPolyElemNumber());

  /* Writing the covariance part */

  for (int ipol = 0; ret && ipol < getPolyElemNumber(); ipol++)
  {
    const PolyElem& polyelem = getPolyElem(ipol);
    ret = ret && polyelem._serialize(os, verbose);
  }
  return ret;
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
 */  VectorDouble _emptyVec; // dummy
 PolyElem     _emptyElem; // dummy
Polygons* Polygons::createFromNF(const String& neutralFilename, bool verbose)
{
  Polygons* polygons = nullptr;
  std::ifstream is;
  polygons = new Polygons();
  bool success = false;
  if (polygons->_fileOpenRead(neutralFilename, is, verbose))
  {
    success = polygons->deserialize(is, verbose);
  }
  if (! success)
  {
    delete polygons;
    polygons = nullptr;
  }
  return polygons;
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

Polygons* Polygons::createFromWKT(const String& filename,
                                  const CSVformat& csv,
                                  int verbose,
                                  int ncol_max,
                                  int nrow_max)
{
  Polygons* polygons = new Polygons();
  if (polygons->resetFromWKT(filename, csv, verbose, ncol_max, nrow_max))
  {
    if (verbose) messerr("Problem reading the CSV File (WKT).");
    delete polygons;
    return nullptr;
  }
  return polygons;
}

Polygons* Polygons::createFromDb(const Db* db, double dilate, bool verbose)
{
  Polygons* polygons = new Polygons();
  if (polygons->resetFromDb(db, dilate, verbose))
  {
    messerr("Problem building Polygons from DB.");
    delete polygons;
    return nullptr;
  }
  return polygons;
}

const PolyElem& Polygons::getPolyElem(int ipol) const
{
  if (! _isValidPolyElemIndex(ipol)) return _emptyElem;
  return _polyelems[ipol];
}

PolyElem Polygons::getClosedPolyElem(int ipol) const
{
  if (! _isValidPolyElemIndex(ipol)) return PolyElem();
  PolyElem polyelem = getPolyElem(ipol);
  polyelem.closePolyElem();
  return polyelem;
}

const VectorDouble& Polygons::getX(int ipol) const
{
  if (! _isValidPolyElemIndex(ipol)) return _emptyVec;
  return _polyelems[ipol].getX();
}

const VectorDouble& Polygons::getY(int ipol) const
{
  if (! _isValidPolyElemIndex(ipol)) return _emptyVec;
  return _polyelems[ipol].getY();
}

void Polygons::setX(int ipol, const VectorDouble& x)
{
  if (! _isValidPolyElemIndex(ipol)) return;
  return _polyelems[ipol].setX(x);
}

void Polygons::setY(int ipol, const VectorDouble& y)
{
  if (! _isValidPolyElemIndex(ipol)) return;
  return _polyelems[ipol].setY(y);
}

bool Polygons::_isValidPolyElemIndex(int ipol) const
{
  int npol = getPolyElemNumber();
  if (ipol < 0 || ipol >= npol)
  {
    messerr("PolyElem Index %d is not valid. It should lie in [0,%d[",
            ipol,npol);
    return false;
  }
  return true;
}

/*****************************************************************************/
/*!
 **  Determine the distance to a polyline
 **
 ** \return  Error returned code
 **
 ** \param[in]  db      Db structure
 ** \param[in]  polygon Polygons structure
 ** \param[in]  dmax    Maximum distance
 ** \param[in]  scale   Scaling option
 **                     0 : no scaling
 **                    >0 : scaling between 0 and 1
 **                    <0 : scaling between 1 and 0
 ** \param[in]  polin   Option for checking against the polygon
 **                     0 : no check
 **                    >0 : if sample is outside polygon, return TEST
 **                    <0 : if sample is inside polygon, return TEST
 ** \param[in] namconv Naming convention
 **
 ** \remarks When patching values with respect to the polygon, when abs(polin):
 ** \remarks 1 : put NA
 ** \remarks 2 : put minimum distance
 ** \remarks 3 : put maximum distance
 **
 *****************************************************************************/
int dbPolygonDistance(Db *db,
                      Polygons *polygon,
                      double dmax,
                      int scale,
                      int polin,
                      const NamingConvention &namconv)
{
  PolyPoint2D pldist;
  VectorDouble target(2);

  // Initializations

  int nech = db->getNSample();

  // Create a new attribute

  int iptr = db->addColumnsByConstant(1, TEST);
  if (iptr < 0) return (1);

  // Loop on the polyelems

  double distmax = 0.;
  double distmin = TEST;
  for (int iset = 0; iset < polygon->getPolyElemNumber(); iset++)
  {
    const PolyElem &polyelem = polygon->getPolyElem(iset);

    // Loop on the samples

    for (int iech = 0; iech < nech; iech++)
    {
      if (!db->isActive(iech)) continue;
      target[0] = db->getCoordinate(iech, 0);
      target[1] = db->getCoordinate(iech, 1);
      PolyLine2D polyline(polyelem.getX(), polyelem.getY());
      pldist = polyline.getPLIndex(target);
      double distloc = pldist.dist;
      if (FFFF(distloc)) continue;
      distmin = db->getArray(iech, iptr);
      if (FFFF(distmin))
        distmin = distloc;
      else
        distmin = MIN(distmin, distloc);
      if (!FFFF(dmax) && distmin > dmax) distmin = dmax;
      db->setArray(iech, iptr, distmin);
    }
  }

  // Calculate the extreme values

  if (scale != 0 || polin != 0)
  {
    distmin = 1.e30;
    distmax = 0.;
    for (int iech = 0; iech < nech; iech++)
    {
      if (!db->isActive(iech)) continue;
      double distloc = db->getArray(iech, iptr);
      if (FFFF(distloc)) continue;
      if (polin != 0)
      {
        int inside = polygon->inside(db->getCoordinate(iech, 0),
                                     db->getCoordinate(iech, 1));
        if (polin > 0)
        {
          if (!inside) continue;
        }
        else
        {
          if (inside) continue;
        }
      }
      if (distloc > distmax) distmax = distloc;
      if (distloc < distmin) distmin = distloc;
    }
  }

  // Scaling option

  if (scale != 0)
  {
    if (scale > 0)
    {
      for (int iech = 0; iech < nech; iech++)
      {
        if (!db->isActive(iech)) continue;
        double distloc = db->getArray(iech, iptr);
        if (FFFF(distloc)) continue;
        double value = (distloc - distmin) / (distmax - distmin);
        db->setArray(iech, iptr, value);
      }
    }
    else
    {
      for (int iech = 0; iech < nech; iech++)
      {
        if (!db->isActive(iech)) continue;
        double distloc = db->getArray(iech, iptr);
        if (FFFF(distloc)) continue;
        double value = (distloc - distmax) / (distmin - distmax);
        db->setArray(iech, iptr, value);
      }
    }
    distmin = 0.;
    distmax = 1.;
  }

  // Check if the sample belongs to the polygon or not

  if (polin != 0)
  {
    double valtest = TEST;
    if (ABS(polin) == 2) valtest = distmin;
    if (ABS(polin) == 3) valtest = distmax;
    for (int iech = 0; iech < nech; iech++)
    {
      int inside = polygon->inside(db->getCoordinate(iech, 0),
                                   db->getCoordinate(iech, 1));
      if (polin > 0 && !inside) db->setArray(iech, iptr, valtest);
      if (polin < 0 &&  inside) db->setArray(iech, iptr, valtest);
    }
  }

  /* Set the error return code */

  namconv.setNamesAndLocators(db, VectorString(), ELoc::Z, -1, db, iptr);

  return (0);
}

/****************************************************************************/
/*!
 **  Check if one point belongs to a Polygons
 **
 ** \return  1 if the point belongs to the Polygons; 0 otherwise
 **
 ** \param[in]  coor        Vector of coordinates
 ** \param[in]  flag_nested Option for nested polyelems (see details)
 **
 ** \remarks When coor is dimensioned to 3, the third dimension test is performed
 ** \remarks If flag_nested=TRUE, a sample is masked off if the number of
 ** \remarks polyelems to which it belongs is odd
 ** \remarks If flag_nested=FALSE, a sample is masked off as soon as it
 ** \remarks belongs to one PolyElem
 **
 *****************************************************************************/
bool Polygons::inside(const VectorDouble& coor, bool flag_nested) const
{
  bool flag3d = (int) coor.size() > 2;
  if (flag_nested)
  {
    int number = 0;
    for (int ipol = 0; ipol < getPolyElemNumber(); ipol++)
    {
      PolyElem polyelem = getClosedPolyElem(ipol);
      if (flag3d)
      {
        if (!polyelem.inside3D(coor[2])) return false;
      }
      if (polyelem.inside(coor)) number++;
    }
    if (number % 2 != 0) return true;
  }
  else
  {
    for (int ipol = 0; ipol < getPolyElemNumber(); ipol++)
    {
      PolyElem polyelem = getClosedPolyElem(ipol);
      if (flag3d)
      {
        if (!polyelem.inside3D(coor[2])) return false;
      }
      if (polyelem.inside(coor)) return true;
    }
  }
  return false;
}

/*****************************************************************************/
/*!
 **  Create a polygon from the convex hull of active samples
 **
 ** \return  Error returned code
 **
 ** \param[in]  x    Vector of X coordinates
 ** \param[in]  y    Vector of Y coordinates
 **
 *****************************************************************************/
VectorInt Polygons::_getHullIndices(const VectorDouble &x,
                                    const VectorDouble &y)
{
  int number = (int) x.size();
  VectorInt index(number + 1);

  /* Calculate the center of gravity and the leftmost point */

  int rank = 0;
  int np = 0;
  double xg = 0.;
  double yg = 0.;
  for (int i = 0; i < number; i++)
  {
    xg += x[i];
    yg += y[i];
    if (x[i] < x[rank]) rank = i;
  }
  xg /= number;
  yg /= number;
  index[0] = rank;
  np++;

  /* Implicit loop for finding the other points of the convex hull */

  for (;;)
  {
    double x2, x3, y2, y3;
    double x1 = x[index[np - 1]];
    double y1 = y[index[np - 1]];
    double x21 = xg - x1;
    double y21 = yg - y1;

    for (int i = 0; i < number; i++)
    {
      if ((x[i] - x1) * y21 - (y[i] - y1) * x21 <= 0.) continue;
      x21 = x[i] - x1;
      y21 = y[i] - y1;
      rank = i;
    }
    if (rank == index[0]) goto label_cont;

    /* Discard the previous point if on the line joining the current point */
    /* to the one before the previous one */

    if (np > 1)
    {
      x1 = x[rank];
      y1 = y[rank];
      x2 = x[index[np - 1]];
      y2 = y[index[np - 1]];
      x3 = x[index[np - 2]];
      y3 = y[index[np - 2]];
      if (ABS((x2-x3)*(y1-y3) - (x1-x3)*(y2-y3)) < EPSILON6) np--;
    }
    index[np] = rank;
    np++;
  }

  label_cont:
  index[np++] = index[0];
  index.resize(np);

  return index;
}

void Polygons::_polygonHullPrint(const VectorInt &index,
                                 const VectorDouble &x,
                                 const VectorDouble &y)
{
  mestitle(1,"Polygon Hull");
  message("Ranks (1-based) and coordinates of the Active Samples included in the Convex Hull\n");
  for (int i = 0; i < (int) index.size(); i++)
  {
    int j = index[i];
    message("%3d : %lf %lf\n",j+1, x[j], y[j]);
  }
}

/**
 * Create a set of fictitious samples obtained by dilating the initial ones
 * @param ext Dilation distance
 * @param x   Vector of X-coordinates or initial samples
 * @param y   Vector of Y-coordinates of initial samples
 * @param nsect Number of discretization points for dilation
 */
void Polygons::_getExtend(double ext,
                          VectorDouble &x,
                          VectorDouble &y,
                          int nsect)
{
  int ninit = (int) x.size();

  x.resize(nsect * ninit);
  y.resize(nsect * ninit);

  for (int j = 0; j < ninit; j++)
  {
    int i = ninit - j - 1;
    double x0 = x[i];
    double y0 = y[i];

    int iad = i * nsect;
    for (int k = 0; k < nsect; k++)
    {
      double angle = 2. * GV_PI * k / (double) nsect;
      x[iad + k] = x0 + ext * cos(angle);
      y[iad + k] = y0 + ext * sin(angle);
    }
  }
}

/*****************************************************************************/
/*!
 **  Create a polygon from the convex hull of active samples
 **
 ** \return  Error returned code
 **
 ** \param[in]  db      descriptor of the Db serving for convex hull calculation
 ** \param[in]  dilate  Radius of the dilation
 ** \param[in]  verbose Verbose flag
 **
 *****************************************************************************/
int Polygons::_buildHull(const Db *db, double dilate, bool verbose)

{
  /* Preliminary check */

  if (db->getNDim() < 2)
  {
    messerr("The input Db must be contain at least 2 coordinates");
    return 1;
  }
  int number = db->getNSample(true);
  if (number <= 0)
  {
    messerr("No active data in the input Db. Convex Hull impossible");
    return 1;
  }

  // Load the vector of coordinates (of active samples)

  VectorDouble xinit = db->getColumnByLocator(ELoc::X, 0, true);
  VectorDouble yinit = db->getColumnByLocator(ELoc::X, 1, true);

  // Calculate the indices that are retained in the convex hull

  VectorInt index = _getHullIndices(xinit, yinit);

  // Optional printout

  if (verbose) _polygonHullPrint(index, xinit, yinit);

  /* Create the polygons */

  int np = (int) index.size();
  VectorDouble xret(np);
  VectorDouble yret(np);
  for (int i = 0; i < np; i++)
  {
    xret[i] = xinit[index[i]];
    yret[i] = yinit[index[i]];
  }

  // Extend to dilated hull (optional)



  {
    xinit = xret;
    yinit = yret;
    _getExtend(dilate, xinit, yinit);
    index = _getHullIndices(xinit, yinit);

    np = (int) index.size();
    xret.resize(np);
    yret.resize(np);
    for (int i = 0; i < np; i++)
    {
      xret[i] = xinit[index[i]];
      yret[i] = yinit[index[i]];
    }
  }

  // Load the information within the polygon structure

  PolyElem polyelem = PolyElem(xret, yret);
  addPolyElem(polyelem);

  return 0;
}

Polygons Polygons::reduceComplexity(double distmin) const
{
  Polygons newpolygon;

  for (int ipol = 0, npol = getPolyElemNumber(); ipol < npol; ipol++)
  {
    PolyElem newpolyelem = getPolyElem(ipol).reduceComplexity(distmin);
    newpolygon.addPolyElem(newpolyelem);
  }
  return newpolygon;
}

/****************************************************************************/
/*!
 **  Create a selection if the samples of a Db are inside Polygons
 **
 ** \param[in]  db          Db structure
 ** \param[in]  polygon     Polygons structure
 ** \param[in]  flag_sel    true if previous selection must be taken into account
 ** \param[in]  flag_period true if first coordinate is longitude (in degree) and
 **                         must be cycled for the check
 ** \param[in]  flag_nested Option for nested polyelems (see details)
 ** \param[in]  namconv     Naming Convention
 **
 ** \remarks If flag_nested=true, a sample is masked off if the number of
 ** \remarks polyelems to which it belongs is odd
 ** \remarks If flag_nested=false, a sample is masked off as soon as it
 ** \remarks belongs to one polyelem
 **
 ** \remark The Naming Convention locator Type is overwritten to ELoc::SEL
 **
 *****************************************************************************/
void db_polygon(Db *db,
                const Polygons *polygon,
                bool flag_sel,
                bool flag_period,
                bool flag_nested,
                const NamingConvention& namconv)
{
  // Adding a new variable

  int iatt = db->addColumnsByConstant(1);

  /* Loop on the samples */

  VectorDouble coor(3, TEST);
  for (int iech = 0; iech < db->getNSample(); iech++)
  {
    mes_process("Checking if sample belongs to a polygon",db->getNSample(),iech);
    int selval = 0;
    if (! flag_sel || db->isActive(iech))
    {
      db->getCoordinatesPerSampleInPlace(iech, coor);
      selval = polygon->inside(coor, flag_nested);

      if (flag_period)
      {
        double xx = coor[0];
        coor[0] = xx - 360;
        selval = selval || polygon->inside(coor, flag_nested);
        coor[0] = xx + 360;
        selval = selval || polygon->inside(coor, flag_nested);
      }
    }
    db->setArray(iech, iatt, (double) selval);
  }

  // Setting the output variable
  namconv.setNamesAndLocators(db, iatt);
}

/*****************************************************************************/
/*!
 **  Select samples from a file according to the 2-D convex hull
 **  computed over the active samples of a second file
 **
 ** \return  Error returned code
 **
 ** \param[in]  db1     descriptor of the Db serving for convex hull calculation
 ** \param[in]  db2     descriptor of the Db where the mask must be stored
 ** \param[in]  dilate  Radius of the dilation
 ** \param[in]  verbose Verbose flag
 ** \param[in]  namconv Naming convention
 **
 ** \remark The Naming Convention locator Type is overwritten to ELoc::SEL
 **
 *****************************************************************************/
int db_selhull(Db *db1,
               Db *db2,
               double dilate,
               bool verbose,
               const NamingConvention &namconv)
{
  /* Create the polygon as the convex hull of first Db */

  Polygons* polygons = Polygons::createFromDb(db1, dilate, verbose);
  if (polygons == nullptr) return 1;

  // Create the variable in the output Db

  int isel = db2->addColumnsByConstant(1, 1.);

  /* Loop on the samples of the second Db */
  // Note that all samples must be checked as a sample, initially masked, can be
  // masked OFF as it belongs to the convex hull.

  int ntotal = db2->getNSample();
  int nactive = 0;
  int nout = 0;
  int nin = 0;

  VectorDouble coor(3, TEST);
  for (int iech = 0; iech < ntotal; iech++)
  {
    db2->getCoordinatesPerSampleInPlace(iech, coor);
    if (!polygons->inside(coor, false))
    {
      db2->setArray(iech, isel, 0.);
      nout++;
    }
    else
    {
      nin++;
    }
  }

  // Verbose optional output
  if (verbose)
  {
    mestitle(1, "Convex Hull calculation");
    message("- Number of target samples = %d\n", ntotal);
    message("- Number of active samples = %d\n", nactive);
    message("- Number of masked samples = %d\n", nout);
    message("- Number of valid samples  = %d\n", nin);
  }

  // Set the Naming Convention
  namconv.setNamesAndLocators(db2, isel);

  delete polygons;
  return 0;
}

