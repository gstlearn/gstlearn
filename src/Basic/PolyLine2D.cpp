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
#include "geoslib_f_private.h"
#include "geoslib_old_f.h"

#include "Basic/PolyLine2D.hpp"
#include "Basic/NamingConvention.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"

PolyLine2D::PolyLine2D(const VectorDouble& x,
               const VectorDouble& y)
    : AStringable(),
      ASerializable(),
      _x(x),
      _y(y)
{
  // Check that both vectors have the same dimension
  if (x.size() != y.size())
  {
    _x.clear();
    _y.clear();
  }
}

PolyLine2D::PolyLine2D(const PolyLine2D &m)
    : AStringable(m),
      ASerializable(m),
      _x(m._x),
      _y(m._y)
{

}

PolyLine2D& PolyLine2D::operator=(const PolyLine2D &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    ASerializable::operator=(m);
    _x = m._x;
    _y = m._y;
  }
  return *this;
}

PolyLine2D::~PolyLine2D()
{

}

PolyLine2D* PolyLine2D::create(const VectorDouble& x, const VectorDouble& y)
{
  return new PolyLine2D(x, y);
}

String PolyLine2D::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  int npoints = getNPoints();

  VectorDouble tab = VectorDouble(2 * npoints);
  for (int i = 0; i < npoints; i++)
  {
    tab[i] = _x[i];
    tab[i + npoints] = _y[i];
  }
  sstr << toMatrix("Line Vertex Coordinates", VectorString(), VectorString(),
                   true, 2, npoints, tab);
  return sstr.str();
}

PolyLine2D* PolyLine2D::createFromNF(const String& neutralFilename, bool verbose)
{
  PolyLine2D* line2D = nullptr;
  std::ifstream is;
  line2D = new PolyLine2D();
  bool success = false;
  if (line2D->_fileOpenRead(neutralFilename, is, verbose))
  {
    success =  line2D->deserialize(is, verbose);
  }
  if (! success)
  {
    delete line2D;
    line2D = nullptr;
  }
  return line2D;
}

/**
 * Serialization (by Point rather than by Coordinate)
 * This is maintained for all classes using this interface for serialization
 * @param os Output Stream
 * @param verbose Verbose flag
 * @return
 */
bool PolyLine2D::_serialize(std::ostream& os, bool /*verbose*/) const
{
  if (getNPoints() <= 0) return false;
  bool ret = true;
  ret = ret && _recordWrite<int>(os, "Number of Points", (int) _x.size());

  VectorDouble buffer(2);
  for (int i = 0; i < (int) _x.size(); i++)
  {
    buffer[0] = _x[i];
    buffer[1] = _y[i];
    ret = ret && _recordWriteVec<double>(os, "", buffer);
  }
  return ret;
}

/**
 * Deserialization (by sample)
 * @param is Input stream
 * @param verbose Verbose flag
 * @return
 */
bool PolyLine2D::_deserialize(std::istream& is, bool /*verbose*/)
{
  int np = 0;
  bool ret = true;
  VectorDouble buffer(2);
  ret = ret && _recordRead<int>(is, "Number of Points", np);
  if (np < 0)
  {
    messerr("Something wrong in your file: the Number of Points read is equal to %d",np);
    return false;
  }
  _x.resize(np);
  _y.resize(np);
  for (int i = 0; i < np; i++)
  {
    ret = ret && _recordReadVec<double>(is, "", buffer, 2);
    _x[i] = buffer[0];
    _y[i] = buffer[1];
  }
  return ret;
}

void PolyLine2D::init(const VectorDouble& x, const VectorDouble& y)
{
  int nvert = static_cast<int> (x.size());

  /* Load the new Line */

  _x.resize(nvert,0);
  _y.resize(nvert,0);

  /* Copy the arrays */

  for (int i=0; i<nvert; i++)
  {
    _x[i] = x[i];
    _y[i] = y[i];
  }
}

void PolyLine2D::addPoint(double x, double y)
{
  int n = getNPoints();
  _x.resize(n+1);
  _y.resize(n+1);
  _x[n] = x;
  _y[n] = y;
}

VectorDouble PolyLine2D::getPoint(int i) const
{
  VectorDouble vec(2);
  vec[0] = getX(i);
  vec[1] = getY(i);
  return vec;
}

/****************************************************************************/
/*!
 **  Returns the point of the PolyLine located at shortest distance
 **  from Target
 **
 ** \return PolyPoint2D structure
 **
 ** \param[in]  xy0      Coordinates of the target point
 **
 *****************************************************************************/
PolyPoint2D PolyLine2D::getPLIndex(const VectorDouble &xy0) const
{

  double xx, yy, dist;
  int nint;
  int nvert = getNPoints();

  // Structure allocation

  PolyPoint2D pldist;
  pldist.coor.resize(2);

  /* Dispatch */

  double dmin = 1.e30;
  for (int i = 0; i < nvert - 1; i++)
  {
    dist = GH::distancePointToSegment(xy0[0], xy0[1], getX(i), getY(i), getX(i + 1),
                                      getY(i + 1), &xx, &yy, &nint);
    if (ABS(dist) > dmin) continue;
    pldist.rank = i;
    pldist.coor[0] = xx;
    pldist.coor[1] = yy;
    pldist.dist = dmin = ABS(dist);
  }
  return pldist;
}

/****************************************************************************/
/*!
 **  Shift a point along a segment
 **
 ** \param[in]  xy1     Coordinates of the first point
 ** \param[in]  xy2     Coordinates of the second point
 ** \param[in]  ratio   Shifting ratio
 **
 ** \param[out] xy0     Shifted point
 **
 ** \remarks 'ratio' varies between 0 and 1
 ** \remarks When 'ratio' =0, (x0,y0) coincides with (x1,y1)
 ** \remarks When 'ratio'>=1, (x0,y0) coincides with (x2,y2)
 **
 *****************************************************************************/
void PolyLine2D::_shiftPoint(const VectorDouble& xy1,
                             const VectorDouble& xy2,
                             double ratio,
                             VectorDouble& xy0) const
{
  if (ratio <= 0.)
  {
    xy0 = xy1;
  }
  else if (ratio >= 1.)
  {
    xy0 = xy2;
  }
  else
  {
    xy0[0] = xy1[0] + ratio * (xy2[0] - xy1[0]);
    xy0[1] = xy1[1] + ratio * (xy2[1] - xy1[1]);
  }
}

/****************************************************************************/
/*!
 **  Find the shortest distance between two points (x1,y1) and (x2,y2)
 **  passing through a polyline
 **
 ** \return Minimum distance
 **
 ** \param[in]  ap      Coefficient applied to the projected distances
 ** \param[in]  al      Coefficient applied to the distance along line
 ** \param[in]  xy1     Coordinates of the first point
 ** \param[in]  xy2     Coordinates of the second point
 **
 *****************************************************************************/
double PolyLine2D::distanceBetweenPoints(double ap,
                                         double al,
                                         const VectorDouble &xy1,
                                         const VectorDouble &xy2) const
{
  double dist, d1, d2, dh, dv, dloc, dmin, dist1, dist2;
  VectorDouble xyp1(2), xyp2(2);

  /* Calculate the projection of each end point */

  PolyPoint2D pldist1 = getPLIndex(xy1);
  PolyPoint2D pldist2 = getPLIndex(xy2);

  /* Calculate the minimum distance */

  dist = 1.e30;
  dh = dv = 0.;
  d1 = pldist1.dist;
  d2 = pldist2.dist;
  dh = ap * ABS(d1 - d2);

  if (al > 0.)
  {
    dv = distanceBetweenPlIndices(pldist1, pldist2);
    d1 = ABS(d1);
    d2 = ABS(d2);
    dmin = MIN(d1, d2);

    dist1 = ut_distance(2, pldist1.coor.data(), pldist2.coor.data());
    if (ABS(d1) > 0.) _shiftPoint(pldist1.coor, xy1, dmin / d1, xyp1);
    if (ABS(d2) > 0.) _shiftPoint(pldist2.coor, xy2, dmin / d2, xyp2);
    dist2 = ut_distance(2, xyp1.data(), xyp2.data());
    dv = (dist1 <= 0.) ? 0. : dv * al * sqrt(dist2 / dist1);
  }
  dloc = sqrt(dh * dh + dv * dv);
  if (dloc < dist) dist = dloc;

  return (dist);
}

/****************************************************************************/
/*!
 **  Find the shortest distance between two points (x1,y1) and (x2,y2)
 **  which belong to the current polyline
 **
 ** \return Minimum distance
 **
 ** \param[in]  pldist1 First PolyPoint2D structure
 ** \param[in]  pldist2 Second PolyPoint2D structure
 **
 *****************************************************************************/
double PolyLine2D::distanceBetweenPlIndices(const PolyPoint2D &pldist1,
                                            const PolyPoint2D &pldist2) const
{
  int i;
  double dist, local1[2], local2[2];
  PolyPoint2D pl1,pl2;

  /* Initializations */

  dist = 0.;
  if (pldist1.rank < pldist2.rank)
  {
    pl1 = pldist1;
    pl2 = pldist2;
  }
  else
  {
    pl1 = pldist2;
    pl2 = pldist1;
  }

  /* If both projected points belong to the same segment */

  if (pl1.rank == pl2.rank)
  {
    dist += ut_distance(2, pl1.coor.data(), pl2.coor.data());
  }
  else
  {

    /* Distance on the first segment */

    local1[0] = getX(pl1.rank + 1);
    local1[1] = getY(pl1.rank + 1);
    dist += ut_distance(2, pl1.coor.data(), local1);

    /* Distance on the last segment */

    local2[0] = getX(pl2.rank + 1);
    local2[1] = getY(pl2.rank + 1);
    dist += ut_distance(2, pl2.coor.data(), local2);

    for (i = pl1.rank + 1; i < pl2.rank; i++)
    {
      local1[0] = getX(i + 1);
      local1[1] = getY(i + 1);
      local2[0] = getX(i);
      local2[1] = getY(i);
      dist += ut_distance(2, local1, local2);
    }
  }
  return (dist);
}

double PolyLine2D::angleAtPolyline(const PolyPoint2D &pldist) const
{
  int r1, r2;
  int npoint = getNPoints();
  int rank = pldist.rank;
  if (rank < npoint - 1)
  {
    r1 = rank;
    r2 = rank + 1;
  }
  else
  {
    r1 = rank - 1;
    r2 = rank;
  }

  double incr_y = getY(r2) - getY(r1);
  double incr_x = getX(r2) - getX(r1);
  double angle  = atan2(incr_y, incr_x) * 180. / GV_PI;
  return angle;
}

/*****************************************************************************/
/*!
 **  Unfold a 2-D Db with respect to a polyline
 **
 ** \return  Error return code
 **
 ** \param[in]  db       Db structure
 ** \param[in]  polyline PolyLine2D structure
 ** \param[in]  namconv  Naming convention
 **
 *****************************************************************************/
int dbUnfoldPolyline(Db *db,
                     const PolyLine2D &polyline,
                     const NamingConvention &namconv)
{
  VectorDouble target(2);

  /* Initializations */

  int nvert = polyline.getNPoints();

  /* Preliminary checks */

  if (db->getNDim() != 2)
  {
    messerr("This function is restricted to 2-D Db");
    return 1;
  }
  if (nvert <= 1)
  {
    messerr("This function requires a Polyline with at least one segment");
    return 1;
  }

  /* Add the variables */

  int iptr = db->addColumnsByConstant(2, 0.);
  if (iptr < 0) return 1;

  /* Project the starting point */

  PolyPoint2D pldist0 = polyline.getPLIndex(polyline.getPoint(0));

  /* Loop on the samples of the Db */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    target[0] = db->getCoordinate(iech, 0);
    target[1] = db->getCoordinate(iech, 1);
    PolyPoint2D pldist = polyline.getPLIndex(target);
    double newx = pldist.dist;
    double newy = polyline.distanceBetweenPlIndices(pldist0, pldist);
    db->setArray(iech, iptr, newx);
    db->setArray(iech, iptr + 1, newy);
  }

  /* Set the error return code */

  namconv.setNamesAndLocators(db, ELoc::Z, -1, db, iptr);
  return 0;
}

/*****************************************************************************/
/*!
 **  Fold an input Db into an output Db with respect to a polyline
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin    Input Db structure
 ** \param[in]  dbout   Output Db structure
 ** \param[in]  cols    Vector of the target variable ranks
 ** \param[in]  polyline PolyLine2D structure
 ** \param[in]  namconv  Naming convention
 **
 *****************************************************************************/
int dbFoldPolyline(DbGrid *dbin,
                   Db *dbout,
                   const VectorInt& cols,
                   const PolyLine2D &polyline,
                   const NamingConvention &namconv)
{
  VectorDouble coor(2);

  /* Initializations */

  int nvert = polyline.getNPoints();

  /* Preliminary checks */

  if (dbin->getNDim() != 2 || ! dbin->isGrid())
  {
    messerr("This function is restricted to 2-D Input Grid Db");
    return 1;
  }
  if (dbout->getNDim() != 2)
  {
    messerr("This function is restricted to 2-D Output Db");
    return 1;
  }
  if (nvert <= 1)
  {
    messerr("This function requires a PolyLine2D with at least one segment");
    return 1;
  }

  /* Add the variables */

  int ncol = (int) cols.size();
  int iptr = dbout->addColumnsByConstant(ncol, TEST);
  if (iptr < 0) return 1;

  /* Project the starting point */

  PolyPoint2D pldist0 = polyline.getPLIndex(polyline.getPoint(0));

  /* Loop on the samples of the output Db */

  VectorDouble target(2);
  for (int iech = 0; iech < dbout->getSampleNumber(); iech++)
  {
    if (!dbout->isActive(iech)) continue;
    target[0] = dbout->getCoordinate(iech, 0);
    target[1] = dbout->getCoordinate(iech, 1);

    /* Project the target point according to the line */

    PolyPoint2D pldist = polyline.getPLIndex(target);
    coor[0] = pldist.dist;
    coor[1] = polyline.distanceBetweenPlIndices(pldist0, pldist);

    /* Locate the sample on the Input Grid */

    int iad = dbin->coordinateToRank(coor);
    if (iad < 0) continue;

    /* Loop on the variables */

    for (int icol = 0; icol < ncol; icol++)
    {
      double value = dbin->getArray(iad, cols[icol]);
      dbout->setArray(iech, iptr + icol, value);
    }
  }

  /* Set the error return code */

  namconv.setNamesAndLocators(dbout, ELoc::Z, -1, dbout, iptr);

  return 0;
}

/**
 * Calculate quantities on the Db by comparison with top and bottom polylines
 * @param db    Pointer to the Db where relevant information will be stored
 * @param top   2-D Polyline defining the Top surface
 * @param bot   2-D Polyline defining the Bottom surface
 * @param namconv Naming convention
 * @return Error return code
 */
int dbFromPolylines(Db* db,
                    const PolyLine2D& top,
                    const PolyLine2D& bot,
                    const NamingConvention &namconv)
{
  VectorDouble target(2);

  if (db == nullptr)
  {
    messerr("You must provide a Data Base for this method");
    return 1;
  }
  if (db->getNDim() != 2)
  {
    messerr("This method is restricted to a 2-D Data Base");
    return 1;
  }

  // Allocate the output variables

  int iptr = db->addColumnsByConstant(2, TEST);
  if (iptr < 0) return 1;

  // Loop on the active samples of the Data Base

  int nech = db->getSampleNumber();
  for (int iech = 0; iech < nech; iech++)
  {
    if (! db->isActive(iech)) continue;
    target[0] = db->getCoordinate(iech, 0);
    target[1] = db->getCoordinate(iech, 1);

    // Project the target point on the two polylines
    PolyPoint2D pl_top = top.getPLIndex(target);
    PolyPoint2D pl_bot = bot.getPLIndex(target);

    // Get distances to the polyline projections
    double dist_top = pl_top.dist;
    double dist_bot = pl_bot.dist;

    // Get the angles at the polyline projections
    double angle_top = top.angleAtPolyline(pl_top);
    double angle_bot = bot.angleAtPolyline(pl_bot);

    // Calculate relevant quantities
    double dist  = dist_bot + dist_top;
    double angle = (angle_top * dist_bot + angle_bot * dist_top) / (dist_bot + dist_top);

    db->setArray(iech, iptr  , dist);
    db->setArray(iech, iptr+1, angle);
  }

  /* Set the error return code */

  namconv.setNamesAndLocators(db, iptr  , "Dist");
  namconv.setNamesAndLocators(db, iptr+1, "Angle");

  return 0;
}
