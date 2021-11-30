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
#include "geoslib_e.h"
#include "Polygon/Polygons.hpp"
#include "Polygon/PolySet.hpp"
#include "Basic/Utilities.hpp"

/****************************************************************************/
/*!
 **  Create the Polygons structure
 **
 ** \return  Pointer to the newly created Polygonsn structure
 **
 *****************************************************************************/
Polygons* polygon_create(void)

{
  Polygons *polygon;

  /* Initializations */

  polygon = nullptr;

  /* Core allocation */

  polygon = new (Polygons);

  return (polygon);
}

/****************************************************************************/
/*!
 **  Free the Polygons structure
 **
 ** \return  Pointer to the newly freed Polygons structure
 **
 *****************************************************************************/
Polygons* polygon_free(Polygons *polygon)

{
  if (polygon == nullptr) return (polygon);
  delete polygon;
  polygon = (Polygons*) nullptr;
  return (polygon);
}

/****************************************************************************/
/*!
 **  Add a PolySet to an existing Polygons structure
 **
 ** \return  Pointer to the Polygons structure
 **
 ** \param[in]  polygon  Polygons structure where the PolySet should be added
 ** \param[in]  x        Array containing the vertices coordinates along X
 ** \param[in]  y        Array containing the vertices coordinates along Y
 ** \param[in]  zmin     Lower bound for the PolySet (or TEST)
 ** \param[in]  zmax     Upper bound for the PolySet (or TEST)
 **
 ** \remark  Polygons are closed (if necessary) when added
 **
 *****************************************************************************/
Polygons* polygon_add(Polygons *polygon,
                                      const VectorDouble &x,
                                      const VectorDouble &y,
                                      double zmin,
                                      double zmax)
{
  if (polygon == nullptr) return (polygon);
  PolySet polyset = PolySet();
  polyset.init(x, y, zmin, zmax);
  polygon->addPolySet(polyset);
  return (polygon);
}

/****************************************************************************/
/*!
 **  Print a Polygons structure
 **
 ** \param[in]  polygon     Polygons structure to be printed
 ** \param[in]  flag_print  Level of detail
 ** \li                      0 : Only the number of polysets
 ** \li                      1 : The number of vertices per polyset
 ** \li                      2 : The vertices for each polyset
 **
 *****************************************************************************/
void polygon_print(Polygons *polygon, int flag_print)

{
  if (polygon == nullptr) return;
  polygon->display(flag_print);
}

/****************************************************************************/
/*!
 **  Check if one point belongs to a vertical interval of a (limited) polyset
 **
 ** \return  1 if the point belongs to the polygon; 0 otherwise
 **
 ** \param[in]  zz   array of point coordinates of the point along Z or TEST
 ** \param[in]  zmin Lower bound
 ** \param[in]  zmax Upper bound
 **
 *****************************************************************************/
static int st_polyset_inside_3D(double zz, double zmin, double zmax)
{
  if (FFFF(zz)) return (1);
  if (!FFFF(zmin) && zz < zmin) return (0);
  if (!FFFF(zmax) && zz > zmax) return (0);
  return (1);
}

/****************************************************************************/
/*!
 **  Check if one point belongs to a 2-D polyset
 **
 ** \return  1 if the point belongs to the polygon; 0 otherwise
 **
 ** \param[in]  xx  array of point coordinates of the point along X
 ** \param[in]  yy  array of point coordinates of the point along Y
 ** \param[in]  np  number of vertices of the polygon
 ** \param[in]  xp  array containing the polygon vertices along X
 ** \param[in]  yp  array containing the polygon vertices along Y
 **
 *****************************************************************************/
static int st_polyset_inside(double xx,
                             double yy,
                             int np,
                             const double *xp,
                             const double *yp)
{
  double dx, dy, xinter;
  int inter, j, sel;

  inter = 0;

  /* Loop on the polygon vertices */

  for (j = 0; j < np - 1; j++)
  {

    /* Horizontal segment */

    dy = yp[j + 1] - yp[j];
    if (dy == 0 && yy == yp[j])
    {
      if (xp[j + 1] > xp[j] && xx > xp[j] && xx < xp[j + 1])
      {
        inter = 1;
        continue;
      }
      if (xp[j + 1] < xp[j] && xx < xp[j] && xx > xp[j + 1])
      {
        inter = 1;
        continue;
      }
    }

    /* One vertex below and one vertex above: point distinct from segment */

    if (dy != 0 && ((yp[j] > yy && yp[j + 1] < yy)
        || (yp[j] < yy && yp[j + 1] > yy)))
    {
      dx = xp[j + 1] - xp[j];
      xinter = (dx * yy + dy * xp[j] - dx * yp[j]) / dy;
      if (xinter > xx) inter++;
    }

    /* One vertex below and one vertex above: point belongs to segment */

    if (dy != 0 && ((yp[j] > yy && yp[j + 1] < yy)
        || (yp[j] < yy && yp[j + 1] > yy)))
    {
      dx = xp[j + 1] - xp[j];
      xinter = (dx * yy + dy * xp[j] - dx * yp[j]) / dy;
      if (xinter == xx)
      {
        inter = 1;
        continue;
      }
    }

    /* Point is in contact with the highest vertex */

    if (yy == yp[j] && yp[j] > yp[j + 1] && xx < xp[j]) inter++;
    if (yy == yp[j + 1] && yp[j + 1] > yp[j] && xx < xp[j + 1]) inter++;

    /* Point coincides with a vertex */

    if (xx == xp[j] && yy == yp[j])
    {
      inter = 1;
      continue;
    }
  }
  sel = ((inter % 2) != 0);
  return (sel);
}

/****************************************************************************/
/*!
 **  Check if one point belongs to a Polygons
 **
 ** \return  1 if the point belongs to the Polygons; 0 otherwise
 **
 ** \param[in]  xx          coordinate of the point along X
 ** \param[in]  yy          coordinate of the point along Y
 ** \param[in]  zz          coordinate of the point along Z (or TEST)
 ** \param[in]  flag_nested Option for nested polysets (see details)
 ** \param[in]  polygon     Polygons structure
 **
 ** \remarks If flag_nested=TRUE, a sample is masked off if the number of
 ** \remarks polysets to which it belongs is odd
 ** \remarks If flag_nested=FALSE, a sample is masked off as soon as it
 ** \remarks belongs to one PolySet
 **
 *****************************************************************************/
int polygon_inside(double xx,
                                   double yy,
                                   double zz,
                                   int flag_nested,
                                   Polygons *polygon)
{
  if (flag_nested)
  {

    /* Loop on the polysets */

    int number = 0;
    for (int ipol = 0; ipol < polygon->getPolySetNumber(); ipol++)
    {
      const PolySet &polyset = polygon->getPolySet(ipol);
      if (st_polyset_inside(xx, yy, polyset.getNVertices(),
                            polyset.getX().data(), polyset.getY().data()))
        number++;
      if (number % 2 != 0 && st_polyset_inside_3D(zz, polyset.getZmin(),
                                                  polyset.getZmax()))
        return (1);
    }
  }
  else
  {
    for (int ipol = 0; ipol < polygon->getPolySetNumber(); ipol++)
    {
      const PolySet &polyset = polygon->getPolySet(ipol);
      if (st_polyset_inside(xx, yy, polyset.getNVertices(),
                            polyset.getX().data(), polyset.getY().data())
          && st_polyset_inside_3D(zz, polyset.getZmin(), polyset.getZmax()))
        return (1);
    }
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Returns the polygon extension
 **
 ** \param[in]  polygon     Polygons structure
 **
 ** \param[out] xmin        Minimum coordinate along X-axis
 ** \param[out] xmax        Maximum coordinate along X-axis
 ** \param[out] ymin        Minimum coordinate along Y-axis
 ** \param[out] ymax        Maximum coordinate along Y-axis
 **
 *****************************************************************************/
void polygon_extension(Polygons *polygon,
                                       double *xmin,
                                       double *xmax,
                                       double *ymin,
                                       double *ymax)
{
  polygon->getExtension(xmin, xmax, ymin, ymax);
}

/****************************************************************************/
/*!
 **  Returns the surface of the Polygons
 **
 ** \return  Surface
 **
 ** \param[in]  polygon  Polygons structure
 **
 *****************************************************************************/
double polygon_surface(Polygons *polygon)

{
  return polygon->getSurface();
}

/****************************************************************************/
/*!
 **  Ask the characteristics of the Polygon structure
 **
 ** \return  Pointer to the newly allocated Polygon
 **
 *****************************************************************************/
Polygons* input_polygon(void)

{
  Polygons *polygon;
  double xval, yval;
  int nvert, nquant, error, flag_stop;
  VectorDouble x, y;
  static int quant = 100;

  /* Initializations */

  error = 1;
  polygon = nullptr;

  /* Core allocation */

  polygon = polygon_create();
  if (polygon == nullptr) goto label_end;

  /* Implicit loop on the PolySets */

  label_loop: nquant = quant;
  x.resize(nquant);
  y.resize(nquant);

  /* Read the next polyset */

  nvert = flag_stop = 0;
  while (!flag_stop)
  {
    if (nvert > 0) flag_stop = !_lire_logical("Add another Vertex", 1, 1);

    if (flag_stop)
    {
      xval = x[0];
      yval = y[0];
    }
    else
    {
      xval = _lire_double("X-vertex", 1, 0, TEST, TEST);
      yval = _lire_double("Y-vertex", 1, 0, TEST, TEST);
    }

    if (nvert >= nquant)
    {
      nquant += quant;
      x.resize(nquant);
      y.resize(nquant);
    }

    x[nvert] = xval;
    y[nvert] = yval;
    nvert++;
  }

  /* Add the new PolySet */

  polygon = polygon_add(polygon, x, y, TEST, TEST);

  /* Continue with other PolySets */

  if (_lire_logical("Add another PolySet", 1, 0)) goto label_loop;

  /* Set the error returned code */

  error = 0;

  label_end: if (error) polygon = polygon_free(polygon);
  return (polygon);
}

/*****************************************************************************/
/*!
 **  Create a polygon from the convex hull of active samples
 **
 ** \return  Error returned code
 **
 ** \param[in]  db    descriptor of the Db serving for convex hull calculation
 **
 *****************************************************************************/
Polygons* polygon_hull(const Db *db)

{
  Polygons *polygons = nullptr;

  /* Preliminary check */

  if (db->getNDim() < 2)
  {
    messerr("The input Db must be contain at least 2 coordinates");
    return polygons;
  }
  int number = db->getSampleNumber();
  if (number <= 0)
  {
    messerr("No active data in the input Db. Convex Hull impossible");
    return polygons;
  }

  /* Core allocation */
  VectorInt index(number + 1);

  /* Calculate the center of gravity and the leftmost point */

  int rank = 0;
  int np = 0;
  double xg = 0.;
  double yg = 0.;
  for (int i = 0; i < number; i++)
  {
    if (!db->isActive(i)) continue;
    xg += db->getCoordinate(i, 0);
    yg += db->getCoordinate(i, 1);
    if (db->getCoordinate(i, 0) < db->getCoordinate(rank, 0)) rank = i;
  }
  xg /= number;
  yg /= number;
  index[0] = rank;
  np++;

  /* Implicit loop for finding the other points of the convex hull */

  for (;;)
  {
    double x2, x3, y2, y3;
    double x1 = db->getCoordinate(index[np - 1], 0);
    double y1 = db->getCoordinate(index[np - 1], 1);
    double x21 = xg - x1;
    double y21 = yg - y1;

    for (int i = 0; i < number; i++)
    {
      if (!db->isActive(i)) continue;
      if ((db->getCoordinate(i, 0) - x1) * y21 - (db->getCoordinate(i, 1) - y1)
          * x21
          <= 0.) continue;
      x21 = db->getCoordinate(i, 0) - x1;
      y21 = db->getCoordinate(i, 1) - y1;
      rank = i;
    }
    if (rank == index[0]) goto label_cont;

    /* Discard the previous point if on the line joining the current point */
    /* to the one before the previous one */

    if (np > 1)
    {
      x1 = db->getCoordinate(rank, 0);
      y1 = db->getCoordinate(rank, 1);
      x2 = db->getCoordinate(index[np - 1], 0);
      y2 = db->getCoordinate(index[np - 1], 1);
      x3 = db->getCoordinate(index[np - 2], 0);
      y3 = db->getCoordinate(index[np - 2], 1);
      if (ABS((x2-x3)*(y1-y3) - (x1-x3)*(y2-y3)) < EPSILON6) np--;
    }
    index[np] = rank;
    np++;
  }

  label_cont: index[np++] = index[0];
  index.resize(np);

  /* Create the polygons */

  polygons = polygon_create();
  if (polygons == nullptr) return polygons;

  VectorDouble x(np);
  VectorDouble y(np);
  for (int i = 0; i < np; i++)
  {
    x[i] = db->getCoordinate(index[i], 0);
    y[i] = db->getCoordinate(index[i], 1);
  }
  polygons = polygon_add(polygons, x, y, TEST, TEST);
  return polygons;
}

/*****************************************************************************/
/*!
 **  Fill the arrays 'x' and 'y' with the coordinates of the polygon
 **  obtained as the convex hull of the active samples of 'db'.
 **
 ** \return  Error returned code
 **
 ** \param[in]  db    descriptor of the Db serving for convex hull calculation
 **
 ** \param[out] x     Array of first coordinates
 ** \param[out] y     Array of second coordinates
 **
 *****************************************************************************/
int polygon_hull(const Db *db, VectorDouble &x, VectorDouble &y)
{
  /* Preliminary check */

  if (db->getNDim() < 2)
  {
    messerr("The input Db must be contain at least 2 coordinates");
    return 1;
  }
  int number = db->getSampleNumber();
  if (number <= 0)
  {
    messerr("No active data in the input Db. Convex Hull impossible");
    return 1;
  }

  /* Core allocation */
  VectorInt index(number + 1);

  /* Calculate the center of gravity and the leftmost point */

  int rank = 0;
  int np = 0;
  double xg = 0.;
  double yg = 0.;
  for (int i = 0; i < number; i++)
  {
    if (!db->isActive(i)) continue;
    xg += db->getCoordinate(i, 0);
    yg += db->getCoordinate(i, 1);
    if (db->getCoordinate(i, 0) < db->getCoordinate(rank, 0)) rank = i;
  }
  xg /= number;
  yg /= number;
  index[0] = rank;
  np++;

  /* Implicit loop for finding the other points of the convex hull */

  for (;;)
  {
    double x2, x3, y2, y3;
    double x1 = db->getCoordinate(index[np - 1], 0);
    double y1 = db->getCoordinate(index[np - 1], 1);
    double x21 = xg - x1;
    double y21 = yg - y1;

    for (int i = 0; i < number; i++)
    {
      if (!db->isActive(i)) continue;
      if ((db->getCoordinate(i, 0) - x1) * y21 - (db->getCoordinate(i, 1) - y1)
          * x21
          <= 0.) continue;
      x21 = db->getCoordinate(i, 0) - x1;
      y21 = db->getCoordinate(i, 1) - y1;
      rank = i;
    }
    if (rank == index[0]) goto label_cont;

    /* Discard the previous point if on the line joining the current point */
    /* to the one before the previous one */

    if (np > 1)
    {
      x1 = db->getCoordinate(rank, 0);
      y1 = db->getCoordinate(rank, 1);
      x2 = db->getCoordinate(index[np - 1], 0);
      y2 = db->getCoordinate(index[np - 1], 1);
      x3 = db->getCoordinate(index[np - 2], 0);
      y3 = db->getCoordinate(index[np - 2], 1);
      if (ABS((x2-x3)*(y1-y3) - (x1-x3)*(y2-y3)) < EPSILON6) np--;
    }
    index[np] = rank;
    np++;
  }

  label_cont: index[np++] = index[0];
  index.resize(np);

  x.resize(np);
  y.resize(np);
  for (int i = 0; i < np; i++)
  {
    x[i] = db->getCoordinate(index[i], 0);
    y[i] = db->getCoordinate(index[i], 1);
  }
  return 0;
}
