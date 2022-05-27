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
#include "Basic/Plane.hpp"
#include "Basic/Law.hpp"
#include "Db/DbGrid.hpp"

#include <math.h>

Plane::Plane()
    : AStringable(),
      _coor(3),
      _intercept(0),
      _value(0.),
      _rndval(0.)
{
}

Plane::Plane(const Plane &m)
    : AStringable(m),
      _coor(m._coor),
      _intercept(m._intercept),
      _value(m._value),
      _rndval(m._rndval)
{
}

Plane& Plane::operator=(const Plane &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _coor = m._coor;
    _intercept = m._intercept;
    _value = m._value;
    _rndval = m._rndval;
  }
  return *this;
}

Plane::~Plane()
{

}

String Plane::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;


  return sstr.str();
}

/****************************************************************************/
/*!
 **  Generate the Poisson planes that cover the grid
 **
 ** \param[in]  dbgrid   Db corresponding to the target grid
 ** \param[in]  np       Number of planes
 **
 ** \remarks  The array 'planes' contains successively a,b,c,d such that
 ** \remarks  ax + by + cz + d = 0
 ** \remarks  The valuation of each line is assigned a uniform value [0,1]
 **
 *****************************************************************************/
std::vector<Plane> Plane::poissonPlanesGenerate(DbGrid *dbgrid, int np)
{
  double ap[3];

  VectorDouble center = dbgrid->getCenter();
  double diagonal = dbgrid->getExtensionDiagonal();
  std::vector<Plane> planes;
  planes.resize(np);

  /* Loop on the planes to be generated */

  for (int ip = 0; ip < np; ip++)
  {
    double d0 = diagonal * law_uniform(-1., 1.) / 2.;
    double u = 0.;
    for (int idim = 0; idim < 3; idim++)
    {
      ap[idim] = law_gaussian();
      u += ap[idim] * ap[idim];
    }
    u = sqrt(u);
    for (int idim = 0; idim < 3; idim++)
    {
      ap[idim] /= u;
      d0 -= ap[idim] * center[idim];
    }
    if (d0 < 0)
    {
      for (int idim = 0; idim < 3; idim++)
        ap[idim] = -ap[idim];
      d0 = -d0;
    }

    /* Storing the plane */

    for (int idim = 0; idim < 3; idim++)
      planes[ip].setCoor(idim, ap[idim]);
    planes[ip].setIntercept(d0);
    planes[ip].setRndval(law_uniform(0., 1.));
  }

  return planes;
}

