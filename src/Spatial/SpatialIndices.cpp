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
#include "Spatial/SpatialIndices.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"

SpatialIndices::SpatialIndices(Db *db)
    : _db(db),
      _center(), _mvalues(), _mvectors(), 
      _inertia(TEST), _wztot(TEST), _iso(TEST), _nvalid(0), _theta(TEST),
      _ra(TEST), _rb(TEST), 
      _totab(TEST), _parea(TEST), _eqarea(TEST)
{
}

SpatialIndices::SpatialIndices(const SpatialIndices &r)
    : AStringable(r), _center(r._center), _mvalues(r._mvalues),
      _mvectors(r._mvectors), _inertia(r._inertia), _wztot(r._wztot),
      _iso(r._iso), _nvalid(r._nvalid), _theta(r._theta), _ra(r._ra),
      _rb(r._rb), 
      _totab(r._totab), _parea(r._parea), _eqarea(r._eqarea)
{
}

SpatialIndices &SpatialIndices::operator=(const SpatialIndices &r) {
  if (this != &r) {
    AStringable::operator=(r);
    _db = r._db;
    _center = r._center;
    _mvalues = r._mvalues;
    _mvectors = r._mvectors;
    _inertia = r._inertia;
    _wztot = r._wztot;
    _iso = r._iso;
    _nvalid = r._nvalid;
    _theta = r._theta;
    _ra = r._ra;
    _rb = r._rb;
    _totab = r._totab;
    _parea = r._parea;
    _eqarea = r._eqarea;
  }
  return *this;
}

SpatialIndices::~SpatialIndices()
{
  
}

/****************************************************************************/
/*!
 **  Load the valid data and calculate its weight
 **
 ** \return  Error returned code
 **
 ** \param[in]  flag_w True if the weight is defined
 ** \param[in]  iech   Rank of the sample
 ** \param[in]  name   Name of the attribute (used if flag_z)
 **
 ** \param[out]  coor   Array of coordinates
 ** \param[out]  wz     Weighted value
 **
 *****************************************************************************/
int SpatialIndices::_getData(bool flag_w, int iech, const String &name,
                             VectorDouble &coor, double *wz)
{
  // Check if the sample is masked off

  if (! _db->isActive(iech)) return 1;

  /* Check if the variable is defined */

  double value = 1.;
  if (! name.empty())
  {
    value = _db->getValue(name, iech);
    if (FFFF(value)) return (1);
    if (value < 0.)
    {
      messerr("The variable cannot be negative (Sample %d = %lf)", iech + 1, value);
      messerr("Procedure is interrupted");
      return 1;
    }
  }

  /* Check if the weight is defined */

  double weight = 1.;
  if (flag_w)
  {
    weight = _db->getWeight(iech);
    if (FFFF(weight)) return (1);
    if (weight < 0.)
    {
      messerr("The weight cannot be negative (Sample %d = %lf)", iech + 1, weight);
      messerr("Procedure is interrupted");
      return 1;
    }
  }

  /* Check if the sample has defined coordinates */

  _db->getCoordinatesPerSampleInPlace(iech, coor);
  for (int idim = 0, ndim = _db->getNDim(); idim < ndim; idim++)
    if (FFFF(coor[idim])) return 1;

  /* Returning argument */

  *wz = value * weight;

  return 0;
}

/****************************************************************************/
/*!
 **  Calculate the Center of Gravity
 **
 ** \return  Error returned code
 **
 ** \param[in]  db   Db structure
 ** \param[in]  name Name of the optional attribute
 **
 *****************************************************************************/
int SpatialIndices::computeCGI(Db *db, const String &name)
{
  // Initializations
  double wz;
  int nech = db->getSampleNumber();
  int ndim = db->getNDim();
  bool flag_w = db->hasLocVariable(ELoc::W);

  /* Calculate the Center of Gravity */

  _wztot = 0.;
  _nvalid = 0;
  _center.resize(ndim, 0.);
  VectorDouble coor(ndim, 0.);
  for (int iech = 0; iech < nech; iech++)
  {
    if (_getData(flag_w, iech, name, coor, &wz)) continue;
    for (int idim = 0; idim < ndim; idim++)
      _center[idim] += wz * coor[idim];
    _wztot += wz;
    _nvalid++;
  }
  if (_wztot <= 0.)
  {
    messerr("The sum of the weights must be positive : %lf", _wztot);
    return 1;
  }
  for (int idim = 0; idim < ndim; idim++)
    _center[idim] /= _wztot;

  /* Calculate the inertia and the weighted PCA */

  _inertia = 0.;
  MatrixSquareSymmetric mm(ndim);
  for (int iech = 0; iech < nech; iech++)
  {
    if (_getData(flag_w, iech, name, coor, &wz)) continue;
    for (int idim = 0; idim < ndim; idim++)
      coor[idim] -= _center[idim];
    for (int idim = 0; idim < ndim; idim++)
    {
      _inertia += wz * coor[idim] * coor[idim];
      for (int jdim = 0; jdim <= idim; jdim++)
        mm.updValue(idim, jdim, EOperator::ADD, wz * coor[idim] * coor[jdim]);
    }
  }

  /* Normation */
  _inertia /= _wztot;
  for (int idim = 0; idim < ndim; idim++)
    for (int jdim = 0; jdim <= idim; jdim++)
      mm.updValue(idim, jdim, EOperator::DIVIDE, _wztot);

  /* Calculate the eigen values and vectors */
  if (mm.computeEigen()) return 1;
  _mvalues = mm.getEigenValues();
  _mvectors = *mm.getEigenVectors();

  double r  = _mvalues[0] / _mvalues[1];
  double e2 = (_mvectors.getValue(1,0) / _mvectors.getValue(0,0));
  e2 = e2 * e2;
  _iso = 1. / sqrt(r);

  MatrixRectangular axes = getMatrixInertia();
  double dx1 = axes.getValue(1, 0) - axes.getValue(0, 0);
  double dy1 = axes.getValue(1, 1) - axes.getValue(0, 1);
  double dx2 = axes.getValue(3, 0) - axes.getValue(2, 0);
  double dy2 = axes.getValue(3, 1) - axes.getValue(2, 1);
  _theta = atan(dy1 / dx1) * 180. / GV_PI;
  _ra = 0.5 * sqrt(dx1*dx1 + dy1*dy1);
  _rb = 0.5 * sqrt(dx2*dx2 + dy2*dy2);

   return 0;
}

/*
* Calculate the axes for representing the Inertia
*
* \return The vector containing successively: xl1, yl1, xl2, yl2
*/
VectorVectorDouble SpatialIndices::getAxes() const
{
  VectorVectorDouble vec(4);
  if (_mvalues.empty())
  {
    messerr("You must use 'computeCGI() beforehand");
    return vec;
  }

  MatrixRectangular axes = getMatrixInertia();
  vec[0] = {axes.getValue(0, 0), axes.getValue(1, 0)};
  vec[1] = {axes.getValue(0, 1), axes.getValue(1, 1)};
  vec[2] = {axes.getValue(2, 0), axes.getValue(3, 0)};
  vec[3] = {axes.getValue(2, 1), axes.getValue(3, 1)};
  return vec;
}

MatrixRectangular SpatialIndices::getMatrixEllipse() const
{
  MatrixRectangular axes(4, 2);
  {
    messerr("You must use 'computeCGI() beforehand");
    return axes;
  }

  double r = _mvalues[0] / _mvalues[1];
  double K = sqrt(r) * (3. + r) / (2. * sqrt(2.) * pow(1. + r, 1.5));
  double e2 = (_mvectors.getValue(1, 0) / _mvectors.getValue(0, 0));
  e2 = e2 * e2;

  double sx1 = _mvectors.getValue(0, 0) / ABS(_mvectors.getValue(0, 0));
  double sy1 = _mvectors.getValue(1, 0) / ABS(_mvectors.getValue(1, 0));
  double sx2 = _mvectors.getValue(0, 1) / ABS(_mvectors.getValue(0, 1));
  double sy2 = _mvectors.getValue(1, 1) / ABS(_mvectors.getValue(1, 1));

  axes.setValue(0, 0, _center[0] + sx1 * sqrt(_inertia / (K * (1. + e2))));
  axes.setValue(0, 1, _center[1] + sy1 * sqrt(_inertia / (K * (1. + (1. / e2)))));
  axes.setValue(1, 0, 2. * _center[0] - axes.getValue(0, 0));
  axes.setValue(1, 1, 2. * _center[1] - axes.getValue(0, 1));
  axes.setValue(2, 0, _center[0] + sx2 * sqrt(_inertia / (K * r * (1. + (1. / e2)))));
  axes.setValue(2, 1, _center[1] + sy2 * sqrt(_inertia / (K * r * (1. + e2))));
  axes.setValue(3, 0, 2. * _center[0] - axes.getValue(2, 0));
  axes.setValue(3, 1, 2. * _center[1] - axes.getValue(2, 1));
  return axes;
}

MatrixRectangular SpatialIndices::getMatrixInertia() const
{
  MatrixRectangular axes(4, 2);
  {
    messerr("You must use 'computeCGI() beforehand");
    return axes;
  }

  double r1 = _mvalues[0] / (_mvalues[0] + _mvalues[1]);
  double e2 = (_mvectors.getValue(1, 0) / _mvectors.getValue(0, 0));
  e2 = e2 * e2;

  double sx1 = _mvectors.getValue(0, 0) / ABS(_mvectors.getValue(0, 0));
  double sy1 = _mvectors.getValue(1, 0) / ABS(_mvectors.getValue(1, 0));
  double sx2 = _mvectors.getValue(0, 1) / ABS(_mvectors.getValue(0, 1));
  double sy2 = _mvectors.getValue(1, 1) / ABS(_mvectors.getValue(1, 1));

  axes.setValue(0, 0, _center[0] + sx1 * sqrt((r1 * _inertia) / (1. + e2)));
  axes.setValue(0, 1, _center[1] + sy1 * sqrt((r1 * _inertia) / (1. + (1. / e2))));
  axes.setValue(1, 0, 2. * _center[0] - axes.getValue(0, 0));
  axes.setValue(1, 1, 2. * _center[1] - axes.getValue(0, 1));
  axes.setValue(2, 0, _center[0] + sx2 * sqrt(((1. - r1) * _inertia) / (1. + (1. / e2))));
  axes.setValue(2, 1, _center[1] + sy2 * sqrt(((1. - r1) * _inertia) / (1. + e2)));
  axes.setValue(3, 0, 2. * _center[0] - axes.getValue(2, 0));
  axes.setValue(3, 1, 2. * _center[1] - axes.getValue(2, 1));
  return axes;
}

/****************************************************************************/
/*!
 **  Calculate several Spatial indices
 **
 ** \param[in]  db   Db structure
 **
 ** \remark This functions have been developped in the scope of the UE
 ** \remark program Fisboat, DG-Fish, STREP #502572
 **
 *****************************************************************************/
void SpatialIndices::spatial(Db *db)
    {
  double top = 0.;
  double bot = 0.;
  double sum = 0.;
  double maille = 1.;
  if (db->isGrid())
  {
    DbGrid *dbgrid = dynamic_cast<DbGrid *>(db);
    maille = dbgrid->getCellSize();
  }

  /* Loop on the samples */

  for (int iech = 0; iech < db->getSampleNumber(); iech++) {
    if (!db->isActive(iech))
      continue;
    double z = db->getLocVariable(ELoc::Z, iech, 0);
    if (FFFF(z))
      continue;
    double w = db->getWeight(iech);
    if (FFFF(w))
      continue;
    if (z > 0)
      sum += w;
    top += w * z;
    bot += w * z * z;
  }
  top *= maille;
  bot *= maille;

  /* Returning arguments */

  _totab = top;
  _parea = sum;
  _eqarea = (bot == 0.) ? TEST : top * top / bot;
}

String SpatialIndices::toString(const AStringFormat *strfmt) const {
  DECLARE_UNUSED(strfmt);
  std::stringstream sstr;

  sstr << toTitle(0, "Spatiale Indices");

  sstr << "Gravity Center" << toVectorDouble(_center) << std::endl;
  sstr << "Number of active samples = " << _nvalid << std::endl;
  sstr << "Inertia = " << _inertia << std::endl;
  sstr << "Isotropy = " << _iso << std::endl;

  return sstr.str();
}
