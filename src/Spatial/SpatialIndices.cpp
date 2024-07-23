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
#include "Db/DbStringFormat.hpp"
#include "Enum/ECalcVario.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Polygon/Polygons.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Space/SpacePoint.hpp"
#include "Variogram/VMap.hpp"

SpatialIndices::SpatialIndices(Db *db)
    : _db(db),
      _center(), _mvalues(), _mvectors(), 
      _inertia(TEST), _wztot(TEST), _iso(TEST), _nvalid(0), _theta(TEST),
      _ra(TEST), _rb(TEST)
{
}

SpatialIndices::SpatialIndices(const SpatialIndices &r)
    : AStringable(r), _center(r._center), _mvalues(r._mvalues),
      _mvectors(r._mvectors), _inertia(r._inertia), _wztot(r._wztot),
      _iso(r._iso), _nvalid(r._nvalid), _theta(r._theta), _ra(r._ra),
      _rb(r._rb)
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
 ** \return  True if the datum must be discarded
 **
 ** \param[in]  flag_w True if the weight is defined
 ** \param[in]  iech   Rank of the sample
 ** \param[in]  name   Name of the attribute (used if flag_z)
 **
 ** \param[out]  coor   Array of coordinates
 ** \param[out]  value  Target value
 ** \param[out]  weight Weight value
 ** \param[out]  wvalue Weighted value
 **
 *****************************************************************************/
bool SpatialIndices::_discardData(bool flag_w, int iech, const String &name,
                                  VectorDouble &coor, double *value,
                                  double *weight, double *wvalue) const
{
  // Check if the sample is masked off

  if (! _db->isActive(iech)) return 1;

  /* Check if the variable is defined */

  *value = 1.;
  if (! name.empty())
  {
    *value = _db->getValue(name, iech);
    if (FFFF(*value)) return true;
    if ((*value) < 0.)
    {
      messerr("The variable cannot be negative (Sample %d = %lf)", iech + 1, *value);
      messerr("Procedure is interrupted");
      return 1;
    }
  }

  /* Check if the weight is defined */

  *weight = 1.;
  if (flag_w)
  {
    *weight = _db->getWeight(iech);
    if (FFFF(*weight)) return true;
    if ((*weight) < 0.)
    {
      messerr("The weight cannot be negative (Sample %d = %lf)", iech + 1, *weight);
      messerr("Procedure is interrupted");
      return 1;
    }
  }

  /* Check if the sample has defined coordinates */

  _db->getCoordinatesPerSampleInPlace(iech, coor);
  for (int idim = 0, ndim = _db->getNDim(); idim < ndim; idim++)
    if (FFFF(coor[idim])) return true;

  /* Returning argument */

  *wvalue = (*value) * (*weight);

  return false;
}

/****************************************************************************/
/*!
 **  Calculate the Center of Gravity
 **
 ** \return  Error returned code
 **
 ** \param[in]  name Name of the optional attribute
 **
 *****************************************************************************/
int SpatialIndices::computeCGI(const String &name)
{
  // Initializations
  double wvalue, value, weight;
  int nech = _db->getSampleNumber();
  int ndim = _db->getNDim();
  bool flag_w = _db->hasLocVariable(ELoc::W);

  /* Calculate the Center of Gravity */

  _wztot = 0.;
  _nvalid = 0;
  _center.resize(ndim, 0.);
  VectorDouble coor(ndim, 0.);
  for (int iech = 0; iech < nech; iech++)
  {
    if (_discardData(flag_w, iech, name, coor, &value, &weight, &wvalue)) continue;
    for (int idim = 0; idim < ndim; idim++)
      _center[idim] += wvalue * coor[idim];
    _wztot += wvalue;
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
    if (_discardData(flag_w, iech, name, coor, &value, &weight, &wvalue)) continue;
    for (int idim = 0; idim < ndim; idim++)
      coor[idim] -= _center[idim];
    for (int idim = 0; idim < ndim; idim++)
    {
      _inertia += wvalue * coor[idim] * coor[idim];
      for (int jdim = 0; jdim <= idim; jdim++)
        mm.updValue(idim, jdim, EOperator::ADD, wvalue * coor[idim] * coor[jdim]);
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

double SpatialIndices::getLIC(const String &name1, const String &name2)
{
  // Initializations
  double wvalue, value1, value2, weight;
  int nech = _db->getSampleNumber();
  int ndim = _db->getNDim();
  bool flag_w = _db->hasLocVariable(ELoc::W);

  /* Calculate the Local Index of Collocation */

  int number = 0;
  double lic_z11 = 0.;
  double lic_z12 = 0.;
  double lic_z22 = 0.;
  VectorDouble coor(ndim, 0.);
  for (int iech = 0; iech < nech; iech++)
  {
    if (_discardData(flag_w, iech, name1, coor, &value1, &weight, &wvalue))
      continue;
    if (_discardData(flag_w, iech, name2, coor, &value2, &weight, &wvalue))
      continue;

    number++;
    lic_z12 += weight * value1 * value2;
    lic_z11 += weight * value1 * value1;
    lic_z22 += weight * value2 * value2;
  }

  if (number <= 0) return TEST;
  double lic = lic_z12 / sqrt(lic_z11 * lic_z22);
  return lic;
}

double SpatialIndices::getGIC(const String &name1, const String &name2)
{
  if (computeCGI(name1) != 0)
    return TEST;
  VectorDouble center1 = getCenter();
  double inertia1 = getInertia();
  if (computeCGI(name2) != 0)
    return TEST;
  VectorDouble center2 = getCenter();
  double inertia2 = getInertia();

  double dx = center1[0] - center2[0];
  double dy = center1[1] - center2[1];
  double d2 = dx * dx + dy * dy;
  double gic = 1. - d2 / (d2 + inertia1 + inertia2);
  return gic;
}

/*
 * Calculate the axes for representing the Inertia
 *
 * \return The vector containing successively: xl1, yl1, xl2, yl2
 */
VectorVectorDouble SpatialIndices::getAxes() const {
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

VectorDouble SpatialIndices::getAxe(int rank) const {
  VectorDouble vec;
  if (rank < 0 || rank > 3)
  {
    messerr("Argument 'rank' should lie between 0 and 3");
    return vec;
  }

  VectorVectorDouble axes = getAxes();
  return axes[rank];
}

MatrixRectangular SpatialIndices::getMatrixEllipse() const
{
  MatrixRectangular axes(4, 2);
  if (_mvalues.empty())
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
  if (_mvalues.empty())
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

/**
 * @brief Calculate several Spatial indices
 *
 * @param name Name of the Target variable
 *
 * \remark This functions have been developped in the scope of the UE
 * \remark program Fisboat, DG-Fish, STREP #502572
 */
void SpatialIndices::spatial(const String &name)
{
  double maille = 1.;
  if (_db->isGrid())
  {
    DbGrid *dbgrid = dynamic_cast<DbGrid *>(_db);
    maille = dbgrid->getCellSize();
  }
  double wvalue, value, weight;
  int ndim = _db->getNDim();
  bool flag_w = _db->hasLocVariable(ELoc::W);
  VectorDouble coor(ndim, 0.);

  /* Loop on the samples */

  double top = 0.;
  double bot = 0.;
  double sum = 0.;
  for (int iech = 0, nech = _db->getSampleNumber(); iech < nech; iech++)
  {
    if (_discardData(flag_w, iech, name, coor, &value, &weight, &wvalue)) continue;
    if (value > 0)
      sum += weight;
    top += weight * value;
    bot += weight * value * value;
  }
  top *= maille;
  bot *= maille;

  /* Returning arguments */

  double eqarea = (bot == 0.) ? TEST : top * top / bot;
  message("Abundance Index = %lf\n", top);
  message("Positive Area   = %lf\n", sum);
  message("Equivalent Area = %lf\n", eqarea);
}

String SpatialIndices::toString(const AStringFormat *strfmt) const {
  DECLARE_UNUSED(strfmt);
  std::stringstream sstr;

  sstr << toTitle(0, "Spatial Indices");

  if (!_center.empty())
    sstr << "Gravity Center" << toVectorDouble(_center) << std::endl;
  if (! FFFF(_inertia))
    sstr << "Inertia = " << _inertia << std::endl;
  if (! FFFF(_iso))
    sstr << "Isotropy = " << _iso << std::endl;
 
  return sstr.str();
}

VectorVectorDouble SpatialIndices::getQT(const String &name) const
{
  VectorVectorDouble vec;
  vec.resize(2);

  // Initializations
  double wvalue, value, weight;
  int ndim = _db->getNDim();
  bool flag_w = _db->hasLocVariable(ELoc::W);
  VectorDouble coor(ndim, 0.);

  /* Calculate the Center of Gravity */

  VectorDouble zz;
  VectorDouble ww;
  VectorDouble wz;
  for (int iech = 0, nech = _db->getSampleNumber(); iech < nech; iech++)
  {
    if (_discardData(flag_w, iech, name, coor, &value, &weight, &wvalue))
      continue;
    zz.push_back(value);
    ww.push_back(weight);
    wz.push_back(wvalue);
  }

  // Sort the data in ascending order

  VectorInt ranks  = VH::orderRanks(zz);
  VectorDouble zzs = VH::reorder(zz, ranks);
  VectorDouble wws = VH::reorder(ww, ranks);
  VectorDouble wzs = VH::reorder(wz, ranks);

  // Calculate the Spreading Area
  double Q = VH::cumul(wzs);
  double SA = 0.;
  VectorDouble QT = VH::cumsum(wzs, true);
  for (int ib = 0, nb = (int)ww.size(); ib < nb; ib++)
    SA += (QT[ib] + QT[ib + 1]) * wws[ib];
  SA /= Q;
  message("Spreading Area  = %lf\n", SA);

  // Calculate the T and Q(T) vectors
  vec[0] = VH::cumsum(wws, true, true);
  vec[1] = VH::cumsum(wzs, true, true);

  // Upgrade Q(T) into (Q-Q(T))/Q
  for (int i = 0, n = (int)vec[1].size(); i < n; i++)
    vec[1][i] = (Q - vec[1][i]) / Q;
  
  return vec;
}

double SpatialIndices::getMicroStructure(const String &name, double h0,
                                         const Polygons *polygon, double dlim,
                                         int ndisc)
{
  // Initializations
  double wvalue, value, weight;
  int ndim = _db->getNDim();
  bool flag_w = _db->hasLocVariable(ELoc::W);
  VectorDouble coor(ndim, 0.);

  // Calculate the Field extension 
  int number = 0;
  double xmin = +1.e30;
  double xmax = -1.e30;
  double ymin = +1.e30;
  double ymax = -1.e30;
  for (int iech = 0, nech = _db->getSampleNumber(); iech < nech; iech++)
  {
    if (_discardData(flag_w, iech, name, coor, &value, &weight, &wvalue))
      continue;
    number++;
    if (coor[0] < xmin)
      xmin = coor[0];
    if (coor[0] > xmax)
      xmax = coor[0];
    if (coor[1] < ymin)
      ymin = coor[1];
    if (coor[1] > ymax)
      ymax = coor[1];
  }

  if (number <= 0)
    return TEST;

  double dx = xmax - xmin;
  double dy = ymax - ymin;
  message("Initial dx=%lf dy=%lf\n", dx, dy);
  double extend = 2. * MAX(dlim / dx, dlim / dy);
  xmin -= dx * extend;
  xmax += dx * extend;
  ymin -= dy * extend;
  ymax += dy * extend;
  dx = (xmax - xmin) / (double) ndisc;
  dy = (ymax - ymin) / (double) ndisc;
  double maille = dx * dy;
  message("xmin=%lf xmax=%lf ymin=%lf ymax=%lf\n", xmin, xmax, ymin, ymax);
  message("dx=%lf dy=%lf nx=%d ny=%d\n", dx,dy,ndisc,ndisc);

  // Create the internal Grid
  DbGrid *grid = DbGrid::create({ndisc, ndisc}, {dx, dy}, {xmin, ymin});

  // Select grid nodes inside the polygon (if defined)
  if (polygon != nullptr)
    db_polygon(grid, polygon);

  // Migrate the Data to the Grid
  migrate(_db, grid, name);
  grid->display();


  // Prepare the variogram map calculation
  int nlag = ceil((3. * h0 / 2.) / MIN(dx, dy));
  nlag=3;
  int nrow = 2 * nlag + 1;
  int ncol = 2 * nlag + 1;
  message("nrow=%d ncol=%d dx=%lf dy=%lf maille=%lf\n", nrow, ncol, dx, dy, maille);
  DbGrid *vmap = db_vmap(grid, ECalcVario::E_COVARIANCE_NC, {nlag, nlag}, {dx, dy});

  // Calculate the Microstructure index
  int icenter = nrow * ncol / 2;
  double g0 = vmap->getValue("VMAP.Migrate.Var", icenter);

  // Delete the internal storage
  delete grid;
  delete vmap;

  return g0;
}

void tot(double a) {
  double b = 12.;
    }

static void
_updateGravityCenter(const VectorDouble &xxs, const VectorDouble &yys,
                     const VectorDouble &zzs, const VectorDouble &wws,
                     std::vector<SpacePoint> &centers, VectorDouble &pa,
                     VectorDouble &pb, VectorInt &ig, int found) {
  // Review the list of SpacePoints assigned to the current target group
  // to update the center of gravity

  double xt = 0.;
  double xb = 0.;
  double yt = 0.;
  double yb = 0.;
  double paval = 0.;
  double pbval = 0.;
  for (int iech = 0, nech = (int)ig.size(); iech < nech; iech++)
  {
    if (ig[iech] != found)
      continue;
    xt += wws[iech] * zzs[iech] * xxs[iech];
    xb += wws[iech] * zzs[iech];
    yt += wws[iech] * zzs[iech] * yys[iech];
    yb += wws[iech] * zzs[iech];
    paval += wws[iech];
    pbval += wws[iech] * zzs[iech];
  }

  // Update the registered center of gravity
  centers[found].setCoord(0, xt / xb);
  centers[found].setCoord(1, yt / yb);
  pa[found] = paval;
  pb[found] = pbval;
}

static SpacePoint _calculateGlobalGravityCenter(const VectorDouble &xxs,
                                                const VectorDouble &yys,
                                                const VectorDouble &zzs,
                                                const VectorDouble &wws,
                                                double *patot, double *pbtot)
{
  double xt = 0.;
  double xb = 0.;
  double yt = 0.;
  double yb = 0.;
  (*patot) = 0.;
  (*pbtot) = 0.;
  for (int iech = 0, nech = (int)xxs.size(); iech < nech; iech++) {
    xt += wws[iech] * zzs[iech] * xxs[iech];
    xb += wws[iech] * zzs[iech];
    yt += wws[iech] * zzs[iech] * yys[iech];
    yb += wws[iech] * zzs[iech];
    (*patot) += wws[iech];
    (*pbtot) += wws[iech] * zzs[iech];
  }

  // Global center of gravity
  SpacePoint centerG;
  centerG.setCoord(0, xt / xb);
  centerG.setCoord(1, yt / yb);
  return centerG;
}

static void _createNewPatch(int iech, const VectorDouble &xxs,
                            const VectorDouble &yys, const VectorDouble &zzs,
                            const VectorDouble &wws,
                            std::vector<SpacePoint> &centers, VectorDouble &pa,
                            VectorDouble &pb)
{
  SpacePoint current;
  current.setCoord(0, xxs[iech]);
  current.setCoord(1, yys[iech]);

  centers.push_back(current);
  pb.push_back(wws[iech] * zzs[iech]);
  pa.push_back(wws[iech]);
}

/**
* Returns the list of center of gravity of the different patches
* The last center of gravity is the global one
*
* @param name Name of the target variable
* @param Dmin Minimum rejection distance between patches
* @param Amin Abundance percentage above which patches are displayed
*/
std::vector<SpacePoint> SpatialIndices::getPatches(const String &name,
                                                   double Dmin,
                                                   double Amin)
{
  std::vector<SpacePoint> centers;
  VectorDouble pa;
  VectorDouble pb;

  // Initializations
  double wvalue, value, weight;
  int ndim = _db->getNDim();
  bool flag_w = _db->hasLocVariable(ELoc::W);
  VectorDouble coor(ndim, 0.);

  // Store the active samples in a list of Space Points
  VectorDouble xx;
  VectorDouble yy;
  VectorDouble ww;
  VectorDouble zz;
  VectorInt origRank;
  for (int iech = 0, nech = _db->getSampleNumber(); iech < nech; iech++)
  {
    if (_discardData(flag_w, iech, name, coor, &value, &weight, &wvalue))
      continue;
    xx.push_back(coor[0]);
    yy.push_back(coor[1]);
    zz.push_back(value);
    ww.push_back(weight);
    origRank.push_back(iech);
  }

  int nech = (int)origRank.size();
  if (nech <= 0)
    return centers;
  VH::normalize(ww, 1); // Normalize the weights

  // Sort the values by decreasing values
  VectorInt ranks = VH::orderRanks(zz, false);
  VectorDouble xxs = VH::reorder(xx, ranks);
  VectorDouble yys = VH::reorder(yy, ranks);
  VectorDouble zzs = VH::reorder(zz, ranks);
  VectorDouble wws = VH::reorder(ww, ranks);
  VectorInt origs = VH::reorder(origRank, ranks);
  VectorInt ig(nech, -1);

   // The first point is automatically assigned to the first center of gravity
  int iech = 0;
  ig[iech] = (int)centers.size();
  _createNewPatch(0, xxs, yys, zzs, wws, centers, pa, pb);

  // Loop on the samples
  SpacePoint current;
  for (int iech = 0; iech < nech; iech++)
  {
    current.setCoord(0, xxs[iech]);
    current.setCoord(1, yys[iech]);

    // Find which gravity center the current point aggregates to
    int found = -1;
    double dmin = 1.e30;
    for (int ic = 0, ncenter = (int)centers.size(); ic < ncenter; ic++)
    {
      double dist = current.getDistance(centers[ic]);
      if (dist > dmin)
        continue;
      found = ic;
      dmin = dist;
    }

    if (dmin < Dmin)
    {
      // The current sample has been assigned to group 'found': update gravity center
      ig[iech] = found; _updateGravityCenter(xxs, yys, zzs, wws, centers, pa, pb, ig, found);
    }
    else
    {
      // Create a new center of gravity
      ig[iech] = (int)centers.size();
      _createNewPatch(iech, xxs, yys, zzs, wws, centers, pa, pb);
    }
  }

  // Store the patch rank in the Data Base
  int iuid = _db->addColumnsByConstant(1, ITEST, "Patch");
  for (int i = 0, n = (int)origRank.size(); i < n; i++)
    _db->setArray(origs[i], iuid, ig[i]);

  // Calculate the global center of gravity

  double patot;
  double pbtot;
  SpacePoint globalCenter = _calculateGlobalGravityCenter(xxs, yys, zzs, wws, &patot, &pbtot);
  VH::divideConstant(pa, patot);
  VH::divideConstant(pb, pbtot);

  // Printout
  int ncenter = (int) centers.size();
  message("Regrouping the information of %s by patches\n", name.c_str());
  message("- Distance to Center of Gravity = %lf\n", Dmin);
  message("- Total Number of patches = %d\n", ncenter);

  int nover = 0;
  for (int ic = 0; ic < ncenter; ic++)
  {
    if (pb[ic] > Amin / 100) nover++;
  }
  message("- Number of patches with abundance > %3.0lf %% = %d\n", Amin, nover);
  message("- Percentage of abundance in these patches = ");
  for (int ic = 0; ic < ncenter; ic++) {
    if (pb[ic] > Amin / 100)
      message(" %6.3lf", 100. * pb[ic]);
  }
  message("\n");
  message("- Percentage of area in these patches = ");
  for (int ic = 0; ic < ncenter; ic++) {
    if (pb[ic] > Amin / 100)
      message(" %6.3lf", 100. * pa[ic]);
  }
  message("\n");

  centers.push_back(globalCenter);
  return centers;
}