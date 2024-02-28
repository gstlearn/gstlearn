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
#include "geoslib_define.h"
#include "geoslib_f_private.h"
#include "geoslib_f.h"

#include "Enum/EAnam.hpp"

#include "Variogram/VCloud.hpp"
#include "Variogram/Vario.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Model/Model.hpp"
#include "Variogram/VarioParam.hpp"
#include "Basic/Limits.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/OptDbg.hpp"
#include "Stats/Classical.hpp"
#include "Anamorphosis/AAnam.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Space/SpacePoint.hpp"
#include "Space/SpaceRN.hpp"
#include "Geometry/BiTargetCheckCode.hpp"
#include "Geometry/BiTargetCheckDate.hpp"
#include "Geometry/BiTargetCheckFaults.hpp"
#include "Geometry/BiTargetCheckGeometry.hpp"
#include "Morpho/Morpho.hpp"
#include "Polygon/Polygons.hpp"

static int IPTR;
static Polygons* POLYGON = nullptr;
static VectorDouble IDS;

VCloud::VCloud(DbGrid* dbcloud,
               const VarioParam* varioparam)
    : AVario(),
      _dbcloud(dbcloud),
      _varioparam(varioparam)
{
}

VCloud::VCloud(const VCloud& r)
    : AVario(r),
      _dbcloud(r._dbcloud),
      _varioparam(r._varioparam)
{
}

VCloud& VCloud::operator=(const VCloud& r)
{
  if (this != &r)
  {
    AVario::operator=(r);
    _dbcloud = r._dbcloud;
    _varioparam = r._varioparam;
  }
  return *this;
}

VCloud::~VCloud()
{
}

double VCloud::_getIVAR(const Db *db, int iech, int ivar) const
{
  return db->getLocVariable(ELoc::Z, iech, ivar);
}

/****************************************************************************/
/*!
 **  Internal function for setting a VCloud value
 **
 ** \param[in]  iech1       Rank of the first sample
 ** \param[in]  iech2       Rank of the second sample
 ** \param[in]  nvar        Number of variables
 ** \param[in]  ipas        Rank of the variogram lag
 ** \param[in]  ivar        Index of the first variable
 ** \param[in]  jvar        Index of the second variable
 ** \param[in]  orient      Orientation
 ** \param[in]  ww          Weight
 ** \param[in]  dist        Distance
 ** \param[in]  value       Variogram value
 **
 *****************************************************************************/
void VCloud::_setResult(int iech1,
                        int iech2,
                        int nvar,
                        int ipas,
                        int ivar,
                        int jvar,
                        int orient,
                        double ww,
                        double dist,
                        double value)
{
  DECLARE_UNUSED(nvar);
  DECLARE_UNUSED(ipas);
  DECLARE_UNUSED(ivar);
  DECLARE_UNUSED(jvar);
  DECLARE_UNUSED(orient);
  DECLARE_UNUSED(ww);

  int igrid = _update_discretization_grid(dist, value);
  if (igrid < 0) return;

  if (POLYGON == nullptr)
  {
    // Store in the output grid
    _dbcloud->updArray(igrid, IPTR, 0, 1.);
  }
  else
  {
    VectorInt indg(2);
    VectorDouble coor(2);
    db_index_sample_to_grid(_dbcloud, igrid, indg.data());
    grid_to_point(_dbcloud, indg.data(), NULL, coor.data());
    if (! POLYGON->inside(coor, false)) return;
    {
      IDS[iech1] += 1.;
      IDS[iech2] += 1.;
    }
  }
}

/****************************************************************************/
/*!
 **  Evaluate the experimental variogram cloud on irregular data
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  namconv      Naming convention
 **
 *****************************************************************************/
int VCloud::compute(Db *db, const NamingConvention &namconv)
{
  if (db == nullptr) return (1);

  /* Preliminary checks */

  if (db->getNDim() != _varioparam->getDimensionNumber())
  {
    messerr("Inconsistent parameters:");
    messerr("Data Base: NDIM=%d", db->getNDim());
    messerr("Variogram: NDIM=%d", _varioparam->getDimensionNumber());
    return (1);
  }
  if (!db->isVariableNumberComparedTo(1)) return 1;
  if (_dbcloud->getNDim() != 2)
  {
    messerr("The output Db for storing the variogram cloud must be 2-D");
    return (1);
  }

  /* Allocate new variables */

  setCalcul(ECalcVario::VARIOGRAM);
  int ndir = _varioparam->getDirectionNumber();
  int iptr = _dbcloud->addColumnsByConstant(ndir, 0.);
  if (iptr < 0) return (1);

  /* Loop on the directions to evaluate */

  for (int idir = 0; idir < ndir; idir++)
  {
    IPTR = iptr + idir;
    _variogram_cloud(db, idir);
    _final_discretization_grid();
  }

  // Naming of the newly created variables

  namconv.setNamesAndLocators(db, VectorString(), ELoc::Z, -1, _dbcloud, iptr,
                              String(), ndir, false);

  return (0);
}

/****************************************************************************/
/*!
 **  Replace zero values by TEST values
 **
 *****************************************************************************/
void VCloud::_final_discretization_grid()
{
  int nech = _dbcloud->getSampleNumber();
  for (int iech = 0; iech < nech; iech++)
  {
    double value = _dbcloud->getArray(iech, IPTR);
    if (value != 0.) continue;
    _dbcloud->setArray(iech, IPTR, TEST);
  }
}

/****************************************************************************/
/*!
 **  Evaluate the variogram cloud
 **
 ** \param[in]  db      Db descriptor
 ** \param[in]  idir    Rank of the Direction
 **
 *****************************************************************************/
void VCloud::_variogram_cloud(Db *db, int idir)
{
  double dist;
  SpaceTarget T1(_varioparam->getSpace());
  SpaceTarget T2(_varioparam->getSpace());

  // Creating a local Vario structure (to constitute the BiTargetCheck list
  Vario* vario = Vario::create(*_varioparam);
  vario->setDb(db);
  if (vario->prepare()) return;

  // Local variables to speed up calculations
  bool hasSel = db->hasLocVariable(ELoc::SEL);
  int nech = db->getSampleNumber();
  int nvar = db->getLocNumber(ELoc::Z);

  /* Loop on the first point */

  for (int iech = 0; iech < nech - 1; iech++)
  {
    if (hasSel && !db->isActive(iech)) continue;
    db->getSampleAsST(iech, T1);

    int ideb = (_varioparam->isDateUsed(db)) ? 0 : iech + 1;
    for (int jech = ideb; jech < nech; jech++)
    {
      if (hasSel && !db->isActive(jech)) continue;
      db->getSampleAsST(jech, T2);

      // Reject the point as soon as one BiTargetChecker is not correct
      if (! vario->keepPair(idir, T1, T2, &dist)) continue;

      evaluate(db, nvar, iech, jech, 0, dist, 0);
    }
  }

  delete vario;
  return;
}

/****************************************************************************/
/*!
 **  Add one pick to the discretization grid
 **
 ** \return  Index of the grid cell (or -1)
 **
 ** \param[in]  x     Coordinate along the first axis
 ** \param[in]  y     Coordinate along the first axis
 **
 *****************************************************************************/
int VCloud::_update_discretization_grid(double x, double y)
{
  int indg[2];

  int ix = (int) floor((x - _dbcloud->getX0(0)) / _dbcloud->getDX(0) + 0.5);
  int iy = (int) floor((y - _dbcloud->getX0(1)) / _dbcloud->getDX(1) + 0.5);
  if (ix < 0 || ix >= _dbcloud->getNX(0)) return (-1);
  if (iy < 0 || iy >= _dbcloud->getNX(1)) return (-1);
  indg[0] = ix;
  indg[1] = iy;
  return db_index_grid_to_sample(_dbcloud, indg);
}

/****************************************************************************/
/*!
 **  Evaluate the experimental variogram cloud
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  varioparam   VarioParam structure
 ** \param[in]  lagmax       Maximum distance
 ** \param[in]  varmax       Maximum Variance value (see remarks)
 ** \param[in]  lagnb        Number of discretization steps along distance axis
 ** \param[in]  varnb        Number of discretization steps along variance axis
 ** \param[in]  namconv      Naming convention
 **
 ** \remarks If 'varmax' is not defined, it is set to 3 times the experimental
 ** variance of the first variable (Z_locator)
 **
 *****************************************************************************/
DbGrid* db_vcloud(Db *db,
                  const VarioParam *varioparam,
                  double lagmax,
                  double varmax,
                  int lagnb,
                  int varnb,
                  const NamingConvention &namconv)
{
  if (FFFF(lagmax)) lagmax = db->getExtensionDiagonal();
  if (FFFF(varmax)) varmax = 3. * db->getVariance(db->getNameByLocator(ELoc::Z));

  // Create a grid as a support for the variogram cloud calculations

  VectorInt nx(2);
  nx[0] = lagnb;
  nx[1] = varnb;
  VectorDouble dx(2);
  dx[0] = lagmax / (double) lagnb;
  dx[1] = varmax / (double) varnb;
  VectorDouble x0(2);
  x0[0] = 0.;
  x0[1] = 0.;
  DbGrid *dbgrid = DbGrid::create(nx, dx, x0);

  // Calling the variogram cloud calculation function

  VCloud vcloud(dbgrid, varioparam);
  int error = vcloud.compute(db, namconv);

  // In case of error, free the newly created structure

  if (error)
  {
    delete dbgrid;
    dbgrid = nullptr;
  }
  return dbgrid;
}

/****************************************************************************/
/*!
 **  Check the samples which are involved in the pairs which are located
 **  within the polygon
 **
 ** \param[in]  db      Db descriptor
 ** \param[in]  polygon Polygons structure
 ** \param[in]  idir    Rank of the direction of itnerest
 **
 *****************************************************************************/
int VCloud::selectFromPolygon(Db *db, Polygons *polygon, int idir)
{
  POLYGON = polygon;
  int nech = db->getSampleNumber();
  IDS.resize(nech, 0.);

  _variogram_cloud(db, idir);

  /* Printout the scores: they are ranked by decreasing number */

  mestitle(0, "Samples in variogram cloud (by decreasing order of occurence)");
  VectorInt rank = VH::sequence(nech);
  ut_sort_double(0, nech, rank.data(), IDS.data());

  for (int iech = 0; iech < nech; iech++)
  {
    int jech = nech - iech - 1;
    if (IDS[jech] <= 0.) break;
    message("Sample #%3d: %d occurence(s)\n", rank[jech] + 1, (int) IDS[jech]);
  }

  POLYGON = nullptr;
  IDS.clear();
  return 0;
}

String VCloud::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  // Print the calculation type

  sstr << _elemString(strfmt) << std::endl;
  if (getCalcul() == ECalcVario::UNDEFINED) return sstr.str();

  // TODO: to be completed

  return sstr.str();
}
