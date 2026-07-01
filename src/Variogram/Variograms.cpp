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
#include "Variogram/Variograms.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Space/SpaceRN.hpp"
#include "Variogram/VCloud.hpp"
#include "Variogram/VMap.hpp"
#include "Variogram/Vario.hpp"
#include "Variogram/VarioParam.hpp"

namespace gstlrn
{
  /**
   * @brief Returns a new Vario pointer initialized and calculated on the Z-variables of the 'db' database.
   *
   * @param db The 'db' database containing the Z-variables to calculate the variogram from.
   * @param calculType Calculation type (variogram, covariance, etc.)
   * @param flag_ergodic Ergodic flag
   * @param nlag Number of lags
   * @param dlag Lag value
   * @param ndir Number of directions
   * @param angles List of calculation angles
   * @param toldis Tolerance on distance
   * @param tolang Tolerance on angle
   * @param verbose Verbose flag
   * @return Vario* Pointer to the calculated Vario object (nullptr if an error occurred)
   *
   * @remarks The decision algorithm is described as follows:
   * - If 'ndir' > 1, multiple regular directions
   * - else if 'angles' is not empty, several 2D directions based on the provided angles
   * - else, omni-directional variogram
   */
  Vario* variogramCalculate(
    Db* db,
    const ECalcVario& calculType,
    bool flag_ergodic,
    Id nlag,
    double dlag,
    Id ndir,
    const VectorDouble& angles,
    double toldis,
    double tolang,
    bool verbose)
  {
    auto space = SpaceRN::create(db->getNDim());
    auto* varioparam =
      VarioParam::createForDb(ndir, nlag, dlag, angles, toldis, tolang, space);

    auto* vario = new Vario(*varioparam);
    if (vario->compute(
          db, calculType, flag_ergodic, false, false, nullptr, 0, verbose))
      return nullptr;
    return vario;
  }

  /**
   * @brief Calculates the variogram on a DbGrid object.
   *
   * @param dbgrid Target Db organized as a Grid
   * @param calculType Calculation type (variogram, covariance, etc.)
   * @param flag_ergodic Ergodic flag
   * @param nlag Number of lags (on each direction)
   * @param flagAllDirections if True, all directions are considered
   * @param dirincr This is the vector of direction increments to consider if flagAllDirections is False
   * @param verbose Verbose flag
   * @return Vario*
   */
  Vario* varioGridCalculate(
    DbGrid* dbgrid,
    const ECalcVario& calculType,
    bool flag_ergodic,
    Id nlag,
    bool flagAllDirections,
    const VectorVectorInt& dirincr,
    bool verbose)
  {
    auto space = SpaceRN::create(dbgrid->getNDim());
    auto* varioparam = VarioParam::createForGrid(
      dbgrid, flagAllDirections, nlag, dirincr, space);

    auto* vario = new Vario(*varioparam);
    if (vario->compute(
          dbgrid, calculType, flag_ergodic, false, false, nullptr, 0, verbose))
      return nullptr;
    return vario;
  }

  VectorInt variogramPerPoint(
    Db* db,
    Id iech0,
    Id ilag0,
    double dlag,
    const VectorDouble& codir,
    double toldis,
    double tolang,
    double bench,
    double cylrad,
    const VectorDouble& benchdir)
  {
    auto space = SpaceRN::create(db->getNDim());

    Id nlag = MAX(1, ilag0 + 1);
    Id opt_code = 0;
    double tolcode = 0.;
    auto dirparam = DirParam(
      nlag, dlag, toldis, tolang, opt_code, 0, bench, cylrad, tolcode,
      VectorDouble(), codir, TEST, benchdir, space);
    auto* varioparam = new VarioParam();
    varioparam->addDir(dirparam);

    auto* vario = new Vario(*varioparam);

    VectorInt pairs = vario->getPairs(db, iech0, ilag0);
    delete varioparam;
    delete vario;
    return pairs;
  }

  /****************************************************************************/
  /*!
   **  Calculate the variogram map (integrated function)
   **
   ** \return  Error return code
   **
   ** \param[in]  db          Db containing the data
   ** \param[in]  radius      Dilation radius (mooth resulting maps) only on points
   ** \param[in]  flag_FFT    Use FFT method (only valid on grid)
   ** \param[in]  calculType  Type of calculation (ECalcVario... only symmetrical ones)
   ** \param[in]  flag_ergodic Use ergodic assumption
   ** \param[in]  nxx         Vector of (Half-) number of nodes for Vmap (def:20)
   ** \param[in]  dxx         Vector of mesh for Vmap (see details)
   ** \param[in]  namconv     Naming convention
   **
   ** \remarks For calculating the default values:
   ** \remarks - for nx: it is set to 20 in all directions
   ** \remarks - for dx:
   ** \remarks   . If 'Db' is a grid, the mesh of the grid is used
   ** \remarks   - Otherwise, the mesh is set to the field extension / nx
   **
   *****************************************************************************/
  DbGrid* vmapFromDb(
    Db* db,
    Id radius,
    bool flag_FFT,
    const ECalcVario& calculType,
    bool flag_ergodic,
    const VectorInt& nxx,
    const VectorDouble& dxx,
    const NamingConvention& namconv)
  {
    if (db == nullptr)
    {
      messerr("You need a Db to compute a Variogram Map");
      return nullptr;
    }
    Id ndim = db->getNDim();
    VectorInt nxloc = nxx;
    if (nxloc.empty()) nxloc.resize(ndim, 20);
    if (ndim != static_cast<Id>(nxloc.size()))
    {
      messerr("Argument 'nxx' should have same Space Dimension as 'db'");
      return nullptr;
    }
    if (!dxx.empty() && ndim != static_cast<Id>(dxx.size()))
    {
      messerr("Argument 'dxx'  should have same Space Dimension as 'db'");
      return nullptr;
    }
    VectorInt nx_map(ndim);
    VectorDouble dx_map(ndim);
    VectorDouble x0_map(ndim);

    for (Id idim = 0; idim < ndim; idim++) nx_map[idim] = 2 * nxloc[idim] + 1;
    if (db->isGrid())
    {
      auto* dbgrid = dynamic_cast<DbGrid*>(db);
      for (Id idim = 0; idim < ndim; idim++) dx_map[idim] = dbgrid->getDX(idim);
    }
    else
    {
      for (Id idim = 0; idim < ndim; idim++)
        dx_map[idim] =
          (!dxx.empty() && !FFFF(dxx[idim]))
            ? dxx[idim]
            : db->getExtension(idim) / static_cast<double>(nxloc[idim]);
    }
    for (Id idim = 0; idim < ndim; idim++)
      x0_map[idim] = -nxloc[idim] * dx_map[idim];

    DbGrid* dbmap = DbGrid::create(nx_map, dx_map, x0_map);

    // Calculating the variogram map in different ways

    VMap vmap(dbmap);
    Id error =
      vmap.compute(db, radius, flag_FFT, calculType, flag_ergodic, namconv);

    // In case of error, free the newly created VMAP structure

    if (error)
    {
      delete dbmap;
      dbmap = nullptr;
    }
    return dbmap;
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
   ** \param[in]  calculType   Calculation type
   ** \param[in]  flag_ergodic Ergodic flag
   ** \param[in]  varmax       Maximum Variance value (see remarks)
   ** \param[in]  lagnb        Number of discretization steps along distance axis
   ** \param[in]  varnb        Number of discretization steps along variance axis
   ** \param[in]  polygon      Optional polygon to select the pairs of points
   ** \param[in]  distmax      Optional maximum distance to select the pairs of points
   ** \param[in]  varmin       Optional minimum variance to select the pairs of points
   ** \param[in]  countmax     Maximum number of outliers to be displayed (ITEST for no limit)
   ** \param[in]  samples      Optional list of samples to consider (otherwise all active samples are considered)
   ** \param[in]  namconv      Naming convention
   **
   ** \remarks If 'lagmax' is not provided, it is set to the diagonal of the
   ** area covered by the active samples within 'db'.
   ** \remarks If 'varmax' is not defined, it is set to 3 times the experimental
   ** variance of the first variable (Z_locator)
   **
   *****************************************************************************/
  DbGrid* vcloudFromDb(
    Db* db,
    const VarioParam* varioparam,
    double lagmax,
    const ECalcVario& calculType,
    bool flag_ergodic,
    double varmax,
    Id lagnb,
    Id varnb,
    const Polygons* polygon,
    double distmax,
    double varmin,
    Id countmax,
    const VectorInt& samples,
    const NamingConvention& namconv)
  {
    // Create a grid as a support for the variogram cloud calculations

    if (FFFF(lagmax)) lagmax = db->getExtensionDiagonal();
    if (FFFF(varmax))
      varmax = 3. * db->getVariance(db->getNameByLocator(ELoc::Z));

    VectorInt nx(2);
    nx[0] = lagnb;
    nx[1] = varnb;
    VectorDouble dx(2);
    dx[0] = lagmax / static_cast<double>(lagnb);
    dx[1] = varmax / static_cast<double>(varnb);
    VectorDouble x0(2);
    x0[0] = 0.;
    x0[1] = 0.;
    DbGrid* dbgrid = DbGrid::create(nx, dx, x0);

    // Create an omnidirectional varioparam (if not provided)

    if (varioparam == nullptr)
    {
      auto lagmax_local = lagmax;
      if (FFFF(lagmax_local)) lagmax_local = db->getExtensionDiagonal();
      // Evaluate the lag so that the maximum distance fit 'lagmax'
      // (given the tolerance on distance which is defaulted to 0.5)
      auto dlag = lagmax_local / 1.5;
      varioparam = VarioParam::createOmniDirection(1, dlag, 0.5);
    }

    // Initialize the variogram cloud structure

    VCloud vcloud(dbgrid, varioparam);

    // Additional selections

    vcloud.setPolygon(polygon);
    vcloud.setIntervals(distmax, varmin, countmax);
    vcloud.setSamples(samples);

    // Calling the variogram cloud calculation function

    if (vcloud.compute(db, calculType, flag_ergodic, namconv))
    {
      delete dbgrid;
      dbgrid = nullptr;
    }
    return dbgrid;
  }

  DbGrid* vcloudFromParam(
    Db* db,
    const ECalcVario& calculType,
    bool flag_ergodic,
    Id nlag,
    double dlag,
    Id ndir,
    const VectorDouble& angles,
    double toldis,
    double tolang,
    Id lagnb,
    Id varnb,
    const Polygons* polygon,
    double distmax,
    double varmin,
    Id countmax,
    const VectorInt& samples)
  {
    VarioParam* varioparam =
      VarioParam::createForDb(ndir, nlag, dlag, angles, toldis, tolang);

    auto* dbgrid = vcloudFromDb(
      db, varioparam, TEST, calculType, flag_ergodic, TEST, lagnb, varnb,
      polygon, distmax, varmin, countmax, samples);

    return dbgrid;
  }

} // namespace gstlrn
