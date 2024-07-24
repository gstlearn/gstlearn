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

#include "Variogram/VMap.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Model/Model.hpp"
#include "Variogram/VarioParam.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AStringable.hpp"
#include "Stats/Classical.hpp"
#include "Anamorphosis/AAnam.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Morpho/Morpho.hpp"

static int IPTV, IPTW;

#define ADD(ix,iy,iz,nx)    ((iz) + nx[2] * ((iy) + nx[1] * (ix)))
#define OPP(idim,i)         (dims[idim] - i - 1)

VMap::VMap(DbGrid* dbmap)
    : AVario(),
      _dbmap(dbmap)
{
}

VMap::VMap(const VMap& r)
    : AVario(r),
      _dbmap(r._dbmap)
{
}

VMap& VMap::operator=(const VMap& r)
{
  if (this != &r)
  {
    AVario::operator=(r);
    _dbmap = r._dbmap;
  }
  return *this;
}

VMap::~VMap()
{
}

double VMap::_getIVAR(const Db *db, int iech, int ivar) const
{
  return db->getLocVariable(ELoc::Z, iech, ivar);
}

/****************************************************************************/
/*!
 **  Internal function for setting a VMAP value
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
void VMap::_setResult(int iech1,
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
  DECLARE_UNUSED(iech1);
  DECLARE_UNUSED(iech2);
  DECLARE_UNUSED(orient);
  DECLARE_UNUSED(dist);
  int ijvar = _get_variable_order(nvar, ivar, jvar);
  _dbmap->updArray(ipas, IPTV + ijvar, EOperator::ADD, ww * value);
  _dbmap->updArray(ipas, IPTW + ijvar, EOperator::ADD, ww);
}

/****************************************************************************/
/*!
 **  Calculate the variogram map
 **
 ** \return  Error return code
 **
 ** \param[in]  db          Db containing the data
 ** \param[in]  calcul_type Type of calculation (ECalcVario)
 ** \param[in]  radius      Dilation radius (mooth resulting maps) only on points
 ** \param[in]  flag_FFT    Use FFT method (only valid on grid)
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int VMap::compute(Db *db,
                  const ECalcVario &calcul_type,
                  int radius,
                  bool flag_FFT,
                  const NamingConvention &namconv)
{
  if (db == nullptr) return 1;
  setCalcul(calcul_type);

  /* Create the variables in the Variogram Map file */

  int nvar = db->getLocNumber(ELoc::Z);
  int nvs2 = nvar * (nvar + 1) / 2;
  IPTV = _dbmap->addColumnsByConstant(nvs2, 0.);
  if (IPTV < 0) return 1;
  IPTW = _dbmap->addColumnsByConstant(nvs2, 0.);
  if (IPTW < 0) return 1;

  // Calculating the variogram map in different ways

  if (db->isGrid())
  {
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(db);

    // Case where Data are on a regular grid

    if (flag_FFT)
    {
      if (_grid_fft(dbgrid, namconv)) return 1;
    }
    else
      if (_vmap_grid(dbgrid, namconv)) return 1;
  }
  else
  {

    // Case where Data are on a set of points

    if (_vmap_general(db, radius, namconv)) return 1;
  }

  if (IPTW >= 0)
    namconv.setNamesAndLocators(db, VectorString(), ELoc::Z, -1,
                                _dbmap, IPTW, "Nb", 1, false);
  if (IPTV >= 0)
    namconv.setNamesAndLocators(db, VectorString(), ELoc::Z, -1,
                                _dbmap, IPTV, "Var");
  return 0;
}

/****************************************************************************/
/*!
 **  Calculate the variogram map (integrated function)
 **
 ** \return  Error return code
 **
 ** \param[in]  db          Db containing the data
 ** \param[in]  calcul_type Type of calculation (ECalcVario)
 ** \param[in]  nxx         Vector of (Half-) number of nodes for Vmap (def:20)
 ** \param[in]  dxx         Vector of mesh for Vmap (see details)
 ** \param[in]  radius      Dilation radius (mooth resulting maps) only on points
 ** \param[in]  flag_FFT    Use FFT method (only valid on grid)
 ** \param[in]  namconv     Naming convention
 **
 ** \remarks For calculating the default values:
 ** \remarks - for nx: it is set to 20 in all directions
 ** \remarks - for dx:
 ** \remarks   . If 'Db' is a grid, the mesh of the grid is used
 ** \remarks   - Otherwise, the mesh is set to the field extension / nx
 **
 *****************************************************************************/
DbGrid* db_vmap(Db *db,
                const ECalcVario &calcul_type,
                const VectorInt &nxx,
                const VectorDouble &dxx,
                int radius,
                bool flag_FFT,
                const NamingConvention &namconv)
{
  int error = 0;

  // Creating the output Variogram Map grid

  int ndim = db->getNDim();
  VectorInt nxloc = nxx;
  if (nxloc.empty()) nxloc.resize(ndim, 20);
  if (ndim != (int) nxloc.size())
  {
    messerr("Argument 'nxx' should have same Space Dimension as 'db'");
    return nullptr;
  }
  if (! dxx.empty() && ndim != (int) dxx.size())
  {
    messerr("Argument 'dxx'  should have same Space Dimension as 'db'");
    return nullptr;
  }
  VectorInt nx_map(ndim);
  VectorDouble dx_map(ndim);
  VectorDouble x0_map(ndim);

  for (int idim = 0; idim<ndim; idim++)
    nx_map[idim] = 2 * nxloc[idim] + 1;
  if (db->isGrid())
  {
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(db);
    for (int idim = 0; idim < ndim; idim++)
      dx_map[idim] = dbgrid->getDX(idim);
  }
  else
  {
    for (int idim = 0; idim < ndim; idim++)
      dx_map[idim] = (! dxx.empty() && !FFFF(dxx[idim])) ?
          dxx[idim] : db->getExtension(idim) / (double) nxloc[idim];
  }
  for (int idim = 0; idim < ndim; idim++)
    x0_map[idim] = -nxloc[idim] * dx_map[idim];

  DbGrid *dbmap = DbGrid::create(nx_map, dx_map, x0_map);

  // Calculating the variogram map in different ways

  VMap vmap(dbmap);
  error = vmap.compute(db, calcul_type, radius, flag_FFT, namconv);

  // In case of error, free the newly created VMAP structure

  if (error)
  {
    delete dbmap;
    dbmap = nullptr;
  }
  return dbmap;
}

 /*****************************************************************************/
 /*!
  **  Calculates the 2-D variogram map on a grid using FFT
  **
  ** \return  Error return code
  **
  ** \param[in]  dbgrid       Db of Grid type containing the data
  ** \param[in]  namconv      Naming convention
  **
  *****************************************************************************/
int VMap::_grid_fft(DbGrid *dbgrid, const NamingConvention &namconv)
{
  DECLARE_UNUSED(namconv);
   int dims[3], dinv[3], nxmap[3], nxgrid[3], sizemap, sizegrid;
   int ndim, nvar, ijvar;
   static bool verbose = false;
   VectorVectorDouble i1i1;
   VectorVectorDouble z1i1;
   VectorVectorDouble i2i2;
   VectorVectorDouble z2i2;
   VectorVectorDouble i1i2;
   VectorVectorDouble z1i2;
   VectorVectorDouble z2i1;
   VectorVectorDouble z2z1;
   VectorVectorDouble ztab;
   VectorDouble res_gg;
   VectorDouble res_nn;
   VectorDouble res_m1;
   VectorDouble res_m2;

   /* Initializations */

   int sizetot = 0;

   /* Preliminary checks */

   if (getCalcul() != ECalcVario::VARIOGRAM &&
       getCalcul() != ECalcVario::COVARIOGRAM &&
       getCalcul() != ECalcVario::COVARIANCE &&
       getCalcul() != ECalcVario::COVARIANCE_NC)
   {
     messerr("This function is limited to the calculation of");
     messerr("Variogram, Covariance (centered or not) or Covariogram");
     return (1);
   }
   if (dbgrid->getNDim() != 2 && dbgrid->getNDim() != 3)
   {
     messerr("The Variogram Map can only be calculated on a grid data set");
     messerr("with dimension equal to 2 or 3");
     return (1);
   }
   if (_dbmap->getNDim() > dbgrid->getNDim())
   {
     messerr("The space dimension of the VMAP (%d)", _dbmap->getNDim());
     messerr(
         "must not be larger than the space dimension of the input Grid (%d)",
         dbgrid->getNDim());
     return (1);
   }
   for (int idim = 0; idim < _dbmap->getNDim(); idim++)
   {
     if (ABS(_dbmap->getDX(idim) - dbgrid->getDX(idim)) > 1.e-03)
     {
       messerr("The grid mesh in the direction %d (dx=%lf)", idim + 1,
               dbgrid->getDX(idim));
       messerr("must match the mesh of the Variogram Map grid (dx=%lf)",
               _dbmap->getDX(idim));
       return (1);
     }
   }

   for (int idim = 0; idim < 3; idim++)
     nxgrid[idim] = nxmap[idim] = 1;
   for (int idim = 0; idim < dbgrid->getNDim(); idim++)
     nxgrid[idim] = dbgrid->getNX(idim);
   for (int idim = 0; idim < _dbmap->getNDim(); idim++)
     nxmap[idim] = _dbmap->getNX(idim);

   /* Preliminary calculations */

   nvar = dbgrid->getLocNumber(ELoc::Z);
   ndim = 0;
   sizetot = sizemap = sizegrid = 1;
   for (int i = 0; i < 3; i++)
   {
     dinv[i] = 1;
     if (nxgrid[i] <= 1)
     {
       dims[i] = 1;
     }
     else
     {
       dims[i] = (int) ceil((double) (nxgrid[i] + nxmap[i] - 1) / 8.) * 8;
       sizegrid *= nxgrid[i];
       sizemap *= nxmap[i];
       sizetot *= dims[i];
       ndim++;
     }
   }
   for (int i = 0; i < ndim; i++)
     dinv[i] = dims[ndim - i - 1];
   if (verbose)
   {
     mestitle(0, "Simulation of a grid using FFT");
     message("Grid dimension:");
     for (int idim = 0; idim < ndim; idim++)
       message(" %4d", nxgrid[idim]);
     message("\n");
     message("Variogram Map :");
     for (int idim = 0; idim < ndim; idim++)
       message(" %4d", nxmap[idim]);
     message("\n");
     message("Working array :");
     for (int idim = 0; idim < ndim; idim++)
       message(" %4d", dinv[idim]);
     message("\n");
   }

   /* Core allocation */

   _complexArrayAlloc(sizetot, ztab);
   if (getCalcul() == ECalcVario::VARIOGRAM)
   {
     _complexArrayAlloc(sizetot, i1i2);
     _complexArrayAlloc(sizetot, z1i2);
     _complexArrayAlloc(sizetot, z2i1);
     _complexArrayAlloc(sizetot, z2z1);
   }
   else
   {
     _complexArrayAlloc(sizetot, i1i1);
     _complexArrayAlloc(sizetot, z1i1);
     _complexArrayAlloc(sizetot, i2i2);
     _complexArrayAlloc(sizetot, z2i2);
   }

   res_nn.resize(sizemap, 0);
   res_gg.resize(sizemap, 0);
   if (getCalcul() == ECalcVario::COVARIANCE ||
       getCalcul() == ECalcVario::COVARIOGRAM)
   {
     res_m1.resize(sizemap, 0);
     res_m2.resize(sizemap, 0);
   }

   /* Loop on the variables */

   ijvar = 0;
   for (int ivar = 0; ivar < nvar; ivar++)
     for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
     {

       /* Calculate the structural function */

       if (getCalcul() == ECalcVario::VARIOGRAM)
       {
         if (_vmap_load_simple(dbgrid, ndim, sizetot, dims, dinv, ivar, jvar,
                               i1i2, z1i2, z2i1, z2z1)) continue;

         /* Calculate the number of pairs */
         _vmap_blank(ztab);
         _product_conjugate(1., i1i2, i1i2, ztab);
         if (fftn(ndim, dinv, ztab[0].data(), ztab[1].data(), -1, -1.)) continue;
         _extract(nxmap, nxgrid, dims, ztab[0], res_nn);

         /* Structural function */
         _vmap_blank(ztab);
         _product_conjugate( 1., z2z1, i1i2, ztab);
         _product_conjugate( 1., i1i2, z2z1, ztab);
         _product_conjugate(-1., z1i2, z2i1, ztab);
         _product_conjugate(-1., z2i1, z1i2, ztab);
         if (fftn(ndim, dinv, ztab[0].data(), ztab[1].data(), -1, -1.)) continue;
         _extract(nxmap, nxgrid, dims, ztab[0], res_gg);
         _vmap_rescale(2., res_gg, res_nn);
       }
       else
       {
         if (_vmap_load_cross(dbgrid, ndim, sizetot, dims, dinv, ivar, jvar,
                              i1i1, z1i1, i2i2, z2i2)) continue;

         /* Calculate the number of pairs */
         _vmap_blank(ztab);
         _product_conjugate(1., i1i1, i2i2, ztab);
         if (fftn(ndim, dinv, ztab[0].data(), ztab[1].data(), -1, -1.)) continue;
         _extract(nxmap, nxgrid, dims, ztab[0], res_nn);

         if (getCalcul() == ECalcVario::COVARIANCE ||
             getCalcul() == ECalcVario::COVARIOGRAM)
        {
          /* Calculate the means */
          _vmap_blank(ztab);
          _product_conjugate(1., z1i1, i2i2, ztab);
          if (fftn(ndim, dinv, ztab[0].data(), ztab[1].data(), -1, -1.)) continue;
          _extract(nxmap, nxgrid, dims, ztab[0], res_m1);
          _vmap_rescale(1., res_m1, res_nn);

          _vmap_blank(ztab);
          _product_conjugate(1., i1i1, z2i2, ztab);
          if (fftn(ndim, dinv, ztab[0].data(), ztab[1].data(), -1, -1.)) continue;
          _extract(nxmap, nxgrid, dims, ztab[0], res_m2);
          _vmap_rescale(1., res_m2, res_nn);
        }

        /* Structural function */
        _vmap_blank(ztab);
        _product_conjugate(1., z1i1, z2i2, ztab);
        if (fftn(ndim, dinv, ztab[0].data(), ztab[1].data(), -1, -1.)) continue;
        _extract(nxmap, nxgrid, dims, ztab[0], res_gg);
        _vmap_rescale(1., res_gg, res_nn);

        if (getCalcul() == ECalcVario::COVARIANCE ||
            getCalcul() == ECalcVario::COVARIOGRAM)
          _vmap_shift(res_gg, res_m1, res_m2);
       }

       /* Store the results */
       _vmap_store(res_nn, IPTW + ijvar);
       _vmap_store(res_gg, IPTV + ijvar);
     }

   return 0;
 }

/****************************************************************************/
/*!
 **  Load the resulting array in the output VMAP grid
 **
 ** \param[in]  nxmap   Dimensions of the Vmap Grid
 ** \param[in]  nxgrid  Dimensions of the Db Grid
 ** \param[in]  dims    Array of dimensions of the extended images
 ** \param[in]  tabin   Array to be stored (in FFT space)
 **
 ** \param[out] tabout  Array containing the resulting VMAP
 **
 *****************************************************************************/
void VMap::_extract(const int *nxmap,
                    const int *nxgrid,
                    int *dims,
                    VectorDouble& tabin,
                    VectorDouble& tabout)
{
  int ix, iy, iz, jx, jy, jz, nxs2, nys2, nzs2, nxloc, nyloc, nzloc;

  /* Initializations */

  nxs2 = nxmap[0] / 2;
  nys2 = nxmap[1] / 2;
  nzs2 = nxmap[2] / 2;
  nxloc = MIN(1 + nxs2, nxgrid[0]);
  nyloc = MIN(1 + nys2, nxgrid[1]);
  nzloc = MIN(1 + nzs2, nxgrid[2]);

  /* Fill the array (IX+,IY+,IZ+) */

  for (ix = 0; ix < nxloc; ix++)
    for (iy = 0; iy < nyloc; iy++)
      for (iz = 0; iz < nzloc; iz++)
      {
        jx = nxs2 + ix;
        jy = nys2 + iy;
        jz = nzs2 + iz;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(ix, iy, iz, dims)];
      }

  /* Fill the array (IX-,IY+,IZ+) */

  for (ix = 0; ix < nxloc - 1; ix++)
    for (iy = 0; iy < nyloc; iy++)
      for (iz = 0; iz < nzloc; iz++)
      {
        jx = nxs2 - ix - 1;
        jy = nys2 + iy;
        jz = nzs2 + iz;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(OPP(0,ix), iy, iz, dims)];
      }

  /* Fill the array (IX+,IY-,IZ+) */

  for (ix = 0; ix < nxloc; ix++)
    for (iy = 0; iy < nyloc - 1; iy++)
      for (iz = 0; iz < nzloc; iz++)
      {
        jx = nxs2 + ix;
        jy = nys2 - iy - 1;
        jz = nzs2 + iz;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(ix, OPP(1,iy), iz, dims)];
      }

  /* Fill the array (IX-,IY-,IZ+) */

  for (ix = 0; ix < nxloc - 1; ix++)
    for (iy = 0; iy < nyloc - 1; iy++)
      for (iz = 0; iz < nzloc; iz++)
      {
        jx = nxs2 - ix - 1;
        jy = nys2 - iy - 1;
        jz = nzs2 + iz;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(OPP(0,ix), OPP(1,iy), iz,dims)];
      }

  /* Fill the array (IX+,IY+,IZ-) */

  for (ix = 0; ix < nxloc; ix++)
    for (iy = 0; iy < nyloc; iy++)
      for (iz = 0; iz < nzloc - 1; iz++)
      {
        jx = nxs2 + ix;
        jy = nys2 + iy;
        jz = nzs2 - iz - 1;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(ix, iy, OPP(2,iz), dims)];
      }

  /* Fill the array (IX-,IY+,IZ-) */

  for (ix = 0; ix < nxloc - 1; ix++)
    for (iy = 0; iy < nyloc; iy++)
      for (iz = 0; iz < nzloc - 1; iz++)
      {
        jx = nxs2 - ix - 1;
        jy = nys2 + iy;
        jz = nzs2 - iz - 1;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(OPP(0,ix), iy, OPP(2,iz),dims)];
      }

  /* Fill the array (IX+,IY-,IZ-) */

  for (ix = 0; ix < nxloc; ix++)
    for (iy = 0; iy < nyloc - 1; iy++)
      for (iz = 0; iz < nzloc - 1; iz++)
      {
        jx = nxs2 + ix;
        jy = nys2 - iy - 1;
        jz = nzs2 - iz - 1;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(ix, OPP(1,iy), OPP(2,iz),dims)];
      }

  /* Fill the array (IX-,IY-,IZ-) */

  for (ix = 0; ix < nxloc - 1; ix++)
    for (iy = 0; iy < nyloc - 1; iy++)
      for (iz = 0; iz < nzloc - 1; iz++)
      {
        jx = nxs2 - ix - 1;
        jy = nys2 - iy - 1;
        jz = nzs2 - iz - 1;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(OPP(0,ix), OPP(1,iy),OPP(2,iz), dims)];
      }
}

/****************************************************************************/
/*!
 **  Calculate the variogram map when data are isolated points
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db containing the data
 ** \param[in]  radius       Dilation radius (used to smooth the resulting maps)
 ** \param[in]  namconv      Naming convention
 **
 *****************************************************************************/
int VMap::_vmap_general(Db *db, int radius, const NamingConvention &namconv)
{
  DECLARE_UNUSED(namconv);
  int error, nvar, nv2, i, idim, flag_out, nbmax;
  int *indg0, *indg1, *ind1, iech0, iech1, iech2, jech1, jech2, nech, ndim;
  double *delta, *mid, *coor, x0;
  VectorInt neigh;

  /* Preliminary checks */

  error = 0;
  indg0 = indg1 = ind1 = nullptr;
  delta = coor = mid = nullptr;

  if (db->getNDim() != 2 && db->getNDim() != 3)
  {
    messerr("The Variogram Map can only be calculated on a grid data set");
    messerr("with dimension equal to 2 or 3");
    return (1);
  }
  if (_dbmap->getNDim() > db->getNDim())
  {
    messerr("The space dimension of the VMAP (%d)", _dbmap->getNDim());
    messerr(
        "must not be larger than the space dimension of the input Grid (%d)",
        db->getNDim());
    return (1);
  }

  /* Initializations */

  ndim = _dbmap->getNDim();
  nvar = db->getLocNumber(ELoc::Z);
  nech = db->getSampleNumber();
  nv2 = nvar * (nvar + 1) / 2;

  /* Core allocation */

  indg0 = db_indg_alloc(_dbmap);
  if (indg0 == nullptr) goto label_end;
  indg1 = db_indg_alloc(_dbmap);
  if (indg1 == nullptr) goto label_end;
  ind1 = (int*) mem_alloc(sizeof(int) * nech, 0);
  if (ind1 == nullptr) goto label_end;
  delta = db_sample_alloc(db, ELoc::X);
  if (delta == nullptr) goto label_end;
  mid = db_sample_alloc(db, ELoc::X);
  if (mid == nullptr) goto label_end;
  coor = db_vector_alloc(db);
  if (coor == nullptr) goto label_end;

  /* Calculate a neighborhood (if radius > 0) */

  neigh = gridcell_neigh(ndim, 1, radius, false, false);
  nbmax = (int) neigh.size() / ndim;

  /* Calculate the VMAP half-extension */

  for (idim = 0; idim < ndim; idim++)
    mid[idim] = _dbmap->getNX(idim) * _dbmap->getDX(idim) / 2;

  /* Sorting the samples according to their first coordinate */

  if (db_coorvec_get(db, 0, coor)) goto label_end;
  for (i = 0; i < nech; i++)
    ind1[i] = i;
  ut_sort_double(1, nech, ind1, coor);

  /* Loop on the first data */

  for (jech1 = 0; jech1 < nech; jech1++)
  {
    iech1 = ind1[jech1];
    if (!db->isActive(iech1)) continue;
    x0 = db->getCoordinate(iech1, 0);

    /* Loop on the second data */

    for (jech2 = jech1; jech2 < nech; jech2++)
    {
      iech2 = ind1[jech2];
      if (!db->isActive(iech2)) continue;
      delta[0] = db->getCoordinate(iech2, 0) - x0;
      if (delta[0] > mid[0]) break;

      for (idim = 1, flag_out = 0; idim < ndim && flag_out == 0; idim++)
      {
        delta[idim] = db->getDistance1D(iech2, iech1, idim);
        if (delta[idim] > mid[idim]) flag_out = 1;
      }
      if (flag_out) continue;

      // Apply to the target cell
      if (point_to_grid(_dbmap, delta, 0, indg0)) continue;
      for (int in = 0; in < nbmax; in++)
      {
        iech0 = _findNeighCell(indg0, neigh.data(), in, indg1);
        if (iech0 < 0) continue;
        evaluate(db, nvar, iech1, iech2, iech0, TEST, 0);
      }

      // Avoid symmetry if point is compared to itself
      if (iech1 == iech2) continue;

      // Apply to the opposite target cell
      for (idim = 0; idim < ndim; idim++)
        delta[idim] = -delta[idim];
      if (point_to_grid(_dbmap, delta, 0, indg0)) continue;
      for (int in = 0; in < nbmax; in++)
      {
        iech0 = _findNeighCell(indg0, neigh.data(), in, indg1);
        if (iech0 < 0) continue;
        evaluate(db, nvar, iech1, iech2, iech0, TEST, 0);
      }
    }
  }

  /* Normalization */

  _vmap_normalize(nv2);

  /* Set the error return code */

  error = 0;

  label_end:
  db_indg_free(indg0);
  db_indg_free(indg1);
  mem_free((char* ) ind1);
  db_sample_free(delta);
  db_sample_free(mid);
  db_vector_free(coor);
  return (error);
}

/****************************************************************************/
/*!
 **  Calculate the variogram map when data are defined on a Grid
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid       Db of Grid type containing the data
 ** \param[in]  namconv      Naming Convention
 **
 *****************************************************************************/
int VMap::_vmap_grid(DbGrid *dbgrid, const NamingConvention &namconv)
{
  DECLARE_UNUSED(namconv);
  int error, nvar, nv2, idim, delta;
  int *ind1, *ind2, *ind0, iech0, iech1, iech2, flag_out, ndim;

  /* Preliminary checks */

  error = 0;
  ind0 = ind1 = ind2 = nullptr;
  if (dbgrid == nullptr) return (1);

  if (dbgrid->getNDim() != 2 && dbgrid->getNDim() != 3)
  {
    messerr("The Variogram Map can only be calculated on a grid data set");
    messerr("with dimension equal to 2 or 3");
    return (1);
  }
  if (_dbmap->getNDim() > dbgrid->getNDim())
  {
    messerr("The space dimension of the VMAP (%d)", _dbmap->getNDim());
    messerr(
        "must not be larger than the space dimension of the input Grid (%d)",
        dbgrid->getNDim());
    return (1);
  }
  for (idim = 0; idim < _dbmap->getNDim(); idim++)
  {
    if (ABS(_dbmap->getDX(idim) - dbgrid->getDX(idim)) > 1.e-03)
    {
      messerr("The grid mesh in the direction %d (dx=%lf)", idim + 1,
              dbgrid->getDX(idim));
      messerr("must match the mesh of the Variogram Map grid (dx=%lf)",
              _dbmap->getDX(idim));
      return (1);
    }
  }

  /* Initializations */

  ndim = _dbmap->getNDim();
  nvar = dbgrid->getLocNumber(ELoc::Z);
  nv2 = nvar * (nvar + 1) / 2;

  /* Core allocation */

  ind0 = db_indg_alloc(_dbmap);
  if (ind0 == nullptr) goto label_end;
  ind1 = db_indg_alloc(dbgrid);
  if (ind1 == nullptr) goto label_end;
  ind2 = db_indg_alloc(dbgrid);
  if (ind2 == nullptr) goto label_end;

  /* Loop on the first data */

  for (iech1 = 0; iech1 < dbgrid->getSampleNumber(); iech1++)
  {
    if (!dbgrid->isActive(iech1)) continue;
    db_index_sample_to_grid(dbgrid, iech1, ind1);

    /* Loop on the second data */

    for (iech2 = 0; iech2 < dbgrid->getSampleNumber(); iech2++)
    {
      if (!dbgrid->isActive(iech2)) continue;
      db_index_sample_to_grid(dbgrid, iech2, ind2);

      for (idim = flag_out = 0; idim < ndim && flag_out == 0; idim++)
      {
        delta = ind1[idim] - ind2[idim];
        int moitie = (_dbmap->getNX(idim) - 1) / 2;
        if (delta < -moitie || delta > moitie) flag_out = 1;
        ind0[idim] = delta + moitie;
      }
      if (flag_out) continue;

      /* Evaluate the variogram map */

      iech0 = db_index_grid_to_sample(_dbmap, ind0);
      evaluate(dbgrid, nvar, iech1, iech2, iech0, TEST, 0);
    }
  }

  /* Normalization */

  _vmap_normalize(nv2);

  /* Set the error return code */

  error = 0;

  label_end:
  db_indg_free(ind0);
  db_indg_free(ind1);
  db_indg_free(ind2);
  return (error);
}

/****************************************************************************/
/*!
 **  Get variable order
 **
 ** \return  Rank of the pair of variables (-1 if incorrect arguments)
 **
 ** \param[in]  nvar      Number of variables
 ** \param[in]  ivar0     Rank of the first variable
 ** \param[in]  jvar0     Rank of the second variable
 **
 *****************************************************************************/
int VMap::_get_variable_order(int nvar, int ivar0, int jvar0)
{
  int rank = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++, rank++)
    {
      if (ivar == ivar0 && jvar == jvar0) return rank;
      if (ivar == jvar0 && jvar == ivar0) return rank;
    }
  return -1;
}

/****************************************************************************/
/*!
 **  Allocate an array of complex values
 **
 ** \param[in]  size    Dimension of the complex array
 **
 ** \param[out] tab     Complex array to be allocated
 **
 *****************************************************************************/
void VMap::_complexArrayAlloc(int size, VectorVectorDouble& tab)
{
  tab.resize(2);
  for (int ic = 0; ic < 2; ic++)
    tab[ic].resize(size);
}

/****************************************************************************/
/*!
 **  Load the data into the FFT arrays
 **
 ** \return  Error returned code
 **
 ** \param[in] dbgrid    Db structure containing the input grid
 ** \param[in] ndim      Space dimension
 ** \param[in] sizetot   Dimension of the vectors
 ** \param[in] dims      Array of dimensions of the extended images
 ** \param[in] dinv      Array of dimensions of the extended images (inverted)
 ** \param[in] ivar      Rank of the first variable
 ** \param[in] jvar      Rank of the second variable
 **
 ** \param[out] i1i2    Pointer to the complex array 1(Z1) * 1(Z2)
 ** \param[out] z1i2    Pointer on the complex array Z1 * 1(Z2)
 ** \param[out] z2i1    Pointer on the complex array Z2 * 1(Z1)
 ** \param[out] z2z1    Pointer on the complex array Z1 * Z2
 **
 ** \remark The arrays are evaluated only if the input pointer is defined
 **
 *****************************************************************************/
int VMap::_vmap_load_simple(DbGrid* dbgrid,
                            int ndim,
                            int sizetot,
                            const int* dims,
                            int* dinv,
                            int ivar,
                            int jvar,
                            VectorVectorDouble& i1i2,
                            VectorVectorDouble& z1i2,
                            VectorVectorDouble& z2i1,
                            VectorVectorDouble& z2z1)
{
  DECLARE_UNUSED(sizetot);
  int ind1, ind2, indice[3];

  /* Initialize the complex array */

  if (! i1i2.empty()) _vmap_blank(i1i2);
  if (! z1i2.empty()) _vmap_blank(z1i2);
  if (! z2i1.empty()) _vmap_blank(z2i1);
  if (! z2z1.empty()) _vmap_blank(z2z1);

  /* Loop on the grid cells */

  int ecr = 0;
  for (int ix = 0; ix < dims[0]; ix++)
    for (int iy = 0; iy < dims[1]; iy++)
      for (int iz = 0; iz < dims[2]; iz++, ecr++)
      {
        if (ndim >= 1 && ix >= dbgrid->getNX(0)) continue;
        if (ndim >= 2 && iy >= dbgrid->getNX(1)) continue;
        if (ndim >= 3 && iz >= dbgrid->getNX(2)) continue;
        indice[0] = ix;
        indice[1] = iy;
        indice[2] = iz;
        int iech = db_index_grid_to_sample(dbgrid, indice);
        if (!dbgrid->getSelection(iech)) continue;
        double val1 = dbgrid->getLocVariable(ELoc::Z, iech, jvar);
        double val2 = dbgrid->getLocVariable(ELoc::Z, iech, ivar);
        ind1 = (!FFFF(val1));
        ind2 = (!FFFF(val2));

        if (! i1i2.empty()) i1i2[0][ecr] = (ind1 && ind2);
        if (! z1i2.empty()) z1i2[0][ecr] = (ind1 && ind2) ? val1 : 0.;
        if (! z2i1.empty()) z2i1[0][ecr] = (ind1 && ind2) ? val2 : 0.;
        if (! z2z1.empty()) z2z1[0][ecr] = (ind1 && ind2) ? val1 * val2 : 0.;
      }

  /* Perform the Fourier transform of the complex arrays defined */

  if (! i1i2.empty() && fftn(ndim, dinv, i1i2[0].data(), i1i2[1].data(), 1, 1.)) return (1);
  if (! z1i2.empty() && fftn(ndim, dinv, z1i2[0].data(), z1i2[1].data(), 1, 1.)) return (1);
  if (! z2i1.empty() && fftn(ndim, dinv, z2i1[0].data(), z2i1[1].data(), 1, 1.)) return (1);
  if (! z2z1.empty() && fftn(ndim, dinv, z2z1[0].data(), z2z1[1].data(), 1, 1.)) return (1);
  return (0);
}

/****************************************************************************/
/*!
 **  Load the data into the FFT arrays
 **
 ** \return  Error returned code
 **
 ** \param[in] dbgrid    Db structure containing the input grid
 ** \param[in] ndim      Space dimension
 ** \param[in] sizetot   Dimension of the vectors
 ** \param[in] dims      Array of dimensions of the extended images
 ** \param[in] dinv      Array of dimensions of the extended images (inverted)
 ** \param[in] ivar      Rank of the first variable
 ** \param[in] jvar      Rank of the second variable
 **
 ** \param[out] i1i1    Pointer to the complex array 1(Z1)
 ** \param[out] z1i1    Pointer to the complex array Z1
 ** \param[out] i2i2    Pointer to the complex array 1(Z2)
 ** \param[out] z2i2    Pointer to the complex array Z2
 **
 ** \remark The arrays are evaluated only if the input pointer is defined
 **
 *****************************************************************************/
int VMap::_vmap_load_cross(DbGrid *dbgrid,
                           int ndim,
                           int sizetot,
                           const int *dims,
                           int *dinv,
                           int ivar,
                           int jvar,
                           VectorVectorDouble &i1i1,
                           VectorVectorDouble &z1i1,
                           VectorVectorDouble &i2i2,
                           VectorVectorDouble &z2i2)
{
  DECLARE_UNUSED(sizetot);
  int ind1, ind2, indice[3];

  /* Initialize the complex array */

  if (! i1i1.empty()) _vmap_blank(i1i1);
  if (! z1i1.empty()) _vmap_blank(z1i1);
  if (! i2i2.empty()) _vmap_blank(i2i2);
  if (! z2i2.empty()) _vmap_blank(z2i2);

  /* Loop on the grid cells */

  int ecr = 0;
  for (int ix = 0; ix < dims[0]; ix++)
    for (int iy = 0; iy < dims[1]; iy++)
      for (int iz = 0; iz < dims[2]; iz++, ecr++)
      {
        if (ndim >= 1 && ix >= dbgrid->getNX(0)) continue;
        if (ndim >= 2 && iy >= dbgrid->getNX(1)) continue;
        if (ndim >= 3 && iz >= dbgrid->getNX(2)) continue;
        indice[0] = ix;
        indice[1] = iy;
        indice[2] = iz;
        int iech = db_index_grid_to_sample(dbgrid, indice);
        if (!dbgrid->getSelection(iech)) continue;
        double val1 = dbgrid->getLocVariable(ELoc::Z, iech, jvar);
        double val2 = dbgrid->getLocVariable(ELoc::Z, iech, ivar);
        ind1 = (!FFFF(val1));
        ind2 = (!FFFF(val2));

        if (! i1i1.empty()) i1i1[0][ecr] = (ind1);
        if (! z1i1.empty()) z1i1[0][ecr] = (ind1) ? val1 : 0.;
        if (! i2i2.empty()) i2i2[0][ecr] = (ind2);
        if (! z2i2.empty()) z2i2[0][ecr] = (ind2) ? val2 : 0.;
      }

  /* Perform the Fourier transform of the complex arrays defined */

  if (! i1i1.empty() && fftn(ndim, dinv, i1i1[0].data(), i1i1[1].data(), 1, 1.)) return (1);
  if (! z1i1.empty() && fftn(ndim, dinv, z1i1[0].data(), z1i1[1].data(), 1, 1.)) return (1);
  if (! i2i2.empty() && fftn(ndim, dinv, i2i2[0].data(), i2i2[1].data(), 1, 1.)) return (1);
  if (! z2i2.empty() && fftn(ndim, dinv, z2i2[0].data(), z2i2[1].data(), 1, 1.)) return (1);
  return (0);
}

/****************************************************************************/
/*!
 **  Blank the FFT arrays
 **
 ** \param[in] tab     Complex array to be blanked out
 **
 *****************************************************************************/
void VMap::_vmap_blank(VectorVectorDouble& tab)
{
  for (int ic = 0; ic < 2; ic++)
    for (int i = 0, size = (int) tab[ic].size(); i < size; i++)
      tab[ic][i] = 0.;
}

/****************************************************************************/
/*!
 **  Product of a vector by its conjugate
 **
 ** \param[in]  coef      Multiplicative coefficient
 ** \param[in]  tab1      First complex array
 ** \param[in]  tab2      Second complex array
 **
 ** \param[out]  tab      Output complex array
 **
 *****************************************************************************/
void VMap::_product_conjugate(double coef,
                              VectorVectorDouble& tab1,
                              VectorVectorDouble& tab2,
                              VectorVectorDouble& tab)
{
  for (int i = 0, size=(int) tab1[0].size(); i < size; i++)
  {
    tab[0][i] += coef * (tab1[0][i] * tab2[0][i] + tab1[1][i] * tab2[1][i]);
    tab[1][i] += coef * (tab1[0][i] * tab2[1][i] - tab1[1][i] * tab2[0][i]);
  }
}

/****************************************************************************/
/*!
 **  Scale the FFT arrays
 **
 ** \param[in] scale    Scaling factor
 ** \param[in,out] tab1 First complex array
 ** \param[in] tab2     Second complex array
 **
 *****************************************************************************/
void VMap::_vmap_rescale(double scale,
                         VectorDouble &tab1,
                         VectorDouble &tab2)
{
  for (int i = 0, size = (int) tab1.size(); i < size; i++)
  {
    double value = tab2[i];
    if (value > EPSILON8) tab1[i] /= (scale * value);
  }
}

/****************************************************************************/
/*!
 **  Shift the product of means for the FFT arrays
 **
 ** \param[in,out] tab Input/Output complex array
 ** \param[in] tabm1   Complex array for mean of variable 1
 ** \param[in] tabm2   Complex array for mean of variable 2
 **
 *****************************************************************************/
void VMap::_vmap_shift(VectorDouble &tab,
                       VectorDouble &tabm1,
                       VectorDouble &tabm2)
{
  for (int i = 0, size = (int) tab.size(); i < size; i++)
    tab[i] -= tabm1[i] * tabm2[i];
}

/****************************************************************************/
/*!
 **  Load the resulting array in the output VMAP grid
 **
 ** \param[in]  tab       Array to be stored (in FFT space)
 ** \param[in]  iptr      Pointer for storage
 **
 *****************************************************************************/
void VMap::_vmap_store(VectorDouble& tab, int iptr)
{
  int indice[3];
  VectorDouble dims(3);
  int ndim = _dbmap->getNDim();

  for (int idim = 0; idim < 3; idim++)
  {
    if (idim < ndim)
      dims[idim] = _dbmap->getNX(idim);
    else
      dims[idim] = 1;
  }

  /* Loop on the sample (supposedly ordered) */

  int ecr = 0;
  for (int ix = 0; ix < dims[0]; ix++)
    for (int iy = 0; iy < dims[1]; iy++)
      for (int iz = 0; iz < dims[2]; iz++, ecr++)
      {
        indice[0] = ix;
        indice[1] = iy;
        indice[2] = iz;
        int iech = db_index_grid_to_sample(_dbmap, indice);
        _dbmap->setArray(iech, iptr, tab[ecr]);
      }
}

/****************************************************************************/
/*!
 **  Get the absolute index of a grid node, shifted within a neighborhood
 **
 ** \return  Returned absolute index (<0 is not relevant)
 **
 ** \param[in]  indg0     Index decomposition for the central cell
 ** \param[in]  neigh     Array describing the neighborhooed
 **                       (Dimension: ndim * nbmax)
 ** \param[in]  rank      Rank of the cell within the neghborhood
 **
 ** \param[out] indg1     Working array for grid indices
 **
 *****************************************************************************/
int VMap::_findNeighCell(const int* indg0,
                         const int* neigh,
                         int rank,
                         int* indg1)
{
  int ndim;

  // Initializations

  ndim = _dbmap->getNDim();

  // Get the indices of the neighboring cell

  for (int idim = 0; idim < ndim; idim++)
    indg1[idim] = indg0[idim] + neigh[rank * ndim + idim];

  return (db_index_grid_to_sample(_dbmap, indg1));
}

/****************************************************************************/
/*!
 **  Scale the variogram maps
 **
 ** \param[in]  nv2       nvar ( nvar + 1) /2
 **
 *****************************************************************************/
void VMap::_vmap_normalize(int nv2)
{
  for (int iech = 0; iech < _dbmap->getSampleNumber(); iech++)
  {
    for (int ijvar = 0; ijvar < nv2; ijvar++)
    {
      double value = _dbmap->getArray(iech, IPTW + ijvar);
      if (value <= 0.)
        _dbmap->setArray(iech, IPTV + ijvar, TEST);
      else
        _dbmap->updArray(iech, IPTV + ijvar, EOperator::DIVIDE, value);
    }
  }
}

String VMap::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  // Print the calculation type

  sstr << _elemString(strfmt) << std::endl;
  if (getCalcul() == ECalcVario::UNDEFINED) return sstr.str();

  // TODO: to be completed

  return sstr.str();
}

