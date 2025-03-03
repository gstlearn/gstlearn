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

#include "Enum/ECst.hpp"
#include "Enum/ELoadBy.hpp"
#include "Enum/ESpaceType.hpp"

#include "API/SPDE.hpp"
#include "API/SPDEParam.hpp"

#include "Mesh/MeshETurbo.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Basic/Law.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/File.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Space/ASpaceObject.hpp"
#include "Covariances/CovContext.hpp"
#include "Model/Model.hpp"
#include "LinearOp/PrecisionOpMatrix.hpp"
#include "Basic/Memory.hpp"

/*****************************************************************************/
/*!
**  Save the resulting vectors in 'dbgrid' and save it in a neutral file
**
** \return Error returned code
**
** \param[in]  colors   Array giving the color
** \param[in]  consmin  Array giving the minimum bounds per pixel (optional)
** \param[in]  consmax  Array giving the maximum bounds per pixel (optional)
** \param[in]  z        Array giving the simulated field
**
*****************************************************************************/
static int st_save(Db    *dbgrid,
                   const VectorInt&    colors,
                   const VectorDouble& consmin,
                   const VectorDouble& consmax,
                   const VectorDouble& z)
{
  int iptr;
  int nech = dbgrid->getNSample();

  /* Add the terms to 'dbgrid' */
  
  if (! consmin.empty())
  {
    if (db_locator_attribute_add(dbgrid,ELoc::L,1,0,0.,&iptr)) return(1);
    for (int i=0; i<nech; i++) dbgrid->setArray(i,iptr,consmin[i]);
  }
  if (! consmax.empty())
  {
    if (db_locator_attribute_add(dbgrid,ELoc::U,1,0,0.,&iptr)) return(1);
    for (int i=0; i<nech; i++) dbgrid->setArray(i,iptr,consmax[i]);
  }
  iptr = dbgrid->addColumnsByConstant(1,0., "Color");
  for (int i=0; i<nech; i++) dbgrid->setArray(i,iptr,colors[i]);
  if (db_locator_attribute_add(dbgrid,ELoc::Z,1,0,0.,&iptr)) return(1);
  for (int i=0; i<nech; i++) dbgrid->setArray(i,iptr,z[i]);

  /* Save the resulting 'dbgrid' in a neutral file */

  (void) dbgrid->dumpToNF("Colored_Gibbs");
  return(0);
}

/*****************************************************************************/
/*!
**  Print all the imported information (used for check only)
**
** \param[in]  colors   Array of colors (same dimension as 'z')
** \param[in]  consmin  Array giving the minimum bounds per pixel (optional)
** \param[in]  consmax  Array giving the maximum bounds per pixel (optional)
** \param[in]  sigma    Array giving the st. dev. per pixel
** \param[in]  Q        Precision matrix (for all colors)
**
*****************************************************************************/
static void st_print_all(const VectorInt& colors,
                         const VectorDouble& consmin,
                         const VectorDouble& consmax,
                         const VectorDouble& sigma,
                         const MatrixSparse *Q)
{
  if (! consmin.empty())
    print_matrix ("consmin",0,0,1,10,NULL,consmin.data());
  if (! consmax.empty())
    print_matrix ("consmax",0,0,1,10,NULL,consmax.data());
  print_matrix ("sigma"  ,0,0,1,10,NULL,sigma.data());
  print_imatrix("colors" ,0,0,1,10,NULL,colors.data());
  Q->display();
}

/*****************************************************************************/
/*!
**  Compressing the input vector 'z' into 'zred' by selecting the only pixels
**  where the value of 'colors' is equal to 'colref'
**
** \return  Number of selected pixels
**
** \param[in]  nvertex  Dimension of the input array 'z'
** \param[in]  colref   Value of the reference color
** \param[in]  z        Array of input values
** \param[in]  colors   Array of colors (same dimension as 'z')
**
** \param[out]  ind     Array which contains address in 'z' of items of 'zred'
** \param[out]  zred    Array of output values retained
**
*****************************************************************************/
static void st_vector_compress(int nvertex,
                               int colref,
                               const VectorDouble &z,
                               const VectorInt &colors,
                               VectorInt &ind,
                               VectorDouble &zred)
{
  ind.clear();
  zred.clear();
  for (int i=0; i<nvertex; i++)
  {
    if (colors[i] != colref) 
      zred.push_back(z[i]);
    else
      ind.push_back(i);
  }
}
  
/*****************************************************************************/
/*!
**  Draw a random value from a Gaussian distribution under the constraint
**  to be larger to a minimum bound
**
** \return  Returned Gaussian value
**
** \param[in]  iter     Rank of the current iteration
** \param[in]  niter    Total count of iterations
** \param[in]  valmin   Value of the minimum bound
** \param[in]  valmax   Value of the maximum bound
** \param[in]  mean     Mean of the Gaussian distribution
** \param[in]  sigma    Standard deviation of the Gaussian distribution
**
*****************************************************************************/
static double st_simcond(int    iter,
                         int    niter,
                         double valmin,
                         double valmax,
                         double mean,
                         double sigma)
{
  double locmin,locmax,x,delta,ratio;

  ratio  = (double) (niter-iter-1) / (double) (iter+1);
  delta  = valmax - valmin;
  if (valmin == -10.)
    locmin = -10.;
  else 
  {    
    locmin = valmin - delta * ratio;
    locmin = (locmin - mean) / sigma;
  }
  if (valmax == +10.)
    locmax = valmax;
  else
  {
    locmax = valmax + delta * ratio;
    locmax = (locmax - mean) / sigma;
  }

  x = law_gaussian_between_bounds(locmin,locmax);
  return(mean + x * sigma);
}

/*****************************************************************************/
/*!
**  Perform the Gibbs sampler algorithm on the contents of the vector 'z'
**  This algorithm treats simultaneously all pixels sharing the same 'color'
**  (due to Markovian property, their influence do not overlap)
**
** \return  Error returned code
**
** \param[in]  niter    Maximum number of Gibbs iterations
** \param[in]  ncolor   Total number of different colors (in 'colors')
** \param[in]  nvertex  Dimension of the arrays 'z', 'colors', 'cons', 'sigma'
** \param[in]  colors   Array giving the colors per pixel
** \param[in]  colref   Array giving the unique colors
** \param[in]  Qcols    Array of 'cs' structures (sparse matrices per color)
** \param[in]  consmin  Array giving the minimum bounds per pixel (optional)
** \param[in]  consmax  Array giving the maximum bounds per pixel (optional)
** \param[in]  sigma    Array giving the st. dev. per pixel
**
** \param[in/out] z     Array of gaussian values simulated
** \param[out] krig     Working array
**
** \remarks Arrays 'colors' is integer
** \remarks        (Dimension 'nvertex')
** \remarks Arrays 'cons', 'sigma', 'z', 'krig' are double precision
** \remarks        (Dimension 'nvertex')
**
*****************************************************************************/
static int st_gibbs(int  niter,
                    int  ncolor,
                    int  nvertex,
                    const VectorInt& colors,
                    const VectorInt& colref,
                    MatrixSparse **Qcols,
                    const VectorDouble& consmin,
                    const VectorDouble& consmax,
                    const VectorDouble& sigma,
                    VectorDouble& z,
                    VectorDouble& krig)
{
  mestitle(1,"Entering in Gibbs algorithm with niter=%d and ncolor=%d",niter,ncolor);
  VectorInt ind;
  VectorDouble zred;

  for (int iter=0; iter<niter; iter++)
  {
    if (iter % 1000 == 0) message("Iteration %d\n",iter);
    for (int icol=0; icol<ncolor; icol++)
    {
      st_vector_compress(nvertex,colref[icol],z,colors,ind,zred);
      Qcols[icol]->prodVecMatInPlacePtr(zred.data(), krig.data(), false);

      for (int ic=0, nc = (int) ind.size(); ic<nc; ic++)
      {
        int i  = ind[ic];
        double valmin = (! consmin.empty()) ? consmin[i] : -10.;
        double valmax = (! consmax.empty()) ? consmax[i] : +10.;
        z[i] = st_simcond(iter, niter, valmin, valmax, krig[ic], sigma[i]);
      }
    }
  }
  return(0);
}

/****************************************************************************/
/*!
** Main Program
**
*****************************************************************************/
int main(int argc, char *argv[])

{
  bool   flag_print =  false;
  bool   flag_save  =  true;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  /***********************/
  /* 1 - Initializations */
  /***********************/
  int ndim     = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);
  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Gibbs-");

  // Setup constants

  OptDbg::reset();
  OptCst::define(ECst::NTCAR,10.);
  OptCst::define(ECst::NTDEC,6.);
  
  setGlobalFlagEigen(true);

  // 2-D grid output file

  VectorInt nx = { 100, 100 };
  VectorDouble x0 = { 0., 0. };
  VectorDouble dx = { 1., 1. };
  DbGrid *dbgrid = DbGrid::create(nx, dx, x0, VectorDouble(), ELoadBy::COLUMN,
                                  VectorDouble(), VectorString(), VectorString(), 1);
    
  // Model for SPDE

  double range_spde =   30.;
  double param_spde =    1.;
  double sill_spde  =    1.;
  Model* model1 = Model::createFromParam(ECov::MATERN,range_spde,sill_spde,param_spde);

  // Model for constraints

  double range_cons =   50.;
  double param_cons =    2.;
  double sill_cons  =    1.;
  Model* model2 = Model::createFromParam(ECov::MATERN,range_cons,sill_cons,param_cons);

  // Creating the meshing for extracting Q

  MeshETurbo mesh(dbgrid);
  auto P = PrecisionOpMatrix(&mesh, model1->getCovAniso(0));
  const MatrixSparse* Qref = P.getQ();
  MatrixSparse* Q = new MatrixSparse(*Qref);
  int nvertex = mesh.getNApices();

  // Coding the various colors

  VectorInt colors = Q->colorCoding();
  VectorInt colref = VH::unique(colors);
  int ncolor = (int) colref.size();

  // Core allocation
  
  VectorDouble z(nvertex,0);
  VectorDouble krig(nvertex,0);
  VectorDouble consmin(nvertex,0);
  VectorDouble consmax(nvertex,0);

  // Creating the constraints

  bool verbose = false;
  int seed = 31415;
  law_set_random_seed(seed);
  int nsimu = 2;
  int useCholesky = 1;
  (void)simulateSPDE(NULL, dbgrid, model2, nullptr, nsimu, NULL, useCholesky,
                     SPDEParam(), verbose);

  int rank = dbgrid->getNColumn();
  for (int i=0; i<nvertex; i++)
  {
    consmin[i] = MIN(dbgrid->getArray(i,rank-1), dbgrid->getArray(i,rank-2));
    consmax[i] = MAX(dbgrid->getArray(i,rank-1), dbgrid->getArray(i,rank-2));
    z[i] = (consmin[i] + consmax[i]) / 2.;
  }
  
  // Creating the variance
  
  VectorDouble sigma = Q->getDiagonal(0);
  VH::transformVD(sigma, -3);

  // Scaling the Q matrix
  
  (void) Q->scaleByDiag();

  // Check the imported information

  if (flag_print) st_print_all(colors,consmin,consmax,sigma,Q);

  //----------------//
  // Main Algorithm //
  //----------------//

  MatrixSparse** Qcols = (MatrixSparse **) mem_alloc(sizeof(MatrixSparse *) * ncolor,1);
  for (int icol=0; icol<ncolor; icol++)
    Qcols[icol] = Q->extractSubmatrixByColor(colors, colref[icol], true, false);

  // Perform the Gibbs sampler
  int niter = 10;
  (void) st_gibbs(niter, ncolor, nvertex, colors, colref, Qcols, consmin, consmax,
                  sigma, z, krig);

  // Add the newly created field to the grid for printout
  if (flag_save)
  {
    (void) st_save(dbgrid,colors,consmin,consmax,z);
    DbStringFormat dbfmt(FLAG_STATS);
    dbgrid->display(&dbfmt);
  }
  
  for (int icol=0; icol<ncolor; icol++)
    delete Qcols[icol];
  delete dbgrid;
  delete model1;
  delete model2;
  delete Q;
  
  return(0);
}
