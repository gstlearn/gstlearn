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
#include "geoslib_f_private.h"

#include "Calculators/ACalcInterpolator.hpp"
#include "Estimation/CalcImage.hpp"
#include "Estimation/KrigingSystem.hpp"
#include "Neigh/NeighImage.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Db/DbGrid.hpp"
#include "Morpho/Morpho.hpp"
#include "Model/Model.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/OptCustom.hpp"
#include "Basic/Convolution.hpp"

#include "geoslib_old_f.h"

CalcImage::CalcImage()
  : ACalcInterpolator()
  , _iattOut(-1)
  , _flagFilter(false)
  , _flagFFT(false)
  , _seed(13242)
  , _flagMorpho(false)
  , _nvarMorpho(1)
  , _oper(EMorpho::UNKNOWN)
  , _vmin(0.5)
  , _vmax(1.5)
  , _option(0)
  , _radius()
  , _distErode(false)
  , _verbose(false)
  , _flagSmooth(false)
  , _smoothType(0)
  , _smoothRange(1.)
{
}

CalcImage::~CalcImage()
{
}

bool CalcImage::_check()
{
  if (! ACalcInterpolator::_check()) return false;

  if (!hasDbin()) return false;
  int nvar = getDbin()->getNLoc(ELoc::Z);
  if (! getDbin()->isGrid())
  {
    messerr("This method requires the Db to be a Grid");
    return false;
  }

  if (_flagFilter)
  {
    if (nvar <= 0)
    {
      messerr("This method requires some Variables to be defined in 'Db'");
      return false;
    }
  }

  if (_flagMorpho)
  {
    if (nvar != 1)
    {
      messerr("This method requires a single Variable to be defined in 'Db'");
      return false;
    }
  }

  if (_flagSmooth)
  {
    if (_smoothType != 1 && _smoothType != 2)
    {
      messerr("Filtering 'type' should be 1 or 2");
      return false;
    }
    if (nvar != 1)
    {
      messerr("This method requires a single Variable to be defined in 'Db'");
      return false;
    }
  }

  return true;
}

bool CalcImage::_preprocess()
{
  if (!ACalcInterpolator::_preprocess()) return false;

  int nvar = _getNVar();
  if (_flagFilter)
    _iattOut = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, nvar, 0.);

  if (_flagMorpho)
    _iattOut = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, _nvarMorpho, 0.);

  if (_flagSmooth)
    _iattOut = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, 1, 0.);

  if (_iattOut < 0) return false;
  return true;
}

bool CalcImage::_postprocess()
{
  /* Free the temporary variables */
  _cleanVariableDb(2);

  if (_flagFilter)
    _renameVariable(2, VectorString(), ELoc::Z, getDbin()->getNLoc(ELoc::Z), _iattOut, String(), 1);

  if (_flagMorpho)
    _renameVariable(2, VectorString(), ELoc::Z, 1, _iattOut, String{_oper.getKey()}, _nvarMorpho);

  if (_flagSmooth)
    _renameVariable(2, VectorString(), ELoc::Z, 1, _iattOut, String(), 1);

  return true;
}

void CalcImage::_rollback()
{
  _cleanVariableDb(1);
}

/**
 * @brief Create the vector of centered grid indices for neighborhood
 * 
 * @param dblocal Neighborhood Grid template
 * @return Returned vector of sample indices 
 */
VectorVectorInt CalcImage::_getActiveRanks(const DbGrid* dblocal)
{
  int ndim = dblocal->getNDim();
  int nech = dblocal->getNSample();

  // Get the indices of the center grid node
  VectorInt center = dblocal->getCenterIndices();

  VectorVectorInt ranks;
  VectorInt local(ndim);
  for (int iech = 0; iech < nech; iech++)
  {
    if (FFFF(dblocal->getZVariable(iech, 0))) continue;

    // The sample is valid, get its indices
    dblocal->rankToIndice(iech, local);

    // Center the indices
    VH::subtractInPlace(local, center);

    // Store these indices to the output vector
    ranks.push_back(local);
  }
  return ranks;
}

bool CalcImage::_filterImage(DbGrid* dbgrid, const ModelGeneric* modelgeneric)
{
  VectorDouble means;
  if (modelgeneric->getNDrift() == 0) means = modelgeneric->getMeans();

  int ndim = dbgrid->getNDim();
  int nvar = _getNVar();

  int optref = OptDbg::getReference();
  OptDbg::setReference(0);
  OptCst::define(ECst::NTROW, -1);

  const NeighImage* neighI = dynamic_cast<const NeighImage*>(getNeigh());
  DbGrid* dblocal          = neighI->buildImageGrid(dbgrid, _seed);
  VectorVectorInt ranks    = _getActiveRanks(dblocal);

  Db* target = Db::createFromOnePoint(VectorDouble(ndim));
  int iuid   = target->addColumnsByConstant(nvar);

  // We perform a Kriging of the center 'dbaux' in Unique Neighborhood
  NeighUnique* neighU = NeighUnique::create();
  KrigingSystem ksys(dblocal, target, modelgeneric, neighU);
  if (ksys.updKrigOptEstim(iuid, -1, -1, true)) return false;
  if (!ksys.isReady()) return false;
  if (ksys.estimate(0)) return false;
  MatrixRectangular wgt = ksys.getWeights();
  ksys.conclusion();

  // Cleaning
  delete target;
  delete neighU;
  OptDbg::setReference(optref);

  // Perform the Sparse convolution
  Convolution conv(dbgrid);

  int retcode = 0;
  if (!_flagFFT)
  {
    retcode = conv.ConvolveSparse(_iattOut, ranks, wgt, means);
  }
  else
  {
    DbGrid* marpat = _buildMarpat(neighI, ranks, wgt);
    retcode        = conv.ConvolveFFT(_iattOut, nvar, marpat, means);
    delete marpat;
  }

  return (retcode == 0);
}

bool CalcImage::_filterImageOld(DbGrid* dbgrid, const ModelGeneric* modelgeneric)
{
  message("on passe page le Krigeage\n");
  KrigingSystem ksys(dbgrid, dbgrid, modelgeneric, getNeigh());
  if (ksys.updKrigOptEstim(_iattOut, -1, -1)) return false;
  if (!ksys.isReady()) return false;

  /* Loop on the targets to be processed */

  for (int iech_out = 0; iech_out < dbgrid->getNSample(); iech_out++)
  {
    mes_process("Image filtering", dbgrid->getNSample(), iech_out);
    if (ksys.estimate(iech_out)) return false;
  }
  ksys.conclusion();
  return true;
}

  /**
   * @brief Construct a regular DbGrid of correct dimension
   * where the weights for the different variables are stored
   *
   * @param neigh NeighImage description
   * @param ranks Vector of Vector of neighborhood ranks
   * @param wgt   Matrix of weights
   * @return DbGrid
   */
  DbGrid* CalcImage::_buildMarpat(const NeighImage* neigh,
                                  const VectorVectorInt& ranks,
                                  const MatrixRectangular& wgt)
{
  int nbneigh = ranks.size();
  int ndim = ranks[0].size();
  int nvar    = wgt.getNCols();
  VectorInt nx(ndim);
  for (int i = 0; i < ndim; i++)
    nx[i] = 2 * neigh->getImageRadius(i)+ 1;

  // Create the relevant DbGrid
  DbGrid* dbgrid   = DbGrid::create(nx);
  int iuid         = dbgrid->addColumnsByConstant(nvar * nvar, 0., "Weights", ELoc::Z);
  VectorInt center = dbgrid->getCenterIndices();

  // Loop on the valid weights
  VectorInt local(ndim);
  for (int ineigh = 0; ineigh < nbneigh; ineigh++)
  {
    local = ranks[ineigh];
    VH::addInPlace(local,  center);
    int iadd = dbgrid->indiceToRank(local);

    // Load the weights as variables
    int ecr = 0;
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++, ecr++)
        dbgrid->setArray(iadd, iuid + ecr, wgt.getValue(nbneigh * jvar + ineigh, ivar));
  }
  return dbgrid;
}

/****************************************************************************/
/*!
 **  Standard Kriging
 **
 ** \return  Error return code
 **
 *****************************************************************************/
bool CalcImage::_run()
{
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbin());

  if (_flagFilter)
  {
    const ModelGeneric* modelgeneric = dynamic_cast<const ModelGeneric*>(getModel());

    bool _oldStyle = OptCustom::query("oldStyle", 1.) == 1.;
    if (!_oldStyle)
    {
      if (!_filterImage(dbgrid, modelgeneric)) return false;
    }
    else
    {
      if (!_filterImageOld(dbgrid, modelgeneric)) return false;
    }
  }

  if (_flagMorpho)
  {
    if (db_morpho_calc(dbgrid, _iattOut, _oper, _vmin, _vmax, _option, _radius,
                       _distErode, _verbose)) return false;
  }

  if (_flagSmooth)
  {
    const NeighImage* neighI = dynamic_cast<const NeighImage*>(getNeigh());
    _image_smoother(dbgrid, neighI, _smoothType, _smoothRange, _iattOut);
  }

  return true;
}

/****************************************************************************/
/*!
 **  Kriging (Factorial) a regular grid
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid     input and output Db grid structure
 ** \param[in]  model      Model structure
 ** \param[in]  neigh      ANeigh structure
 ** \param[in]  flagFFT    True if the FFT version is to be used
 ** \param[in]  seed       Seed used for random number generation
 ** \param[in]  namconv    Naming Convention
 **
 *****************************************************************************/
int krimage(DbGrid* dbgrid,
            Model* model,
            ANeigh* neigh,
            bool flagFFT,
            int seed,
            const NamingConvention& namconv)
{
  CalcImage image;

  image.setDbin(dbgrid);
  image.setDbout(dbgrid);
  image.setModel(model);
  image.setNeigh(neigh);
  image.setFlagFFT(flagFFT);
  image.setSeed(seed);
  image.setNamingConvention(namconv);

  image.setFlagFilter(true);

  // Run the calculator
  int error = (image.run()) ? 0 : 1;
  return error;
}

/****************************************************************************/
/*!
 **  Smooth a regular grid
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid     input and output Db grid structure
 ** \param[in]  neigh      ANeigh structure
 ** \param[in]  type       1 for Uniform; 2 for Gaussian
 ** \param[in]  range      Range (used for Gaussian only)
 ** \param[in]  namconv    Naming Convention
 **
 *****************************************************************************/
int dbSmoother(DbGrid *dbgrid,
               ANeigh *neigh,
               int type,
               double range,
               const NamingConvention &namconv)
{
  CalcImage image;

  image.setDbin(dbgrid);
  image.setDbout(dbgrid);
  image.setNeigh(neigh);
  image.setNamingConvention(namconv);

  image.setFlagSmooth(true);
  image.setSmoothType(type);
  image.setSmoothRange(range);

  // Run the calculator
  int error = (image.run()) ? 0 : 1;
  return error;
}

/**
 * Perform a Morphological operation on an image stored in Db
 * @param dbgrid  Target IN/OUT Db (must be a Grid)
 * @param oper    Type of morphological operation
 * @param vmin    Minimum threshold value
 * @param vmax    Maximum threshold value
 * @param option  Option
 * @param radius  Radius
 * @param verbose Verbose option
 * @param flagDistErode True: Inflate the grain; False: Reduce the grain
 * @param namconv Naming convention
 * @return
 */
GSTLEARN_EXPORT int dbMorpho(DbGrid *dbgrid,
                             const EMorpho &oper,
                             double vmin,
                             double vmax,
                             int option,
                             const VectorInt &radius,
                             bool flagDistErode,
                             bool verbose,
                             const NamingConvention &namconv)
{
  CalcImage image;

  image.setDbin(dbgrid);
  image.setDbout(dbgrid);
  image.setNamingConvention(namconv);

  image.setFlagMorpho(true);
  image.setOper(oper);
  image.setVmin(vmin);
  image.setVmax(vmax);
  image.setOption(option);
  image.setRadius(radius);
  image.setDistErode(flagDistErode);
  image.setVerbose(verbose);

  // Particular case of the number of output variables
  int nvar = 1;
  if (oper == EMorpho::GRADIENT) nvar = dbgrid->getNDim();
  image.setNvarMorpho(nvar);

  // Run the calculator
  int error = (image.run()) ? 0 : 1;
  return error;
}
