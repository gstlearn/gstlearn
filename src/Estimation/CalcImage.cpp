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
#include "Estimation/CalcImage.hpp"
#include "Basic/Convolution.hpp"
#include "Basic/Law.hpp"
#include "Basic/NamingConvention.hpp"
#include "Calculators/ACalcInterpolator.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Enum/EStatOption.hpp"
#include "Estimation/KrigingSystem.hpp"
#include "Model/Model.hpp"
#include "Morpho/Morpho.hpp"
#include "Neigh/NeighImage.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Stats/Classical.hpp"
#include "geoslib_f_private.h"
#include "geoslib_old_f.h"

namespace gstlrn
{
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
  if (!ACalcInterpolator::_check()) return false;

  if (!hasDbin()) return false;
  auto nvar = getDbin()->getNLoc(ELoc::Z);
  if (!getDbin()->isGrid())
  {
    messerr("This method requires the Db to be a Grid");
    return false;
  }

  if (_flagFilter)
  {
    const auto* model = dynamic_cast<const ModelCovList*>(getModel());
    if (model == nullptr)
    {
      messerr("Model should be a ModelCovList");
      return false;
    }
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

  auto nvar = _getNVar();
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
    _renameVariable(2, VectorString(), ELoc::Z, 1, _iattOut, String {_oper.getKey()}, _nvarMorpho);

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
  Id ndim = dblocal->getNDim();
  Id nech = dblocal->getNSample();

  // Get the indices of the center grid node
  VectorInt center = dblocal->getCenterIndices();

  VectorVectorInt ranks;
  VectorInt local(ndim);
  for (Id iech = 0; iech < nech; iech++)
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

bool CalcImage::_filterImage(DbGrid* dbgrid, const ModelCovList* model)
{
  VectorDouble means;
  if (model->getNDrift() == 0) means = model->getMeans();

  Id ndim = dbgrid->getNDim();
  auto nvar = _getNVar();

  const auto* neighI = dynamic_cast<const NeighImage*>(getNeigh());
  DbGrid* dblocal          = neighI->buildImageGrid(dbgrid, _seed);
  VectorVectorInt ranks    = _getActiveRanks(dblocal);

  Db* target = Db::createFromOnePoint(VectorDouble(ndim));
  auto iuid  = target->addColumnsByConstant(nvar);

  // We perform a Kriging of the center 'dbaux' in Unique Neighborhood
  NeighUnique* neighU = NeighUnique::create();
  KrigingSystem ksys(dblocal, target, model, neighU, getKrigopt());
  if (ksys.updKrigOptEstim(iuid, -1, -1, true)) return false;
  if (!ksys.isReady()) return false;
  if (ksys.estimate(0)) return false;
  MatrixDense wgt = ksys.getWeights();
  ksys.conclusion();

  // Cleaning
  delete target;
  delete neighU;

  // Perform the Sparse convolution
  Convolution conv(dbgrid);

  Id retcode = 0;
  if (!_flagFFT)
  {
    retcode = conv.ConvolveSparse(_iattOut, ranks, wgt, means, static_cast<Id>(_verbose));
  }
  else
  {
    DbGrid* marpat = _buildMarpat(neighI, ranks, wgt, static_cast<Id>(_verbose));
    retcode        = conv.ConvolveFFT(_iattOut, nvar, marpat, means);
    delete marpat;
  }
  return (retcode == 0);
}

/**
 * @brief Construct a regular DbGrid of correct dimension
 * where the weights for the different variables are stored
 *
 * @param neigh NeighImage description
 * @param ranks Vector of Vector of neighborhood ranks
 * @param wgt   Matrix of weights
 * @param optionVerbose Verbose option (0: silent; 1: statistics; 2: whole display)
 * @return DbGrid
 */
DbGrid* CalcImage::_buildMarpat(const NeighImage* neigh,
                                const VectorVectorInt& ranks,
                                const MatrixDense& wgt,
                                Id optionVerbose)
{
  Id nbneigh = static_cast<Id>(ranks.size());
  Id ndim    = static_cast<Id>(ranks[0].size());
  auto nvar   = wgt.getNCols();
  VectorInt nx(ndim);
  for (Id i = 0; i < ndim; i++)
    nx[i] = 2 * neigh->getImageRadius(i) + 1;

  // Create the relevant DbGrid
  DbGrid* dbgrid   = DbGrid::create(nx);
  Id iuid         = dbgrid->addColumnsByConstant(nvar * nvar, 0., "Weights", ELoc::Z);
  VectorInt center = dbgrid->getCenterIndices();

  // Loop on the valid weights
  VectorInt local(ndim);
  for (Id ineigh = 0; ineigh < nbneigh; ineigh++)
  {
    local = ranks[ineigh];
    VH::addInPlace(local, center);
    Id iadd = dbgrid->indiceToRank(local);

    // Load the weights as variables
    Id ecr = 0;
    for (Id ivar = 0; ivar < nvar; ivar++)
      for (Id jvar = 0; jvar < nvar; jvar++, ecr++)
        dbgrid->setArray(iadd, iuid + ecr, wgt.getValue(nbneigh * jvar + ineigh, ivar));
  }

  // Optional printout
  if (optionVerbose)
  {
    mestitle(1, "Convolution Pattern");
    if (optionVerbose == 1)
    {
      Table table = dbStatisticsMono(dbgrid, {"Weights*"},
                                     {EStatOption::MINI, EStatOption::MAXI});
      for (Id ivar = 0, irow = 0; ivar < nvar; ivar++)
        for (Id jvar = 0; jvar < nvar; jvar++, irow++)
          table.setRowName(irow, "Weight of Z" + std::to_string(jvar + 1) +
                                   " for Z*" + std::to_string(ivar + 1));
      table.display();
    }
    else
    {
      DbStringFormat* dbfmt = DbStringFormat::createFromFlags(false, false, false,
                                                              false, true, false,
                                                              {"Weights*"});
      dbgrid->display(dbfmt);
    }
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
    const auto* model = dynamic_cast<const ModelCovList*>(getModel());
    if (!_filterImage(dbgrid, model)) return false;
  }

  if (_flagMorpho)
  {
    if (db_morpho_calc(dbgrid, _iattOut, _oper, _vmin, _vmax, _option, _radius,
                       _distErode, _verbose)) return false;
  }

  if (_flagSmooth)
  {
    const auto* neighI = dynamic_cast<const NeighImage*>(getNeigh());
    _image_smoother(dbgrid, neighI, _smoothType, _smoothRange, _iattOut);
  }

  return true;
}

/****************************************************************************/
/*!
**  Smooth a regular grid
**
** \param[in]  dbgrid    input and output Db grid structure
** \param[in]  neigh     Neigh structure
** \param[in]  type      1 for Uniform; 2 for Gaussian
** \param[in]  range     Range (used for Gaussian only)
** \param[in]  iptr0     Storage address
**
** \remarks Limited to the monovariate case
**
*****************************************************************************/
void CalcImage::_image_smoother(DbGrid* dbgrid,
                                const NeighImage* neigh,
                                Id type,
                                double range,
                                Id iptr0)
{
  Id ndim  = dbgrid->getNDim();
  double r2 = (type == 1) ? 1. : range * range;

  /* Core allocation */

  VectorInt indg0(ndim);
  VectorInt indgl(ndim);
  VectorInt indn0(ndim);
  VectorInt indnl(ndim);

  /* Create the secondary grid for image processing */

  VectorInt nx(ndim);
  Id nech = 1;
  for (Id idim = 0; idim < ndim; idim++)
  {
    nx[idim] = 2 * neigh->getImageRadius(idim) + 1;
    nech *= nx[idim];
  }

  law_set_random_seed(12345);
  double seuil = 1. / neigh->getSkip();
  VectorDouble tab(nech);
  for (Id iech = 0; iech < nech; iech++)
    tab[iech] = (law_uniform(0., 1.) < seuil) ? 0. : TEST;

  DbGrid* dbaux =
    DbGrid::create(nx, dbgrid->getDXs(), dbgrid->getX0s(), dbgrid->getAngles(),
                   ELoadBy::COLUMN, tab, {"test"}, {String {ELoc::Z.getKey()}}, 1);

  Id nb_neigh = dbaux->getNSample(true);
  dbaux->rankToIndice(nb_neigh / 2, indn0);

  /* Loop on the targets to be processed */

  for (Id iech_out = 0; iech_out < dbgrid->getNSample(); iech_out++)
  {
    if (!dbgrid->isActive(iech_out)) continue;
    dbgrid->rankToIndice(iech_out, indg0);

    /* Loop on the neighboring points */

    double estim = 0.;
    double total = 0.;
    for (Id iech = 0; iech < nb_neigh; iech++)
    {
      if (FFFF(dbaux->getZVariable(iech, 0))) continue;
      dbaux->rankToIndice(iech, indnl);
      double d2 = 0.;
      for (Id i = 0; i < ndim; i++)
      {
        Id idelta   = (indnl[i] - indn0[i]);
        double delta = idelta * dbgrid->getDX(i);
        d2 += delta * delta;
        indgl[i] = indg0[i] + idelta;
        indgl[i] = dbgrid->getMirrorIndex(i, indgl[i]);
      }

      Id jech    = dbgrid->indiceToRank(indgl);
      double data = dbgrid->getZVariable(jech, 0);
      if (!FFFF(data))
      {
        double weight = (type == 1) ? 1. : exp(-d2 / r2);
        estim += data * weight;
        total += weight;
      }
    }
    estim = (total <= 0.) ? TEST : estim / total;
    dbgrid->setArray(iech_out, iptr0, estim);
  }
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
 ** \param[in]  verbose    Verbose flag
 ** \param[in]  seed       Seed used for random number generation
 ** \param[in]  namconv    Naming Convention
 **
 *****************************************************************************/
Id krimage(DbGrid* dbgrid,
            Model* model,
            ANeigh* neigh,
            bool flagFFT,
            bool verbose,
            Id seed,
            const NamingConvention& namconv)
{
  CalcImage image;

  image.setDbin(dbgrid);
  image.setDbout(dbgrid);
  image.setModel(model);
  image.setNeigh(neigh);
  image.setFlagFFT(flagFFT);
  image.setSeed(seed);
  image.setVerbose(verbose);
  image.setNamingConvention(namconv);

  image.setFlagFilter(true);

  // Run the calculator
  Id error = (image.run()) ? 0 : 1;
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
Id dbSmoother(DbGrid* dbgrid,
               ANeigh* neigh,
               Id type,
               double range,
               const NamingConvention& namconv)
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
  Id error = (image.run()) ? 0 : 1;
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
GSTLEARN_EXPORT Id dbMorpho(DbGrid* dbgrid,
                             const EMorpho& oper,
                             double vmin,
                             double vmax,
                             Id option,
                             const VectorInt& radius,
                             bool flagDistErode,
                             bool verbose,
                             const NamingConvention& namconv)
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
  Id nvar = 1;
  if (oper == EMorpho::GRADIENT) nvar = dbgrid->getNDim();
  image.setNvarMorpho(nvar);

  // Run the calculator
  Id error = (image.run()) ? 0 : 1;
  return error;
}

} // namespace gstlrn