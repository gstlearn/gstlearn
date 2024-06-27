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
#include "Simulation/SimuSpectral.hpp"
#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Model/Model.hpp"

#include <math.h>

  SimuSpectral::  SimuSpectral(const Model* model)
    : _ndim(0),
      _nb(0),
      _isPrepared(false),
      _phi(),
      _gamma(),
      _omega(),
      _model(model)
{
}

  SimuSpectral::  SimuSpectral(const   SimuSpectral &r)
    : _ndim(r._ndim),
      _nb(r._nb),
      _isPrepared(r._isPrepared),
      _phi(r._phi),
      _gamma(r._gamma),
      _omega(r._omega),
      _model(r._model)
{
}

  SimuSpectral&   SimuSpectral::operator=(const   SimuSpectral &r)
{
  if (this != &r)
  {
    _ndim = r._ndim;
    _nb = r._nb;
    _isPrepared = r._isPrepared;
    _phi = r._phi;
    _gamma = r._gamma;
    _omega = r._omega;
    _model = r._model;
  }
  return *this;
}

SimuSpectral::~SimuSpectral()
{
}

int SimuSpectral::simulate(int nb, int seed)
{
  if (_model == nullptr)
  {
    messerr("A Model should be attached beforehand");
    return 1;
  }
  if (! isValidForSpectral(_model)) return 1;
  if (nb <= 0)
  {
    messerr("The number of spectral components should be positive");
    return 1;
  }

  _ndim = _model->getDimensionNumber();
  _nb = nb;

  law_set_random_seed(seed);

  _phi = VectorDouble(_nb);
  for (int ib = 0; ib < _nb; ib++)
    _phi[ib] = 2. * GV_PI * law_uniform();
  _gamma = VectorDouble(_nb);
  for (int ib = 0; ib < _nb; ib++)
    _gamma[ib] = sqrt(-log(law_uniform()));
  _omega = _model->getCova(0)->simulateSpectralOmega(_nb);

  _isPrepared = true;
  return 0;
}

int SimuSpectral::simulateOnSphere(int nb, int seed)
{
  if (_model == nullptr)
  {
    messerr("A Model should be attached beforehand");
    return 1;
  }
  if (! isValidForSpectral(_model)) return 1;
  if (nb <= 0)
  {
    messerr("The number of spectral components should be positive");
    return 1;
  }

  _ndim = _model->getDimensionNumber();
  _nb = nb;

  //  NK_  = simulate_spectrum_V2(ns = ns, spectrum = sp, seed = seed)

  law_set_random_seed(seed);

  _phi = VectorDouble(_nb);
  for (int ib = 0; ib < _nb; ib++)
    _phi[ib] = 2. * GV_PI * law_uniform();

  _isPrepared = true;
  return 0;
}

/**
 * Simulation of the spectral components (N,K) from spectrum values (version 2)
 *
 * @param verbose Verbose flag
 *
 * @return It returns the list with two vectors of length _nb
 * @return N contains the simulated degrees, K contains the simulated order -N <= K <= N)
 */
VectorVectorInt SimuSpectral::_simulateOnSphereV0()
{
  return VectorVectorInt();
}


VectorVectorInt SimuSpectral::_simulateOnSphere()
{
  return VectorVectorInt();
}

int SimuSpectral::compute(Db *dbout,
                          const VectorDouble &xref,
                          bool verbose,
                          const NamingConvention &namconv)
{
  int nech = dbout->getSampleNumber(true);
  int ndim = dbout->getNDim();

  // Preliminary checks

  if (ndim != _ndim)
  {
    messerr("The Space dimension of 'dbout'(%d) should match the one of Model(%d)",
            ndim, _ndim);
    return 1;
  }
  if (nech <= 0)
  {
    messerr("'dbout' must have a positive number of active samples");
    return 1;
  }
  if (! _isPrepared)
  {
    messerr("You should run 'simulate' beforehand");
    return 1;
  }

  // Preparation
  bool flagCenter = ((int) xref.size() == ndim);
  MatrixSquareGeneral tensor = _model->getCova(0)->getAniso().getTensorInverse();
  double scale = sqrt(2. / _nb);
  AMatrix* res = MatrixFactory::prodMatMat(&_omega, &tensor);

  // Create the variable

  int iuid = dbout->addColumnsByConstant(1, 0., String(), ELoc::Z);
  if (iuid < 0) return 1;

  // Optional printout
  if (verbose)
  {
    message("Spectral Simulation on a set of Isolated Points\n");
    message("- Number of samples = %d\n", nech);
    message("- Space dimension   = %d\n", ndim);
    message("- Number of spectral components = %d\n", _nb);
  }

  // Loop on the active samples
  VectorInt ranks = dbout->getRanksActive();
  VectorDouble coor(ndim);
  for (int jech = 0; jech < nech; jech++)
  {
    int iech = ranks[jech];
    dbout->getSampleCoordinatesInPlace(iech, coor);
    if (flagCenter) VH::subtractInPlace(coor, xref);
    VectorDouble u = res->prodMatVec(coor);

    double value = 0.;
    for (int ib = 0; ib < _nb; ib++)
      value += _gamma[ib] * cos(u[ib] + _phi[ib]);
    value *= scale;

    dbout->setArray(iech, iuid, value);
  }

  // Delete temporary matrix
  delete res;

  // Modify the name of the output
  namconv.setNamesAndLocators(dbout, iuid);
  return 0;
}

int SimuSpectral::computeOnSphere(Db *dbout,
                                  bool verbose,
                                  const NamingConvention &namconv)
{
  return 0;
}

/****************************************************************************/
/*!
 **  Check if the Model can be simulated using Spectral Method
 **
 ** \return  True if the Model is valid; 0 otherwise
 **
 ** \param[in]  model    Model structure
 **
 *****************************************************************************/
bool SimuSpectral::isValidForSpectral(const Model* model)
{
  ESpaceType type = getDefaultSpaceType();
  if (model->getCovaNumber() != 1)
  {
    messerr("This method only considers Model with a single covariance structure");
    return false;
  }

  /* Loop on the structures */

  for (int is = 0; is < model->getCovaNumber(); is++)
  {
    if (! model->getCova(is)->isValidForSpectral())
    {
      messerr("The current structure is not valid for Spectral Simulation");
      return false;
    }
    // Check that the model is valid for the current space type
  }
  return true;
}
