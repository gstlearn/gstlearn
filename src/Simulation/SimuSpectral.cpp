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
#include "Stats/Classical.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Model/Model.hpp"

#include <math.h>

  SimuSpectral::  SimuSpectral(const Model* model)
    : _ndim(0),
      _nb(0),
      _nd(0),
      _ns(0),
      _isPrepared(false),
      _phi(),
      _gamma(),
      _omega(),
      _spSims(),
      _model(model)
{
}

  SimuSpectral::  SimuSpectral(const   SimuSpectral &r)
    : _ndim(r._ndim),
      _nb(r._nb),
      _nd(r._nd),
      _ns(r._ns),
      _isPrepared(r._isPrepared),
      _phi(r._phi),
      _gamma(r._gamma),
      _omega(r._omega),
      _spSims(r._spSims),
      _model(r._model)
{
}

  SimuSpectral&   SimuSpectral::operator=(const   SimuSpectral &r)
{
  if (this != &r)
  {
    _ndim = r._ndim;
    _nb = r._nb;
    _nd = r._nd;
    _ns = r._ns;
    _isPrepared = r._isPrepared;
    _phi = r._phi;
    _gamma = r._gamma;
    _omega = r._omega;
    _spSims = r._spSims;
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

int SimuSpectral::simulateOnSphere(int ns, int nd, int seed, bool verbose)
{
  if (_model == nullptr)
  {
    messerr("A Model should be attached beforehand");
    return 1;
  }
  if (! isValidForSpectral(_model)) return 1;
  if (ns <= 0)
  {
    messerr("The number of simulated harmonic components should be positive");
    return 1;
  }
  if (nd <= 0)
  {
    messerr("The number of degrees considered in the spectrum should be positive");
    return 1;
  }

  _ndim = _model->getDimensionNumber();
  _ns = ns;
  _nd = nd;

  law_set_random_seed(seed);

  _phi = VectorDouble(_ns);
  for (int is = 0; is < _ns; is++)
    _phi[is] = 2. * GV_PI * law_uniform();

  _simulateOnSphere(verbose);

  _isPrepared = true;
  return 0;
}

/**
 * Simulation of the spectral components (N,K) from spectrum values
 *
 * @param verbose Verbose flag
 *
 * @return It returns the list with two vectors of length _nb
 * @return N contains the simulated degrees, K contains the simulated order -N <= K <= N)
 */
void SimuSpectral::_simulateOnSphere(bool verbose)
{
  // Simulation of the spectrum
  VectorDouble U = VH::simulateUniform(_ns);
  VH::sortInPlace(U);
  double maxU = VH::maximum(U);

  VectorDouble spectrum = _model->getCova(0)->evalSpectrumOnSphere(_nd);

  // Simulate vector N
  int n = 0;
  double p = 0.;
  VectorInt N = VectorInt(_ns, 0);
  while (p < maxU && n < _ns)
  {
    p += spectrum[n++];
    for (int is = 0; is < _ns; is++)
    {
      if (U[is] > p) N[is]++;
    }
  }

  // Simulate vector K
  VectorInt K = VectorInt(_ns, ITEST);
  for (int is = 0; is < _ns; is++)
    K[is] = sampleInteger(-N[is], N[is]);

  // Derive the vector of orders
  VectorInt Kabs = K;
  for (int is = 0; is < _ns; is++) Kabs[is] = ABS(Kabs[is]);
  VectorInt orders = VH::unique(Kabs);
  int order_size = (int) orders.size();

  // Loop on the orders
  _spSims.resize(order_size);
  for (int kk = 0; kk < order_size; kk++)
  {
    int k = orders[kk];
    VectorInt Is;
    VectorInt Ns;
    int countP = 0;
    int countM = 0;
    for (int is = 0; is < _ns; is++)
    {
      int ki = K[is];
      if (ABS(ki) != k) continue;

      // Derive the restricted list of indices
      if (ki >= 0)
      {
        countP++;
        Is.push_back(1);
      }
      else
      {
        countM++;
        Is.push_back(-1);
      }
      Ns.push_back(N[is]);
    }

    // Create the table of contingency
    _spSims[kk]._k = kk;
    _spSims[kk]._countP = countP;
    _spSims[kk]._countM = countM;
    _spSims[kk]._tab = contingencyTable2(Ns, Is);
  }

  // Optional printout
  if (verbose) _printSpSims(1);
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

VectorInt SimuSpectral::_getKeys1(const spSim& spsim) const
{
  VectorInt keys;
  for (auto e: spsim._tab)
    keys.push_back(e.first);
  return keys;
}

int SimuSpectral::_getKey1Maximum(const spSim& spsim) const
{
  double maxi = -9999;
  for (auto e: spsim._tab)
    if (e.first > maxi) maxi = e.first;
  return maxi;
}

int SimuSpectral::_getSumValue(const spSim& spsim) const
{
  double sum = 0;
  for (const auto e1: spsim._tab)
  {
    for (const auto e2: e1.second)
      sum += e2.second;
  }
  return sum;
}

VectorInt SimuSpectral::_getKeys2(const spSim& spsim, int key1) const
{
  VectorInt keys;
  for (auto e1: spsim._tab)
  {
    if (e1.first != key1) continue;

    // The target key has been encountered
    for (const auto& e2: e1.second)
      keys.push_back(e2.first);
    return keys;
  }
  return keys;
}

VectorInt SimuSpectral::_getValues2(const spSim& spsim, int key1) const
{
  VectorInt keys;
  for (auto e1: spsim._tab)
  {
    if (e1.first != key1) continue;

    // The target key has been encountered
    for (const auto& e2: e1.second)
      keys.push_back(e2.second);
    return keys;
  }
  return keys;
}

void SimuSpectral::_printSpSim(const spSim& spsim, int status) const
{
  message("Component %2d (%2d / %2d)\n", spsim._k, spsim._countP, spsim._countM);
  if (status == 0) return;

  for (const auto& e1: spsim._tab)
  {
    message(" Key=%2d", e1.first);
    for (const auto& e2: e1.second)
    {
      message(" %2d", e2.second);
      if (e2.first > 0)
        message(" (+)");
      else
        message(" (-)");
    }
    message("\n");
  }
}

void SimuSpectral::_printSpSims(int status)
{
  int totalP = 0;
  int totalM = 0;
  int ns = (int) _spSims.size();
  mestitle(1,"List of Orders");
  for (int is = 0; is < ns; is++)
  {
    _printSpSim(_spSims[is], status);
    totalP += _spSims[is]._countP;
    totalM += _spSims[is]._countM;
  }

  message("\n");
  message("Summary:\n");
  message("- Number of Orders         = %d\n", ns);
  message("- Number of components (+) = %d\n", totalP);
  message("- Number of components (-) = %d\n", totalM);
}

int SimuSpectral::computeOnSphere(Db *dbout,
                                  bool verbose,
                                  const NamingConvention &namconv)
{
  int np   = dbout->getSampleNumber(true);
  int ndim = dbout->getNDim();

  // Preliminary checks

  if (ndim != _ndim)
  {
    messerr("The Space dimension of 'dbout'(%d) should match the one of Model(%d)",
            ndim, _ndim);
    return 1;
  }
  if (np <= 0)
  {
    messerr("'dbout' must have a positive number of active samples");
    return 1;
  }
  if (! _isPrepared)
  {
    messerr("You should run 'simulate' beforehand");
    return 1;
  }

  // Create the variable

  int iuid = dbout->addColumnsByConstant(1, 0., String(), ELoc::Z);
  if (iuid < 0) return 1;

  int nb = 0;
  int N_max = -9999;
  VectorInt K_list;
  for (int is = 0, size = (int) _spSims.size(); is < size; is++)
  {
    nb += _spSims[is]._countP + _spSims[is]._countM;
    K_list.push_back(_spSims[is]._k);
    int nmax = _getKey1Maximum(_spSims[is]);
    if (nmax > N_max) N_max = nmax;
  }
  int K_max = VH::maximum(K_list);

  // Optional printout
  if (verbose)
  {
    mestitle(1,">>> simulation on Sphere");
    message(">>> point number    : %d\n", np);
    message(">>> component number: %d\n", nb);
    message(">>> Maximum order   : %d\n", K_max);
    message(">>> Maximum degree  : %d\n", N_max);
  }

  // Simulation

  VectorDouble phi   = dbout->getCoordinates(0);
  VectorDouble theta = dbout->getCoordinates(1);
  VectorDouble sim(np,0.);
  VectorDouble x(np);
  VectorDouble w(np);
  for (int i = 0; i < np; i++)
  {
    double cosval = cos(theta[i]);
    x[i] = cosval;
    w[i] = sqrt(1 - cosval * cosval);
  }

  int K_idx = 0; // Index running in spectrum list
  int jk = 0;    // Index running in components
  int cumComp = 0;
  VectorDouble val(np, 0.);
  VectorDouble Pmm(np, 0.);
  VectorDouble Plm(np, 0.);
  VectorDouble P1(np, 0.);
  VectorDouble P2(np, 0.);

  for (int m = 0; m <= K_max; m++)
  {
    // From m-1 to m
    if (m == 0)
      Pmm.fill(0);
    else
    {
      double scale = sqrt((2.*m+1.)/(2.*m));
      for (int ip = 0; ip < np; ip++)
        Pmm[ip] = -scale * w[ip] * Pmm[ip];
    }

    if (VH::whereElement(K_list, m) >= 0)
    {
      const spSim& spsimK = _spSims[K_idx++];
      VectorInt N_list = _getKeys1(spsimK);

      if (verbose)
        message(">>> Simulating order K = %d: component number = %d\n",
                m, _getSumValue(spsimK));

      // From n-1 to n
      int NK_max = VH::maximum(N_list);
      for (int n = m; n <= NK_max; n++)
      {
        if (n == m)
        {
          Plm = Pmm;
          P2.fill(0);
          P1 = Plm;
        }
        else
        {
          double a = sqrt((2.*n+1.)*(2.*n-1.)/(n-m)/(n+m));
          double b = sqrt((2.*n+1.)/(2.*n-3.)*(n-1.-m)/(n-m)*(n-1.+m)/(n+m));
          for (int ip = 0; ip < np; ip++)
            Plm[ip] = a * x[ip] * P1[ip] - b * P2[ip];
          P2 = P1;
          P1 = Plm;
        }

        // Simulation
        if (VH::whereElement(N_list, n) >= 0)
        {
          VectorInt valComp = _getKeys2(spsimK, n);
          VectorInt nbrComp = _getValues2(spsimK, n);

          if (verbose)
          {
            int sumComp = VH::cumul(nbrComp);
            cumComp += sumComp;
            message("K = %d and N = %d : %d / %d  jk = %d\n",
                    m, n, sumComp, cumComp, jk);
          }

          for (int ii = 0, ncomp = (int) valComp.size(); ii < ncomp; ii++)
          {
            if (nbrComp[ii] > 0)
            {
              double s = valComp[ii];
              double fac = (s > 0) ? 1. : pow(-1., m);

              for (int ns = 0; ns < nbrComp[ii]; ns++)
              {

                for (int ip = 0; ip < np; ip++)
                  val[ip] += fac * Plm[ip] * cos(s * m * phi[ip] + _phi[jk]);
                jk++;
              }
            }
          }
        }
      }
    }
  }

  // Save the resulting array
  dbout->setColumnByUID(val, iuid);

  // Modify the name of the output
  namconv.setNamesAndLocators(dbout, iuid);
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
