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
#include "Covariances/ACov.hpp"
#include "Covariances/CovBase.hpp"
#include "Db/Db.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Model/Model.hpp"
#include "Model/ModelGeneric.hpp"
#include "Stats/Classical.hpp"

#include <cmath>

namespace gstlrn
{
SimuSpectral::SimuSpectral(const ACov* cova)
  : _ndim(0)
  , _ns(0)
  , _isPrepared(false)
  , _phi()
  , _gamma()
  , _omega()
  , _spSims()
  , _cova(cova)
{
}

SimuSpectral::SimuSpectral(const SimuSpectral& r)
  : _ndim(r._ndim)
  , _ns(r._ns)
  , _isPrepared(r._isPrepared)
  , _phi(r._phi)
  , _gamma(r._gamma)
  , _omega(r._omega)
  , _spSims(r._spSims)
  , _cova(r._cova)
{
}

SimuSpectral& SimuSpectral::operator=(const SimuSpectral& r)
{
  if (this != &r)
  {
    _ndim       = r._ndim;
    _ns         = r._ns;
    _isPrepared = r._isPrepared;
    _phi        = r._phi;
    _gamma      = r._gamma;
    _omega      = r._omega;
    _spSims     = r._spSims;
    _cova       = r._cova;
  }
  return *this;
}

SimuSpectral::~SimuSpectral()
{
}

/**
 * Perpare the simulation context
 * (This method works for Rn and Sphere)
 *
 * @param ns Number of components
 * @param seed Seed for random number generation 5à: avoid setting the seed)
 * @param verbose Verbose flag
 * @param nd Number of discretization of the spectrum
 */
Id SimuSpectral::simulate(Id ns, Id seed, bool verbose, Id nd)
{
  if (_cova == nullptr)
  {
    messerr("A Covariance should be attached beforehand");
    return 1;
  }
  if (!isValidForSpectral(_cova)) return 1;
  if (ns <= 0)
  {
    messerr("The number of simulated harmonic components should be positive");
    return 1;
  }
  if (nd <= 0)
  {
    messerr("The number of degrees considered in the spectrum should be positive (Sphere only)");
    return 1;
  }

  _ndim = static_cast<Id>(_cova->getNDim());
  _ns   = ns;

  // Cleaning any previously allocated memory
  _phi.clear();
  _gamma.clear();
  _omega.reset(0, 0);
  _spSims.clear();

  if (seed > 0) law_set_random_seed(seed);

  _phi       = VectorDouble(_ns);
  double pi2 = 2. * GV_PI;
  for (Id is = 0; is < _ns; is++)
    _phi[is] = pi2 * law_uniform();

  if (getDefaultSpaceType() == ESpaceType::RN)
    _simulateOnRn();
  else
    _simulateOnSphere(nd, verbose);

  _isPrepared = true;
  return 0;
}

void SimuSpectral::_simulateOnRn()
{
  _gamma = VectorDouble(_ns);
  for (Id ib = 0; ib < _ns; ib++)
    _gamma[ib] = sqrt(-log(law_uniform()));
  _omega = _cova->simulateSpectralOmega(_ns);
}

/**
 * Simulation of the spectral components (N,K) from spectrum values
 *
 * @param nd Number of discretization steps for Spectrum calculation
 * @param verbose Verbose flag
 */
void SimuSpectral::_simulateOnSphere(Id nd, bool verbose)
{
  // Simulation of the spectrum
  VectorDouble U = VH::simulateUniform(_ns);
  VH::sortInPlace(U);
  double maxU = VH::maximum(U);

  const auto* cl = dynamic_cast<const CovList*>(_cova);
  VectorDouble spectrum = cl->getCov(0)->evalSpectrumOnSphere(nd);

  // Simulate vector N
  Id n     = 0;
  double p = 0.;
  VectorInt N(_ns, 0);
  while (p < maxU && n < _ns)
  {
    p += spectrum[n++];
    for (Id is = 0; is < _ns; is++)
    {
      if (U[is] > p) N[is]++;
    }
  }

  // Simulate vector K
  VectorInt K(_ns, ITEST);
  for (Id is = 0; is < _ns; is++)
    K[is] = sampleInteger(-N[is], N[is]);

  // Derive the vector of orders
  VectorInt Kabs = K;
  for (Id is = 0; is < _ns; is++) Kabs[is] = ABS(Kabs[is]);
  VectorInt orders = VH::unique(Kabs);
  Id order_size    = static_cast<Id>(orders.size());

  // Loop on the orders
  _spSims.resize(order_size);
  for (Id kk = 0; kk < order_size; kk++)
  {
    Id k = orders[kk];
    VectorInt Is;
    VectorInt Ns;
    Id countP = 0;
    Id countM = 0;
    for (Id is = 0; is < _ns; is++)
    {
      Id ki = K[is];
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
    _spSims[kk]._k      = kk;
    _spSims[kk]._countP = countP;
    _spSims[kk]._countM = countM;
    _spSims[kk]._tab    = contingencyTable2(Ns, Is);
  }

  // Optional printout
  if (verbose) _printSpSims(1);
}

void SimuSpectral::_computeOnRn(Db* dbout, Id iuid, bool verbose)
{
  Id nech      = dbout->getNSample(true);
  double scale = sqrt(2. / _ns);

  // Optional printout
  if (verbose)
  {
    message("Spectral Simulation on a set of Isolated Points\n");
    message("- Number of samples = %d\n", nech);
    message("- Space dimension   = %d\n", _ndim);
    message("- Number of spectral components = %d\n", _ns);
  }

  // Loop on the active samples
  VectorInt ranks;
  dbout->getSampleRanksPerVariable(ranks);
  VectorDouble coor(_ndim);
  for (Id jech = 0; jech < nech; jech++)
  {
    Id iech = ranks[jech];
    dbout->getCoordinatesInPlace(coor, iech);
    VectorDouble u = _omega.prodMatVec(coor);

    double value = 0.;
    for (Id ib = 0; ib < _ns; ib++)
      value += _gamma[ib] * cos(u[ib] + _phi[ib]);
    value *= scale;

    dbout->setArray(iech, iuid, value);
  }
}

/**
 * Compute the simulation on Dbout using Spectral Method
 *
 * @param dbout Db containing the results
 * @param iuid  Address for storage (or 0 if the variable must be created locally)
 * @param verbose Verbose flag
 * @param namconv Naming convention (only used when 'iuid' == 0)
 */
Id SimuSpectral::compute(Db* dbout,
                         Id iuid,
                         bool verbose,
                         const NamingConvention& namconv)
{
  Id nech              = dbout->getNSample(true);
  Id ndim              = dbout->getNDim();
  bool flagNewVariable = (iuid <= 0);

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
  if (!_isPrepared)
  {
    messerr("You should run 'simulate' beforehand");
    return 1;
  }

  // Create the variable
  if (flagNewVariable)
  {
    iuid = dbout->addColumnsByConstant(1, 0., String(), ELoc::Z);
    if (iuid < 0) return 1;
  }

  if (getDefaultSpaceType() == ESpaceType::RN)
    _computeOnRn(dbout, iuid, verbose);
  else
    _computeOnSphere(dbout, iuid, verbose);

  // Modify the name of the output
  if (flagNewVariable) namconv.setNamesAndLocators(dbout, iuid);
  return 0;
}

VectorInt SimuSpectral::_getKeys1(const spSim& spsim)
{
  VectorInt keys;
  for (const auto& e: spsim._tab)
    keys.push_back(e.first);
  return keys;
}

Id SimuSpectral::_getKey1Maximum(const spSim& spsim)
{
  double maxi = -9999;
  for (const auto& e: spsim._tab)
    if (e.first > maxi) maxi = e.first;
  return maxi;
}

Id SimuSpectral::_getSumValue(const spSim& spsim)
{
  double sum = 0;
  for (const auto& e1: spsim._tab)
  {
    for (const auto& e2: e1.second)
      sum += e2.second;
  }
  return sum;
}

VectorInt SimuSpectral::_getKeys2(const spSim& spsim, Id key1)
{
  VectorInt keys;
  for (const auto& e1: spsim._tab)
  {
    if (e1.first != key1) continue;

    // The target key has been encountered
    for (const auto& e2: e1.second)
      keys.push_back(e2.first);
    return keys;
  }
  return keys;
}

VectorInt SimuSpectral::_getValues2(const spSim& spsim, Id key1)
{
  VectorInt keys;
  for (const auto& e1: spsim._tab)
  {
    if (e1.first != key1) continue;

    // The target key has been encountered
    for (const auto& e2: e1.second)
      keys.push_back(e2.second);
    return keys;
  }
  return keys;
}

void SimuSpectral::_printSpSim(const spSim& spsim, Id status)
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

void SimuSpectral::_printSpSims(Id status)
{
  Id totalP = 0;
  Id totalM = 0;
  Id ns     = static_cast<Id>(_spSims.size());
  mestitle(1, "List of Orders");
  for (Id is = 0; is < ns; is++)
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

void SimuSpectral::_computeOnSphere(Db* dbout, Id iuid, bool verbose)
{
  Id np = dbout->getNSample(true);

  Id nb    = 0;
  Id N_max = -9999;
  VectorInt K_list;
  for (Id is = 0, size = static_cast<Id>(_spSims.size()); is < size; is++)
  {
    nb += _spSims[is]._countP + _spSims[is]._countM;
    K_list.push_back(_spSims[is]._k);
    auto nmax = _getKey1Maximum(_spSims[is]);
    if (nmax > N_max) N_max = nmax;
  }
  Id K_max = VH::maximumVI(K_list);

  // Optional printout
  if (verbose)
  {
    mestitle(1, ">>> simulation on Sphere");
    message(">>> point number    : %d\n", np);
    message(">>> component number: %d\n", nb);
    message(">>> Maximum order   : %d\n", K_max);
    message(">>> Maximum degree  : %d\n", N_max);
  }

  // Simulation
  VectorDouble phi   = dbout->getOneCoordinate(0);
  VectorDouble theta = dbout->getOneCoordinate(1);
  VectorDouble sim(np, 0.);
  VectorDouble x(np);
  VectorDouble w(np);
  for (Id i = 0; i < np; i++)
  {
    double cosval = cos(theta[i]);
    x[i]          = cosval;
    w[i]          = sqrt(1 - cosval * cosval);
  }

  Id K_idx   = 0; // Index running in spectrum list
  Id jk      = 0; // Index running in components
  Id cumComp = 0;
  VectorDouble val(np, 0.);
  VectorDouble Pmm(np, 0.);
  VectorDouble Plm(np, 0.);
  VectorDouble P1(np, 0.);
  VectorDouble P2(np, 0.);

  for (Id m = 0; m <= K_max; m++)
  {
    // From m-1 to m
    if (m == 0)
      Pmm.fill(1.);
    else
    {
      double scale = sqrt((2. * m + 1.) / (2. * m));
      for (Id ip = 0; ip < np; ip++)
        Pmm[ip] = -scale * w[ip] * Pmm[ip];
    }

    if (VH::whereElement(K_list, m) >= 0)
    {
      const spSim& spsimK = _spSims[K_idx++];
      VectorInt N_list    = _getKeys1(spsimK);

      if (verbose)
        message(">>> Simulating order K = %d: component number = %d\n", m,
                _getSumValue(spsimK));

      // From n-1 to n
      Id NK_max = VH::maximumVI(N_list);
      for (Id n = m; n <= NK_max; n++)
      {
        if (n == m)
        {
          Plm = Pmm;
          P2.fill(0);
          P1 = Plm;
        }
        else
        {
          double a = sqrt((2. * n + 1.) * (2. * n - 1.) / (n - m) / (n + m));
          double b = sqrt((2. * n + 1.) / (2. * n - 3.) * (n - 1. - m) / (n - m) * (n - 1. + m) / (n + m));
          for (Id ip = 0; ip < np; ip++)
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
            Id sumComp = VH::cumul(nbrComp);
            cumComp += sumComp;
            message("K = %d and N = %d : %d / %d  jk = %d\n", m, n, sumComp,
                    cumComp, jk);
          }

          for (Id ii = 0, ncomp = static_cast<Id>(valComp.size()); ii < ncomp; ii++)
          {
            if (nbrComp[ii] > 0)
            {
              double s   = valComp[ii];
              double fac = (s > 0) ? 1. : pow(-1., m);

              for (Id ns = 0; ns < nbrComp[ii]; ns++)
              {
                for (Id ip = 0; ip < np; ip++)
                  val[ip] += fac * Plm[ip] * cos(s * m * phi[ip] + _phi[jk]);
                jk++;
              }
            }
          }
        }
      }
    }
  }

  // Normalize
  VH::multiplyConstant(val, sqrt(2. / nb));

  // Save the resulting array
  dbout->setColumnByUID(val, iuid);
}

/****************************************************************************/
/*!
 **  Check if the Model can be simulated using Spectral Method
 **
 ** \return  True if the covariance is valid; 0 otherwise
 **
 ** \param[in]  cov    Covariance structure
 **
 *****************************************************************************/
bool SimuSpectral::isValidForSpectral(const ACov* cov)
{
  return cov->isValidForSpectral();
}

/**
 * Perform a series of simulations (on Rn or on the Spehere) using Spectral Method
 *
 * @param dbin Input Db where the conditioning data are read
 * @param dbout Output Db where the results are stored
 * @param model Model (should only contain covariances that can cope with spectral method)
 * @param nbsimu Number of simulations processed simultaneously
 * @param seed Seed used for the Random number generator
 * @param ns Number of spectral components
 * @param nd Number of discretization steps (used for the Spectrum on Sphere)
 * @param verbose Verbose flag
 * @param namconv Naming Convention
 *
 * @note The conditional version is not yet available
 */
Id simuSpectral(Db* dbin,
                Db* dbout,
                ModelGeneric* model,
                Id nbsimu,
                Id seed,
                Id ns,
                Id nd,
                bool verbose,
                const NamingConvention& namconv)
{
  if (dbin != nullptr)
  {
    messerr("The current version does not allow Conditional Simulations");
    return 1;
  }
  if (nbsimu <= 0)
  {
    messerr("You must provide a positive number of simulations");
    return 1;
  }
  if (dbout->getNDim() != static_cast<Id>(model->getNDim()))
  {
    messerr("The Space dimension of 'dbout'(%d) should match the one of Model(%d)",
            dbout->getNDim(), model->getNDim());
    return 1;
  }

  // Instantiate the SimuSpectral class
  SimuSpectral simsph = SimuSpectral(model->getCov());

  // Set the seed for Random number generator (for all simulations)
  law_set_random_seed(seed);

  // Creating the output variables
  Id iuid = dbout->addColumnsByConstant(nbsimu, 0., String(), ELoc::Z);
  if (iuid < 0) return 1;

  // Loop on the simulations
  for (Id isimu = 0; isimu < nbsimu; isimu++)
  {
    if (simsph.simulate(ns, 0, verbose, nd)) return 1;
    if (simsph.compute(dbout, iuid + isimu, verbose)) return 1;
  }

  // Modify the name of the output
  namconv.setNamesAndLocators(dbin, VectorString(), ELoc::Z, 1, dbout, iuid, "", nbsimu);

  return 0;
}
} // namespace gstlrn