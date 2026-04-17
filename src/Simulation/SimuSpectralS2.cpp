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
#include "Simulation/SimuSpectralS2.hpp"
#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Simulation/CalcSimuSpectral.hpp"
#include "Stats/Classical.hpp"
#include "geoslib_define.h"

#include <cmath>

namespace gstlrn
{
  /*
   * -----------------------------------------
   * Simulation on the sphere S2
   * -----------------------------------------
   */

  SimuSpectralS2::SimuSpectralS2(Id nbsimu, Id ns, Id nd, Id seed, bool verbose)
    : CalcSimuSpectral(nbsimu, ns, nd, seed, verbose)
    , _spSims()
    , _phi()
  {
  }

  SimuSpectralS2::~SimuSpectralS2() {}

  /**
   * Simulate the spectrum components for Rn
   */
  Id SimuSpectralS2::_simulate()
  {
    const ACov* cov = getModelGeneric()->getCov();
    if (cov == nullptr) return -1;

    Id ns = _getNs();
    Id nd = _getNd();

    // Optional printout
    if (getVerbose())
    {
      message("Simulation of the spectrum\n");
      message("- Space dimension   = S%d\n", _getNDim());
      message("- Number of variables  = %d\n", _getNVar());
      message("- Number of spectral components = %d\n", ns);
    }

    // Cleaning any previously allocated memory
    _spSims.clear();

    // Simulation of the spectrum
    VectorDouble U = VH::simulateUniform(ns);
    VH::sortInPlace(U);
    double maxU = U.maximum();

    VectorDouble spectrum = cov->evalSpectrumOnSphere(nd);
    if (spectrum.empty()) return 1;

    // Simulate vector N
    Id n = 0;
    double p = 0.;
    VectorInt N(ns, 0);
    while (p < maxU && n < ns)
    {
      p += spectrum[n++];
      for (Id is = 0; is < ns; is++)
      {
        if (U[is] > p) N[is]++;
      }
    }

    // simulation of random phases, uniform on [0,2 pi]
    _phi.resize(ns);
    double pi2 = 2. * GV_PI;
    for (Id is = 0; is < ns; is++) _phi[is] = pi2 * law_uniform();

    // Simulate vector K
    VectorInt K(ns, ITEST);
    for (Id is = 0; is < ns; is++) K[is] = sampleInteger(-N[is], N[is]);

    // Derive the vector of orders
    VectorInt Kabs = K;
    for (Id is = 0; is < ns; is++) Kabs[is] = ABS(Kabs[is]);
    VectorInt orders = VH::unique(Kabs);
    Id order_size = static_cast<Id>(orders.size());

    // Loop on the orders
    _spSims.resize(order_size);
    for (Id kk = 0; kk < order_size; kk++)
    {
      Id k = orders[kk];
      VectorInt Is;
      VectorInt Ns;
      Id countP = 0;
      Id countM = 0;
      for (Id is = 0; is < ns; is++)
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
      _spSims[kk]._k = kk;
      _spSims[kk]._countP = countP;
      _spSims[kk]._countM = countM;
      _spSims[kk]._tab = contingencyTable2(Ns, Is);
    }

    // Optional printout
    if (getVerbose()) _printSpSims(1);

    return 0;
  }

  VectorInt SimuSpectralS2::_getKeys1(const spSim& spsim)
  {
    VectorInt keys;
    for (const auto& e: spsim._tab) keys.push_back(e.first);
    return keys;
  }

  Id SimuSpectralS2::_getKey1Maximum(const spSim& spsim)
  {
    double maxi = -9999;
    for (const auto& e: spsim._tab)
      if (e.first > maxi) maxi = e.first;
    return maxi;
  }

  Id SimuSpectralS2::_getSumValue(const spSim& spsim)
  {
    double sum = 0;
    for (const auto& e1: spsim._tab)
    {
      for (const auto& e2: e1.second) sum += e2.second;
    }
    return sum;
  }

  VectorInt SimuSpectralS2::_getKeys2(const spSim& spsim, Id key1)
  {
    VectorInt keys;
    for (const auto& e1: spsim._tab)
    {
      if (e1.first != key1) continue;

      // The target key has been encountered
      for (const auto& e2: e1.second) keys.push_back(e2.first);
      return keys;
    }
    return keys;
  }

  VectorInt SimuSpectralS2::_getValues2(const spSim& spsim, Id key1)
  {
    VectorInt keys;
    for (const auto& e1: spsim._tab)
    {
      if (e1.first != key1) continue;

      // The target key has been encountered
      for (const auto& e2: e1.second) keys.push_back(e2.second);
      return keys;
    }
    return keys;
  }

  void SimuSpectralS2::_printSpSim(const spSim& spsim, Id status)
  {
    message(
      "Component %2d (%2d / %2d)\n", spsim._k, spsim._countP, spsim._countM);
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

  void SimuSpectralS2::_printSpSims(Id status)
  {
    Id totalP = 0;
    Id totalM = 0;
    Id ns = static_cast<Id>(_spSims.size());
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

  Id SimuSpectralS2::_compute(
    Db* db,
    Id isimu,
    const VectorBool& activeArray,
    VectorVectorDouble& tab)
  {
    DECLARE_UNUSED(isimu);
    auto nech = db->getNSample();

    Id nb = 0;
    Id N_max = -9999;
    VectorInt K_list;
    for (Id is = 0, size = static_cast<Id>(_spSims.size()); is < size; is++)
    {
      nb += _spSims[is]._countP + _spSims[is]._countM;
      K_list.push_back(_spSims[is]._k);
      auto nmax = _getKey1Maximum(_spSims[is]);
      if (nmax > N_max) N_max = nmax;
    }
    Id K_max = K_list.maximum();

    // Optional printout
    if (getVerbose())
    {
      mestitle(1, ">>> Simulation on Sphere");
      message(">>> point number    : %d\n", nech);
      message(">>> component number: %d\n", nb);
      message(">>> Maximum order   : %d\n", K_max);
      message(">>> Maximum degree  : %d\n", N_max);
    }

    // Simulation
    VectorDouble phi = db->getOneCoordinate(0);
    VectorDouble theta = db->getOneCoordinate(1);
    VectorDouble sim(nech, 0.);
    VectorDouble x(nech);
    VectorDouble w(nech);
    for (Id iech = 0; iech < nech; iech++)
    {
      double cosval = cos(theta[iech]);
      x[iech] = cosval;
      w[iech] = sqrt(1 - cosval * cosval);
    }

    Id K_idx = 0; // Index running in spectrum list
    Id jk = 0; // Index running in components
    Id cumComp = 0;
    VectorDouble Pmm(nech, 0.);
    VectorDouble Plm(nech, 0.);
    VectorDouble P1(nech, 0.);
    VectorDouble P2(nech, 0.);

    for (Id m = 0; m <= K_max; m++)
    {
      // From m-1 to m
      if (m == 0)
        Pmm.fill(1.);
      else
      {
        double scale = sqrt((2. * m + 1.) / (2. * m));
        for (Id iech = 0; iech < nech; iech++)
          Pmm[iech] = -scale * w[iech] * Pmm[iech];
      }

      if (VH::whereElement(K_list, m) >= 0)
      {
        const spSim& spsimK = _spSims[K_idx++];
        VectorInt N_list = _getKeys1(spsimK);

        if (getVerbose())
          message(
            ">>> Simulating order K = %d: component number = %d\n", m,
            _getSumValue(spsimK));

        // From n-1 to n
        Id NK_max = N_list.maximum();
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
            double b = sqrt(
              (2. * n + 1.) / (2. * n - 3.) * (n - 1. - m) / (n - m)
              * (n - 1. + m) / (n + m));
            for (Id iech = 0; iech < nech; iech++)
              Plm[iech] = a * x[iech] * P1[iech] - b * P2[iech];
            P2 = P1;
            P1 = Plm;
          }

          // Simulation
          if (VH::whereElement(N_list, n) >= 0)
          {
            VectorInt valComp = _getKeys2(spsimK, n);
            VectorInt nbrComp = _getValues2(spsimK, n);

            if (getVerbose())
            {
              Id sumComp = nbrComp.sum();
              cumComp += sumComp;
              message(
                "K = %d and N = %d : %d / %d  jk = %d\n", m, n, sumComp,
                cumComp, jk);
            }

            for (Id ii = 0, ncomp = static_cast<Id>(valComp.size()); ii < ncomp;
                 ii++)
            {
              if (nbrComp[ii] > 0)
              {
                double s = valComp[ii];
                double fac = (s > 0) ? 1. : pow(-1., m);

                for (Id ns = 0; ns < nbrComp[ii]; ns++)
                {
                  for (Id iech = 0; iech < nech; iech++)
                    if (activeArray[iech])
                      tab[0][iech] +=
                        fac * Plm[iech] * cos(s * m * phi[iech] + getPhi(jk));
                  jk++;
                }
              }
            }
          }
        }
      }
    }

    // Normalize
    tab[0] *= sqrt(2. / nb);
    return 0;
  }

} // namespace gstlrn
