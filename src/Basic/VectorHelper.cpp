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
#include "Basic/VectorHelper.hpp"
#include "Basic/AException.hpp"
#include "Basic/Law.hpp"
#include "Basic/Message.hpp"
#include "Basic/Utilities.hpp"
#include "geoslib_define.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>

namespace gstlrn
{
  VectorInt VectorHelper::initVInt(Id nval, Id value)
  {
    VectorInt tab(nval, value);
    return tab;
  }

  VectorDouble VectorHelper::initVDouble(Id nval, double value)
  {
    VectorDouble tab(nval, value);
    return tab;
  }

  VectorVectorDouble
    VectorHelper::initVVDouble(Id nval1, Id nval2, double value)
  {
    VectorVectorDouble tab(nval1, VectorDouble(nval2, value));
    return tab;
  }

  VectorVectorInt VectorHelper::initVVInt(Id nval1, Id nval2, Id value)
  {
    VectorVectorInt tab(nval1, VectorInt(nval2, value));
    return tab;
  }

  VectorInt VectorHelper::initVInt(const Id* values, Id number)
  {
    if (values == nullptr) return VectorInt();
    VectorInt vec(number);
    for (Id i = 0; i < number; i++) vec[i] = values[i];
    return vec;
  }

  VectorInt VectorHelper::initVInt(const I32* values, Id number)
  {
    if (values == nullptr) return VectorInt();
    VectorInt vec(number);
    for (Id i = 0; i < number; i++) vec[i] = values[i];
    return vec;
  }

  VectorDouble VectorHelper::initVDouble(const double* values, Id number)
  {
    if (values == nullptr) return VectorDouble();
    VectorDouble vec(number);
    for (Id i = 0; i < number; i++) vec[i] = values[i];
    return vec;
  }

  VectorVectorDouble
    VectorHelper::initVVDouble(const double* value, Id n1, Id n2)
  {
    if (value == nullptr) return VectorVectorDouble();
    VectorVectorDouble vec;
    vec.resize(n1);
    for (Id i1 = 0; i1 < n1; i1++) vec[i1].resize(n2);

    Id lec = 0;
    for (Id i1 = 0; i1 < n1; i1++)
      for (Id i2 = 0; i2 < n2; i2++) vec[i1][i2] = value[lec++];
    return vec;
  }

  VectorString VectorHelper::initVString(Id ntab, char** names)
  {
    VectorString rettab(ntab);
    if (names == nullptr) return rettab;
    for (Id i = 0; i < ntab; i++) rettab[i] = names[i];
    return rettab;
  }

  void VectorHelper::dumpStats(const String& title, constvect vec, Id nmax)
  {
    Id ntotal = static_cast<Id>(vec.size());
    if (nmax > 0 && nmax < ntotal) ntotal = nmax;
    Id number = 0;
    double mean = 0.;
    double stdv = 0.;
    double mini = MAXIMUM_BIG;
    double maxi = MINIMUM_BIG;

    for (Id i = 0; i < ntotal; i++)
    {
      double value = vec[i];
      if (FFFF(value)) continue;
      number++;
      mean += value;
      stdv += value * value;
      if (value < mini) mini = value;
      if (value > maxi) maxi = value;
    }

    if (!title.empty()) message("%s\n", title.c_str());
    if (number > 0)
    {
      mean /= static_cast<double>(number);
      stdv = stdv / static_cast<double>(number) - mean * mean;
      stdv = (stdv > 0.) ? sqrt(stdv) : 0.;

      message("- Number of samples = %d / %d\n", number, ntotal);
      message("- Minimum  = %s\n", toStr(mini).c_str());
      message("- Maximum  = %s\n", toStr(maxi).c_str());
      message("- Mean     = %s\n", toStr(mean).c_str());
      message("- St. Dev. = %s\n", toStr(stdv).c_str());
    }
    else
    {
      message("No value defined\n");
    }
  }

  void VectorHelper::dumpStats(
    const String& title,
    const VectorDouble& vectin,
    Id nmax)
  {
    constvect vec(vectin);
    dumpStats(title, vec, nmax);
  }

  void VectorHelper::dumpRange(
    const String& title,
    const VectorDouble& vectin,
    Id nmax)
  {
    constvect vec(vectin);
    dumpRange(title, vec, nmax);
  }

  void VectorHelper::dumpRange(const String& title, constvect vec, Id nmax)
  {
    Id ntotal = static_cast<Id>(vec.size());
    if (nmax > 0 && nmax < ntotal) ntotal = nmax;
    Id number = 0;
    double mini = MAXIMUM_BIG;
    double maxi = MINIMUM_BIG;

    for (Id i = 0; i < ntotal; i++)
    {
      double value = vec[i];
      if (FFFF(value)) continue;
      number++;
      if (value < mini) mini = value;
      if (value > maxi) maxi = value;
    }

    if (!title.empty()) message("%s\n", title.c_str());
    if (number > 0)
    {
      message("- Number of samples = %d / %d\n", number, ntotal);
      message("- Minimum  = %lf\n", mini);
      message("- Maximum  = %lf\n", maxi);
    }
    else
    {
      message("No value defined\n");
    }
  }

  void VectorHelper::dumpRange(const String& title, const VectorInt& vec)
  {
    Id ntotal = static_cast<Id>(vec.size());
    Id number = 0;
    Id mini = 100000000;
    Id maxi = -100000000;

    for (Id i = 0; i < ntotal; i++)
    {
      Id value = vec[i];
      if (FFFF(value)) continue;
      number++;
      if (value < mini) mini = value;
      if (value > maxi) maxi = value;
    }

    if (!title.empty()) message("%s\n", title.c_str());
    if (number > 0)
    {
      message("- Number of samples = %d / %d\n", number, ntotal);
      message("- Minimum  = %d\n", mini);
      message("- Maximum  = %d\n", maxi);
    }
    else
    {
      message("No value defined\n");
    }
  }

  void VectorHelper::dumpNNZ(
    const String& title,
    const VectorDouble& vec,
    Id nclass)
  {
    Id ntotal = static_cast<Id>(vec.size());
    VectorInt total(nclass);
    for (Id ic = 0; ic < nclass; ic++) total[ic] = 0.;

    for (Id i = 0; i < ntotal; i++)
    {
      double value = ABS(vec[i]);
      double tol = 1.;
      for (Id ic = 0; ic < nclass; ic++)
      {
        tol /= 10.;
        if (value > tol) break;
        total[ic] += 1;
      }
    }

    if (!title.empty()) message("%s\n", title.c_str());
    for (Id ic = 0; ic < nclass; ic++)
      message("Count below 10.e-%d = %d\n", ic + 1, total[ic]);
  }

  Id VectorHelper::count(const VectorVectorInt& vec)
  {
    Id total = 0.;
    for (const auto& v: vec)
    {
      total += static_cast<Id>(v.size());
    }
    return total;
  }

  double VectorHelper::cumulLog(const VectorDouble& vec)
  {
    double total = 0.;
    for (const auto& v: vec)
    {
      if (!FFFF(v)) total += log(v);
    }
    return total;
  }

  VectorInt VectorHelper::cumulIncrement(const VectorVectorInt& vec)
  {
    Id nvar = static_cast<Id>(vec.size());
    VectorInt cumul(nvar, 0);
    Id number = 0;
    for (Id ivar = 0; ivar < nvar; ivar++)
    {
      cumul[ivar] = number;
      number += static_cast<Id>(vec[ivar].size());
    }
    return cumul;
  }

  /**
   * Calculate the quantiles
   * @param vec    Array of values
   * @param probas Array of probabilities (sorted by ascending order)
   * @return Vector of data values for the different probabilities
   */
  VectorDouble
    VectorHelper::quantiles(const VectorDouble& vec, const VectorDouble& probas)
  {
    Id nproba = static_cast<Id>(probas.size());
    Id nech = static_cast<Id>(vec.size());
    if (nech <= 0 || nproba <= 0) return VectorDouble();

    VectorDouble retval(nproba, TEST);

    // Sort the data in ascending order
    VectorDouble sorted = VH::sort(vec, true);

    for (Id ip = 0; ip < nproba; ip++)
    {
      double proba = probas[ip];
      Id rank = static_cast<Id>(proba * static_cast<double>(nech));

      double value = TEST;
      if (rank < 0)
      {
        value = TEST;
      }
      else if (rank < nech - 1)
      {
        double v1 = sorted[rank];
        double v2 = sorted[rank + 1];
        double p1 = static_cast<double>(rank) / static_cast<double>(nech);
        double p2 = static_cast<double>(1 + rank) / static_cast<double>(nech);
        value = v1 + (proba - p1) * (v2 - v1) / (p2 - p1);
      }
      else
      {
        value = sorted[nech - 1];
      }
      retval[ip] = value;
    }
    return retval;
  }

  void VectorHelper::normalize(double* tab, Id ntab)
  {
    Id i;
    double norme;

    norme = 0.;
    for (i = 0; i < ntab; i++) norme += tab[i] * tab[i];
    norme = sqrt(norme);

    if (norme <= 0.) return;
    for (i = 0; i < ntab; i++) tab[i] /= norme;
  }

  void VectorHelper::normalizeFromGaussianDistribution(
    VectorDouble& vec,
    double mini,
    double maxi)
  {
    double* iptr = vec.data();
    for (Id i = 0, n = static_cast<Id>(vec.size()); i < n; i++)
    {
      if (!FFFF(*iptr))
        (*iptr) = mini + (maxi - mini) * law_cdf_gaussian(*iptr);
      iptr++;
    }
  }

  VectorDouble VectorHelper::qnormVec(const VectorDouble& vec)
  {
    Id number = static_cast<Id>(vec.size());
    VectorDouble retvec(number, TEST);
    for (Id i = 0; i < number; i++) retvec[i] = law_invcdf_gaussian(vec[i]);
    return retvec;
  }

  VectorDouble VectorHelper::pnormVec(const VectorDouble& vec)
  {
    Id number = static_cast<Id>(vec.size());
    VectorDouble retvec(number, TEST);
    for (Id i = 0; i < number; i++) retvec[i] = law_cdf_gaussian(vec[i]);
    return retvec;
  }

  VectorDouble
    VectorHelper::normalScore(const VectorDouble& data, const VectorDouble& wt)
  {
    Id nech = static_cast<Id>(data.size());
    VectorDouble vec(nech, TEST);
    if (nech <= 0) return vec;

    // Check dimension of vector
    if (!wt.empty() && nech != static_cast<Id>(wt.size()))
    {
      messerr("Arguments 'data' and 'wt' should have the same dimension");
      return VectorDouble();
    }

    // Check that weights of active samples are positive
    double wtotal = 0.;
    double nechtot = 0.;
    for (Id iech = 0; iech < nech; iech++)
    {
      if (!FFFF(data[iech]))
      {
        double wtloc = 1.;
        if (!wt.empty()) wtloc = wt[iech];
        if (wtloc < 0.)
        {
          messerr(
            "The weight of sample (%d) is negative (%lf)", iech + 1, wtloc);
          return VectorDouble();
        }
        wtotal += wtloc;
        nechtot += 1;
      }
    }
    if (wtotal <= 0.)
    {
      messerr("The sum of weights of active samples is not positive");
      return VectorDouble();
    }
    wtotal *= (1. + nechtot) / nechtot;

    // Get the list of indices sorted by increasing values of data
    VectorInt idx = VH::orderRanks(data);

    // Loop on the samples
    double wpartial = 0.;
    for (Id iech = 0; iech < nech; iech++)
    {
      Id jech = idx[iech];
      if (!FFFF(data[jech]))
      {
        double wtloc = (wt.empty()) ? 1. : wt[jech];
        wpartial += wtloc;
        vec[jech] = law_invcdf_gaussian(wpartial / wtotal);
      }
      else
      {
        vec[jech] = TEST;
      }
    }
    return vec;
  }

  bool VectorHelper::isEqual(
    const VectorDouble& v1,
    const VectorDouble& v2,
    double eps)
  {
    if (v1.size() != v2.size()) return false;
    auto it1(v1.begin());
    auto it2(v2.begin());
    while (it1 < v1.end())
    {
      if (ABS(*it1 - *it2) > eps) return false;
      it1++;
      it2++;
    }
    return true;
  }

  void VectorHelper::fillUndef(VectorDouble& vec, double repl)
  {
    auto it(vec.begin());
    while (it < vec.end())
    {
      if (FFFF(*it)) *it = repl;
      it++;
    }
  }

  void VectorHelper::sequenceInPlace(Id n, VectorInt& vec)
  {
    vec.resize(n);
    for (Id i = 0; i < n; i++) vec[i] = i;
  }

  /**
   * Create an output vector containing the 'number' consecutive numbers starting from 'ideb'
   *
   * @param number  Length of the output vector
   * @param ideb    Index of the first element of the output vector
   * @param step    Step between two consecutive values
   */
  VectorInt VectorHelper::sequence(Id number, Id ideb, Id step)
  {
    VectorInt vec(number);

    Id jdeb = ideb;
    for (Id i = 0; i < number; i++)
    {
      vec[i] = jdeb;
      jdeb += step;
    }
    return vec;
  }

  /**
   * Create an output vector going from 'valFrom' to 'ValTo' by step of 'valStep'
   */

  /**
   * Create a vector containing the a sequence of numbers
   * @param valFrom Starting value
   * @param valTo   Ending value
   * @param valStep Step
   * @param ratio   The whole sequence can be ultimately scaled by 'ratio'
   * @return
   */
  VectorDouble VectorHelper::sequenceVD(
    double valFrom,
    double valTo,
    double valStep,
    double ratio)
  {
    VectorDouble vec;

    double value = valFrom;
    while (value <= valTo)
    {
      vec.push_back(value / ratio);
      value = value + valStep;
    }
    return vec;
  }

  VectorBool VectorHelper::simulateBoolean(Id n, double probaTrue)
  {
    VectorBool vec(n);
    for (auto& el: vec) el = law_uniform(0., 1.) < probaTrue;
    return vec;
  }

  VectorInt VectorHelper::simulateInteger(Id n, const VectorDouble& probas)
  {
    VectorInt vec(n);

    // Normalize the probabilities
    double total = 0.;
    for (const auto& p: probas)
    {
      if (p < 0.)
      {
        messerr("Probabilities should be positive. Nothing is done.");
        return VectorInt();
      }
      total += p;
    }
    if (total <= 0.)
    {
      messerr("The sum of probabilities should be positive. Nothing is done.");
      return VectorInt();
    }
    auto nproba = static_cast<Id>(probas.size());
    VectorDouble normprobas(nproba);
    for (Id i = 0; i < nproba; i++) normprobas[i] = probas[i] / total;

    for (auto& el: vec)
    {
      double rand = law_uniform(0., 1.);
      Id ic = 0;
      double tol = normprobas[0];
      while (rand > tol && ic < nproba - 1)
      {
        ic++;
        tol += normprobas[ic];
      }
      el = ic;
    }
    return vec;
  }

  VectorDouble VectorHelper::simulateUniform(
    Id n,
    double mini,
    double maxi,
    double normSum)
  {
    VectorDouble vec(n);
    double total = 0.;
    for (auto& el: vec)
    {
      el = law_uniform(mini, maxi);
      total += el;
    }
    if (normSum != TEST && total != 0.)
    {
      for (auto& el: vec) el /= total;
    }
    return vec;
  }

  VectorDouble VectorHelper::simulateBernoulli(
    Id n,
    double proba,
    double vone,
    double velse)
  {
    VectorDouble vec(n);
    for (auto& el: vec)
    {
      double rand = law_uniform(0., 1.);
      if (rand < proba)
        el = vone;
      else
        el = velse;
    }
    return vec;
  }

  VectorDouble VectorHelper::simulateGaussian(Id n, double mean, double sigma)
  {
    VectorDouble vec(n);
    simulateGaussianInPlace(vec, mean, sigma);
    return vec;
  }

  void VectorHelper::simulateGaussianInPlace(
    VectorDouble& vec,
    double mean,
    double sigma)
  {
    for (auto& el: vec)
    {
      el = mean + sigma * law_gaussian();
    }
  }

  VectorDouble VectorHelper::concatenate(
    const VectorDouble& veca,
    const VectorDouble& vecb)
  {
    VectorDouble res = veca;
    for (const auto& e: vecb) res.push_back(e);
    return res;
  }

  void VectorHelper::concatenateInPlace(
    VectorDouble& veca,
    const VectorDouble& vecb)
  {
    for (const auto& e: vecb) veca.push_back(e);
  }

  void VectorHelper::cumulateInPlace(VectorDouble& vec)
  {
    auto it(vec.begin());
    double old = 0.;
    while (it < vec.end())
    {
      *it += old;
      old = *it;
      it++;
    }
  }

  VectorDouble VectorHelper::cumsum(
    const VectorDouble& vecin,
    bool flagAddZero,
    bool revert)
  {
    VectorDouble vecout;
    if (flagAddZero) vecout.push_back(0.);

    double total = 0.;
    for (Id i = 0, n = static_cast<Id>(vecin.size()); i < n; i++)
    {
      total += vecin[i];
      vecout.push_back(total);
    }

    if (revert)
    {
      Id size = static_cast<Id>(vecout.size());
      double lastval = vecout[size - 1];
      for (Id i = 0; i < size; i++) vecout[i] = lastval - vecout[i];
    }
    return vecout;
  }

  void VectorHelper::cumulate(
    VectorDouble& veca,
    const VectorDouble& vecb,
    double coeff,
    double addval)
  {
    if (veca.size() != vecb.size())
    {
      messerr(
        "Arguments 'veca' and 'vecb' should have the same dimension. Nothing "
        "is done");
      return;
    }

    auto ita(veca.begin());
    auto itb(vecb.begin());
    while (ita < veca.end())
    {
      *ita += coeff * (*itb) + addval;
      ita++;
      itb++;
    }
  }

  /**
   * Display the first significant values of the input vector.
   * A "significant" value is a value larger than 'tol' in absolute value
   * Values are listed by decreasing importance.
   * @param vec  Input Vector
   * @param tol  Tolerance above which a value is significant (in absolute value)
   * @param nmax Limit on the number of values printed (-1: no limit)
   */
  void VectorHelper::getMostSignificant(
    const VectorDouble& vec,
    double tol,
    Id nmax)
  {
    Id nsize = static_cast<Id>(vec.size());
    VectorDouble absval(nsize, 0.);
    Id ninvalid = 0;
    for (Id i = 0; i < nsize; i++)
    {
      double value = vec[i];
      if (FFFF(value)) continue;
      value = ABS(value);
      if (value <= tol) continue;
      absval[i] = value;
      ninvalid++;
    }

    if (ninvalid <= 0) return;

    VectorInt ranks = orderRanks(absval, false);
    Id nend = ninvalid;
    if (nmax > 0) nend = MIN(ninvalid, nmax);
    for (Id i = 0; i < nend; i++)
    {
      Id j = ranks[i];
      message("Sample %d - Value = %lf\n", j, vec[j]);
    }
    if (nmax > 0 && ninvalid > nmax)
      message(
        "Found %d (out of %d) samples. Print limited to the %d most important "
        "ones.\n",
        ninvalid, nsize, nmax);
  }

  /**
   * Sample a set of 'ntotal' ranks (unique occurrence)
   * @param ntotal      Dimension to be sampled
   * @param proportion  Proportion of elected samples (in [0,1])
   * @param number      Number of elected samples
   * @param seed        Seed used for the random number generator
   * @param optSort     Sorting: 0 for None; 1 for Ascending; -1 for Descending
   * @return A vector of indices lying between 0 and ntotal-1. No duplicate.
   *
   * @remark If 'proportion' and 'number' are not specified,
   * @remark the output vector has dimension equal to 'ntotal'
   */
  VectorInt VectorHelper::sampleRanks(
    Id ntotal,
    double proportion,
    Id number,
    Id seed,
    Id optSort)
  {
    if (proportion <= 0. && number <= 0) return VectorInt();
    law_set_random_seed(seed);

    // Find the number of expected values
    Id count;
    if (proportion <= 0. && number <= 0)
      count = ntotal;
    else if (proportion > 0.)
      count = static_cast<Id>(ntotal * proportion);
    else
      count = number;
    count = MIN(ntotal, MAX(1, count));

    VectorInt ranks = law_random_path(ntotal);

    ranks.resize(count);

    // Sort them out
    if (optSort > 0)
      ranks = sort(ranks, true);
    else if (optSort < 0)
      ranks = sort(ranks, false);

    VectorInt::iterator it;
    it = std::unique(ranks.begin(), ranks.end());
    ranks.resize(distance(ranks.begin(), it));

    return ranks;
  }

  void VectorHelper::addInPlace(constvect in, vect dest)
  {
    const double* inp = in.data();
    double* outp = dest.data();
    for (Id i = 0, n = static_cast<Id>(in.size()); i < n; i++)
    {
      *(outp++) += *(inp++);
    }
  }

  /**
   * Performs: veca += vecb**2
   * @param dest Input/Output vector
   * @param src Auxiliary vector
   */
  void
    VectorHelper::addSquareInPlace(VectorDouble& dest, const VectorDouble& src)
  {
    if (dest.size() != src.size())
    {
      messerr(
        "Arguments 'dest' and 'src' should have the same dimension. Nothing is "
        "done");
      return;
    }

    auto itd(dest.begin());
    auto its(src.begin());
    while (itd < dest.end())
    {
      *itd += (*its) * (*its);
      itd++;
      its++;
    }
  }

  void VectorHelper::addInPlace(
    const double* veca,
    const double* vecb,
    double* res,
    Id size)
  {
    for (Id i = 0; i < size; i++) res[i] = veca[i] + vecb[i];
  }

  /**
   * Return a vector containing vecb - veca
   * @param veca Input Vector
   * @param vecb Input Vector
   * @return
   */
  VectorDouble VectorHelper::subtract(constvect veca, constvect vecb)
  {
    VectorDouble res(veca.size());
    for (Id i = 0; i < static_cast<Id>(veca.size()); i++)
    {
      res[i] = vecb[i] - veca[i];
    }
    return res;
  }

  /**
   * Performs: outv = in2 - in1
   * @param in1 Input vector
   * @param in2 Input vector
   * @param outv Output vector
   */
  void VectorHelper::subtractInPlace(
    const constvect in1,
    const constvect in2,
    vect outv)
  {
    for (Id is = 0, ns = static_cast<Id>(in1.size()); is < ns; is++)
    {
      outv[is] = in2[is] - in1[is];
    }
  }

  void VectorHelper::multiplyComplexInPlace(
    const VectorDouble& vecaRe,
    const VectorDouble& vecaIm,
    const VectorDouble& vecbRe,
    const VectorDouble& vecbIm,
    VectorDouble& resRe,
    VectorDouble& resIm)
  {
    VectorDouble temp(vecaRe);

    multiply(resRe, vecaRe, vecbRe);
    multiply(temp, vecaIm, vecbIm);
    resRe -= temp;

    multiply(resIm, vecaRe, vecbIm);
    multiply(temp, vecaIm, vecbRe);
    resIm += temp;
  }

  void VectorHelper::addMultiplyConstantInPlace(
    double val1,
    const VectorDouble& in,
    VectorDouble& out,
    Id iad)
  {
    double* outp = out.data() + iad;
    const double* inp = in.data();
    for (Id i = 0; i < static_cast<Id>(in.size()); i++)
    {
      *(outp++) += val1 * *(inp++);
    }
  }

  void VectorHelper::addMultiplyVectVectInPlace(
    const constvect in1,
    const constvect in2,
    vect out,
    Id iad)
  { // TODO check if one can use eigen operators
    double* outp = out.data() + iad;
    const double* inp1 = in1.data();
    const double* inp2 = in2.data();
    for (Id i = 0; i < static_cast<Id>(in1.size()); i++)
    {
      *(outp++) += *(inp1++) * *(inp2++);
    }
  }

  void VectorHelper::addMultiplyConstantInPlace(
    double val1,
    const constvect in,
    vect out,
    Id iad)
  {
    double* outp = out.data() + iad;
    const double* inp = in.data();
    for (Id i = 0; i < static_cast<Id>(in.size()); i++)
    {
      *(outp++) += val1 * *(inp++);
    }
  }

  void VectorHelper::addMultiplyConstantInPlace(
    double val1,
    const VectorVectorDouble& in1,
    VectorVectorDouble& outv)
  {
    for (Id is = 0, ns = static_cast<Id>(in1.size()); is < ns; is++)
    {
      for (Id i = 0, n = static_cast<Id>(in1[is].size()); i < n; i++)
      {
        outv[is][i] += val1 * in1[is][i];
      }
    }
  }

  void
    VectorHelper::copy(const VectorDouble& vecin, VectorDouble& vecout, Id size)
  {
    if (size < 0) size = static_cast<Id>(vecin.size());
    if (size > static_cast<Id>(vecout.size())) my_throw("Wrong size");

    auto itout(vecout.begin());
    auto itin(vecin.begin());
    for (Id i = 0; i < size; i++)
    {
      (*itout) = (*itin);
      itin++;
      itout++;
    }
  }

  void VectorHelper::copy(const VectorInt& vecin, VectorInt& vecout, Id size)
  {
    if (size < 0) size = static_cast<Id>(vecin.size());
    if (size > static_cast<Id>(vecout.size())) my_throw("Wrong size");

    auto itout(vecout.begin());
    auto itin(vecin.begin());
    for (Id i = 0; i < size; i++)
    {
      (*itout) = (*itin);
      itin++;
      itout++;
    }
  }

  void
    VectorHelper::copy(const VectorVectorDouble& inv, VectorVectorDouble& outv)
  {
    for (Id is = 0, ns = static_cast<Id>(inv.size()); is < ns; is++)
    {
      for (Id i = 0, n = static_cast<Id>(inv[is].size()); i < n; i++)
      {
        outv[is][i] = inv[is][i];
      }
    }
  }

  void VectorHelper::mean1AndMean2ToStdev(
    const VectorDouble& mean1,
    const VectorDouble& mean2,
    VectorDouble& std,
    Id number)
  {
    auto dnumber = static_cast<double>(number);
    Id size = static_cast<Id>(mean1.size());
    if (static_cast<Id>(mean2.size()) != size)
    {
      messerr(
        "Arguments 'mean1'(%d) and 'mean2'(%d) should have same dimension",
        size, static_cast<Id>(mean2.size()));
      return;
    }
    if (static_cast<Id>(std.size()) != size)
    {
      messerr(
        "Arguments 'mean1'(%d) and 'std'(%d) should have same dimension", size,
        static_cast<Id>(std.size()));
      return;
    }

    for (Id i = 0; i < size; i++)
    {
      if (FFFF(mean1[i]) || FFFF(mean2[i]))
        std[i] = TEST;
      else
      {
        double dmean1 = mean1[i] / dnumber;
        double dmean2 = mean2[i] / dnumber;
        double value = dmean2 - dmean1 * dmean1;
        std[i] = (value > 0) ? sqrt(value) : 0.;
      }
    }
  }

  VectorDouble VectorHelper::power(const VectorDouble& vec, double power)
  {
    VectorDouble res;
    VectorHelper::power(res, vec, power);
    return res;
  }

  VectorDouble VectorHelper::inverse(const VectorDouble& vec)
  {
    VectorDouble res;
    VectorHelper::inverse(res, vec);
    return res;
  }

  void VectorHelper::power(VectorDouble& res, const constvect vec, double power)
  {
    res.resize(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
    {
      res[i] = pow(vec[i], power);
    }
  }

  void VectorHelper::inverse(VectorDouble& res, const constvect vec)
  {
    res.resize(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
    {
      res[i] = 1.0 / vec[i];
    }
  }

  /**
   * Calculate the diagonal of the box extension
   * @param mini Array of lower coordinates of the box
   * @param maxi Array of upper coordinates of the box
   * @return
   * @remark If one coordinate is undefined, TEST is returned.
   */
  double VectorHelper::extensionDiagonal(
    const VectorDouble& mini,
    const VectorDouble& maxi)
  {
    double diag = 0.;
    VectorDouble delta = maxi - mini;
    Id ndim = static_cast<Id>(delta.size());
    for (Id idim = 0; idim < ndim; idim++)
    {
      double dval = delta[idim];
      if (FFFF(dval)) return TEST;
      diag += dval * dval;
    }
    diag = sqrt(diag);
    return diag;
  }

  VectorInt VectorHelper::unique(const VectorInt& vecin, Id size)
  {
    if (size < 0) size = static_cast<Id>(vecin.size());

    VectorInt vecout = vecin;
    vecout.resize(size);
    std::sort(vecout.begin(), vecout.end());
    auto last = std::unique(vecout.begin(), vecout.end());
    vecout.erase(last, vecout.end());
    return vecout;
  }

  VectorDouble VectorHelper::unique(const VectorDouble& vecin, Id size)
  {
    if (size < 0) size = static_cast<Id>(vecin.size());

    VectorDouble vecout = vecin;
    vecout.resize(size);
    std::sort(vecout.begin(), vecout.end());
    auto last = std::unique(vecout.begin(), vecout.end());
    vecout.erase(last, vecout.end());
    return vecout;
  }

  bool VectorHelper::isInList(const VectorInt& vec, Id item)
  {
    return std::count(vec.begin(), vec.end(), item);
  }

  VectorInt VectorHelper::sort(const VectorInt& vecin, bool ascending, Id size)
  {
    if (vecin.empty()) return VectorInt();
    if (size < 0) size = static_cast<Id>(vecin.size());
    VectorInt vecout = vecin;
    vecout.resize(size);
    std::sort(vecout.begin(), vecout.end());
    if (!ascending) std::reverse(vecout.begin(), vecout.end());
    return vecout;
  }

  VectorDouble
    VectorHelper::sort(const VectorDouble& vecin, bool ascending, Id size)
  {
    if (vecin.empty()) return VectorDouble();
    if (size < 0) size = static_cast<Id>(vecin.size());
    VectorDouble vecout = vecin;
    vecout.resize(size);
    std::sort(vecout.begin(), vecout.end());
    if (!ascending) std::reverse(vecout.begin(), vecout.end());
    return vecout;
  }

  void VectorHelper::sortInPlace(VectorInt& vecin, bool ascending, Id size)
  {
    if (vecin.empty()) return;
    VectorInt vecout = sort(vecin, ascending, size);
    copy(vecout, vecin, size);
  }

  void VectorHelper::sortInPlace(VectorDouble& vecin, bool ascending, Id size)
  {
    if (vecin.empty()) return;
    VectorDouble vecout = sort(vecin, ascending, size);
    copy(vecout, vecin, size);
  }

  bool VectorHelper::isSorted(const VectorDouble& vec, bool ascending)
  {
    Id nval = static_cast<Id>(vec.size());

    if (ascending)
    {
      // Ascending order
      for (Id i = 1; i < nval; i++)
      {
        if (vec[i] > vec[i - 1]) continue;
        return false;
      }
    }
    else
    {
      // Descending order
      for (Id i = 1; i < nval; i++)
      {
        if (vec[i] < vec[i - 1]) continue;
        return false;
      }
    }
    return true;
  }

  /**
   * From an input list, filter out all the elements which do no lie within [vmin, vmax[,
   * suppress double occurrences and sort them out (ascending or descending)
   * @param vecin Input array (integer)
   * @param vmin  lower bound included (or ITEST)
   * @param vmax  upper bound excluded (or ITEST)
   * @param ascending True for ascending order; False for descending order
   * @return Output array (integers)
   */
  VectorInt VectorHelper::filter(
    const VectorInt& vecin,
    Id vmin,
    Id vmax,
    bool ascending)
  {
    VectorInt vecout = vecin;

    // Sort the vector
    std::sort(vecout.begin(), vecout.end());
    if (!ascending) std::reverse(vecout.begin(), vecout.end());

    // Unique occurrence
    VectorInt::iterator it;
    it = std::unique(vecout.begin(), vecout.end());
    vecout.resize(distance(vecout.begin(), it));

    // Filter out the irrelevant values
    Id nech = static_cast<Id>(vecout.size());
    for (Id j = 0; j < nech; j++)
    {
      Id i = nech - j - 1;
      if (!isNA(vmin))
      {
        if (vecout[i] < vmin)
        {
          vecout.erase(vecout.begin() + i);
          continue;
        }
      }
      if (!isNA(vmax))
      {
        if (vecout[i] >= vmax)
        {
          vecout.erase(vecout.begin() + i);
          continue;
        }
      }
    }
    return vecout;
  }

  /**
   * Returns the list complementary to 'sel' within 'vecin'
   * @param vec Initial list
   * @param sel Vector of forbidden elements
   * @return Complementary list
   */
  VectorInt VectorHelper::complement(const VectorInt& vec, const VectorInt& sel)
  {
    VectorInt rest;
    if (vec.empty()) return rest;
    if (sel.empty()) return vec;

    // Sort

    VectorInt allVec = vec;
    std::sort(allVec.begin(), allVec.end());

    VectorInt offVec = sel;
    std::sort(offVec.begin(), offVec.end());

    Id j, k, idx;
    Id nvec = static_cast<Id>(allVec.size());
    Id noff = static_cast<Id>(offVec.size());
    for (Id i = 0; i < nvec; i++)
    {
      j = allVec.at(i);

      // I go through offVec as long as element is strictly less than j
      k = 0;
      idx = offVec.at(k);
      while (idx < j && k < noff)
      {
        idx = offVec.at(k++);
      }

      if (idx != j) // idx not in offElemsVec
      {
        rest.push_back(j);
      }
    }
    return rest;
  }

  /**
   * Returns the permutation which rearranges the input vector into any order
   * @param vecin Input vector
   * @param ascending True for ascending order; False otherwise
   * @param size Optional dimension of the input vector
   * @return Vector of orders
   */
  VectorInt
    VectorHelper::orderRanks(const VectorInt& vecin, bool ascending, Id size)
  {
    if (vecin.empty()) return VectorInt();

    if (size < 0) size = static_cast<Id>(vecin.size());
    VectorInt idx(size);
    for (Id i = 0; i < size; i++) idx[i] = i;

    // sort indexes based on comparing values in v using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings when v contains elements of equal values
    auto last = idx.begin() + size;
    if (ascending)
    {
      stable_sort(
        idx.begin(), last,
        [&vecin](size_t i1, size_t i2) { return vecin[i1] < vecin[i2]; });
    }
    else
    {
      stable_sort(
        idx.begin(), last,
        [&vecin](size_t i1, size_t i2) { return vecin[i1] > vecin[i2]; });
    }

    return idx;
  }

  /**
   * @overload
   */
  VectorInt
    VectorHelper::orderRanks(const VectorDouble& vecin, bool ascending, Id size)
  {
    if (vecin.empty()) return VectorInt();
    if (size < 0) size = static_cast<Id>(vecin.size());
    VectorInt idx(size);
    for (Id i = 0; i < size; i++) idx[i] = i;

    // sort indexes based on comparing values in v using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings when v contains elements of equal values
    auto last = idx.begin() + size;
    if (ascending)
    {
      stable_sort(
        idx.begin(), last,
        [&vecin](size_t i1, size_t i2) { return vecin[i1] < vecin[i2]; });
    }
    else
    {
      stable_sort(
        idx.begin(), last,
        [&vecin](size_t i1, size_t i2) { return vecin[i1] > vecin[i2]; });
    }

    return idx;
  }

  VectorInt
    VectorHelper::sortRanks(const VectorDouble& vecin, bool ascending, Id size)
  {
    if (vecin.empty()) return VectorInt();
    if (size < 0) size = static_cast<Id>(vecin.size());
    VectorInt order = orderRanks(vecin, ascending, size);
    VectorInt idx(size);
    for (Id i = 0; i < size; i++) idx[order[i]] = i;

    return idx;
  }

  VectorInt
    VectorHelper::sortRanks(const VectorInt& vecin, bool ascending, Id size)
  {
    if (vecin.empty()) return VectorInt();
    if (size < 0) size = static_cast<Id>(vecin.size());
    VectorInt order = orderRanks(vecin, ascending, size);
    VectorInt idx(size);
    for (Id i = 0; i < size; i++) idx[order[i]] = i;

    return idx;
  }

  VectorInt VectorHelper::reorder(
    const VectorInt& vecin,
    const VectorInt& order,
    Id size)
  {
    if (size < 0) size = static_cast<Id>(vecin.size());
    VectorInt vecout(size);
    for (Id i = 0; i < size; i++) vecout[i] = vecin[order[i]];
    return vecout;
  }

  VectorDouble VectorHelper::reorder(
    const VectorDouble& vecin,
    const VectorInt& order,
    Id size)
  {
    if (size < 0) size = static_cast<Id>(vecin.size());
    VectorDouble vecout(size);
    for (Id i = 0; i < size; i++) vecout[i] = vecin[order[i]];
    return vecout;
  }

  VectorDouble VectorHelper::revert(const VectorDouble& vecin)
  {
    Id nech = static_cast<Id>(vecin.size());
    VectorDouble vecout(nech);
    for (Id iech = 0; iech < nech; iech++)
      vecout[nech - 1 - iech] = vecin[iech];
    return vecout;
  }

  VectorInt VectorHelper::revert(const VectorInt& vecin)
  {
    Id nech = static_cast<Id>(vecin.size());
    VectorInt vecout(nech);
    for (Id iech = 0; iech < nech; iech++)
      vecout[nech - 1 - iech] = vecin[iech];
    return vecout;
  }

  /*****************************************************************************/
  /*!
   ** Sorts the (double) array value() and the array ranks() if provided
   *  @brief Arrange values in place (VectorDouble version)
   ** \param[in]  safe   1 if the value array if preserved
   **                    0 if the value array is also sorted
   ** \param[in]  ascending Sorting order
   ** \param[in]  size   Optional vector dimension
   **
   ** \param[out] ranks  input and output Id array
   ** \param[out] values input and output double array
   **
   ** \remark  If ranks = NULL, ranks is ignored
   ** \remark  When using 'size', the remaining part of arrays is unchanged
   **
   *****************************************************************************/
  void VectorHelper::arrangeInPlace(
    Id safe,
    VectorInt& ranks,
    VectorDouble& values,
    bool ascending,
    Id size)
  {
    VectorInt order = orderRanks(values, ascending, size);

    if (!ranks.empty())
    {
      VectorInt newranks = reorder(ranks, order, size);
      copy(newranks, ranks, size);
    }

    if (safe == 0)
    {
      VectorDouble newvals = reorder(values, order, size);
      copy(newvals, values, size);
    }
  }

  /*****************************************************************************/
  /*!
   ** Sorts the static_cast<Id>(array) value() and the array ranks() if provided
   ** @overload
   **/
  void VectorHelper::arrangeInPlace(
    Id safe,
    VectorInt& ranks,
    VectorInt& values,
    bool ascending,
    Id size)
  {
    VectorInt order = orderRanks(values, ascending, size);

    if (!ranks.empty())
    {
      VectorInt newranks = reorder(ranks, order, size);
      ranks = newranks;
    }

    if (safe == 0)
    {
      VectorInt newval = reorder(values, order, size);
      values = newval;
    }
  }

  std::pair<double, double> VectorHelper::rangeVals(const VectorDouble& vec)
  {
    std::pair<double, double> res(vec[0], vec[0]);
    for (Id i = 1; i < static_cast<Id>(vec.size()); i++)
    {
      res.first = MIN(res.first, vec[i]);
      res.second = MAX(res.second, vec[i]);
    }
    return res;
  }

  double
    VectorHelper::innerProductCV(const constvect veca, const constvect vecb)
  {
    return innerProduct(veca.data(), vecb.data(), static_cast<Id>(veca.size()));
  }

  double
    VectorHelper::innerProduct(const double* veca, const double* vecb, Id size)
  {
    double prod = 0.;
    const double* ptra = &veca[0];
    const double* ptrb = &vecb[0];
    for (Id i = 0; i < size; i++)
    {
      prod += (*ptra) * (*ptrb);
      ptra++;
      ptrb++;
    }
    return prod;
  }

  /**
   * Cross product (limited to 3D)
   * @param veca First vector
   * @param vecb Second Vector
   * @return
   */
  VectorDouble VectorHelper::crossProduct3D(
    const VectorDouble& veca,
    const VectorDouble& vecb)
  {
    if (veca.size() != vecb.size()) my_throw("Wrong size");
    VectorDouble res(3);
    crossProduct3DInPlace(veca.data(), vecb.data(), res.data());
    return res;
  }

  void VectorHelper::crossProduct3DInPlace(
    const double* a,
    const double* b,
    double* v)
  {
    v[0] = a[1] * b[2] - a[2] * b[1];
    v[1] = a[2] * b[0] - a[0] * b[2];
    v[2] = a[0] * b[1] - a[1] * b[0];
  }

  /**
   * Method which flattens a VectorVectorDouble into a VectorDouble
   * @param vvd Input VectorVectorDouble
   * @return Returned VectorDouble
   */
  VectorDouble VectorHelper::flatten(const VectorVectorDouble& vvd)
  {
    VectorDouble vd;

    for (Id i = 0; i < static_cast<Id>(vvd.size()); i++)
      for (Id j = 0; j < static_cast<Id>(vvd[i].size()); j++)
        vd.push_back(vvd[i][j]);

    return vd;
  }

  void VectorHelper::flattenInPlace(
    const VectorVectorDouble& vvd,
    VectorDouble& vd)
  {
    Id ecr = 0;
    for (Id i = 0; i < static_cast<Id>(vvd.size()); i++)
      for (Id j = 0; j < static_cast<Id>(vvd[i].size()); j++)
        vd[ecr++] = (vvd[i][j]);
  }

  void VectorHelper::unflattenInPlace(
    const VectorDouble& vd,
    VectorVectorDouble& vvd)
  {
    Id lec = 0;
    for (Id i = 0, n = static_cast<Id>(vvd.size()); i < n; i++)
      for (Id j = 0; j < static_cast<Id>(vvd[i].size()); j++)
        vvd[i][j] = vd[lec++];
  }

  VectorVectorDouble
    VectorHelper::unflatten(const VectorDouble& vd, const VectorInt& sizes)
  {
    VectorVectorDouble vvd;

    Id lec = 0;
    for (Id i = 0, n = static_cast<Id>(sizes.size()); i < n; i++)
    {
      Id lng = sizes[i];
      VectorDouble local(lng);
      for (Id j = 0; j < lng; j++) local[j] = vd[lec++];
      vvd.push_back(local);
    }
    return vvd;
  }

  VectorDouble VectorHelper::suppressTest(const VectorDouble& vecin)
  {
    VectorDouble vecout;
    for (Id i = 0, n = static_cast<Id>(vecin.size()); i < n; i++)
    {
      if (!FFFF(vecin[i])) vecout.push_back(vecin[i]);
    }
    return vecout;
  }

  void VectorHelper::linearCombinationInPlace(
    double val1,
    const VectorDouble& vd1,
    double val2,
    const VectorDouble& vd2,
    VectorDouble& outv)
  {
    if (vd1.empty() || vd2.empty()) return;
    for (Id i = 0, n = static_cast<Id>(vd1.size()); i < n; i++)
    {
      double value = 0.;
      if (val1 != 0. && !vd1.empty()) value += val1 * vd1[i];
      if (val2 != 0. && !vd2.empty()) value += val2 * vd2[i];
      outv[i] = value;
    }
  }

  void VectorHelper::linearCombinationVVDInPlace(
    double val1,
    const VectorVectorDouble& vvd1,
    double val2,
    const VectorVectorDouble& vvd2,
    VectorVectorDouble& outv)
  {
    if (vvd1.empty() || vvd2.empty()) return;

    for (Id is = 0, ns = static_cast<Id>(vvd1.size()); is < ns; is++)
    {
      for (Id i = 0, n = static_cast<Id>(vvd1[is].size()); i < n; i++)
      {
        double value = 0.;
        if (val1 != 0. && !vvd1.empty()) value += val1 * vvd1[is][i];
        if (val2 != 0. && !vvd2.empty()) value += val2 * vvd2[is][i];
        outv[is][i] = value;
      }
    }
  }

  /**
   * Extract the part of a vector 'vecin' (dimensioned to dimension of 'vecout')
   * starting at address 'istart' and copy it into 'vecout'
   * @param vecin  Initial vector
   * @param vecout Resulting vector (already allocated)
   * @param start  Starting address (within 'vecin')
   */
  void VectorHelper::extractInPlace(
    const VectorDouble& vecin,
    VectorDouble& vecout,
    Id start)
  {
    std::copy(
      vecin.begin() + start, vecin.begin() + start + vecout.size(),
      vecout.begin());
  }

  /**
   * Merge 'vecin' into 'vecout' starting at address 'istart'
   * @param vecin  Initial vector
   * @param vecout Vector where 'vecin' should be copied
   * @param start  Starting address (in 'vecout')
   */
  void VectorHelper::mergeInPlace(
    const VectorDouble& vecin,
    VectorDouble& vecout,
    Id start)
  {
    std::copy(vecin.begin(), vecin.end(), vecout.begin() + start);
  }

  /**
   * Transform a vector of double values as follows
   * @param tab  Vector of double values
   * @param oper_choice Operation on the diagonal term (see Utilities::operate_XXX)
   */
  void VectorHelper::transformVD(VectorDouble& tab, Id oper_choice)
  {
    operate_function oper_func = operate_Identify(oper_choice);
    Id number = static_cast<Id>(tab.size());
    for (Id i = 0; i < number; i++) tab[i] = oper_func(tab[i]);
  }

  /****************************************************************************/
  /*!
   **  Fix plausible values for the Direction coefficients.
   **  They must be defined and with norm equal to 1
   **
   ** \param[in]  ndim      Space dimension
   ** \param[in,out]  codir Input/Output Direction coefficients
   **
   *****************************************************************************/
  void VectorHelper::normalizeCodir(Id ndim, VectorDouble& codir)
  {
    double norme;

    if (codir.empty()) return;
    norme = codir.norm2();
    if (norme <= 0.)
    {
      for (Id idim = 0; idim < ndim; idim++) codir[idim] = 0.;
      codir[0] = 1.;
    }
    else
    {
      norme = sqrt(norme);
      for (Id idim = 0; idim < ndim; idim++) codir[idim] /= norme;
    }
  }

  /**
   * Operate the squeeze-and-stretch algorithm forward (see remarks)
   * @param vecin  Input vector (in structural system)
   * @param vecout Output vector (in sugar box system)
   * @param origin Origin of the vertical axis (structural system)
   * @param mesh   Mesh of the vertical axis (structural system)
   * @param top    Elevation of the Top surface
   * @param bot    Elevation of the Bottom surface
   *
   * @remarks The information is contained in 'vecin' which is defined on a regular 1D grid
   * @remarks in the structural system. The purpose is to sample the relevant sub-information
   * @remarks (between 'top' and 'bot') densely in 'vecout'
   */
  void VectorHelper::squeezeAndStretchInPlaceForward(
    const VectorDouble& vecin,
    VectorDouble& vecout,
    double origin,
    double mesh,
    double top,
    double bot)
  {
    Id nzin = static_cast<Id>(vecin.size());
    Id nzout = static_cast<Id>(vecout.size());
    double thick = top - bot;
    double ratio = thick / nzout;

    // Loop on the positions of the pile in the sugar box system
    for (Id iz = 0; iz < nzout; iz++)
    {
      // Corresponding coordinate of the sample in the structural system
      double zzin = bot + static_cast<double>(iz) * ratio;

      // Find the index in the input vector
      Id izin = static_cast<Id>((zzin - origin) / mesh);
      if (izin < 0 || izin >= nzin) continue;

      // Assign the value
      vecout[iz] = vecin[izin];
    }
  }

  /**
   * Operate the squeeze-and-stretch algorithm backward (see remarks)
   * @param vecin  Input vector (in sugar box system)
   * @param vecout Output vector (in structural system)
   * @param origin Origin of the vertical axis (structural system)
   * @param mesh   Mesh of the vertical axis (structural system)
   * @param top    Elevation of the Top surface
   * @param bot    Elevation of the Bottom surface
   *
   * @remarks The information is contained in 'vecin' which is defined on a regular 1D grid
   * @remarks (characterized by 'base' and 'mesh')
   * @remarks Extend the relevant information, lying between 'bot' and 'top' in order to fill
   * @remarks the whole vector 'vecout'
   */
  void VectorHelper::squeezeAndStretchInPlaceBackward(
    const VectorDouble& vecin,
    VectorDouble& vecout,
    double origin,
    double mesh,
    double top,
    double bot)
  {
    Id nzin = static_cast<Id>(vecin.size());
    Id nzout = static_cast<Id>(vecout.size());

    // Blank out the output vector
    vecout.fill(TEST);
    double thick = top - bot;
    if (thick <= 0) return;

    // Get the top and bottom indices in the output vector
    Id indbot = floor((bot - origin) / mesh);
    if (indbot < 0) indbot = 0;
    Id indtop = ceil((top - origin) / mesh);
    if (indtop >= nzout) indtop = nzout - 1;

    double ratio = static_cast<double>(nzin) / thick;

    // Loop on the positions of the pile in the structural system
    for (Id izout = indbot; izout <= indtop; izout++)
    {
      // Get the location of the sample in the structural system
      double zzout = origin + izout * mesh;

      // Find the index in the input vector (sugar box)
      Id izin = ratio * (zzout - bot);
      if (izin < 0 || izin >= nzin) continue;

      // Assign the value
      vecout[izout] = vecin[izin];
    }
  }

  /*****************************************************************************/
  /*!
   **  Find the location of the minimum value within a vector
   **
   ** \return Rank of the minimum value
   **
   ** \param[in]  tab  Vector of values
   **
   *****************************************************************************/
  Id VectorHelper::whereMinimum(const VectorDouble& tab)
  {
    Id ibest = -1;
    double vbest = MAXIMUM_BIG;
    for (Id i = 0, ntab = static_cast<Id>(tab.size()); i < ntab; i++)
    {
      if (FFFF(tab[i])) continue;
      if (tab[i] > vbest) continue;
      vbest = tab[i];
      ibest = i;
    }
    return ibest;
  }

  /*****************************************************************************/
  /*!
   **  Find the location of the maximum value within a vector
   **
   ** \return Rank of the maximum value
   **
   ** \param[in]  tab  Vector of values
   **
   *****************************************************************************/
  Id VectorHelper::whereMaximum(const VectorDouble& tab)
  {
    Id ibest = -1;
    double vbest = MINIMUM_BIG;
    for (Id i = 0, ntab = static_cast<Id>(tab.size()); i < ntab; i++)
    {
      if (FFFF(tab[i])) continue;
      if (tab[i] < vbest) continue;
      vbest = tab[i];
      ibest = i;
    }
    return ibest;
  }

  /*
   * Returns the rank where 'target' has been found within 'tab'
   *
   * @param tab Vector of integer values to be searched
   * @param target Target value to be searched for
   *
   * @return Rank at which the target value has been found (-1 if not found)
   */
  Id VectorHelper::whereElement(const VectorInt& tab, Id target)
  {
    for (Id i = 0, ntab = static_cast<Id>(tab.size()); i < ntab; i++)
    {
      if (tab[i] == target) return i;
    }
    return -1;
  }

  /*
   * Returns the rank where 'target' has been found within 'tab' starting from 'start'
   * This optimized version is useful when the vector is sorted and we're searching
   * for increasing values sequentially.
   *
   * @param tab Vector of integer values to be searched (assumed to be sorted)
   * @param target Target value to be searched for
   * @param start Starting position for the search (default: 0)
   *
   * @return Rank at which the target value has been found (-1 if not found)
   */
  Id VectorHelper::whereElement(const VectorInt& tab, Id target, Id start)
  {
    Id ntab = static_cast<Id>(tab.size());
    // Search from start to end
    for (Id i = start; i < ntab; i++)
    {
      if (tab[i] == target) return i;
    }
    return -1;
  }

  /**
   * Reduce the input vector 'vecin' by suppressing the element referred by 'index'
   *
   * @param vecin Input vector (double)
   * @param index Index to be suppressed
   */
  VectorDouble VectorHelper::reduceOne(const VectorDouble& vecin, Id index)
  {
    VectorInt vindex(1);
    vindex[0] = index;
    return reduce(vecin, vindex);
  }

  /**
   * Reduce the input vector 'vecin' by suppressing the elements referred by 'index'
   *
   * @param vecin Input vector (double)
   * @param vindex Vector of indices to be suppressed
   */
  VectorDouble
    VectorHelper::reduce(const VectorDouble& vecin, const VectorInt& vindex)
  {
    VectorDouble vecout = vecin;

    // Sort the indices to be removed in ascending order
    VectorInt indexLocal = vindex;
    std::sort(indexLocal.begin(), indexLocal.end());

    Id nsel = static_cast<Id>(indexLocal.size());
    for (Id j = 0; j < nsel; j++)
    {
      Id i = indexLocal[nsel - j - 1];
      vecout.erase(vecout.begin() + i);
    }
    return vecout;
  }

  /**
   * Reduce the input vector 'vecin' by returning the only elements referred by 'index'
   *
   * @param vecin Input vector (double)
   * @param vindex Vector of indices to be kept
   */
  VectorDouble
    VectorHelper::compress(const VectorDouble& vecin, const VectorInt& vindex)
  {
    VectorDouble vecout;
    for (Id j = 0, nsel = static_cast<Id>(vindex.size()); j < nsel; j++)
    {
      Id i = vindex[j];
      vecout.push_back(vecin[i]);
    }
    return vecout;
  }

  void VectorHelper::truncateDecimalsInPlace(VectorDouble& vec, Id ndec)
  {
    for (Id i = 0, n = static_cast<Id>(vec.size()); i < n; i++)
    {
      if (FFFF(vec[i])) continue;
      vec[i] = truncateDecimals(vec[i], ndec);
    }
  }

  void VectorHelper::truncateDigitsInPlace(VectorDouble& vec, Id ndec)
  {
    for (Id i = 0, n = static_cast<Id>(vec.size()); i < n; i++)
    {
      if (FFFF(vec[i])) continue;
      vec[i] = truncateDigits(vec[i], ndec);
    }
  }

  /**
   * @brief Create an output VectorDouble by selecting some indices
   *        of the Input VectorDouble 'vecin'
   *
   * @param vecin    Input Rectangular Matrix
   * @param indKeep  Set of Indices to be kept (all if not defined)
   */
  VectorDouble
    VectorHelper::sample(const VectorDouble& vecin, const VectorInt& indKeep)
  {
    VectorDouble vecout;

    VectorInt indices = indKeep;
    if (indices.empty()) indices = VH::sequence(static_cast<Id>(vecin.size()));

    Id nindices = static_cast<Id>(indices.size());
    if (nindices <= 0) return vecout;

    for (Id i = 0; i < nindices; i++)
    {
      if (!checkArg(
            "Selected index", indices[i], static_cast<Id>(vecin.size())))
        return vecout;
    }

    vecout.resize(nindices);
    for (Id i = 0; i < nindices; i++) vecout[i] = vecin[indices[i]];
    return vecout;
  }

  /**
   * @brief Function checking that two values are equal
   * This verbose option is essentially used in tests
   *
   * @param v1 First value to be compared
   * @param v2  Second value to be compared
   * @param eps Tolerance used for comparison
   * @param flagRelative when True, the values are compared without paying
   * attention to their sign
   * @param flagAbsolute when True, test is run on absolute difference
   * @param string Message to be displayed when the vectors are not similar
   * @return Boolean
   *
   * @note: When the two vectors do not share the same dimension, the test is not
   * performed and a message is printed.
   */
  bool VectorHelper::isEqualExtended(
    const VectorDouble& v1,
    const VectorDouble& v2,
    double eps,
    bool flagRelative,
    bool flagAbsolute,
    const String& string)
  {
    // Check that the two vectors have the same dimension
    if (v1.size() != v2.size())
    {
      if (!string.empty()) message("%s : ", string.c_str());
      message("Impossible to compare vectors of different dimensions\n");
      return false;
    }
    Id size = static_cast<Id>(v1.size());
    VectorDouble vec1 = v1;
    VectorDouble vec2 = v2;

    // Check is performed on the absolute value of each term of each vector
    if (flagAbsolute)
    {
      for (Id i = 0; i < size; i++)
      {
        vec1[i] = ABS(vec1[i]);
        vec2[i] = ABS(vec2[i]);
      }
    }

    // Evaluate the comparison test
    double diff = 0.;
    for (Id i = 0; i < size; i++)
    {
      double value = (vec1[i] - vec2[i]);
      if (flagRelative) value /= (vec1[i] + vec2[i] + eps);
      diff += value * value;
    }
    diff = sqrt(diff) / size;

    if (diff >= eps)
    {
      if (!string.empty()) message("%s : ", string.c_str());
      message(
        "Experimental value = %lf is larger than tolerance (%lf)\n", diff, eps);
      return false;
    }
    return true;
  }

  bool VectorHelper::isIsotropic(const VectorVectorInt& sampleRanks)
  {
    Id nvar = static_cast<Id>(sampleRanks.size());
    if (nvar <= 0) return true;

    Id refSize = static_cast<Id>(sampleRanks[0].size());
    for (Id ivar = 1; ivar < nvar; ivar++)
      if (refSize != static_cast<Id>(sampleRanks[ivar].size())) return false;
    return true;
  }

  void VectorHelper::capInPlace(VectorDouble& vec, double vmin, double vmax)
  {
    for (Id i = 0, n = static_cast<Id>(vec.size()); i < n; i++)
    {
      double value = vec[i];
      if (!FFFF(vmin))
      {
        if (FFFF(value) || value < vmin) value = vmin;
      }
      if (!FFFF(vmax))
      {
        if (FFFF(value) || value > vmax) value = vmax;
      }
      vec[i] = value;
    }
  }

  void VectorHelper::capInPlaceVVD(
    VectorVectorDouble& vec,
    double vmin,
    double vmax)
  {
    for (Id i = 0, n = static_cast<Id>(vec.size()); i < n; i++)
      VH::capInPlace(vec[i], vmin, vmax);
  }

} // namespace gstlrn
