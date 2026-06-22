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
#include "Variogram/VarioOrder.hpp"
#include "Basic/VectorHelper.hpp"

#include <cmath>
#include <cstring>

#define QUANT_DIR 10000

namespace gstlrn
{
  VarioOrder::VarioOrder(Id flag_dist, Id size_aux)
    : _nalloc(0)
    , _npair(0)
    , _sizeAux(size_aux)
    , _flagDist(flag_dist)
    , _tabIech()
    , _tabJech()
    , _tabIpas()
    , _tabSort()
    , _tabAuxIech()
    , _tabAuxJech()
    , _tabDist()
  {
  }

  Id VarioOrder::final()
  {
    _nalloc = _npair;
    if (_npair > 0)
    {
      _tabIech.resize(_npair);
      _tabJech.resize(_npair);
      _tabIpas.resize(_npair);
      _tabSort.resize(_npair);
      if (_flagDist) _tabDist.resize(_npair);

      if (_sizeAux > 0)
      {
        _tabAuxIech.resize(_npair * _sizeAux);
        _tabAuxJech.resize(_npair * _sizeAux);
      }
      for (Id i = 0; i < _npair; i++) _tabSort[i] = i;
      VectorHelper::arrangeInPlace(1, _tabSort, _tabIpas, true, _npair);
    }
    return _npair;
  }

  void VarioOrder::clear()
  {
    _nalloc = 0;
    _npair = 0;
    _tabIech.clear();
    _tabJech.clear();
    _tabIpas.clear();
    _tabSort.clear();
    _tabAuxIech.clear();
    _tabAuxJech.clear();
    _tabDist.clear();
  }

  /****************************************************************************/
  /*!
   **  Print the VarioOrder structure
   **
   ** \param[in]  idir_target Rank of the target direction (starting from 0) or -1
   ** \param[in]  ipas_target Rank of the target lag (starting from 0) or -1
   ** \param[in]  verbose     1 for a complete printout
   **
   *****************************************************************************/
  void VarioOrder::printout(Id idir_target, Id ipas_target, Id verbose)
  {
    mestitle(0, "Variogram Order structure");
    message("Allocated size    = %d\n", _nalloc);
    message("Number of pairs   = %d\n", _npair);
    if (!verbose) return;

    bool flagFirst = true;

    for (Id i = 0; i < _npair; i++)
    {
      Id j = (_tabSort.empty()) ? i : _tabSort[i];
      Id ilag = _tabIpas[j];
      Id idir = ilag / QUANT_DIR;
      ilag = ilag - QUANT_DIR * idir;
      if (idir_target >= 0 && idir != idir_target) continue;
      if (ipas_target >= 0 && ilag != ipas_target) continue;

      if (flagFirst)
      {
        if (!_flagDist)
          message("Rank - Dir - Lag - I - J\n");
        else
          message("Rank - Dir - Lag - I - J - Dist\n");
        flagFirst = false;
      }

      message("%5d", i + 1);
      message(" %5d", idir + 1);
      message(" %5d", ilag + 1);
      message(" %5d", _tabIech[j] + 1);
      message(" %5d", _tabJech[j] + 1);
      if (_flagDist) message(" %lf", _tabDist[j]);
      message("\n");
    }
  }

  /****************************************************************************/
  /*!
   **  Returns the first and last indices matching a target lag
   **
   ** \param[in]  idir        Rank of the target direction
   ** \param[in]  ilag        Rank of the target lag
   **
   ** \param[out] ifirst      Rank of the first sample of the lag (included)
   ** \param[out] ilast       Rank of the last sample of the lag (excluded)
   **
   *****************************************************************************/
  void VarioOrder::getBounds(Id idir, Id ilag, Id* ifirst, Id* ilast) const
  {
    Id ipair, jpair;

    Id ival = ilag + idir * QUANT_DIR;
    if (_npair > 0 && _tabSort.empty()) messageAbort("getBounds");
    *ifirst = _npair;
    *ilast = -1;
    for (ipair = 0; ipair < _npair; ipair++)
    {
      jpair = _tabSort[ipair];
      if (_tabIpas[jpair] == ival)
      {
        if (ipair < *ifirst) *ifirst = ipair;
      }
      else
      {
        if (ipair > *ifirst)
        {
          *ilast = ipair;
          return;
        }
      }
    }

    /* Particular case of the last lag */

    if (*ifirst < _npair)
    {
      *ilast = _npair;
      return;
    }
  }

  /****************************************************************************/
  /*!
   **  Returns the two samples for a given (ordered) pair
   **
   ** \param[in]  ipair       Rank of the sorted pair
   **
   ** \param[out] iech        Rank of the first sample
   ** \param[out] jech        Rank of the second sample
   ** \param[out] dist        Calculated distance or TEST (if flag_dist == 0)
   **
   *****************************************************************************/
  void VarioOrder::getIndices(Id ipair, Id* iech, Id* jech, double* dist) const
  {
    if (_tabSort.empty()) messageAbort("getIndices");
    Id jpair = _tabSort[ipair];
    *iech = _tabIech[jpair];
    *jech = _tabJech[jpair];
    *dist = (_flagDist) ? _tabDist[jpair] : TEST;
  }

  /****************************************************************************/
  /*!
   **  Returns the two auxiliary arrays for a given (ordered) pair
   **
   ** \param[in]  ipair       Rank of the sorted pair
   **
   ** \param[out] aux_iech    Array to auxiliary information for sample 'iech'
   ** \param[out] aux_jech    Array to auxiliary information for sample 'jech'
   **
   *****************************************************************************/
  void VarioOrder::getAuxiliary(Id ipair, char* aux_iech, char* aux_jech) const
  {
    if (_tabSort.empty()) messageAbort("getAuxiliary");
    Id jpair = _tabSort[ipair];
    Id iad = _sizeAux * jpair;
    (void)memcpy(aux_iech, &_tabAuxIech[iad], _sizeAux);
    (void)memcpy(aux_jech, &_tabAuxJech[iad], _sizeAux);
  }

  /****************************************************************************/
  /*!
   **  Add a record to the Variogram Order structure
   **
   ** \return Error return code
   **
   ** \param[in]  iech        Rank of the first sample
   ** \param[in]  jech        Rank of the second sample
   ** \param[in]  aux_iech    Auxiliary array for sample 'iech' (or NULL)
   ** \param[in]  aux_jech    Auxiliary array for sample 'jech' (or NULL)
   ** \param[in]  ilag        Rank of the lag
   ** \param[in]  idir        Rank of the direction (or 0)
   ** \param[in]  dist        Calculated distance (only stored if flag_dist == 1)
   **
   *****************************************************************************/
  Id VarioOrder::add(
    Id iech,
    Id jech,
    void* aux_iech,
    void* aux_jech,
    Id ilag,
    Id idir,
    double dist)
  {
    static Id VARIOORDER_QUANT = 1000;

    /* Resize the array */

    if (_npair >= _nalloc)
    {
      _nalloc += VARIOORDER_QUANT;
      _tabIech.resize(_nalloc);
      _tabJech.resize(_nalloc);
      _tabIpas.resize(_nalloc);
      _tabSort.resize(_nalloc);
      if (_sizeAux > 0)
      {
        _tabAuxIech.resize(_nalloc * _sizeAux);
        _tabAuxJech.resize(_nalloc * _sizeAux);
      }
      if (_flagDist) _tabDist.resize(_nalloc);
    }

    /* Add the new information */

    _tabIech[_npair] = (dist > 0) ? iech : jech;
    _tabJech[_npair] = (dist > 0) ? jech : iech;
    _tabIpas[_npair] = ilag + idir * QUANT_DIR;
    if (_flagDist) _tabDist[_npair] = dist;
    if (_sizeAux > 0)
    {
      Id iad = _npair * _sizeAux;
      if (aux_iech != nullptr)
        (void)memcpy(&_tabAuxIech[iad], aux_iech, _sizeAux);
      if (aux_jech != nullptr)
        (void)memcpy(&_tabAuxJech[iad], aux_jech, _sizeAux);
    }
    _npair++;
    return (0);
  }

} // namespace gstlrn
