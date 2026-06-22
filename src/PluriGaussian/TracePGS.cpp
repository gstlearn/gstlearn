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
#include "PluriGaussian/TracePGS.hpp"
#include "Basic/Message.hpp"
#include "Basic/VectorT.hpp"
#include "Core/Keypair.hpp"
#include "geoslib_define.h"

namespace gstlrn
{
  TracePGS::TracePGS(Id ngrf, Id npar, bool flag_rho, bool flag_correl)
    : _flagTrace(false)
    , _nrow(0)
    , _ncol(0)
    , _trace()
  {
    define(ngrf, npar, flag_rho, flag_correl);
  }

  /****************************************************************************/
  /*!
   **  Extract the information of the Trace
   **
   ** \return The calculated score
   **
   *****************************************************************************/
  double TracePGS::extractTrace()
  {
    Id nrow = _nrow;
    Id ncol = _ncol;
    if (nrow <= 0 || ncol <= 0) return TEST;

    /* Evaluate the sum of the score */

    double totsum = 0.;
    for (Id irow = 0; irow < nrow; irow++)
    {
      totsum += _trace[ncol * irow + 2];
      if (ncol >= 5) totsum += _trace[ncol * irow + 4];
    }
    set_keypair("vario.pgs_score", 1, 1, 1, &totsum);

    if (_flagTrace)
      set_keypair("vario.pgs_stack", 1, nrow, ncol, _trace.data());

    return (totsum);
  }

  /****************************************************************************/
  /*!
   **  Add a new row to the trace
   **
   *****************************************************************************/
  void TracePGS::addRow()
  {
    if (!_flagTrace) return;
    Id nrow = _nrow;
    Id iad = _ncol * nrow;

    nrow++;
    _trace.resize(nrow * _ncol);
    for (Id icol = 0; icol < _ncol; icol++) _trace[iad + icol] = TEST;
    _nrow = nrow;
  }

  /****************************************************************************/
  /*!
   **  Update the Trace array
   **
   ** \param[in]  value0     First value in a Trace row
   ** \param[in]  value1     Second value in a Trace row
   ** \param[in]  origin     Origin for values in record (after 2 heading values)
   ** \param[in]  number     Number of values assigned
   ** \param[in]  values     Array of values assigned
   **
   *****************************************************************************/
  void TracePGS::update(
    double value0,
    double value1,
    Id origin,
    Id number,
    const double* values)
  {
    if (!_flagTrace) return;
    Id iad = _ncol * (_nrow - 1);
    if (2 + origin + number > _ncol)
      messageAbort(
        "Error in Trace dimension (ncol=%d origin=%d number=%d)", _ncol, origin,
        number);

    /* Store the information */

    _trace[iad] = value0;
    _trace[iad + 1] = value1;

    for (Id i = 0; i < number; i++) _trace[iad + 2 + origin + i] = values[i];
  }

  void TracePGS::define(Id ngrf, Id npar, Id flag_rho, bool flag_correl)
  {
    _flagTrace = !flag_rho;
    if (!_flagTrace) return;

    if (!flag_correl)
      _ncol = 2 + 2 * ngrf;
    else
      _ncol = 2 + 2 + npar;

    _nrow = 0;
    _trace.clear();
  }

} // namespace gstlrn
