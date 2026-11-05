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
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/AStringable.hpp"

namespace gstlrn
{
  /**
   * @brief Set of parameters used by the Simulation using FFT method
   *
   * The parameters are:
   * - flag_aliasing: this option allows subdiscretizing the lag for a better covariance rendition
   *   This option is now deprecated: a message is issued if the user switches it ON
   * - percent: this is quality control parameter, defined as a percentage (see Note)
   *
   * Note: the FFT method requires to fill a grid (similar to the target grid) with the covariance function.
   * In fact, the work is done in a Dilated Grid in order to avoid aliasing effects.
   * The extension is done until the covariance function is below a certain threshold,
   * defined as a percentage of the maximum value of the covariance function.
   * A value of 5 is considered as reasonable, 1 is a better choice to reproduce the covariance better.
   * It can even be set to a smaller number, inducing a larger dilated grid and therefore more memory consumption and CPU time.
   */
  class GSTLEARN_EXPORT SimuFFTParam: public AStringable
  {
  public:
    SimuFFTParam(bool flag_aliasing = false, double percent = 1.);
    SimuFFTParam(const SimuFFTParam& r);
    SimuFFTParam& operator=(const SimuFFTParam& r);
    virtual ~SimuFFTParam();

    /// Interface to AStringable
    String toString(const AStringFormat* strfmt = nullptr) const override;

    double getPercent() const { return _percent; }

    static void _checkAliasing(bool aliasing);

  private:
    bool _flagAliasing;
    double _percent;
  };
} // namespace gstlrn
