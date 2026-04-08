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

#include "Basic/NamingConvention.hpp"
#include "Basic/VectorNumT.hpp"
#include "Simulation/CalcSimuSpectral.hpp"

namespace gstlrn
{
  typedef struct
  {
    Id _k;
    Id _countP;
    Id _countM;
    std::map<Id, std::map<Id, Id>> _tab;
  } spSim;

  class ACov;
  class ModelGeneric;

  /**
   * Class for operating the Spectral simulations on S2
   */
  class GSTLEARN_EXPORT SimuSpectralS2: public CalcSimuSpectral
  {
  public:
    SimuSpectralS2(
      Id nbsimu = 1,
      Id ns = 10000,
      Id nd = 100,
      Id seed = 4324324,
      bool verbose = false);
    SimuSpectralS2(const SimuSpectralS2& r) = delete;
    SimuSpectralS2& operator=(const SimuSpectralS2& r) = delete;
    virtual ~SimuSpectralS2();

  protected:
    Id _simulate() override;
    Id _compute(
      Db* dbout,
      const VectorBool& activeArray,
      VectorVectorDouble& tab) override;

  private:
    static void _printSpSim(const spSim& spsim, Id status = 0);
    void _printSpSims(Id status = 0);

    double getPhi(Id i) { return _phi[i]; };

    static Id _getKey1Maximum(const spSim& spsim);
    static Id _getSumValue(const spSim& spsim);
    static VectorInt _getKeys1(const spSim& spsim);
    static VectorInt _getKeys2(const spSim& spsim, Id key1);
    static VectorInt _getValues2(const spSim& spsim, Id key1);

  private:
    // spectrum for S2
    std::vector<spSim> _spSims;
    VectorDouble _phi; // Vector length = _ns
  };

} // namespace gstlrn
