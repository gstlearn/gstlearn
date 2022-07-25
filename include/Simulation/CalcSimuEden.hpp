/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "geoslib_define.h"

#include "Db/DbGrid.hpp"
#include "Skin/ISkinFunctions.hpp"
#include "Basic/AStringable.hpp"
#include "Simulation/ACalcSimulation.hpp"

class Skin;
class MatrixRectangular;

class GSTLEARN_EXPORT CalcSimuEden: public ACalcSimulation, public AStringable, public ISkinFunctions
{
public:
  CalcSimuEden(int nfacies = 0,
           int nfluids = 0,
           int niter = 1,
           int nbsimu = 0,
           int seed = 4324324,
           bool verbose = false);
  CalcSimuEden(const CalcSimuEden &r) = delete;
  CalcSimuEden& operator=(const CalcSimuEden &r) = delete;
  virtual ~CalcSimuEden();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface to ISkinFunctions
  int    isAlreadyFilled(int ipos) const override;
  int    isToBeFilled(int ipos) const override;
  double getWeight(int ipos, int idir) const override;

  void setIndFacies(int indFacies) { _indFacies = indFacies; }
  void setIndFluid(int indFluid) { _indFluid = indFluid; }
  void setIndPerm(int indPerm) { _indPerm = indPerm; }
  void setIndPoro(int indPoro) { _indPoro = indPoro; }
  void setSpeeds(const VectorInt &speeds) { _speeds = speeds; }
  void setNumberMax(double numberMax) { _numberMax = numberMax; }
  void setShowFluid(bool showFluid) { _showFluid = showFluid; }
  void setVolumeMax(double volumeMax) { _volumeMax = volumeMax; }

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;

  bool _simulate();
  bool _fluid_check(void);
  int  _getWT(int ifacies, int ifluid, int perm, int idir);
  int  _getFACIES(int iech) const;
  int  _getFLUID(int iech) const;
  int  _getFLUID_OLD(int iech) const;
  int  _getPERM(int iech) const;
  double _getDATE(int iech);
  double _getPORO(int iech) const;
  void _setFLUID(int iech, int ifluid);
  void _setFACIES(int iech, int ifacies);
  void _setFACIES_CORK(int iech);
  void _setDATE(int iech, int idate);
  void _printParams(bool verbose);
  void _statsDefine(void);
  void _statsReset();
  void _statsInit();
  void _setStatNumber(int ifacies, int ifluid, int value);
  void _setStatVolume(int ifacies, int ifluid, double value);
  void _addStatNumber(int ifacies, int ifluid, int value);
  void _addStatVolume(int ifacies, int ifluid, double value);
  void _checkInconsistency(bool verbose);
  int _getStatNumber(int ifacies, int ifluid) const;
  double _getStatVolume(int ifacies, int ifluid) const;
  int _checkMax(double number_max, double volume_max);
  int _fluidModify(Skin *skin, int ipos, int *ref_fluid_loc);
  void _statsPrint(const char *title);
  void _statsEmpty(const char *title);
  void _calculateCumul(void);
  void _updateResults(int reset_facies, int show_fluid);
  void _normalizeCumul(int niter);
  int  _countAlreadyFilled() const;
  int  _countIsToBeFilled() const;

private:
  bool _verbose;
  bool _showFluid;
  int _iptrStatFluid;
  int _iptrStatCork;
  int _iptrFluid;
  int _iptrDate;
  int _niter;
  int _nfacies;
  int _nfluids;
  VectorInt _speeds;
  double _numberMax;
  double _volumeMax;

  int _indFacies;
  int _indFluid;
  int _indPerm;
  int _indPoro;
  int _indDate;

  int _nxyz;
  int _ncork;
  VectorInt    _numbers;
  VectorDouble _volumes;
};

GSTLEARN_EXPORT int fluid_propagation(DbGrid *dbgrid,
                                      const String& name_facies,
                                      const String& name_fluid,
                                      const String& name_perm,
                                      const String& name_poro,
                                      int nfacies,
                                      int nfluids,
                                      int niter = 1,
                                      const VectorInt& speeds = VectorInt(),
                                      bool show_fluid = false,
                                      double number_max = TEST,
                                      double volume_max = TEST,
                                      int seed = 321321,
                                      bool verbose = false,
                                      const NamingConvention& namconv = NamingConvention("Eden"));
