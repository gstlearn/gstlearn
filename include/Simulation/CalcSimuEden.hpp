/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
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

/**
 * Multivariate multiphase propagation into a set of components
 * constrained by initial conditions and fluid densities
 *
 * \remark  Directions are ordered as follows :
 * \remark  0: +X; 1: -X; 2: +Y; 3: -Y; 4: +Z(up); 5: -Z(down)
 * \remark  The coding of the matrix is:
 * \remark              facies + nfacies * fluid
 * \remark  Facies: 0 (Shale), 1 to nfacies, -1 (Cork)
 * \remark  Fluids: 0 (undefined), 1 to nfluids, -1 (No Fluid)
 * \remark  Fluids should be ordered by increasing weight
 * \remark  A Permeability variable is a value (>=1) which divides
 * \remark  the velocities. This variable is optional.
 * \remark  A Porosity variable is a value (in [0,1]) which multiplies
 * \remark  the volumes. This variable is optional.
 * \remark  Volume_max represents the volumic part of the invaded area:
 * \remark  it is always <= number of cells invaded.
 */
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
  /// 1 for modifying the value of the cells to show
  ///\li                       the initial valid fluid information
  ///\li                       the cork (different from shale)
  bool _showFluid;
  int _iptrStatFluid;
  int _iptrStatCork;
  int _iptrFluid;
  int _iptrDate;
  int _niter; /// Number of iterations
  int _nfacies; /// number of facies (facies 0 excluded)
  int _nfluids; /// number of fluids
  VectorInt _speeds; /// array containing the travel speeds
  double _numberMax; /// Maximum count of cells invaded (or TEST)
  double _volumeMax; /// Maximum volume invaded (or TEST)

  int _indFacies; /// Rank of the variable containing the Facies
  int _indFluid; /// Rank of the variable containing the Fluid
  int _indPerm; /// Rank of the variable containing the Permeability
  int _indPoro; /// Rank of the variable containing the Porosity
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
