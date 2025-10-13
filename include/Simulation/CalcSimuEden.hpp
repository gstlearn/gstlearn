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

#include "geoslib_define.h"

#include "Basic/AStringable.hpp"
#include "Db/DbGrid.hpp"
#include "Simulation/ACalcSimulation.hpp"
#include "Skin/ISkinFunctions.hpp"

namespace gstlrn
{
class Skin;
class MatrixDense;

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
  CalcSimuEden(Id nfacies  = 0,
               Id nfluids  = 0,
               Id niter    = 1,
               Id nbsimu   = 0,
               Id seed     = 4324324,
               bool verbose = false);
  CalcSimuEden(const CalcSimuEden& r)            = delete;
  CalcSimuEden& operator=(const CalcSimuEden& r) = delete;
  virtual ~CalcSimuEden();

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface to ISkinFunctions
  Id isAlreadyFilled(Id ipos) const override;
  Id isToBeFilled(Id ipos) const override;
  double getWeight(Id ipos, Id idir) const override;

  void setIndFacies(Id indFacies) { _indFacies = indFacies; }
  void setIndFluid(Id indFluid) { _indFluid = indFluid; }
  void setIndPerm(Id indPerm) { _indPerm = indPerm; }
  void setIndPoro(Id indPoro) { _indPoro = indPoro; }
  void setSpeeds(const VectorInt& speeds) { _speeds = speeds; }
  void setNMax(double numberMax) { _numberMax = numberMax; }
  void setShowFluid(bool showFluid) { _showFluid = showFluid; }
  void setVolumeMax(double volumeMax) { _volumeMax = volumeMax; }

private:
  bool _check() override;
  bool _preprocess() override;
  bool _run() override;
  bool _postprocess() override;
  void _rollback() override;

  bool _simulate();
  bool _fluid_check(void);
  Id _getWT(Id ifacies, Id ifluid, Id perm, Id idir);
  Id _getFACIES(Id iech) const;
  Id _getFLUID(Id iech) const;
  Id _getFLUID_OLD(Id iech) const;
  Id _getPERM(Id iech) const;
  double _getDATE(Id iech);
  double _getPORO(Id iech) const;
  void _setFLUID(Id iech, Id ifluid);
  void _setFACIES(Id iech, Id ifacies);
  void _setFACIES_CORK(Id iech);
  void _setDATE(Id iech, Id idate);
  void _printParams(bool verbose);
  void _statsDefine(void);
  void _statsReset();
  void _statsInit();
  void _setStatCount(Id ifacies, Id ifluid, Id value);
  void _setStatVolume(Id ifacies, Id ifluid, double value);
  void _addStatCount(Id ifacies, Id ifluid, Id value);
  void _addStatVolume(Id ifacies, Id ifluid, double value);
  void _checkInconsistency(bool verbose);
  Id _getStatCount(Id ifacies, Id ifluid) const;
  double _getStatVolume(Id ifacies, Id ifluid) const;
  Id _checkMax(double number_max, double volume_max);
  Id _fluidModify(Skin* skin, Id ipos, Id* ref_fluid_loc);
  void _statsPrint(const char* title);
  void _statsEmpty(const char* title);
  void _calculateCumul(void);
  void _updateResults(Id reset_facies, Id show_fluid);
  void _normalizeCumul(Id niter);
  Id _countAlreadyFilled() const;
  Id _countIsToBeFilled() const;

private:
  bool _verbose;
  /// 1 for modifying the value of the cells to show
  ///\li                       the initial valid fluid information
  ///\li                       the cork (different from shale)
  bool _showFluid;
  Id _iptrStatFluid;
  Id _iptrStatCork;
  Id _iptrFluid;
  Id _iptrDate;
  Id _niter;        /// Number of iterations
  Id _nfacies;      /// number of facies (facies 0 excluded)
  Id _nfluids;      /// number of fluids
  VectorInt _speeds; /// array containing the travel speeds
  double _numberMax; /// Maximum count of cells invaded (or TEST)
  double _volumeMax; /// Maximum volume invaded (or TEST)

  Id _indFacies; /// Rank of the variable containing the Facies
  Id _indFluid;  /// Rank of the variable containing the Fluid
  Id _indPerm;   /// Rank of the variable containing the Permeability
  Id _indPoro;   /// Rank of the variable containing the Porosity
  Id _indDate;

  Id _nxyz;
  Id _ncork;
  VectorInt _numbers;
  VectorDouble _volumes;
};

GSTLEARN_EXPORT Id fluid_propagation(DbGrid* dbgrid,
                                      const String& name_facies,
                                      const String& name_fluid,
                                      const String& name_perm,
                                      const String& name_poro,
                                      Id nfacies,
                                      Id nfluids,
                                      Id niter                       = 1,
                                      const VectorInt& speeds         = VectorInt(),
                                      bool show_fluid                 = false,
                                      double number_max               = TEST,
                                      double volume_max               = TEST,
                                      Id seed                        = 321321,
                                      bool verbose                    = false,
                                      const NamingConvention& namconv = NamingConvention("Eden"));
} // namespace gstlrn