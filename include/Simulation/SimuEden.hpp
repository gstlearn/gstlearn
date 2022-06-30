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

#include "Simulation/ASimulation.hpp"
#include "Db/DbGrid.hpp"
#include "Skin/ISkinFunctions.hpp"
#include "Basic/AStringable.hpp"

class Skin;
class MatrixRectangular;

class GSTLEARN_EXPORT SimuEden: public ASimulation, public AStringable, public ISkinFunctions
{
public:
  SimuEden(int nbsimu = 0, int seed = 4324324);
  SimuEden(const SimuEden &r);
  SimuEden& operator=(const SimuEden &r);
  virtual ~SimuEden();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface to ISkinFunctions
  int    isAlreadyFilled(int ipos) const override;
  int    isToBeFilled(int ipos) const override;
  double getWeight(int ipos, int idir) const override;

  int simulate(DbGrid *dbgrid,
               int ind_facies,
               int ind_fluid,
               int ind_perm,
               int ind_poro,
               int nfacies,
               int nfluids,
               int niter,
               int iptr_fluid,
               int iptr_date,
               int iptr_stat_fluid,
               int iptr_stat_cork,
               const VectorInt& speeds,
               bool verbose = false,
               bool show_fluid = false,
               double number_max = TEST,
               double volume_max = TEST);
  int getTimeInterval(double date, int ntime, double time0, double dtime);
  MatrixRectangular fluidExtract(DbGrid* dbgrid,
                                 int ind_facies,
                                 int ind_fluid,
                                 int ind_poro,
                                 int ind_date,
                                 int nfacies,
                                 int nfluids,
                                 int facies0,
                                 int fluid0,
                                 int ntime,
                                 double time0,
                                 double dtime,
                                 bool verbose = false);

private:
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
  int _nxyz;
  int _nfluids;
  int _nfacies;
  int _indFacies;
  int _indFluid;
  int _indPerm;
  int _indPoro;
  int _indDate;
  int _iptrFluid;
  int _iptrStatFluid;
  int _iptrStatCork;
  int _iptrDate;
  DbGrid* _dbgrid;
  VectorInt _speeds;

  int _ncork;
  VectorInt _number;
  VectorDouble _volume;
};
