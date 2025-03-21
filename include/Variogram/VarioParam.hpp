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

#include "Space/ASpace.hpp"
#include "gstlearn_export.hpp"

#include "Variogram/DirParam.hpp"
#include "Faults/Faults.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/AStringable.hpp"

class Db;
class Model;

/**
 * \brief
 * Class containing the definition of the criteria for calculating the Spatial (and Temporal) Characteristics
 * from samples contained in a Db.
 *
 * These criteria consist in:
 * - some criteria based on the **dates**: this information will is used for calculating the Temporal Characteristics
 * - a collection of definitions of **Calculation Directions** for Spatial Characteristics.
 * For more information on a Direction definition, please refer to DirParam.hpp
 *
 * Note that this class also stores a pointer to any Faults definition, if to be used during the
 * calculation of the Spatial Characteristics.
 */
class GSTLEARN_EXPORT VarioParam : public AStringable, public ICloneable
{
public:
  VarioParam(double scale = 0.,
             const VectorDouble& dates = VectorDouble(),
             const Faults* faults = nullptr);
  VarioParam(const VarioParam& VarioParam,
             const VectorInt& dircols,
             const Faults* faults = nullptr);
  VarioParam(const VarioParam& r);
  VarioParam& operator=(const VarioParam& r);
  virtual ~VarioParam();

public:
  /// ICloneable interface
  IMPLEMENT_CLONING(VarioParam)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Shortcuts
  static VarioParam* createOmniDirection(int nlag = 10,
                                         double dlag = 1.,
                                         double toldis = 0.5,
                                         int opt_code = 0,
                                         int idate = 0,
                                         double bench = TEST,
                                         double cylrad = TEST,
                                         double tolcode = 0.,
                                         const VectorDouble& breaks = VectorDouble(),
                                         double scale = 0.,
                                         const VectorDouble& dates = VectorDouble(),
                                         const ASpaceSharedPtr& space = ASpaceSharedPtr());
  static VarioParam* createMultiple(int ndir,
                                    int nlag = 10,
                                    double dlag = 1.,
                                    double toldis = 0.5,
                                    double angref = 0.,
                                    double scale = 0.,
                                    const VectorDouble& dates = VectorDouble(),
                                    const ASpaceSharedPtr& space = ASpaceSharedPtr());
  static VarioParam* createMultipleFromGrid(const DbGrid* dbgrid,
                         int nlag,
                         double scale              = 0.,
                         const VectorDouble& dates = VectorDouble(),
                         const ASpaceSharedPtr& space       = ASpaceSharedPtr(),
                         int ndimax = 0);
  static VarioParam* createFromSpaceDimension(int nlag = 10,
                                              double dlag = 1.,
                                              double toldis = 0.5,
                                              double tolang = 45.,
                                              double scale = 0.,
                                              const VectorDouble &dates = VectorDouble(),
                                              const ASpaceSharedPtr& space = ASpaceSharedPtr());
  static VarioParam* createSeveral2D(const VectorDouble &angles,
                                     int nlag = 10,
                                     double dlag = 1.,
                                     double toldis = 0.5,
                                     double tolang = TEST,
                                     double scale = 0.,
                                     const VectorDouble& dates = VectorDouble(),
                                     const ASpaceSharedPtr& space = ASpaceSharedPtr());

  void addDir(const DirParam& dirparam);
  void addMultiDirs(const std::vector<DirParam>& dirparams);
  void delDir(int rank);
  void delAllDirs();

  ASpaceSharedPtr getSpace() const { return _dirparams[0].getSpace(); }
  double getScale() const { return _scale; }
  int    getNDate() const { return (int) _dates.size() / 2; }
  int    getNDir() const { return (int) _dirparams.size(); }
  const VectorDouble& getDates() const { return _dates; }
  double getDate(int idate, int icas) const;
  int getNLag(int idir) const;
  VectorDouble getCodirs(int idir = 0) const;
  const std::vector<DirParam>& getDirParams() const { return _dirparams; }
  const DirParam& getDirParam(int idir) const { return _dirparams[idir]; }
  int getNDim() const;
  bool isDefinedForGrid() const;

  int hasDate() const { return (getNDate() > 0 && (_dates[0] > -1.e30 || _dates[1] < 1.e30)); }
  bool isDateUsed(const Db *db1, const Db *db2 = nullptr) const;

  void setScale(double scale) { _scale = scale; }
  void setDates(const VectorDouble& dates) { _dates = dates; }
  void setDPas(int idir,const DbGrid* db);
  void setGrincr(int idir, const VectorInt& grincr);

  String toStringMain(const AStringFormat* strfmt = nullptr) const;

  const Faults* getFaults() const { return _faults; }
  bool hasFaults() const { return _faults != nullptr; }
  void addFaults(const Faults* faults) { _faults = faults; }

private:
  int  _getAddress(int ivar, int jvar) const;
  bool _isVariableValid(int ivar) const;
  bool _isDirectionValid(int idir) const;
  bool _isBivariableValid(int i) const;
  bool _isDateValid(int idate) const;
  void _initMeans();
  void _initVars();
  VectorDouble _getDirectionInterval(int idir) const;
  bool _validDefinedFromGrid(const DirParam& dirparam) const;

private:
  double                _scale;
  VectorDouble          _dates;
  std::vector<DirParam> _dirparams;
  const Faults*         _faults; // Pointer copy (not to be deleted)
};

GSTLEARN_EXPORT Db* buildDbFromVarioParam(Db *db, const VarioParam& varioparam);
