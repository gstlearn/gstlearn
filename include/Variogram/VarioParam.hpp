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

#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Faults/Faults.hpp"
#include "Variogram/DirParam.hpp"

namespace gstlrn
{
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
class GSTLEARN_EXPORT VarioParam: public AStringable, public ICloneable
{
public:
  VarioParam(double scale              = 0.,
             const VectorDouble& dates = VectorDouble(),
             const Faults* faults      = nullptr);
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
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Shortcuts
  static VarioParam* createOmniDirection(Id nlag                      = 10,
                                         double dlag                  = 1.,
                                         double toldis                = 0.5,
                                         Id opt_code                  = 0,
                                         Id idate                     = 0,
                                         double bench                 = TEST,
                                         double cylrad                = TEST,
                                         double tolcode               = 0.,
                                         const VectorDouble& breaks   = VectorDouble(),
                                         double scale                 = 0.,
                                         const VectorDouble& dates    = VectorDouble(),
                                         const ASpaceSharedPtr& space = ASpaceSharedPtr());
  static VarioParam* createMultiple(Id ndir,
                                    Id nlag                      = 10,
                                    double dlag                  = 1.,
                                    double toldis                = 0.5,
                                    double angref                = 0.,
                                    double scale                 = 0.,
                                    const VectorDouble& dates    = VectorDouble(),
                                    const ASpaceSharedPtr& space = ASpaceSharedPtr());
  static VarioParam* createMultipleFromGrid(const DbGrid* dbgrid,
                                            Id nlag,
                                            double scale                 = 0.,
                                            const VectorDouble& dates    = VectorDouble(),
                                            const ASpaceSharedPtr& space = ASpaceSharedPtr(),
                                            Id ndimax                    = 0);
  static VarioParam* createFromSpaceDimension(Id nlag                      = 10,
                                              double dlag                  = 1.,
                                              double toldis                = 0.5,
                                              double tolang                = 45.,
                                              double scale                 = 0.,
                                              const VectorDouble& dates    = VectorDouble(),
                                              const ASpaceSharedPtr& space = ASpaceSharedPtr());
  static VarioParam* createSeveral2D(const VectorDouble& angles,
                                     Id nlag                      = 10,
                                     double dlag                  = 1.,
                                     double toldis                = 0.5,
                                     double tolang                = TEST,
                                     double scale                 = 0.,
                                     const VectorDouble& dates    = VectorDouble(),
                                     const ASpaceSharedPtr& space = ASpaceSharedPtr());

  void addDir(const DirParam& dirparam);
  void addMultiDirs(const std::vector<DirParam>& dirparams);
  void delDir(Id rank);
  void delAllDirs();

  ASpaceSharedPtr getSpace() const { return _dirparams[0].getSpace(); }
  double getScale() const { return _scale; }
  Id getNDate() const { return static_cast<Id>(_dates.size() / 2); }
  Id getNDir() const { return static_cast<Id>(_dirparams.size()); }
  const VectorDouble& getDates() const { return _dates; }
  double getDate(Id idate, Id icas) const;
  Id getNLag(Id idir) const;
  VectorDouble getCodirs(Id idir = 0) const;
  const std::vector<DirParam>& getDirParams() const { return _dirparams; }
  const DirParam& getDirParam(Id idir) const { return _dirparams[idir]; }
  Id getNDim() const;
  bool isDefinedForGrid() const;

  Id hasDate() const { return (getNDate() > 0 && (_dates[0] > MINIMUM_BIG || _dates[1] < MAXIMUM_BIG)); }
  bool isDateUsed(const Db* db1, const Db* db2 = nullptr) const;

  void setScale(double scale) { _scale = scale; }
  void setDates(const VectorDouble& dates) { _dates = dates; }
  void setDPas(Id idir, const DbGrid* db);
  void setGrincr(Id idir, const VectorInt& grincr);

  String toStringMain(const AStringFormat* strfmt = nullptr) const;

  const Faults* getFaults() const { return _faults; }
  bool hasFaults() const { return _faults != nullptr; }
  void addFaults(const Faults* faults) { _faults = faults; }

private:
  Id _getAddress(Id ivar, Id jvar) const;
  bool _isVariableValid(Id ivar) const;
  bool _isDirectionValid(Id idir) const;
  bool _isBivariableValid(Id i) const;
  bool _isDateValid(Id idate) const;
  void _initMeans();
  void _initVars();
  VectorDouble _getDirectionInterval(Id idir) const;
  bool _validDefinedFromGrid(const DirParam& dirparam) const;

private:
  double _scale;
  VectorDouble _dates;
  std::vector<DirParam> _dirparams;
  const Faults* _faults; // Pointer copy (not to be deleted)
};

GSTLEARN_EXPORT Db* buildDbFromVarioParam(Db* db, const VarioParam& varioparam);
} // namespace gstlrn