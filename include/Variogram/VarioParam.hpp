/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Variogram/DirParam.hpp"
#include "Faults/Faults.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"

class Db;
class Model;

/**
 * Experimental Variogram calculation parameters TODO : to be improved
 */
// TODO : Inherits from ASpaceParam which inherits from ASPaceObject and AParam, which inherits from ASerializable, AStringable, IClonable
class GSTLEARN_EXPORT VarioParam : public AStringable, public ICloneable
{
public:
  VarioParam(double scale = 0.,
             const VectorDouble& dates = VectorDouble(),
             const Faults* faults = nullptr);
  VarioParam(const VarioParam& VarioParam,
             const VectorInt& rankdirs,
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
  static VarioParam* createOmniDirection(int npas = 10,
                                         double dpas = 1.,
                                         double toldis = 0.5,
                                         int opt_code = 0,
                                         int idate = 0,
                                         double bench = TEST,
                                         double cylrad = TEST,
                                         double tolcode = 0.,
                                         const VectorDouble& breaks = VectorDouble(),
                                         double scale = 0.,
                                         const VectorDouble& dates = VectorDouble(),
                                         const ASpace* space = nullptr);
  static VarioParam* createMultiple(int ndir,
                                    int npas = 10,
                                    double dpas = 1.,
                                    double toldis = 0.5,
                                    double scale = 0.,
                                    const VectorDouble& dates = VectorDouble(),
                                    const ASpace* space = nullptr);
  static VarioParam* createMultipleFromGrid(int npas,
                                            double scale = 0.,
                                            const VectorDouble& dates = VectorDouble(),
                                            const ASpace* space = nullptr);
  static VarioParam* createFromSpaceDimension(int npas = 10,
                                              double dpas = 1.,
                                              double toldis = 0.5,
                                              double tolang = 45.,
                                              double scale = 0.,
                                              const VectorDouble &dates = VectorDouble(),
                                              const ASpace *space = nullptr);

  void addDir(const DirParam& dirparam);
  void addMultiDirs(const std::vector<DirParam>& dirparams);
  void delDir(int rank);
  void delAllDirs();

  double getScale() const { return _scale; }
  int    getDateNumber() const { return (int) _dates.size() / 2; }
  int    getDirectionNumber() const { return (int) _dirparams.size(); }
  const VectorDouble& getDates() const { return _dates; }
  double getDate(int idate, int icas) const;
  int getLagNumber(int idir) const;
  VectorDouble getCodirs(int idir = 0) const;
  const std::vector<DirParam>& getDirParams() const { return _dirparams; }
  const DirParam& getDirParam(int idir) const { return _dirparams[idir]; }
  int getDimensionNumber() const;
  bool isDefinedForGrid() const;

  int hasDate() const { return (getDateNumber() > 0 && (_dates[0] > -1.e30 || _dates[1] < 1.e30)); }

  void setScale(double scale) { _scale = scale; }
  void setDates(VectorDouble dates) { _dates = dates; }
  void setDPas(int idir,const DbGrid* db);
  void setGrincr(int idir, const VectorInt& grincr);

  String toStringMain(const AStringFormat* strfmt) const;

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
  const Faults*         _faults;
};
