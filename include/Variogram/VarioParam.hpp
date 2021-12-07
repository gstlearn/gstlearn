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
#include "Variogram/DirParam.hpp"
#include "Basic/Vector.hpp"
#include "Basic/IClonable.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"

class Db;
class Model;

class GSTLEARN_EXPORT VarioParam : public AStringable, public IClonable
{
public:
  VarioParam(double scale = 0.,
             VectorDouble dates = VectorDouble());
  VarioParam(const VarioParam& VarioParam,
             const VectorInt& dircols);
  VarioParam(const VarioParam& r);
  VarioParam& operator=(const VarioParam& r);
  virtual ~VarioParam();

public:
  virtual String toString(int level = 0) const override;
  virtual IClonable* clone() const override;

  void addDirs(const DirParam& dirparam);
  void addDirs(const std::vector<DirParam>& dirparams);
  void delDir(int rank);
  void delAllDirs();

  double getScale() const { return _scale; }

  int    getDateNumber() const { return (int) _dates.size() / 2; }
  int    getDirectionNumber() const { return (int) _dirparams.size(); }

  const VectorDouble& getDates() const { return _dates; }
  double getDates(int idate, int icas) const;

  int hasDate() const { return (getDateNumber() > 0 && (_dates[0] > -1.e30 || _dates[1] < 1.e30)); }

  void setScale(double scale) { _scale = scale; }

  int getLagNumber(int idir) const;
  VectorDouble getCodir(int idir = 0) const;

  void setDates(VectorDouble dates) { _dates = dates; }

  const std::vector<DirParam>& getDirParam() const { return _dirparams; }
  const DirParam& getDirParam(int idir) const { return _dirparams[idir]; }

  void setDPas(int idir,const Db* db);
  void setGrincr(int idir, const VectorInt& grincr);

  int getDimensionNumber() const { return _dirparams[0].getDimensionNumber(); }
  String toStringMain(int level) const;

private:
  int  _getAddress(int ivar, int jvar) const;
  bool _isVariableValid(int ivar) const;
  bool _isDirectionValid(int idir) const;
  bool _isBivariableValid(int i) const;
  bool _isDateValid(int idate) const;
  void _initMeans();
  void _initVars();
  VectorDouble _getDirectionInterval(int idir) const;

private:
  double                _scale;
  VectorDouble          _dates;
  std::vector<DirParam> _dirparams;
};
