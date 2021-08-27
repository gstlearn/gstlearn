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

#include "Variogram/DirParam.hpp"
#include "Basic/Vector.hpp"
#include "Basic/IClonable.hpp"
#include "Basic/AStringable.hpp"

int identifyCalculTypeN(const String& calcul_name);
int identifyFlagAsymN(const String& calcul_name);

class Db;
class Model;

class VarioParam : public AStringable, public IClonable
{
public:
  VarioParam(double scale = 0.,
             bool flagSample = false,
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
  void addDirs(const std::vector<DirParam> dirparams);
  void delDir(int rank);
  void delAllDirs();

  const  std::string& getCalculName() const { return _calculName; }
  int    getCalculType() const { return identifyCalculTypeN(_calculName); }
  int    getFlagAsym() const { return identifyFlagAsymN(_calculName); }
  int    getFlagSample() const { return _flagSample; }
  double getScale() const { return _scale; }

  int    getDateNumber() const { return (int) _dates.size() / 2; }
  int    getDimensionNumber() const { return _nDim; }
  int    getDirectionNumber() const { return (int) _dirparams.size(); }

  const VectorDouble& getDates() const { return _dates; }
  double getDates(int idate, int icas) const;

  int hasDate() const { return (getDateNumber() > 0 && (_dates[0] > -1.e30 || _dates[1] < 1.e30)); }

  void setCalculName(std::string calcul_name) { _calculName = calcul_name; }
  void setFlagSample(int flag_sample) { _flagSample = flag_sample; }
  void setDimensionNumber(int ndim) { _nDim = ndim; }
  void setScale(double scale) { _scale = scale; }

  int getLagNumber(int idir) const;
  int getLagTotalNumber(int idir) const;
  VectorDouble getCodir(int idir = 0) const;

  void setDates(VectorDouble dates) { _dates = dates; }

  const std::vector<DirParam>& getDirParam() const { return _dirparams; }
  const DirParam& getDirParam(int idir) const { return _dirparams[idir]; }

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
  String       _calculName;
  int          _nDim;
  bool         _flagSample;
  double       _scale;
  VectorDouble _dates;
  std::vector<DirParam> _dirparams;
};
