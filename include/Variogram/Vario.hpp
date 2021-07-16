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

#include "Variogram/Dir.hpp"
#include "Basic/Vector.hpp"
#include "Basic/IClonable.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"

int identifyCalculType(const String& calcul_name);
int identifyFlagAsym(const String& calcul_name);

class Db;
class Model;

class Vario : public AStringable, ASerializable, IClonable
{
public:
  Vario(double scale = 0.,
        bool flagSample = false,
        VectorDouble dates = VectorDouble());
  Vario(const Vario& vario, const VectorInt& varcols, const VectorInt& dircols, bool flagVario);
  Vario(const String& neutralFileName, bool verbose);
  Vario(const Vario& r);
  Vario& operator=(const Vario& r);
  virtual ~Vario();

public:
  virtual String toString(int level = 0) const override;
  int deSerialize(const String& filename, bool verbose = false) override;
  int serialize(const String& filename, bool verbose = false) override;
  virtual IClonable* clone() const override;

  void addDirs(const Dir& dir);
  void addDirs(const std::vector<Dir> dirs);
  void delDir(int rank);
  void delAllDirs();

  const  std::string& getCalculName() const { return _calculName; }
  int    getCalculType() const { return identifyCalculType(_calculName); }
  int    getFlagAsym() const { return identifyFlagAsym(_calculName); }
  int    getFlagSample() const { return _flagSample; }

  int    getDateNumber() const { return (int) _dates.size() / 2; }
  int    getDimensionNumber() const { return _nDim; }
  int    getDirectionNumber() const { return (int) _dirs.size(); }
  int    getVariableNumber() const { return _nVar; }
  double getScale() const { return _scale; }

  const VectorDouble& getDates() const { return _dates; }
  double getDates(int idate, int icas) const;
  const VectorDouble& getMeans() const { return _means; }
  double getMeans(int ivar) const;
  const VectorDouble& getVars()  const { return _vars; }
  double getVars(int ivar, int jvar) const;
  double getVars(int i) const;

  int hasDate() const { return (getDateNumber() > 0 && (_dates[0] > -1.e30 || _dates[1] < 1.e30)); }

  void setCalculName(std::string calcul_name) { _calculName = calcul_name; }
  void setScale(double scale) { _scale = scale; }
  void setFlagSample(int flag_sample) { _flagSample = flag_sample; }
  void setMeans(const VectorDouble& means);
  void setMeans(int ivar, double mean);
  void setVars(const VectorDouble& vars);
  void setVars(int i, double value);
  void setVars(int ivar, int jvar, double value);

  int getLagNumber(int idir) const;
  VectorDouble getGg(int ivar = 0, int jvar = 0, int idir = 0) const;
  VectorDouble getHh(int ivar = 0, int jvar = 0, int idir = 0) const;
  VectorDouble getSw(int ivar = 0, int jvar = 0, int idir = 0) const;
  VectorDouble getCodir(int idir = 0) const;

  void setDates(VectorDouble dates) { _dates = dates; }

  const std::vector<Dir>& getDirs() const { return _dirs; }
  const Dir& getDirs(int idir) const { return _dirs[idir]; }

  void internalResize(int ndim, int nvar, const String& calculName);

  double getHmax(int ivar, int jvar=0) const;
  double getHmax() const;
  double getGmax(int ivar, int jvar=0, bool flagAbs=false) const;
  double getGmax(bool flagAbs=false) const;

  int compute(Db *db,
              const String& calculName = "vg",
              const VectorDouble& means = VectorDouble(),
              const VectorDouble& vars = VectorDouble(),
              bool flag_grid = false,
              bool flag_gen = false,
              bool flag_sample = false,
              bool verr_mode = false,
              bool flag_model = false,
              Model *model = nullptr,
              bool verbose = false);
  int computeIndic(Db *db,
                   const String& calculName = "vg",
                   const VectorDouble& means = VectorDouble(),
                   const VectorDouble& vars = VectorDouble(),
                   bool flag_grid = false,
                   bool flag_gen = false,
                   bool flag_sample = false,
                   bool verr_mode = false,
                   bool flag_model = false,
                   Model *model = nullptr,
                   bool verbose = false,
                   int nfacmax = 10);
  bool isCalculated() const;

  // Pipe the Dir methods for setting values
  void setSw(int idir, int iad, double sw) { _dirs[idir].setSw(iad, sw); }
  void setHh(int idir, int iad, double hh) { _dirs[idir].setHh(iad, hh); }
  void setGg(int idir, int iad, double gg) { _dirs[idir].setGg(iad, gg); }
  void setSw(int idir, int ivar, int jvar, int ipas, double sw) { return _dirs[idir].setSw(ivar, jvar, ipas, sw); }
  void setHh(int idir, int ivar, int jvar, int ipas, double hh) { return _dirs[idir].setHh(ivar, jvar, ipas, hh); }
  void setGg(int idir, int ivar, int jvar, int ipas, double gg) { return _dirs[idir].setGg(ivar, jvar, ipas, gg); }
  void updSw(int idir, int iad, double sw) { return _dirs[idir].updSw(iad, sw); }
  void updHh(int idir, int iad, double hh) { return _dirs[idir].updHh(iad, hh); }
  void updGg(int idir, int iad, double gg) { return _dirs[idir].updGg(iad, gg); }
  void setUtilize(int idir, int iad, double val) { return _dirs[idir].setUtilize(iad, val); }
  void clean(int idir) { return _dirs[idir].clean(); }

private:
  int  _getAddress(int ivar, int jvar) const;
  bool _isVariableValid(int ivar) const;
  bool _isDirectionValid(int idir) const;
  bool _isBivariableValid(int i) const;
  bool _isDateValid(int idate) const;
  void _initMeans();
  void _initVars();

private:
  String _calculName;
  int _nDim;
  int _nVar;
  bool _flagSample;
  double _scale;
  VectorDouble _means;
  VectorDouble _vars;
  VectorDouble _dates;
  std::vector<Dir> _dirs;
};
