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

#include "Variogram/VarioParam.hpp"
#include "Basic/Vector.hpp"
#include "Basic/IClonable.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"

int identifyCalculTypeC(const String& calcul_name);
int identifyFlagAsymC(const String& calcul_name);

class Db;
class Model;

class VarioC : public AStringable, public ASerializable, public IClonable
{
public:
  VarioC(bool flagSample = false, VectorDouble dates = VectorDouble());
  VarioC(const VarioC& VarioC,
         const VectorInt& varcols,
         const VectorInt& dircols);
  VarioC(const String& neutralFileName, bool verbose);
  VarioC(const VarioC& r);
  VarioC& operator=(const VarioC& r);
  virtual ~VarioC();

public:
  virtual String toString(int level = 0) const override;
  int deSerialize(const String& filename, bool verbose = false) override;
  int serialize(const String& filename, bool verbose = false) const override;
  virtual IClonable* clone() const override;

  int getDirectionNumber() const { return _varioparam.getDirectionNumber(); }

  int    getVariableNumber() const { return _nVar; }
  const  VectorDouble& getMeans() const { return _means; }
  double getMeans(int ivar) const;
  const  VectorDouble& getVars()  const { return _vars; }
  double getVars(int ivar, int jvar) const;
  double getVars(int ijvar) const;

  void setMeans(const VectorDouble& means);
  void setMeans(int ivar, double mean);
  void setVars(const VectorDouble& vars);
  void setVars(int ijvar, double value);
  void setVars(int ivar, int jvar, double value);

  int getDirSize(int idir) const;

  double getGg(int idir, int i) const;
  double getHh(int idir, int i) const;
  double getSw(int idir, int i) const;
  double getUtilize(int idir, int i) const;

  double getGg(int idir, int ivar, int jvar, int ipas) const;
  double getHh(int idir, int ivar, int jvar, int ipas) const;
  double getSw(int idir, int ivar, int jvar, int ipas) const;
  double getUtilize(int idir, int ivar, int jvar, int ipas) const;

  VectorDouble getGgVec(int idir, int ivar, int jvar) const;
  VectorDouble getHhVec(int idir, int ivar, int jvar) const;
  VectorDouble getSwVec(int idir, int ivar, int jvar) const;
  VectorDouble getUtilizeVec(int idir, int ivar, int jvar) const;

  const VectorDouble& getGgVec(int idir) const;
  const VectorDouble& getHhVec(int idir) const;
  const VectorDouble& getSwVec(int idir) const;
  const VectorDouble& getUtilizeVec(int idir) const;

  void setGg(int idir, int i, double gg);
  void setHh(int idir, int i, double hh);
  void setSw(int idir, int i, double sw);
  void setUtilize(int idir, int i, double utilize);

  void setSw(int idir, int ivar, int jvar, int ipas, double sw);
  void setHh(int idir, int ivar, int jvar, int ipas, double hh);
  void setGg(int idir, int ivar, int jvar, int ipas, double gg);
  void setUtilize(int idir, int ivar, int jvar, int ipas, double utilize);

  void updSw(int idir, int i, double sw);
  void updHh(int idir, int i, double hh);
  void updGg(int idir, int i, double gg);
  void updUtilize(int idir, int i, double utilize);

  int getCenter(int ivar = 0, int jvar = 0, int idir = 0) const;

  void directionResize(int idir);
  void internalResize(int nvar);
  int  getCalculType() const;
  int  getFlagAsym() const;

  double getHmax(int ivar=-1, int jvar=-1, int idir=-1) const;
  VectorDouble getHRange(int ivar=-1, int jvar=-1, int idir=-1) const;
  double getGmax(int ivar = -1,
                 int jvar = -1,
                 int idir = -1,
                 bool flagAbs = false,
                 bool flagSill = false) const;
  VectorDouble getGRange(int ivar=-1, int jvar=-1, int idir=-1, bool flagSill = false) const;

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

  void patchCenter(int idir, int nech, double rho);

  int fill(int idir,
           const VectorDouble& sw,
           const VectorDouble& gg,
           const VectorDouble& hh);

private:
  int _getVarAddress(int ivar, int jvar) const;
  int _getDirAddress(int idir,
                     int ivar,
                     int jvar,
                     int ipas,
                     bool flag_abs,
                     int sens) const;
  bool _isVariableValid(int ivar) const;
  bool _isDirectionValid(int idir) const;
  bool _isBivariableValid(int ijvar) const;
  bool _isAddressValid(int i, int idir) const;
  void _initMeans();
  void _initVars();
  VectorDouble _getVariableInterval(int ivar) const;
  VectorDouble _getDirectionInterval(int idir) const;

private:
  int _nVar;
  VarioParam   _varioparam;
  VectorDouble _means;
  VectorDouble _vars;
  std::vector<VectorDouble> _sw;      /* Array for number of lags */
  std::vector<VectorDouble> _gg;      /* Array for average variogram values */
  std::vector<VectorDouble> _hh;      /* Array for average distance values */
  std::vector<VectorDouble> _utilize; /* Array to mention if a lag is used or not */
};
