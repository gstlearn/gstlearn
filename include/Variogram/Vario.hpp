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

// Enums
#include "Variogram/ECalcVario.hpp"

#include "Variogram/VarioParam.hpp"

#include "Basic/IClonable.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"

class Db;
class Model;
class DirParam;

/**
 * Experimental Variogram (not only): TODO : to be improved
 */
class GSTLEARN_EXPORT Vario : public AStringable, public ASerializable, public IClonable
{
public:
  Vario(const VarioParam* varioparam,
        Db* db = nullptr,
        const VectorDouble& means = VectorDouble(),
        const VectorDouble& vars = VectorDouble());
  Vario(const Vario& r);
  Vario& operator=(const Vario& r);
  virtual ~Vario();

public:
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;
  virtual IClonable* clone() const override { return new Vario(*this); };

  int dumpToNF(const String& neutralFilename, bool verbose = false) const;
  int dumpToNF2(const String& neutralFilename, bool verbose = false) const;
  static Vario* create(const VarioParam* varioparam,
                       Db* db = nullptr,
                       const VectorDouble& means = VectorDouble(),
                       const VectorDouble& vars = VectorDouble());
  static Vario* createFromNF(const String& neutralFilename, bool verbose = false);
  static Vario* createFromNF2(const String& neutralFilename, bool verbose = false);

  void reduce(const VectorInt& varcols,
              const VectorInt& dircols,
              bool asSymmetric = false);

  const ECalcVario& getCalcul() const { return _calcul; }
  ECalcVario    getCalculType(const String& calcul_name) const;
  bool          getFlagAsym() const { return _flagAsym; }

  int    getVariableNumber() const { return _nVar; }
  const  VectorDouble& getMeans() const { return _means; }
  double getMean(int ivar) const;

  double getVar(int ivar, int jvar) const;
  double getVarIndex(int ijvar) const;
  const  VectorDouble& getVars()  const { return _vars; }
  void setMeans(const VectorDouble& means);
  void setMean(int ivar, double mean);
  void setVar(int ivar, int jvar, double value);
  void setVars(const VectorDouble& vars);
  void setVarIndex(int ijvar, double value);

  int getDirSize(int idir) const;

  double getGgByIndex(int idir, int i) const;
  double getHhByIndex(int idir, int i) const;
  double getSwByIndex(int idir, int i) const;
  double getUtilizeByIndex(int idir, int i) const;

  double getGg(int idir,
               int ivar,
               int jvar,
               int ipas,
               bool asCov = false,
               bool flagNormalized = false) const;
  double getHh(int idir, int ivar, int jvar, int ipas) const;
  double getSw(int idir, int ivar, int jvar, int ipas) const;
  double getUtilize(int idir, int ivar, int jvar, int ipas) const;

  VectorDouble getGgVec(int idir,
                        int ivar,
                        int jvar,
                        bool asCov = false,
                        bool flagNormalized = false) const;
  VectorDouble getHhVec(int idir, int ivar, int jvar) const;
  VectorDouble getSwVec(int idir, int ivar, int jvar) const;
  VectorDouble getUtilizeVec(int idir, int ivar, int jvar) const;

  const VectorDouble& getAllGg(int idir) const;
  const VectorDouble& getAllHh(int idir) const;
  const VectorDouble& getAllSw(int idir) const;
  const VectorDouble& getAllUtilize(int idir) const;

  void setGgByIndex(int idir, int i, double gg);
  void setHhByIndex(int idir, int i, double hh);
  void setSwByIndex(int idir, int i, double sw);
  void setUtilizeByIndex(int idir, int i, double utilize);

  void setSw(int idir, int ivar, int jvar, int ipas, double sw);
  void setHh(int idir, int ivar, int jvar, int ipas, double hh);
  void setGg(int idir, int ivar, int jvar, int ipas, double gg);
  void setUtilize(int idir, int ivar, int jvar, int ipas, double utilize);

  void updateSwByIndex(int idir, int i, double sw);
  void updateHhByIndex(int idir, int i, double hh);
  void updateGgByIndex(int idir, int i, double gg);

  int getCenter(int ivar = 0, int jvar = 0, int idir = 0) const;

  int  internalVariableResize();
  void internalDirectionResize(int ndir = 0, bool flagDirs = true);

  double getHmax(int ivar=-1, int jvar=-1, int idir=-1) const;
  VectorDouble getHRange(int ivar=-1, int jvar=-1, int idir=-1) const;
  double getGmax(int ivar = -1,
                 int jvar = -1,
                 int idir = -1,
                 bool flagAbs = false,
                 bool flagSill = false) const;
  VectorDouble getGRange(int ivar = -1,
                         int jvar = -1,
                         int idir = -1,
                         bool flagSill = false) const;

  void patchCenter(int idir, int nech, double rho);

  int fill(int idir,
           const VectorDouble& sw,
           const VectorDouble& gg,
           const VectorDouble& hh);

  int getDirAddress(int idir,
                     int ivar,
                     int jvar,
                     int ipas,
                     bool flag_abs = false,
                     int sens = 0) const;
  int getVarAddress(int ivar, int jvar) const;
  int getLagTotalNumber(int idir) const;

  int attachDb(Db* db,
               const VectorDouble& vars = VectorDouble(),
               const VectorDouble& means = VectorDouble());
  int computeByKey(const String& calcul_name = "vg",
                   bool flag_grid = false,
                   bool flag_gen = false,
                   bool flag_sample = false,
                   bool verr_mode = false,
                   Model *model = nullptr,
                   bool verbose = false);
  int computeIndicByKey(const String& calcul_name = "vg",
                        bool flag_grid = false,
                        bool flag_gen = false,
                        bool flag_sample = false,
                        bool verr_mode = false,
                        Model *model = nullptr,
                        bool verbose = false,
                        int nfacmax = -1);
  int compute(const ECalcVario& calcul = ECalcVario::VARIOGRAM,
              bool flag_grid = false,
              bool flag_gen = false,
              bool flag_sample = false,
              bool verr_mode = false,
              Model *model = nullptr,
              bool verbose = false);
  int computeIndic(const ECalcVario& calcul = ECalcVario::VARIOGRAM,
                   bool flag_grid = false,
                   bool flag_gen = false,
                   bool flag_sample = false,
                   bool verr_mode = false,
                   Model *model = nullptr,
                   bool verbose = false,
                   int nfacmax = -1);

  // Pipe to the DirParam
  const DirParam& getDirParam(int idir) const { return _varioparam.getDirParam(idir); }
  int getDirectionNumber() const { return _varioparam.getDirectionNumber(); }
  const VectorDouble& getDates() const { return _varioparam.getDates(); }
  bool hasDate() const { return _varioparam.hasDate(); }
  double getDates(int idate, int icas) const { return _varioparam.getDate(idate, icas); }
  int getDateNumber() const { return _varioparam.getDateNumber(); }
  double getScale() const { return _varioparam.getScale(); }
  int getDimensionNumber() const { return getDirParam(0).getDimensionNumber(); }

  void setScale(double scale) { _varioparam.setScale(scale); }
  void addDirs(const DirParam& dirparam) { _varioparam.addDirs(dirparam); }

  int getLagNumber(int idir) const { return getDirParam(idir).getLagNumber(); }
  double getDPas(int idir) const { return getDirParam(idir).getDPas(); }
  int getDimensionNumber(int idir) const { return getDirParam(idir).getDimensionNumber(); }
  const VectorDouble& getCodir(int idir) const { return getDirParam(idir).getCodir(); }
  double getCodir(int idir, int idim) const { return getDirParam(idir).getCodir(idim); }
  double getMaximumDistance(int idir) const { return getDirParam(idir).getMaximumDistance(); }
  int getIdate(int idir) const { return getDirParam(idir).getIdate(); }
  double getGrincr(int idir, int idim) { return getDirParam(idir).getGrincr(idim); }

  void setNVar(int nvar) { _nVar = nvar; }
  void setCalculName(const String calcul_name);

  const VarioParam& getVarioParam() const { return _varioparam; }

protected:
  virtual int _deserialize(FILE* file, bool verbose = false);
  virtual int _serialize(FILE* file, bool verbose = false) const;

  virtual int _deserialize2(std::istream& is, bool verbose = false) override;
  virtual int _serialize2(std::ostream& os, bool verbose = false) const override;

private:
  bool _isVariableValid(int ivar) const;
  bool _isDirectionValid(int idir) const;
  bool _isBivariableValid(int ijvar) const;
  bool _isAddressValid(int i, int idir) const;
  void _initMeans();
  void _initVars();
  int  _getNVar(const Db* db);
  VectorInt _getVariableInterval(int ivar) const;
  VectorInt _getDirectionInterval(int idir) const;
  String _toStringByDirection(const AStringFormat* strfmt, int idir) const;
  void _directionResize(int idir);
  void _setDPasFromGrid(bool flag_grid);
  void _setFlagAsym();
  VectorDouble _varsFromProportions(VectorDouble props);

private:
  int                _nVar;
  VarioParam         _varioparam;
  VectorDouble       _means;
  VectorDouble       _vars;
  ECalcVario         _calcul;
  bool               _flagSample;
  Db*                _db;
  VectorVectorDouble _sw;      /* Array for number of lags */
  VectorVectorDouble _gg;      /* Array for average variogram values */
  VectorVectorDouble _hh;      /* Array for average distance values */
  VectorVectorDouble _utilize; /* Array to mention if a lag is used or not */
  mutable bool       _flagAsym;
};
