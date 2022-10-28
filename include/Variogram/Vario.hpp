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

#include "Enum/ECalcVario.hpp"

#include "Variogram/VarioParam.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"

class Db;
class Model;
class DirParam;
class AAnam;

/**
 * Experimental Variogram (not only): TODO : to be improved
 */
class GSTLEARN_EXPORT Vario : public AStringable, public ASerializable, public ICloneable
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
  /// ICloneable interface
  IMPLEMENT_CLONING(Vario)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static Vario* create(const VarioParam* varioparam,
                       Db* db = nullptr,
                       const VectorDouble& means = VectorDouble(),
                       const VectorDouble& vars = VectorDouble());
  static Vario* createFromNF(const String& neutralFilename, bool verbose = true);
  static Vario* computeFromDb(const VarioParam* varioparam,
                              Db* db,
                              const ECalcVario& calcul = ECalcVario::VARIOGRAM,
                              bool flag_gen = false,
                              bool flag_sample = false,
                              bool verr_mode = false,
                              Model *model = nullptr,
                              bool verbose = false);
  static Vario* createRegularizeFromModel(const Model* model,
                                          const VarioParam* varioparam,
                                          const VectorDouble& ext,
                                          const VectorInt& ndisc,
                                          const VectorDouble& angles);
  static Vario* createTransformZToY(const Vario* varioZ,
                                    const AAnam* anam,
                                    double cvv);
  static Vario* createTransformYToZ(const Vario* varioY,
                                    const AAnam* anam,
                                    const Model* model);
  static Vario* createReduce(const Vario *varioIn,
                             const VectorInt &varcols,
                             const VectorInt &dircols,
                             bool asSymmetric = false);
  void reduce(const VectorInt& varcols,
              const VectorInt& dircols,
              bool asSymmetric = false);

  const ECalcVario& getCalcul() const { return _calcul; }
  ECalcVario    getCalculType(const String& calcul_name) const;
  bool          getFlagAsym() const { return _flagAsym; }
  bool          drawOnlyPositiveX(int ivar, int jvar) const;
  bool          drawOnlyPositiveY(int ivar, int jvar) const;

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

  double getGgByIndex(int idir = 0, int i = 0) const;
  double getHhByIndex(int idir = 0, int i = 0) const;
  double getSwByIndex(int idir = 0, int i = 0) const;
  double getUtilizeByIndex(int idir = 0, int i = 0) const;

  double getGg(int idir = 0,
               int ivar = 0,
               int jvar = 0,
               int ipas = 0,
               bool asCov = false,
               bool flagNormalized = false) const;
  double getHh(int idir = 0, int ivar = 0, int jvar = 0, int ipas = 0) const;
  double getSw(int idir = 0, int ivar = 0, int jvar = 0, int ipas = 0) const;
  double getUtilize(int idir = 0, int ivar = 0, int jvar = 0, int ipas = 0) const;

  VectorVectorDouble getVec(int idir = 0, int ivar = 0, int jvar = 0) const;
  VectorDouble getGgVec(int idir = 0,
                        int ivar = 0,
                        int jvar = 0,
                        bool asCov = false,
                        bool flagNormalized = false) const;
  VectorDouble getHhVec(int idir = 0, int ivar = 0, int jvar = 0) const;
  VectorDouble getSwVec(int idir = 0, int ivar = 0, int jvar = 0) const;
  VectorDouble getUtilizeVec(int idir = 0, int ivar = 0, int jvar = 0) const;

  const VectorDouble& getAllGg(int idir = 0) const;
  const VectorDouble& getAllHh(int idir = 0) const;
  const VectorDouble& getAllSw(int idir = 0) const;
  const VectorDouble& getAllUtilize(int idir = 0) const;

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
                   bool flag_gen = false,
                   bool flag_sample = false,
                   bool verr_mode = false,
                   Model *model = nullptr,
                   bool verbose = false);
  int computeIndicByKey(const String& calcul_name = "vg",
                        bool flag_gen = false,
                        bool flag_sample = false,
                        bool verr_mode = false,
                        Model *model = nullptr,
                        bool verbose = false,
                        int nfacmax = -1);
  int compute(const ECalcVario& calcul = ECalcVario::VARIOGRAM,
              bool flag_gen = false,
              bool flag_sample = false,
              bool verr_mode = false,
              Model *model = nullptr,
              bool verbose = false);
  int computeIndic(const ECalcVario& calcul = ECalcVario::VARIOGRAM,
                   bool flag_gen = false,
                   bool flag_sample = false,
                   bool verr_mode = false,
                   Model *model = nullptr,
                   bool verbose = false,
                   int nfacmax = -1);
  int transformZToY(const AAnam *anam, double cvv);
  int transformYToZ(const AAnam *anam, const Model *model);
  int modelRegularize(const Model* model,
                      const VectorDouble& ext,
                      const VectorInt& ndisc,
                      const VectorDouble& angles = VectorDouble(),
                      const CovCalcMode& mode = CovCalcMode());

  // Pipe to the DirParam
  const DirParam& getDirParam(int idir) const { return _varioparam.getDirParam(idir); }
  int getDirectionNumber() const { return _varioparam.getDirectionNumber(); }
  const VectorDouble& getDates() const { return _varioparam.getDates(); }
  bool hasDate() const { return _varioparam.hasDate(); }
  double getDates(int idate, int icas) const { return _varioparam.getDate(idate, icas); }
  int getDateNumber() const { return _varioparam.getDateNumber(); }
  double getScale() const { return _varioparam.getScale(); }
  int getDimensionNumber() const { return getDirParam(0).getNDim(); }

  void setScale(double scale) { _varioparam.setScale(scale); }
  void addDirs(const DirParam& dirparam) { _varioparam.addDir(dirparam); }

  int getLagNumber(int idir) const { return getDirParam(idir).getLagNumber(); }
  double getDPas(int idir) const { return getDirParam(idir).getDPas(); }
  int getDimensionNumber(int idir) const { return getDirParam(idir).getNDim(); }
  const VectorDouble& getCodir(int idir) const { return getDirParam(idir).getCodir(); }
  double getCodir(int idir, int idim) const { return getDirParam(idir).getCodir(idim); }
  double getMaximumDistance(int idir) const { return getDirParam(idir).getMaximumDistance(); }
  int getIdate(int idir) const { return getDirParam(idir).getIdate(); }
  double getGrincr(int idir, int idim) { return getDirParam(idir).getGrincr(idim); }
  bool isDefinedForGrid() const { return _varioparam.isDefinedForGrid(); }

  void setNVar(int nvar) { _nVar = nvar; }
  void setCalculName(const String calcul_name);

  const VarioParam& getVarioParam() const { return _varioparam; }

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "Vario"; }

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
