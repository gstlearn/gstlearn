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

#include "gstlearn_export.hpp"
#include "geoslib_define.h"
#include "geoslib_d.h"

#include "Variogram/AVario.hpp"
#include "Variogram/VarioParam.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Basic/ASerializable.hpp"
#include "Geometry/ABiTargetCheck.hpp"

class Db;
class Model;
class DirParam;
class AAnam;

/**
 * \brief
 * Class containing the Spatial Characteristics as calculated experimentally from the data (contained in a Db).
 *
 * The experimental Spatial Characteristics is usually referred to as the experimental **variogram**.
 * However, note that it can rather calculate other results such as a Covariance or a Madogram. All these
 * quantities can be regrouped by considering them as **two-points** statistics.
 * For a complete list of calculation methods, please refer to ECalcVario.hpp.
 *
 * This class is composed of two parts:
 * - the first part describes the rule when comparing two samples from the Db. They are defined by:
 *
 *    - the definition of the **Geometry**: e.g. definition of calculation direction, tolerances.
 * For more information, please refer to VarioParam.hpp
 *    - the definition of the calculations **Options**: e.g. calculation method.
 *    - some additional **Conditions** used during calculations: e.g. usage of *Faults*.
 * For more information, please refer to ABiTargetCheck.hpp.
 *
 * - the second part are the results of the calculations
 *
 * **Results**
 *
 * All the Spatial Characteristics are calculated:
 * - from the sample values of active samples contained in a Db,
 * - for all the variables (defined with the locator ELoc.Z): in the multivariate case, simple and
 * cross-variograms are calculated
 * - for a series of distance lags.
 *
 * They are always expressed as a table with one row per distance lag and three columns containing:
 * - the number of pairs
 * - the average value of the distance
 * - the average value of the two-points statistics
 *
 * Note that:
 * - the lags for which no pair is found are skipped.
 * - some methods correspond to an **even** function (values are equal whether the distance between
 * the two end-points is counted positively or negatively: then only one-sided results are stored.
 * For **odd**, the results of both sides are stored.
 * - for a number of lags equal to *N*, the number of rows is {N+1} when the function is even and
 * {2N+1} when the function is odd.
 * - in the multivariate case (NV variables), the number of rows is multiplied by NV*(NV+1)/2.
 * In order to avoid any indexing problem, the user should use the assessors provided in order to access to the information
 * relative to the target pair of variables.
 .
 *
 */
class GSTLEARN_EXPORT Vario : public AVario, public ASerializable
{
public:
  Vario(const VarioParam& varioparam);
  Vario(const Vario& r);
  Vario& operator=(const Vario& r);
  virtual ~Vario();

public:
  /// ICloneable interface
  IMPLEMENT_CLONING(Vario)

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// AVario Interface
  double _getIVAR(const Db *db, int iech, int ivar) const override;
  void _setResult(int iech1,
                  int iech2,
                  int nvar,
                  int ipas,
                  int ivar,
                  int jvar,
                  int orient,
                  double ww,
                  double dist,
                  double value) override;

  static Vario* create(const VarioParam& varioparam);
  static Vario* createFromNF(const String& neutralFilename, bool verbose = true);
  static Vario* createRegularizeFromModel(const Model& model,
                                          const VarioParam& varioparam,
                                          const VectorDouble& ext,
                                          const VectorInt& ndisc,
                                          const VectorDouble& angles,
                                          bool asCov = false);
  static Vario* createTransformZToY(const Vario& varioZ,
                                    const AAnam* anam);
  static Vario* createTransformYToZ(const Vario& varioY,
                                    const AAnam* anam);
  static Vario* createReduce(const Vario& varioIn,
                             const VectorInt &varcols,
                             const VectorInt &dircols,
                             bool asSymmetric = false);
  static Vario* computeFromDb(const VarioParam& varioparam,
                              Db* db,
                              const ECalcVario& calcul = ECalcVario::fromKey("VARIOGRAM"),
                              bool flag_sample = false,
                              bool verr_mode = false,
                              Model *model = nullptr,
                              int niter_UK = 0,
                              bool verbose = false);

  void resetReduce(const VectorInt &varcols,
                   const VectorInt &dircols,
                   bool asSymmetric = false);

  bool              getFlagAsym() const { return _flagAsym; }
  bool              drawOnlyPositiveX(int ivar, int jvar) const;
  bool              drawOnlyPositiveY(int ivar, int jvar) const;

  int    getVariableNumber() const { return _nVar; }
  const  VectorDouble& getMeans() const { return _means; }
  double getMean(int ivar) const;

  double getVar(int ivar, int jvar) const;
  double getVarIndex(int ijvar) const;
  const  VectorDouble& getVars()  const { return _vars; }
  void setMeans(const VectorDouble& means);
  void setMean(double mean, int ivar=0);
  void setVar(double value, int ivar=0, int jvar=0);
  void setVars(const VectorDouble& vars);
  void setVarIndex(int ijvar, double value);
  void setDb(Db* db);

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
                        bool flagNormalized = false,
                        bool compress = true) const;
  VectorDouble getHhVec(int idir = 0,
                        int ivar = 0,
                        int jvar = 0,
                        bool compress = true) const;
  VectorDouble getSwVec(int idir = 0,
                        int ivar = 0,
                        int jvar = 0,
                        bool compress = true) const;
  VectorDouble getUtilizeVec(int idir = 0,
                             int ivar = 0,
                             int jvar = 0,
                             bool compress = true) const;

  void setSwVec(int idir, int ivar, int jvar, const VectorDouble& sw);
  void setHhVec(int idir, int ivar, int jvar, const VectorDouble& hh);
  void setGgVec(int idir, int ivar, int jvar, const VectorDouble& gg);

  VectorDouble getGgs(int idir = 0,
                      int ivar = 0,
                      int jvar = 0,
                      const VectorInt &ipas = VectorInt()) const;
  VectorDouble setGgs(int idir, int ivar, int jvar, const VectorInt& ipas, const VectorDouble& values);

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
  int getNext(int ivar = 0, int jvar = 0, int idir = 0, int shift = 1) const;

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

  int compute(Db* db,
              const ECalcVario& calcul = ECalcVario::fromKey("VARIOGRAM"),
              bool flag_sample = false,
              bool verr_mode = false,
              const Model* model = nullptr,
              int niter_UK = 0,
              bool verbose = false);
  int computeIndic(Db* db,
                   const ECalcVario& calcul = ECalcVario::fromKey("VARIOGRAM"),
                   bool flag_sample = false,
                   bool verr_mode = false,
                   const Model* model = nullptr,
                   int niter_UK = 0,
                   bool verbose = false,
                   int nfacmax = -1);
  int computeGeometry(Db *db, Vario_Order *vorder, int *npair);
  int computeVarioVect(Db *db, int ncomp);
  int computeGeometryMLayers(Db *db, VectorInt& seltab, Vario_Order *vorder) const;

  int regularizeFromModel(const Model &model,
                          const VectorDouble &ext,
                          const VectorInt &ndisc,
                          const VectorDouble &angles = VectorDouble(),
                          const CovCalcMode *mode = nullptr,
                          bool asCov = false);
  int regularizeFromDbGrid(Model* model,
                           const Db& db,
                           const CovCalcMode *mode = nullptr);
  void getExtension(int ivar,
                    int jvar,
                    int idir0,
                    int flag_norm,
                    int flag_vars,
                    double distmin,
                    double distmax,
                    double varmin,
                    double varmax,
                    int *flag_hneg,
                    int *flag_gneg,
                    double *c0,
                    double *hmin,
                    double *hmax,
                    double *gmin,
                    double *gmax);
  int sampleModel(Model *model, const CovCalcMode*  mode = nullptr);

  // Pipe to the DirParam
  const DirParam& getDirParam(int idir) const { return _varioparam.getDirParam(idir); }
  int getDirectionNumber() const { return _varioparam.getDirectionNumber(); }
  const VectorDouble& getDates() const { return _varioparam.getDates(); }
  bool hasDate() const { return _varioparam.hasDate(); }
  double getDates(int idate, int icas) const { return _varioparam.getDate(idate, icas); }
  int getDateNumber() const { return _varioparam.getDateNumber(); }
  double getScale() const { return _varioparam.getScale(); }
  int getDimensionNumber() const { return getDirParam(0).getNDim(); }
  const ASpace* getSpace() const { return getDirParam(0).getSpace(); }

  void setScale(double scale) { _varioparam.setScale(scale); }
  void addDirs(const DirParam& dirparam) { _varioparam.addDir(dirparam); }

  int getLagNumber(int idir) const { return getDirParam(idir).getLagNumber(); }
  double getDPas(int idir) const { return getDirParam(idir).getDPas(); }
  int getDimensionNumber(int idir) const { return getDirParam(idir).getNDim(); }
  VectorDouble getCodirs(int idir) const;
  double getCodir(int idir, int idim) const;
  double getMaximumDistance(int idir) const { return getDirParam(idir).getMaximumDistance(); }
  int getIdate(int idir) const { return getDirParam(idir).getIdate(); }
  VectorInt getGrincrs(int idir) const { return getDirParam(idir).getGrincrs(); }
  double getGrincr(int idir, int idim) const { return getDirParam(idir).getGrincr(idim); }
  bool isDefinedForGrid() const { return _varioparam.isDefinedForGrid(); }
  void setNVar(int nvar) { _nVar = nvar; }
  void setCalculName(const String& calcul_name);
  void setVariableNames(const VectorString &variableNames) { _variableNames = variableNames; }
  void setVariableName(int ivar, const String &variableName);

  int  prepare(const ECalcVario &calcul = ECalcVario::fromKey("VARIOGRAM"), bool defineList = true);

  const VarioParam& getVarioParam() const { return _varioparam; }
  int getBiPtsNumberPerDirection() const { return _biPtsPerDirection; }
  const ABiTargetCheck* getBipts(int idir, int rank) const { return _bipts[_getBiPtsRank(idir, rank)]; }
  bool keepPair(int idir, SpaceTarget &T1, SpaceTarget &T2, double *dist) const;
  int getRankFromDirAndDate(int idir, int idate) const;
  const VectorString& getVariableNames() const { return _variableNames; }
  String getVariableName(int ivar) const;

  int transformCut(int nh, double ycut);
  int transformZToY(const AAnam *anam);
  int transformYToZ(const AAnam *anam);

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
  static VectorDouble _varsFromProportions(VectorDouble props);
  void _clearBiTargetCheck();
  void _addBiTargetCheck(ABiTargetCheck* abpc);
  void _setListBiTargetCheck();
  int  _getBiPtsNumber() const { return (int) _bipts.size(); }
  int  _getBiPtsRank(int idir, int rank) const;

  int _compute(Db *db,
               int flag_sample,
               int verr_mode,
               const Model* model,
               int niter_UK,
               bool verbose);
  int _calculateGeneral(Db *db,
                        int flag_sample,
                        int verr_mode);
  int  _calculateGenOnLine(Db *db, int norder);
  int  _calculateGenOnGrid(DbGrid *db, int norder);
  int  _calculateOnGrid(DbGrid *db);

  static int  _getRelativeSampleRank(Db *db, int iech0);
  int  _updateUK(Db *db, Vario_Order *vorder);
  void _patchC00(Db *db, int idir);
  int  _get_generalized_variogram_order();
  void _getStatistics(Db *db);
  int  _updateVerr(Db *db, int idir, Vario_Order *vorder, int verr_mode);
  static double _s(Db *db, int iech, int jech);
  double _g(Db *db, int iech, int jech) const;
  void _calculateBiasLocal(Db *db,
                           int idir,
                           int ipas,
                           Vario_Order *vorder,
                           int ifirst,
                           int ilast);
  void _calculateBiasGlobal(Db *db);
  double _getBias(int iiech, int jjech);

  void _calculateFromGeometry(Db *db, int idir, Vario_Order *vorder);
  int  _calculateGeneralSolution1(Db *db, int idir, const int *rindex, Vario_Order *vorder);
  int  _calculateGeneralSolution2(Db *db, int idir, const int *rindex);
  int  _calculateOnGridSolution(DbGrid *db, int idir);
  int  _calculateGenOnGridSolution(DbGrid *db, int idir, int norder);
  int  _calculateVarioVectSolution(Db *db, int idir, int ncomp, const int *rindex);
  void _calculateOnLineSolution(Db *db, int idir, int norder);

  void _driftManage(Db *db);
  int  _driftEstimateCoefficients(Db *db);

  static void _printDebug(int iech1,
                          int iech2,
                          int ivar,
                          int jvar,
                          int ilag,
                          double scale,
                          double value);
  void _centerCovariance(Db *db, int idir);
  void _getVarioVectStatistics(Db *db, int ncomp);
  void _rescale(int idir);
  bool _isCompatible(const Db *db) const;
  static double _linear_interpolate(int n,
                                    const VectorDouble& x,
                                    const VectorDouble& y,
                                    double x0);
  MatrixSquareGeneral _evalAverageDbIncr(Model *model,
                                         const Db &db,
                                         const VectorDouble &incr = VectorDouble(),
                                         const CovCalcMode *mode = nullptr) const;

private:
  int                _nVar;
  VarioParam         _varioparam;
  VectorDouble       _means;
  VectorDouble       _vars;
  bool               _flagSample;
  Db*                _db;
  VectorVectorDouble _sw;      /* Array for number of lags */
  VectorVectorDouble _gg;      /* Array for average variogram values */
  VectorVectorDouble _hh;      /* Array for average distance values */
  VectorVectorDouble _utilize; /* Array to mention if a lag is used or not */

  int _biPtsPerDirection;
  std::vector<ABiTargetCheck*> _bipts;
  mutable bool       _flagAsym;

  bool _verbose;
  bool _flag_UK;
  int  _niter_UK;

  VectorString                _variableNames;
  mutable Model*              _model;  // Model pointer (not to be deleted) for drift removal
  mutable VectorDouble        _BETA;
  mutable VectorDouble        _DRFDIAG;
  mutable MatrixRectangular   _DRFXA;
  mutable MatrixRectangular   _DRFGX;
  mutable MatrixRectangular   _DRFTAB;
  mutable MatrixSquareGeneral _DRFXGX;
};
