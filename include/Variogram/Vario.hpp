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

#include "Anamorphosis/AAnam.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Basic/ASerializable.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Db/DbGrid.hpp"
#include "Geometry/ABiTargetCheck.hpp"
#include "Variogram/AVario.hpp"
#include "Variogram/VarioParam.hpp"

namespace gstlrn
{
typedef struct
{
  Id nalloc;
  Id npair;
  Id size_aux;
  Id flag_dist;
  VectorInt tab_iech;
  VectorInt tab_jech;
  VectorInt tab_ipas;
  VectorInt tab_sort;
  std::vector<char> tab_aux_iech;
  std::vector<char> tab_aux_jech;
  VectorDouble tab_dist;
} Vario_Order;

class Db;
class Model;
class DirParam;

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
 *
 */

class GSTLEARN_EXPORT Vario: public AVario, public ASerializable
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
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// AVario Interface
  double _getIVAR(const Db* db, Id iech, Id ivar) const override;
  void _setResult(Id iech1,
                  Id iech2,
                  Id nvar,
                  Id ilag,
                  Id ivar,
                  Id jvar,
                  Id orient,
                  double ww,
                  double dist,
                  double value) override;

  static Vario* create(const VarioParam& varioparam);
  static Vario* createFromNF(const String& NFFilename, bool verbose = true);
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
                             const VectorInt& varcols,
                             const VectorInt& dircols,
                             bool asSymmetric = false);
  static Vario* computeFromDb(const VarioParam& varioparam,
                              Db* db,
                              const ECalcVario& calcul = ECalcVario::fromKey("VARIOGRAM"),
                              bool flag_sample         = false,
                              bool verr_mode           = false,
                              Model* model             = nullptr,
                              Id niter_UK              = 0,
                              bool verbose             = false);

  void resetReduce(const VectorInt& varcols,
                   const VectorInt& dircols,
                   bool asSymmetric = false);

  bool getFlagAsym() const { return _flagAsym; }
  bool drawOnlyPositiveX(Id ivar, Id jvar) const;
  bool drawOnlyPositiveY(Id ivar, Id jvar) const;

  Id getNVar() const { return _nVar; }
  const VectorDouble& getMeans() const { return _means; }
  double getMean(Id ivar) const;

  double getVar(Id ivar, Id jvar) const;
  MatrixSymmetric getVarMatrix() const;
  double getVarIndex(Id ijvar) const;
  const VectorDouble& getVars() const { return _vars; }
  void setMeans(const VectorDouble& means);
  void setMean(double mean, Id ivar = 0);
  void setVar(double value, Id ivar = 0, Id jvar = 0);
  void setVars(const VectorDouble& vars);
  void setVarIndex(Id ijvar, double value);
  void setDb(Db* db);

  Id getDirSize(Id idir) const;

  double getGgByIndex(Id idir = 0, Id i = 0) const;
  double getHhByIndex(Id idir = 0, Id i = 0) const;
  double getSwByIndex(Id idir = 0, Id i = 0) const;
  double getUtilizeByIndex(Id idir = 0, Id i = 0) const;

  double getGg(Id idir             = 0,
               Id ivar             = 0,
               Id jvar             = 0,
               Id ilag             = 0,
               bool asCov          = false,
               bool flagNormalized = false) const;
  double getHh(Id idir = 0, Id ivar = 0, Id jvar = 0, Id ilag = 0) const;
  double getSw(Id idir = 0, Id ivar = 0, Id jvar = 0, Id ilag = 0) const;
  double getUtilize(Id idir = 0, Id ivar = 0, Id jvar = 0, Id ilag = 0) const;

  VectorVectorDouble getVec(Id idir = 0, Id ivar = 0, Id jvar = 0) const;
  VectorDouble getGgVec(Id idir             = 0,
                        Id ivar             = 0,
                        Id jvar             = 0,
                        bool asCov          = false,
                        bool flagNormalized = false,
                        bool compress       = true) const;
  VectorDouble getHhVec(Id idir       = 0,
                        Id ivar       = 0,
                        Id jvar       = 0,
                        bool compress = true) const;
  VectorDouble getSwVec(Id idir       = 0,
                        Id ivar       = 0,
                        Id jvar       = 0,
                        bool compress = true) const;
  VectorDouble getUtilizeVec(Id idir       = 0,
                             Id ivar       = 0,
                             Id jvar       = 0,
                             bool compress = true) const;

  void setSwVec(Id idir, Id ivar, Id jvar, const VectorDouble& sw);
  void setHhVec(Id idir, Id ivar, Id jvar, const VectorDouble& hh);
  void setGgVec(Id idir, Id ivar, Id jvar, const VectorDouble& gg);

  VectorDouble getGgs(Id idir               = 0,
                      Id ivar               = 0,
                      Id jvar               = 0,
                      const VectorInt& ilag = VectorInt()) const;
  VectorDouble setGgs(Id idir, Id ivar, Id jvar, const VectorInt& ilag, const VectorDouble& values);

  const VectorDouble& getAllGg(Id idir = 0) const;
  const VectorDouble& getAllHh(Id idir = 0) const;
  const VectorDouble& getAllSw(Id idir = 0) const;
  const VectorDouble& getAllUtilize(Id idir = 0) const;

  void setGgByIndex(Id idir, Id i, double gg, bool flagCheck = true);
  void setHhByIndex(Id idir, Id i, double hh, bool flagCheck = true);
  void setSwByIndex(Id idir, Id i, double sw, bool flagCheck = true);
  void setUtilizeByIndex(Id idir, Id i, double utilize, bool flagCheck = true);

  void setSw(Id idir, Id ivar, Id jvar, Id ilag, double sw, bool flagCheck = true);
  void setHh(Id idir, Id ivar, Id jvar, Id ilag, double hh, bool flagCheck = true);
  void setGg(Id idir, Id ivar, Id jvar, Id ilag, double gg, bool flagCheck = true);
  void setUtilize(Id idir,
                  Id ivar,
                  Id jvar,
                  Id ilag,
                  double utilize,
                  bool flagCheck = true);

  void updateSwByIndex(Id idir, Id i, double sw, bool flagCheck = true);
  void updateHhByIndex(Id idir, Id i, double hh, bool flagCheck = true);
  void updateGgByIndex(Id idir, Id i, double gg, bool flagCheck = true);

  Id getCenter(Id ivar = 0, Id jvar = 0, Id idir = 0) const;
  Id getNext(Id ivar = 0, Id jvar = 0, Id idir = 0, Id shift = 1) const;

  Id internalVariableResize();
  void internalDirectionResize(Id ndir = 0, bool flagDirs = true);

  double getHmax(Id ivar = -1, Id jvar = -1, Id idir = -1) const;
  VectorDouble getHRange(Id ivar = -1, Id jvar = -1, Id idir = -1) const;
  double getGmax(Id ivar       = -1,
                 Id jvar       = -1,
                 Id idir       = -1,
                 bool flagAbs  = false,
                 bool flagSill = false) const;
  VectorDouble getGRange(Id ivar       = -1,
                         Id jvar       = -1,
                         Id idir       = -1,
                         bool flagSill = false) const;

  void patchCenter(Id idir, Id nech, double rho);

  Id fill(Id idir,
          const VectorDouble& sw,
          const VectorDouble& gg,
          const VectorDouble& hh);

  Id getDirAddress(Id idir,
                   Id ivar,
                   Id jvar,
                   Id ilag,
                   bool flag_abs  = false,
                   Id sens        = 0,
                   bool flagCheck = true) const;
  Id getVarAddress(Id ivar, Id jvar) const;
  Id getNLagTotal(Id idir) const;

  Id compute(Db* db,
             const ECalcVario& calcul = ECalcVario::fromKey("VARIOGRAM"),
             bool flag_sample         = false,
             bool verr_mode           = false,
             const Model* model       = nullptr,
             Id niter_UK              = 0,
             bool verbose             = false);
  Id computeIndic(Db* db,
                  const ECalcVario& calcul = ECalcVario::fromKey("VARIOGRAM"),
                  bool flag_sample         = false,
                  bool verr_mode           = false,
                  const Model* model       = nullptr,
                  Id niter_UK              = 0,
                  bool verbose             = false,
                  Id nfacmax               = -1);
  Id computeGeometry(Db* db, Vario_Order* vorder, Id* npair);
  Id computeVarioVect(Db* db, Id ncomp);
  Id computeGeometryMLayers(Db* db, VectorInt& seltab, Vario_Order* vorder) const;

  Id regularizeFromModel(const Model& model,
                         const VectorDouble& ext,
                         const VectorInt& ndisc,
                         const VectorDouble& angles = VectorDouble(),
                         const CovCalcMode* mode    = nullptr,
                         bool asCov                 = false);
  Id regularizeFromDbGrid(Model* model,
                          const Db& db,
                          const CovCalcMode* mode = nullptr);
  void getExtension(Id ivar,
                    Id jvar,
                    Id idir0,
                    Id flag_norm,
                    Id flag_vars,
                    double distmin,
                    double distmax,
                    double varmin,
                    double varmax,
                    Id* flag_hneg,
                    Id* flag_gneg,
                    double* c0,
                    double* hmin,
                    double* hmax,
                    double* gmin,
                    double* gmax);
  Id sampleModel(Model* model, const CovCalcMode* mode = nullptr);

  // Pipe to the DirParam
  const DirParam& getDirParam(Id idir) const { return _varioparam.getDirParam(idir); }
  Id getNDir() const { return _varioparam.getNDir(); }
  const VectorDouble& getDates() const { return _varioparam.getDates(); }
  bool hasDate() const { return _varioparam.hasDate(); }
  double getDates(Id idate, Id icas) const { return _varioparam.getDate(idate, icas); }
  Id getNDate() const { return _varioparam.getNDate(); }
  double getScale() const { return _varioparam.getScale(); }
  Id getNDim() const { return static_cast<Id>(getDirParam(0).getNDim()); }
  ASpaceSharedPtr getSpace() const { return getDirParam(0).getSpace(); }

  void setScale(double scale) { _varioparam.setScale(scale); }
  void addDirs(const DirParam& dirparam) { _varioparam.addDir(dirparam); }

  Id getNLag(Id idir) const { return getDirParam(idir).getNLag(); }
  double getDPas(Id idir) const { return getDirParam(idir).getDPas(); }
  Id getNDim(Id idir) const { return static_cast<Id>(getDirParam(idir).getNDim()); }
  VectorDouble getCodirs(Id idir) const;
  double getCodir(Id idir, Id idim) const;
  double getMaximumDistance(Id idir) const { return getDirParam(idir).getMaximumDistance(); }
  double getMaximumDistance() const;
  Id getIdate(Id idir) const { return getDirParam(idir).getIdate(); }
  VectorInt getGrincrs(Id idir) const { return getDirParam(idir).getGrincrs(); }
  double getGrincr(Id idir, Id idim) const { return getDirParam(idir).getGrincr(idim); }
  bool isDefinedForGrid() const { return _varioparam.isDefinedForGrid(); }
  void setNVar(Id nvar) { _nVar = nvar; }
  void setCalculByName(const String& calcul_name);
  void setVariableNames(const VectorString& variableNames) { _variableNames = variableNames; }
  void setVariableName(Id ivar, const String& variableName);

  Id prepare(const ECalcVario& calcul = ECalcVario::fromKey("VARIOGRAM"), bool defineList = true);

  const VarioParam& getVarioParam() const { return _varioparam; }
  Id getNBiPtsPerDir() const { return _biPtsPerDirection; }
  const ABiTargetCheck* getBipts(Id idir, Id rank) const { return _bipts[_getBiPtsRank(idir, rank)]; }
  bool keepPair(Id idir, SpaceTarget& T1, SpaceTarget& T2, double* dist) const;
  Id getRankFromDirAndDate(Id idir, Id idate) const;
  const VectorString& getVariableNames() const { return _variableNames; }
  String getVariableName(Id ivar) const;

  Id transformCut(Id nh, double ycut);
  Id transformZToY(const AAnam* anam);
  Id transformYToZ(const AAnam* anam);

  bool isLagCorrect(Id idir, Id k) const;
  double getC00(Id idir, Id ivar, Id jvar) const;
  VectorDouble computeWeightPerDirection() const;
  Id getTotalLagsPerDirection() const;
  VectorDouble computeWeightsFromVario(Id wmode) const;

protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "Vario"; }

private:
  bool _isVariableValid(Id ivar, bool flagCheck = true) const;
  bool _isDirectionValid(Id idir, bool flagCheck = true) const;
  bool _isBivariableValid(Id ijvar, bool flagCheck = true) const;
  bool _isAddressValid(Id i, Id idir, bool flagCheck = true) const;
  void _initMeans();
  void _initVars();
  Id _getNVar(const Db* db);
  VectorInt _getVariableInterval(Id ivar) const;
  VectorInt _getDirectionInterval(Id idir) const;
  String _toStringByDirection(const AStringFormat* strfmt, Id idir) const;
  void _directionResize(Id idir);
  void _setDPasFromGrid(bool flag_grid);
  void _setFlagAsym();
  static VectorDouble _varsFromProportions(VectorDouble props);
  void _clearBiTargetCheck();
  void _addBiTargetCheck(ABiTargetCheck* abpc);
  void _setListBiTargetCheck();
  Id _getNBiPts() const { return static_cast<Id>(_bipts.size()); }
  Id _getBiPtsRank(Id idir, Id rank) const;

  Id _compute(Db* db,
              Id flag_sample,
              Id verr_mode,
              const Model* model,
              Id niter_UK,
              bool verbose);
  Id _calculateGeneral(Db* db,
                       Id flag_sample,
                       Id verr_mode);
  Id _calculateGenOnLine(Db* db, Id norder);
  Id _calculateGenOnGrid(DbGrid* db, Id norder);
  Id _calculateOnGrid(DbGrid* db);

  static Id _getRelativeSampleRank(Db* db, Id iech0);
  Id _updateUK(Db* db, Vario_Order* vorder);
  void _patchC00(Db* db, Id idir);
  Id _get_generalized_variogram_order();
  void _getStatistics(Db* db);
  Id _updateVerr(Db* db, Id idir, Vario_Order* vorder, Id verr_mode);
  static double _s(Db* db, Id iech, Id jech);
  double _g(Db* db, Id iech, Id jech) const;
  void _calculateBiasLocal(Db* db,
                           Id idir,
                           Id ilag,
                           Vario_Order* vorder,
                           Id ifirst,
                           Id ilast);
  void _calculateBiasGlobal(Db* db);
  double _getBias(Id iiech, Id jjech);

  void _calculateFromGeometry(Db* db, Id idir, Vario_Order* vorder);
  Id _calculateGeneralSolution1(Db* db, Id idir, const Id* rindex, Vario_Order* vorder);
  Id _calculateGeneralSolution2(Db* db, Id idir, const Id* rindex);
  Id _calculateOnGridSolution(DbGrid* db, Id idir);
  Id _calculateGenOnGridSolution(DbGrid* db, Id idir, Id norder);
  Id _calculateVarioVectSolution(Db* db, Id idir, Id ncomp, const Id* rindex);
  void _calculateOnLineSolution(Db* db, Id idir, Id norder);

  void _driftManage(Db* db);
  Id _driftEstimateCoefficients(Db* db);

  static void _printDebug(Id iech1,
                          Id iech2,
                          Id ivar,
                          Id jvar,
                          Id ilag,
                          double scale,
                          double value);
  void _centerCovariance(Db* db, Id idir);
  void _getVarioVectStatistics(Db* db, Id ncomp);
  void _rescale(Id idir);
  bool _isCompatible(const Db* db) const;
  static double _linear_interpolate(Id n,
                                    const VectorDouble& x,
                                    const VectorDouble& y,
                                    double x0);
  MatrixSquare _evalAverageDbIncr(Model* model,
                                  const Db& db,
                                  const VectorDouble& incr = VectorDouble(),
                                  const CovCalcMode* mode  = nullptr) const;

private:
  Id _nVar;
  VarioParam _varioparam;
  VectorDouble _means;
  VectorDouble _vars;
  bool _flagSample;
  Db* _db;
  VectorVectorDouble _sw;      /* Array for number of lags */
  VectorVectorDouble _gg;      /* Array for average variogram values */
  VectorVectorDouble _hh;      /* Array for average distance values */
  VectorVectorDouble _utilize; /* Array to mention if a lag is used or not */

  Id _biPtsPerDirection;
  std::vector<ABiTargetCheck*> _bipts;
  mutable bool _flagAsym;

  bool _verbose;
  bool _flag_UK;
  Id _niter_UK;

  VectorString _variableNames;
  mutable Model* _model; // Model pointer (not to be deleted) for drift removal
  mutable VectorDouble _BETA;
  mutable VectorDouble _DRFDIAG;
  mutable MatrixDense _DRFXA;
  mutable MatrixDense _DRFGX;
  mutable MatrixDense _DRFTAB;
  mutable MatrixSquare _DRFXGX;
};

GSTLEARN_EXPORT Vario_Order*
vario_order_manage(Id mode, Id flag_dist, Id size_aux, Vario_Order* vorder);

GSTLEARN_EXPORT Vario_Order* vario_order_final(Vario_Order* vorder, Id* npair);
GSTLEARN_EXPORT void vario_order_print(Vario_Order* vorder,
                                       Id idir_target,
                                       Id ipas_target,
                                       Id verbose);
GSTLEARN_EXPORT void vario_order_get_bounds(Vario_Order* vorder,
                                            Id idir,
                                            Id ilag,
                                            Id* ifirst,
                                            Id* ilast);
GSTLEARN_EXPORT void vario_order_get_indices(Vario_Order* vorder,
                                             Id ipair,
                                             Id* iech,
                                             Id* jech,
                                             double* dist);
GSTLEARN_EXPORT void vario_order_get_auxiliary(Vario_Order* vorder,
                                               Id ipair,
                                               char* aux_iech,
                                               char* aux_jech);
GSTLEARN_EXPORT Id vario_order_add(Vario_Order* vorder,
                                   Id iech,
                                   Id jech,
                                   void* aux_iech,
                                   void* aux_jech,
                                   Id ilag,
                                   Id idir,
                                   double dist);
} // namespace gstlrn
