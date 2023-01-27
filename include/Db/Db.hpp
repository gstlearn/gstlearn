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
#include "geoslib_d.h"

#include "Enum/ELoadBy.hpp"
#include "Enum/EStatOption.hpp"

#include "Db/PtrGeos.hpp"
#include "Basic/Grid.hpp"
#include "Basic/Limits.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/CSVformat.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/ICloneable.hpp"

class DbGrid;
class Polygons;

/**
 * Class containing a Data Set organized as a set of Isolated Points.
 */
class GSTLEARN_EXPORT Db: public AStringable, public ASerializable, public ICloneable
{
public:
  Db();
  Db(const Db& r);
  Db& operator=(const Db& r);
  virtual ~Db();

public:
  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// ICloneable interface
  IMPLEMENT_CLONING(Db)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface for Db
  virtual bool isGrid() const { return false; }
  virtual double getCoordinate(int iech, int idim, bool flag_rotate=true) const;
  virtual double getUnit(int idim = 0) const;
  virtual int getNDim() const;
  virtual bool mayChangeSampleNumber() const { return true; }
  virtual void resetDims(int ncol, int nech);

  static Db* createFromNF(const String& neutralFilename,
                           bool verbose = true);
  int resetFromSamples(int nech,
                       const ELoadBy& order = ELoadBy::fromKey("SAMPLE"),
                       const VectorDouble& tab = VectorDouble(),
                       const VectorString& names = VectorString(),
                       const VectorString& locatorNames = VectorString(),
                       int flag_add_rank = 1);
  int resetFromCSV(const String& filename,
                   bool verbose,
                   const CSVformat& csvfmt,
                   int ncol_max = -1,
                   int nrow_max = -1,
                   int flag_add_rank = 1);
  int resetFromBox(int nech,
                   const VectorDouble& coormin,
                   const VectorDouble& coormax,
                   int ndim = 2,
                   int seed = 321415,
                   int flag_add_rank = 1);
  int resetFromOnePoint(const VectorDouble &tab = VectorDouble(),
                        int flag_add_rank = 1);
  int resetSamplingDb(const Db* dbin,
                      double proportion = 0,
                      int number = 0,
                      const VectorString& names = VectorString(),
                      int seed = 23241,
                      bool verbose = false);

  static Db* create();
  static Db* createFromSamples(int nech,
                               const ELoadBy& order = ELoadBy::fromKey("SAMPLE"),
                               const VectorDouble& tab = VectorDouble(),
                               const VectorString& names = VectorString(),
                               const VectorString& locatorNames = VectorString(),
                               int flag_add_rank = 1);
  static Db* createFromCSV(const String& filename,
                           const CSVformat& csv = CSVformat(),
                           bool verbose = false,
                           int ncol_max = -1,
                           int nrow_max = -1,
                           int flag_add_rank = 1);
  static Db* createFromBox(int nech,
                           const VectorDouble& coormin,
                           const VectorDouble& coormax,
                           int seed = 43241,
                           bool flag_exact = true,
                           bool flag_repulsion = false,
                           double range = 0.,
                           double beta = 0.,
                           int flag_add_rank = 1);
  static Db* createFromOnePoint(const VectorDouble &tab = VectorDouble(),
                                int flag_add_rank = 1);
  static Db* createSamplingDb(const Db* dbin,
                              double proportion = 0.,
                              int number = 0,
                              const VectorString& names = VectorString(),
                              int seed = 23241,
                              bool verbose = false);
  static Db* createFromDbGrid(int nech,
                              DbGrid* dbgrid,
                              int seed = 432423,
                              bool flag_exact = true,
                              bool flag_repulsion = false,
                              double range = 0.,
                              double beta = 0.,
                              int flag_add_rank = 1);

  DbGrid* coveringDb(const VectorInt& nodes = VectorInt(),
                     const VectorDouble& dcell = VectorDouble(),
                     const VectorDouble& origin = VectorDouble(),
                     const VectorDouble& margin = VectorDouble()) const;

  const VectorDouble& getArrays() const { return _array; }

  String getNameByLocator(const ELoc& locatorType, int locatorIndex=0) const;
  String getNameByColIdx(int icol) const;
  String getNameByUID(int iuid) const;

  VectorString getName(const String& name) const;
  VectorString getNames(const VectorString& names) const;
  VectorString getNamesByLocator(const ELoc& locatorType) const;
  VectorString getNamesByColIdx(const VectorInt& icols) const;
  VectorString getNamesByUID(const VectorInt& iuids) const;
  VectorString getAllNames() const;

  void setName(const String& old_name, const String& name);
  void setName(const VectorString list, const String& name);
  void setNameByUID(int iuid, const String& name);
  void setNameByColIdx(int icol, const String& name);
  void setNameByLocator(const ELoc& locatorType, const String& name);

  inline int getUIDMaxNumber() const { return (int) _uidcol.size(); }
  inline int getColumnNumber() const { return _ncol; }
  double getColumnSize(bool useSel = false) const;
  int getSampleNumber(bool useSel = false) const;
  int getActiveSampleNumber() const;
  int getRankRelativeToAbsolute(int irel) const;
  int getRankAbsoluteToRelative(int iabs) const;

  VectorString expandNameList(const VectorString& names) const;
  VectorString expandNameList(const String& names) const;

  // Locator and UID methods

  void clearLocators(const ELoc& locatorType);
  void setLocatorByUID(int iuid,
                       const ELoc& locatorType = ELoc::fromKey("UNKNOWN"),
                       int locatorIndex = 0,
                       bool cleanSameLocator = false);
  void setLocatorByColIdx(int icol,
                          const ELoc& locatorType = ELoc::fromKey("UNKNOWN"),
                          int locatorIndex = 0,
                          bool cleanSameLocator = false);
  void setLocator(const String& names,
                  const ELoc& locatorType = ELoc::fromKey("UNKNOWN"),
                  int locatorIndex = 0,
                  bool cleanSameLocator = false);
  void setLocators(const VectorString &names,
                   const ELoc &locatorType = ELoc::fromKey("UNKNOWN"),
                   int locatorIndex = 0,
                   bool cleanSameLocator = false);
  void setLocatorsByUID(int number,
                        int iuid,
                        const ELoc& locatorType = ELoc::fromKey("UNKNOWN"),
                        int locatorIndex = 0,
                        bool cleanSameLocator = false);
  void setLocatorsByUID(const VectorInt& iuids,
                        const ELoc& locatorType = ELoc::fromKey("UNKNOWN"),
                        int locatorIndex = 0,
                        bool cleanSameLocator = false);
  void setLocatorsByColIdx(const VectorInt& icols,
                           const ELoc& locatorType = ELoc::fromKey("UNKNOWN"),
                           int locatorIndex = 0,
                           bool cleanSameLocator = false);
  void addColumnsByVVD(const VectorVectorDouble tab,
                       const String &radix,
                       const ELoc& locatorType,
                       int locatorIndex,
                       bool useSel = false,
                       double valinit = 0.,
                       int nvar = 1);
  int addColumns(const VectorDouble& tab,
                 const String& radix = "New",
                 const ELoc& locatorType = ELoc::fromKey("UNKNOWN"),
                 int locatorIndex = 0,
                 bool useSel = false,
                 double valinit = 0.,
                 int nvar = 1);
  int addColumnsByConstant(int nadd,
                           double valinit = 0.,
                           const String& radix = "New",
                           const ELoc& locatorType = ELoc::fromKey("UNKNOWN"),
                           int locatorIndex = 0,
                           int nechInit = 0);
  int addSelection(const VectorDouble& tab = VectorDouble(),
                   const String& name = "NewSel",
                   const String& combine = "set");
  int addSelectionByRanks(const VectorInt &ranks,
                          const String &name = "NewSel",
                          const String &combine = "set");
  int addSelectionByLimit(const String& testvar,
                          const Limits& limits = Limits(),
                          const String& name = "NewSel",
                          const String& combine = "set");
  int addSamples(int nadd, double valinit);
  int deleteSample(int e_del);
  void switchLocator(const ELoc& locatorTypein, const ELoc& locatorTypeout);
  int  getLastUID(int number = 0) const;
  String getLastName(int number = 0) const;

  int getColIdx(const String& name) const;
  int getColIdxByUID(int iuid) const;
  int getColIdxByLocator(const ELoc& locatorType, int locatorIndex=0) const;
  VectorInt getColIdxs(const String& name) const;
  VectorInt getColIdxs(const VectorString& names) const;
  VectorInt getColIdxsByUID(const VectorInt iuids) const;
  VectorInt getColIdxsByLocator(const ELoc& locatorType) const;

  void setColumn(const VectorDouble &tab,
                 const String &name,
                 const ELoc& locatorType = ELoc::fromKey("UNKNOWN"),
                 int locatorIndex = 0,
                 bool useSel = false);
  void setColumnByUIDOldStyle(const double* tab, int iuid, bool useSel = false);
  void setColumnByUID(const VectorDouble& tab, int iuid, bool useSel = false);
  void setColumnByColIdx(const VectorDouble& tab, int icol, bool useSel = false);
  void setColumnsByColIdx(const VectorDouble& tabs, const VectorInt& icols, bool useSel = false);
  void setColumnByColIdxOldStyle(const double* tab, int icol, bool useSel = false);
  void duplicateColumnByUID(int iuid_in, int iuid_out);


  VectorVectorDouble getItem(const VectorInt& rows,
                             const VectorString& colnames,
                             bool useSel = false) const;
  VectorVectorDouble getItem(const VectorInt& rows,
                             const String& colname,
                             bool useSel = false) const;
  VectorVectorDouble getItem(const VectorInt& rows,
                             const ELoc& locatorType,
                             bool useSel = false) const;
  VectorVectorDouble getItem(const VectorString& colnames,
                             bool useSel = false) const;
  VectorVectorDouble getItem(const String& colname,
                             bool useSel = false) const;
  VectorVectorDouble getItem(const ELoc& locatorType,
                             bool useSel = false) const;
  VectorString getItemNames(const VectorString& colnames);
  VectorString getItemNames(const String& colname);
  VectorString getItemNames(const ELoc& locatorType);

  int setItem(const VectorInt& rows,
              const VectorString& colnames,
              const VectorVectorDouble& values,
              bool useSel = false);
  int setItem(const VectorInt& rows,
              const ELoc& locatorType,
              const VectorVectorDouble& values,
              bool useSel = false);
  int setItem(const VectorString& colnames,
              const VectorVectorDouble& values,
              bool useSel = false);
  int setItem(const ELoc& locatorType,
              const VectorVectorDouble& values,
              bool useSel = false);
  int setItem(const VectorInt& rows,
              const String& colname,
              const VectorDouble& values,
              bool useSel = false);
  int setItem(const String& colname,
              const VectorDouble& values,
              bool useSel = false);

  bool getLocator(const String& name,
                  ELoc* ret_locatorType,
                  int* ret_locatorIndex) const;
  bool getLocatorByColIdx(int icol,
                          ELoc* ret_locatorType,
                          int* ret_locatorIndex) const;
  bool getLocatorByUID(int iuid,
                       ELoc* ret_locatorType,
                       int* ret_locatorIndex) const;
  VectorString getLocators(bool anyLocator = true,
                           const ELoc& locatorType = ELoc::fromKey("UNKNOWN")) const;
  int getLocatorNumber(const ELoc& locatorType) const;
  bool isUIDDefined(int iuid) const;

  int getUID(const String &name) const;
  int getUIDByColIdx(int icol) const;
  int getUIDByLocator(const ELoc& locatorType, int locatorIndex=0) const;

  VectorInt getUIDs(const VectorString& names) const;
  VectorInt getUIDsByLocator(const ELoc& locatorType) const;
  VectorInt getUIDsByColIdx(const VectorInt& icols) const;
  VectorInt getAllUIDs() const;

  int getFaciesNumber(void) const;
  bool hasLocatorDefined(const String& name, const ELoc& locatorType, int locatorIndex=0) const;

  // Accessing elements of the contents

  VectorDouble getSampleCoordinates(int iech) const;
  void getSampleCoordinates(int iech, VectorDouble& coor) const;
  VectorDouble getSampleLocators(const ELoc& locatorType, int iech) const;

  void   getCoordinatesInPlace(int iech, VectorDouble& coor, bool flag_rotate = true) const;
  VectorDouble getCoordinates(int idim, bool useSel = false, bool flag_rotate = true) const;
  VectorVectorDouble getAllCoordinates(bool useSel = false) const;
  void   setCoordinate(int iech, int idim, double value);

  double getDistance1D(int iech, int jech, int idim, bool flagAbs = false) const;
  double getDistance(int iech, int jech) const;
  int    getDistanceVec(int iech, int jech, VectorDouble& dd, const Db* db2 = nullptr) const;

  double getValue(const String& name, int iech) const;
  void   setValue(const String& name, int iech, double value);

  double getArray(int iech, int iuid) const;
  void   setArray(int iech, int iuid, double value);
  void   updArray(int iech, int iuid, int oper, double value);
  VectorDouble getArray(int iuid, bool useSel = false) const;
  VectorDouble getArrayBySample(int iech) const;
  void setArrayBySample(int iech, const VectorDouble& vec);

  int    getFromLocatorNumber(const ELoc& locatorType) const;
  double getFromLocator(const ELoc& locatorType, int iech, int locatorIndex=0) const;
  void   setFromLocator(const ELoc& locatorType,
                        int iech,
                        int locatorIndex,
                        double value);

  double getValueByColIdx(int iech, int icol) const;
  VectorDouble getValuesByColIdx(const VectorInt &iechs,
                                 const VectorInt &icols,
                                 bool bySample = false) const;
  void   setValueByColIdx(int iech, int icol, double value);
  void   setValuesByColIdx(const VectorInt &iechs,
                           const VectorInt &icols,
                           const VectorDouble &values,
                           bool bySample = false);
  int    getVariableNumber() const;
  bool   hasVariable() const;
  double getVariable(int iech, int item) const;
  void   setVariable(int iech, int item, double value);
  void   updVariable(int iech, int item, int oper, double value);
  bool   isVariableNumberComparedTo(int nvar, int compare = 0) const;
  bool   isIsotopic(int iech, int nvar_max = -1) const;
  bool   isAllUndefined(int iech) const;

  int    getLowerIntervalNumber() const;
  bool   hasLowerInterval() const;
  double getLowerInterval(int iech, int item) const;
  void   setLowerInterval(int iech, int item, double rklow);

  int    getUpperIntervalNumber() const;
  bool   hasUpperInterval() const;
  double getUpperInterval(int iech, int item) const;
  void   setUpperInterval(int iech, int item, double rkup);
  void   setIntervals(int iech, int item, double rklow, double rkup);

  int    getIntervalNumber() const;

  int    getLowerBoundNumber() const;
  bool   hasLowerBound() const;
  double getLowerBound(int iech, int item) const;
  void   setLowerBound(int iech, int item, double lower);

  int    getUpperBoundNumber() const;
  bool   hasUpperBound() const;
  double getUpperBound(int iech, int item) const;
  void   setUpperBound(int iech, int item, double upper);
  void   setBounds(int iech, int item, double lower, double upper);
  VectorDouble getWithinBounds(int item, bool useSel = false) const;

  int    getGradientNumber() const;
  bool   hasGradient() const;
  double getGradient(int iech, int item) const;
  VectorDouble getGradients(int item, bool useSel = false) const;
  void   setGradient(int iech, int item, double value);

  int    getTangentNumber() const;
  bool   hasTangent() const;
  double getTangent(int iech, int item) const;
  VectorDouble getTangents(int item, bool useSel = false) const;
  void   setTangent(int iech, int item, double value);

  int    getProportionNumber() const;
  bool   hasProportion() const;
  double getProportion(int iech, int item) const;
  void   setProportion(int iech, int item, double value);

  bool   hasSelection() const;
  int    getSelection(int iech) const;
  void   setSelection(int iech, int value);
  VectorDouble getSelection(void) const;

  bool   hasWeight() const;
  double getWeight(int iech) const;
  void   setWeight(int iech, double value);
  VectorDouble getWeight(bool useSel = false) const;

  int    getExternalDriftNumber() const;
  bool   hasExternalDrift() const;
  double getExternalDrift(int iech, int item) const;
  void   setExternalDrift(int iech, int item, double value);

  int    getBlockExtensionNumber() const;
  bool   hasBlockExtension() const;
  double getBlockExtension(int iech, int item) const;
  void   setBlockExtension(int iech, int item, double value);

  int    getCodeNumber() const;
  bool   hasCode() const;
  double getCode(int iech) const;
  void   setCode(int iech, double value);
  VectorDouble getCodeList(void);

  int    getVarianceErrorNumber() const;
  bool   hasVarianceError() const;
  double getVarianceError(int iech, int item) const;
  void   setVarianceError(int iech, int item, double value);

  bool   hasDomain() const;
  int    getDomain(int iech) const;
  void   setDomain(int iech, int value);

  int    getDipDirectionNumber() const;
  bool   hasDipDirection() const;
  double getDipDirection(int iech) const;
  void   setDipDirection(int iech, double value);

  int    getDipAngleNumber() const;
  bool   hasDipAngle() const;
  double getDipAngle(int iech) const;
  void   setDipAngle(int iech, double value);

  int    getObjectSizeNumber() const;
  bool   hasObjectSize() const;
  double getObjectSize(int iech) const;
  void   setObjectSize(int iech, double value);

  int    getBorderUpNumber() const;
  bool   hasBorderUp() const;
  double getBorderUp(int iech) const;
  void   setBorderUp(int iech, double value);

  int    getBorderDownNumber() const;
  bool   hasBorderDown() const;
  double getBorderDown(int iech) const;
  void   setBorderDown(int iech, double value);

  int    getDateNumber() const;
  bool   hasDate() const;
  double getDate(int iech) const;
  void   setDate(int iech, double value);

  int getSimvarRank(int isimu, int ivar, int icase, int nbsimu, int nvar);
  double getSimvar(const ELoc& locatorType,
                   int iech,
                   int isimu,
                   int ivar,
                   int icase,
                   int nbsimu,
                   int nvar) const;
  void setSimvar(const ELoc& locatorType,
                 int iech,
                 int isimu,
                 int ivar,
                 int icase,
                 int nbsimu,
                 int nvar,
                 double value);
  void updSimvar(const ELoc& locatorType,
                 int iech,
                 int isimu,
                 int ivar,
                 int icase,
                 int nbsimu,
                 int nvar,
                 int oper,
                 double value);

  bool isActive(int iech) const;
  bool isActiveAndDefined(int iech, int item) const;
  int  getActiveAndDefinedNumber(int item) const;
  int  getActiveAndDefinedNumber(const String& name) const;
  VectorBool getMaskArray() const;

  VectorInt getSortArray() const;
  double getCosineToDirection(int iech1, int iech2, const VectorDouble& codir) const;

  VectorDouble getColumn(const String &name,
                         bool useSel = false,
                         bool flagCompress = true) const;
  VectorDouble getColumnByUID(int iuid,
                              bool useSel = false,
                              bool flagCompress = true) const;
  VectorDouble getColumnByLocator(const ELoc &locatorType,
                                  int locatorIndex = 0,
                                  bool useSel = false,
                                  bool flagCompress = true) const;
  VectorDouble getColumnByColIdx(int icol,
                                 bool useSel = false,
                                 bool flagCompress = true) const;

  VectorDouble getAllColumns(bool useSel = false,
                             bool flagCompress = true) const;
  VectorDouble getColumns(const VectorString &names = VectorString(),
                          bool useSel = false,
                          bool flagCompress = true) const;
  VectorVectorDouble getColumnsAsVVD(const VectorString &names = VectorString(),
                                     bool useSel = false,
                                     bool flagCompress = true) const;
  MatrixRectangular getColumnsAsMatrix(const VectorString &names,
                                       bool useSel = false,
                                       bool flagCompress = true) const;
  VectorDouble getColumnsByColIdx(const VectorInt &icols = VectorInt(),
                                  bool useSel = false,
                                  bool flagCompress = true) const;
  VectorDouble getColumnsByColIdxInterval(int icol_beg,
                                          int icol_end,
                                          bool useSel = false,
                                          bool flagCompress = true) const;

  VectorDouble getColumnsByLocator(const ELoc &locatorType,
                                   bool useSel = false,
                                   bool flagCompress = true) const;
  VectorDouble getColumnsByUID(const VectorInt &iuids,
                               bool useSel = false,
                               bool flagCompress = true) const;
  VectorDouble getColumnsByUIDRange(int iuid_beg,
                                    int iuid_end,
                                    bool useSel = false,
                                    bool flagCompress = true) const;
  void setAllColumns(const VectorVectorDouble& tabs, bool useSel = false);

  VectorDouble getFFFFs(const VectorString &names = VectorString(),
                        bool useSel = false) const;

  void deleteColumn(const String& name);
  void deleteColumnByUID(int iuid_del);
  void deleteColumnByColIdx(int icol_del);

  void deleteColumns(const VectorString& names);
  void deleteColumnsByLocator(const ELoc& locatorType);
  void deleteColumnsByUID(const VectorInt& iuids);
  void deleteColumnsByColIdx(const VectorInt& icols);

  VectorDouble getExtrema(int idim, bool useSel = false) const;
  VectorVectorDouble getExtremas(bool useSel = false) const;
  VectorDouble getCoorMinimum(bool useSel = false) const;
  VectorDouble getCoorMaximum(bool useSel = false) const;
  double getExtension(int idim, bool useSel = false) const;
  double getExtensionDiagonal(bool useSel = false) const;
  double getCenter(int idim, bool useSel = false) const;
  VectorDouble getCenters(bool useSel = false) const;
  void getExtensionInPlace(VectorDouble &mini, VectorDouble &maxi);

  double getMinimum(const String& name, bool useSel = false) const;
  double getMaximum(const String& name, bool useSel = false) const;
  VectorDouble getRange(const String& name, bool useSel = false) const;
  double getMean(const String& name, bool useSel = false) const;
  double getVariance(const String& name, bool useSel = false) const;
  double getStdv(const String& name, bool useSel = false) const;
  double getCorrelation(const String& name1, const String& name2,bool useSel = false) const;

  bool hasSameDimension(const Db* dbaux) const;
  bool hasLargerDimension(const Db* dbaux) const;

  // Functions for checking validity of parameters

  bool isColIdxValid(int icol) const;
  bool isUIDValid(int iuid) const;
  bool isSampleIndexValid(int iech) const;
  bool isSampleIndicesValid(const VectorInt& iechs, bool useSel = false) const;
  bool isLocatorIndexValid(const ELoc& locatorType, int locatorIndex) const;
  bool isDimensionIndexValid(int idim) const;

  void combineSelection(VectorDouble& sel, const String& combine = "set") const;

  void generateRank(const String& radix = "rank");

  VectorInt shrinkToValidRows(const VectorInt& rows);
  VectorInt shrinkToValidCols(const VectorInt& cols);

  // Statistics
  VectorDouble statistics(const VectorString& names,
                          const std::vector<EStatOption>& opers = EStatOption::fromKeys({"MEAN"}),
                          bool flagIso = true,
                          bool flagStoreInDb = false,
                          bool verbose = true,
                          double vmin = TEST,
                          double vmax = TEST,
                          double proba = TEST,
                          const String& title = "",
                          const NamingConvention& namconv = NamingConvention("Stats"));
  VectorDouble statisticsByLocator(const ELoc& locatorType,
                                   const std::vector<EStatOption>& opers = EStatOption::fromKeys({"MEAN"}),
                                   bool flagIso = true,
                                   bool flagStoreInDb = false,
                                   bool verbose = true,
                                   double vmin = TEST,
                                   double vmax = TEST,
                                   double proba = TEST,
                                   const String& title = "",
                                   const NamingConvention& namconv = NamingConvention("Stats"));
  VectorDouble statisticsByUID(const VectorInt& iuids,
                               const std::vector<EStatOption>& opers = EStatOption::fromKeys({"MEAN"}),
                               bool flagIso = true,
                               bool flagStoreInDb = false,
                               bool verbose = true,
                               double proba = TEST,
                               double vmin = TEST,
                               double vmax = TEST,
                               const String& title = "",
                               const NamingConvention& namconv = NamingConvention("Stats"));
  VectorDouble statisticsMulti(const VectorString& names,
                               bool flagIso = true,
                               bool verbose = false,
                               const String& title = "");
  VectorDouble statisticsMultiByUID(const VectorInt& iuids,
                                    bool flagIso = true,
                                    bool verbose = false,
                                    const String& title = "");
  bool areSame(const String& name1,
               const String& name2,
               double eps = EPSILON3,
               bool useSel = true,
               bool verbose = false);

  VectorInt filter(const String& name,
                   const Interval& interval,
                   int belowRow = ITEST,
                   int aboveRow = ITEST) const;

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os,bool verbose = false) const override;
  String _getNFName() const override { return "Db"; }

  void _clear();
  void _createRank(int icol = 0);
  void _loadData(const VectorDouble& tab,
                 const VectorString& names,
                 const VectorString& locatorNames,
                 const ELoadBy& order,
                 int shift);
  void _loadData(const ELoadBy& order, int flag_add_rank, const VectorDouble& tab);
  void _defineDefaultNames(int shift, const VectorString& names);
  void _defineDefaultLocators(int shift, const VectorString& locatorNames);
  void _setNameByColIdx(int icol, const String& name);
  String _toStringCommon(const AStringFormat *strfmt) const;
  String _summaryString(void) const;

private:
  const VectorInt& _getUIDcol() const { return _uidcol; }
  const VectorString _getNames() const { return _colNames; }
  int _getUIDcol(int iuid) const;
  int _getAddress(int iech, int icol) const;
  void _columnInit(int ncol, int icol0, double valinit);
  double _updateValue(int oper, double oldval, double value);
  String _summaryVariables(void) const;
  String _summaryExtensions(void) const;
  String _summaryStats(VectorInt cols, int mode = 1, int maxNClass = 50) const;
  String _summaryLocators(void) const;
  String _summaryUIDs(void) const;
  String _summaryArrays(VectorInt cols, bool useSel = true) const;

  void _defineDefaultLocatorsByNames(int shift, const VectorString& names);
  int  _getSimrank(int isimu, int ivar, int icase, int nbsimu, int nvar) const;
  VectorInt _getUIDsBasic(const VectorString& names) const;

  int _getLastColumn(int number = 0) const;

  int _findColumnInLocator(const ELoc& locatorType, int icol) const;
  int _findUIDInLocator(const ELoc& locatorType, int iuid) const;
  String _getLocatorNameByColIdx(int icol) const;
  VectorInt _ids(const String& name, bool flagOne, bool verbose = true) const;
  VectorInt _ids(const VectorString& names, bool flagOne, bool verbose = true) const;
  VectorInt _ids(const ELoc& locatorType, bool flagOne, bool verbose = true) const;
  VectorInt _ids(const VectorInt& iuids, bool flagOne, bool verbose = true) const;

  VectorDouble _getItem(const String& exp_name,
                        bool useSel,
                        const VectorInt& rows) const;
  void _setItem(const String& name,
                const VectorInt& rows,
                const VectorDouble& values);
  void _setItem(const String& name, bool useSel, const VectorDouble& values);
  bool _isValidCountRows(const VectorInt& rows,
                         bool useSel,
                         const VectorDouble& values) const;
  bool _isValidCountRows(bool useSel, const VectorDouble& values) const;
  VectorString _getVarNames(const VectorString& colnames,
                            int expectedVarCount);

  // Higher level methods
  bool _isCountValid(const VectorInt iuds, bool flagOne, bool verbose = true) const;

private:
  int _ncol;                 //!< Number of Columns of data
  int _nech;                 //!< Number of samples
  VectorDouble _array;       //!< Array of values
  VectorInt _uidcol;         //!< UID to Column
  VectorString _colNames;    //!< Names of the variables
  std::map<ELoc,PtrGeos> _p; //!< Locator characteristics
};
