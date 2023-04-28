/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_d.h"

#include "Enum/ELoadBy.hpp"
#include "Enum/EStatOption.hpp"

#include "Db/PtrGeos.hpp"
#include "Basic/Grid.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/CSVformat.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/Limits.hpp"

class DbGrid;
class Interval;

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
  virtual void getCoordinatesPerSampleInPlace(int iech, VectorDouble& coor, bool flag_rotate = true) const;
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
                   double extend = 0.,
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
  int resetReduce(const Db *dbin,
                  const VectorString &names = VectorString(),
                  const VectorInt &ranks = VectorInt(),
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
                           double extend = 0.,
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
  static Db* createReduce(const Db *dbin,
                          const VectorString &names = VectorString(),
                          const VectorInt &ranks = VectorInt(),
                          bool verbose = false);
  static Db* createFillRandom(int ndat,
                              int ndim = 2,
                              int nvar = 1,
                              int nfex = 0,
                              double varmax = 0.,
                              double selRatio = 0.,
                              const VectorDouble heteroRatio = VectorDouble(),
                              int seed = 124234);

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

  int getSampleNumber(bool useSel = false) const;
  int getActiveSampleNumber() const;
  int getRankRelativeToAbsolute(int irel) const;
  int getRankAbsoluteToRelative(int iabs) const;

  VectorString expandNameList(const VectorString& names) const;
  VectorString expandNameList(const String& names) const;
  VectorString identifyNames(const VectorString& names) const;

  // Locator and UID methods

  void clearLocators(const ELoc& locatorType);
  void clearSelection() { clearLocators(ELoc::SEL); }
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
                       int locatorIndex = 0,
                       bool useSel = false);
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
  int addSelectionFromDbByConvexHull(Db* db,
                                     double dilate = 0.,
                                     bool verbose = false,
                                     const NamingConvention &namconv = NamingConvention("Hull", true, true, true,
                                                                                        ELoc::fromKey("SEL")));

  int addSamples(int nadd, double valinit);
  int deleteSample(int e_del);
  int deleteSamples(const VectorInt& e_dels);
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
  VectorDouble getArrayByUID(int iuid, bool useSel = false) const;
  VectorDouble getArrayBySample(int iech) const;
  void setArrayBySample(int iech, const VectorDouble& vec);

  int    getFromLocatorNumber(const ELoc& locatorType) const;
  double getFromLocator(const ELoc& locatorType, int iech, int locatorIndex=0) const;
  void   setFromLocator(const ELoc& locatorType,
                        int iech,
                        int locatorIndex,
                        double value);

  double getValueByColIdx(int iech, int icol) const;
  void   setValueByColIdx(int iech, int icol, double value);
  VectorDouble getValuesByNames(const VectorInt& iechs,
                                const VectorString& names,
                                bool bySample = false) const;
  VectorDouble getValuesByColIdx(const VectorInt &iechs,
                                 const VectorInt &icols,
                                 bool bySample = false) const;
  void   setValuesByNames(const VectorInt &iechs,
                          const VectorString &names,
                          const VectorDouble &values,
                          bool bySample = false);
  void   setValuesByColIdx(const VectorInt &iechs,
                           const VectorInt &icols,
                           const VectorDouble &values,
                           bool bySample = false);

  /**
   * \defgroup DB Db: Numerical Data Base
   *
   * Here are the implementation of several functions regarding the manipulation of the Numerical Data Base
   *
   * These operations are generic: they are available for any class derived from the Db one (such as DbGrid).
   *
   **/

  /** @addtogroup DB_0 Getting and Setting functions by Locator
   * \ingroup DB
   *
   * Various functions for accessing fields of the Db using the **locator** designation.
   * They use the argument 'loctype' which refers to the Locator type (see ELoc enumeration).
   * In most cases, they also refer to 'item' i.e. the rank (0 based) for the target locator.
   *
   * @param loctype Target locator
   * @param iech    Target sample (0 based)
   * @param item    Rank of the 'loctype' locator (0 based)
   * @param oper    Type of operation
   * \li                 0 : New = New + Old
   * \li                 1 : New = New * Old
   * \li                 2 : New = New - Old
   * \li                 3 : New = Old / New
   * \li                 4 : New = New (only if old is defined)
   * \li                 5 : New = MAX(New, Old)
   * \li                 6 : New = MIN(New, Old)
   * @param value   Assigned value
   *  @{
   */
  int    getLocNumber(const ELoc& loctype) const;
  bool   hasLocVariable(const ELoc& loctype) const;
  double getLocVariable(const ELoc& loctype, int iech, int item) const;
  void   setLocVariable(const ELoc& loctype, int iech, int item, double value);
  void   updLocVariable(const ELoc& loctype, int iech, int item, int oper, double value);
  /**@}*/

  bool   isVariableNumberComparedTo(int nvar, int compare = 0) const;
  bool   isIsotopic(int iech, int nvar_max = -1) const;
  bool   isAllUndefined(int iech) const;

  void   setInterval(int iech, int item, double rklow = TEST, double rkup = TEST);
  int    getIntervalNumber() const;
  void   setBound(int iech, int item, double lower = TEST, double upper = TEST);
  VectorDouble getWithinBounds(int item, bool useSel = false) const;
  VectorDouble getGradient(int item, bool useSel = false) const;
  VectorDouble getTangent(int item, bool useSel = false) const;
  VectorDouble getCodeList(void);

  int    getSelection(int iech) const;
  VectorDouble getSelections(void) const;
  VectorInt getSelectionRanks() const;

  double getWeight(int iech) const;
  VectorDouble getWeights(bool useSel = false) const;

  /** @addtogroup DB_1 Variable designation (used for simulations in particular)
   * \ingroup DB
   *
   * These functions allow designation of columns which contain the results of one simulation
   * for one variable in particular.
   *
   * @param locatorType Target locator type
   * @param iech Rank of the target sample
   * @param isimu Rank of the simulation (0-based)
   * @param ivar Rank of the variable (0-based)
   * @param icase Rank of the GRF / PGS
   * @param nbsimu Number of simulations
   * @param nvar Number of variables
   * @param value Value to be assigned
   *  @{
   */
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
  /**@}*/

  bool isActive(int iech) const;
  bool isActiveDomain(int iech) const;
  bool isActiveAndDefined(int iech, int item) const;
  int  getActiveAndDefinedNumber(int item) const;
  int  getActiveAndDefinedNumber(const String& name) const;
  VectorBool getMaskArray() const;

  VectorInt getSortArray() const;
  double getCosineToDirection(int iech1, int iech2, const VectorDouble& codir) const;

  /** @addtogroup DB_2 Reading one or several Columns
   * \ingroup DB
   *
   * The **column** refers to one element of the Db (which can be viewed as an Excel spread sheet).
   * Each variable stands as a column of this table: it is also attached a 'name' (which will serve
   * as the name of the variable) and a possible 'locator' (which characterizes the role of the
   * variable, e.g; coordinate, variable, code, ...).
   * These functions can refer to a single column or to several of them.
   * The columns can be referred to by the variable name, the column index, the internal Id (UID) or the locator.
   * @param useSel Option when reading a masked sample:
   * \li TRUE: the contents of the masked samples is set to TEST
   * \li FALSE: the masked samples are returned with no impact of the selection
   * @param flagCompress Option when reading a masked sample:
   * \li TRUE: the returned array is compressed to the only non-masked samples
   * \li FALSE: the returned array is not compressed
   *
   * @param name Name of the target column
   * @param names Vector of target variable names
   * @param locatorType Type of target locator
   * @param iuids Vector of target user-identified ranks
   * @param icols Vector of Column ranks
   * @param icol_beg Lower bound of the rank interval (included)
   * @param icol_end Upper bound of the rank interval (excluded)
   * @param iuid_beg Lower bound of the user-identification interval (included)
   * @param iuid_end Upper bound of the user-identification interval (excluded)
   * @param locatorType Type of the target locator
   * @param locatorIndex Rank of the item (0-based) for the target locator
   * @param flagCompress When True, the masked values are skipped
   *  @{
   */
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
  /**@}*/

  void setAllColumns(const VectorVectorDouble& tabs);

  /** @addtogroup DB_3 Deleting one or several Columns
   * \ingroup DB
   *
   * These Columns are defined by their names, column number of user-identification rank
   *
   * @param name Name of the variable to be deleted
   * @param names Vector of variable names to be deleted
   * @param icol_del Column number of the variable to be deleted
   * @param icols Vector of Column ranks for the variables to be deleted
   * @param iuid_del User-identification rank for the variable to be deleted
   * @param iuids Vector of user-identification ranks for variables to be deleted
   * @param locatorType Locator of the variables to be deleted
   * @{
   */
  void deleteColumn(const String& name);
  void deleteColumnByUID(int iuid_del);
  void deleteColumnByColIdx(int icol_del);

  void deleteColumns(const VectorString& names);
  void deleteColumnsByLocator(const ELoc& locatorType);
  void deleteColumnsByUID(const VectorInt& iuids);
  void deleteColumnsByColIdx(const VectorInt& icols);
  /**@}*/

  /** @addtogroup DB_4 Calculating Spatial characteristics on the Db
   * \ingroup DB
   *
   * @param idim Rank of the target space dimension (0 based)
   * @param useSel When TRUE, the characteristics are derived from the only active samples
   * @param mini Vector of minimum values (modified by this function)
   * @param maxi Vector of maximum values (modified by this function)
   *
   *  @{
   */
  VectorDouble getExtrema(int idim, bool useSel = false) const;
  VectorVectorDouble getExtremas(bool useSel = false) const;
  VectorDouble getCoorMinimum(bool useSel = false) const;
  VectorDouble getCoorMaximum(bool useSel = false) const;
  double getExtension(int idim, bool useSel = false) const;
  double getExtensionDiagonal(bool useSel = false) const;
  double getCenter(int idim, bool useSel = false) const;
  VectorDouble getCenters(bool useSel = false) const;
  void getExtensionInPlace(VectorDouble &mini, VectorDouble &maxi, bool useSel = false);
  /**@}*/

  /** @addtogroup DB_5 Calculating basic Statistics
   * \ingroup DB
   *
   * Calculate some basic statistics on the active samples of variables stored in a Db.
   *
   * @param name Target variable name
   * @param name1 First  target variable name
   * @param name2 Second  target variable name
   *
   * @param useSel When TRUE, the statistics are derived from the only active samples
   *
   *  @{
   */
  double getMinimum(const String& name, bool useSel = false) const;
  double getMaximum(const String& name, bool useSel = false) const;
  VectorDouble getRange(const String& name, bool useSel = false) const;
  double getMean(const String& name, bool useSel = false) const;
  double getVariance(const String& name, bool useSel = false) const;
  double getStdv(const String& name, bool useSel = false) const;
  double getCorrelation(const String& name1, const String& name2,bool useSel = false) const;
  /**@}*/

  bool hasSameDimension(const Db* dbaux) const;
  bool hasLargerDimension(const Db* dbaux) const;

  /** @addtogroup DB_6 Checking validity for various parameters
   * \ingroup DB
   *
   * These functions are used in order to check that the arguments are valid
   * (such as the sample rank, the locator type, the user-designation rank)
   *
   * @param icol Column rank to be checked
   * @param iuid User-designated rank
   * @param iech Sample rank to be checked
   * @param idim Space rank to be checked
   * @param iechs Vector of sample ranks to be checked
   * @param useSel When TRUE, the rank corresponds to the *active* sample
   * @param locatorType Type of the Locator
   * @param locatorIndex Rank of the locator (0-based)
   *
   *  @{
   */
  bool isColIdxValid(int icol) const;
  bool isUIDValid(int iuid) const;
  bool isSampleIndexValid(int iech) const;
  bool isSampleIndicesValid(const VectorInt& iechs, bool useSel = false) const;
  bool isLocatorIndexValid(const ELoc& locatorType, int locatorIndex) const;
  bool isDimensionIndexValid(int idim) const;
  /**@}*/

  void combineSelection(VectorDouble& sel, const String& combine = "set") const;

  void generateRank(const String& radix = "rank");

  VectorInt shrinkToValidRows(const VectorInt& rows);
  VectorInt shrinkToValidCols(const VectorInt& cols);

  /** @addtogroup DB_7 Calculating several statistics in Db
   * \ingroup DB
   *
   * These functions are meant to calculate several statistics on a set of target variables per sample.
   * The resulting values are stored in variables newly created in the same Db.
   *
   * @param names Vector of target variable names
   * @param iuids Vector of user-designation ranks
   * @param locatorType Target Locator
   * @param opers Vector of operations to be performed
   * @param flagIso The statistics are calculated only for samples where all target variables have defined values
   * @param flagStoreInDb When TRUE, the results are stored in the Db; otherwise the statistics are returned
   * @param verbose Verbose flag
   * @param proba              For 'quant': the quantile for this probability is calculated
   * @param vmin               For 'prop', 'T', 'Q', 'M', 'B': defines the lower bound of the interval to work in
   * @param vmax               For 'prop', 'T', 'Q', 'M', 'B': defines the upper bound of the interval to work in
   * @param title              If verbose, the title of the printed statistics.
   * @param namconv            Naming Convention used as a radix for the variables newly created in the Db
   * (only used when 'flagStoreInDb' is TRUE)
   *
   * @return If 'flagStoreInDb' is FALSE, the function returns a vector containing the statistics.
   * @return If there is more than one operator and more than one variable, the statistics are ordered first by variables
   * (all the statistics of the first variable, then all the statistics of the second variable...).
   *
   *  @{
   */
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
  /**@}*/

  /** @addtogroup DB_8 Calculating correlations on variables of a Db
   * \ingroup DB
   *
   * These functions calculate the correlation matrix based on a set of variables contained in a Db.
   * Although the result stands as a matrix, they are returned as a Vector.
   *
   * @param names Vector of target variable names
   * @param iuids Vector of user-designation ranks
   * @param flagIso The statistics are calculated only for samples where all target variables have defined values
   * @param verbose Verbose flag
   * @param title If verbose, the title of the printed statistics.
   *
   * @return These functions return a vector containing the correlation matrix.
   *  @{
   */
  VectorDouble statisticsMulti(const VectorString& names,
                               bool flagIso = true,
                               bool verbose = false,
                               const String& title = "");
  VectorDouble statisticsMultiByUID(const VectorInt& iuids,
                                    bool flagIso = true,
                                    bool verbose = false,
                                    const String& title = "");
  /**@}*/

  bool areSame(const String& name1,
               const String& name2,
               double eps = EPSILON3,
               bool useSel = true,
               bool verbose = false);

  VectorInt filter(const String& name,
                   const Interval& interval,
                   int belowRow = ITEST,
                   int aboveRow = ITEST) const;

  VectorInt getSampleRanks() const;

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

  inline const PtrGeos& _getPtrGeos(const ELoc& locatorType) const { return _p.at(locatorType.getValue()); }

protected:
  void _defineVariableAndLocators(const Db* dbin, const VectorString& names);
  void _loadValues(const Db* db, const VectorString& names, const VectorInt& ranks);

private:
  int _ncol;                 //!< Number of Columns of data
  int _nech;                 //!< Number of samples
  VectorDouble _array;       //!< Array of values
  VectorInt _uidcol;         //!< UID to Column
  VectorString _colNames;    //!< Names of the variables
  std::map<int,PtrGeos> _p; //!< Locator characteristics
};
