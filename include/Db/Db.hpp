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

#include "Space/ASpace.hpp"
#include "gstlearn_export.hpp"

#include "Enum/ELoadBy.hpp"
#include "Enum/EStatOption.hpp"

#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/CSVformat.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/Limits.hpp"
#include "Basic/NamingConvention.hpp"
#include "Db/PtrGeos.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/Table.hpp"

namespace gstlrn
{
  class ASpace;
  class DbGrid;
  class Interval;
  class SpacePoint;
  class SpaceTarget;
  class ModelGeneric;

  /**
   * \brief
   * Class containing the Data Information.
   *
   * This Data Set benefits from the comparison to an Excel spread sheet:
   * it can be considered as a Data Frame with a number of rows and a number of columns.
   * The columns will correspond to **variables** and the rows to **samples**.
   *
   * Notes:
   * - For short, this Data Base organization is referred to as **Db**,
   * - At any time, variables can be added, renamed or deleted,
   * - In the current version, this Data Set is currently limited to numerical contents.
   * The data frame is necessarily completely filled with values: therefore, a specific value
   * stands for a missing value (printed as **NA**);
   *
   * Moreover each variable may be assigned a *role* (or functionality) in the rest of a script.
   * This role is defined by a **locator**: for example, a variable can serve as a coordinate,
   * or a target variable. The locators can be viewed in the following list (see ELoc.hpp).
   * Note that a variable which does not play any role in particular may be assigned to an idle
   * role (locator NA);
   * This locator may be modified or cancelled (by the user) at any time.
   *
   * There are two status for these locators: **unique** or **multiple**.
   *
   * *Unique*. Only one variable can be assigned this unique locator: this is the
   * case for the selection (*SEL*) as there may be only one current selection activated at a time.
   *
   * *Multiple*. Several variable may be assigned the same locator: this is the case for the
   * coordinates (*X*). In that case, the locator name is followed by its rank (1-based).
   * Then the first coordinate will correspond to the locator *X1*, the second coordinate to *X2*, ...
   * Note that:
   * - for a given locator name, the ranks are always consecutive between 1 and N (if it
   * happens that you delete X3, the locator X4 is automatically modified into X3, X5 into X4, and
   * so on up to X{N} into X{N-1}.
   * - there is no limitation in the number of ranks for a given locator name.
   *
   * Each variable (or column) can be designated:
   * - by its name (unique in the Data Base) or
   * - by its locator (name and rank) or
   * - by its column index (0-based): this designation mode is dangerous (and not recommended) as the index may change over time.
   */
  class GSTLEARN_EXPORT Db: public AStringable,
                            public ASerializable,
                            public ICloneable
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
    String toString(const AStringFormat* strfmt = nullptr) const override;

    /// ASerializable interface
    String getNFName() const override { return "Db"; }
#ifdef HDF5
    bool deserializeH5(H5::Group& grp) override;
    bool serializeH5(H5::Group& grp) const override;
#endif

    /// Interface for Db
    virtual bool isGrid() const { return false; }

    virtual bool isLine() const { return false; }

    virtual bool isMesh() const { return false; }

    virtual double
      getCoordinate(Id iech, Id idim, bool flag_rotate = true) const;
    virtual void getCoordinatesInPlace(
      VectorDouble& coor,
      Id iech,
      bool flag_rotate = true) const;

    virtual double getUnit(Id idim = 0) const;
    virtual Id getNDim() const;

    virtual bool mayChangeSampleNumber() const { return true; }

    virtual void resetDims(Id ncol, Id nech);

    virtual bool isConsistent() const { return true; };

    /**
     * \defgroup DB Db: Numerical Data Base
     *
     **/

    /** @addtogroup DB_Reset Reset the contents of an already existing Db
     * \ingroup DB
     *
     * All methods enabling to Reset the contents of an already existing Db.
     *
     * They clean the initial contents and replace it by the new one.
     *  @{
     */
    Id resetFromSamples(
      Id nech,
      const ELoadBy& order = ELoadBy::fromKey("SAMPLE"),
      const VectorDouble& tab = VectorDouble(),
      const VectorString& names = VectorString(),
      const VectorString& locatorNames = VectorString(),
      bool flagAddSampleRank = true);
    Id resetFromCSV(
      const String& filename,
      bool verbose,
      const CSVformat& csvfmt,
      Id ncol_max = -1,
      Id nrow_max = -1,
      bool flagAddSampleRank = true);
    Id resetFromBox(
      Id nech,
      const VectorDouble& coormin,
      const VectorDouble& coormax,
      Id ndim = 2,
      double extend = 0.,
      Id seed = 321415,
      bool flagAddSampleRank = true);
    Id resetFromOnePoint(
      const VectorDouble& tab = VectorDouble(),
      bool flagAddSampleRank = true);
    Id resetSamplingDb(
      const Db* dbin,
      double proportion = 0,
      Id number = 0,
      const VectorString& names = VectorString(),
      Id seed = 23241,
      bool verbose = false,
      bool flagAddSampleRank = true);
    Id resetReduce(
      const Db* dbin,
      const VectorString& names = VectorString(),
      const VectorInt& ranks = VectorInt(),
      bool flagIsotopic = false,
      bool verbose = false);
    Id resetFromGridRandomized(
      const DbGrid* dbin,
      double randperc = 0.,
      bool flagAddSampleRank = true);
    /**@}*/

    /** @addtogroup DB_Creators Creating a Db structure
     * \ingroup DB
     *
     * All methods enabling to Create a new Db in various conditions.
     *
     * They all return a pointer to the newly created Db structure.
     *  @{
     */
    static Db* create();
    static Db* createFromNF(const String& NFFilename, bool verbose = true);
    static Db* createFromSamples(
      Id nech,
      const ELoadBy& order = ELoadBy::fromKey("SAMPLE"),
      const VectorDouble& tab = VectorDouble(),
      const VectorString& names = VectorString(),
      const VectorString& locatorNames = VectorString(),
      bool flagAddSampleRank = true);
    static Db* createFromCSV(
      const String& filename,
      const CSVformat& csv = CSVformat(),
      bool verbose = false,
      Id ncol_max = -1,
      Id nrow_max = -1,
      bool flagAddSampleRank = true);
    static Db* createFromBox(
      Id nech,
      const VectorDouble& coormin,
      const VectorDouble& coormax,
      Id seed = 43241,
      bool flag_exact = true,
      bool flag_repulsion = false,
      double range = 0.,
      double beta = 0.,
      double extend = 0.,
      bool flagAddSampleRank = true);
    static Db* createFromOnePoint(
      const VectorDouble& tab = VectorDouble(),
      bool flagAddSampleRank = true);
    static Db* createSamplingDb(
      const Db* dbin,
      double proportion = 0.,
      Id number = 0,
      const VectorString& names = VectorString(),
      Id seed = 23241,
      bool verbose = false,
      bool flagAddSampleRank = true);
    static Db* createFromDbGrid(
      Id nech,
      DbGrid* dbgrid,
      Id seed = 432423,
      bool flag_exact = true,
      bool flag_repulsion = false,
      double range = 0.,
      double beta = 0.,
      bool flagAddSampleRank = true);
    static Db* createReduce(
      const Db* dbin,
      const VectorString& names = VectorString(),
      const VectorInt& ranks = VectorInt(),
      bool flagIsotopic = false,
      bool verbose = false);
    static Db* createFillRandom(
      Id ndat = 10,
      Id ndim = 2,
      Id nvar = 1,
      Id nfex = 0,
      Id ncode = 0,
      double varmax = 0.,
      double selRatio = 0.,
      const VectorDouble& heteroRatio = VectorDouble(),
      const VectorDouble& coormin = VectorDouble(),
      const VectorDouble& coormax = VectorDouble(),
      Id seed = 124234,
      bool flagAddSampleRank = true);
    static Db* createEmpty(
      Id ndat,
      Id ndim = 2,
      Id nvar = 1,
      Id nfex = 0,
      Id ncode = 0,
      bool flagVerr = false,
      bool flagSel = false,
      bool flagAddSampleRank = true);
    static Db* createFromGridRandomized(
      DbGrid* dbgrid,
      double randperc = 0.,
      bool flagAddSampleRank = true);

    /**@}*/

    const VectorDouble& getArrays() const { return _array; }

    /** @addtogroup DB_Names Manipulating Names of the variables contained in a Db
     * \ingroup DB
     *
     * All methods used to manipulated Names of one or several Variables
     * contained in a Db.
     *  @{
     */
    String getNameByLocator(const ELoc& locatorType, Id locatorIndex = 0) const;
    String getNameByColIdx(Id icol) const;
    String getNameByUID(Id iuid) const;

    virtual void initThread() const {}

    VectorString getName(const String& name) const;
    VectorString getNames(const VectorString& names) const;
    VectorString getNamesByLocator(const ELoc& locatorType) const;
    VectorString getNamesByColIdx(const VectorInt& icols) const;
    VectorString getNamesByUID(const VectorInt& iuids) const;
    VectorString getAllNames(
      bool excludeRankAndCoordinates = false,
      bool verbose = false) const;

    void setName(const String& old_name, const String& name);
    void setName(const VectorString& list, const String& name);
    void setNameByUID(Id iuid, const String& name);
    void setNameByColIdx(Id icol, const String& name);
    void setNameByLocator(const ELoc& locatorType, const String& name);

    VectorString expandNameList(const VectorString& names) const;
    VectorString expandNameList(const String& names) const;
    VectorString identifyNames(const VectorString& names) const;

    /**@}*/

    inline Id getNUIDMax() const { return static_cast<Id>(_uidcol.size()); }

    inline Id getNColumn() const { return _ncol; }

    static Id getNEloc();
    Id getNSample(bool useSel = false) const;
    Id getNSampleActiveAndDefined(Id item) const;
    Id getNSampleActiveAndDefined(const String& name) const;
    Id getNSampleActive() const;
    Id getRankRelativeToAbsolute(Id irel) const;
    Id getRankAbsoluteToRelative(Id iabs) const;
    VectorInt getRankRelativeToAbsoluteVec() const;
    VectorInt getRankAbsoluteToRelativeVec() const;

    void clearLocators(const ELoc& locatorType);

    void clearSelection() { clearLocators(ELoc::SEL); }

    void setLocatorByUID(
      Id iuid,
      const ELoc& locatorType = ELoc::fromKey("UNDEFINED"),
      Id locatorIndex = 0,
      bool cleanSameLocator = false);
    void setLocatorByColIdx(
      Id icol,
      const ELoc& locatorType = ELoc::fromKey("UNDEFINED"),
      Id locatorIndex = 0,
      bool cleanSameLocator = false);
    void setLocator(
      const String& name,
      const ELoc& locatorType = ELoc::fromKey("UNDEFINED"),
      Id locatorIndex = 0,
      bool cleanSameLocator = false);
    void setLocators(
      const VectorString& names,
      const ELoc& locatorType = ELoc::fromKey("UNDEFINED"),
      Id locatorIndex = 0,
      bool cleanSameLocator = false);
    void setLocatorsByUID(
      Id number,
      Id iuid,
      const ELoc& locatorType = ELoc::fromKey("UNDEFINED"),
      Id locatorIndex = 0,
      bool cleanSameLocator = false);
    void setLocatorsByUID(
      const VectorInt& iuids,
      const ELoc& locatorType = ELoc::fromKey("UNDEFINED"),
      Id locatorIndex = 0,
      bool cleanSameLocator = false);
    void setLocatorsByColIdx(
      const VectorInt& icols,
      const ELoc& locatorType = ELoc::fromKey("UNDEFINED"),
      Id locatorIndex = 0,
      bool cleanSameLocator = false);
    void addColumnsByVVD(
      const VectorVectorDouble& tab,
      const String& radix,
      const ELoc& locatorType,
      Id locatorIndex = 0,
      bool useSel = false);
    Id addColumns(
      const VectorDouble& tab,
      const String& radix = "New",
      const ELoc& locatorType = ELoc::fromKey("UNDEFINED"),
      Id locatorIndex = 0,
      bool useSel = false,
      double valinit = 0.,
      Id nvar = 1);
    Id addColumnsByConstant(
      Id nadd = 1,
      double valinit = 0.,
      const String& radix = "New",
      const ELoc& locatorType = ELoc::fromKey("UNDEFINED"),
      Id locatorIndex = 0,
      Id nechInit = 0);
    Id addColumnsRandom(
      Id nadd,
      const String& radix = "New",
      const ELoc& locatorType = ELoc::fromKey("Z"),
      Id locatorIndex = 0,
      Id seed = 1352,
      Id nechInit = 0);

    Id addSelection(
      const VectorDouble& tab = VectorDouble(),
      const String& name = "NewSel",
      const String& combine = "set");
    Id addSelectionByVariable(
      const String& varname,
      double lower = TEST,
      double upper = TEST,
      const String& name = "NewSel",
      const String& combine = "set");
    Id addSelectionByRanks(
      const VectorInt& ranks,
      const String& name = "NewSel",
      const String& combine = "set");
    Id addSelectionByLimit(
      const String& testvar,
      const Limits& limits = Limits(),
      const String& name = "NewSel",
      const String& combine = "set");
    Id addSelectionFromDbByConvexHull(
      Db* db,
      double dilate = 0.,
      bool verbose = false,
      const NamingConvention& namconv =
        NamingConvention("Hull", true, true, true, ELoc::fromKey("SEL")));
    Id addSelectionRandom(
      double prop,
      Id seed = 138213,
      const String& name = "NewSel",
      const String& combine = "set");

    Id addSamples(Id nadd, double valinit = TEST);
    Id deleteSample(Id e_del);
    Id deleteSamples(const VectorInt& e_dels);
    void switchLocator(const ELoc& locatorType_in, const ELoc& locatorType_out);
    Id getLastUID(Id number = 0) const;
    String getLastName(Id number = 0) const;

    Id getColIdx(const String& name) const;
    Id getColIdxByUID(Id iuid) const;
    Id getColIdxByLocator(const ELoc& locatorType, Id locatorIndex = 0) const;
    VectorInt getColIdxs(const String& name) const;
    VectorInt getColIdxs(const VectorString& names) const;
    VectorInt getColIdxsByUID(const VectorInt& iuids) const;
    VectorInt getColIdxsByLocator(const ELoc& locatorType) const;

    void setColumn(
      const VectorDouble& tab,
      const String& name,
      const ELoc& locatorType = ELoc::fromKey("UNDEFINED"),
      Id locatorIndex = 0,
      bool useSel = false);
    void
      setColumnByUIDOldStyle(const double* tab, Id iuid, bool useSel = false);
    void setColumnByUID(const VectorDouble& tab, Id iuid, bool useSel = false);
    void
      setColumnByColIdx(const VectorDouble& tab, Id icol, bool useSel = false);
    void setColumnsByColIdx(
      const VectorDouble& tabs,
      const VectorInt& icols,
      bool useSel = false);
    void setColumnByColIdxOldStyle(
      const double* tab,
      Id icol,
      bool useSel = false);
    void duplicateColumnByUID(Id iuid_in, Id iuid_out);

    const double* getColumnPtr(const ELoc& locatorType, Id locatorIndex = 0);
    VectorVectorDouble getItem(
      const VectorInt& rows,
      const VectorString& colnames,
      bool useSel = false) const;
    VectorVectorDouble
      getItem(const VectorInt& rows, const String& colname, bool useSel = false)
        const;
    VectorVectorDouble getItem(
      const VectorInt& rows,
      const ELoc& locatorType,
      bool useSel = false) const;
    VectorVectorDouble
      getItem(const VectorString& colnames, bool useSel = false) const;
    VectorVectorDouble
      getItem(const String& colname, bool useSel = false) const;
    VectorVectorDouble
      getItem(const ELoc& locatorType, bool useSel = false) const;
    VectorString getItemNames(const VectorString& colnames) const;
    VectorString getItemNames(const String& colname) const;
    VectorString getItemNames(const ELoc& locatorType) const;

    Id setItem(
      const VectorInt& rows,
      const VectorString& colnames,
      const VectorVectorDouble& values,
      bool useSel = false);
    Id setItem(
      const VectorInt& rows,
      const ELoc& locatorType,
      const VectorVectorDouble& values,
      bool useSel = false);
    Id setItem(
      const VectorString& colnames,
      const VectorVectorDouble& values,
      bool useSel = false);
    Id setItem(
      const ELoc& locatorType,
      const VectorVectorDouble& values,
      bool useSel = false);
    Id setItem(
      const VectorInt& rows,
      const String& colname,
      const VectorDouble& values,
      bool useSel = false);
    Id setItem(
      const String& colname,
      const VectorDouble& values,
      bool useSel = false);

    bool getLocator(
      const String& name,
      ELoc* ret_locatorType,
      Id* ret_locatorIndex) const;
    bool
      getLocatorByColIdx(Id icol, ELoc* ret_locatorType, Id* ret_locatorIndex)
        const;
    bool getLocatorByUID(Id iuid, ELoc* ret_locatorType, Id* ret_locatorIndex)
      const;
    VectorString getLocators(
      bool anyLocator = true,
      const ELoc& locatorType = ELoc::fromKey("UNDEFINED")) const;
    bool isUIDDefined(Id iuid) const;

    Id getUID(const String& name) const;
    Id getUIDByColIdx(Id icol) const;
    Id getUIDByLocator(const ELoc& locatorType, Id locatorIndex = 0) const;

    VectorInt getUIDs(const VectorString& names) const;
    VectorInt getUIDsByLocator(const ELoc& locatorType) const;
    VectorInt getUIDsByColIdx(const VectorInt& icols) const;
    VectorInt getAllUIDs() const;
    void getAllUIDs(VectorInt& iuids) const;

    void copyByUID(Id iuidIn, Id iuidOut);
    void copyByCol(Id icolIn, Id icolOut);

    Id getNFacies(void) const;
    bool hasLocatorDefined(
      const String& name,
      const ELoc& locatorType,
      Id locatorIndex = 0) const;

    // Accessing elements of the contents

    VectorDouble getSampleCoordinates(Id iech) const;
    void getSampleAsSPInPlace(SpacePoint& P, Id iabs) const;
    virtual void getSampleAsSTInPlace(Id iech, SpaceTarget& P) const;
    VectorDouble getSampleLocators(const ELoc& locatorType, Id iech) const;
    VectorVectorDouble
      getIncrements(const VectorInt& iechs, const VectorInt& jechs) const;
    VectorDouble
      getSamplesOneCoordinate(const VectorInt& iechs, Id idim = 0) const;
    Id getSampleClosestTo(const VectorDouble& coor, bool useSel = false) const;

    VectorDouble
      getOneCoordinate(Id idim, bool useSel = false, bool flag_rotate = true)
        const;
    VectorVectorDouble getAllCoordinates(bool useSel = false) const;
    MatrixDense
      getAllCoordinatesMat(const MatrixDense& box = MatrixDense()) const;
    void setCoordinate(Id iech, Id idim, double value);
    void setCoordinates(Id idim, const VectorDouble& coor, bool useSel = false);
    void setSampleCoordinates(Id iech, const VectorDouble& coor);

    double
      getDistance1D(Id iech, Id jech, Id idim = 0, bool flagAbs = false) const;
    double getDistance(Id iech, Id jech) const;
    Id getDistanceVecInPlace(
      Id iech,
      Id jech,
      VectorDouble& dd,
      const Db* db2 = nullptr) const;

    double getValue(const String& name, Id iech) const;
    void setValue(const String& name, Id iech, double value);

    double getArray(Id iech, Id iuid) const;
    void
      getArrayVec(const VectorInt& iechs, Id iuid, VectorDouble& values) const;
    void setArray(Id iech, Id iuid, double value);
    void
      setArrayVec(const VectorInt& iechs, Id iuid, const VectorDouble& values);
    void updArray(Id iech, Id iuid, const EOperator& oper, double value);
    void updArrayVec(
      const VectorInt& iechs,
      Id iuid,
      const EOperator& oper,
      VectorDouble& values);
    VectorDouble getArrayByUID(Id iuid, bool useSel = false) const;
    void setArrayByUID(const VectorDouble& tab, Id iuid);
    void getArrayBySample(VectorDouble& vals, Id iech) const;
    void setArrayBySample(Id iech, const VectorDouble& vec);

    void getSamplesAsSP(
      std::vector<SpacePoint>& pvec,
      const ASpaceSharedPtr& space,
      bool useSel = false) const;
    void getSamplesFromNbghAsSP(
      std::vector<SpacePoint>& pvec,
      const ASpaceSharedPtr& space,
      const VectorInt& nbgh) const;

    bool hasLocator(const ELoc& locatorType) const;
    double getFromLocator(const ELoc& locatorType, Id iech, Id locatorIndex = 0)
      const;
    void setFromLocator(
      const ELoc& locatorType,
      Id iech,
      Id locatorIndex,
      double value);

    double getValueByColIdx(Id iech, Id icol, bool flagCheck = true) const;
    const double* getColAdressByColIdx(Id icol) const;

    void
      setValueByColIdx(Id iech, Id icol, double value, bool flagCheck = true);
    VectorDouble getValuesByNames(
      const VectorInt& iechs,
      const VectorString& names,
      bool bySample = false) const;
    VectorDouble getValuesByColIdx(
      const VectorInt& iechs,
      const VectorInt& icols,
      bool bySample = false) const;
    void setValuesByNames(
      const VectorInt& iechs,
      const VectorString& names,
      const VectorDouble& values,
      bool bySample = false);
    void setValuesByColIdx(
      const VectorInt& iechs,
      const VectorInt& icols,
      const VectorDouble& values,
      bool bySample = false);

    Table getStatsAsTable(
      const VectorString& names = VectorString(),
      const std::vector<EStatOption>& opers = EStatOption::fromKeys(
        {"NUM", "MINI", "MAXI", "MEAN", "STDV", "VAR"})) const;
    Table getCorrelationAsTable(const VectorString& names) const;
    Table getStatsByCategoryAsTable(
      const String& name,
      const String& category,
      const std::vector<EStatOption>& opers =
        EStatOption::fromKeys({"NUM", "MINI", "MAXI", "MEAN", "STDV", "VAR"}),
      double eps = EPSILON6);
    Table getContentsAsTable(
      const VectorString& names = VectorString(),
      bool useSel = false) const;

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
    Id getNLoc(const ELoc& loctype) const;
    bool hasLocVariable(const ELoc& loctype) const;
    double getLocVariable(const ELoc& loctype, Id iech, Id item) const;
    void setLocVariable(const ELoc& loctype, Id iech, Id item, double value);
    void updLocVariable(
      const ELoc& loctype,
      Id iech,
      Id item,
      const EOperator& oper,
      double value);
    /**@}*/

    Id getNZValues() const;
    bool hasZVariable() const;
    double getZVariable(Id iech, Id item) const;
    void setZVariable(Id iech, Id item, double value);
    void updZVariable(Id iech, Id item, const EOperator& oper, double value);

    VectorDouble
      getLocVariables(const ELoc& loctype, Id iech, Id nitemax = 0) const;
    void
      setLocVariables(const ELoc& loctype, Id iech, const VectorDouble& values);

    bool isNVarComparedTo(Id nvar, Id compare = 0) const;
    bool isIsotopic(Id iech, Id nvar_max = -1) const;
    bool isAllUndefined(Id iech) const;
    bool isAllUndefinedByType(const ELoc& loctype, Id iech) const;
    bool isAllIsotopic() const;

    void setInterval(Id iech, Id item, double rklow = TEST, double rkup = TEST);
    Id getNInterval() const;
    Id getNBounds() const;
    void setBound(Id iech, Id item, double lower = TEST, double upper = TEST);
    VectorDouble getWithinBounds(Id item, bool useSel = false) const;
    VectorDouble getGradient(Id item, bool useSel = false) const;
    VectorDouble getTangent(Id item, bool useSel = false) const;
    VectorDouble getCodeList(void) const;

    Id getSelection(Id iech) const;
    VectorDouble getSelections(void) const;
    void getSampleRanksPerVariable(
      VectorInt& ranks,
      const VectorInt& nbgh = VectorInt(),
      Id ivar = -1,
      bool useSel = true,
      bool useZ = true,
      bool useVerr = false,
      bool useExtD = true) const;
    VectorVectorInt getSampleRanks(
      const VectorInt& ivars = VectorInt(),
      const VectorInt& nbgh = VectorInt(),
      bool useSel = true,
      bool useZ = true,
      bool useVerr = false,
      bool useExtD = true) const;
    void getSampleRanksInPlace(
      VectorVectorInt& sampleRanks,
      const VectorInt& ivars = VectorInt(),
      const VectorInt& nbgh = VectorInt(),
      bool useSel = true,
      bool useZ = true,
      bool useVerr = false,
      bool useExtD = true) const;
    VectorDouble getValuesByRanks(
      const VectorVectorInt& sampleRanks,
      const VectorDouble& means = VectorDouble(),
      bool subtractMean = true) const;
    void getValuesByRanksInPlace(
      VectorDouble* values,
      const VectorVectorInt& sampleRanks,
      const VectorDouble& means = VectorDouble(),
      bool subtractMean = true) const;
    static VectorInt getMultipleSelectedRanks(
      const VectorVectorInt& index,
      const VectorInt& ivars = VectorInt(),
      const VectorInt& nbgh = VectorInt());
    static VectorInt getMultipleSelectedVariables(
      const VectorVectorInt& index,
      const VectorInt& ivars = VectorInt(),
      const VectorInt& nbgh = VectorInt());

    Id getListOfSampleIndicesInPlace(
      Id nvar,
      VectorInt& cumul,
      VectorVectorInt& ranks,
      bool useSel = true) const;
    double getWeight(Id iech) const;
    VectorDouble getWeights(bool useSel = false) const;

    /** @addtogroup DB_1 Designating Variables (used for simulations in particular)
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
    static Id getSimRank(Id isimu, Id ivar, Id icase, Id nbsimu, Id nvar);
    double getSimvar(
      const ELoc& locatorType,
      Id iech,
      Id isimu,
      Id ivar,
      Id icase,
      Id nbsimu,
      Id nvar) const;
    void setSimvar(
      const ELoc& locatorType,
      Id iech,
      Id isimu,
      Id ivar,
      Id icase,
      Id nbsimu,
      Id nvar,
      double value);
    void updSimvar(
      const ELoc& locatorType,
      Id iech,
      Id isimu,
      Id ivar,
      Id icase,
      Id nbsimu,
      Id nvar,
      const EOperator& oper,
      double value);
    /**@}*/

    bool isActive(Id iech) const;
    bool isActiveDomain(Id iech) const;
    bool isActiveAndDefined(Id iech, Id item) const;
    VectorBool getActiveArray() const;

    VectorInt getSortArray() const;
    double
      getCosineToDirection(Id iech1, Id iech2, const VectorDouble& codir) const;

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
    VectorDouble getColumn(
      const String& name,
      bool useSel = false,
      bool flagCompress = true) const;
    VectorDouble
      getColumnByUID(Id iuid, bool useSel = false, bool flagCompress = true)
        const;
    VectorDouble getColumnByLocator(
      const ELoc& locatorType,
      Id locatorIndex = 0,
      bool useSel = false,
      bool flagCompress = true) const;
    VectorDouble
      getColumnByColIdx(Id icol, bool useSel = false, bool flagCompress = true)
        const;

    VectorDouble
      getAllColumns(bool useSel = false, bool flagCompress = true) const;
    VectorDouble getColumns(
      const VectorString& names = VectorString(),
      bool useSel = false,
      bool flagCompress = true,
      const VectorDouble& origins = VectorDouble()) const;
    VectorVectorDouble getColumnsAsVVD(
      const VectorString& names = VectorString(),
      bool useSel = false,
      bool flagCompress = true) const;
    MatrixDense getColumnsAsMatrix(
      const VectorString& names,
      bool useSel = false,
      bool flagCompress = true) const;
    VectorDouble getColumnsByColIdx(
      const VectorInt& icols = VectorInt(),
      bool useSel = false,
      bool flagCompress = true,
      const VectorDouble& origins = VectorDouble()) const;
    VectorDouble getColumnsByColIdxInterval(
      Id icol_beg,
      Id icol_end,
      bool useSel = false,
      bool flagCompress = true) const;
    VectorDouble getColumnsActiveAndDefined(
      const ELoc& locatorType,
      const VectorDouble& origins = VectorDouble()) const;
    VectorDouble getColumnsByLocator(
      const ELoc& locatorType,
      bool useSel = false,
      bool flagCompress = true,
      const VectorDouble& origins = VectorDouble()) const;
    VectorDouble getColumnsByUID(
      const VectorInt& iuids,
      bool useSel = false,
      bool flagCompress = true,
      const VectorDouble& origins = VectorDouble()) const;
    VectorDouble getColumnsByUIDInterval(
      Id iuid_beg,
      Id iuid_end,
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
    void deleteColumnByUID(Id iuid_del);
    void deleteColumnByColIdx(Id icol_del);

    void deleteColumns(const VectorString& names);
    void deleteColumnsByLocator(const ELoc& locatorType);
    void deleteColumnsByUID(const VectorInt& iuids);
    void deleteColumnsByColIdx(const VectorInt& icols);
    void deleteColumnsByUIDRange(Id i_del, Id n_del);
    /**@}*/

    /** @addtogroup DB_4 Calculating Spatial characteristics on the Db
     * \ingroup DB
     *
     * @param idim Rank of the target space dimension (0 based)
     * @param useSel When TRUE, the characteristics are derived from the only
     * active samples
     * @param mini Vector of minimum values (modified by this function)
     * @param maxi Vector of maximum values (modified by this function)
     *
     *  @{
     */
    VectorDouble getExtrema(Id idim, bool useSel = false) const;
    VectorVectorDouble getExtremas(bool useSel = false) const;
    VectorDouble getExtends(bool useSel = false) const;
    VectorDouble getCoorMinimum(bool useSel = false) const;
    VectorDouble getCoorMaximum(bool useSel = false) const;
    double getExtension(Id idim, bool useSel = false) const;
    double getExtensionDiagonal(bool useSel = false) const;
    double getCenter(Id idim, bool useSel = false) const;
    VectorDouble getCenters(bool useSel = false) const;
    void getExtensionInPlace(
      VectorDouble& mini,
      VectorDouble& maxi,
      bool flagPreserve = false,
      bool useSel = false) const;
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
    double getCorrelation(
      const String& name1,
      const String& name2,
      bool useSel = false) const;
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
    bool isColIdxValid(Id icol) const;
    bool isUIDValid(Id iuid) const;
    bool isSampleIndexValid(Id iech) const;
    bool
      isSampleIndicesValid(const VectorInt& iechs, bool useSel = false) const;
    bool isLocatorIndexValid(const ELoc& locatorType, Id locatorIndex) const;
    bool isDimensionIndexValid(Id idim) const;
    /**@}*/

    void
      _combineSelection(VectorDouble& sel, const String& combine = "set") const;

    void generateRank(const String& radix = "rank");

    VectorInt shrinkToValidRows(const VectorInt& rows) const;
    VectorInt shrinkToValidCols(const VectorInt& cols) const;

    static const Db*
      coverSeveralDbs(const Db* db1, const Db* db2, bool* isBuilt);

    /** @addtogroup DB_7 Calculating several statistics in Db
     * \ingroup DB
     *
     * These functions are meant to calculate several statistics on a set of target variables per sample.
     * The resulting values are stored in variables newly created in the same Db.
     *
     * @param names Vector of target variable names
     * @param opers Vector of operations to be performed
     * @param flagIso The statistics are calculated only for samples where all target variables have defined values
     * @param proba              For 'quant': the quantile for this probability is calculated
     * @param vmin               For 'prop', 'T', 'Q', 'M', 'B': defines the lower bound of the interval to work in
     * @param vmax               For 'prop', 'T', 'Q', 'M', 'B': defines the upper bound of the interval to work in
     * @param namconv            Naming Convention used as a radix for the variables newly created in the Db
     * (only used when 'flagStoreInDb' is TRUE)
     *
     * @return If there is more than one operator and more than one variable, the statistics are ordered first by variables
     * (all the statistics of the first variable, then all the statistics of the second variable...).
     *
     *  @{
     */
    void statisticsBySample(
      const VectorString& names,
      const std::vector<EStatOption>& opers = EStatOption::fromKeys({"MEAN"}),
      bool flagIso = true,
      double proba = TEST,
      double vmin = TEST,
      double vmax = TEST,
      const NamingConvention& namconv = NamingConvention("Stats"));
    /**@}*/

    /** @addtogroup DB_8 Calculating correlations on variables of a Db
     * \ingroup DB
     *
     * These functions calculate the correlation matrix based on a set of variables contained in a Db.
     * Although the result stands as a matrix, they are returned as a Vector.
     *
     * @param names Vector of target variable names
     * @param flagIso The statistics are calculated only for samples where all target variables have defined values
     * @param verbose Verbose flag
     * @param title If verbose, the title of the printed statistics.
     *
     * @return These functions return a vector containing the correlation matrix.
     *  @{
     */
    VectorDouble statisticsMulti(
      const VectorString& names,
      bool flagIso = true,
      bool verbose = false,
      const String& title = "") const;
    /**@}*/

    bool areSame(
      const String& name1,
      const String& name2,
      double eps = EPSILON3,
      bool useSel = true,
      bool verbose = false) const;

    VectorInt filter(
      const String& name,
      const Interval& interval,
      Id belowRow = ITEST,
      Id aboveRow = ITEST) const;

    Table printOneSample(
      Id iech,
      const VectorString& names = VectorString(),
      bool excludeCoordinates = true,
      bool skipTitle = false) const;

    void dumpGeometry(Id iech, Id jech) const;

    // Operator overload
    double& operator()(Id iech, const String& name)
    {
      static double dummy = std::numeric_limits<double>::quiet_NaN();
      auto iuid = getUID(name);
      if (!isUIDValid(iuid)) return dummy;
      auto icol = getColIdxByUID(iuid);
      if (!isColIdxValid(icol)) return dummy;
      if (!isSampleIndexValid(iech)) return dummy;
      auto iad = _getAddress(iech, icol);
      return _array[iad];
    }

    const double& operator()(Id iech, const String& name) const
    {
      static const double dummy = std::numeric_limits<double>::quiet_NaN();
      auto iuid = getUID(name);
      if (!isUIDValid(iuid)) return dummy;
      auto icol = getColIdxByUID(iuid);
      if (!isColIdxValid(icol)) return dummy;
      if (!isSampleIndexValid(iech)) return dummy;
      auto iad = _getAddress(iech, icol);
      return _array[iad];
    }

  protected:
    bool _deserializeAscii(std::istream& is) override;
    bool _serializeAscii(std::ostream& os) const override;

    void _clear();
    void _createRank(Id icol = 0);
    void _addRank(Id nech);
    void _loadData(
      const VectorDouble& tab,
      const VectorString& names,
      const VectorString& locatorNames,
      const ELoadBy& order,
      Id shift);
    void _loadData(
      const ELoadBy& order,
      bool flagAddSampleRank,
      const VectorDouble& tab);
    void _defineDefaultNames(Id shift, const VectorString& names);
    void _defineDefaultLocators(Id shift, const VectorString& locatorNames);
    void _setNameByColIdx(Id icol, const String& name);
    String _toStringCommon(const AStringFormat* strfmt) const;
    String _summaryString(void) const;

  private:
    Id _getNextLocator(const ELoc& locatorType) const;

    const VectorInt& _getUIDcol() const { return _uidcol; }

    VectorString _getNames() const { return _colNames; }

    Id _getUIDcol(Id iuid) const;
    Id _getAddress(Id iech, Id icol) const;
    void _columnInit(
      Id ncol,
      Id icol0,
      bool flagCst = true,
      double valinit = TEST);
    String _summaryVariables(void) const;
    String _summaryExtensions(void) const;
    String _summaryStats(VectorInt cols, Id mode = 1, Id maxNClass = 50) const;
    String _summaryLocators(void) const;
    String _summaryUIDs(void) const;
    String _summaryArrays(VectorInt cols, bool useSel = true) const;

    void _defineDefaultLocatorsByNames(Id shift, const VectorString& names);
    VectorInt _getUIDsBasic(const VectorString& names) const;

    Id _getLastColumn(Id number = 0) const;

    Id _findColumnInLocator(const ELoc& locatorType, Id icol) const;
    Id _findUIDInLocator(const ELoc& locatorType, Id iuid) const;
    String _getLocatorNameByColIdx(Id icol) const;
    VectorInt _ids(const String& name, bool flagOne, bool verbose = true) const;
    VectorInt
      _ids(const VectorString& names, bool flagOne, bool verbose = true) const;
    VectorInt
      _ids(const ELoc& locatorType, bool flagOne, bool verbose = true) const;
    VectorInt
      _ids(const VectorInt& iuids, bool flagOne, bool verbose = true) const;

    VectorDouble
      _getItem(const String& exp_name, bool useSel, const VectorInt& rows)
        const;
    void _setItem(
      const String& name,
      const VectorInt& rows,
      const VectorDouble& values);
    void _setItem(const String& name, bool useSel, const VectorDouble& values);
    bool _isValidCountRows(
      const VectorInt& rows,
      bool useSel,
      const VectorDouble& values) const;
    bool _isValidCountRows(bool useSel, const VectorDouble& values) const;
    VectorString
      _getVarNames(const VectorString& colnames, Id expectedVarCount);
    Id _getListOfSampleIndicesPerVariableInPlace(
      VectorInt& ranks,
      Id ivar = 0,
      bool useSel = true) const;

    // Higher level methods
    bool
      _isCountValid(const VectorInt& iuids, bool flagOne, bool verbose = true)
        const;

  protected:
    void _defineVariableAndLocators(
      const Db* dbin,
      const VectorString& names,
      Id shift = 0);
    void _loadValues(
      const Db* db,
      const VectorString& names,
      const VectorInt& ranks,
      Id shift = 0);

  private:
    Id _ncol; //!< Number of Columns of data
    Id _nech; //!< Number of samples
    VectorDouble _array; //!< Array of values
    VectorInt _uidcol; //!< UID to Column
    VectorString _colNames; //!< Names of the variables
    std::vector<PtrGeos> _p; //!< Locator characteristics

    /// factor allocations
    mutable VectorInt _uids;
  };

  GSTLEARN_EXPORT bool haveSameNDim(
    const Db* db1,
    const Db* db2,
    const ModelGeneric* model,
    Id* ndim);
  GSTLEARN_EXPORT bool haveCompatibleNVar(
    const Db* db1,
    const Db* db2,
    const ModelGeneric* model,
    Id* nvar);

} // namespace gstlrn
