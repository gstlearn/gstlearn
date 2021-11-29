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

#include "Basic/Vector.hpp"
#include "Db/PtrGeos.hpp"
#include "Db/ELoadBy.hpp"
#include "Polygon/Polygons.hpp"
#include "Basic/Limits.hpp"
#include "Basic/GridC.hpp"
#include "Basic/Vector.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/CSVformat.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"

// Do not convert to AEnum (mask combination is not used as enum)
typedef enum
{
  FLAG_RESUME = 1,    //!< Print the Db summary
  FLAG_VARS = 2,      //!< Print the Field names
  FLAG_EXTEND = 4,    //!< Print the Db extension
  FLAG_STATS = 8,     //!< Print the variable statistics
  FLAG_ARRAY = 16,    //!< Print the variable contents
} DISPLAY_PARAMS;

class Limits;

/**
 * Class containing the Data Set.
 * It can be organized as a set of Isolated Points or as a regular Grid
 */
class Db: public AStringable, public ASerializable
{
public:
  Db();
  Db(int nech,
     const ELoadBy& order = ELoadBy::SAMPLE,
     const VectorDouble& tab = VectorDouble(),
     const VectorString& names = VectorString(),
     const VectorString& locatorNames = VectorString(),
     int flag_add_rank = 1);
  Db(const VectorInt& nx,
     const VectorDouble& dx = VectorDouble(),
     const VectorDouble& x0 = VectorDouble(),
     const VectorDouble& angles = VectorDouble(),
     const ELoadBy& order = ELoadBy::SAMPLE,
     const VectorDouble& tab = VectorDouble(),
     const VectorString& names = VectorString(),
     const VectorString& locatorNames = VectorString(),
     int flag_add_rank = 1);
  Db(const String& filename,
     bool verbose,
     const CSVformat& csv,
     int ncol_max = -1,
     int nrow_max = -1,
     int flag_add_rank = 1);
  Db(Db* db,
     const VectorInt&    nodes  = VectorInt(),
     const VectorDouble& dcell  = VectorDouble(),
     const VectorDouble& origin = VectorDouble(),
     const VectorDouble& margin = VectorDouble(),
     int flag_add_rank = 1);
  Db(const String& neutralFileName, bool verbose = false);
  Db(Polygons* polygon,
     const VectorInt& nodes,
     const VectorDouble& dcell,
     int flag_add_rank = 1);
  Db(const Db* dbin,
     double proportion,
     const VectorString& names = VectorString(),
     int seed = 23241,
     bool verbose = false);
  Db(int nech,
     const VectorDouble& coormin,
     const VectorDouble& coormax,
     int ndim = 2,
     int seed = 321415,
     int flag_add_rank = 1);
  Db(const VectorDouble& tab, int flag_add_rank = 1);
  Db(const Db& r);
  Db& operator=(const Db& r);
  virtual ~Db();

public:
  virtual String toString(int level = 0) const override;
  int deSerialize(const String& filename, bool verbose = false) override;
  int serialize(const String& filename, bool verbose = false) const override;

  const VectorDouble& getArrays() const { return _array; }
  String getNameByColumn(int icol) const { return _colNames[icol]; }
  String getName(int iatt) const;
  String getName(const ELoc& locatorType, int locatorIndex=0) const;
  VectorString getNames(const VectorString& names) const;
  VectorString getNames(const String& name) const;
  VectorString getNames(const ELoc& locatorType) const;
  VectorString getNames(const VectorInt& iatts) const;
  VectorString getNames() const;
  void setName(const String& old_name, const String& name);
  void setName(const VectorString list, const String& name);
  void setName(int iatt, const String& name);
  void setName(const ELoc& locatorType, const String& name);
  const GridC& getGrid() const { return _grid; }
  int getAttributeMaxNumber() const { return static_cast<int>(_attcol.size()); }
  int getFieldNumber() const { return _ncol; }
  double getFieldSize(bool useSel = false) const;
  int getSampleNumber() const { return _nech; }
  int getActiveSampleNumber() const;
  int isGrid() const { return _isGrid; }

  VectorString expandNameList(const VectorString& names) const;
  VectorString expandNameList(const String& names) const;
  VectorInt ids(const String& name, bool flagOne) const;
  VectorInt ids(const VectorString& names, bool flagOne) const;
  VectorInt ids(const ELoc& locatorType, bool flagOne) const;
  VectorInt ids(const VectorInt& iatts, bool flagOne) const;

  void reset(int ncol, int nech);

  // Locator and Attribute methods

  void clearLocators(const ELoc& locatorType);
  void setLocatorByAttribute(int iatt,
                             const ELoc& locatorType,
                             int locatorIndex = 0);
  void setLocator(const VectorString& names,
                  const ELoc& locatorType = ELoc::UNKNOWN,
                  int locatorIndex = 0);
  void setLocator(const String& names,
                  const ELoc& locatorType = ELoc::UNKNOWN,
                  int locatorIndex = 0);
  void setLocatorsByAttribute(int number,
                              int iatt,
                              const ELoc& locatorType,
                              int locatorIndex = 0);
  int addFields(const VectorDouble& tab,
                const String& radix = "New",
                const ELoc& locatorType = ELoc::UNKNOWN,
                int locatorIndex = 0,
                bool useSel = false,
                double valinit = 0.,
                int nvar = 1);
  int addFields(int nadd,
                double valinit = 0.,
                const String& radix = "New",
                const ELoc& locatorType = ELoc::UNKNOWN,
                int locatorIndex = 0,
                int nechInit = 0);
  int addSelection(const VectorDouble& tab,
                   const String& name = "NewSel");
  int addSelection(const String& testvar,
                   const Limits& limits = Limits(),
                   const String& name = "NewSel");
  int addSamples(int nadd, double valinit);
  void deleteSample(int e_del);
  void switchLocator(const ELoc& locatorTypein, const ELoc& locatorTypeout);
  void printLocators(void) const;
  void printAttributes(void) const;
  int  getLastAttribute(int number = 0) const;
  String getLastName(int number = 0) const;

  int getColumn(const String& name) const;
  int getColumnByAttribute(int iatt) const;
  VectorInt getColumnByAttribute(const VectorInt iatts) const;
  int getColumnByLocator(const ELoc& locatorType, int locatorIndex=0) const;
  VectorDouble getColumnByRank(int icol, bool useSel = false) const;
  void setColumnByRank(const VectorDouble& tab, int icol, bool useSel = false);
  void setColumnByRank(const double* tab, int icol, bool useSel = false);
  void duplicateColumnByAttribute(int iatt_in, int iatt_out);

  VectorInt getColumns(const VectorString& names) const;
  VectorInt getColumnsByAttribute(const ELoc& locatorType) const;
  VectorDouble getColumnsByRank(const VectorInt& icols = VectorInt(),
                                bool useSel = false) const;
  VectorDouble getColumnsByRank(int icol_beg,
                                int icol_end,
                                bool useSel = false) const;

  int getLocatorByColumn(int icol,
                         ELoc* ret_locatorType,
                         int* ret_locatorIndex) const;
  int getLocator(int iatt,
                 ELoc* ret_locatorType,
                 int* ret_locatorIndex) const;
  int getLocator(const String& name,
                 ELoc* ret_locatorType,
                 int* ret_locatorIndex) const;
  VectorString getLocators(bool anyLocator = true,
                           const ELoc& locatorType = ELoc::UNKNOWN) const;
  int getLocatorNumber(const ELoc& locatorType) const;
  bool isAttributeDefined(int iatt) const;

  int getAttribute(const ELoc& locatorType, int locatorIndex=0) const;
  int getAttribute(const String &name) const;
  VectorInt getAttributes(const VectorString& names) const;
  VectorInt getAttributes(const ELoc& locatorType) const;
  VectorInt getAttributes() const;
  VectorInt getAttributesBasic(const VectorString& names) const;
  int getFaciesNumber(void) const;

  // Accessing elements of the contents

  VectorDouble getSampleCoordinates(int iech) const;
  void   getSampleCoordinates(int iech, VectorDouble& coor) const;
  VectorDouble getSampleAttributes(const ELoc& locatorType, int iech) const;

  double getCoordinate(int iech, int idim, bool flag_rotate=true) const;
  void   getCoordinate(int iech, VectorDouble& coor, bool flag_rotate=true) const;
  void   setCoordinate(int iech, int idim, double value);
  VectorDouble getCoordinate(int idim, bool useSel = false, bool flag_rotate = true) const;
  double getDistance1D(int iech, int jech, int idim, bool flagAbs = false) const;
  double getDistance(int iech, int jech) const;

  VectorVectorDouble getCoordinates(bool useSel = false) const;

  double getValue(const String& name, int iech) const;
  void   setValue(const String& name, int iech, double value);

  double getArray(int iech, int iatt) const;
  void   setArray(int iech, int iatt, double value);
  void   updArray(int iech, int iatt, int oper, double value);
  VectorDouble getArray(int iatt, bool useSel = false) const;

  int    getFromLocatorNumber(const ELoc& locatorType) const;
  double getFromLocator(const ELoc& locatorType, int iech, int locatorIndex=0) const;
  void   setFromLocator(const ELoc& locatorType,
                        int iech,
                        int locatorIndex,
                        double value);

  double getByColumn(int iech, int icol) const;
  void   setByColumn(int iech, int icol, double value);

  int    getVariableNumber() const;
  bool   hasVariable() const;
  double getVariable(int iech, int item) const;
  void   setVariable(int iech, int item, double value);
  void   updVariable(int iech, int item, int oper, double value);
  bool   isVariableNumberComparedTo(int nvar, int compare = 0) const;

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
  void   setGradient(int iech, int item, double value);

  int    getTangentNumber() const;
  bool   hasTangent() const;
  double getTangent(int iech, int item) const;
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

  int getVarianceErrorNumber() const;
  bool hasVarianceError() const;
  double getVarianceError(int iech, int item) const;
  void setVarianceError(int iech, int item, double value);

  bool hasDomain() const;
  int getDomain(int iech) const;
  void setDomain(int iech, int value);

  int getDipDirectionNumber() const;
  bool hasDipDirection() const;
  double getDipDirection(int iech) const;
  void setDipDirection(int iech, double value);

  int getDipAngleNumber() const;
  bool hasDipAngle() const;
  double getDipAngle(int iech) const;
  void setDipAngle(int iech, double value);

  int getObjectSizeNumber() const;
  bool hasObjectSize() const;
  double getObjectSize(int iech) const;
  void setObjectSize(int iech, double value);

  int getBorderUpNumber() const;
  bool hasBorderUp() const;
  double getBorderUp(int iech) const;
  void setBorderUp(int iech, double value);

  int getBorderDownNumber() const;
  bool hasBorderDown() const;
  double getBorderDown(int iech) const;
  void setBorderDown(int iech, double value);

  int getDateNumber() const;
  bool hasDate() const;
  double getDate(int iech) const;
  void setDate(int iech, double value);

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
  int getActiveAndDefinedNumber(int item) const;
  int getActiveAndDefinedNumber(const String& name) const;
  VectorInt getSortArray() const;
  double getCosineToDirection(int iech1, int iech2, const VectorDouble& codir) const;

  VectorDouble getField(const String& name, bool useSel = false) const;
  VectorDouble getFieldByAttribute(int iatt, bool useSel = false) const;
  VectorDouble getFieldByLocator(const ELoc& locatorType,
                                 int locatorIndex=0,
                                 bool useSel = false) const;
  void setFieldByAttribute(const double* tab, int iatt, bool useSel = false);
  void setFieldByAttribute(const VectorDouble& tab, int iatt, bool useSel = false);
  void setField(const VectorDouble& tab,
                const String& name,
                bool useSel = false);

  VectorDouble getFields(const VectorString& names = VectorString(),
                         bool useSel = false) const;
  VectorDouble getFieldsByLocator(const ELoc& locatorType,
                                  bool useSel = false) const;
  VectorDouble getFieldsByAttribute(const VectorInt& iatts,
                                    bool useSel = false) const;
  VectorDouble getFieldsByAttribute(int iatt_beg,
                                    int iatt_end,
                                    bool useSel = false) const;

  void deleteField(const String& names);
  void deleteField(const VectorString& names);
  void deleteFieldByAttribute(int iatt_del);
  void deleteFieldByLocator(const ELoc& locatorType);
  void deleteFieldByRank(int rank_del);

  VectorDouble getExtrema(int idim, bool useSel = false) const;
  double getExtension(int idim, bool useSel = false) const;
  double getMinimum(const String& name, bool useSel = false) const;
  double getMaximum(const String& name, bool useSel = false) const;
  double getMean(const String& name, bool useSel = false) const;
  double getVariance(const String& name, bool useSel = false) const;
  double getStdv(const String& name, bool useSel = false) const;

  // Statistics

  VectorDouble statistics(const VectorString& names,
                          const VectorString& opers = { "Mean" },
                          bool flagIso = true,
                          bool flagVariableWise = true,
                          bool flagPrint = true,
                          double vmin = TEST,
                          double vmax = TEST,
                          double proba = TEST,
                          const String& title = "",
                          NamingConvention namconv = NamingConvention("Stats"));

  VectorDouble statisticsMulti(const VectorString& names,
                               bool flagIso = true,
                               bool flagPrint = false,
                               const String& title = "");

  // Pipe to the GridC class

  /**
   * Definition of the grid. Function to access the class Grid
   * @param nx     Array giving the number of nodes
   * @param dx     Array giving the mesh dimensions along each space dimension
   * @param x0     Array giving the origin of the Grid (lower left corner)
   * @param angles Array of rotation angles
   */
  int gridDefine(const VectorInt& nx,
                  const VectorDouble& dx = VectorDouble(),
                  const VectorDouble& x0 = VectorDouble(),
                  const VectorDouble& angles = VectorDouble());
  void gridCopyParams(int mode, const GridC& gridaux);
  bool isSameGrid(const GridC& grid) const;
  bool isSameGridMesh(const Db& dbaux) const;
  bool isSameGridRotation(const Db& dbaux) const;
  bool isGridRotated() const;
  int  getNDim() const; /// TODO : rename to getDimensionNumber etc...
  int  getNX(int idim) const;
  int  getNTotal() const { return _grid.getNTotal(); }
  double getCellSize() const { return _grid.getCellSize(); }
  bool hasSameDimension(const Db* dbaux) const;
  bool hasLargerDimension(const Db* dbaux) const;
  VectorInt getNX() const { return _grid.getNX(); }
  double getDX(int idim) const;
  VectorDouble getDX() const { return _grid.getDX(); }
  double getX0(int idim) const;
  VectorDouble getX0() const { return _grid.getX0(); }
  double getAngles(int idim) const;
  VectorDouble getAngles() const { return _grid.getRotAngles(); }
  VectorDouble getRotMat() const { return _grid.getRotMat(); }
  void setNX(int idim, int value) { _grid.setNX(idim, value); }
  void setX0(int idim, double value) { _grid.setX0(idim, value); }
  void setDX(int idim, double value) { _grid.setDX(idim, value); }
  VectorDouble getGridAxis(int idim) const { return _grid.getAxis(idim); }
  VectorDouble getCoordinateFromCorner(const VectorInt& icorner) const
  {
    return _grid.getCoordinateFromCorner(icorner);
  }
  int coordinateToRank(const VectorDouble& coor, double eps = EPSILON6) const
  {
    return _grid.coordinateToRank(coor,eps);
  }
  void rankToCoordinate(int rank,
                        VectorDouble& coor,
                        const VectorDouble& percent = VectorDouble()) const
  {
    _grid.rankToCoordinate(rank, coor, percent);
  }

  // Functions for checking validity of parameters

  bool isColumnIndexValid(int icol) const;
  bool isAttributeIndexValid(int iatt) const;
  bool isSampleIndexValid(int iech) const;
  bool isLocatorIndexValid(const ELoc& locatorType, int locatorIndex) const;
  bool isDimensionIndexValid(int idim) const;

  using AStringable::display; // https://stackoverflow.com/questions/18515183/c-overloaded-virtual-function-warning-by-clang
  /**
   * Display the contents of the Db
   * @param params  Mask defining the printout (see remarks)
   * @param cols    Vector of Column indices on which Stats or Array is applied (optional)
   * @param flagSel Take the selection into account
   * @param mode    Way to consider the variable for Stats (1: Real; 2: Categorical)
   *
   * @remark The Mask is a combination of DISPLAY_PARAMS, i.e.:
   * @remark - FLAG_RESUME: for a Summary of the contents
   * @remark - FLAG_VARS:   for the Field Names and Locators
   * @remark - FLAG_EXTEND: for the area covered by the Db
   * @remark - FLAG_STATS:  for Basic Statistics on the variables
   * @remark - FLAG_ARRAY:  for the extensive printout of the variables
   * @remark For flags, FLAG_STATS and FLAG_ARRAY, you can use the 'cols' argument
   * @remark to restrain the variables of interest
   */
  void displayMore(unsigned char params,
                   const VectorInt& cols = VectorInt(),
                   bool flagSel = true,
                   int mode = 1) const;
  void displayMore(unsigned char params,
                   const VectorString& names,
                   bool flagSel = true,
                   int mode = 1) const;

private:
  void  _initP();
  const VectorInt& _getAttcol() const { return _attcol; }
  const VectorString _getNames() const { return _colNames; }
  double _getAttcol(int icol) const { return _attcol[icol]; }
  int _getAddress(int iech, int icol) const;
  void _columnInit(int ncol, int icol0, double valinit);
  double _updateValue(int oper, double oldval, double value);
  String _summaryString(void) const;
  String _summaryVariableString(void) const;
  String _summaryExtensionString(void) const;
  String _summaryVariableStat(VectorInt cols,
                              int mode = 1,
                              int maxNClass = 50) const;
  String _summaryArrayString(VectorInt cols, bool flagSel = true) const;
  String _display(unsigned char params,
                  const VectorInt& cols = VectorInt(),
                  bool flagSel = true,
                  int mode = 1) const;
  void _loadData(const VectorDouble& tab,
                 const VectorString& names,
                 const VectorString& locatorNames,
                 const ELoadBy& order,
                 int shift);
  void _createRank(int shift = 0);
  void _createGridCoordinates(int shift);
  void _defineDefaultNames(int shift, const VectorString& names);
  void _defineDefaultLocators(int shift, const VectorString& locatorNames);
  void _defineDefaultLocatorsByNames(int shift, const VectorString& names);
  int  _getSimrank(int isimu, int ivar, int icase, int nbsimu, int nvar) const;

  // Column dependent indexing
  void _setNameByColumn(int icol, const String& name);
  int _getLastColumn(int number = 0) const;
  int _getAttributeByColumn(int icol) const;
  int _findColumnInLocator(const ELoc& locatorType, int icol) const;
  int _findAttributeInLocator(const ELoc& locatorType, int iatt) const;
  String _getLocatorNameByColumn(int icol) const;

  // Higher level methods
  VectorDouble _statistics(const VectorInt& iatts,
                          const VectorString& opers = { "Mean" },
                          bool flagIso = true,
                          bool flagVariableWise = true,
                          bool flagPrint = true,
                          double proba = TEST,
                          double vmin = TEST,
                          double vmax = TEST,
                          const String& title = "",
                          NamingConvention namconv = NamingConvention("Stats"));
  VectorDouble _statisticsMulti(const VectorInt& iatts,
                               bool flagIso = true,
                               bool flagPrint = false,
                               const String& title = "");
  int  _variableWrite(bool flag_grid, bool onlyLocator=false,
                      bool writeCoorForGrid=false) const;
  void _variableRead(int *natt_r,
                     int *ndim_r,
                     int *nech_r,
                     std::vector<ELoc>& tabatt,
                     VectorInt& tabnum,
                     VectorString& tabnam,
                     VectorDouble& tab);
  void _loadData(const ELoadBy& order, int flag_add_rank, const VectorDouble& tab);
  bool _isCountValid(const VectorInt iatts, bool flagOne) const;

private:
  int _ncol;                 //!< Number of Columns of data
  int _nech;                 //!< Number of samples
  int _isGrid;               //!< Is the Data Base organized as a Grid
  VectorDouble _array;       //!< Array of values
  VectorInt _attcol;         //!< Attribute to Column
  VectorString _colNames;    //!< Names of the variables
  std::map<ELoc,PtrGeos> _p; //!< Locator characteristics
  GridC _grid;               //!< Grid characteristics
};
