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
#include "Basic/AStringFormat.hpp"
#include "Basic/Vector.hpp"

#include "geoslib_define.h"

// Do not convert to AEnum (mask combination is not used as enum)
//typedef enum
//{
//  FLAG_RESUME = 1,    //!< Print the Db summary
//  FLAG_VARS = 2,      //!< Print the Field names
//  FLAG_EXTEND = 4,    //!< Print the Db extension
//  FLAG_STATS = 8,     //!< Print the variable statistics
//  FLAG_ARRAY = 16,    //!< Print the variable contents
//  FLAG_LOCATOR = 32,  //!< Print the locators
//} DISPLAY_PARAMS;
//
 static unsigned char FLAG_RESUME = 1;  //!< Print the Db summary
 static unsigned char FLAG_VARS = 2;   //!< Print the Field names
 static unsigned char FLAG_EXTEND = 4; //!< Print the Db extension
 static unsigned char FLAG_STATS = 8; //!< Print the variable statistics
 static unsigned char FLAG_ARRAY = 16; //!< Print the variable contents
 static unsigned char FLAG_LOCATOR = 32;  //!< Print the locators

class GSTLEARN_EXPORT DbStringFormat: public AStringFormat
{
public:
  DbStringFormat(unsigned char params = FLAG_RESUME | FLAG_VARS,
                 const VectorString& names = VectorString(),
                 const VectorInt& cols = VectorInt(),
                 bool useSel = true);
  DbStringFormat(const DbStringFormat& r);
  DbStringFormat& operator=(const DbStringFormat& r);
  virtual ~DbStringFormat();

  static DbStringFormat* create(unsigned char params = FLAG_RESUME | FLAG_VARS,
                                const VectorString& names = VectorString(),
                                const VectorInt& cols = VectorInt(),
                                bool useSel = true);

  DbStringFormat* createFromKeys(unsigned char params,
                                 const VectorString &names,
                                 const VectorInt &cols,
                                 bool useSel);

  const VectorInt& getCols() const { return _cols; }
  bool getUseSel() const { return _useSel; }
  int getMode() const { return _mode; }
  unsigned char getParams() const { return _params; }
  const VectorString& getNames() const { return _names; }

  /**
   * Reduce the set of variables for which the print is provided
   * @param cols Vector of Column indices on which Stats or Array is applied (optional)
   * @remark This selection is performed by Column Rank. It invalidates any possible selection already performed.
   */
  void setCols(const VectorInt& cols) { _names.clear(); _cols = cols; }
  /**
   * Using the current Selection or Not
   * @param useSel Take the selection into account when true
   */
  void setUseSel(bool useSel) { _useSel = useSel; }
  /**
   * @param mode Way to consider the variable for Stats (1: Real; 2: Categorical)
   */
  void setMode(int mode) { _mode = mode; }
  /**
   * Reduce the set of variables for which the print is provided
   * @param names Vector of Column names
   * @remark This selection is performed by Variable Name. It invalidates any possible selection already performed.
   */
  void setNames(const VectorString& names) { _names = names; }
  /**
   * Set the String Format parameters
   * @param params Mask defining the printout
   *
   * @remark The Mask is a combination of DISPLAY_PARAMS, i.e.:
   * @remark - FLAG_RESUME: for a Summary of the contents
   * @remark - FLAG_VARS:   for the Field Names and Locators
   * @remark - FLAG_EXTEND: for the area covered by the Db
   * @remark - FLAG_STATS:  for Basic Statistics on the variables
   * @remark - FLAG_ARRAY:  for the extensive printout of the variables
   */
  void setParams(unsigned char params) { _params = params; }
  void setFlags(bool flag_resume = true,
                bool flag_vars = true,
                bool flag_extend = false,
                bool flag_stats = false,
                bool flag_array = false,
                bool flag_locator = false,
                const VectorString& names = VectorString(),
                const VectorInt& cols = VectorInt(),
                bool useSel = true);

  bool matchResume()  const { return _matchFlag(FLAG_RESUME); }
  bool matchVars()    const { return _matchFlag(FLAG_VARS); }
  bool matchExtend()  const { return _matchFlag(FLAG_EXTEND); }
  bool matchStats()   const { return _matchFlag(FLAG_STATS); }
  bool matchArray()   const { return _matchFlag(FLAG_ARRAY); }
  bool matchLocator() const { return _matchFlag(FLAG_LOCATOR); }

private:
  bool _matchFlag(int flag) const;

private:
  unsigned char _params;
  VectorInt _cols;
  VectorString _names;
  bool _useSel;
  int _mode;
};
