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
#include "Basic/StringHelper.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/String.hpp"
#include "Enum/ECst.hpp"

#include <cctype>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <locale>
#include <sstream>

#define CASE_DOUBLE 0
#define CASE_REAL 1
#define CASE_INT 2
#define CASE_COL 3
#define CASE_ROW 4

namespace gstlrn
{
  Id _getColumnRank()
  {
    return static_cast<Id>(OptCst::query(ECst::NTRANK));
  }

  Id _getColumnName()
  {
    return static_cast<Id>(OptCst::query(ECst::NTNAME));
  }

  Id _getColumnSize(Id localSize, Id nColumns)
  {
    if (localSize > 0) return localSize;
    Id size = 1 + static_cast<Id>(OptCst::query(ECst::NTCAR));
    if (nColumns == 1) return size;
    return size * nColumns; // account for spaces between columns
  }

  Id _getMaxNCols()
  {
    return static_cast<Id>(OptCst::query(ECst::NTCOL));
  }

  Id _getMaxNRows()
  {
    return static_cast<Id>(OptCst::query(ECst::NTROW));
  }

  Id _getNBatch()
  {
    return static_cast<Id>(OptCst::query(ECst::NTBATCH));
  }

  Id _getDecimalNumber()
  {
    return static_cast<Id>(OptCst::query(ECst::NTDEC));
  }

  double _getThresh()
  {
    Id ndec = static_cast<Id>(OptCst::query(ECst::NTDEC));
    double thresh = (0.5 * pow(10, -ndec));
    return thresh;
  }

  String _toStrRowColumn(Id icase, Id value, Id flagAdd)
  {
    std::stringstream sstr;
    I32 rank = static_cast<I32>(_getColumnRank());
    I32 width = static_cast<I32>(_getColumnSize() - _getColumnRank() - 1);
    sstr << std::setw(width) << std::right;
    if (icase == CASE_ROW)
    {
      if (!flagAdd)
        sstr << "[" << std::setw(rank) << value << ",]";
      else
        sstr << "[" << std::setw(rank) << value << "+]";
    }
    else
    {
      if (!flagAdd)
        sstr << "[," << std::setw(rank) << value << "]";
      else
        sstr << "[ " << std::setw(rank) << value << "]";
    }
    return sstr.str();
  }

  // Functions for beautifying a suite of values
  String _toStrRowHeader(const VectorString& rownames, Id iy, Id rowSize)
  {
    std::stringstream sstr;
    if (!rownames.empty())
      sstr << toStr(rownames[iy], -1, rowSize);
    else
      sstr << _toStrRowColumn(CASE_ROW, iy, false);
    return sstr.str();
  }

  String _toStrColumnHeaders(const VectorString& colnames,
                             Id colfrom,
                             Id colto,
                             Id colSize)
  {
    std::stringstream sstr;
    if (!colnames.empty())
    {
      // By Names
      sstr << toStr(" ", 1) << " ";
      for (Id ix = colfrom; ix < colto; ix++)
        sstr << toStr(colnames[ix], 1, colSize);
      sstr << std::endl;
    }
    else
    {
      // By Numbers
      sstr << toStr(" ", 1) << " ";
      for (Id ix = colfrom; ix < colto; ix++)
        sstr << _toStrRowColumn(CASE_COL, ix, false);
      sstr << std::endl;
    }
    return sstr.str();
  }

  String _toStrTrailer(Id ncols, Id nrows, Id ncols_util, Id nrows_util)
  {
    std::stringstream sstr;

    bool one_used = (ncols != ncols_util || nrows != nrows_util);
    bool all_used = (ncols != ncols_util && nrows != nrows_util);

    if (one_used) sstr << "(";

    if (ncols != ncols_util)
    {
      if (ncols == ncols_util)
        sstr << "Ncols=" << ncols;
      else
        sstr << "Ncols=" << ncols_util << "[from " << ncols << "]";
    }

    if (all_used) sstr << ",";

    if (nrows != nrows_util)
    {
      if (nrows == nrows_util)
        sstr << "Nrows=" << nrows;
      else
        sstr << "Nrows=" << nrows_util << "[from " << nrows << "]";
    }

    if (one_used) sstr << ")" << std::endl;
    return sstr.str();
  }

  /**
   * Decode an double from a string. Returns TEST if impossible
   * @param v String to be decoded
   * @param dec Decimal separator character
   * @return The double value or TEST (in case of failure)
   */
  template<typename T>
  class dec_separator: public std::numpunct<T>
  {
  public:
    dec_separator(char dec = ',')
      : _dec(dec)
    {
    }

  private:
    typename std::numpunct<T>::char_type do_decimal_point() const override
    {
      return _dec;
    }

    char _dec;
  };

  /**
   * Decode an integer from a string. Returns ITEST when impossible
   * @param v String to be decoded
   * @param dec Symbol used for decimal point (unused here)
   * @return The integer value or ITEST (in case of failure)
   */
  Id _convertToId(const String& v, char dec)
  {
    DECLARE_UNUSED(dec);
    std::istringstream iss(v);
    Id number;
    iss >> number;
    if (iss.fail()) return ITEST;
    return number;
  }

  /**
   * Decode a double value from a string. Returns TEST when impossible
   * @param v String to be decoded
   * @param dec Symbol used for decimal point
   * @return The double value or TEST (in case of failure)
   */
  double _convertToDouble(const String& v, char dec)
  {
    std::istringstream iss(v);
    double number;
    iss.imbue(std::locale(iss.getloc(), new dec_separator<char>(dec)));
    iss >> number;
    if (iss.fail()) return TEST;
    return number;
  }

  /**
   * Ask interactively for the value of one integer
   * @param v Text of the question
   * @param defval Default value (or IFFFF)
   * @param authTest True if TEST value is authorized (TEST)
   */
  Id _askInt(const String& v, Id defval, bool authTest)
  {
    bool hasDefault = !isNA(defval) || authTest;
    Id answer = defval;
    std::cin.exceptions(std::istream::failbit | std::istream::badbit);

    try
    {
      while (true)
      {
        // Display the question
        if (hasDefault)
        {
          if (isNA(defval))
            std::cout << v << " (Default = TEST) : ";
          else
            std::cout << v << " (Default = " << defval << ") : ";
        }
        else
          std::cout << v << " : ";

        // Read the answer
        String str;
        std::getline(std::cin, str);

        // Check for empty line: set to default value
        if (str.empty() && hasDefault)
        {
          answer = defval;
          break;
        }

        // Check the TEST answer

        if (authTest && str == "TEST")
        {
          answer = ITEST;
          break;
        }

        // Try casting in integer
        std::stringstream ss(str);
        if (ss >> answer) break;

        std::cout << "The answer is not a valid integer!" << std::endl;
      }
    }
    catch (std::istream::failure& e)
    {
      std::cerr << "Problem when reading integer:" << e.what() << std::endl;
    }
    return answer;
  }

  /**
   * Ask interactively for the value of one Real (Double value)
   * @param v Text of the question
   * @param defval Default value (or IFFFF)
   * @param authTest True if a TEST answer is authorized (TEST)
   */
  double _askDouble(const String& v, double defval, bool authTest)
  {
    bool hasDefault = !FFFF(defval) || authTest;
    double answer = defval;
    std::cin.exceptions(std::istream::failbit | std::istream::badbit);

    try
    {
      while (true)
      {
        // Display the question
        if (hasDefault)
        {
          if (FFFF(defval))
            std::cout << v << " (Default = TEST) : ";
          else
            std::cout << v << " (Default = " << defval << ") : ";
        }
        else
          std::cout << v << " : ";

        // Read the answer
        String str;
        std::getline(std::cin, str);

        // Check for empty line: set to default value
        if (str.empty() && hasDefault)
        {
          answer = defval;
          break;
        }

        // Catch the TEST answer
        if (authTest && str == "TEST")
        {
          answer = TEST;
          break;
        }

        // Try casting in integer
        std::stringstream ss(str);
        if (ss >> answer) break;

        std::cout << "The answer is not a valid double!" << std::endl;
      }
    }
    catch (std::istream::failure& e)
    {
      std::cerr << "Problem when reading double:" << e.what() << std::endl;
    }
    return answer;
  }

  /**
   * Ask interactively for the value of one boolean
   * @param v Text of the question
   * @param defval Default value
   * @param authTest True if a TEST answer is authorized (TEST)
   */
  bool _askBool(const String& v, bool defval, bool authTest)
  {
    DECLARE_UNUSED(authTest);
    bool hasDefault = !isNA(defval);
    bool answer = defval;
    std::cin.exceptions(std::istream::failbit | std::istream::badbit);

    try
    {
      while (true)
      {
        // Display the question
        if (hasDefault)
        {
          String defstr;
          if (defval)
            defstr = "Y";
          else
            defstr = "N";
          std::cout << v << " (Default = " << defstr << ") : ";
        }
        else
          std::cout << v << " : ";

        // Read the answer
        String str;
        std::getline(std::cin, str);

        // Check for empty line: set to default value
        if (str.empty() && hasDefault)
        {
          answer = defval;
          break;
        }

        // Try checking authorized answer
        if (str == "Y")
        {
          answer = true;
          break;
        }
        if (str == "N")
        {
          answer = false;
          break;
        }

        std::cout << "The answer is not a valid bool!" << std::endl;
      }
    }
    catch (std::istream::failure& e)
    {
      std::cerr << "Problem when reading bool:" << e.what() << std::endl;
    }
    return answer;
  }

} // namespace gstlrn
