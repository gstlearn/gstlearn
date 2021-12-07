#pragma once

#include "gstlearn_export.hpp"

// - Standard includes --------------------------------
//#include <time.h>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <string>

// -Class definition ----------------------------------

class GSTLEARN_EXPORT TimeUtility {
public :
  static bool is_expired(const std::string& expiration_date);
  static void shift_one_year(std::string & expiration_date);

private:
  static void get_current_date(struct tm* timeinfo);
  static bool convert_to_tm(const std::string& date, struct tm* time);
};

