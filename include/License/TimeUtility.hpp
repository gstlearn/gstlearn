#ifndef TIME_UTILITY_HPP
#define TIME_UTILITY_HPP

// - Standard includes --------------------------------
#include <time.h>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <string>

// -Class definition ----------------------------------

class TimeUtility {
public :
  static bool is_expired(const std::string& expiration_date);
  static void shift_one_year(std::string & expiration_date);

private:
  static void get_current_date(struct tm* timeinfo);
  static bool convert_to_tm(const std::string& date, struct tm* time);
};

#endif // TIME_UTILITY_HPP
