#include "License/TimeUtility.hpp"
#include "Basic/AStringable.hpp"

#include <string.h>
#include <sstream>

const std::string day_keys = "AZE9RTY5UPM8LKJ4HGF2D6S3W7XCVBNQ";

//
// Check if date is expired
//
bool TimeUtility::is_expired(const std::string &expiration_date) 
{
  struct tm exp_time;
  if(convert_to_tm(expiration_date, &exp_time)) 
  {
    // Check if expired comparing to current time
    struct tm cur_time;
    get_current_date(&cur_time);
    if( cur_time.tm_year < exp_time.tm_year      || // year before
        (cur_time.tm_mon  < exp_time.tm_mon  && 
         cur_time.tm_year == exp_time.tm_year   ) || // same year
        (cur_time.tm_mday < exp_time.tm_mday && 
         cur_time.tm_mon  == exp_time.tm_mon    ) || // same month
        (true                                && 
         cur_time.tm_mday == exp_time.tm_mday &&     // same day
         cur_time.tm_mon  == exp_time.tm_mon  &&
         cur_time.tm_year == exp_time.tm_year))
      return false;
  }
  return true;
}

bool TimeUtility::convert_to_tm(const std::string& date, 
                                struct tm* time) 
{
  // Convert string to integers
  std::string in_date = date;
  std::replace(in_date.begin(), in_date.end(), '/', ' ');
  int day, month, year;
  std::istringstream strdate(in_date);
  strdate >> day >> month >> year;

  // Fill tm structure
  get_current_date(time);
  time->tm_mday = day;
  time->tm_mon  = month - 1;
  time->tm_year = year  - 1900;

  // Check date validity
  if(mktime(time) == -1) {
    messerr("Invalid date");
    return false;
  }
  return true;
}

void TimeUtility::get_current_date(struct tm* timeinfo) {
  time_t cur = time(NULL);
  struct tm* cur_time;

#ifdef _WIN32
  locatime_s(cur_time, &cur);
#else
  cur_time = localtime(&cur);
#endif
  memcpy(timeinfo, cur_time, sizeof(struct tm));
}

void TimeUtility::shift_one_year(std::string & date) 
{
  struct tm cur_time;

  get_current_date(&cur_time);

  std::stringstream out;
  out << cur_time.tm_mday         << "/" 
      << cur_time.tm_mon+1        << "/" 
      << cur_time.tm_year+1900+1;
  date = out.str();
}
