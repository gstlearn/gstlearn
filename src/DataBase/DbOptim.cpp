/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2026) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "DataBase/DbOptim.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <utility>
#include <vector>

namespace gstlrn
{
  std::unordered_map<std::string, long long> DbOptim::_counters{};
  bool DbOptim::_atExitRegistered = false;

  void DbOptim::incrementCounter(const String& name)
  {
    if constexpr (areCountersActive)
    {
      // Register automatic summary output once on process exit
      if (!_atExitRegistered)
      {
        _atExitRegistered = true;
        std::atexit([]() { DbOptim::displayCounters(); });
      }
      _counters[name] += 1;
    }
  }

  void DbOptim::resetCounters()
  {
    if constexpr (areCountersActive)
    {
      _counters.clear();
    }
  }

  void DbOptim::displayCounters()
  {
    if constexpr (areCountersActive)
    {
      if (_counters.empty())
      {
        std::fprintf(stderr, "No counter hits recorded.\n");
        std::fflush(stderr);
        return;
      }

      // Sort entries by descending hit count to highlight hotspots first
      std::vector<std::pair<std::string, long long>> sorted(
        _counters.begin(), _counters.end());
      std::sort(
        sorted.begin(), sorted.end(),
        [](const auto& a, const auto& b) { return a.second > b.second; });

      std::fprintf(stderr, "\n=== Loop / Node Pass Counters ===\n");
      for (const auto& [name, count]: sorted)
      {
        std::fprintf(stderr, "  %-65s : %lld\n", name.c_str(), count);
      }
      std::fflush(stderr);
    }
  }

} // namespace gstlrn
