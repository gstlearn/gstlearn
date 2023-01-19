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

// WARNING: Make this include list as small as possible!
#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Basic/VectorT.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/RepeatMacro.hpp"
#include "Basic/WarningMacro.hpp"
#include "Basic/String.hpp"

#include <iostream>
#include <string>
#include <iomanip>
#include <map>
#include <vector>

class GSTLEARN_EXPORT AEnum
{
public:
  //! Return the enum key as a string (max 10 characters)
  const String& getKey() const { return _key; }

  //! Return enum value as an integer value (max 32 enum)
  int getValue() const { return _value; }

  //! Return the enum description as a string
  const String& getDescr() const { return _descr; }

#ifndef SWIG
  // Remove this: too much dangerous (implicit casts)
  // => Force compilation error where enum were used as integer
  //operator int() const { return _value; }
#endif

  bool operator< (const AEnum& e) const { return getValue() <  e.getValue(); }
  bool operator<=(const AEnum& e) const { return getValue() <= e.getValue(); }
  bool operator> (const AEnum& e) const { return getValue() >  e.getValue(); }
  bool operator>=(const AEnum& e) const { return getValue() >= e.getValue(); }
  bool operator==(const AEnum& e) const { return getValue() == e.getValue(); }
  bool operator!=(const AEnum& e) const { return !operator==(e); }

  bool isSmaller       (const AEnum& e) const { return e <  *this; }  
  bool isSmallerOrEqual(const AEnum& e) const { return e <= *this; }  
  bool isGreater       (const AEnum& e) const { return e >  *this; }
  bool isGreaterOrEqual(const AEnum& e) const { return e >= *this; }
  bool isEqual         (const AEnum& e) const { return e == *this; }
  bool isDifferent     (const AEnum& e) const { return e != *this; }

  void printEnum() const;

protected:
  AEnum(const String& key, int value, const String& descr)
  : _key(key), _value(value), _descr(descr)
  {
  }
  AEnum(const AEnum&) = default;
  ~AEnum() = default;
  AEnum& operator=(const AEnum&) = default;

  template<typename ... Args>
  static void _printMsg(const char *format, Args... args);
  // TODO: Should be used for the error messages (not possible as this function
  // (template), as it is present in a macro (thus expanded)
  // would be embarked in SWIG... who does not know it.
 
private:
  String _key;
  int    _value;
  String _descr;
};

#define ENUM_ITEM(NAME, x,y,z) E_ ## x = y,
#define ENUM_ITEMS(NAME, ...) EXPAND(REPEAT3(ENUM_ITEM, NAME, __VA_ARGS__))

#define ENUM_DECL(NAME, x,y,z) static const NAME x;
#define ENUM_DECLS(NAME, ...) EXPAND(REPEAT3(ENUM_DECL, NAME, __VA_ARGS__))

#define ENUM_IMPL(NAME, x,y,z) const NAME NAME::x = NAME(#x, y, z);
#define ENUM_IMPLS(NAME, ...) EXPAND(REPEAT3(ENUM_IMPL, NAME, __VA_ARGS__))

// ######################
//      ENUM DECLARE
// ######################
#define ENUM_DECLARE_(NAME, DEFAULT, ...)\
class NAME;\
\
typedef std::map<int, NAME*> NAME ## Map;\
\
class GSTLEARN_EXPORT NAME ## Iterator\
{\
  friend class NAME;\
\
  NAME ## Iterator() = delete;\
  NAME ## Iterator(NAME ## Map* map);\
public:\
  ~NAME ## Iterator() = default;\
  NAME ## Iterator(const NAME ## Iterator&) = default;\
  NAME ## Iterator& operator=(const NAME ## Iterator&) = default;\
\
  const NAME& operator*() const;\
  bool hasNext() const;\
  const NAME& toNext();\
  const NAME& toFront();\
  const NAME& getEnum() const;\
  int getValue() const;\
  const String& getKey() const;\
  const String& getDescr() const;\
\
private:\
  NAME ## Map::iterator _stditer;\
  NAME ## Map*          _refmap;\
};\
\
class GSTLEARN_EXPORT NAME : public AEnum\
{\
\
public:\
  NAME();\
  ~NAME();\
  NAME(const NAME&) = default;\
  NAME(int value);\
  NAME(const String& key);\
  NAME& operator=(const NAME&) = default;\
\
  static size_t getSize();\
  static NAME ## Iterator getIterator();\
  static void printAll();\
\
  static bool existsKey(const String& key);\
  static bool existsValue(int value);\
  static const NAME& fromKey(const String& key);\
  static const NAME& fromValue(int value);\
  static std::vector<NAME> fromKeys(const VectorString& keys);\
  static std::vector<NAME> fromValues(const VectorInt& values);\
\
private:\
  NAME(const String& key, int value, const String& descr);\
\
  static NAME ## Map      _map;\
  static NAME ## Iterator _iterator;\
  static const NAME*      _default;\
\
public:\
  enum E ## NAME\
  {\
    EXPAND(ENUM_ITEMS(NAME, __VA_ARGS__))\
  };\
  E ## NAME toEnum() const;\
\
  EXPAND(ENUM_DECLS(NAME, __VA_ARGS__))\
};\
\


// ######################
//       ENUM DEFINE
// ######################
#define ENUM_DEFINE_(NAME, DEFAULT, ...)\
NAME ## Map NAME::_map = NAME ## Map();\
NAME ## Iterator NAME::_iterator = NAME ## Iterator(&NAME::_map);\
\
const NAME* NAME::_default = &NAME::DEFAULT;\
\
NAME::NAME()\
: AEnum(*_default)\
{\
}\
\
NAME::NAME(int value)\
: AEnum(fromValue(value))\
{\
}\
\
NAME::NAME(const String& key)\
: AEnum(fromKey(key))\
{\
}\
\
NAME::NAME(const String& key, int value, const String& descr)\
: AEnum(key, value, descr)\
{\
  if (_map.find(value) != _map.end())\
    throw("Duplicated item");\
  _map[value] = this;\
}\
\
NAME::~NAME()\
{\
}\
\
size_t NAME::getSize()\
{\
  return _map.size();\
}\
\
NAME ## Iterator NAME::getIterator()\
{\
  auto it(_iterator);\
  it.toFront();\
  return it;\
}\
\
void NAME::printAll()\
{\
  auto it(getIterator());\
  while (it.hasNext())\
  {\
    (*it).printEnum();\
    it.toNext();\
  }\
}\
\
bool NAME::existsKey(const String& key)\
{\
  auto it = _map.begin();\
  while (it != _map.end())\
  {\
    if (it->second->getKey() == key)\
      return true;\
    it++;\
  }\
  return false;\
}\
\
bool NAME::existsValue(int value)\
{\
  return (_map.find(value) != _map.end());\
}\
\
const NAME& NAME::fromKey(const String& key)\
{\
  auto it = _map.begin();\
  while (it != _map.end())\
  {\
    if (it->second->getKey() == toUpper(key))\
      return (*(it->second));\
    it++;\
  }\
  std::cout << "Unknown key " << key << " for enum " << #NAME << std::endl;\
  return *_default;\
}\
\
const NAME& NAME::fromValue(int value)\
{\
  if (existsValue(value))\
    return (*(_map[value]));\
  std::cout << "Unknown value " << value << " for enum " << #NAME << std::endl;\
  return *_default;\
}\
\
std::vector<NAME> NAME::fromKeys(const VectorString& keys)\
{\
  std::vector<NAME> vec;\
  for (auto v : keys)\
    vec.push_back(fromKey(v));\
  return vec;\
}\
\
std::vector<NAME> NAME::fromValues(const VectorInt& values)\
{\
  std::vector<NAME> vec;\
  for (auto v : values)\
    vec.push_back(fromValue(v));\
  return vec;\
}\
\
NAME::E ## NAME NAME::toEnum() const\
{\
  return static_cast<E ## NAME>(getValue());\
}\
\
EXPAND(ENUM_IMPLS(NAME, __VA_ARGS__))\
\
NAME ## Iterator::NAME ## Iterator(NAME ## Map* map) \
: _stditer(map->begin())\
, _refmap(map)\
{\
}\
\
const NAME& NAME ## Iterator::operator*() const\
{\
  return (*(_stditer->second));\
}\
\
bool NAME ## Iterator::hasNext() const\
{\
  return (_stditer != _refmap->end());\
}\
\
const NAME& NAME ## Iterator::toNext()\
{\
  return (*((_stditer++)->second));\
}\
\
const NAME& NAME ## Iterator::toFront()\
{\
  _stditer = _refmap->begin();\
  return (*(_stditer->second));\
}\
\
const NAME& NAME ## Iterator::getEnum() const\
{\
  return (*(_stditer->second));\
}\
\
int NAME ## Iterator::getValue() const\
{\
  return (_stditer->second->getValue());\
}\
\
const String& NAME ## Iterator::getKey() const\
{\
  return (_stditer->second->getKey());\
}\
\
const String& NAME ## Iterator::getDescr() const\
{\
  return (_stditer->second->getDescr());\
}\
\

// Top level macros
#define ENUM_DECLARE(...) EXPAND(ENUM_DECLARE_(__VA_ARGS__))
#define ENUM_DEFINE(...)  EXPAND(ENUM_DEFINE_(__VA_ARGS__))
