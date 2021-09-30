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

#include <iostream>
#include <string>
#include <map>

#include "Basic/RepeatMacro.hpp"

class AEnum
{
public:
  //! Return the enum key as a string
  const std::string& getKey() const { return _key; }

  //! Return enum value as an integer value
  int getValue() const { return _value; }

  //! Return the enum description as a string
  const std::string& getDescr() const { return _descr; }

#ifndef SWIG
  //! Cast to constexpr integer (for switch usage) (-std=c++11)
  operator int() const { return _value; }
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

protected:
  AEnum(const std::string& key, int value, const std::string& descr)
  : _key(key), _value(value), _descr(descr)
  {
  }
  AEnum(const AEnum&) = default;
  ~AEnum() = default;
  AEnum& operator=(const AEnum&) = default;
 
private:
  std::string _key;
  int         _value;
  std::string _descr;
};

#define ENUM_ITEM(NAME, X,Y,Z) E_ ## X = Y,
#define ENUM_ITEMS(NAME, ...) EXPAND(REPEAT3(ENUM_ITEM, NAME, __VA_ARGS__))

#define ENUM_DECL(NAME, X,Y,Z) static const NAME X;
#define ENUM_DECLS(NAME, ...) EXPAND(REPEAT3(ENUM_DECL, NAME, __VA_ARGS__))

#define ENUM_IMPL(NAME, X,Y,Z) const NAME NAME::X = NAME(#X, Y, Z);
#define ENUM_IMPLS(NAME, ...) EXPAND(REPEAT3(ENUM_IMPL, NAME, __VA_ARGS__))

#define ENUM_DECLARE_(NAME, DEFAULT, ...)\
class NAME ## Iterator;\
class NAME;\
\
typedef std::map<int, NAME*> NAME ## Map;\
\
class NAME : public AEnum\
{\
\
public:\
  NAME() : AEnum(*_default) {}\
  ~NAME();\
  NAME(const NAME&) = default;\
  NAME& operator=(const NAME&) = default;\
\
  static size_t getSize();\
  static NAME ## Iterator* getIterator();\
  \
  static bool existsKey(const std::string& key);\
  static bool existsValue(int value);\
  static const NAME& fromKey(const std::string& key);\
  static const NAME& fromValue(int value);\
\
private:\
  NAME(const std::string& key, int value, const std::string& descr);\
\
  static const NAME*   _default;\
  static NAME ## Map       _map;\
  static NAME ## Iterator _iterator;\
\
public:\
  enum E ## NAME\
  {\
    EXPAND(ENUM_ITEMS(NAME, __VA_ARGS__))\
  };\
  E ## NAME toEnum() const { return static_cast<E ## NAME>(getValue()); }\
\
  EXPAND(ENUM_DECLS(NAME, __VA_ARGS__))\
};\
\
class NAME ## Iterator\
{\
  friend class NAME;\
\
  NAME ## Iterator() = delete;\
  NAME ## Iterator(NAME ## Map& map) : _iterator(map.begin()), _refmap(map) {}\
public:\
  ~NAME ## Iterator() = default;\
  NAME ## Iterator(const NAME ## Iterator&) = default;\
  NAME ## Iterator& operator=(const NAME ## Iterator&) = default;\
\
  const NAME& operator*() const { return (*(_iterator->second)); }\
  bool hasNext() const { return (_iterator != _refmap.end()); }\
  const NAME& toNext() { return (*(_iterator++)->second); }\
  const NAME& toFront() { _iterator = _refmap.begin(); return (*(_iterator->second));}\
  int getValue() const { return (_iterator->second->getValue()); }\
  const std::string& getKey() const { return (_iterator->second->getKey()); }\
  const std::string& getDescr() const { return (_iterator->second->getDescr()); }\
\
private:\
  NAME ## Map::iterator _iterator;\
  NAME ## Map&          _refmap;\
};\
\

#define ENUM_DEFINE_(NAME, DEFAULT, ...)\
NAME ## Map NAME::_map = NAME ## Map();\
NAME ## Iterator NAME::_iterator = NAME ## Iterator(NAME::_map);\
\
NAME::NAME(const std::string& key, int value, const std::string& descr)\
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
NAME ## Iterator* NAME::getIterator()\
{\
  _iterator.toFront();\
  return &_iterator;\
}\
\
bool NAME::existsKey(const std::string& key)\
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
const NAME& NAME::fromKey(const std::string& key)\
{\
  auto it = _map.begin();\
  while (it != _map.end())\
  {\
    if (it->second->getKey() == key)\
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
EXPAND(ENUM_IMPLS(NAME, __VA_ARGS__))\
\
const NAME* NAME::_default = &NAME::DEFAULT;\
\

// Top level macros
#define ENUM_DECLARE(...) EXPAND(ENUM_DECLARE_(__VA_ARGS__))
#define ENUM_DEFINE(...) EXPAND(ENUM_DEFINE_(__VA_ARGS__))

