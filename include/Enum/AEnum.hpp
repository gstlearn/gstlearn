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

#include <map>
#include <string>

/****************************************************************************/
/*!
 * Base class for all enums.
 *
 ****************************************************************************/
class AEnum
{
public:
  //! Return the enum key
  const std::string& getKey() const { return _key; }

  //! Return enum as an integer value
  int getValue() const { return _value; }

  //! Return the enum description
  const std::string& getDescr() const { return _descr; }

  //! Cast to constexpr integer (for switch usage)
  constexpr operator int() const { return _value; }

protected:
  AEnum(const std::string& key, int value, const std::string& descr);
  AEnum(const AEnum&) = default;
  ~AEnum() = default;
  AEnum& operator=(const AEnum&) = default;

private:
  std::string _key;
  int         _value;
  std::string _descr;
};

#ifndef SWIG

//********************************
#define ENUM_DECLARE(NAME)\
public:\
  NAME() : AEnum(*_default) {}\
  ~NAME() = default;\
  NAME(const NAME&) = default;\
  NAME& operator=(const NAME&) = default;\
\
  bool operator< (const NAME& e) const { return getValue() <  e.getValue(); }\
  bool operator<=(const NAME& e) const { return getValue() <= e.getValue(); }\
  bool operator> (const NAME& e) const { return getValue() >  e.getValue(); }\
  bool operator>=(const NAME& e) const { return getValue() >= e.getValue(); }\
  bool operator==(const NAME& e) const { return (getValue() == e.getValue()); }\
  bool operator!=(const NAME& e) const { return !operator==(e); }\
\
  static bool existsKey(const std::string& key);\
  static bool existsValue(int value);\
\
  static const NAME& fromKey(const std::string& key);\
  static const NAME& fromValue(int value);\
\
  static NAME fromKey(const std::string& key, const NAME& def);\
  static NAME fromValue(int value, const NAME& def);\
\
  static size_t size() { return _map.size(); }\
\
  class Iterator\
  {\
    friend class NAME;\
\
    Iterator() = delete;\
    Iterator(std::map<int, NAME*>& map) : _iterator(map.begin()), _refmap(map) {}\
  public:\
    ~Iterator() = default;\
    Iterator(const Iterator&) = default;\
    Iterator& operator=(const Iterator&) = default;\
\
    bool isEnd() const { return _iterator == _refmap.end(); }\
    const NAME& next() { return *((_iterator++)->second); }\
    void toFront() { _iterator = _refmap.begin(); }\
    int getValue() { return _iterator->second->getValue(); }\
    const std::string& getKey() { return _iterator->second->getKey(); }\
    const std::string& getDescr() { return _iterator->second->getDescr(); }\
    int first() { return _iterator->first; }\
    NAME* second() { return _iterator->second; }\
\
  private:\
    std::map<int, NAME*>::iterator _iterator;\
    std::map<int, NAME*>&          _refmap;\
  };\
  static Iterator createIterator() { return Iterator(_map); }\
\
private:\
  NAME(const std::string& key, int value, const std::string& descr);\
\
  static const NAME*          _default;\
  static std::map<int, NAME*> _map;\

//********************************
#define ENUM_INIT(NAME)\
std::map<int, NAME*> NAME::_map = std::map<int, NAME*>();\
NAME::NAME(const std::string& key, int value, const std::string& descr) : AEnum(key, value, descr)\
{\
  if (_map.find(value) != _map.end())\
    throw("Duplicated item");\
  _map[value] = this;\
}\
\
bool NAME::existsKey(const std::string& key)\
{\
  auto it = NAME::createIterator();\
  while (!it.isEnd())\
  {\
    if (it.getKey() == key)\
      return true;\
    it.next();\
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
  auto it = NAME::createIterator();\
  while (!it.isEnd())\
  {\
      if (it.getKey() == key)\
      return *(it.second());\
    it.next();\
  }\
  throw("Unknown key");\
}\
\
const NAME& NAME::fromValue(int value)\
{\
  if (existsValue(value))\
      return *_map[value];\
  throw("unknown value");\
}\
\

//*******************************
#define ENUM_DEFINE(NAME, KEY, VALUE, DESCR) NAME NAME::KEY(#KEY, (VALUE), DESCR);

//*******************************
#define ENUM_DEFAULT(NAME, KEY) const NAME* NAME::_default = &NAME::KEY;


#endif // SWIG
