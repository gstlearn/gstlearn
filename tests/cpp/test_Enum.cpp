/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "Basic/RepeatMacro.hpp"
#include "Basic/File.hpp"
#include "Enum/AEnum.hpp"

#include <string>
#include <iostream>
#include <sstream>

#define PRINT_THAT(X) std::cout << #X << std::endl; 
#define PRINT_ALL_THAT(...) REPEAT(PRINT_THAT, __VA_ARGS__)
#define PRINT2_THAT(X,Y) std::cout << #X << "|" << #Y << std::endl; 
#define PRINT2_ALL_THAT(...) REPEAT2(PRINT2_THAT, __VA_ARGS__)

// Under Windows, prevent from declspec(import) next enums
#ifdef GSTLEARN_EXPORT
  #undef GSTLEARN_EXPORT
  #define GSTLEARN_EXPORT
#endif

// This must be put in header files (hpp)
#define ENUM_FRUIT Fruit, APPLE,\
                   APPLE,  1, "Apple is the best",\
                   PEAR,   2, "Pear is juicy",\
                   BANANA, 4, "Banana is strange"
ENUM_DECLARE(ENUM_FRUIT)

#define ENUM_DAY EDay, MONDAY,\
                 MONDAY   , 12, "Monday is the first day",\
                 TUESDAY  , 13, "Tuesday",\
                 WEDNESDAY, 36, "Wednesday",\
                 THURSDAY , 15, "Thursday",\
                 FRIDAY   , 46, "Friday",\
                 SATURDAY , 49, "Saturday",\
                 SUNDAY   , 47, "Sunday"
ENUM_DECLARE(ENUM_DAY)


// This must be put in body source files (cpp)
ENUM_DEFINE(ENUM_FRUIT)

ENUM_DEFINE(ENUM_DAY)



int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  EDay d1 = EDay::SUNDAY;
  std::cout << "d1=" << "Enum#" << d1.getValue() << ": " << d1.getKey() << "|" << d1.getDescr() << std::endl;
  EDay d2;
  std::cout << "d2=" << "Enum#" << d2.getValue() << ": " << d2.getKey() << "|" << d2.getDescr() << std::endl;
  EDay d3 = EDay::fromValue(5);
  std::cout << "d3=" << "Enum#" << d3.getValue() << ": " << d3.getKey() << "|" << d3.getDescr() << std::endl;
  
  auto it = EDay::getIterator();
  while (it.hasNext())
  {
    EDay d = *it;
    std::cout << "Enum#" << it.getValue() << ": " << d.getKey() << "|" << (*it).getDescr() << std::endl;
    it.toNext();
  }

  std::string day("TOTODAY");
  if (EDay::existsKey(day))
    std::cout << day << " exists" << std::endl;
  else
    std::cout << day << " doesn't exists" << std::endl;

  if (d2 == EDay::MONDAY)
    std::cout << "I hate d2!" << std::endl;

  switch(d1.toEnum())
  {
    case EDay::E_SUNDAY:
    case EDay::E_SATURDAY:
    {
      std::cout << "d1 is amazing!" << std::endl;
      break;
    }
    default:
    {
      std::cout << "d1 is boring!" << std::endl;
      break;
    }
  }
  
  std::cout << "Test Banana : " << Fruit::BANANA.getDescr() << std::endl;
  
  std::cout << "NARG=" << NARG3(a,b,c,d,e,f,g,h,i) << std::endl;
  
  PRINT_ALL_THAT(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,Bonjour,19,20,21,22,23,24,25,26,27,28,29,30,31,32)
  PRINT2_ALL_THAT(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,Bonjour,19,20,21,22,23,24,25,26,27,28,29,30,31,32,
                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,Bonjour,19,20,21,22,23,24,25,26,27,28,29,30,31,32)
  
  return 0;
}
