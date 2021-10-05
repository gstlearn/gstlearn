#ifndef  INTERFACE_D_HPP
#define  INTERFACE_D_HPP

#include "geoslib_d.h"
#include "Enum/AEnum.hpp"
#include <vector>

#define UNDEF_DOUBLE 1.234e30
#define UNDEF_STRING "??"
#define UNDEF_INT -2147483648

#define UNDEF_CAT_VAL UNDEF_INT
#define UNDEF_CAT_LABEL ""

/*
typedef enum 
{
  ROLE_NOROLE = -1,
  ROLE_COORD  = 0,
  ROLE_Z      = 1,
  ROLE_CODE   = 9,  
  ROLE_SEL    = 10,
  ROLE_MAX
} Roles;

typedef enum
{
  CALCUL_BY_LAG = 0,
  CALCUL_BY_SAMPLE = 1,
} CalculRules;
*/

#define ENUM_ROLES ERoles, NOROLE, \
                   NOROLE, -1, "Equivalent to ELoc::UNKNOWN", \
                   COORD,   0, "Equivalent to ELoc::X", \
                   Z,       1, "Equivalent to ELoc::Z", \
                   CODE,    9, "Equivalent to ELoc::C", \
                   SEL,    10, "Equivalent to ELoc::SEL"

ENUM_DECLARE(ENUM_ROLES)

#define ENUM_CALC_RULES ECalcRules, CALCUL_BY_LAG, \
                        CALCUL_BY_LAG,    0, "Calculation by lag", \
                        CALCUL_BY_SAMPLE, 1, "Calculation by sample"

ENUM_DECLARE(ENUM_CALC_RULES)

#endif
