#ifndef  INTERFACE_D_HPP
#define  INTERFACE_D_HPP

#include "geoslib_d.h"
#include <vector>

#define UNDEF_DOUBLE 1.234e30
#define UNDEF_STRING "??"
#define UNDEF_INT -2147483648

#define UNDEF_CAT_VAL UNDEF_INT
#define UNDEF_CAT_LABEL ""

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

#endif
