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
#include "Neigh/ENeigh.hpp"
#include "Db/ELoadBy.hpp"
#include "Db/ELoc.hpp"
#include "Covariances/ECov.hpp"
#include "Covariances/ECalcMember.hpp"
#include "Drifts/EDrift.hpp"
#include "Model/EModelProperty.hpp"
#include "Model/EConsElem.hpp"
#include "Model/EConsType.hpp"
#include "Variogram/ECalcVario.hpp"
#include "LithoRule/ERule.hpp"
#include "Enum/EKrigOpt.hpp"
#include "Anamorphosis/EAnam.hpp"
#include "Basic/EJustify.hpp"
#include "LithoRule/EProcessOper.hpp"

ENUM_DEFINE(ENUM_LOAD_BY)
ENUM_DEFINE(ENUM_NEIGH)
ENUM_DEFINE(ENUM_LOC)
ENUM_DEFINE(ENUM_COV)
ENUM_DEFINE(ENUM_DRIFT)
ENUM_DEFINE(ENUM_CALC_VARIO)
ENUM_DEFINE(ENUM_MODEL_PROPERTY)
ENUM_DEFINE(ENUM_RULE)
ENUM_DEFINE(ENUM_CALC_MEMBER)
ENUM_DEFINE(ENUM_KRIG_OPT)
ENUM_DEFINE(ENUM_ANAM)
ENUM_DEFINE(ENUM_CONS_ELEM)
ENUM_DEFINE(ENUM_CONS_TYPE)
ENUM_DEFINE(ENUM_JUSTIFY)
ENUM_DEFINE(ENUM_PROCESS_OPER)
