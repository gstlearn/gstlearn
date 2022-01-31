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
#include "Basic/ECst.hpp"
#include "Basic/EDbg.hpp"
#include "Neigh/ENeigh.hpp"
#include "Db/ELoadBy.hpp"
#include "Db/ELoc.hpp"
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
#include "LinearOp/EPowerPT.hpp"
#include "API/ESPDECalcMode.hpp"
#include "Covariances/ECov.hpp"
#include "Covariances/ECalcMember.hpp"
#include "Covariances/ETape.hpp"
#include "Covariances/EConvType.hpp"
#include "Covariances/EConvDir.hpp"

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
ENUM_DEFINE(ENUM_DEBUG)
ENUM_DEFINE(ENUM_CST)
ENUM_DEFINE(ENUM_PROCESS_OPER)
ENUM_DEFINE(ENUM_POWER_PT)
ENUM_DEFINE(ENUM_SPDE_CALC_MODE)
ENUM_DEFINE(ENUM_TAPE)
ENUM_DEFINE(ENUM_CONVTYPE)
ENUM_DEFINE(ENUM_CONVDIR)

