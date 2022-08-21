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
#include "geoslib_f.h"
#include "geoslib_old_f.h"

#include "Basic/NamingConvention.hpp"
#include "Morpho/Morpho.hpp"
#include "Estimation/CalcImage.hpp"
#include "Estimation/KrigingSystem.hpp"
#include "Db/DbGrid.hpp"

CalcImage::CalcImage()
    : ACalcInterpolator(),
      _iattOut(-1),
      _flagFilter(false),
      _flagMorpho(false),
      _oper(EMorpho::UNKNOWN),
      _vmin(0.5),
      _vmax(1.5),
      _option(0),
      _radius(),
      _verbose(false)
{
}

CalcImage::~CalcImage()
{
}

int CalcImage::_getNVar() const
{
  return getDbin()->getVariableNumber();
}

bool CalcImage::_check()
{
  if (! ACalcInterpolator::_check()) return false;

  if (! hasDbin()) return false;
  if (! getDbin()->isGrid())
  {
    messerr("This method requires the Db to be a Grid");
    return false;
  }

  if (_flagFilter)
  {
    if (_getNVar() <= 0)
    {
      messerr("This method requires some Variables to be defined in 'Db'");
      return false;
    }
  }

  if (_flagMorpho)
  {
    if (_getNVar() != 1)
    {
      messerr("This method requires a single Variable to be defined in 'Db'");
      return false;
    }
  }

  return true;
}

bool CalcImage::_preprocess()
{
  _iattOut = _addVariableDb(2, 1, ELoc::UNKNOWN, 1, 0.);
  if (_iattOut < 0) return false;
  return true;
}

bool CalcImage::_postprocess()
{
  if (_flagFilter)
    _renameVariable(2, 1, _iattOut, String(), 1);

  if (_flagMorpho)
    _renameVariable(2, 1, _iattOut, _oper.getKey(), 1);

  return true;
}

void CalcImage::_rollback()
{
  _cleanVariableDb(1);
}

/****************************************************************************/
/*!
 **  Standard Kriging
 **
 ** \return  Error return code
 **
 *****************************************************************************/
bool CalcImage::_run()
{
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbin());

  if (_flagFilter)
  {
    KrigingSystem ksys(dbgrid, dbgrid, getModel(), getNeighparam());
    if (ksys.updKrigOptEstim(_iattOut, -1, -1)) return false;
    if (! ksys.isReady()) return false;

    /* Loop on the targets to be processed */

    for (int iech_out = 0; iech_out < dbgrid->getSampleNumber(); iech_out++)
    {
      mes_process("Image filtering", dbgrid->getSampleNumber(), iech_out);
      if (ksys.estimate(iech_out)) return false;
    }
  }

  if (_flagMorpho)
  {
    if (db_morpho_calc(dbgrid, _iattOut, _oper, _vmin, _vmax, _option, _radius,
                       _verbose)) return false;
  }
  return true;
}

/****************************************************************************/
/*!
 **  Kriging (Factorial) a regular grid
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid     input and output Db grid structure
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure
 ** \param[in]  namconv    Naming Convention
 **
 *****************************************************************************/
int krimage(DbGrid *dbgrid,
            Model *model,
            NeighImage *neighparam,
            const NamingConvention& namconv)
{
  CalcImage image;

  image.setDbin(dbgrid);
  image.setDbout(dbgrid);
  image.setModel(model);
  image.setNeighparam(neighparam);
  image.setNamingConvention(namconv);

  image.setFlagFilter(true);

  // Run the calculator
  int error = (image.run()) ? 0 : 1;
  return error;
}

/**
 * Perform a Morphological operation on an image stored in Db
 * @param dbgrid  Target IN/OUT Db (must be a Grid)
 * @param oper    Type of morphological operation
 * @param vmin    Minimum threshold value
 * @param vmax    Maximum threshold value
 * @param option  Option
 * @param radius  Radius
 * @param verbose Verbose option
 * @param namconv Naming convention
 * @return
 */
GSTLEARN_EXPORT int morpho(DbGrid *dbgrid,
                           const EMorpho &oper,
                           double vmin,
                           double vmax,
                           int option,
                           const VectorInt &radius,
                           bool verbose,
                           const NamingConvention &namconv)
{
  CalcImage image;

  image.setDbin(dbgrid);
  image.setDbout(dbgrid);
  image.setNamingConvention(namconv);

  image.setFlagMorpho(true);
  image.setOper(oper);
  image.setVmin(vmin);
  image.setVmax(vmax);
  image.setOption(option);
  image.setRadius(radius);
  image.setVerbose(verbose);

  // Run the calculator
  int error = (image.run()) ? 0 : 1;
  return error;
}
