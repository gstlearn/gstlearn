/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f_private.h"
#include "geoslib_old_f.h"

#include "Basic/NamingConvention.hpp"
#include "Morpho/Morpho.hpp"
#include "Estimation/CalcImage.hpp"
#include "Estimation/KrigingSystem.hpp"
#include "Neigh/NeighImage.hpp"
#include "Db/DbGrid.hpp"

CalcImage::CalcImage()
    : ACalcInterpolator(),
      _iattOut(-1),
      _flagFilter(false),
      _flagMorpho(false),
      _nvarMorpho(1),
      _oper(EMorpho::UNKNOWN),
      _vmin(0.5),
      _vmax(1.5),
      _option(0),
      _radius(),
      _distErode(false),
      _verbose(false),
      _flagSmooth(false),
      _smoothType(0),
      _smoothRange(1.)
{
}

CalcImage::~CalcImage()
{
}

int CalcImage::_getNVar() const
{
  return getDbin()->getLocNumber(ELoc::Z);
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

  if (_flagSmooth)
  {
    if (_smoothType != 1 && _smoothType != 2)
    {
      messerr("Filtering 'type' should be 1 or 2");
      return false;
    }
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
  int nvar = _getNVar();
  if (_flagFilter)
    _iattOut = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, nvar, 0.);

  if (_flagMorpho)
    _iattOut = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, _nvarMorpho, 0.);

  if (_flagSmooth)
    _iattOut = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, 1, 0.);

  if (_iattOut < 0) return false;
  return true;
}

bool CalcImage::_postprocess()
{
  /* Free the temporary variables */
  _cleanVariableDb(2);

  if (_flagFilter)
    _renameVariable(2, getDbin()->getLocNumber(ELoc::Z), _iattOut, String(), 1);

  if (_flagMorpho)
    _renameVariable(2, 1, _iattOut, _oper.getKey(), _nvarMorpho);

  if (_flagSmooth)
    _renameVariable(2, 1, _iattOut, String(), 1);

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
                       _distErode, _verbose)) return false;
  }

  if (_flagSmooth)
  {
    NeighImage* neighI = dynamic_cast<NeighImage*>(getNeighparam());
    _image_smoother(dbgrid, neighI, _smoothType, _smoothRange, _iattOut);
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

/****************************************************************************/
/*!
 **  Smooth a regular grid
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid     input and output Db grid structure
 ** \param[in]  neighparam Neigh structure
 ** \param[in]  type       1 for Uniform; 2 for Gaussian
 ** \param[in]  range      Range (used for Gaussian only)
 ** \param[in]  namconv    Naming Convention
 **
 *****************************************************************************/
int dbSmoother(DbGrid *dbgrid,
               NeighImage *neighparam,
               int type,
               double range,
               const NamingConvention &namconv)
{
  CalcImage image;

  image.setDbin(dbgrid);
  image.setDbout(dbgrid);
  image.setNeighparam(neighparam);
  image.setNamingConvention(namconv);

  image.setFlagSmooth(true);
  image.setSmoothType(type);
  image.setSmoothRange(range);

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
 * @param flagDistErode True: Inflate the grain; False: Reduce the grain
 * @param namconv Naming convention
 * @return
 */
GSTLEARN_EXPORT int dbMorpho(DbGrid *dbgrid,
                             const EMorpho &oper,
                             double vmin,
                             double vmax,
                             int option,
                             const VectorInt &radius,
                             bool flagDistErode,
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
  image.setDistErode(flagDistErode);
  image.setVerbose(verbose);

  // Particular case of the number of output variables
  int nvar = 1;
  if (oper == EMorpho::GRADIENT) nvar = dbgrid->getNDim();
  image.setNvarMorpho(nvar);

  // Run the calculator
  int error = (image.run()) ? 0 : 1;
  return error;
}
