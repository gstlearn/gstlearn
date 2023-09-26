/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "geoslib_define.h"

#include "Db/Db.hpp"
#include "Basic/Grid.hpp"
#include "Calculators/CalcGridToGrid.hpp"

CalcGridToGrid::CalcGridToGrid()
    : ACalcDbToDb(),
      _iattOut(-1),
      _flagCopy(false),
      _flagExpand(false),
      _flagShrink(false),
      _iattAux(-1),
      _flagInter(false),
      _nameTops(),
      _nameBots()
{
}

CalcGridToGrid::~CalcGridToGrid()
{
}

int CalcGridToGrid::_getNVar() const
{
  int nvar = 0;
  if (getDbin() != nullptr) nvar = getDbin()->getLocNumber(ELoc::Z);
  return nvar;
}

bool CalcGridToGrid::_check()
{
  /*************************************/
  /* Both Files are compulsory as Grid */
  /*************************************/

  if (! hasDbin())  return false;
  if (! hasDbout()) return false;
  if (! isGridIn()) return false;
  if (! isGridOut()) return false;

  /**************************************************/
  /* Cross-checking the Space Dimension consistency */
  /**************************************************/

  if (! getGridin()->getGrid().isSame(getGridout()->getGrid()))
  {
    messerr("The two Grids do not share the same common dimensions");
    return false;
  }

  /**************************************************/
  /* Cross-Checking the Variable Number consistency */
  /**************************************************/

  int nvar = _getNVar();
  int nvar_requested = (_flagInter) ? 2 : 1;
  if (nvar != nvar_requested)
  {
    messerr("This application requires %d variable(s) to be defined",nvar_requested);
    return false;
  }

  /******************/
  /* Specific check */
  /******************/

  if (_flagCopy)
  {
    if (_compareInMinusOut() != 0)
    {
      messerr("The two files should have the same Space Dimension");
      return false;
    }
  }

  if (_flagExpand)
  {
    if (_compareInMinusOut() >= 0)
    {
      messerr("The Space dimension of 'dbout' should be larger then the one of 'dbin'");
      return false;
    }
  }

  if (_flagShrink)
  {
    if (_compareInMinusOut() <= 0)
    {
      messerr("The Space dimension of 'dbout' should be smaller then the one of 'dbin'");
      return false;
    }
  }

  if (_flagInter)
  {
    if (_compareInMinusOut() >= 0)
    {
      messerr("The Space dimension of 'dbout' should be larger then the one of 'dbin'");
      return false;
    }
    int ndiff = -_compareInMinusOut();
    if (ndiff != (int) _nameBots.size())
    {
      messerr("The argument 'nameBots' (%d) should be dimensioned to %d",
              (int) _nameBots.size(), ndiff);
      return false;
    }
    if (ndiff != (int) _nameTops.size())
    {
      messerr("The argument 'nameTops' (%d) should be dimensioned to %d",
              (int) _nameTops.size(), ndiff);
      return false;
    }
  }
  return true;
}

int CalcGridToGrid::_compareInMinusOut() const
{
  return getDbin()->getNDim() - getDbout()->getNDim();
}

bool CalcGridToGrid::_preprocess()
{
  _iattOut = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, 1, 0.);
  if (_iattOut < 0) return false;

  if (_flagShrink)
  {
    _iattAux = _addVariableDb(2, 2, ELoc::UNKNOWN, 0, 1, 0.);
     if (_iattAux < 0) return false;
  }

  return true;
}

bool CalcGridToGrid::_postprocess()
{
  /* Free the temporary variables */
  _cleanVariableDb(2);

  _renameVariable(2, VectorString(), ELoc::Z, 1, _iattOut, String(), 1);
  return true;
}

void CalcGridToGrid::_rollback()
{
  _cleanVariableDb(1);
}

bool CalcGridToGrid::_run()
{
  if (_flagCopy)
    return _g2gCopy();

  if (_flagExpand)
    return _g2gExpand();

  if (_flagShrink)
    return _g2gShrink();

  if (_flagInter)
    return _g2gInter();

  return false;
}

int dbg2gCopy(DbGrid *dbin,
              DbGrid *dbout,
              const NamingConvention &namconv)
{
  CalcGridToGrid calcul;
  calcul.setDbin(dbin);
  calcul.setDbout(dbout);
  calcul.setNamingConvention(namconv);

  calcul.setFlagCopy(true);

  // Run the calculator
  int error = (calcul.run()) ? 0 : 1;
  return error;
}

int dbg2gExpand(DbGrid *dbin,
                DbGrid *dbout,
                const NamingConvention &namconv)
{
  CalcGridToGrid calcul;
  calcul.setDbin(dbin);
  calcul.setDbout(dbout);
  calcul.setNamingConvention(namconv);

  calcul.setFlagExpand(true);

  // Run the calculator
  int error = (calcul.run()) ? 0 : 1;
  return error;
}

int dbg2gInterpolate(DbGrid *dbin,
                     DbGrid *dbout,
                     const VectorString& tops,
                     const VectorString& bots,
                     const NamingConvention &namconv)
{
  CalcGridToGrid calcul;
  calcul.setDbin(dbin);
  calcul.setDbout(dbout);
  calcul.setNamingConvention(namconv);

  calcul.setFlagInter(true);
  calcul.setNameBots(bots);
  calcul.setNameTops(tops);

  // Run the calculator
  int error = (calcul.run()) ? 0 : 1;
  return error;
}

int dbg2gShrink(DbGrid *dbin,
                DbGrid *dbout,
                const NamingConvention &namconv)
{
  CalcGridToGrid calcul;
  calcul.setDbin(dbin);
  calcul.setDbout(dbout);
  calcul.setNamingConvention(namconv);

  calcul.setFlagShrink(true);

  // Run the calculator
  int error = (calcul.run()) ? 0 : 1;
  return error;
}

/**
 * Copy the cell indices of 'indgIn' into 'indgOut'.
 * It is assumed that the dimension of 'indgOut' is always smaller than the one of 'indgIn'
 * No test if performed for efficiency reason.
 * @param indgIn  Input vector of indices
 * @param indgOut Output vector of indices
 */
void CalcGridToGrid::_reduceIndices(const VectorInt& indgIn, VectorInt& indgOut)
{
  int ndimOut = (int) indgOut.size();

  for (int i = 0; i < ndimOut; i++)
    indgOut[i] = indgIn[i];
}

bool CalcGridToGrid::_g2gCopy()
{
  int nech = getDbin()->getSampleNumber();

  for (int iech = 0; iech < nech; iech++)
  {
    if (! getDbin()->isActive(iech)) continue;
    double value = getDbin()->getLocVariable(ELoc::Z,iech, 0);
    getDbout()->setArray(iech, _iattOut, value);
  }
  return true;
}

bool CalcGridToGrid::_g2gExpand()
{
  int ndim_in  = getDbin()->getNDim();
  int ndim_out = getDbout()->getNDim();
  VectorInt indgIn(ndim_in);
  VectorInt indgOut(ndim_out);

  // Loop on the output file
  for (int iech_out = 0; iech_out < getDbout()->getActiveSampleNumber(); iech_out++)
  {
    if (! getDbout()->isActive(iech_out)) continue;
    getGridout()->rankToIndice(iech_out, indgOut);
    _reduceIndices(indgOut, indgIn);
    int iech_in = getGridin()->indiceToRank(indgIn);
    double value = getDbin()->getLocVariable(ELoc::Z,iech_in, 0);
    getDbout()->setArray(iech_out, _iattOut, value);
  }
  return true;
}

bool CalcGridToGrid::_g2gShrink()
{
  int ndim_in  = getDbin()->getNDim();
  int ndim_out = getDbout()->getNDim();
  VectorInt indgIn(ndim_in);
  VectorInt indgOut(ndim_out);

  // Loop on the input file
  for (int iech_in = 0; iech_in < getDbin()->getActiveSampleNumber(); iech_in++)
  {
    if (! getDbin()->isActive(iech_in)) continue;
    getGridin()->rankToIndice(iech_in, indgIn);
    _reduceIndices(indgIn, indgOut);
    int iech_out = getGridout()->indiceToRank(indgOut);
    double value = getDbout()->getLocVariable(ELoc::Z,iech_out, 0);
    if (! FFFF(value))
    {
      value += getDbout()->getArray(iech_out, _iattOut);
      getDbout()->setArray(iech_out, _iattOut, value);
      double count = getDbout()->getArray(iech_out, _iattAux);
      getDbout()->setArray(iech_out, _iattAux, count + 1.);
    }
  }

  // Normalization
  for (int iech_out = 0; iech_out < getDbout()->getActiveSampleNumber(); iech_out++)
  {
    double ratio = getDbout()->getArray(iech_out, _iattAux);
    if (ratio > 0.)
    {
      double value = getDbout()->getArray(iech_out, _iattOut) / ratio;
      getDbout()->setArray(iech_out, _iattOut, value);
    }
    else
      getDbout()->setArray(iech_out, _iattOut, TEST);
  }
  return true;
}

bool CalcGridToGrid::_g2gInter()
{
  int ndim_in  = getDbin()->getNDim();
  int ndim_out = getDbout()->getNDim();
  int nvar = ndim_out - ndim_in;
  VectorInt indgIn(ndim_in);
  VectorInt indgOut(ndim_out);
  VectorDouble coorTop(nvar);
  VectorDouble coorBot(nvar);
  VectorDouble coorOut(ndim_out);

  // Identify the coordinate variables
  VectorInt iuidTop = getDbin()->getUIDs(_nameTops);
  VectorInt iuidBot = getDbin()->getUIDs(_nameBots);

  // Loop on the output file
  for (int iech_out = 0; iech_out < getDbout()->getActiveSampleNumber(); iech_out++)
  {
    if (! getDbout()->isActive(iech_out)) continue;
    getGridout()->rankToIndice(iech_out, indgOut);
    _reduceIndices(indgOut, indgIn);
    int iech_in = getGridin()->indiceToRank(indgIn);

    bool ret = true;
    ret = ret && getDbin()->isActive(iech_in);
    ret = ret &&_loadExtrema(nvar, iech_in, iuidTop, coorTop);
    ret = ret &&_loadExtrema(nvar, iech_in, iuidBot, coorBot);
    double value = TEST;
    if (ret)
    {
      getGridout()->rankToCoordinatesInPlace(iech_out, coorOut);
      double valTop = getDbin()->getLocVariable(ELoc::Z,iech_in, 0);
      double valBot = getDbin()->getLocVariable(ELoc::Z,iech_in, 1);
      value = _interpolate(nvar, valTop, valBot, coorTop, coorBot, coorOut);
    }
    getDbout()->setArray(iech_out, _iattOut, value);
  }
  return true;
}

/** Perform the linear interpolation
 *
 * @param nvar     Number of relevant coordinates
 * @param valTop   Value at the Top
 * @param valBot   Value at the Bottom
 * @param coorTop  Coordinates of the Top
 * @param coorBot  Coordinates of the Bottom
 * @param coorOut  Coordinates at Target (see remarks)
 * @return
 *
 * @remarks The vectors 'coorTop' and 'coorBot' are dimensioned to nvar
 * @remarks The vector 'coorOut' contains the ndimIn first coordinates followed
 * @remarks by the actual 'nvar' coordinates of the target
 */
double CalcGridToGrid::_interpolate(int nvar,
                                    double valTop,
                                    double valBot,
                                    const VectorDouble& coorTop,
                                    const VectorDouble& coorBot,
                                    const VectorDouble& coorOut)
{
  double delta;
  if (FFFF(valTop) || FFFF(valBot)) return TEST;

  int shift = getDbin()->getNDim();
  double dT0 = 0.;
  double dB0 = 0.;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    double zTop = coorTop[ivar];
    double zBot = coorBot[ivar];
    double zMin = MIN(zBot, zTop);
    double zMax = MAX(zBot, zTop);

    double zCur = coorOut[ivar + shift];
    if (zCur < zMin || zCur > zMax) return TEST;
    delta = zMax - zCur;
    dT0 += delta * delta;
    delta = zCur - zMin;
    dB0 += delta * delta;
  }

  double value = (dT0 * valTop + dB0 * valBot) / (dT0 + dB0);
  return value;
}

bool CalcGridToGrid::_loadExtrema(int nvar,
                                  int iech,
                                  const VectorInt &iuids,
                                  VectorDouble &coor)
{
  for (int ivar=0; ivar < nvar; ivar++)
  {
    coor[ivar] = getDbin()->getArray(iech, iuids[ivar]);
    if (FFFF(coor[ivar])) return false;
  }
  return true;
}
