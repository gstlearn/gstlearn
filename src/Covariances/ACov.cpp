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
#include "Covariances/ACov.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Basic/AException.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"

ACov::ACov(const ASpace* space)
: ASpaceObject(space)
{
}

ACov::ACov(const ACov &r)
: ASpaceObject(r)
{
}

ACov& ACov::operator=(const ACov &r)
{
  if (this != &r)
  {
    ASpaceObject::operator=(r);
  }
  return *this;
}

ACov::~ACov()
{
}

MatrixSquareGeneral ACov::eval0Nvar(const CovCalcMode& mode) const
{
  int nvar = getNVariables();
  MatrixSquareGeneral mat(nvar);
  for (int ivar=0; ivar<nvar; ivar++)
    for (int jvar=0; jvar<nvar; jvar++)
      mat.setValue(ivar, jvar, eval0(ivar, jvar, mode)); // pure virtual method
  return mat;
}

VectorDouble ACov::eval(int ivar,
                        int jvar,
                        const std::vector<SpacePoint>& vec_p1,
                        const std::vector<SpacePoint>& vec_p2,
                        const CovCalcMode& mode) const
{
  VectorDouble vec;
  if (vec_p1.size() != vec_p2.size())
    my_throw ("Error: 'p1' and 'p2' should have same dimension");
  for (int i=0, n=static_cast<int> (vec_p1.size()); i < n; i++)
    vec.push_back(eval(ivar, jvar, vec_p1[i], vec_p2[i], mode)); // pure virtual method
  return vec;
}

MatrixSquareGeneral ACov::eval(const SpacePoint& p1,
                               const SpacePoint& p2,
                               const CovCalcMode& mode) const
{
  int nvar = getNVariables();
  MatrixSquareGeneral mat(nvar);
  for (int ivar=0; ivar<nvar; ivar++)
    for (int jvar=0; jvar<nvar; jvar++)
      mat.setValue(ivar, jvar, eval(ivar, jvar, p1, p2, mode)); // pure virtual method
  return mat;
}

/**
 * Covariance from a given point (center) in a given direction (dir *step)
 * @param ivar   Rank of the first variable
 * @param jvar   Rank of the second variable
 * @param step   Step value
 * @param dir    Direction definition
 * @param center Coordinates of the center
 * @param mode   CovCalcMode structure
 * @return
 */
double ACov::evalIvarIpas(int ivar,
                          int jvar,
                          double step,
                          const VectorDouble& dir,
                          const VectorDouble& center,
                          const CovCalcMode& mode) const
{
  // Define the point in the ACov space (center will be checked)
  SpacePoint p1(center,getSpace());
  SpacePoint p2(center,getSpace());
  VectorDouble dirloc(dir);
  if (dirloc.empty())
  {
    dirloc.resize(getNDim(),0.);
    dirloc[0] = 1.;
  }
  VectorDouble vec(dirloc);
  ut_vector_multiply_inplace(vec, step);
  p2.move(vec);
  return eval(ivar, jvar, p1, p2, mode); // pure virtual method
}

/**
 * Covariance vector from a given point (center) in a given direction (dir * steps)
 * for a pair of variables and a set of steps
 * @param ivar      Rank of the first variable
 * @param jvar      Rank of the second variable
 * @param vec_step  Vector of step values
 * @param dir       Direction definition
 * @param center    Coordinates of the Center
 * @param mode      CovCalcMode structure
 * @return
 */
VectorDouble ACov::evalIvarNpas(int ivar,
                                int jvar,
                                const VectorDouble& vec_step,
                                const VectorDouble& dir,
                                const VectorDouble& center,
                                const CovCalcMode& mode) const
{
  VectorDouble vec;
  for (int i=0, n=static_cast<int> (vec_step.size()); i < n; i++)
    vec.push_back(evalIvarIpas(ivar, jvar, vec_step[i], dir, center, mode));
  return vec;
}

/**
 * Covariance Matrix from a given point (center) in a given direction (dir * step)
 * for a set of variables and a given step
 * @param step   Step value
 * @param dir    Direction definition
 * @param center Coordinates of the center
 * @param mode   CovCalcMode structure
 * @return
 */
MatrixSquareGeneral ACov::evalNvarIpas(double step,
                                       const VectorDouble& dir,
                                       const VectorDouble& center,
                                       const CovCalcMode& mode) const
{
  int nvar = getNVariables();
  MatrixSquareGeneral mat(nvar);
  for (int ivar=0; ivar<nvar; ivar++)
    for (int jvar=0; jvar<nvar; jvar++)
      mat.setValue(ivar, jvar, evalIvarIpas(ivar, jvar, step, dir, center, mode));
  return mat;
}

/**
 * Covariance for a given unit global distance (without anisotropy)
 * for a pair of variables and a single step
 * @param ivar Rank of the first variable
 * @param jvar Rank of the second variable
 * @param step Step value
 * @param mode CovCalcMode structure
 * @return
 */
double ACov::evalIsoIvarIpas(int ivar,
                             int jvar,
                             double step,
                             const CovCalcMode& mode) const
{
  /// TODO : Not true whatever the space
  VectorDouble center = getOrigin();
  VectorDouble dir = getUnitaryVector();
  return evalIvarIpas(ivar, jvar, step, dir, center, mode);
}

/**
 * Covariance for a given unit global distance (without anisotropy)
 * for a pair of variables and a set of steps
 * @param ivar
 * @param jvar
 * @param vec_step
 * @param mode
 * @return
 */
VectorDouble ACov::evalIsoIvarNpas(int ivar,
                                   int jvar,
                                   const VectorDouble& vec_step,
                                   const CovCalcMode& mode) const
{
  VectorDouble vec;
  /// TODO : Not true whatever the space
  VectorDouble center = getOrigin();
  VectorDouble dir = getUnitaryVector();
  for (const auto& h : vec_step)
    vec.push_back(evalIvarIpas(ivar, jvar, h, dir, center, mode));
  return vec;
}

/**
 * Covariance for a given unit global distance (without anisotropy)
 * for a set of variables and a single step
 * @param step Step value
 * @param mode CovCalcMode structure
 * @return
 */
MatrixSquareGeneral ACov::evalIsoNvarIpas(double step,
                                          const CovCalcMode& mode) const
{
  int nvar = getNVariables();
  /// TODO : Not true whatever the space
  VectorDouble center = getOrigin();
  VectorDouble dir = getUnitaryVector();
  MatrixSquareGeneral mat(nvar);
  for (int ivar=0; ivar<nvar; ivar++)
    for (int jvar=0; jvar<nvar; jvar++)
      mat.setValue(ivar, jvar, evalIvarIpas(ivar, jvar, step, dir, center, mode));
  return mat;
}
