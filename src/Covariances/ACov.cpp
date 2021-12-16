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
#include "Basic/AException.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
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

bool ACov::isGradientCompatible() const
{
  messerr("This covariance is not compatible with Gradient calculations");
  return false;
}

MatrixSquareGeneral ACov::eval0(const CovCalcMode& mode) const
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

double ACov::eval(int ivar,
                  int jvar,
                  double step,
                  const VectorDouble& dir,
                  const VectorDouble& center,
                  const CovCalcMode& mode) const
{
  // Define the point in the ACov space (center will be checked)
  SpacePoint p1(center,getSpace());
  SpacePoint p2(center,getSpace());
  if (! dir.empty())
  {
    VectorDouble vec(dir);
    ut_vector_multiply_inplace(vec, step);
    p2.move(vec);
  }
  return eval(ivar, jvar, p1, p2, mode); // pure virtual method
}

VectorDouble ACov::eval(int ivar,
                        int jvar,
                        const VectorDouble& vec_step,
                        const VectorDouble& dir,
                        const VectorDouble& center,
                        const CovCalcMode& mode) const
{
  VectorDouble vec;
  for (int i=0, n=static_cast<int> (vec_step.size()); i < n; i++)
    vec.push_back(eval(ivar, jvar, vec_step[i], dir, center, mode));
  return vec;
}

MatrixSquareGeneral ACov::eval(double step,
                               const VectorDouble& dir,
                               const VectorDouble& center,
                               const CovCalcMode& mode) const
{
  int nvar = getNVariables();
  MatrixSquareGeneral mat(nvar);
  for (int ivar=0; ivar<nvar; ivar++)
    for (int jvar=0; jvar<nvar; jvar++)
      mat.setValue(ivar, jvar, eval(ivar, jvar, step, dir, center, mode));
  return mat;
}

double ACov::eval(int ivar,
                  int jvar,
                  double step,
                  const CovCalcMode& mode) const
{
  /// TODO : Not true whatever the space
  VectorDouble center = getOrigin();
  VectorDouble dir = getUnitaryVector();
  return eval(ivar, jvar, step, dir, center, mode);
}

VectorDouble ACov::eval(int ivar,
                        int jvar,
                        const VectorDouble& vec_step,
                        const CovCalcMode& mode) const
{
  VectorDouble vec;
  /// TODO : Not true whatever the space
  VectorDouble center = getOrigin();
  VectorDouble dir = getUnitaryVector();
  for (const auto& h : vec_step)
    vec.push_back(eval(ivar, jvar, h, dir, center, mode));
  return vec;
}

MatrixSquareGeneral ACov::eval(double step,
                       const CovCalcMode& mode) const
{
  int nvar = getNVariables();
  /// TODO : Not true whatever the space
  VectorDouble center = getOrigin();
  VectorDouble dir = getUnitaryVector();
  MatrixSquareGeneral mat(nvar);
  for (int ivar=0; ivar<nvar; ivar++)
    for (int jvar=0; jvar<nvar; jvar++)
      mat.setValue(ivar, jvar, eval(ivar, jvar, step, dir, center, mode));
  return mat;
}
