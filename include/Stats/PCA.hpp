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
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/Vector.hpp"

class GSTLEARN_EXPORT PCA
{
private:
  int          _nVar;
  VectorDouble _mean;
  VectorDouble _sigma;
  VectorDouble _eigen;
  VectorDouble _Z2F;
  VectorDouble _F2Z;

public:
  PCA();
  PCA(const PCA &m);
  PCA& operator= (const PCA &m);
  virtual ~PCA();

  void init(int nvar);
  void clean();
  int calculateEigen(int nvar, VectorDouble& c0);
  void display(int flag_center, int flag_stats);

  const VectorDouble& getEigen() const
  {
    return _eigen;
  }
  double getEigen(int ivar) const
  {
    return _eigen[ivar];
  }

  const VectorDouble& getMean() const
  {
    return _mean;
  }
  double getMean(int ivar) const
  {
    return _mean[ivar];
  }

  int getNVar() const
  {
    return _nVar;
  }

  const VectorDouble& getF2Z() const
  {
    return _F2Z;
  }
  double getF2Z(int ivar, int jvar) const
  {
    return _F2Z[_getAddress(ivar,jvar)];
  }

  const VectorDouble& getZ2F() const
  {
    return _Z2F;
  }
  double getZ2F(int ivar, int jvar) const
  {
    return _Z2F[_getAddress(ivar, jvar)];
  }

  const VectorDouble& getSigma() const
  {
    return _sigma;
  }
  double getSigma(int ivar) const
  {
    return _sigma[ivar];
  }

  void setPcaZ2F(VectorDouble& pcaz2f)
  {
    _setPcaZ2F(pcaz2f);
  }
  void setPcaZ2F(int ivar, int jvar, double pcaz2f)
  {
    _setPcaZ2F(ivar, jvar, pcaz2f);
  }

  void setPcaF2Z(VectorDouble& pcaf2z)
  {
    _setPcaF2Z(pcaf2z);
  }
  void setEigen(VectorDouble& eigen)
  {
    _setEigen(eigen);
  }
  void setEigen(int ivar, double eigen)
  {
    _setEigen(ivar,eigen);
  }

  void setMean(VectorDouble& mean) { _setMean(mean); }
  void setSigma(VectorDouble& sigma) { _setSigma(sigma); }

private:
  int _getAddress(int ivar, int jvar) const
  {
    return (ivar * _nVar + jvar);
  }
  void _setPcaZ2F(VectorDouble& pcaz2f)
  {
    _Z2F = pcaz2f;
  }
  void _setPcaZ2F(int ivar, int jvar, double pcaz2f)
  {
    _Z2F[_getAddress(ivar,jvar)] = pcaz2f;
  }

  void _setPcaF2Z(VectorDouble& pcaf2z)
  {
    _F2Z = pcaf2z;
  }

  void _setEigen(VectorDouble& eigen)
  {
    _eigen = eigen;
  }
  void _setEigen(int ivar, double eigen)
  {
    _eigen[ivar] = eigen;
  }

  void _setMean(VectorDouble& mean)
  {
    _mean = mean;
  }
  void _setSigma(VectorDouble& sigma)
  {
    _sigma = sigma;
  }
};
