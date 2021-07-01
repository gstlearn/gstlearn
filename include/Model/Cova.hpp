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

#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"
class Cova;

typedef struct
{
  int ndim;
  int rank_db1;
  int rank_db2;
  int rank_ech1;
  int rank_ech2;
  VectorDouble x1;
  VectorDouble x2;
} External_Cov;

/* Prototyping the internal covariance functions */
double _st_cov_func(Cova *, int, double, VectorDouble, double);

class Cova: public AStringable
{
private:
  int _type;
  int _nVar;
  int _nDim;
  int _flagFilter;
  int _flagAniso;
  int _flagRotation;
  int _covChoice;
  double _range;  // Practical range; Theoretical = range / _scafac
  double _param;
  VectorDouble _sill; // TODO should be a square symmetric matrix (nvar * nvar)
  VectorDouble _anisoAngles;
  VectorDouble _anisoCoeffs;
  VectorDouble _anisoRotMat; // TODO Should become a square matrix (ndim * ndim)

  double (*_st_cov_external)(int, int, int, int, int,
                             double, double *, double *, double *);

public:
  Cova();
  Cova(const Cova &m);
  Cova& operator=(const Cova &m);
  virtual ~Cova();

  void init(int ndim, int nvar);
  void init(int ndim,
            int nvar,
            int type,
            double range,
            double param,
            int flag_aniso = 0,
            int flag_rotation = 0,
            const VectorDouble& aniso_ranges = VectorDouble(),
            const VectorDouble& aniso_rotmat = VectorDouble(),
            const VectorDouble& sill = VectorDouble());

  int hasExternalCov() const
  {
    return (_st_cov_external != NULL);
  }
  double evaluateExternalCov(External_Cov& E_Cov, double h, VectorDouble d)
  {
    return _st_cov_external(E_Cov.ndim, E_Cov.rank_db1, E_Cov.rank_ech1,
                            E_Cov.rank_db2, E_Cov.rank_ech2, h, d.data(),
                            E_Cov.x1.data(), E_Cov.x2.data());
  }

  /****************************************************************************/
  /*!
   **  Define the external covariance function
   **
   ** \param[in]  cov_func   External Covariance function
   **
   ** \remarks  This external function has a given prototype
   **
   *****************************************************************************/
  void setCovExternal(double (*cov_func)(int,
                                         int,
                                         int,
                                         int,
                                         int,
                                         double,
                                         double *,
                                         double *,
                                         double *))
  {
    _st_cov_external = cov_func;
  }
};
