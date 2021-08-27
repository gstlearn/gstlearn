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
#include "Variogram/DirParam.hpp"

class Db;

class DirC : public AStringable
{
public:
  DirC(int nvar, const DirParam dirparam);
  DirC(int nvar, const DirParam dirparam,
       const VectorDouble& gg,
       const VectorDouble& hh = VectorDouble(),
       const VectorDouble& sw = VectorDouble());
  DirC(const DirC& r);
  DirC& operator=(const DirC& r);
  virtual ~DirC();

public:
  virtual String toString(int level = 0) const override;

  bool isCalculated() const;
  void copy(const DirC &DirC);
  void clean(void);
  int  getAddress(int ivar, int jvar, int ipas, bool flag_abs, int sens) const;

  const  VectorDouble& getUtilize() const { return _utilize; }
  double getUtilize(int i) const { return _utilize[i]; }
  int    getVariableNumber() const { return _nvar; }

  const  VectorDouble& getGg() const { return _gg; }
  double getGg(int iad) const;
  double getGg(int ivar, int jvar, int ipas) const;
  VectorDouble getGg(int ivar, int jvar) const;

  const  VectorDouble& getHh() const { return _hh; }
  double getHh(int iad) const;
  double getHh(int ivar, int jvar, int ipas) const;
  VectorDouble getHh(int ivar, int jvar) const;

  const  VectorDouble& getSw() const { return _sw; }
  double getSw(int iad) const;
  double getSw(int ivar, int jvar, int ipas) const;
  VectorDouble getSw(int ivar, int jvar) const;

  int getCenter(int ivar, int jvar) const;

  int  getSize() const;

  void setSw(int iad, double sw);
  void setHh(int iad, double hh);
  void setGg(int iad, double gg);
  void updSw(int iad, double sw);
  void updHh(int iad, double hh);
  void updGg(int iad, double gg);

  void setSw(int ivar, int jvar, int ipas, double sw);
  void setHh(int ivar, int jvar, int ipas, double hh);
  void setGg(int ivar, int jvar, int ipas, double gg);
  void setUtilize(int iad, double val);

  double getHmax(int ivar=0, int jvar=0) const;
  double getGmax(int ivar=0, int jvar=0, bool flagAbs = false) const;

  void patchCenter(int nech, double rho);

  int fill(int nvar,
           const VectorDouble& sw,
           const VectorDouble& gg,
           const VectorDouble& hh);

  const DirParam& getDirParam() const { return _dirparam; }

private:
  void _internalResize(int nvar);
  void _completeDefinition();
  bool _isVariableValid(int ivar) const;
  bool _isAddressValid(int iad) const;
  bool _isLagValid(int ilag) const;
  bool _isDimensionValid(int idim) const;

private:
  DirParam _dirparam;
  int _nvar;
  VectorDouble _sw;      /* Array for number of lags */
  VectorDouble _gg;      /* Array for average variogram values */
  VectorDouble _hh;      /* Array for average distance values */
  VectorDouble _utilize; /* Array to mention if a lag is used or not */
};
