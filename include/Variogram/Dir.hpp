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

class Dir : public AStringable
{
public:
  Dir(int ndim      = 2,
      int npas      = 0,
      double dpas   = 0.,
      double toldis = 0.5,
      double tolang = 90.);
  Dir(int ndim,
      int nvar,
      int npas,
      int flagAsym,
      double dpas,
      const VectorDouble& gg,
      const VectorDouble& hh = VectorDouble(),
      const VectorDouble& sw = VectorDouble());
  Dir(const Dir& r);
  Dir& operator=(const Dir& r);
  virtual ~Dir();

public:
  virtual String toString(int level = 0) const override;

  void init(int ndim,
            int npas,
            double dpas,
            double toldis,
            double tolang,
            int flag_asym = 0,
            int opt_code = 0,
            int idate = 0,
            double bench = TEST,
            double cylrad = TEST,
            double tolcode = 0,
            VectorDouble breaks = VectorDouble(),
            VectorDouble codir = VectorDouble(),
            VectorDouble grincr = VectorDouble());
  bool isCalculated() const;
  void copy(const Dir &dir);
  void clean(void);
  int  getAddress(int ivar, int jvar, int ipas, bool flag_abs, int sens) const;

  double getBench() const { return _bench; }
  const  VectorDouble& getBreaks() const { return _breaks; }
  const  double getBreaks(int i) const { return _breaks[i]; }
  const  VectorDouble& getCodir() const { return _codir; }
  double getCodir(int i) const { return _codir[i]; }
  double getCylRad() const { return _cylRad; }
  double getDPas() const { return _dPas; }
  double getLag() const { return _dPas; }
  int    getFlagAsym() const { return _flagAsym; }
  int    getIdate() const { return _idate; }
  int    getNPas() const { return _nPas; }
  int    getLagNumber() const { return _nPas; }
  int    getOptionCode() const { return _optionCode; }
  double getTolAngle() const { return _tolAngle; }
  double getTolCode() const { return _tolCode; }
  double getTolDist() const { return _tolDist; }
  const  VectorDouble& getUtilize() const { return _utilize; }
  double getUtilize(int i) const { return _utilize[i]; }
  int    getVariableNumber() const { return _nvar; }
  int    getDimensionNumber() const { return _ndim; }

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

  const  VectorDouble& getGrincr() const { return _grincr; }
  double getGrincr(int i) const;

  int  getBreakNumber() const { return ((int) _breaks.size() / 2); }
  bool getLagRegular() const { return (getBreakNumber() <= 0); }

  int  getLagTotalNumber() const;
  int  getSize() const;

  void setLagNumber(int npas) {_nPas = npas; }
  void setOptionCode(int option_code) {_optionCode = option_code; }
  void setIdate(int idate) {_idate = idate; }
  void setDPas(double dpas) {_dPas = dpas; }
  void setDLag(double dlag) {_dPas = dlag; }
  void setBench(double bench) {_bench = bench; }
  void setCylRad(double cylrad) {_cylRad = cylrad; }
  void setTolDist(double toldist) {_tolDist = toldist; }
  void setTolAngle(double tolang) {_tolAngle = tolang; }
  void setTolCode(double tolcode) {_tolCode = tolcode; }
  void setBreaks(VectorDouble breaks) {_breaks = breaks; }
  void setCodir(VectorDouble codir) {_codir = codir; }
  void setGrincr(VectorDouble grincr) {_grincr = grincr; }

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

  void resize(int nvar, int flagAsym);

  double getHmax(int ivar=0, int jvar=0) const;

private:
  void _completeDefinition();
  bool _isLagValid(int ilag) const;
  bool _isVariableValid(int ivar) const;
  bool _isDimensionValid(int idim) const;
  bool _isAddressValid(int iad) const;

private:
  int _ndim;
  int _nvar;
  int _flagAsym;
  int _nPas;
  int _optionCode;
  int _idate;
  double _dPas;
  double _bench;
  double _cylRad;
  double _tolDist;
  double _tolAngle;
  double _tolCode;
  VectorDouble _breaks;
  VectorDouble _codir;
  VectorDouble _grincr;
  VectorDouble _sw;      /* Array for number of lags */
  VectorDouble _gg;      /* Array for average variogram values */
  VectorDouble _hh;      /* Array for average distance values */
  VectorDouble _utilize; /* Array to mention if a lag is used or not */
};

std::vector<Dir> generateMultipleDirs(int ndim,
                                      int ndir,
                                      int npas = 0,
                                      double dpas = 0.,
                                      double toldis = 0.5);
