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

class Db;

class DirParam : public AStringable
{
public:
  DirParam(int ndim = 2,
           int npas = 0,
           double dpas = 0.,
           double toldis = 0.5,
           double tolang = 90.,
           int opt_code = 0,
           int idate = 0,
           double bench = TEST,
           double cylrad = TEST,
           double tolcode = 0,
           VectorDouble breaks = VectorDouble(),
           VectorDouble codir  = VectorDouble(),
           VectorInt grincr    = VectorInt());
  DirParam(int ndim, int npas, const VectorInt& grincr);
  DirParam(const DirParam& r);
  DirParam& operator=(const DirParam& r);
  virtual ~DirParam();

public:
  virtual String toString(int level = 0) const override;

  void init(int ndim,
            int npas,
            double dpas,
            double toldis,
            double tolang,
            int opt_code = 0,
            int idate = 0,
            double bench = TEST,
            double cylrad = TEST,
            double tolcode = 0,
            VectorDouble breaks = VectorDouble(),
            VectorDouble codir  = VectorDouble(),
            VectorInt grincr    = VectorInt());
  bool isCalculated() const;

  double getBench() const { return _bench; }
  const  VectorDouble& getBreaks() const { return _breaks; }
  const  double getBreaks(int i) const { return _breaks[i]; }
  const  VectorDouble& getCodir() const { return _codir; }
  double getCodir(int i) const { return _codir[i]; }
  double getCylRad() const { return _cylRad; }
  double getDPas() const { return _dPas; }
  double getLag() const { return _dPas; }
  int    getIdate() const { return _idate; }
  int    getLagNumber() const { return _nPas; }
  int    getOptionCode() const { return _optionCode; }
  double getTolAngle() const { return _tolAngle; }
  double getTolCode() const { return _tolCode; }
  double getTolDist() const { return _tolDist; }
  int    getDimensionNumber() const { return _ndim; }

  const  VectorInt& getGrincr() const { return _grincr; }
  int getGrincr(int i) const;
  double getMaximumDistance() const;

  int  getBreakNumber() const { return ((int) _breaks.size() / 2); }
  bool getFlagRegular() const { return (getBreakNumber() <= 0); }

  void setLagNumber(int npas) {_nPas = npas; }
  void setOptionCode(int option_code) {_optionCode = option_code; }
  void setIdate(int idate) {_idate = idate; }
  void setDPas(double dpas) {_dPas = dpas; }
  void setDLag(double dlag) {_dPas = dlag; }
  void setDPas(const Db* db);
  void setBench(double bench) {_bench = bench; }
  void setCylRad(double cylrad) {_cylRad = cylrad; }
  void setTolDist(double toldist) {_tolDist = toldist; }
  void setTolAngle(double tolang) {_tolAngle = tolang; }

  void setTolCode(double tolcode) {_tolCode = tolcode; }
  void setBreaks(VectorDouble breaks) {_breaks = breaks; }
  void setCodir(VectorDouble codir) {_codir = codir; }
  void setGrincr(VectorInt grincr) {_grincr = grincr; }

  bool isLagValid(int ilag) const;
  bool isDimensionValid(int idim) const;

  void setDimensionNumber(int ndim) { _ndim = ndim; }

private:
  void _completeDefinition();

private:
  int _ndim;
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
  VectorInt    _grincr;
};

std::vector<DirParam> generateMultipleDirs(int ndim,
                                            int ndir,
                                            int npas = 0,
                                            double dpas = 0.,
                                            double toldis = 0.5);
std::vector<DirParam> generateMultipleGridDirs(int ndim, int npas);
