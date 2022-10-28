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
#include "Basic/AStringable.hpp"
#include "geoslib_define.h"

class Db;
class DbGrid;

/**
 * Experimental Variogram calculation direction parameters TODO : to be improved
 */
class GSTLEARN_EXPORT DirParam : public AStringable // TODO : Should inherit from ASpaceObject
{
public:
  DirParam(int ndim = 2,
           int npas = 10,
           double dpas = 1.,
           double toldis = 0.5,
           double tolang = 90.,
           int opt_code = 0,
           int idate = 0,
           double bench = TEST,
           double cylrad = TEST,
           double tolcode = 0.,
           const VectorDouble& breaks = VectorDouble(),
           const VectorDouble& codir  = VectorDouble(),
           const VectorInt& grincr    = VectorInt());
  DirParam(const DirParam& r);
  DirParam& operator=(const DirParam& r);
  virtual ~DirParam();

  static DirParam* create(int ndim = 2,
                          int npas = 10,
                          double dpas = 1.,
                          double toldis = 0.5,
                          double tolang = 90.,
                          int opt_code = 0,
                          int idate = 0,
                          double bench = TEST,
                          double cylrad = TEST,
                          double tolcode = 0.,
                          const VectorDouble& breaks = VectorDouble(),
                          const VectorDouble& codir = VectorDouble());
  static DirParam* createOmniDirection(int ndim = 2,
                                       int npas = 10,
                                       double dpas = 1.,
                                       double toldis = 0.5,
                                       int opt_code = 0,
                                       int idate = 0,
                                       double bench = TEST,
                                       double cylrad = TEST,
                                       double tolcode = 0.,
                                       const VectorDouble& breaks = VectorDouble());
  static DirParam* createFromGrid(int ndim = 2,
                                  int npas = 10,
                                  const VectorInt& grincr = VectorInt());
  static std::vector<DirParam> createMultiple(int ndim,
                                              int ndir,
                                              int npas = 10,
                                              double dpas = 1.,
                                              double toldis = 0.5);
  static std::vector<DirParam> createMultipleFromGrid(int ndim, int npas);

public:
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;


  double getBench() const { return _bench; }
  const  VectorDouble& getBreaks() const { return _breaks; }
  double getBreak(int i) const;
  const  VectorDouble& getCodir() const { return _codir; }
  double getCodir(int i) const;
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
  void setDPas(const DbGrid* db);
  void setBench(double bench) {_bench = bench; }
  void setCylRad(double cylrad) {_cylRad = cylrad; }
  void setTolDist(double toldist) {_tolDist = toldist; }
  void setTolAngle(double tolang) {_tolAngle = tolang; }

  void setTolCode(double tolcode) {_tolCode = tolcode; }
  void setBreaks(VectorDouble breaks) {_breaks = breaks; }
  void setCodir(VectorDouble codir) {_codir = codir; }
  void setGrincr(VectorInt grincr) {_grincr = grincr; }

  bool isLagValid(int ilag, bool flagAsym = false) const;
  bool isDimensionValid(int idim) const;
  bool isDefinedForGrid() const { return _definedForGrid; }

private:
  void _completeDefinition();

private:
  int    _ndim;  // TODO : Should be stored by ASpaceObject upper class
  int    _nPas;
  int    _optionCode;
  int    _idate;
  bool   _definedForGrid;
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

