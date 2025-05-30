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
#pragma once

#include "Space/ASpace.hpp"
#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include "Space/ASpaceObject.hpp"
#include "geoslib_define.h"

class Db;
class DbGrid;

/**
 * \brief
 * Class containing the definition of a **Direction** used for the calculation of the experimental Spatial Characteristics
 * as calculated experimentally from the data (contained in a Db).
 * This class corresponds to one item of the list of criteria stored in VarioParam and use for the calculation of Vario.
 *
 * The **Direction** consists in a series of rules (some are optional) for comparing two active samples of the Db:
 * - their distance must be assigned to a lag: i.e. the distance must correspond to a multiple of the lag (**lag**),
 * up to a tolerance (**tolerance**) expressed as a percentage of the lag. The rank of this multiple must be smaller
 * than the number of lags (**nlag**).
 * - the lag definition can be replaced by a series of intervals (**breaks**): the pair is selected if the distance
 * belongs to one of these intervals.
 * - the orientation of the segment joining the two points must be assigned to the current direction
 * characterized by its angle (expressed by its direction coefficients **codir**), up to a tolerance on angle (**tolangle**)
 * given in degrees.
 * - the distance between the two points (measured along the axis perpendicular to the direction)
 * must be smaller than a maximum cylinder distance (**cylrad**).
 * - the distance between the two points (measured along the highest space dimension) must be smaller than a bench height (**bench**)
 * - the difference between the code values (locator ELoc::CODE) defined at both samples must be either smaller or larger
 * than the tolerance on the code (**tolcode**).
 * - the two saples must share the same data (ELoc::DATE)
 *
 * In the case, the Db correspond to a grid, the lag is defined as an increment on the grid meshes (**grincr**)
 */
class GSTLEARN_EXPORT DirParam : public ASpaceObject
{
public:
  DirParam(int nlag = 10,
           double dlag = 1.,
           double toldis = 0.5,
           double tolang = 90.,
           int opt_code = 0,
           int idate = 0,
           double bench = TEST,
           double cylrad = TEST,
           double tolcode = 0.,
           const VectorDouble& breaks = VectorDouble(),
           const VectorDouble& codir  = VectorDouble(),
           double angle2D = TEST,
           const ASpaceSharedPtr& space = ASpaceSharedPtr());
  DirParam(const DbGrid *dbgrid,
           int nlag,
           const VectorInt &grincr,
           const ASpaceSharedPtr& space);
  DirParam(const DirParam& r);
  DirParam& operator=(const DirParam& r);
  virtual ~DirParam();

  static DirParam* create(int nlag = 10,
                          double dlag = 1.,
                          double toldis = 0.5,
                          double tolang = 90.,
                          int opt_code = 0,
                          int idate = 0,
                          double bench = TEST,
                          double cylrad = TEST,
                          double tolcode = 0.,
                          const VectorDouble& breaks = VectorDouble(),
                          const VectorDouble& codir = VectorDouble(),
                          double angle2D = TEST,
                          const ASpaceSharedPtr& space = ASpaceSharedPtr());
  static DirParam* createOmniDirection(int nlag = 10,
                                       double dlag = 1.,
                                       double toldis = 0.5,
                                       int opt_code = 0,
                                       int idate = 0,
                                       double bench = TEST,
                                       double cylrad = TEST,
                                       double tolcode = 0.,
                                       const VectorDouble& breaks = VectorDouble(),
                                       const ASpaceSharedPtr& space = ASpaceSharedPtr());
  static DirParam* createFromGrid(const DbGrid* dbgrid,
                                  int nlag = 10,
                                  const VectorInt& grincr = VectorInt(),
                                  const ASpaceSharedPtr& space = ASpaceSharedPtr());
  static std::vector<DirParam> createMultiple(int ndir,
                                              int nlag = 10,
                                              double dlag = 1.,
                                              double toldis = 0.5,
                                              double angref = 0.,
                                              const ASpaceSharedPtr& space = ASpaceSharedPtr());
  static std::vector<DirParam> createSeveral2D(const VectorDouble& angles,
                                               int nlag = 10,
                                               double dlag = 1.,
                                               double toldis = 0.5,
                                               double tolang = TEST,
                                               const ASpaceSharedPtr& space = ASpaceSharedPtr());
  static std::vector<DirParam> createMultipleInSpace(int nlag,
                                                     double dlag = 1.,
                                                     const ASpaceSharedPtr& space = ASpaceSharedPtr());

public:
  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ASpaceObject Interface
  virtual bool isConsistent(const ASpace* space) const override;

  double getBench() const { return _bench; }
  const  VectorDouble& getBreaks() const { return _breaks; }
  double getBreak(int i) const;
  const  VectorDouble& getCodirs() const { return _codir; }
  double getCodir(int i) const;
  double getCylRad() const { return _cylRad; }
  double getDPas() const { return _dLag; }
  double getLag() const { return _dLag; }
  int    getIdate() const { return _idate; }
  int    getNLag() const { return _nLag; }
  int    getOptionCode() const { return _optionCode; }
  double getTolAngle() const { return _tolAngle; }
  double getTolCode() const { return _tolCode; }
  double getTolDist() const { return _tolDist; }

  VectorInt getGrincrs() const { return _grincr; }
  int    getGrincr(int i) const;
  double getMaximumDistance() const;

  int  getNBreak() const { return ((int) _breaks.size() / 2); }
  bool getFlagRegular() const { return (getNBreak() <= 0); }

  void setNLag(int nlag) {_nLag = nlag; }
  void setOptionCode(int option_code) {_optionCode = option_code; }
  void setIdate(int idate) {_idate = idate; }
  void setDPas(double dlag) {_dLag = dlag; }
  void setDLag(double dlag) {_dLag = dlag; }
  void setDPas(const DbGrid* db);
  void setBench(double bench) {_bench = bench; }
  void setCylRad(double cylrad) {_cylRad = cylrad; }
  void setTolDist(double toldist) {_tolDist = toldist; }
  void setTolAngle(double tolang);

  void setTolCode(double tolcode) {_tolCode = tolcode; }
  void setBreaks(const VectorDouble& breaks) {_breaks = breaks; }
  void setCodir(const VectorDouble& codir) {_codir = codir; }
  void setGrincr(const VectorInt &grincr) { _grincr = grincr; }

  bool isLagValid(int ilag, bool flagAsym = false, bool flagCheck = true) const;
  bool isDimensionValid(int idim) const;
  bool isDefinedForGrid() const { return ! _grincr.empty(); }

  int getLagRank(double dist) const;

private:
  void _completeDefinition(double angle2D = TEST);

private:
  int    _nLag;
  int    _optionCode;
  int    _idate;
  double _dLag;
  double _bench;
  double _cylRad;
  double _tolDist;
  double _tolAngle;
  double _tolCode;
  VectorDouble _breaks;
  VectorDouble _codir;
  VectorInt    _grincr;
};

