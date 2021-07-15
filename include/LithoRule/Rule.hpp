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

#include "LithoRule/Node.hpp"
#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "geoslib_enum.h"


class Db;
class Model;

class Rule: public AStringable, ASerializable
{
public:
  Rule(int mode_rule = RULE_STD,
       double rho = 0.,
       double slope = 0.,
       double sh_down = 0.,
       double sh_dsup = 0.,
       const VectorDouble& shift = VectorDouble());
  // Constructor for standard option using the String definition
  Rule(const VectorString& nodnames,double rho = 0.);
  // Constructor for Standard option
  Rule(const VectorInt& n_type,
       const VectorInt& n_facs,
       double rho = 0.);
  // Constructor for Shift option
  Rule(const VectorDouble& shift);
  // Constructor for Shadow option
  Rule(double slope, double sh_dsup, double sh_down, const VectorDouble& shift);
  Rule(const String& neutralFileName, bool verbose);

  Rule(const Rule& r);
  Rule& operator=(const Rule& r);
  virtual ~Rule();

  virtual std::string toString(int level = 0) const override;
  int deSerialize(const String& filename, bool verbose = false) override;
  int serialize(const String& filename, bool verbose = false) override;

  double getDMax() const { return _dMax; }
  int    getFlagProp() const { return _flagProp; }
  double getIncr() const { return _incr; }
  int    getModeRule() const { return _modeRule; }
  double getRho() const { return _rho; }
  double getShDown() const { return _shDown; }
  double getShDsup() const { return _shDsup; }
  const VectorDouble& getShift() const { return _shift; }
  double getShift(int idim) const { return _shift[idim]; }
  double getSlope() const { return _slope; }
  double getTgte() const { return _tgte; }
  Node*  getMainNode() { return _mainNode; }

  void   setFlagProp(int flagProp) { _flagProp = flagProp; }
  void   setRho(double rho) { _rho = rho; }
  int    setProportions(const VectorDouble& proportions = VectorDouble());

  int init(const VectorInt& nodes);
  int statistics(int  verbose,
                 int *node_tot,
                 int *nfac_tot,
                 int *nmax_tot,
                 int *ny1_tot,
                 int *ny2_tot,
                 double *prop_tot) const;

  using AStringable::display; // https://stackoverflow.com/questions/18515183/c-overloaded-virtual-function-warning-by-clang
  void display(bool flagProp, bool flagThresh = false) const;

  int  getFaciesNumber() const;
  int  getGRFNumber() const;
  int  getY1Number() const;
  int  getY2Number() const;
  bool isYUsed(int igrf) const;
  VectorInt whichGRFUsed() const;
  double getProportion(int facies);
  VectorDouble getThresh(int facies);
  VectorDouble getThreshFromRectangle(int rect, int *facies);
  int getFaciesFromGaussian(double y1, double y2);
  int particularities(Db *db,
                      Db *dbprop,
                      Model *model,
                      int flag_grid_check,
                      int flag_stat);
  int particularities_shadow(Db *db,
                             Db *dbprop,
                             Model *model,
                             int flag_grid_check,
                             int flag_stat);
  double st_grid_eval(Db *dbgrid,
                      int isimu,
                      int icase,
                      int nbsimu,
                      VectorDouble& xyz0);
  void updateShift();

private:
  String _display(bool flagProp, bool flagThresh) const;
  void _st_shadow_max(Db *dbprop,
                      int flag_stat,
                      double *sh_dsup_max,
                      double *sh_down_max);
  int _st_shift_on_grid(Db *db, int ndim, int flag_grid_check);
  void _nodNamesToIds(const VectorString& nodes, VectorInt &n_type, VectorInt& n_facs);
  void _ruleDefine(Node *node,
                   int from_type,
                   int from_rank,
                   int from_vers,
                   int *rank);

private:
  int    _modeRule;    /* Type of usage */
  int    _flagProp;    /* 1 if proportions are defined; 0 otherwise */
  double _rho;         /* Correlation between GRFs */
  double _shDsup;      /* Upper limit */
  double _shDown;      /* Downwards limit */
  double _slope;       /* Slope used for shadow option */
  double _dMax;        /* Longest distance for shadow */
  double _tgte;        /* Tangent of the slope */
  double _incr;        /* Increments used for creating replicates */
  VectorDouble _shift; /* Shadow or translation orientation */
  VectorDouble _xyz;
  VectorInt    _ind1;
  VectorInt    _ind2;
  Node*        _mainNode;
};

void   set_rule_mode(int rule_mode);
int    get_rule_mode(void);
double get_rule_extreme(int mode);
