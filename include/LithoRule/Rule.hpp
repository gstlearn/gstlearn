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
#include "Basic/IClonable.hpp"
#include "geoslib_enum.h"

class Db;
class Model;
class PropDef;

class Rule: public AStringable, public ASerializable
{
public:
  Rule(double rho = 0.);
  Rule(const VectorString& nodnames,double rho = 0.);
  Rule(const VectorInt& nodes,double rho = 0.);
  Rule(const VectorInt& n_type, const VectorInt& n_facs, double rho = 0.);
  Rule(int nfacies, double rho = 0.);
  Rule(const String& neutralFileName, bool verbose = false);

  Rule(const Rule& r);
  Rule& operator=(const Rule& r);
  virtual ~Rule();

  virtual std::string toString(int level = 0) const override;
  virtual int deSerialize(const String& filename, bool verbose = false) override;
  virtual int serialize(const String& filename, bool verbose = false) const override;
  virtual int deSerializeSpecific() { return 0; }
  virtual void serializeSpecific() const { return; }
  virtual String displaySpecific(int flagProp, int flagThresh) const;

  virtual int particularities(Db *db,
                              const Db *dbprop,
                              Model *model,
                              int flag_grid_check,
                              int flag_stat);
  virtual bool checkModel(const Model* model, int nvar = 0) const;
  virtual int gaus2facData(PropDef *propdef,
                           Db *dbin,
                           Db *dbout,
                           int *flag_used,
                           int ipgs,
                           int isimu,
                           int nbsimu);
  virtual int gaus2facResult(PropDef *propdef,
                             Db *dbout,
                             int *flag_used,
                             int ipgs,
                             int isimu,
                             int nbsimu);
  virtual int evaluateBounds(PropDef *propdef,
                             Db *dbin,
                             Db *dbout,
                             int isimu,
                             int igrf,
                             int ipgs,
                             int nbsimu);
  virtual double getShift(int idim) const { return TEST; }

  double getDMax() const { return _dMax; }
  int    getFlagProp() const { return _flagProp; }
  double getIncr() const { return _incr; }
  int    getModeRule() const { return _modeRule; }
  double getRho() const { return _rho; }
  double getTgte() const { return _tgte; }
  const Node*  getMainNode() const { return _mainNode; }

  void   setFlagProp(int flagProp) { _flagProp = flagProp; }
  void   setRho(double rho) { _rho = rho; }
  int    setProportions(const VectorDouble& proportions = VectorDouble());

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
  VectorDouble getThresh(int facies) const;
  VectorDouble getThreshFromRectangle(int rect, int *facies);
  int getFaciesFromGaussian(double y1, double y2) const;

  void updateShift();

  void setModeRule(int modeRule) { _modeRule = modeRule; }

protected:
  void setMainNodeFromNodNames(const VectorInt& n_type,
                               const VectorInt& n_facs);
  void setMainNodeFromNodNames(const VectorString& nodnames);
  int  setMainNodeFromNodNames(const VectorInt& nodes);
  int replicateInvalid(Db *dbin, Db *dbout, int jech);
  VectorString buildNodNames(int nfacies);

private:
  String _display(bool flagProp, bool flagThresh) const;
  void _ruleDefine(const Node *node,
                   int from_type,
                   int from_rank,
                   int from_vers,
                   int *rank) const;
  void _nodNamesToIds(const VectorString& nodes,
                      VectorInt &n_type,
                      VectorInt& n_facs);

private:
  int    _modeRule;    /* Type of usage */
  int    _flagProp;    /* 1 if proportions are defined; 0 otherwise */
  double _rho;         /* Correlation between GRFs */
  double _slope;       /* Slope used for shadow option */
  double _dMax;        /* Longest distance for shadow */
  double _tgte;        /* Tangent of the slope */
  double _incr;        /* Increments used for creating replicates */
  VectorDouble _shift; /* Shadow or translation orientation */
  mutable VectorDouble _xyz;
  VectorInt    _ind1;
  VectorInt    _ind2;
  Node*        _mainNode;
};

void   set_rule_mode(int rule_mode);
int    get_rule_mode(void);
double get_rule_extreme(int mode);
