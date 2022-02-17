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
#include "LithoRule/Node.hpp"
#include "LithoRule/ERule.hpp"
#include "RuleStringFormat.hpp"

#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/IClonable.hpp"

class Db;
class Model;
class PropDef;

class GSTLEARN_EXPORT Rule: public AStringable, public ASerializable
{
public:
  Rule(double rho = 0.);
  Rule(const Rule& r);
  Rule& operator=(const Rule& r);
  virtual ~Rule();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int resetFromNames(const VectorString& nodnames,double rho = 0.);
  int resetFromCodes(const VectorInt& nodes,double rho = 0.);
  int resetFromNumericalCoding(const VectorInt& n_type, const VectorInt& n_facs, double rho = 0.);
  int resetFromFaciesCount(int nfacies, double rho = 0.);

  int dumpToNF(const String& neutralFilename, bool verbose = false) const;
  int dumpToNF2(const String& neutralFilename, bool verbose = false) const;
  static Rule* create(double rho = 0.);
  static Rule* createFromNF(const String& neutralFilename, bool verbose = false);
  static Rule* createFromNF2(const String& neutralFilename, bool verbose = false);
  static Rule* createFromNames(const VectorString& nodnames,double rho = 0.);
  static Rule* createFromCodes(const VectorInt& nodes,double rho = 0.);
  static Rule* createFromNumericalCoding(const VectorInt& n_type,
                                         const VectorInt& n_facs,
                                         double rho = 0.);
  static Rule* createFromFaciesCount(int nfacies, double rho = 0.);

  virtual String displaySpecific() const;

  virtual int particularities(Db *db,
                              const Db *dbprop,
                              Model *model,
                              int flag_grid_check,
                              int flag_stat) const;
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
                             int nbsimu) const;
  virtual int evaluateBounds(PropDef *propdef,
                             Db *dbin,
                             Db *dbout,
                             int isimu,
                             int igrf,
                             int ipgs,
                             int nbsimu) const;

  int          getFlagProp() const { return _flagProp; }
  const ERule& getModeRule() const { return _modeRule; }
  double       getRho()      const { return _rho; }
  const Node*  getMainNode() const { return _mainNode; }

  void setFlagProp(int flagProp)          { _flagProp = flagProp; }
  void setRho(double rho) const           { _rho = rho; } /// TODO : Check if mutable is really necessary
  void setModeRule(const ERule& modeRule) { _modeRule = modeRule; }

  int setProportions(const VectorDouble& proportions = VectorDouble()) const;

  int statistics(int  verbose,
                 int *node_tot,
                 int *nfac_tot,
                 int *nmax_tot,
                 int *ny1_tot,
                 int *ny2_tot,
                 double *prop_tot) const;

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

  void updateShift() const;

protected:
  void setMainNodeFromNodNames(const VectorInt& n_type,
                               const VectorInt& n_facs);
  void setMainNodeFromNodNames(const VectorString& nodnames);
  int  setMainNodeFromNodNames(const VectorInt& nodes);
  int replicateInvalid(Db *dbin, Db *dbout, int jech) const;
  VectorString buildNodNames(int nfacies);
  virtual int _deserializeSpecific(FILE* /*file*/) { return 0; }
  virtual void _serializeSpecific(FILE* /*file*/) const { return; }
  virtual int _deserialize(FILE* file, bool verbose = false) ;
  virtual int _serialize(FILE* file, bool verbose = false) const;

  virtual int _serialize2(std::ostream& os, bool verbose = false) const override;
  virtual int _deserializeSpecific2(std::istream& /*is*/) { return 0; }
  virtual int _deserialize2(std::istream& is, bool verbose = false) override;
  virtual void _serializeSpecific2(std::ostream& /*file*/) const { return; }

private:
  void _ruleDefine(FILE* file,
                   const Node *node,
                   int from_type,
                   int from_rank,
                   int from_vers,
                   int *rank) const;
  void _ruleDefine2(std::ostream& os,
                   const Node *node,
                   int from_type,
                   int from_rank,
                   int from_vers,
                   int *rank) const;
  void _nodNamesToIds(const VectorString& nodes,
                      VectorInt &n_type,
                      VectorInt& n_facs);
  void _clear();

private:
  ERule          _modeRule;  /* Type of usage (ERule) */
  mutable int    _flagProp;  /* 1 if proportions are defined; 0 otherwise */
  mutable double _rho;       /* Correlation between GRFs */
  Node*          _mainNode;
};

GSTLEARN_EXPORT void   set_rule_mode(int rule_mode);
GSTLEARN_EXPORT int    get_rule_mode(void);
GSTLEARN_EXPORT double get_rule_extreme(int mode);
