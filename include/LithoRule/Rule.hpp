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

#include "gstlearn_export.hpp"

#include "Enum/ERule.hpp"

#include "LithoRule/Node.hpp"
#include "LithoRule/PropDef.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"

namespace gstlrn
{
class Db;
class Model;
class PropDef;

class GSTLEARN_EXPORT Rule: public AStringable, public ASerializable
{
public:
  Rule(double rho = 0.);
  Rule(const Rule& m);
  Rule& operator=(const Rule& m);
  virtual ~Rule();

  String toString(const AStringFormat* strfmt = nullptr) const override;

  Id resetFromNames(const VectorString& nodnames,double rho = 0.);
  Id resetFromCodes(const VectorInt& nodes,double rho = 0.);
  Id resetFromNumericalCoding(const VectorInt& n_type, const VectorInt& n_facs, double rho = 0.);
  Id resetFromFaciesCount(Id nfacies, double rho = 0.);

  static Rule* create(double rho = 0.);
  static Rule* createFromNF(const String& NFFilename, bool verbose = true);
  static Rule* createFromNames(const VectorString& nodnames, double rho = 0.);
  static Rule* createFromCodes(const VectorInt& nodes, double rho = 0.);
  static Rule* createFromNumericalCoding(const VectorInt& n_type,
                                         const VectorInt& n_facs,
                                         double rho = 0.);
  static Rule* createFromFaciesCount(Id nfacies, double rho = 0.);

  virtual String displaySpecific() const;

  virtual Id particularities(Db *db,
                              const Db *dbprop,
                              Model *model,
                              Id flag_grid_check,
                              Id flag_stat) const;
  virtual bool checkModel(const Model* model, Id nvar = 0) const;
  virtual Id gaus2facData(PropDef *propdef,
                           Db *dbin,
                           Db *dbout,
                           Id *flag_used,
                           Id ipgs,
                           Id isimu,
                           Id nbsimu);
  virtual Id gaus2facResult(PropDef *propdef,
                             Db *dbout,
                             Id *flag_used,
                             Id ipgs,
                             Id isimu,
                             Id nbsimu) const;
  virtual Id evaluateBounds(PropDef *propdef,
                             Db *dbin,
                             Db *dbout,
                             Id isimu,
                             Id igrf,
                             Id ipgs,
                             Id nbsimu) const;

  Id          getFlagProp() const { return _flagProp; }
  const ERule& getModeRule() const { return _modeRule; }
  double       getRho()      const { return _rho; }
  const Node*  getMainNode() const { return _mainNode; }

  void setFlagProp(Id flagProp)          { _flagProp = flagProp; }
  void setRho(double rho) const           { _rho = rho; } /// TODO : Check if mutable is really necessary
  void setModeRule(const ERule& modeRule) { _modeRule = modeRule; }

  Id setProportions(const VectorDouble& proportions = VectorDouble()) const;

  Id statistics(Id  verbose,
                 Id *node_tot,
                 Id *nfac_tot,
                 Id *nmax_tot,
                 Id *ny1_tot,
                 Id *ny2_tot,
                 double *prop_tot) const;

  Id  getNFacies() const;
  Id  getNGRF() const;
  Id  getNY1() const;
  Id  getNY2() const;
  bool isYUsed(Id igrf) const;
  VectorInt whichGRFUsed() const;
  double getProportion(Id facies);
#ifndef SWIG
  std::array<double, 4> getThresh(Id facies) const;
#endif
  VectorDouble getThreshFromRectangle(Id rect, Id *facies);
  Id getFaciesFromGaussian(double y1, double y2) const;
  VectorInt getNodes() const;

  void updateShift() const;

protected:
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "Rule"; } 

  void setMainNodeFromNodNames(const VectorInt& n_type,
                               const VectorInt& n_facs);
  void setMainNodeFromNodNames(const VectorString& nodnames);
  Id  setMainNodeFromNodNames(const VectorInt& nodes);
  static Id replicateInvalid(Db *dbin, Db *dbout, Id jech);
  static VectorString buildNodNames(Id nfacies);

private:
  void _ruleDefine(std::ostream& os,
                   const Node *node,
                   Id from_type,
                   Id from_rank,
                   Id from_vers,
                   Id *rank) const;
  static void _nodNamesToIds(const VectorString& nodes,
                             VectorInt& n_type,
                             VectorInt& n_facs);
  void _clear();

private:
  ERule          _modeRule;  /* Type of usage (ERule) */
  mutable Id    _flagProp;  /* 1 if proportions are defined; 0 otherwise */
  mutable double _rho;       /* Correlation between GRFs */
  Node*          _mainNode;

  mutable VectorInt _facies;
  mutable VectorDouble _props;
};

GSTLEARN_EXPORT void   set_rule_mode(Id rule_mode);
GSTLEARN_EXPORT Id    get_rule_mode(void);
GSTLEARN_EXPORT double get_rule_extreme(Id mode);
GSTLEARN_EXPORT Rule* rule_free(const Rule* rule);
GSTLEARN_EXPORT Model* model_rule_combine(const Model* model1, const Model* model2, const Rule* rule);
GSTLEARN_EXPORT Id db_rule_shadow(Db* db,
                                   Db* dbprop,
                                   RuleShadow* rule,
                                   Model* model1,
                                   const VectorDouble& props,
                                   Id flag_stat,
                                   Id nfacies);
}