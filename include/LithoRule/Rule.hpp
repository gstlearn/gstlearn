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

#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/VectorNumT.hpp"
#include "LithoRule/Node.hpp"
#include "LithoRule/PropDef.hpp"

namespace gstlrn
{
  class Db;
  class Model;
  class PropDef;
  class RuleShadow;

  /**
   * @brief Lithological Rule which describes the arrangement of facies in the MonoGaussian
   * or the biGaussian space.
   *
   * Let us first describe the structure of the Rule in the MonoGaussian case.
   * The Rule is defined by a set of codes that constitue a grammar, explained hereafter.
   * The idea is to consider an axis (representing the Gaussian scale) and to define thresholds
   * along this axis.
   *
   * Each threshold is coded using "S" symbol. It is followed by the description of the Rule contents
   * on the left side of the threshold first, then on its right side.
   * On any side, we can imagine to define another threshold (using the code "S" again) or a Facies
   * (using the code "F" followed by the Facies rank).
   *
   * In the following example, we consider the following string (set of codes):
   * "S1", "S2", "F4", "F1", "S3", "F3", "F2"
   * (Note that the symbol "S" has been 'decorated' by an optional rank).
   * This corresponds to the following arrangement:
   * - Define Threshold S1 first
   * - On the left of S1, define Threshold S2
   * - On the left of S2, define Facies F4
   * - On the right of S2 (and left of S1), define Facies F1
   * - On the right of S1, define Threshold S3
   * - On the left of S3 (and the right of S1), define Facies F3
   * - On the right of S3, define Facies F2
   * This is one way to define following arrangement of Facies and Thresholds along the Gaussian axis:
   * F4 | S2 | F1 | S1 | F3 | S3 | F2
   * Note that the same result would be obtained with the following string:
   * "S1", "S3", "F2", "F3", "S2", "F1", "F4"
   *
   * In the BiGaussian case, the Rule is defined with a string (set of codes) that is similar to the one
   * used in the MonoGaussian case, but with the addition of a new symbol "T" that is used to define
   * thresholds along the second Gaussian axis.
   * When a "T" threshold is used, it must be followed by the description of the Rule contents
   * on the left side of the threshold first, then on its right side.
   *
   * Note that a convenient representation of the Rule is to consider it as a Square with
   * the first Gaussian axis as the horizontal axis and the second Gaussian axis as the vertical axis.
   * Due to this representation, we must switch from "left and right" to "below and above" when
   * describing the Rule contents on the two sides of a "T" threshold.
   *
   * The following simple BiGaussian example is considered:
   * "S1", "T1", "F4", "F1", "T2", "F3", "F2"
   * This corresponds to the following arrangement:
   * - Define Threshold S1 first (vertical split of the Square)
   * - On the left of S1, define Threshold T1 (horizontal split of the left part of the Square)
   * - On the bottom of T1 (and left of S1), define Facies F4
   * - On the top of T1 (and left of S1), define Facies F1
   * - On the right of S1, define Threshold T2 (horizontal split of the right part of the Square)
   * - On the bottom of T2 (and right of S1), define Facies F3
   * - On the top of T2 (and right of S1), define Facies F2
   * This is one way to define following arrangement of Facies and Thresholds along the two Gaussian axes:
   * F1 | S1 | F3
   * ---+----+---
   * T1 |    | T2
   * ---+----+---
   * F4 | S1 | F2
   * ---+----+---
   *
   */
  class GSTLEARN_EXPORT Rule: public AStringable,
                              public ASerializable,
                              public ICloneable
  {
  public:
    Rule(double rho = 0.);
#ifndef SWIG
    Rule(const VectorInt& icodes);
#endif
    Rule(const Rule& m);
    Rule& operator=(const Rule& m);
    virtual ~Rule();

    /// ICloneable Interface
    IMPLEMENT_CLONING(Rule)

    /// AStringable Interface
    String toString(const AStringFormat* strfmt = nullptr) const override;

    /// ASerializable Interface
    String getNFName() const override { return "Rule"; }
#ifdef HDF5
    bool deserializeH5(H5::Group& grp) override;
    bool serializeH5(H5::Group& grp) const override;
#endif

    static Rule* create(double rho = 0.);
    static Rule* createFromNF(const String& NFFilename, bool verbose = true);
    static Rule* createFromNames(const VectorString& nodnames, double rho = 0.);
    static Rule* createFromFaciesCount(Id nfacies, double rho = 0.);

    Id resetFromNames(const VectorString& nodnames, double rho = 0.);
    Id resetFromFaciesCount(Id nfacies, double rho = 0.);

    virtual String displaySpecific() const;

    virtual Id particularities(
      Db* db,
      const Db* dbprop,
      Model* model,
      Id flag_grid_check,
      Id flag_stat) const;
    virtual bool checkModel(const Model* model, Id nvar = 0) const;
    virtual Id gaus2facData(
      const PropDef& propdef,
      Db* dbin,
      Db* dbout,
      const VectorBool& flag_used,
      Id ipgs,
      Id isimu,
      Id nbsimu) const;
    virtual Id gaus2facResult(
      const PropDef& propdef,
      Db* dbout,
      const VectorBool& flag_used,
      Id ipgs,
      Id isimu,
      Id nbsimu) const;
    virtual Id evaluateBounds(
      const PropDef& propdef,
      Db* dbin,
      Db* dbout,
      Id isimu,
      Id igrf,
      Id ipgs,
      Id nbsimu) const;

    Id getFlagProp() const { return _flagProp; }

    const ERule& getModeRule() const { return _modeRule; }

    double getRho() const { return _rho; }

    const Node* getMainNode() const { return _mainNode; }

    void setFlagProp(Id flagProp) { _flagProp = flagProp; }

    void setRho(double rho) const
    {
      _rho = rho;
    } /// TODO : Check if mutable is really necessary

    void setModeRule(const ERule& modeRule) { _modeRule = modeRule; }

    Id setProportions(
      const VectorDouble& proportions = VectorDouble(),
      bool flagGaussian = true) const;

    Id statistics(
      Id verbose,
      Id* node_tot,
      Id* nfac_tot,
      Id* nmax_tot,
      Id* ny1_tot,
      Id* ny2_tot,
      double* prop_tot) const;

    Id getNFacies() const;
    Id getNGRF() const;
    Id getNY1() const;
    Id getNY2() const;
    bool isYUsed(Id igrf) const;
    VectorBool whichGRFUsed() const;
    double getProportion(Id facies);
#ifndef SWIG
    std::array<double, 4> getThresh(Id facies) const;
#endif
    VectorDouble getThreshFromRectangle(Id rect, Id* facies);
    Id getFaciesFromGaussian(double y1, double y2) const;
    VectorInt getNodes() const;

    void updateShift() const;

    VectorString getFaciesNames() const { return _facnames; }

    VectorInt getFaciesColors() const { return _faccols; }

    VectorInt getFaciesValues() const { return _facvalues; }

    String getFaciesName(Id facies) const;
    Id getFaciesColor(Id facies) const;
    Id getFaciesValue(Id facies) const;
    String getFaciesColorName(Id facies) const;

    double getScore() const { return _score; }

    void setFaciesNames(const VectorString& facnames) { _facnames = facnames; }

    void setFaciesColors(const VectorInt& faccols) { _faccols = faccols; }

    void setFaciesValues(const VectorInt& facvalues) { _facvalues = facvalues; }

    void setFaciesName(Id facies, const String& name);
    void setFaciesValue(Id facies, Id value);
    void setFaciesColorByName(Id facies, const String& color);
    void setFaciesColorByHexa(Id facies, const String& hexa);
    void setFaciesColorByInt(Id facies, const Id& value);
    void setCharacteristics(
      Id facies,
      const String& name = String(),
      const String& color = String(),
      Id value = ITEST);

    void setScore(double score) { _score = score; }

    String namesPrint() const;

  protected:
    Rule(const VectorInt& n_type, const VectorInt& n_facs, double rho = 0.);
    Id _resetFromNumericalCoding(
      const VectorInt& n_type,
      const VectorInt& n_facs,
      double rho = 0.);
    bool _serializeAscii(std::ostream& os) const override;
    bool _deserializeAscii(std::istream& is) override;

    void
      setMainNodeFromNodNames(const VectorInt& n_type, const VectorInt& n_facs);
    void setMainNodeFromNodNames(const VectorString& nodnames);
    Id setMainNodeFromNodNames(const VectorInt& nodes);
    static Id replicateInvalid(Db* dbin, Db* dbout, Id jech);
    static VectorString buildNodNames(Id nfacies);

  private:
    void _ruleDefine(
      std::ostream& os,
      const Node* node,
      Id from_type,
      Id from_rank,
      Id from_vers,
      Id* rank) const;
    static void _nodNamesToIds(
      const VectorString& nodes,
      VectorInt& n_type,
      VectorInt& n_facs);
    void _clear();
    void _initCharacteristics();

  private:
    ERule _modeRule; /* Type of usage (ERule) */
    mutable Id _flagProp; /* 1 if proportions are defined; 0 otherwise */
    mutable double _rho; /* Correlation between GRFs */
    Node* _mainNode;

    mutable double _score;
    mutable VectorInt _facies;
    mutable VectorDouble _props; // Constant proportion per facies
    mutable VectorString _facnames; // Name of each facies
    mutable VectorInt _faccols; // Color of each facies (packed RGB coding)
    mutable VectorInt _facvalues; // Value attached to each facies
  };

  GSTLEARN_EXPORT void set_rule_mode(Id rule_mode);
  GSTLEARN_EXPORT Id get_rule_mode(void);
  GSTLEARN_EXPORT double get_rule_extreme(Id mode);
  GSTLEARN_EXPORT Model* model_rule_combine(
    const Model* model1,
    const Model* model2,
    const Rule* rule);
  GSTLEARN_EXPORT Id db_rule_shadow(
    Db* db,
    Db* dbprop,
    RuleShadow* rule,
    Model* model1,
    const VectorDouble& props,
    bool flag_stat,
    Id nfacies);
} // namespace gstlrn
