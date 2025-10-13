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

#include "Basic/VectorNumT.hpp"
#include "Basic/AStringable.hpp"


namespace gstlrn
{
class GSTLEARN_EXPORT Node: public AStringable
{
class Db;
class Model;
public:
  Node(const String& nodnam, Id orient, Id facies);
  Node(const String& nodnam,
       const VectorInt& n_type,
       const VectorInt& n_facs,
       Id *ipos,
       Id *n_fac,
       Id *n_y1,
       Id *n_y2);
  Node(bool flagShadow = true);
  Node(const Node& m);
  Node& operator=(const Node& m);
  virtual ~Node();

  String toString(const AStringFormat* strfmt = nullptr) const override;

  void getStatistics(Id *node_tot,
                     Id *nfac_tot,
                     Id *ny1_tot,
                     Id *ny2_tot,
                     double *prop_tot);
  Id  isValid(VectorInt& facies);
  void scaleProp(double scale);
  Id  proportionDefine(const VectorDouble& props);
  Id  getProportion(Id facies, double *prop);
  Id getThresh(Id mode,
                Id istop,
                Id *rank,
                Id *facies,
                double *t1min,
                double *t1max,
                double *t2min,
                double *t2max);
  void proportionToThresh(double rho,
                          double t1min,
                          double t1max,
                          double t2min,
                          double t2max);
  Id  gaussianToFacies(double y1, double y2, double *facies);
  void getInfo(Id *nodes) const;

  String nodePrint(bool flagProp, bool flagThresh) const;
  String nodePrintShadow(bool flagProp, bool flagThresh) const;

  double getT1max() const         { return _t1max; }
  double getT1min() const         { return _t1min; }
  double getT2max() const         { return _t2max; }
  double getT2min() const         { return _t2min; }
  void setT1max(double t1max)     { _t1max = t1max; }
  void setT1min(double t1min)     { _t1min = t1min; }
  void setT2max(double t2max)     { _t2max = t2max; }
  void setT2min(double t2min)     { _t2min = t2min; }
  double getCdf1max() const       { return _cdf1max; }
  double getCdf1min() const       { return _cdf1min; }
  double getCdf2max() const       { return _cdf2max; }
  double getCdf2min() const       { return _cdf2min; }
  Id getFacies() const           { return _facies; }
  const String& getNodnam() const { return _nodnam; }
  Id getOrient() const           { return _orient; }
  double getP1() const            { return _p1; }
  double getP2() const            { return _p2; }
  double getProp() const          { return _prop; }
  double getAllThresh() const     { return _thresh; }
  void setProp(double prop)       { _prop = prop; }
  void setCdf1max(double cdf1max) { _cdf1max = cdf1max; }
  void setCdf1min(double cdf1min) { _cdf1min = cdf1min; }
  void setCdf2max(double cdf2max) { _cdf2max = cdf2max; }
  void setCdf2min(double cdf2min) { _cdf2min = cdf2min; }
  void setAllThresh(double thresh){ _thresh = thresh; }
  Node* getR1() const { return _r1; }
  void setR1(Node* r1) { _r1 = r1; }
  Node* getR2() const { return _r2; }
  void setR2(Node* r2) { _r2 = r2; }

private:
  void _getStatistics(Id *node_tot,
                      Id *nfac_tot,
                      Id *ny1_tot,
                      Id *ny2_tot,
                      double *prop_tot);
  void _getInfo(Id *nodes,
                Id parent_type,
                Id parent_rank,
                Id parent_vers,
                Id *rank,
                Id *n_fac,
                Id *n_y1,
                Id *n_y2) const;
  static double _transform(Id mode, double value);
  double _threshFromPropcum(double rho);
  double _threshDichotomy(double rho) const;

private:
  String _nodnam;  /* Name of the node */
  Node  *_r1;      /* Pointer to the left-side */
  Node  *_r2;      /* Pointer to the right-side */
  Id    _orient;  /* Orientation */
  Id    _facies;  /* Facies number */
  double _prop;    /* Proportion */
  double _thresh;  /* Threshold value */
  double _p1;      /* Cumulative proportion on the left-side */
  double _p2;      /* Cumulative proportion on the right-side */
  double _t1min;   /* Lower bound along first gaussian */
  double _t1max;   /* Upper bound along first gaussian */
  double _t2min;   /* Lower bound along second gaussian */
  double _t2max;   /* Upper bound along second gaussian */
  double _cdf1min; /* CDF for Lower bound along first gaussian */
  double _cdf1max; /* CDF for Upper bound along first gaussian */
  double _cdf2min; /* CDF for Lower bound along second gaussian */
  double _cdf2max; /* CDF for Upper bound along second gaussian */
};
}