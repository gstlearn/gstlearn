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

#include "Basic/VectorNumT.hpp"
#include "Basic/AStringable.hpp"

class Db;
class Model;

class GSTLEARN_EXPORT Node: public AStringable
{
public:
  Node(const String& nodnam, int orient, int facies);
  Node(const String& nodnam,
       const VectorInt& n_type,
       const VectorInt& n_facs,
       int *ipos,
       int *n_fac,
       int *n_y1,
       int *n_y2);
  Node(bool flagShadow = true);
  Node(const Node& r);
  Node& operator=(const Node& r);
  virtual ~Node();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void getStatistics(int *node_tot,
                     int *nfac_tot,
                     int *ny1_tot,
                     int *ny2_tot,
                     double *prop_tot);
  int  isValid(VectorInt& facies);
  void scaleProp(double scale);
  int  proportionDefine(const VectorDouble& props);
  int  getProportion(int facies, double *prop);
  int getThresh(int mode,
                int istop,
                int *rank,
                int *facies,
                double *t1min,
                double *t1max,
                double *t2min,
                double *t2max);
  void proportionToThresh(double rho,
                          double t1min,
                          double t1max,
                          double t2min,
                          double t2max);
  int  gaussianToFacies(double y1, double y2, double *facies);
  void getInfo(int *nodes) const;

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
  int getFacies() const           { return _facies; }
  const String& getNodnam() const { return _nodnam; }
  int getOrient() const           { return _orient; }
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
  void _getStatistics(int *node_tot,
                      int *nfac_tot,
                      int *ny1_tot,
                      int *ny2_tot,
                      double *prop_tot);
  void _getInfo(int *nodes,
                int parent_type,
                int parent_rank,
                int parent_vers,
                int *rank,
                int *n_fac,
                int *n_y1,
                int *n_y2) const;
  double _transform(int mode, double value);
  double _threshFromPropcum(double rho);
  double _threshDichotomy(double rho);

private:
  String _nodnam;  /* Name of the node */
  Node  *_r1;      /* Pointer to the left-side */
  Node  *_r2;      /* Pointer to the right-side */
  int    _orient;  /* Orientation */
  int    _facies;  /* Facies number */
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
