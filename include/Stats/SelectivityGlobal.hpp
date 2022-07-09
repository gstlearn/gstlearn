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
#include "Stats/Selectivity.hpp"
#include "Basic/Vector.hpp"

class Db;

class GSTLEARN_EXPORT SelectivityGlobal: public Selectivity
{
public:
  SelectivityGlobal(int nclass = 0);
  SelectivityGlobal(const SelectivityGlobal &m);
  SelectivityGlobal& operator= (const SelectivityGlobal &m);
  virtual ~SelectivityGlobal();

  static SelectivityGlobal* createFromDb(const VectorDouble& zcuts, const Db* db);
  static SelectivityGlobal* createFromTab(const VectorDouble& zcuts,
                                    const VectorDouble& tab,
                                    const VectorDouble& weights = VectorDouble());
  static SelectivityGlobal* createInterpolation(const VectorDouble& zcuts,
                                          const SelectivityGlobal& calest,
                                          bool verbose);

  void setBest(int iclass, double best);
  void setMest(int iclass, double mest);
  void setQest(int iclass, double qest);
  void setQstd(int iclass, double qstd);
  void setTest(int iclass, double test);
  void setTstd(int iclass, double tstd);

  double getBest(int iclass) const;
  double getMest(int iclass) const;
  double getQest(int iclass) const;
  double getQstd(int iclass) const;
  double getTest(int iclass) const;
  double getTstd(int iclass) const;

  void calculateBenefitAndGrade();
  void dumpGini();
  void correctTonnageOrder();

private:
  bool _isValid(int iclass) const;
  void _interpolateInterval(double zval,
                            double zi0,
                            double zi1,
                            double ti0,
                            double ti1,
                            double qi0,
                            double qi1,
                            double *tval,
                            double *qval,
                            double tol = EPSILON3);

private:
  VectorDouble _Test;
  VectorDouble _Qest;
  VectorDouble _Best;
  VectorDouble _Mest;
  VectorDouble _Tstd;
  VectorDouble _Qstd;
};
